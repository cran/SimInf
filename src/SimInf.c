/*
 *  SimInf, a framework for stochastic disease spread simulations
 *  Copyright (C) 2015  Pavol Bauer
 *  Copyright (C) 2015 - 2017  Stefan Engblom
 *  Copyright (C) 2015 - 2017  Stefan Widgren
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <Rdefines.h>

#include "SimInf_arg.h"
#include "SimInf_solver.h"

/**
 * Initiate and run the simulation
 *
 * @param model The siminf_model
 * @param threads Number of threads
 * @param seed Random number seed.
 * @param tr_fun Vector of function pointers to transition rate functions.
 * @param pts_fun Function pointer to callback after each time step
 *        e.g. update infectious pressure.
 */
SEXP SimInf_run(
    SEXP model,
    SEXP threads,
    SEXP seed,
    TRFun *tr_fun,
    PTSFun pts_fun)
{
    int i, j, err = 0, nprotect = 0, n_threads;
    SEXP result = R_NilValue;
    SEXP ext_events, E, G, N, S, prS;
    SEXP tspan, rownames, colnames;
    SEXP U_dimnames, U_rownames, V_dimnames;
    SEXP U, V, U_sparse, V_sparse;
    int *U_ptr = NULL, *irU = NULL, *jcU = NULL;
    double *prU = NULL;
    int *irV = NULL, *jcV = NULL;
    double *V_ptr = NULL, *prV = NULL;
    int Nn, Nc, Nt, Nd, Nld, tlen;
    unsigned long int s;

    if (SimInf_arg_check_model(model)) {
        err = SIMINF_ERR_INVALID_MODEL;
        goto cleanup;
    }

    /* number of threads */
    err = SimInf_get_threads(&n_threads, threads);
    if (err)
        goto cleanup;

    /* seed */
    err =  SimInf_get_seed(&s, seed);
    if (err)
        goto cleanup;

    /* Duplicate model. */
    PROTECT(result = duplicate(model));
    nprotect++;

    /* SimInf model */
    PROTECT(G = GET_SLOT(result, Rf_install("G")));
    nprotect++;
    PROTECT(S = GET_SLOT(result, Rf_install("S")));
    nprotect++;
    PROTECT(prS = coerceVector(GET_SLOT(S, Rf_install("x")), INTSXP));
    nprotect++;

    /* Dimnames */
    rownames = VECTOR_ELT(GET_SLOT(S, Rf_install("Dimnames")), 0);
    PROTECT(tspan = GET_SLOT(result, Rf_install("tspan")));
    nprotect++;
    PROTECT(colnames = Rf_getAttrib(tspan , R_NamesSymbol));
    nprotect++;

    /* Scheduled events */
    PROTECT(ext_events = GET_SLOT(result, Rf_install("events")));
    nprotect++;
    PROTECT(E = GET_SLOT(ext_events, Rf_install("E")));
    nprotect++;
    PROTECT(N = GET_SLOT(ext_events, Rf_install("N")));
    nprotect++;

    /* Constants */
    Nn   = INTEGER(GET_SLOT(GET_SLOT(result, Rf_install("u0")), R_DimSymbol))[1];
    Nc   = INTEGER(GET_SLOT(S, Rf_install("Dim")))[0];
    Nt   = INTEGER(GET_SLOT(S, Rf_install("Dim")))[1];
    Nd   = INTEGER(GET_SLOT(GET_SLOT(result, Rf_install("v0")), R_DimSymbol))[0];
    Nld  = INTEGER(GET_SLOT(GET_SLOT(result, Rf_install("ldata")), R_DimSymbol))[0];
    tlen = LENGTH(GET_SLOT(result, Rf_install("tspan")));

    /* Output array (to hold a single trajectory) */
    PROTECT(U_sparse = GET_SLOT(result, Rf_install("U_sparse")));
    nprotect++;
    if ((INTEGER(GET_SLOT(U_sparse, Rf_install("Dim")))[0] == (Nn * Nc)) &&
        (INTEGER(GET_SLOT(U_sparse, Rf_install("Dim")))[1] == tlen))
    {
        irU = INTEGER(GET_SLOT(U_sparse, Rf_install("i")));
        jcU = INTEGER(GET_SLOT(U_sparse, Rf_install("p")));
        prU = REAL(GET_SLOT(U_sparse, Rf_install("x")));

        PROTECT(U_dimnames = GET_SLOT(U_sparse, Rf_install("Dimnames")));
        nprotect++;
        PROTECT(U_rownames = allocVector(STRSXP, Nn * Nc));
        nprotect++;
        SET_VECTOR_ELT(U_dimnames, 0, U_rownames);
    } else {
        PROTECT(U = allocMatrix(INTSXP, Nn * Nc, tlen));
        nprotect++;
        SET_SLOT(result, Rf_install("U"), U);
        U_ptr = INTEGER(GET_SLOT(result, Rf_install("U")));

        PROTECT(U_dimnames = allocVector(VECSXP, 2));
        nprotect++;
        setAttrib(GET_SLOT(result, Rf_install("U")),
                  R_DimNamesSymbol, U_dimnames);
        PROTECT(U_rownames = allocVector(STRSXP, Nn * Nc));
        nprotect++;
        SET_VECTOR_ELT(U_dimnames, 0, U_rownames);
    }

    /* Add rownames to U */
    for (i = 0; i < Nn; i++)
        for (j = 0; j < Nc; j++)
            SET_STRING_ELT(U_rownames, i * Nc + j, STRING_ELT(rownames, j));

    /* Add colnames to U. Use the the values of 'tspan' if the
     * colnames of 'tspan' is null. */
    if (isNull(colnames))
        SET_VECTOR_ELT(U_dimnames, 1, coerceVector(tspan, STRSXP));
    else
        SET_VECTOR_ELT(U_dimnames, 1, duplicate(colnames));

    /* Output array (to hold a single trajectory) */
    PROTECT(V_sparse = GET_SLOT(result, Rf_install("V_sparse")));
    nprotect++;
    if ((INTEGER(GET_SLOT(V_sparse, Rf_install("Dim")))[0] == (Nn * Nd)) &&
        (INTEGER(GET_SLOT(V_sparse, Rf_install("Dim")))[1] == tlen))
    {
        irV = INTEGER(GET_SLOT(V_sparse, Rf_install("i")));
        jcV = INTEGER(GET_SLOT(V_sparse, Rf_install("p")));
        prV = REAL(GET_SLOT(V_sparse, Rf_install("x")));

        V_dimnames = GET_SLOT(V_sparse, Rf_install("Dimnames"));
    } else {
        PROTECT(V = allocMatrix(REALSXP, Nn * Nd, tlen));
        nprotect++;
        SET_SLOT(result, Rf_install("V"), V);
        V_ptr = REAL(GET_SLOT(result, Rf_install("V")));

        PROTECT(V_dimnames = allocVector(VECSXP, 2));
        nprotect++;
        setAttrib(GET_SLOT(result, Rf_install("V")),
                  R_DimNamesSymbol, V_dimnames);
    }

    /* Add colnames to V. Use the the values of 'tspan' if the
     * colnames of 'tspan' is null. */
    if (isNull(colnames))
        SET_VECTOR_ELT(V_dimnames, 1, coerceVector(tspan, STRSXP));
    else
        SET_VECTOR_ELT(V_dimnames, 1, duplicate(colnames));

    /* Run simulation solver. */
    err = SimInf_run_solver(
        INTEGER(GET_SLOT(result, Rf_install("u0"))),
        REAL(GET_SLOT(result, Rf_install("v0"))),
        INTEGER(GET_SLOT(G, Rf_install("i"))),
        INTEGER(GET_SLOT(G, Rf_install("p"))),
        INTEGER(GET_SLOT(S, Rf_install("i"))),
        INTEGER(GET_SLOT(S, Rf_install("p"))),
        INTEGER(prS),
        REAL(GET_SLOT(result, Rf_install("tspan"))),
        tlen,
        U_ptr, irU, jcU, prU,
        V_ptr, irV, jcV, prV,
        REAL(GET_SLOT(result, Rf_install("ldata"))),
        REAL(GET_SLOT(result, Rf_install("gdata"))),
        Nn, Nc, Nt, Nd, Nld,
        INTEGER(GET_SLOT(E, Rf_install("i"))),
        INTEGER(GET_SLOT(E, Rf_install("p"))),
        INTEGER(N),
        LENGTH(GET_SLOT(ext_events, Rf_install("event"))),
        INTEGER(GET_SLOT(ext_events, Rf_install("event"))),
        INTEGER(GET_SLOT(ext_events, Rf_install("time"))),
        INTEGER(GET_SLOT(ext_events, Rf_install("node"))),
        INTEGER(GET_SLOT(ext_events, Rf_install("dest"))),
        INTEGER(GET_SLOT(ext_events, Rf_install("n"))),
        REAL(GET_SLOT(ext_events,    Rf_install("proportion"))),
        INTEGER(GET_SLOT(ext_events, Rf_install("select"))),
        INTEGER(GET_SLOT(ext_events, Rf_install("shift"))),
        n_threads, s, tr_fun, pts_fun);

cleanup:
    if (err) {
        switch (err) {
        case SIMINF_ERR_NEGATIVE_STATE:
            Rf_error("Negative state detected.");
            break;
        case SIMINF_ERR_ALLOC_MEMORY_BUFFER:
            Rf_error("Unable to allocate memory buffer");
            break;
        case SIMINF_UNDEFINED_EVENT:
            Rf_error("Undefined event type.");
            break;
        case SIMINF_INVALID_SEED_VALUE:
            Rf_error("Invalid 'seed' value");
            break;
        case SIMINF_INVALID_THREADS_VALUE:
            Rf_error("Invalid 'threads' value");
            break;
        case SIMINF_ERR_V_IS_NOT_FINITE:
            Rf_error("The continuous state 'v' is not finite.");
            break;
        case SIMINF_ERR_SAMPLE_SELECT:
            Rf_error("Unable to sample individuals for event.");
            break;
        case SIMINF_ERR_INVALID_MODEL:
            Rf_error("Invalid model.");
            break;
        case SIMINF_ERR_V_IS_NEGATIVE:
            Rf_error("The continuous state 'v' is negative.");
            break;
        case SIMINF_ERR_INVALID_RATE:
            Rf_error("Invalid rate detected (non-finite or < 0.0)");
            break;
        default:
            Rf_error("Unknown error code: %i", err);
            break;
        }
    }

    if (nprotect)
        UNPROTECT(nprotect);

    return result;
}
