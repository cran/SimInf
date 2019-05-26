/*
 *  SimInf, a framework for stochastic disease spread simulations
 *  Copyright (C) 2015 Pavol Bauer
 *  Copyright (C) 2017 - 2019 Robin Eriksson
 *  Copyright (C) 2015 - 2019 Stefan Engblom
 *  Copyright (C) 2015 - 2019 Stefan Widgren
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
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <Rdefines.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "SimInf_arg.h"
#include "solvers/SimInf_solver.h"
#include "solvers/ssm/SimInf_solver_ssm.h"
#include "solvers/aem/SimInf_solver_aem.h"

/**
 * Initiate and run the simulation
 *
 * @param model The SimInf_model
 * @param threads Number of threads
 * @param solver The numerical solver.
 * @param tr_fun Vector of function pointers to transition rate functions.
 * @param pts_fun Function pointer to callback after each time step
 *        e.g. update infectious pressure.
 */
SEXP SimInf_run(
    SEXP model,
    SEXP threads,
    SEXP solver,
    TRFun *tr_fun,
    PTSFun pts_fun)
{
    int i, j, error = 0, nprotect = 0;
    SEXP result = R_NilValue;
    SEXP ext_events, E, G, N, S, prS;
    SEXP tspan, rownames, colnames;
    SEXP U_dimnames, U_rownames, V_dimnames;
    SEXP U, V, U_sparse, V_sparse;
    SimInf_solver_args args = {NULL};

    /* If the model ldata is a 0x0 matrix, i.e. Nld == 0, then use
     * ldata_tmp in the transition rate functions. This is to make
     * &ldata[node * Nld] work in the solvers. The reason for INFINITY
     * is to facilitate for the solvers to detect and raise an error
     * if a model C code uses ldata[0] in the transition rate
     * functions. */
    const double ldata_tmp[1] = {INFINITY};

    if (SimInf_arg_check_model(model)) {
        error = SIMINF_ERR_INVALID_MODEL;
        goto cleanup;
    }

    /* Check solver argument */
    if (!Rf_isNull(solver)) {
        if (!Rf_isString(solver)) {
            error = SIMINF_ERR_UNKNOWN_SOLVER;
            goto cleanup;
        }

        if (Rf_length(solver) != 1 || STRING_ELT(solver, 0) == NA_STRING) {
            error = SIMINF_ERR_UNKNOWN_SOLVER;
            goto cleanup;
        }
    }

    /* number of threads */
    error = SimInf_get_threads(&(args.Nthread), threads);
    if (error)
        goto cleanup;

    /* seed */
    GetRNGstate();
    args.seed = (unsigned long int)(unif_rand() * UINT_MAX);
    PutRNGstate();

    /* Duplicate model. */
    PROTECT(result = Rf_duplicate(model));
    nprotect++;

    /* Dependency graph */
    PROTECT(G = GET_SLOT(result, Rf_install("G")));
    nprotect++;
    args.irG = INTEGER(GET_SLOT(G, Rf_install("i")));
    args.jcG = INTEGER(GET_SLOT(G, Rf_install("p")));

    /* State change matrix */
    PROTECT(S = GET_SLOT(result, Rf_install("S")));
    nprotect++;
    PROTECT(prS = Rf_coerceVector(GET_SLOT(S, Rf_install("x")), INTSXP));
    nprotect++;
    args.irS = INTEGER(GET_SLOT(S, Rf_install("i")));
    args.jcS = INTEGER(GET_SLOT(S, Rf_install("p")));
    args.prS = INTEGER(prS);

    /* tspan */
    PROTECT(tspan = GET_SLOT(result, Rf_install("tspan")));
    nprotect++;
    args.tspan = REAL(tspan);

    /* Dimnames */
    rownames = VECTOR_ELT(GET_SLOT(S, Rf_install("Dimnames")), 0);
    PROTECT(colnames = Rf_getAttrib(tspan , R_NamesSymbol));
    nprotect++;

    /* Scheduled events */
    PROTECT(ext_events = GET_SLOT(result, Rf_install("events")));
    nprotect++;
    args.len = LENGTH(GET_SLOT(ext_events, Rf_install("event")));
    args.event = INTEGER(GET_SLOT(ext_events, Rf_install("event")));
    args.time = INTEGER(GET_SLOT(ext_events, Rf_install("time")));
    args.node = INTEGER(GET_SLOT(ext_events, Rf_install("node")));
    args.dest = INTEGER(GET_SLOT(ext_events, Rf_install("dest")));
    args.n = INTEGER(GET_SLOT(ext_events, Rf_install("n")));
    args.proportion = REAL(GET_SLOT(ext_events, Rf_install("proportion")));
    args.select = INTEGER(GET_SLOT(ext_events, Rf_install("select")));
    args.shift = INTEGER(GET_SLOT(ext_events, Rf_install("shift")));
    PROTECT(E = GET_SLOT(ext_events, Rf_install("E")));
    nprotect++;
    args.irE = INTEGER(GET_SLOT(E, Rf_install("i")));
    args.jcE = INTEGER(GET_SLOT(E, Rf_install("p")));
    PROTECT(N = GET_SLOT(ext_events, Rf_install("N")));
    nprotect++;
    args.N = INTEGER(N);

    /* Constants */
    args.Nn = INTEGER(GET_SLOT(GET_SLOT(result, Rf_install("u0")), R_DimSymbol))[1];
    args.Nc = INTEGER(GET_SLOT(S, Rf_install("Dim")))[0];
    args.Nt = INTEGER(GET_SLOT(S, Rf_install("Dim")))[1];
    args.Nd = INTEGER(GET_SLOT(GET_SLOT(result, Rf_install("v0")), R_DimSymbol))[0];
    args.Nld = INTEGER(GET_SLOT(GET_SLOT(result, Rf_install("ldata")), R_DimSymbol))[0];
    args.tlen = LENGTH(GET_SLOT(result, Rf_install("tspan")));

    /* Output array (to hold a single trajectory) */
    PROTECT(U_sparse = GET_SLOT(result, Rf_install("U_sparse")));
    nprotect++;
    if ((INTEGER(GET_SLOT(U_sparse, Rf_install("Dim")))[0] == (args.Nn * args.Nc)) &&
        (INTEGER(GET_SLOT(U_sparse, Rf_install("Dim")))[1] == args.tlen))
    {
        args.irU = INTEGER(GET_SLOT(U_sparse, Rf_install("i")));
        args.jcU = INTEGER(GET_SLOT(U_sparse, Rf_install("p")));
        args.prU = REAL(GET_SLOT(U_sparse, Rf_install("x")));

        PROTECT(U_dimnames = GET_SLOT(U_sparse, Rf_install("Dimnames")));
        nprotect++;
        PROTECT(U_rownames = Rf_allocVector(STRSXP, args.Nn * args.Nc));
        nprotect++;
        SET_VECTOR_ELT(U_dimnames, 0, U_rownames);
    } else {
        PROTECT(U = Rf_allocMatrix(INTSXP, args.Nn * args.Nc, args.tlen));
        nprotect++;
        SET_SLOT(result, Rf_install("U"), U);
        args.U = INTEGER(GET_SLOT(result, Rf_install("U")));

        PROTECT(U_dimnames = Rf_allocVector(VECSXP, 2));
        nprotect++;
        Rf_setAttrib(GET_SLOT(result, Rf_install("U")),
                     R_DimNamesSymbol, U_dimnames);
        PROTECT(U_rownames = Rf_allocVector(STRSXP, args.Nn * args.Nc));
        nprotect++;
        SET_VECTOR_ELT(U_dimnames, 0, U_rownames);
    }

    /* Add rownames to U */
    for (i = 0; i < args.Nn; i++)
        for (j = 0; j < args.Nc; j++)
            SET_STRING_ELT(U_rownames, i * args.Nc + j, STRING_ELT(rownames, j));

    /* Add colnames to U. Use the the values of 'tspan' if the
     * colnames of 'tspan' is null. */
    if (Rf_isNull(colnames))
        SET_VECTOR_ELT(U_dimnames, 1, Rf_coerceVector(tspan, STRSXP));
    else
        SET_VECTOR_ELT(U_dimnames, 1, Rf_duplicate(colnames));

    /* Output array (to hold a single trajectory) */
    PROTECT(V_sparse = GET_SLOT(result, Rf_install("V_sparse")));
    nprotect++;
    if ((INTEGER(GET_SLOT(V_sparse, Rf_install("Dim")))[0] == (args.Nn * args.Nd)) &&
        (INTEGER(GET_SLOT(V_sparse, Rf_install("Dim")))[1] == args.tlen))
    {
        args.irV = INTEGER(GET_SLOT(V_sparse, Rf_install("i")));
        args.jcV = INTEGER(GET_SLOT(V_sparse, Rf_install("p")));
        args.prV = REAL(GET_SLOT(V_sparse, Rf_install("x")));

        V_dimnames = GET_SLOT(V_sparse, Rf_install("Dimnames"));
    } else {
        PROTECT(V = Rf_allocMatrix(REALSXP, args.Nn * args.Nd, args.tlen));
        nprotect++;
        SET_SLOT(result, Rf_install("V"), V);
        args.V = REAL(GET_SLOT(result, Rf_install("V")));

        PROTECT(V_dimnames = Rf_allocVector(VECSXP, 2));
        nprotect++;
        Rf_setAttrib(GET_SLOT(result, Rf_install("V")),
                     R_DimNamesSymbol, V_dimnames);
    }

    /* Add colnames to V. Use the the values of 'tspan' if the
     * colnames of 'tspan' is null. */
    if (Rf_isNull(colnames))
        SET_VECTOR_ELT(V_dimnames, 1, Rf_coerceVector(tspan, STRSXP));
    else
        SET_VECTOR_ELT(V_dimnames, 1, Rf_duplicate(colnames));

    /* Initial state. Copy u0 to either U[, 1] or U_sparse[, 1] */
    args.u0 = INTEGER(GET_SLOT(result, Rf_install("u0")));
    if (args.U) {
        memcpy(args.U, args.u0, args.Nn * args.Nc * sizeof(int));
    } else {
        for (i = args.jcU[0]; i < args.jcU[1]; i++)
            args.prU[i] = args.u0[args.irU[i]];
    }

    /* Initial state. Copy v0 to either V[, 1] or V_sparse[, 1] */
    args.v0 = REAL(GET_SLOT(result, Rf_install("v0")));
    if (args.V) {
        memcpy(args.V, args.v0, args.Nn * args.Nd * sizeof(double));
    } else {
        for (i = args.jcV[0]; i < args.jcV[1]; i++)
            args.prV[i] = args.v0[args.irV[i]];
    }

    /* Local data */
    if (args.Nld > 0)
        args.ldata = REAL(GET_SLOT(result, Rf_install("ldata")));
    else
        args.ldata = ldata_tmp;

    /* Global data */
    args.gdata = REAL(GET_SLOT(result, Rf_install("gdata")));

    /* Function pointers */
    args.tr_fun = tr_fun;
    args.pts_fun = pts_fun;

    /* Specify the number of threads to use. Make sure to not use more
     * threads than the number of nodes in the model. */
#ifdef _OPENMP
    if (args.Nthread < 1)
        args.Nthread = omp_get_num_procs();
#else
    args.Nthread = 1;
#endif
    if (args.Nn < args.Nthread)
        args.Nthread = args.Nn;
#ifdef _OPENMP
    omp_set_num_threads(args.Nthread);
#endif

    /* Run the simulation solver. */
    if (Rf_isNull(solver) || (strcmp(CHAR(STRING_ELT(solver, 0)), "ssm") == 0))
        error = SimInf_run_solver_ssm(&args);
    else if (strcmp(CHAR(STRING_ELT(solver, 0)), "aem") == 0)
        error = SimInf_run_solver_aem(&args);
    else
        error = SIMINF_ERR_UNKNOWN_SOLVER;

cleanup:
    if (error) {
        switch (error) {
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
        case SIMINF_ERR_UNKNOWN_SOLVER:
            Rf_error("Invalid 'solver' value.");
            break;
        case SIMINF_ERR_DEST_OUT_OF_BOUNDS:
            Rf_error("'dest' is out of bounds.");
            break;
        case SIMINF_ERR_NODE_OUT_OF_BOUNDS:
            Rf_error("'node' is out of bounds.");
            break;
        case SIMINF_ERR_EVENTS_N:
            Rf_error("'N' is invalid.");
            break;
        case SIMINF_ERR_EVENT_SHIFT:
            Rf_error("'shift' is invalid.");
            break;
        case SIMINF_ERR_SHIFT_OUT_OF_BOUNDS:
            Rf_error("'shift' is out of bounds.");
            break;
        default:
            Rf_error("Unknown error code: %i", error);
            break;
        }
    }

    if (nprotect)
        UNPROTECT(nprotect);

    return result;
}
