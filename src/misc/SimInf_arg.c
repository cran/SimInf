/*
 *  SimInf, a framework for stochastic disease spread simulations
 *  Copyright (C) 2015 Pavol Bauer
 *  Copyright (C) 2017 - 2018 Robin Eriksson
 *  Copyright (C) 2015 - 2018 Stefan Engblom
 *  Copyright (C) 2015 - 2018 Stefan Widgren
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

#include "SimInf.h"

/**
 * Check dgCMatrix argument
 *
 * @param arg The arg to check
 * @return 0 if OK, else -1
 */
int SimInf_arg_check_dgCMatrix(SEXP arg)
{
    SEXP class_name;

    if (!Rf_isS4(arg))
        return -1;
    class_name = Rf_getAttrib(arg, R_ClassSymbol);
    if (0 != strcmp(CHAR(STRING_ELT(class_name, 0)), "dgCMatrix"))
        return -1;
    return 0;
}

/**
 * Check integer argument
 *
 * @param arg The arg to check
 * @return 0 if OK, else -1
 */
int SimInf_arg_check_integer(SEXP arg)
{
    if (!Rf_isInteger(arg) || Rf_length(arg) != 1 || NA_INTEGER == INTEGER(arg)[0])
        return -1;
    return 0;
}

/**
 * Check matrix argument
 *
 * @param arg The arg to check
 * @return 0 if OK, else -1
 */
int SimInf_arg_check_matrix(SEXP arg)
{
    if (!Rf_isMatrix(arg))
        return -1;
    return 0;
}

/**
 * Check model argument
 *
 * @param arg The arg to check
 * @return 0 if OK, else -1
 */
int SimInf_arg_check_model(SEXP arg)
{
    static const char *valid[] = {"SimInf_model", ""};

    if (!Rf_isS4(arg) || R_check_class_etc(arg, valid) < 0)
        return -1;

    return 0;
}

/**
 * Get number of threads
 *
 * @param out Number of threads
 * @param threads Number of threads from R
 * @return 0 if Ok, else error code.
 */
int SimInf_get_threads(int *out, SEXP threads)
{
    int error = 0;

    if (Rf_isNull(threads)) {
        *out = 0;
    } else if (Rf_isInteger(threads)) {
        if (LENGTH(threads) != 1)
            error = SIMINF_INVALID_THREADS_VALUE;
        else if (INTEGER(threads)[0] == NA_INTEGER)
            error = SIMINF_INVALID_THREADS_VALUE;
        else if (INTEGER(threads)[0] < 0)
            error = SIMINF_INVALID_THREADS_VALUE;
        else
            *out = INTEGER(threads)[0];
    } else if (Rf_isReal(threads)) {
        if (LENGTH(threads) != 1)
            error = SIMINF_INVALID_THREADS_VALUE;
        else if (!R_finite(REAL(threads)[0]))
            error = SIMINF_INVALID_THREADS_VALUE;
        else if ((int)(REAL(threads)[0] < 0))
            error = SIMINF_INVALID_THREADS_VALUE;
        else
            *out = (int)(REAL(threads)[0]);
    } else {
        error = SIMINF_INVALID_THREADS_VALUE;
    }

    return error;
}
