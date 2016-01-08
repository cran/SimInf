/*
 *  siminf, a framework for stochastic disease spread simulations
 *  Copyright (C) 2015  Pavol Bauer
 *  Copyright (C) 2015  Stefan Engblom
 *  Copyright (C) 2015  Stefan Widgren
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
#include <time.h>

#include "siminf_error.h"

/**
 * Check integer argument
 *
 * @param arg The arg to check
 * @return 0 if OK, else -1
 */
int siminf_arg_check_integer(SEXP arg)
{
    if (arg == R_NilValue || !isInteger(arg) ||
        length(arg) != 1  || NA_INTEGER == INTEGER(arg)[0])
        return -1;
    return 0;
}

/**
 * Check vector argument of type real
 *
 * @param arg The arg to check
 * @param size The length of the vector
 * @return 0 if OK, else -1
 */
int siminf_arg_check_real_vec(SEXP arg, size_t size)
{
    size_t i = 0;

    if (arg == R_NilValue || !isReal(arg) || length(arg) != size)
        return -1;
    for (; i < size; i++) {
        if (!R_finite(REAL(arg)[i]))
            return -1;
    }

    return 0;
}

/**
 * Get seed value
 *
 * @param out The seed value.
 * @param seed Random number seed from R.
 * @return 0 if Ok, else error code.
 */
int siminf_get_seed(unsigned long int *out, SEXP seed)
{
    int err = 0;

    if (seed != R_NilValue) {
        if (isInteger(seed) || isReal(seed)) {
            switch (LENGTH(seed)) {
            case 0:
                *out = (unsigned long int)time(NULL);
                break;
            case 1:
                if (isInteger(seed)) {
                    if (INTEGER(seed)[0] == NA_INTEGER)
                        err = SIMINF_INVALID_SEED_VALUE;
                    else
                        *out = (unsigned long int)INTEGER(seed)[0];
                } else if (!R_finite(REAL(seed)[0])) {
                    err = SIMINF_INVALID_SEED_VALUE;
                } else {
                    *out = (unsigned long int)REAL(seed)[0];
                }
                break;
            default:
                err = SIMINF_INVALID_SEED_VALUE;
                break;
            }
        } else {
            err = SIMINF_INVALID_SEED_VALUE;
        }
    } else {
        *out = (unsigned long int)time(NULL);
    }

    return err;
}

/**
 * Get number of threads
 *
 * @param out Number of threads
 * @param threads Number of threads from R
 * @return 0 if Ok, else error code.
 */
int siminf_get_threads(int *out, SEXP threads)
{
    int err = 0;

    if (threads == R_NilValue) {
        *out = 0;
    } else if (isInteger(threads)) {
        if (LENGTH(threads) != 1)
            err = SIMINF_INVALID_THREADS_VALUE;
        else if (INTEGER(threads)[0] == NA_INTEGER)
            err = SIMINF_INVALID_THREADS_VALUE;
        else if (INTEGER(threads)[0] < 0)
            err = SIMINF_INVALID_THREADS_VALUE;
        else
            *out = INTEGER(threads)[0];
    } else if (isReal(threads)) {
        if (LENGTH(threads) != 1)
            err = SIMINF_INVALID_THREADS_VALUE;
        else if (!R_finite(REAL(threads)[0]))
            err = SIMINF_INVALID_THREADS_VALUE;
        else if ((int)(REAL(threads)[0] < 0))
            err = SIMINF_INVALID_THREADS_VALUE;
        else
            *out = (int)(REAL(threads)[0]);
    } else {
        err = SIMINF_INVALID_THREADS_VALUE;
    }

    return err;
}
