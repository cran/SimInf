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

#include "SimInf.h"
#include "SimInf_forward_euler_linear_decay.h"

/* Offset in integer compartment state vector */
enum {S, I};

/* Offset in real-valued continuous state vector */
enum {PHI};

/* Offsets in node local data (ldata) to parameters in the model */
enum {END_T1, END_T2, END_T3, END_T4};

/* Offsets in global data (gdata) to parameters in the model */
enum {UPSILON, GAMMA, ALPHA, BETA_T1, BETA_T2, BETA_T3, BETA_T4, EPSILON};

/**
 * susceptible to infected: S -> I
 *
 * @param u The compartment state vector in node.
 * @param v The continuous state vector in node.
 * @param ldata The local data vector for the node.
 * @param gdata The global data vector.
 * @param t Current time.
 * @return propensity.
 */
double SISe_S_to_I(
    const int *u,
    const double *v,
    const double *ldata,
    const double *gdata,
    double t)
{
    return gdata[UPSILON] * v[PHI] * u[S];
}

/**
 *  infected to susceptible: I -> S
 *
 * @param u The compartment state vector in node.
 * @param v The continuous state vector in node.
 * @param ldata The local data vector for node.
 * @param gdata The global data vector.
 * @param t Current time.
 * @return propensity.
 */
double SISe_I_to_S(
    const int *u,
    const double *v,
    const double *ldata,
    const double *gdata,
    double t)
{
    return gdata[GAMMA] * u[I];
}

/**
 * Update environmental infectious pressure phi
 *
 * @param v_new The continuous state vector in the node after the post
 * time step
 * @param u The compartment state vector in the node.
 * @param v The current continuous state vector in the node.
 * @param ldata The local data vector for the node.
 * @param gdata The global data vector.
 * @param node The node.
 * @param t The current time.
 * @param rng The random number generator.
 * @return error code (<0), or 1 if node needs to update the
 * transition rates, or 0 when it doesn't need to update the
 * transition rates.
 */
int SISe_post_time_step(
    double *v_new,
    const int *u,
    const double *v,
    const double *ldata,
    const double *gdata,
    int node,
    double t)
{
    const int day = (int)t % 365;
    const double I_n = u[I];
    const double n = u[S] + I_n;
    const double phi = v[PHI];

    /* Time dependent beta in each of the four intervals of the
     * year. Forward Euler step. */
    v_new[PHI] = SimInf_forward_euler_linear_decay(
        phi, day,
        ldata[END_T1], ldata[END_T2], ldata[END_T3], ldata[END_T4],
        gdata[BETA_T1], gdata[BETA_T2], gdata[BETA_T3], gdata[BETA_T4]);

    if (n > 0.0)
        v_new[PHI] += gdata[ALPHA] * I_n / n + gdata[EPSILON];
    else
        v_new[PHI] += gdata[EPSILON];

    if (!isfinite(v_new[PHI]))
        return SIMINF_ERR_V_IS_NOT_FINITE;
    if (v_new[PHI] < 0.0)
        return SIMINF_ERR_V_IS_NEGATIVE;
    return phi != v_new[PHI]; /* 1 if needs update */
}

/**
 * Run simulation with the SISe model
 *
 * @param model The SISe model.
 * @param threads Number of threads.
 * @param solver The numerical solver.
 * @return The simulated trajectory.
 */
SEXP SISe_run(SEXP model, SEXP threads, SEXP solver)
{
    TRFun tr_fun[] = {&SISe_S_to_I, &SISe_I_to_S};

    return SimInf_run(model, threads, solver, tr_fun, &SISe_post_time_step);
}
