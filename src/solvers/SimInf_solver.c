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

#include <string.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "SimInf.h"
#include "SimInf_solver.h"

/**
 * Sample individuals from a node
 *
 * Individuals are sampled from the states determined by select.
 *
 * @param irE Select matrix for events. irE[k] is the row of E[k].
 * @param jcE Select matrix for events. Index to data of first
 *        non-zero element in row k.
 * @param Nc Number of compartments in each node.
 * @param u The state vector with number of individuals in each
 *        compartment at each node. The current state in each node is
 *        offset by node * Nc.
 * @param node The node to sample.
 * @param select Column j in the Select matrix that determines the
 *        states to sample from.
 * @param n The number of individuals to sample. n >= 0.
 * @param proportion If n equals zero, then the number of individuals
 *        to sample is calculated by summing the number of individuals
 *        in the states determined by select and multiplying with the
 *        proportion. 0 <= proportion <= 1.
 * @param individuals The result of the sampling is stored in the
 *        individuals vector.
 * @param u_tmp Help vector for sampling individuals.
 * @param rng Random number generator.
 * @return 0 if Ok, else error code.
 */
static int SimInf_sample_select(
    const int *irE, const int *jcE, int Nc, const int *u,
    int node, int select, int n, double proportion,
    int *individuals, int *u_tmp, gsl_rng *rng)
{
    int i, Nstates, Nindividuals = 0, Nkinds = 0;

    /* Clear vector with number of sampled individuals */
    memset(individuals, 0, Nc * sizeof(int));

    /* 1) Count number of states with individuals */
    /* 2) Count total number of individuals       */
    for (i = jcE[select]; i < jcE[select + 1]; i++) {
        const int nk = u[node * Nc + irE[i]];
        if (nk > 0)
            Nkinds++;
        Nindividuals += nk;
    }

    /* Number of states */
    Nstates = jcE[select + 1] - jcE[select];

    /* If n == 0, use the proportion of Nindividuals, else use n as */
    /* the number of individuals to sample                          */
    if (n == 0)
        n = round(proportion * Nindividuals);

    /* Error checking. */
    if (Nstates <= 0 ||     /* No states to sample from, we shouldn't be here. */
        n > Nindividuals || /* Cannot sample this number of individuals.       */
        n < 0)              /* Cannot sample negative number of individuals.   */
        return SIMINF_ERR_SAMPLE_SELECT;

    /* Handle cases that require no random sampling */
    if (n == 0) {
        /* We are done */
        return 0;
    } else if (Nindividuals == n) {
        /* Include all individuals */
        for (i = jcE[select]; i < jcE[select + 1]; i++)
            individuals[irE[i]] = u[node * Nc + irE[i]];
        return 0;
    } else if (Nstates == 1) {
        /* Only individuals from one state to select from. */
        individuals[irE[jcE[select]]] = n;
        return 0;
    } else if (Nkinds == 1) {
        /* All individuals to choose from in one state */
        for (i = jcE[select]; i < jcE[select + 1]; i++) {
            if (u[node * Nc + irE[i]] > 0) {
                individuals[irE[i]] = n;
                break;
            }
        }
        return 0;
    }

    /* Handle cases that require random sampling */
    if (Nstates == 2) {
        /* Sample from the hypergeometric distribution */
        i = jcE[select];
        individuals[irE[i]] = gsl_ran_hypergeometric(
            rng,
            u[node * Nc + irE[i]],
            u[node * Nc + irE[i+1]],
            n);
        individuals[irE[i+1]] = n - individuals[irE[i]];
    } else {
        /* Randomly sample n individuals from Nindividudals in
         * the Nstates */
        memcpy(u_tmp, &u[node * Nc], Nc * sizeof(int));
        while (n > 0) {
            double cum, rand = gsl_rng_uniform_pos(rng) * Nindividuals;

            /* Determine from which compartment the individual was
             * sampled from */
            for (i = jcE[select], cum = u_tmp[irE[i]];
                 i < jcE[select + 1] && rand > cum;
                 i++, cum += u_tmp[irE[i]]);

            /* Update sampled individual */
            u_tmp[irE[i]]--;
            individuals[irE[i]]++;

            Nindividuals--;
            n--;
        }
    }

    return 0;
}

/**
 * Split scheduled events to E1 and E2 events by number of threads
 * used during simulation
 *
 * Thread id 0 is the main thread. All E2 events are assigned to
 * thread id 0.
 *
 * All E1 events for a node are assigned to the same thread.
 *
 * @param len Number of scheduled events.
 * @param event The type of event i.
 * @param time The time of event i.
 * @param node The source node index (one based) of event i.
 * @param dest The dest node index (one-based) of event i.
 * @param n The number of individuals in event i. n[i] >= 0.
 * @param proportion If n[i] equals zero, then the number of
 *        individuals to sample is calculated by summing the number of
 *        individuals in the states determined by select[i] and
 *        multiplying with the proportion. 0 <= p[i] <= 1.
 * @param select Column j (one-based) in the event matrix that
 *        determines the states to sample from.
 * @param shift Column j (one-based) in the shift matrix S that
 *        determines the shift of the internal and external
 *        transfer event.
 * @param Nn Total number of nodes.
 * @param Nthread Number of threads to use during simulation.
 * @return 0 if Ok, else error code.
 */
static int SimInf_split_events(
    SimInf_scheduled_events *out,
    int len, const int *event, const int *time, const int *node,
    const int *dest, const int *n, const double *proportion,
    const int *select, const int *shift, int Nn, int Nthread)
{
    int i;
    int chunk_size = Nn / Nthread;

    for (i = 0; i < len; i++) {
        int j;
        const SimInf_scheduled_event e = {event[i], time[i], node[i] - 1,
                                          dest[i] - 1, n[i], proportion[i],
                                          select[i] - 1, shift[i] - 1};

        switch (event[i]) {
        case EXIT_EVENT:
        case ENTER_EVENT:
        case INTERNAL_TRANSFER_EVENT:
            j = (node[i] - 1) / chunk_size;
            if (j >= Nthread)
                j = Nthread - 1;
            kv_push(SimInf_scheduled_event, out[j].events, e);
            break;
        case EXTERNAL_TRANSFER_EVENT:
            kv_push(SimInf_scheduled_event, out[0].events, e);
            break;
        default:
            return SIMINF_UNDEFINED_EVENT;
        }
    }

    return 0;
}

/**
 * Create and initialize data to process scheduled events. The
 * generated data structure must be freed by the user.
 *
 * @param out the resulting data structure.
 * @param args structure with data for the solver.
 * @param rng random number generator
 * @return 0 or an error code
 */
int SimInf_scheduled_events_create(
    SimInf_scheduled_events **out, SimInf_solver_args *args, gsl_rng *rng)
{
    int error = SIMINF_ERR_ALLOC_MEMORY_BUFFER, i;
    SimInf_scheduled_events *events = NULL;

    events = calloc(args->Nthread, sizeof(SimInf_scheduled_events));
    if (!events)
        goto on_error;

    for (i = 0; i < args->Nthread; i++) {
        /*** Constants ***/
        events[i].Nthread = args->Nthread;

        /* Matrices to process events */
        events[i].irE = args->irE;
        events[i].jcE = args->jcE;
        events[i].N = args->N;

        /* Scheduled events */
	kv_init(events[i].events);

        events[i].individuals = calloc(args->Nc, sizeof(int));
        if (!events[i].individuals)
            goto on_error;

        events[i].u_tmp = calloc(args->Nc, sizeof(int));
        if (!events[i].u_tmp)
            goto on_error;

        /* Random number generator */
        events[i].rng = gsl_rng_alloc(gsl_rng_mt19937);
        if (!events[i].rng)
            goto on_error;
        gsl_rng_set(events[i].rng, gsl_rng_uniform_int(rng, gsl_rng_max(rng)));
    }

    /* Split scheduled events into E1 and E2 events. */
    error = SimInf_split_events(
        events, args->len, args->event, args->time, args->node,
        args->dest, args->n, args->proportion, args->select,
        args->shift, args->Nn, args->Nthread);
    if (error)
        goto on_error;

    *out = events;
    return 0;

on_error:
    SimInf_scheduled_events_free(events);
    return error;
}

/**
 * Free allocated memory to process events
 *
 * @param events SimInf_scheduled_events to free
 * @param Nthread number of threads that was used during simulation.
 */
void SimInf_scheduled_events_free(
    SimInf_scheduled_events *events)
{
    if (events) {
        int i;

        for (i = 0; i < events->Nthread; i++) {
            SimInf_scheduled_events *e = &events[i];

            if (e) {
                kv_destroy(e->events);
                free(e->individuals);
                e->individuals = NULL;
                free(e->u_tmp);
                e->u_tmp = NULL;
                gsl_rng_free(e->rng);
                e->rng = NULL;
            }
        }

        free(events);
    }
}

/**
 * Print event information to facilitate debugging.
 *
 * @param e The SimInf_scheduled_events object to print.
 * @param irE Select matrix for events. irE[k] is the row of E[k].
 * @param jcE Select matrix for events. Index to data of first
 *        non-zero element in row k.
 * @param Nc Number of compartments in each node.
 * @param u The state vector with number of individuals in each
 *        compartment at each node. The current state in each node is
 *        offset by node * Nc.
 * @param node The node in u.
 */
static void SimInf_print_event(
    const SimInf_scheduled_event *e, const int *irE, const int *jcE,
    const int Nc, const int *u, const int node, const int dest)
{
    #pragma omp critical
    {
        int i;

        if (irE && jcE && u) {
            int n = e->n, Nindividuals = 0, Nkinds = 0;

            /* 1) Count number of states with individuals */
            /* 2) Count total number of individuals       */
            for (i = jcE[e->select]; i < jcE[e->select + 1]; i++) {
                const int nk = u[node * Nc + irE[i]];
                if (nk > 0)
                    Nkinds++;
                Nindividuals += nk;
            }

            /* Number of states */
            if ((jcE[e->select + 1] - jcE[e->select]) <= 0)
                Rprintf("No states to sample from.\n");

            /* If n == 0, use the proportion of Nindividuals, else use
             * n as the number of individuals to sample */
            if (n == 0)
                n = round(e->proportion * Nindividuals);

            if (n > Nindividuals)
                Rprintf("Cannot sample %i for event from %i individuals.\n",
                        n, Nindividuals);

            if (n < 0)
                Rprintf("Cannot sample %i individuals for event.\n", n);

            Rprintf("\n");
        }

        if (u && (node >= 0)) {
            Rprintf("Current state in node\n");
            Rprintf("---------------------\n");

            Rprintf("{");
            for (i = 0; i < Nc; i++) {
                Rprintf("%i", u[node * Nc + i]);
                if (i < (Nc - 1))
                    Rprintf(", ");
            }
            Rprintf("}\n\n");
        }

        if (u && (dest >= 0)) {
            Rprintf("Current state in dest\n");
            Rprintf("---------------------\n");

            Rprintf("{");
            for (i = 0; i < Nc; i++) {
                Rprintf("%i", u[dest * Nc + i]);
                if (i < (Nc - 1))
                    Rprintf(", ");
            }
            Rprintf("}\n\n");
        }

        Rprintf("Scheduled event\n");
        Rprintf("---------------\n");

        switch (e->event) {
        case EXIT_EVENT:
            Rprintf("event: %i (exit event)\n", e->event);
            break;
        case ENTER_EVENT:
            Rprintf("event: %i (enter event)\n", e->event);
            break;
        case INTERNAL_TRANSFER_EVENT:
            Rprintf("event: %i (internal transfer event)\n", e->event);
            break;
        case EXTERNAL_TRANSFER_EVENT:
            Rprintf("event: %i (external transfer event)\n", e->event);
            break;
        default:
            Rprintf("event: %i (undefined event)\n", e->event);
            break;
        }

        Rprintf("time: %i\n", e->time);
        Rprintf("node: %i\n", e->node + 1); /* One based in events data */
        Rprintf("dest: %i\n", e->dest + 1); /* One based in events data */
        Rprintf("n: %i\n", e->n);
        Rprintf("proportion: %g\n", e->proportion);
        Rprintf("select: %i\n", e->select + 1); /* One based in events data */
        Rprintf("shift: %i\n\n", e->shift + 1); /* One based in events data */
    }
}

/**
 * Process all scheduled E1 and E2 events where time is less or equal
 * to the global time in the simulation.
 *
 * @param model The compartment model with information for each node
 * and the global time.
 * @param events Data with events to process.
 * @param process_E2 Process only E1 events (process_E2 = 0), else
 * process both E1 and E2 events.
 * @return 0 if Ok, else error code.
 */
void SimInf_process_events(
    SimInf_compartment_model *model,
    SimInf_scheduled_events *events,
    int process_E2)
{
    SimInf_compartment_model m = *&model[0];
    SimInf_scheduled_events e = *&events[0];

    /* Process events */
    while (e.events_index < kv_size(e.events) && !m.error) {
        int i;
        const SimInf_scheduled_event ee = kv_A(e.events, e.events_index);

        if (ee.time > m.tt)
            goto done;

        if (ee.node < 0 || ee.node >= m.Ntot) {
            SimInf_print_event(&ee, NULL, NULL, 0, NULL, -1, -1);
            m.error = SIMINF_ERR_NODE_OUT_OF_BOUNDS;
            goto done;
        }

        switch (ee.event) {
        case EXIT_EVENT:
            m.error = SimInf_sample_select(
                e.irE, e.jcE, m.Nc, m.u, ee.node - m.Ni, ee.select,
                ee.n, ee.proportion, e.individuals, e.u_tmp, e.rng);

            if (m.error) {
                SimInf_print_event(&ee, e.irE, e.jcE, m.Nc,
                                   m.u, ee.node - m.Ni, -1);
                goto done;
            }

            for (i = e.jcE[ee.select]; i < e.jcE[ee.select + 1]; i++) {
                const int jj = e.irE[i];
                const int kn = (ee.node - m.Ni) * m.Nc + jj;

                /* Remove individuals from node */
                m.u[kn] -= e.individuals[jj];
                if (m.u[kn] < 0) {
                    SimInf_print_event(&ee, NULL, NULL, m.Nc,
                                       m.u, ee.node - m.Ni, -1);
                    m.error = SIMINF_ERR_NEGATIVE_STATE;
                    goto done;
                }
            }
            break;

        case ENTER_EVENT:
            /* All individuals enter first non-zero compartment,
             * i.e. a non-zero entry in element in the select
             * column. */
            if (e.jcE[ee.select] < e.jcE[ee.select + 1]) {
                m.u[(ee.node - m.Ni) * m.Nc + e.irE[e.jcE[ee.select]]] += ee.n;
                if (m.u[(ee.node - m.Ni) * m.Nc + e.irE[e.jcE[ee.select]]] < 0) {
                    SimInf_print_event(&ee, NULL, NULL, m.Nc,
                                       m.u, ee.node - m.Ni, -1);
                    m.error = SIMINF_ERR_NEGATIVE_STATE;
                }
            }
            break;

        case INTERNAL_TRANSFER_EVENT:
            if (!e.N) {
                /* Not possible to shift when N is not defined. */
                SimInf_print_event(&ee, NULL, NULL, 0, NULL, -1, -1);
                m.error = SIMINF_ERR_EVENTS_N;
                goto done;
            }

            if (ee.shift < 0) {
                /* Invalid shift parameter. */
                SimInf_print_event(&ee, NULL, NULL, 0, NULL, -1, -1);
                m.error = SIMINF_ERR_EVENT_SHIFT;
                goto done;
            }

            m.error = SimInf_sample_select(
                e.irE, e.jcE, m.Nc, m.u, ee.node - m.Ni, ee.select,
                ee.n, ee.proportion, e.individuals, e.u_tmp, e.rng);

            if (m.error) {
                SimInf_print_event(&ee, e.irE, e.jcE, m.Nc,
                                   m.u, ee.node - m.Ni, -1);
                goto done;
            }

            for (i = e.jcE[ee.select]; i < e.jcE[ee.select + 1]; i++) {
                const int jj = e.irE[i];
                const int kn = (ee.node - m.Ni) * m.Nc + jj;
                const int ll = e.N[ee.shift * m.Nc + jj];

                /* Check that the index to the new compartment is not
                 * out of bounds. */
                if (jj + ll < 0 || jj + ll >= m.Nc) {
                    SimInf_print_event(&ee, NULL, NULL, 0, NULL, -1, -1);
                    m.error = SIMINF_ERR_SHIFT_OUT_OF_BOUNDS;
                    goto done;
                }

                /* Add individuals to new compartments in node */
                m.u[kn + ll] += e.individuals[jj];
                if (m.u[kn + ll] < 0) {
                    SimInf_print_event(&ee, NULL, NULL, m.Nc,
                                       m.u, ee.node - m.Ni, -1);
                    m.error = SIMINF_ERR_NEGATIVE_STATE;
                    goto done;
                }

                /* Remove individuals from previous compartments in
                 * node */
                m.u[kn] -= e.individuals[jj];
                if (m.u[kn] < 0) {
                    SimInf_print_event(&ee, NULL, NULL, m.Nc,
                                       m.u, ee.node - m.Ni, -1);
                    m.error = SIMINF_ERR_NEGATIVE_STATE;
                    goto done;
                }
            }
            break;

        case EXTERNAL_TRANSFER_EVENT:
            /* Check if we are done because we only want to process E1
             * events. */
            if (!process_E2)
                goto done;

            if (ee.dest < 0 || ee.dest >= m.Ntot) {
                SimInf_print_event(&ee, NULL, NULL, 0, NULL, -1, -1);
                m.error = SIMINF_ERR_DEST_OUT_OF_BOUNDS;
                goto done;
            }

            m.error = SimInf_sample_select(
                e.irE, e.jcE, m.Nc, m.u, ee.node, ee.select, ee.n,
                ee.proportion, e.individuals, e.u_tmp, e.rng);

            if (m.error) {
                SimInf_print_event(&ee, e.irE, e.jcE, m.Nc,
                                   m.u, ee.node, ee.dest);
                goto done;
            }

            for (i = e.jcE[ee.select]; i < e.jcE[ee.select + 1]; i++) {
                const int jj = e.irE[i];
                const int kd = ee.dest * m.Nc + jj;
                const int kn = ee.node * m.Nc + jj;

                if (ee.shift < 0) {
                    /* Add individuals to dest without shifting
                     * compartments */
                    m.u[kd] += e.individuals[jj];
                    if (m.u[kd] < 0) {
                        SimInf_print_event(&ee, NULL, NULL, m.Nc,
                                           m.u, ee.node, ee.dest);
                        m.error = SIMINF_ERR_NEGATIVE_STATE;
                        goto done;
                    }
                } else if (!e.N) {
                    /* Not possible to shift when N is not defined. */
                    SimInf_print_event(&ee, NULL, NULL, 0, NULL, -1, -1);
                    m.error = SIMINF_ERR_EVENTS_N;
                    goto done;
                } else {
                    /* Process a movement event that also involves a
                     * shift between compartments. */
                    const int ll = e.N[ee.shift * m.Nc + jj];

                    /* Check that the index to the new compartment is
                     * not out of bounds. */
                    if (jj + ll < 0 || jj + ll >= m.Nc) {
                        SimInf_print_event(&ee, NULL, NULL, 0, NULL, -1, -1);
                        m.error = SIMINF_ERR_SHIFT_OUT_OF_BOUNDS;
                        goto done;
                    }

                    /* Add individuals to dest */
                    m.u[kd + ll] += e.individuals[jj];
                    if (m.u[kd + ll] < 0) {
                        SimInf_print_event(&ee, NULL, NULL, m.Nc,
                                           m.u, ee.node, ee.dest);
                        m.error = SIMINF_ERR_NEGATIVE_STATE;
                        goto done;
                    }
                }

                /* Remove individuals from node */
                m.u[kn] -= e.individuals[jj];
                if (m.u[kn] < 0) {
                    SimInf_print_event(&ee, NULL, NULL, m.Nc,
                                       m.u, ee.node, ee.dest);
                    m.error = SIMINF_ERR_NEGATIVE_STATE;
                    goto done;
                }
            }

            /* Indicate dest for update */
            m.update_node[ee.dest] = 1;
            break;

        default:
            SimInf_print_event(&ee, NULL, NULL, 0, NULL, -1, -1);
            m.error = SIMINF_UNDEFINED_EVENT;
            break;
        }

        /* Indicate node for update */
        m.update_node[ee.node - m.Ni] = 1;

        e.events_index++;
    }

done:
    *&events[0] = e;
    *&model[0] = m;
}

/**
 * Handle the case where the solution is stored in a sparse matrix
 *
 * Store solution if tt has passed the next time in tspan. Report
 * solution up to, but not including tt.
 *
 * @param SimInf_compartment_model *model data to store.
 */
void SimInf_store_solution_sparse(SimInf_compartment_model *model)
{
    while (!model[0].U && model[0].U_it < model[0].tlen &&
           model[0].tt > model[0].tspan[model[0].U_it]) {
        int j;

        /* Copy compartment state to U_sparse */
        for (j = model[0].jcU[model[0].U_it];
             j < model[0].jcU[model[0].U_it + 1]; j++)
            model[0].prU[j] = model[0].u[model[0].irU[j]];
        model[0].U_it++;
    }

    while (!model[0].V && model[0].V_it < model[0].tlen &&
           model[0].tt > model[0].tspan[model[0].V_it]) {
        int j;

        /* Copy continuous state to V_sparse */
        for (j = model[0].jcV[model[0].V_it];
             j < model[0].jcV[model[0].V_it + 1]; j++)
            model[0].prV[j] = model[0].v_new[model[0].irV[j]];
        model[0].V_it++;
    }
}

/**
 * Free allocated memory for an epidemiological compartment
 * model.
 *
 * @param model the data structure to free.
 * @param Nthread number of threads that was used during simulation.
 */
void SimInf_compartment_model_free(SimInf_compartment_model *model)
{
    if (model) {
        int i;

        for (i = 0; i < model->Nthread; i++) {
            SimInf_compartment_model *m = &model[i];

            if (m) {
                free(m->t_rate);
                m->t_rate = NULL;
                free(m->sum_t_rate);
                m->sum_t_rate = NULL;
                free(m->t_time);
                m->t_time = NULL;
            }
        }

        free(model[0].u);
        model[0].u = NULL;
        free(model[0].v);
        model[0].v = NULL;
        free(model[0].v_new);
        model[0].v_new = NULL;
        free(model[0].update_node);
        model[0].update_node = NULL;
        free(model);
    }
}

/**
 * Create and initialize data for an epidemiological compartment
 * model. The generated model must be freed by the user.
 *
 * @param out the resulting data structure.
 * @param args structure with data for the solver.
 * @return 0 or SIMINF_ERR_ALLOC_MEMORY_BUFFER
 */
int SimInf_compartment_model_create(
    SimInf_compartment_model **out, SimInf_solver_args *args)
{
    int i;
    SimInf_compartment_model *model = NULL;

    /* Allocate memory for the compartment model. */
    model = calloc(args->Nthread, sizeof(SimInf_compartment_model));
    if (!model)
        goto on_error;

    /* Allocate memory to keep track of the continuous state in each
     * node. */
    model[0].v = malloc(args->Nn * args->Nd * sizeof(double));
    if (!model[0].v)
        goto on_error;
    model[0].v_new = malloc(args->Nn * args->Nd * sizeof(double));
    if (!model[0].v_new)
        goto on_error;

    /* Set continuous state to the initial state in each node. */
    memcpy(model[0].v, args->v0, args->Nn * args->Nd * sizeof(double));

    /* Setup vector to keep track of nodes that must be updated due to
     * scheduled events */
    model[0].update_node = calloc(args->Nn, sizeof(int));
    if (!model[0].update_node)
        goto on_error;

    /* Allocate memory for compartment state and set compartment state
     * to the initial state. */
    model[0].u = malloc(args->Nn * args->Nc * sizeof(int));
    if (!model[0].u)
        goto on_error;
    memcpy(model[0].u, args->u0, args->Nn * args->Nc * sizeof(int));

    for (i = 0; i < args->Nthread; i++) {
        /* Constants */
        model[i].Nthread = args->Nthread;
        model[i].Ntot = args->Nn;
        model[i].Ni = i * (args->Nn / args->Nthread);
        model[i].Nn = args->Nn / args->Nthread;
        if (i == (args->Nthread - 1))
            model[i].Nn += (args->Nn % args->Nthread);
        model[i].Nt = args->Nt;
        model[i].Nc = args->Nc;
        model[i].Nd = args->Nd;
        model[i].Nld = args->Nld;

        /* Sparse matrices */
        model[i].irG = args->irG;
        model[i].jcG = args->jcG;
        model[i].irS = args->irS;
        model[i].jcS = args->jcS;
        model[i].prS = args->prS;

        /* Callbacks */
        model[i].tr_fun = args->tr_fun;
        model[i].pts_fun = args->pts_fun;

        /* Keep track of time */
        model[i].tt = args->tspan[0];
        model[i].next_unit_of_time = floor(model[i].tt) + 1.0;
        model[i].tspan = args->tspan;
        model[i].tlen = args->tlen;
        model[i].U_it = 1;
        model[i].V_it = 1;

        /* Data vectors */
        if (args->U) {
            model[i].U = args->U;
        } else if (i == 0) {
            model[i].irU = args->irU;
            model[i].jcU = args->jcU;
            model[i].prU = args->prU;
        }

        if (args->V) {
            model[i].V = args->V;
        } else if (i == 0) {
            model[i].irV = args->irV;
            model[i].jcV = args->jcV;
            model[i].prV = args->prV;
        }

        if (i > 0) {
            model[i].u = &model[0].u[model[i].Ni * args->Nc];
            model[i].v = &model[0].v[model[i].Ni * args->Nd];
            model[i].v_new = &model[0].v_new[model[i].Ni * args->Nd];
            model[i].update_node = &model[0].update_node[model[i].Ni];
        }

        model[i].ldata = &(args->ldata[model[i].Ni * model[i].Nld]);
        model[i].gdata = args->gdata;

        /* Create transition rate matrix (Nt X Nn) and total rate
         * vector. In t_rate we store all propensities for state
         * transitions, and in sum_t_rate the sum of propensities in
         * every node. */
        model[i].t_rate = malloc(args->Nt * model[i].Nn * sizeof(double));
        if (!model[i].t_rate)
            goto on_error;
        model[i].sum_t_rate = malloc(model[i].Nn * sizeof(double));
        if (!model[i].sum_t_rate)
            goto on_error;
        model[i].t_time = malloc(model[i].Nn * sizeof(double));
        if (!model[i].t_time)
            goto on_error;
    }

    *out = model;
    return 0;

on_error:
    SimInf_compartment_model_free(model);
    return SIMINF_ERR_ALLOC_MEMORY_BUFFER;
}
