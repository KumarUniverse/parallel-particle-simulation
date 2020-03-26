/*
 * Filename: openmp.cpp
 * Contributors: Brandon Barker, Akash Kumar, Andrew Valenci
 */
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <iostream>
#include <vector>
#include "common.h"
#include "omp.h"

//
//  benchmarking program
//
int main( int argc, char **argv )
{
    int navg,nabsavg=0,numthreads;
    double dmin, absmin=1.0,davg,absavg=0.0;

    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        printf( "-s <filename> to specify a summary file name\n" );
        printf( "-no turns off all correctness checks and particle output\n");
        return 0;
    }

    int n = read_int(argc, argv, "-n", 1000);
    char *savename = read_string(argc, argv, "-o", NULL);
    char *sumname = read_string(argc, argv, "-s", NULL);

    FILE *fsave = savename ? fopen(savename, "w") : NULL;
    FILE *fsum = sumname ? fopen(sumname, "a") : NULL;

    particle_t*particles = (particle_t*) malloc(n * sizeof(particle_t));
    set_size(n);
    init_particles(n, particles);

    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer();

    {
        for (int step = 0; step < NSTEPS; step++) {
            navg = 0;
            davg = 0.0;
            dmin = 1.0;

            double density = .0005;
            double size = sqrt(density * n);
            int grid_num = 16;
            std::vector<particle_t*> b[grid_num][grid_num];

            for (int i = 0; i < grid_num; i++) {
                for (int j = 0; j < grid_num; j++) {
                    b[i][j] = std::vector<particle_t*>();
                }
            }
            double cell_size = size / grid_num;

            // Fill matrix with particles
            for (int i = 0; i < n; i++) {
                double current_x = particles[i].x;
                int bin_x = floor(current_x / cell_size);

                double current_y = particles[i].y;
                int bin_y = floor(current_y / cell_size);

                b[bin_y][bin_x].push_back(&particles[i]);
            }

            // Compute forces


            #pragma omp parallel
            {
                #pragma omp for collapse(2) reduction(+:navg) reduction(+:davg) firstprivate(dmin) schedule(dynamic)

                for (int grid_row = 0; grid_row < grid_num; grid_row++) {
                    for (int grid_col = 0; grid_col < grid_num; grid_col++) {

                        numthreads = omp_get_num_threads();
                        std::vector<std::vector<particle_t*> > changed;
                        for (int i = grid_row - 1; i <= grid_row + 1; i++) {
                            for (int j = grid_col - 1; j <= grid_col + 1; j++) {
                                if (i >= 0 && i < grid_num && j >= 0 && j < grid_num) {
                                    changed.push_back(b[i][j]);
                                }
                            }
                        }

                        std::vector<particle_t*> particle_b = b[grid_row][grid_col];

                        for (int p = 0; p < particle_b.size(); p++) {
                            particle_b[p]->ax = particle_b[p]->ay = 0;
                            for (int neighbor_index = 0; neighbor_index < changed.size(); neighbor_index++) {
                                std::vector<particle_t*> current_neighbor = changed[neighbor_index];
                                for (int np = 0; np < current_neighbor.size(); np++) {
                                    apply_force(*particle_b[p], *current_neighbor[np], &dmin, &davg, &navg);
                                }
                            }
                        }
                    }
                }



                //
                //  move particles
                //
                #pragma omp for
                for (int i = 0; i < n; i++)
                    move(particles[i]);

                if (find_option(argc, argv, "-no") == -1) {
                    //
                    //  compute statistical data
                    //
                    #pragma omp master
                    if (navg) {
                        absavg += davg / navg;
                        nabsavg++;
                    }

                    #pragma omp critical
                    if (dmin < absmin) absmin = dmin;

                    //
                    //  save if necessary
                    //
                    #pragma omp master
                    if (fsave && (step % SAVEFREQ) == 0)
                        save(fsave, n, particles);

                }
            }
        }
    }
    simulation_time = read_timer() - simulation_time;

    printf("n = %d,threads = %d, simulation time = %g seconds", n, numthreads, simulation_time);

    if (find_option(argc, argv, "-no") == -1) {
        if (nabsavg) absavg /= nabsavg;
        //
        //  -The minimum distance absmin between 2 particles during the run of the simulation
        //  -A Correct simulation will have particles stay at greater than 0.4 (of cutoff) with typical values between .7-.8
        //  -A simulation where particles don't interact correctly will be less than 0.4 (of cutoff) with typical values between .01-.05
        //
        //  -The average distance absavg is ~.95 when most particles are interacting correctly and ~.66 when no particles are interacting
        //
        printf(", absmin = %lf, absavg = %lf", absmin, absavg);
        if (absmin < 0.4) printf("\nThe minimum distance is below 0.4 meaning that some particle is not interacting");
        if (absavg < 0.8) printf("\nThe average distance is below 0.8 meaning that most particles are not interacting");
    }
    printf("\n");

    //
    // Printing summary data
    //
    if (fsum)
        fprintf(fsum, "%d %d %g\n", n, numthreads, simulation_time);

    //
    // Clearing space
    //
    if (fsum)
        fclose(fsum);

    free(particles);
    if (fsave)
        fclose(fsave);

    return 0;
}
