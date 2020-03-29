/*
 * Filename: openmp.cpp
 * Contributors: Brandon Barker, Akash Kumar, Andrew Valenci
 */
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
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

    int n = read_int( argc, argv, "-n", 1000 );
    char *savename = read_string( argc, argv, "-o", NULL );
    char *sumname = read_string( argc, argv, "-s", NULL );

    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname ? fopen ( sumname, "a" ) : NULL;

    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    //bin_t& changed; // To keep track of particles which change bins.
    std::vector<bin_t> thread_bins; // All the changed particle bins. One changed bin per thread.

    set_size( n );
    init_particles( n, particles );
    make_spatial_hash( n , particles );

    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );

    #pragma omp parallel private(dmin)
    {
        #pragma omp master // This code is only executed once by the master thread.
        {
            numthreads = omp_get_num_threads();
            thread_bins.resize(numthreads);
        }

        for( int step = 0; step < NSTEPS; step++ )
        {
            navg = 0;
            davg = 0.0;
            dmin = 1.0;
            //
            //  compute all forces
            //
            // Each thread performs partial reduction on addition operations on navg and davg.
            #pragma omp for reduction (+:navg) reduction(+:davg)
            for (int grid_row = 0; grid_row < bin_count; grid_row++)
            {
                for (int grid_col = 0; grid_col < bin_count; grid_col++)
                {
                    compute_bin_forces(grid_row, grid_col, dmin, davg, navg);
                }
            }

            //
            //  move particles
            //
            int id = omp_get_thread_num();
            bin_t& changed = thread_bins[id];
            changed.clear(); // No particles changed bins yet.
            #pragma omp for
            for (int grid_row = 0; grid_row < bin_count; grid_row++)
            {
                for (int grid_col = 0; grid_col < bin_count; grid_col++)
                {
                    bin_t& bin = particle_bins[grid_row*bin_count + grid_col];
                    int tail = bin.size(), i = 0;
                    while (i < tail)
                    {
                        move(bin[i]);                     // Move particle.
                        int x = int(bin[i].x / bin_size); // Check the position.
                        int y = int(bin[i].y / bin_size);
                        if (x == grid_row && y == grid_col) // Still inside original bin.
                            i++;
                        else
                        {
                            changed.push_back(bin[i]); // Store particles that have changed bins.
                            bin[i] = bin[--tail]; // Remove the particle from the current bin.
                        }
                    }
                    bin.resize(i); // Remove duplicates and shrink the bin.
                }
            }

            #pragma omp master // This code is only executed once by the master thread.
            {
                for (int i = 0; i < numthreads; i++)
                {
                    changed = thread_bins[i];
                    // Put the particles that have changed bins
                    // into their respective bins.
                    for (int j = 0; j < changed.size(); j++)
                    {
                        int x = int(changed[j].x / bin_size);
                        int y = int(changed[j].y / bin_size);
                        particle_bins[x*bin_count + y].push_back(changed[j]);
                    }
                }
            }

            if( find_option( argc, argv, "-no" ) == -1 )
            {
                //
                //  compute statistical data
                //
                #pragma omp master
                if (navg) {
                    absavg += davg/navg;
                    nabsavg++;
                }

                #pragma omp critical
                if (dmin < absmin) absmin = dmin;

                //
                //  save if necessary
                //
                #pragma omp master
                if( fsave && (step%SAVEFREQ) == 0 )
                    save( fsave, n, particles );
            }

            #pragma omp barrier // Wait for all the threads to finish.
        }
    }
    simulation_time = read_timer( ) - simulation_time;

    printf( "n = %d,threads = %d, simulation time = %g seconds", n,numthreads, simulation_time);

    if( find_option( argc, argv, "-no" ) == -1 )
    {
      if (nabsavg) absavg /= nabsavg;
    //
    //  -the minimum distance absmin between 2 particles during the run of the simulation
    //  -A Correct simulation will have particles stay at greater than 0.4 (of cutoff) with typical values between .7-.8
    //  -A simulation were particles don't interact correctly will be less than 0.4 (of cutoff) with typical values between .01-.05
    //
    //  -The average distance absavg is ~.95 when most particles are interacting correctly and ~.66 when no particles are interacting
    //
    printf( ", absmin = %lf, absavg = %lf", absmin, absavg);
    if (absmin < 0.4) printf ("\nThe minimum distance is below 0.4 meaning that some particle is not interacting");
    if (absavg < 0.8) printf ("\nThe average distance is below 0.8 meaning that most particles are not interacting");
    }
    printf("\n");

    //
    // Printing summary data
    //
    if( fsum)
        fprintf(fsum,"%d %d %g\n",n,numthreads,simulation_time);

    //
    // Clearing space
    //
    if( fsum )
        fclose( fsum );

    free( particles );
    if( fsave )
        fclose( fsave );

    return 0;
}
