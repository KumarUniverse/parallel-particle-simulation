/*
 * Filename: mpi.cpp
 * Contributors: Brandon Barker, Akash Kumar, Andrew Valenci
 */
#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <vector>
#include "common.h"

//
//  benchmarking program
//
int main( int argc, char **argv )
{
    int navg, nabsavg=0;
    double dmin, absmin=1.0,davg,absavg=0.0;
    double rdavg,rdmin;
    int rnavg;

    //
    //  process command line parameters
    //
    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        printf( "-s <filename> to specify a summary file name\n" );
        printf( "-no turns off all correctness checks and particle output\n");
        return 0;
    }

    int n = read_int( argc, argv, "-n", 1000 );
    char *savename = read_string( argc, argv, "-o", NULL );
    char *sumname = read_string( argc, argv, "-s", NULL );

    //
    //  set up MPI
    //
    int n_proc, // Size of the group associated with MPI_COMM_WORLD communicator.
        rank;   // Rank of the process in MPI_COMM_WORLD.
    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &n_proc );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );

    //
    //  allocate generic resources
    //
    FILE *fsave = savename && rank == 0 ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname && rank == 0 ? fopen( sumname, "a" ) : NULL;


    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );

    MPI_Datatype PARTICLE;
    MPI_Type_contiguous( 6, MPI_DOUBLE, &PARTICLE );
    MPI_Type_commit( &PARTICLE );

    //
    //  initialize and distribute the particles (that's fine to leave it unoptimized)
    //
    set_size( n );
    if( rank == 0 )
        init_particles( n, particles );
    // MPI_Scatterv( particles, partition_sizes, partition_offsets, PARTICLE, local, nlocal, PARTICLE, 0, MPI_COMM_WORLD );
    MPI_Bcast(particles, n, PARTICLE, 0, MPI_COMM_WORLD);

    // Create spatial hash.
    make_spatial_hash( n, particles );


    // Divide the particle bins evenly between the MPI processes.
    // Each MPI process only performs computations on all bins
    // starting from the start_bin (inclusive) and ending at
    // the end_bin (exclusive).

    int bins_per_proc = bin_count / n_proc;
    int start_bin = bins_per_proc * rank;
    int end_bin = bins_per_proc * (rank + 1);

    // Special case for last MPI process.
    if (rank == n_proc - 1)
        end_bin = bin_count;

    //
    //  Simulate a number of time steps
    //
    double simulation_time = read_timer();
    for( int step = 0; step < NSTEPS; step++ )
    {
        navg = 0;
        dmin = 1.0;
        davg = 0.0;

        // Compute forces
        for (int grid_row = start_bin; grid_row < end_bin; grid_row++)
        {
            for (int grid_col = 0; grid_col < bin_count; grid_col++)
            {
                compute_bin_forces(grid_row, grid_col, dmin, davg, navg);
            }
        }

        if( find_option( argc, argv, "-no" ) == -1 )
        {
            MPI_Reduce(&davg,&rdavg,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
            MPI_Reduce(&navg,&rnavg,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
            MPI_Reduce(&dmin,&rdmin,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);

            if (rank == 0)
            {
                //
                // Computing statistical data
                //
                if (rnavg)
                {
                   absavg +=  rdavg/rnavg;
                   nabsavg++;
                }
                if (rdmin < absmin) absmin = rdmin;
            }
        }

        //
        // Move particles
        //
        bin_t remotely_changed;  // To keep track of particles which change bins but are still within the local process.
        bin_t locally_changed; // To keep track of particles which change bins which are out of the local process's control.

        for (int grid_row = start_bin; grid_row < end_bin; ++grid_row)
        {
            for (int grid_col = 0; grid_col < bin_count; ++grid_col)
            {
                bin_t& bin = particle_bins[grid_row*bin_count + grid_col];
                int tail = bin.size(), i = 0;
                while (i < tail)
                {
                    move(bin[i]);                     // Move particle.
                    int x = int(bin[i].x / bin_size); // Check the position.
                    int y = int(bin[i].y / bin_size);
                    if (x <= start_bin && x < end_bin)
                    {
                        if (x == grid_row && y == grid_col) // Still inside original bin.
                            i++;
                        else
                        {
                            locally_changed.push_back(bin[i]); // Particle moves within current process's control.
                            bin[i] = bin[--tail]; // Remove the particle from the current bin.
                        }
                    }
                    else
                    {
                        remotely_changed.push_back(bin[i]); // Particle moves out of current process's control.
                        bin[i] = bin[--tail]; // Remove the particle from the current bin.
                    }
                }
                bin.resize(i); // Remove duplicates and shrink the bin.
            }
        }

        // Put all the locally changed particles into their correct bins.
        for (int i = 0; i < locally_changed.size(); i++)
        {
            bin_particle(locally_changed[i]);
        }

        bin_t incoming_move;
        int sendcount = remotely_changed.size();
        int recvbuf[n_proc];

        // Have each process send their size of their remotely changed bin to the root process.
        MPI_Gather(&sendcount, 1, MPI_INT, recvbuf, 1, MPI_INT, 0, MPI_COMM_WORLD);

        int displs[n_proc];
        int total_num = 0;

        if (rank == 0)
        {
            displs[0] = 0;
            for (int i = 1; i < n_proc; ++i) {
                displs[i] = displs[i-1] + recvbuf[i-1];
            }
            total_num = displs[n_proc-1] + recvbuf[n_proc-1];
            incoming_move.resize(total_num);
        }

        // Have each process send their remotely changed data to the root process.
        MPI_Gatherv(remotely_changed.data(), sendcount, PARTICLE,
                    incoming_move.data(), recvbuf, displs, PARTICLE,
                    0, MPI_COMM_WORLD);

        std::vector<bin_t> scatter_particles;
        scatter_particles.resize(n_proc);

        // Root process scatters the particles to their respective processes.
        if (rank == 0)
        {
            for (int i = 0; i < incoming_move.size(); i++)
            {
                int x = int(incoming_move[i].x / bin_size);

                assert(incoming_move[i].x >= 0 && incoming_move[i].y >= 0 &&
                       incoming_move[i].x <= size && incoming_move[i].y <= size);

                int who = min(x / bins_per_proc, n_proc-1);
                scatter_particles[who].push_back(incoming_move[i]);

                int row = x % bins_per_proc;
                if (row == 0 && who != 0)
                    scatter_particles[who - 1].push_back(incoming_move[i]);
                if (row == bins_per_proc-1 && who != n_proc-1)
                    scatter_particles[who + 1].push_back(incoming_move[i]);
            }
            for (int i = 0; i < n_proc; i++)
            {
                recvbuf[i] = scatter_particles[i].size();
            }
            displs[0] = 0;
            for (int i = 1; i < n_proc; i++)
            {
                displs[i] = displs[i-1] + recvbuf[i-1];
            }
        }

        sendcount = 0;
        // Scatter the size of their scatter bin to all processes.
        MPI_Scatter(recvbuf, 1, MPI_INT, &sendcount, 1, MPI_INT, 0, MPI_COMM_WORLD);

        bin_t outgoing_move;
        outgoing_move.resize(sendcount);

        bin_t scatter_particles_flatten;
        for (int i = 0; i < scatter_particles.size(); ++i)
        {
            scatter_particles_flatten.insert(scatter_particles_flatten.end(),
                                             scatter_particles[i].begin(), scatter_particles[i].end());
        }

        // Scatter the scatter particles to all processes.
        MPI_Scatterv(scatter_particles_flatten.data(), recvbuf, displs, PARTICLE,
            outgoing_move.data(), sendcount, PARTICLE, 0, MPI_COMM_WORLD);

        // Bin the outgoing particles.
        for (int i = 0; i < sendcount; i++)
        {
            particle_t &p = outgoing_move[i];
            assert(p.x >= 0 && p.y >= 0 && p.x <= size && p.y <= size);
            bin_particle(p);
        }
    }
    simulation_time = read_timer() - simulation_time;

    if (rank == 0)
    {
        printf( "n = %d, simulation time = %g seconds", n, simulation_time);

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
        if(fsum)
            fprintf(fsum,"%d %d %g\n",n,n_proc,simulation_time);
    }

    //
    //  release resources
    //
    if ( fsum )
        fclose( fsum );
    // free( partition_offsets );
    // free( partition_sizes );
    // free( local );
    free( particles );
    if( fsave )
        fclose( fsave );

    MPI_Finalize( );

    return 0;
}
