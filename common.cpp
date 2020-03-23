/*
 * Filename: common.cpp
 * Contributors: Brandon Barker, Akash Kumar, Andrew Valenci
 */
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <float.h>
#include <string.h>
#include <math.h>
#include <vector>
#include <time.h>
#include <sys/time.h>
#include "common.h"

double size;

//
//  tuned constants
//
#define density 0.0005
#define mass    0.01
#define cutoff  0.01         // Range of interaction forces.
#define min_r   (cutoff/100)
#define dt      0.0005

//
//  timer
//
double read_timer( )
{
    static bool initialized = false;
    static struct timeval start;
    struct timeval end;
    if( !initialized )
    {
        gettimeofday( &start, NULL );
        initialized = true;
    }
    gettimeofday( &end, NULL );
    return (end.tv_sec - start.tv_sec) + 1.0e-6 * (end.tv_usec - start.tv_usec);
}

//
//  keep density constant
//
void set_size( int n )
{
    size = sqrt( density * n );
}

//
//  Initialize the particle positions and velocities
//
void init_particles( int n, particle_t *p )
{
    srand48( time( NULL ) );

    int sx = (int)ceil(sqrt((double)n));
    int sy = (n+sx-1)/sx;

    int *shuffle = (int*)malloc( n * sizeof(int) );
    for( int i = 0; i < n; i++ )
        shuffle[i] = i;

    for( int i = 0; i < n; i++ )
    {
        //
        //  make sure particles are not spatially sorted
        //
        int j = lrand48()%(n-i);
        int k = shuffle[j];
        shuffle[j] = shuffle[n-i-1];

        //
        //  distribute particles evenly to ensure proper spacing
        //
        p[i].x = size*(1.+(k%sx))/(1+sx);
        p[i].y = size*(1.+(k/sx))/(1+sy);

        //
        //  assign random velocities within a bound
        //
        p[i].vx = drand48()*2-1;
        p[i].vy = drand48()*2-1;
    }
    free( shuffle );
}

//
//  interact two particles
//
void apply_force( particle_t &particle, particle_t &neighbor , double *dmin, double *davg, int *navg)
{

    double dx = neighbor.x - particle.x;
    double dy = neighbor.y - particle.y;
    double r2 = dx * dx + dy * dy;
    if( r2 > cutoff*cutoff )
        return;
	if (r2 != 0)
        {
	   if (r2/(cutoff*cutoff) < *dmin * (*dmin))
	      *dmin = sqrt(r2)/cutoff;
           (*davg) += sqrt(r2)/cutoff;
           (*navg) ++;
        }

    r2 = fmax( r2, min_r*min_r );
    double r = sqrt( r2 );



    //
    //  very simple short-range repulsive force
    //
    double coef = ( 1 - cutoff / r ) / r2 / mass;
    particle.ax += coef * dx;
    particle.ay += coef * dy;
}

//
//  integrate the ODE
//
void move( particle_t &p )
{
    //
    //  slightly simplified Velocity Verlet integration
    //  conserves energy better than explicit Euler method
    //
    p.vx += p.ax * dt;
    p.vy += p.ay * dt;
    p.x  += p.vx * dt;
    p.y  += p.vy * dt;

    //
    //  bounce from walls
    //
    while( p.x < 0 || p.x > size )
    {
        p.x  = p.x < 0 ? -p.x : 2*size-p.x;
        p.vx = -p.vx;
    }
    while( p.y < 0 || p.y > size )
    {
        p.y  = p.y < 0 ? -p.y : 2*size-p.y;
        p.vy = -p.vy;
    }
}

//
//  I/O routines
//
void save( FILE *f, int n, particle_t *p )
{
    static bool first = true;
    if( first )
    {
        fprintf( f, "%d %g\n", n, size );
        first = false;
    }
    for( int i = 0; i < n; i++ )
        fprintf( f, "%g %g\n", p[i].x, p[i].y );
}

//
//  command line option processing
//
int find_option( int argc, char **argv, const char *option )
{
    for( int i = 1; i < argc; i++ )
        if( strcmp( argv[i], option ) == 0 )
            return i;
    return -1;
}

int read_int( int argc, char **argv, const char *option, int default_value )
{
    int iplace = find_option( argc, argv, option );
    if( iplace >= 0 && iplace < argc-1 )
        return atoi( argv[iplace+1] );
    return default_value;
}

char *read_string( int argc, char **argv, const char *option, char *default_value )
{
    int iplace = find_option( argc, argv, option );
    if( iplace >= 0 && iplace < argc-1 )
        return argv[iplace+1];
    return default_value;
}

//
// Student additions
//
double bin_size;
int bin_count;   // Number of bins per row/col in spatial hash.
// Represent a 2D particle grid as a 1D vector of particle bins.
std::vector<bin_t> particle_bins;

// A function that is used to create a spatial hash
// given particles of type particle_t.
void make_spatial_hash(int n, particle_t* particles)
{
    bin_size = cutoff * 2;
    bin_count = int(size / bin_size) + 1;

    particle_bins.resize(bin_count * bin_count);

    // Put each particle into its appropriate bin
    // in the 1D particle_bins vector.
    for (int i = 0; i < n; i++)
    {
        bin_particle(particles[i]);
        // int x = int(particles[i].x / bin_size);
        // int y = int(particles[i].y / bin_size);
        // particle_bins[x*bin_count + y].push_back(particles[i]);
    }
}

// A function that is used to add a particle to the spatial hash.
void bin_particle(particle_t& particle)
{
    int x = particle.x / size;
    int y = particle.y / size;
    particle_bins[x*bin_count + y].push_back(particle);
}

// A function that is used to compute all particle forces in a spatial hash bin.
void compute_bin_forces(int grid_row, int grid_col,
                        double& dmin, double& davg, int& navg)
{
    bin_t& bin = particle_bins[grid_row*bin_count + grid_col];
    // Reset acceleration of all particles in the bin.
    for (int k = 0; k < bin.size(); k++)
        bin[k].ax = bin[k].ay = 0;
    // Each non-edge particle bin is surround by 8 bins.
    // Apply forces between particles in the current bin
    // and particles in the surrounding 8 bins as well as
    // between particles in the same bin.
    for (int dx = -1; dx <= 1; dx++)
    {
        for (int dy = -1; dy <= 1; dy++)
        {
            if (grid_row + dx >= 0 && grid_row + dx < bin_count &&
                grid_col + dy >= 0 && grid_col + dy < bin_count)
            {   // bin2 represents one of the surrounding 8 bins.
                bin_t& bin2 = particle_bins[(grid_row+dx)*bin_count + grid_col + dy];
                // Apply particle forces between the two bins.
                for (int i = 0; i < bin.size(); i++)
                    for (int j = 0; j < bin2.size(); j++)
                        apply_force(bin[i], bin2[j], &dmin, &davg, &navg);
            }
        }
    }
}
