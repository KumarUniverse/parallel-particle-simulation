/*
 * Filename: common.h
 * Contributors: Brandon Barker, Akash Kumar, Andrew Valenci
 */
#ifndef __CS267_COMMON_H__
#define __CS267_COMMON_H__

#include <vector>

inline int min( int a, int b ) { return a < b ? a : b; }
inline int max( int a, int b ) { return a > b ? a : b; }

//
//  saving parameters
//
const int NSTEPS = 1000;
const int SAVEFREQ = 10;

//
// particle data structure
//
typedef struct
{
  double x;  // x position
  double y;  // y position
  double vx; // velocity in the x direction
  double vy; // velocity in the y direction
  double ax; // acceleration in the x direction
  double ay; // acceleration in the y direction
} particle_t;

//
//  timing routines
//
double read_timer( );

//
//  simulation routines
//
void set_size( int n );
void init_particles( int n, particle_t *p );
void apply_force( particle_t &particle, particle_t &neighbor , double *dmin, double *davg, int *navg);
void move( particle_t &p );


//
//  I/O routines
//
FILE *open_save( char *filename, int n );
void save( FILE *f, int n, particle_t *p );

//
//  argument processing routines
//
int find_option( int argc, char **argv, const char *option );
int read_int( int argc, char **argv, const char *option, int default_value );
char *read_string( int argc, char **argv, const char *option, char *default_value );

//
//  our additions
//
typedef std::vector<particle_t> bin_t; // spatial hash bin data structure
extern double grid_size, bin_size;
extern int bin_count;
// Represent a 2D particle grid as a 1D vector of particle bins.
extern std::vector<bin_t> particle_bins;

// A function that is used to create a spatial hash
// given particles of type particle_t.
void make_spatial_hash(int n, particle_t* particles);

#endif
