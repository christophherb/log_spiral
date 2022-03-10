/* Automatically generated file. Do not edit. 
 * Format:     ANSI C source code
 * Creator:    McStas <http://www.mcstas.org>
 * Instrument: reverse_test.instr (logspir_test)
 * Date:       Thu Mar 10 10:18:10 2022
 * File:       ./reverse_test.c
 * Compile:    cc -o logspir_test.out ./reverse_test.c 
 * CFLAGS=
 */


#define MCCODE_STRING "McStas 2.7 - Nov. 27, 2020"
#define FLAVOR "mcstas"
#define FLAVOR_UPPER "MCSTAS"
#define MC_USE_DEFAULT_MAIN
#define MC_TRACE_ENABLED
#define MC_EMBEDDED_RUNTIME

#line 1 "mccode-r.h"
/*******************************************************************************
*
* McCode, neutron/xray ray-tracing package
*         Copyright (C) 1997-2009, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Runtime: share/mccode-r.h
*
* %Identification
* Written by: KN
* Date:    Aug 29, 1997
* Release: McStas 2.7
* Version: $Revision$
*
* Runtime system header for McStas/McXtrace.
*
* In order to use this library as an external library, the following variables
* and macros must be declared (see details in the code)
*
*   struct mcinputtable_struct mcinputtable[];
*   int mcnumipar;
*   char mcinstrument_name[], mcinstrument_source[];
*   int mctraceenabled, mcdefaultmain;
*   extern MCNUM  mccomp_storein[];
*   extern MCNUM  mcAbsorbProp[];
*   extern MCNUM  mcScattered;
*   #define MCCODE_STRING "the McStas/McXtrace version"
*
* Usage: Automatically embbeded in the c code.
*
* $Id$
*
*******************************************************************************/

#ifndef MCCODE_R_H
#define MCCODE_R_H "$Revision$"

#include <math.h>
#include <string.h>
#include <strings.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <limits.h>
#include <errno.h>
#include <time.h>
#include <float.h>
#include <inttypes.h>

/* If the runtime is embedded in the simulation program, some definitions can
   be made static. */

#ifdef MC_EMBEDDED_RUNTIME
#define mcstatic static
#else
#define mcstatic
#endif

#ifdef __dest_os
#if (__dest_os == __mac_os)
#define MAC
#endif
#endif

#ifdef __FreeBSD__
#define NEED_STAT_H
#endif

#if defined(__APPLE__) && defined(__GNUC__)
#define NEED_STAT_H
#endif

#ifdef NEED_STAT_H
#include <sys/stat.h>
#endif

#ifndef MC_PATHSEP_C
#ifdef WIN32
#define MC_PATHSEP_C '\\'
#define MC_PATHSEP_S "\\"
#else  /* !WIN32 */
#define MC_PATHSEP_C '/'
#define MC_PATHSEP_S "/"
#endif /* !WIN32 */
#endif /* MC_PATHSEP_C */

#ifndef WIN32
#ifndef _POSIX_SOURCE
#define _POSIX_SOURCE 1
#endif
#endif

/* the version string is replaced when building distribution with mkdist */
#ifndef MCCODE_STRING
#define MCCODE_STRING "McStas 2.7 - Nov. 27, 2020"
#endif

#ifndef MCCODE_DATE
#define MCCODE_DATE "Nov. 27, 2020"
#endif

#ifndef MCCODE_VERSION
#define MCCODE_VERSION "2.7"
#endif

#ifndef MCCODE_NAME
#define MCCODE_NAME "McStas"
#endif

#ifndef MCCODE_PARTICLE
#define MCCODE_PARTICLE "neutron"
#endif

#ifndef MCCODE_LIBENV
#define MCCODE_LIBENV "MCSTAS"
#endif

#ifndef FLAVOR_UPPER
#define FLAVOR_UPPER MCCODE_NAME
#endif

#ifdef MC_PORTABLE
#ifndef NOSIGNALS
#define NOSIGNALS 1
#endif
#endif

#ifdef MAC
#ifndef NOSIGNALS
#define NOSIGNALS 1
#endif
#endif

#if (USE_MPI == 0)
#undef USE_MPI
#endif

#ifdef USE_MPI  /* default is to disable signals with MPI, as MPICH uses them to communicate */
#ifndef NOSIGNALS
#define NOSIGNALS 1
#endif
#endif

#if (NOSIGNALS == 0)
#undef NOSIGNALS
#endif

/* Note: the enum instr_formal_types definition MUST be kept
   synchronized with the one in mccode.h and with the
   instr_formal_type_names array in cogen.c. */
enum instr_formal_types
  {
    instr_type_double, instr_type_int, instr_type_string
  };
struct mcinputtable_struct { /* defines instrument parameters */
  char *name; /* name of parameter */
  void *par;  /* pointer to instrument parameter (variable) */
  enum instr_formal_types type;
  char *val;  /* default value */
};

typedef double MCNUM;
typedef struct {MCNUM x, y, z;} Coords;
typedef MCNUM Rotation[3][3];

/* the following variables are defined in the McStas generated C code
   but should be defined externally in case of independent library usage */
#ifndef DANSE
extern struct mcinputtable_struct mcinputtable[]; /* list of instrument parameters */
extern int    mcnumipar;                          /* number of instrument parameters */
extern char   mcinstrument_name[], mcinstrument_source[]; /* instrument name and filename */
extern char  *mcinstrument_exe;                           /* executable path = argv[0] or NULL */
extern MCNUM  mccomp_storein[]; /* 11 coords * number of components in instrument */
extern MCNUM  mcAbsorbProp[];
extern MCNUM  mcScattered;      /* number of SCATTER calls in current component */
extern MCNUM  mcRestore;        /* Flag to indicate if neutron needs to be restored */
#ifndef MC_ANCIENT_COMPATIBILITY
extern int mctraceenabled, mcdefaultmain;
#endif
#endif


/* Useful macros ============================================================ */

/* MPI stuff */

#ifdef USE_MPI
#include "mpi.h"

#ifdef OMPI_MPI_H  /* openmpi does not use signals: we may install our sighandler */
#undef NOSIGNALS
#endif

/*
 * MPI_MASTER(i):
 * execution of i only on master node
 */
#define MPI_MASTER(statement) { \
  if(mpi_node_rank == mpi_node_root)\
  { statement; } \
}

#ifndef MPI_REDUCE_BLOCKSIZE
#define MPI_REDUCE_BLOCKSIZE 1000
#endif

int mc_MPI_Sum(double* buf, long count);
int mc_MPI_Send(void *sbuf, long count, MPI_Datatype dtype, int dest);
int mc_MPI_Recv(void *rbuf, long count, MPI_Datatype dtype, int source);

/* MPI_Finalize exits gracefully and should be preferred to MPI_Abort */
#define exit(code) do {                                   \
    MPI_Finalize();                                       \
    exit(code);                                           \
  } while(0)

#else /* !USE_MPI */
#define MPI_MASTER(instr) instr
#endif /* USE_MPI */

#ifdef USE_MPI
static int mpi_node_count;
#endif

#ifdef USE_THREADS  /* user want threads */
#error Threading (USE_THREADS) support has been removed for very poor efficiency. Use MPI/SSH grid instead.
#endif


void   mcset_ncount(unsigned long long count);    /* wrapper to get mcncount */
unsigned long long int mcget_ncount(void);            /* wrapper to set mcncount */
unsigned long long mcget_run_num(void);           /* wrapper to get mcrun_num=0:mcncount */


/* Following part is only embedded when not redundant with mccode.h ========= */

#ifndef MCCODE_H

#ifndef NOSIGNALS
#include <signal.h>
#define SIG_MESSAGE(msg) strcpy(mcsig_message, msg);
#else
#define SIG_MESSAGE(msg)
#endif /* !NOSIGNALS */

/* Useful macros and constants ============================================== */

#ifndef FLT_MAX
#define FLT_MAX         3.40282347E+38F /* max decimal value of a "float" */
#endif

#ifndef MIN
#define MIN(a, b)  (((a) < (b)) ? (a) : (b))
#endif
#ifndef MAX
#define MAX(a, b)  (((a) > (b)) ? (a) : (b))
#endif
#ifndef SQR
#define SQR(x) ( (x) * (x) )
#endif
#ifndef SIGN
#define SIGN(x) (((x)>0.0)?(1):(-1))
#endif

#ifndef PI
# ifdef M_PI
#  define PI M_PI
# else
/* When using c99 in the CFLAGS, some of these consts
   are lost... Perhaps we should in fact include everything from
   https://www.gnu.org/software/libc/manual/html_node/Mathematical-Constants.html
*/
#  define PI 3.14159265358979323846
#  define M_PI PI
#  define M_PI_2 M_PI/2.0
#  define M_PI_4 M_PI/4.0
#  define M_1_PI 1.0/M_PI
#  define M_2_PI 2*M_1_PI
#  define M_2_SQRTPI 2/sqrt(M_PI)
#  define M_SQRT2 sqrt(2)
#  define M_SQRT1_2 sqrt(1/2)
# endif
#endif

#define RAD2MIN  ((180*60)/PI)
#define MIN2RAD  (PI/(180*60))
#define DEG2RAD  (PI/180)
#define RAD2DEG  (180/PI)
#define FWHM2RMS 0.424660900144    /* Convert between full-width-half-max and */
#define RMS2FWHM 2.35482004503     /* root-mean-square (standard deviation) */
#define HBAR     1.05457168e-34    /* [Js] h bar Planck constant CODATA 2002 */
#define MNEUTRON 1.67492728e-27    /* [kg] mass of neutron CODATA 2002 */
#define GRAVITY  9.81              /* [m/s^2] gravitational acceleration */
#define NA       6.02214179e23     /* [#atoms/g .mole] Avogadro's number*/


/* wrapper to get absolute and relative position of comp */
/* mccomp_posa and mccomp_posr are defined in McStas generated C code */
#define POS_A_COMP_INDEX(index) \
    (mccomp_posa[index])
#define POS_R_COMP_INDEX(index) \
    (mccomp_posr[index])
/* number of SCATTER calls in current comp: mcScattered defined in generated C code */
#define SCATTERED mcScattered
/* Flag to indicate if neutron needs to be restored: mcRestore defined in generated C code */
#define RESTORE mcRestore


/* Retrieve component information from the kernel */
/* Name, position and orientation (both absolute and relative)  */
/* Any component: For "redundancy", see comment by KN */
#define tmp_name_comp(comp) #comp
#define NAME_COMP(comp) tmp_name_comp(comp)
#define tmp_pos_a_comp(comp) (mcposa ## comp)
#define POS_A_COMP(comp) tmp_pos_a_comp(comp)
#define tmp_pos_r_comp(comp) (mcposr ## comp)
#define POS_R_COMP(comp) tmp_pos_r_comp(comp)
#define tmp_rot_a_comp(comp) (mcrota ## comp)
#define ROT_A_COMP(comp) tmp_rot_a_comp(comp)
#define tmp_rot_r_comp(comp) (mcrotr ## comp)
#define ROT_R_COMP(comp) tmp_rot_r_comp(comp)

/* Current component name, index, position and orientation */
#define NAME_CURRENT_COMP  NAME_COMP(mccompcurname)
#define INDEX_CURRENT_COMP mccompcurindex
#define POS_A_CURRENT_COMP POS_A_COMP(mccompcurname)
#define POS_R_CURRENT_COMP POS_R_COMP(mccompcurname)
#define ROT_A_CURRENT_COMP ROT_A_COMP(mccompcurname)
#define ROT_R_CURRENT_COMP ROT_R_COMP(mccompcurname)

/* Note: The two-stage approach to MC_GETPAR is NOT redundant; without it,
* after #define C sample, MC_GETPAR(C,x) would refer to component C, not to
* component sample. Such are the joys of ANSI C.

* Anyway the usage of MCGETPAR requires that we use sometimes bare names...
*/
#define MC_GETPAR2(comp, par) (mcc ## comp ## _ ## par)
#define MC_GETPAR(comp, par) MC_GETPAR2(comp,par)

/* MCDISPLAY/trace and debugging message sent to stdout */
#ifdef MC_TRACE_ENABLED
#define DEBUG
#endif

#ifdef DEBUG
#define mcDEBUG_INSTR() if(!mcdotrace); else { printf("\nINSTRUMENT:\n"); printf("Instrument '%s' (%s)\n", mcinstrument_name, mcinstrument_source); }
#define mcDEBUG_COMPONENT(name,c,t) if(!mcdotrace); else {\
  printf("COMPONENT: \"%s\"\n" \
         "POS: %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n", \
         name, c.x, c.y, c.z, t[0][0], t[0][1], t[0][2], \
         t[1][0], t[1][1], t[1][2], t[2][0], t[2][1], t[2][2]); \
  mcAccumulatedILength += coords_len(coords_sub(mcLastComp,c)); \
  printf("Component %30s AT (%g,%g,%g)    %g m from origin\n", name, c.x, c.y, c.z, mcAccumulatedILength); \
  mcLastComp=c;\
  }
#define mcDEBUG_INSTR_END() if(!mcdotrace); else printf("INSTRUMENT END:\n");
#define mcDEBUG_ENTER() if(!mcdotrace); else printf("ENTER:\n");
#define mcDEBUG_COMP(c) if(!mcdotrace); else printf("COMP: \"%s\"\n", c);
#define mcDEBUG_LEAVE() if(!mcdotrace); else printf("LEAVE:\n");
#define mcDEBUG_ABSORB() if(!mcdotrace); else printf("ABSORB:\n");
#else
#define mcDEBUG_INSTR()
#define mcDEBUG_COMPONENT(name,c,t)
#define mcDEBUG_INSTR_END()
#define mcDEBUG_ENTER()
#define mcDEBUG_COMP(c)
#define mcDEBUG_LEAVE()
#define mcDEBUG_ABSORB()
#endif

// mcDEBUG_STATE and mcDEBUG_SCATTER are defined by mcstas-r.h and mcxtrace-r.h



#ifdef TEST
#define test_printf printf
#else
#define test_printf while(0) printf
#endif

/* send MCDISPLAY message to stdout to show gemoetry */
void mcdis_magnify(char *what);
void mcdis_line(double x1, double y1, double z1,
                double x2, double y2, double z2);
void mcdis_dashed_linemcdis_dashed_line(double x1, double y1, double z1,
		       double x2, double y2, double z2, int n);
void mcdis_multiline(int count, ...);
void mcdis_rectangle(char* plane, double x, double y, double z,
		     double width, double height);
void mcdis_box(double x, double y, double z,
	       double width, double height, double length);
void mcdis_circle(char *plane, double x, double y, double z, double r);
void mcdis_Circle(double x, double y, double z, double r, double nx, double ny, double nz);
void mcdis_cylinder( double x, double y, double z,
        double r, double height, int N, double nx, double ny, double nz);
void mcdis_sphere(double x, double y, double z, double r, int N);

/* selection of random number generator. default is MT */
#ifndef MC_RAND_ALG
#define MC_RAND_ALG 1
#endif

#if MC_RAND_ALG == 0
   /* Use system random() (not recommended). */
#  define MC_RAND_MAX RAND_MAX
#elif MC_RAND_ALG == 1
   /* "Mersenne Twister", by Makoto Matsumoto and Takuji Nishimura. */
#  define MC_RAND_MAX ((unsigned long)0xffffffff)
#  define random mt_random
#  define srandom mt_srandom
#elif MC_RAND_ALG == 2
   /* Algorithm used in McStas CVS-080208 and earlier (not recommended). */
#  define MC_RAND_MAX 0x7fffffff
#  define random mc_random
#  define srandom mc_srandom
#else
#  error "Bad value for random number generator choice."
#endif

typedef int mc_int32_t;
mc_int32_t mc_random(void);
void mc_srandom (unsigned int x);
unsigned long mt_random(void);
void mt_srandom (unsigned long x);

double rand01();
double randpm1();
double rand0max(double max);
double randminmax(double min, double max);

double randnorm(void);
double randtriangle(void);

#ifndef DANSE
void mcinit(void);
void mcraytrace(void);
void mcsave(FILE *);
void mcfinally(void);
void mcdisplay(void);
#endif

/* simple vector algebra ==================================================== */
#define vec_prod(x, y, z, x1, y1, z1, x2, y2, z2) \
	vec_prod_func(&x, &y, &z, x1, y1, z1, x2, y2, z2)
mcstatic void vec_prod_func(double *x, double *y, double *z,
		double x1, double y1, double z1, double x2, double y2, double z2);

mcstatic double scalar_prod(
		double x1, double y1, double z1, double x2, double y2, double z2);

#define NORM(x,y,z) \
	norm_func(&x, &y, &z)
mcstatic void norm_func(double *x, double *y, double *z) {
	double temp = (*x * *x) + (*y * *y) + (*z * *z);
	if (temp != 0) {
		temp = sqrt(temp);
		*x /= temp;
		*y /= temp;
		*z /= temp;
	}
}
#define normal_vec(nx, ny, nz, x, y, z) \
    normal_vec_func(&(nx), &(ny), &(nz), x, y, z)
mcstatic void normal_vec_func(double *nx, double *ny, double *nz,
    double x, double y, double z);

/**
 * Rotate the vector vx,vy,vz psi radians around the vector ax,ay,az
 * and put the result in x,y,z.
 */
#define rotate(x, y, z, vx, vy, vz, phi, ax, ay, az) \
  do { \
    double mcrt_tmpx = (ax), mcrt_tmpy = (ay), mcrt_tmpz = (az); \
    double mcrt_vp, mcrt_vpx, mcrt_vpy, mcrt_vpz; \
    double mcrt_vnx, mcrt_vny, mcrt_vnz, mcrt_vn1x, mcrt_vn1y, mcrt_vn1z; \
    double mcrt_bx, mcrt_by, mcrt_bz; \
    double mcrt_cos, mcrt_sin; \
    NORM(mcrt_tmpx, mcrt_tmpy, mcrt_tmpz); \
    mcrt_vp = scalar_prod((vx), (vy), (vz), mcrt_tmpx, mcrt_tmpy, mcrt_tmpz); \
    mcrt_vpx = mcrt_vp*mcrt_tmpx; \
    mcrt_vpy = mcrt_vp*mcrt_tmpy; \
    mcrt_vpz = mcrt_vp*mcrt_tmpz; \
    mcrt_vnx = (vx) - mcrt_vpx; \
    mcrt_vny = (vy) - mcrt_vpy; \
    mcrt_vnz = (vz) - mcrt_vpz; \
    vec_prod(mcrt_bx, mcrt_by, mcrt_bz, \
             mcrt_tmpx, mcrt_tmpy, mcrt_tmpz, mcrt_vnx, mcrt_vny, mcrt_vnz); \
    mcrt_cos = cos((phi)); mcrt_sin = sin((phi)); \
    mcrt_vn1x = mcrt_vnx*mcrt_cos + mcrt_bx*mcrt_sin; \
    mcrt_vn1y = mcrt_vny*mcrt_cos + mcrt_by*mcrt_sin; \
    mcrt_vn1z = mcrt_vnz*mcrt_cos + mcrt_bz*mcrt_sin; \
    (x) = mcrt_vpx + mcrt_vn1x; \
    (y) = mcrt_vpy + mcrt_vn1y; \
    (z) = mcrt_vpz + mcrt_vn1z; \
  } while(0)

/**
 * Mirror (xyz) in the plane given by the point (rx,ry,rz) and normal (nx,ny,nz)
 *
 * TODO: This define is seemingly never used...
 */
#define mirror(x,y,z,rx,ry,rz,nx,ny,nz) \
  do { \
    double mcrt_tmpx= (nx), mcrt_tmpy = (ny), mcrt_tmpz = (nz); \
    double mcrt_tmpt; \
    NORM(mcrt_tmpx, mcrt_tmpy, mcrt_tmpz); \
    mcrt_tmpt=scalar_prod((rx),(ry),(rz),mcrt_tmpx,mcrt_tmpy,mcrt_tmpz); \
    (x) = rx -2 * mcrt_tmpt*mcrt_rmpx; \
    (y) = ry -2 * mcrt_tmpt*mcrt_rmpy; \
    (z) = rz -2 * mcrt_tmpt*mcrt_rmpz; \
  } while (0)

Coords coords_set(MCNUM x, MCNUM y, MCNUM z);
Coords coords_get(Coords a, MCNUM *x, MCNUM *y, MCNUM *z);
Coords coords_add(Coords a, Coords b);
Coords coords_sub(Coords a, Coords b);
Coords coords_neg(Coords a);
Coords coords_scale(Coords b, double scale);
double coords_sp(Coords a, Coords b);
Coords coords_xp(Coords b, Coords c);
double coords_len(Coords a);
void   coords_print(Coords a);
mcstatic void coords_norm(Coords* c);

void rot_set_rotation(Rotation t, double phx, double phy, double phz);
int  rot_test_identity(Rotation t);
void rot_mul(Rotation t1, Rotation t2, Rotation t3);
void rot_copy(Rotation dest, Rotation src);
void rot_transpose(Rotation src, Rotation dst);
Coords rot_apply(Rotation t, Coords a);

void mccoordschange(Coords a, Rotation t, double *x, double *y, double *z,
    double *vx, double *vy, double *vz, double *sx, double *sy, double *sz);
void
mccoordschange_polarisation(Rotation t, double *sx, double *sy, double *sz);

double mcestimate_error(double N, double p1, double p2);
void mcreadparams(void);

/* this is now in mcstas-r.h and mcxtrace-r.h as the number of state parameters is no longer equal*/
/* void mcsetstate(double x, double y, double z, double vx, double vy, double vz,
                double t, double sx, double sy, double sz, double p);
*/
void mcgenstate(void);

/* trajectory/shape intersection routines */
int inside_rectangle(double, double, double, double);
int box_intersect(double *dt_in, double *dt_out, double x, double y, double z,
    double vx, double vy, double vz, double dx, double dy, double dz);
int cylinder_intersect(double *t0, double *t1, double x, double y, double z,
    double vx, double vy, double vz, double r, double h);
int sphere_intersect(double *t0, double *t1, double x, double y, double z,
                 double vx, double vy, double vz, double r);
/* second order equation roots */
int solve_2nd_order(double *t1, double *t2,
    double A,  double B,  double C);

/* random vector generation to shape */
void randvec_target_circle(double *xo, double *yo, double *zo,
    double *solid_angle, double xi, double yi, double zi, double radius);
#define randvec_target_sphere randvec_target_circle
void randvec_target_rect_angular(double *xo, double *yo, double *zo,
    double *solid_angle,
               double xi, double yi, double zi, double height, double width, Rotation A);
#define randvec_target_rect(p0,p1,p2,p3,p4,p5,p6,p7,p8,p9)  randvec_target_rect_real(p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,0,0,0,1)
void randvec_target_rect_real(double *xo, double *yo, double *zo,
    double *solid_angle,
	       double xi, double yi, double zi, double height, double width, Rotation A,
			 double lx, double ly, double lz, int order);

/* this is the main() */
int mccode_main(int argc, char *argv[]);


#endif /* !MCCODE_H */

#ifndef MCCODE_R_IO_H
#define MCCODE_R_IO_H "$Revision$"

#if (USE_NEXUS == 0)
#undef USE_NEXUS
#endif

#ifndef CHAR_BUF_LENGTH
#define CHAR_BUF_LENGTH 1024
#endif

/* I/O section part ========================================================= */

/* ========================================================================== */

/*                               MCCODE_R_IO_C                                */

/* ========================================================================== */


/* main DETECTOR structure which stores most information to write to data files */
struct mcdetector_struct {
  char   filename[CHAR_BUF_LENGTH];   /* file name of monitor */
  char   position[CHAR_BUF_LENGTH];   /* position of detector component */
  char   component[CHAR_BUF_LENGTH];  /* component instance name */
  char   instrument[CHAR_BUF_LENGTH]; /* instrument name */
  char   type[CHAR_BUF_LENGTH];       /* data type, e.g. 0d, 1d, 2d, 3d */
  char   user[CHAR_BUF_LENGTH];       /* user name, e.g. HOME */
  char   date[CHAR_BUF_LENGTH];       /* date of simulation end/write time */
  char   title[CHAR_BUF_LENGTH];      /* title of detector */
  char   xlabel[CHAR_BUF_LENGTH];     /* X axis label */
  char   ylabel[CHAR_BUF_LENGTH];     /* Y axis label */
  char   zlabel[CHAR_BUF_LENGTH];     /* Z axis label */
  char   xvar[CHAR_BUF_LENGTH];       /* X variable name */
  char   yvar[CHAR_BUF_LENGTH];       /* Y variable name */
  char   zvar[CHAR_BUF_LENGTH];       /* Z variable name */
  char   ncount[CHAR_BUF_LENGTH];     /* number of events initially generated */
  char   limits[CHAR_BUF_LENGTH];     /* X Y Z limits, e.g. [xmin xmax ymin ymax zmin zmax] */
  char   variables[CHAR_BUF_LENGTH];  /* variables written into data block */
  char   statistics[CHAR_BUF_LENGTH]; /* center, mean and half width along axis */
  char   signal[CHAR_BUF_LENGTH];     /* min max and mean of signal (data block) */
  char   values[CHAR_BUF_LENGTH];     /* integrated values e.g. [I I_err N] */
  double xmin,xmax;                   /* min max of axes */
  double ymin,ymax;
  double zmin,zmax;
  double intensity;                   /* integrated values for data block */
  double error;
  double events;
  double min;                         /* statistics for data block */
  double max;
  double mean;
  double centerX;                     /* statistics for axes */
  double halfwidthX;
  double centerY;
  double halfwidthY;
  int    rank;                        /* dimensionaly of monitor, e.g. 0 1 2 3 */
  char   istransposed;                /* flag to transpose matrix for some formats */

  long   m,n,p;                       /* dimensions of data block and along axes */
  long   date_l;                      /* same as date, but in sec since 1970 */

  double *p0, *p1, *p2;               /* pointers to saved data, NULL when freed */
  char   format[CHAR_BUF_LENGTH];    /* format for file generation */
};

typedef struct mcdetector_struct MCDETECTOR;

static   char *mcdirname             = NULL;      /* name of output directory */
static   char *mcsiminfo_name        = "mccode";  /* default output sim file name */
char    *mcformat                    = NULL;      /* NULL (default) or a specific format */

/* file I/O definitions and function prototypes */

#ifndef MC_EMBEDDED_RUNTIME /* the mcstatic variables (from mccode-r.c) */
extern FILE * mcsiminfo_file;     /* handle to the output siminfo file */
extern int    mcgravitation;      /* flag to enable gravitation */
extern int    mcdotrace;          /* flag to print MCDISPLAY messages */
#else
mcstatic FILE *mcsiminfo_file        = NULL;
#endif

/* I/O function prototypes ================================================== */

/* output functions */
MCDETECTOR mcdetector_out_0D(char *t, double p0, double p1, double p2, char *c, Coords pos);
MCDETECTOR mcdetector_out_1D(char *t, char *xl, char *yl,
                  char *xvar, double x1, double x2, long n,
                  double *p0, double *p1, double *p2, char *f, char *c, Coords pos);
MCDETECTOR mcdetector_out_2D(char *t, char *xl, char *yl,
                  double x1, double x2, double y1, double y2, long m,
                  long n, double *p0, double *p1, double *p2, char *f,
                  char *c, Coords pos);
MCDETECTOR mcdetector_out_list(char *t, char *xl, char *yl,
                  long m, long n,
                  double *p1, char *f,
                  char *c, Coords posa);

/* wrappers to output functions, that automatically set NAME and POSITION */
#define DETECTOR_OUT(p0,p1,p2) mcdetector_out_0D(NAME_CURRENT_COMP,p0,p1,p2,NAME_CURRENT_COMP,POS_A_CURRENT_COMP)
#define DETECTOR_OUT_0D(t,p0,p1,p2) mcdetector_out_0D(t,p0,p1,p2,NAME_CURRENT_COMP,POS_A_CURRENT_COMP)
#define DETECTOR_OUT_1D(t,xl,yl,xvar,x1,x2,n,p0,p1,p2,f) \
     mcdetector_out_1D(t,xl,yl,xvar,x1,x2,n,p0,p1,p2,f,NAME_CURRENT_COMP,POS_A_CURRENT_COMP)
#define DETECTOR_OUT_2D(t,xl,yl,x1,x2,y1,y2,m,n,p0,p1,p2,f) \
     mcdetector_out_2D(t,xl,yl,x1,x2,y1,y2,m,n,p0,p1,p2,f,NAME_CURRENT_COMP,POS_A_CURRENT_COMP)

#ifdef USE_NEXUS
#include "napi.h"
NXhandle nxhandle;
#endif

#endif /* ndef MCCODE_R_IO_H */

#endif /* MCCODE_R_H */
/* End of file "mccode-r.h". */

#line 712 "./reverse_test.c"

#line 1 "mcstas-r.h"
/*******************************************************************************
*
* McStas, neutron ray-tracing package
*         Copyright (C) 1997-2009, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Runtime: share/mcstas-r.h
*
* %Identification
* Written by: KN
* Date:    Aug 29, 1997
* Release: McStas X.Y
* Version: $Revision$
*
* Runtime system header for McStas.
*
* In order to use this library as an external library, the following variables
* and macros must be declared (see details in the code)
*
*   struct mcinputtable_struct mcinputtable[];
*   int mcnumipar;
*   char mcinstrument_name[], mcinstrument_source[];
*   int mctraceenabled, mcdefaultmain;
*   extern MCNUM  mccomp_storein[];
*   extern MCNUM  mcAbsorbProp[];
*   extern MCNUM  mcScattered;
*   #define MCCODE_STRING "the McStas version"
*
* Usage: Automatically embbeded in the c code.
*
* $Id$
*
*******************************************************************************/

#ifndef MCSTAS_R_H
#define MCSTAS_R_H "$Revision$"

/* Following part is only embedded when not redundant with mcstas.h ========= */

#ifndef MCCODE_H

#define AA2MS    629.622368        /* Convert k[1/AA] to v[m/s] */
#define MS2AA    1.58825361e-3     /* Convert v[m/s] to k[1/AA] */
#define K2V      AA2MS
#define V2K      MS2AA
#define Q2V      AA2MS
#define V2Q      MS2AA
#define SE2V     437.393377        /* Convert sqrt(E)[meV] to v[m/s] */
#define VS2E     5.22703725e-6     /* Convert (v[m/s])**2 to E[meV] */

#define SCATTER do {mcDEBUG_SCATTER(mcnlx, mcnly, mcnlz, mcnlvx, mcnlvy, mcnlvz, \
        mcnlt,mcnlsx,mcnlsy,mcnlsz, mcnlp); mcScattered++;} while(0)
#define ABSORB do {mcDEBUG_STATE(mcnlx, mcnly, mcnlz, mcnlvx, mcnlvy, mcnlvz, \
        mcnlt,mcnlsx,mcnlsy,mcnlsz, mcnlp); mcDEBUG_ABSORB(); MAGNET_OFF; goto mcabsorb;} while(0)

#define STORE_NEUTRON(index, x, y, z, vx, vy, vz, t, sx, sy, sz, p) \
  mcstore_neutron(mccomp_storein,index, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
#define RESTORE_NEUTRON(index, x, y, z, vx, vy, vz, t, sx, sy, sz, p) \
  mcrestore_neutron(mccomp_storein,index, &x, &y, &z, &vx, &vy, &vz, &t, &sx, &sy, &sz, &p);

#define MAGNET_ON \
  do { \
    mcMagnet = 1; \
  } while(0)

#define MAGNET_OFF \
  do { \
    mcMagnet = 0; \
  } while(0)

#define ALLOW_BACKPROP \
  do { \
    mcallowbackprop = 1; \
  } while(0)

#define DISALLOW_BACKPROP \
  do { \
    mcallowbackprop = 0; \
  } while(0)

#define PROP_MAGNET(dt) \
  do { \
  }while (0)
    /* change coordinates from local system to magnet system */
/*    Rotation rotLM, rotTemp; \
      Coords   posLM = coords_sub(POS_A_CURRENT_COMP, mcMagnetPos); \
      rot_transpose(ROT_A_CURRENT_COMP, rotTemp); \
      rot_mul(rotTemp, mcMagnetRot, rotLM); \
      mcMagnetPrecession(mcnlx, mcnly, mcnlz, mcnlt, mcnlvx, mcnlvy, mcnlvz, \
               &mcnlsx, &mcnlsy, &mcnlsz, dt, posLM, rotLM); \
      } while(0)
*/

#define mcPROP_DT(dt) \
  do { \
    if (mcMagnet && dt > 0) PROP_MAGNET(dt);\
    mcnlx += mcnlvx*(dt); \
    mcnly += mcnlvy*(dt); \
    mcnlz += mcnlvz*(dt); \
    mcnlt += (dt); \
    if (isnan(p) || isinf(p)) { mcAbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }\
  } while(0)

/* ADD: E. Farhi, Aug 6th, 2001 PROP_GRAV_DT propagation with acceleration */
#define PROP_GRAV_DT(dt, Ax, Ay, Az) \
  do { \
    if(dt < 0 && mcallowbackprop == 0) { mcAbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }\
    if (mcMagnet) printf("Spin precession gravity\n"); \
    mcnlx  += mcnlvx*(dt) + (Ax)*(dt)*(dt)/2; \
    mcnly  += mcnlvy*(dt) + (Ay)*(dt)*(dt)/2; \
    mcnlz  += mcnlvz*(dt) + (Az)*(dt)*(dt)/2; \
    mcnlvx += (Ax)*(dt); \
    mcnlvy += (Ay)*(dt); \
    mcnlvz += (Az)*(dt); \
    mcnlt  += (dt); \
    DISALLOW_BACKPROP;\
  } while(0)


#define PROP_DT(dt) \
  do { \
    if(dt < 0) { RESTORE=1; goto mcabsorbComp; }; \
    if (mcgravitation) { Coords mcLocG; double mc_gx, mc_gy, mc_gz; \
    mcLocG = rot_apply(ROT_A_CURRENT_COMP, coords_set(0,-GRAVITY,0)); \
    coords_get(mcLocG, &mc_gx, &mc_gy, &mc_gz); \
    PROP_GRAV_DT(dt, mc_gx, mc_gy, mc_gz); } \
    else mcPROP_DT(dt); \
    DISALLOW_BACKPROP;\
  } while(0)


#define PROP_Z0 \
  do { \
    if (mcgravitation) { Coords mcLocG; int mc_ret; \
    double mc_dt, mc_gx, mc_gy, mc_gz; \
    mcLocG = rot_apply(ROT_A_CURRENT_COMP, coords_set(0,-GRAVITY,0)); \
    coords_get(mcLocG, &mc_gx, &mc_gy, &mc_gz); \
    mc_ret = solve_2nd_order(&mc_dt, NULL, -mc_gz/2, -mcnlvz, -mcnlz); \
    if (mc_ret && mc_dt>=0) {PROP_GRAV_DT(mc_dt, mc_gx, mc_gy, mc_gz); mcnlz=0;}\
    else { if (mcallowbackprop ==0) {mcAbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }}; }\
    else mcPROP_Z0; \
    DISALLOW_BACKPROP;\
  } while(0)

#define mcPROP_Z0 \
  do { \
    double mc_dt; \
    if(mcnlvz == 0) { mcAbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }; \
    mc_dt = -mcnlz/mcnlvz; \
    if(mc_dt < 0 && mcallowbackprop == 0) { mcAbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }; \
    mcPROP_DT(mc_dt); \
    mcnlz = 0; \
    DISALLOW_BACKPROP;\
  } while(0)

#define PROP_X0 \
  do { \
    if (mcgravitation) { Coords mcLocG; int mc_ret; \
    double mc_dt, mc_gx, mc_gy, mc_gz; \
    mcLocG = rot_apply(ROT_A_CURRENT_COMP, coords_set(0,-GRAVITY,0)); \
    coords_get(mcLocG, &mc_gx, &mc_gy, &mc_gz); \
    mc_ret = solve_2nd_order(&mc_dt, NULL, -mc_gx/2, -mcnlvx, -mcnlx); \
    if (mc_ret && mc_dt>=0) PROP_GRAV_DT(mc_dt, mc_gx, mc_gy, mc_gz); \
    else { if (mcallowbackprop ==0) {mcAbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }}; }\
    else mcPROP_X0; \
    DISALLOW_BACKPROP;\
  } while(0)

#define mcPROP_X0 \
  do { \
    double mc_dt; \
    if(mcnlvx == 0) { mcAbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }; \
    mc_dt = -mcnlx/mcnlvx; \
    if(mc_dt < 0 && mcallowbackprop == 0) { mcAbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }; \
    mcPROP_DT(mc_dt); \
    mcnlx = 0; \
    DISALLOW_BACKPROP;\
  } while(0)

#define PROP_Y0 \
  do { \
    if (mcgravitation) { Coords mcLocG; int mc_ret; \
    double mc_dt, mc_gx, mc_gy, mc_gz; \
    mcLocG = rot_apply(ROT_A_CURRENT_COMP, coords_set(0,-GRAVITY,0)); \
    coords_get(mcLocG, &mc_gx, &mc_gy, &mc_gz); \
    mc_ret = solve_2nd_order(&mc_dt, NULL, -mc_gy/2, -mcnlvy, -mcnly); \
    if (mc_ret && mc_dt>=0) PROP_GRAV_DT(mc_dt, mc_gx, mc_gy, mc_gz); \
    else { if (mcallowbackprop ==0) {mcAbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }}; }\
    else mcPROP_Y0; \
    DISALLOW_BACKPROP;\
  } while(0)


#define mcPROP_Y0 \
  do { \
    double mc_dt; \
    if(mcnlvy == 0) { mcAbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }; \
    mc_dt = -mcnly/mcnlvy; \
    if(mc_dt < 0 && mcallowbackprop == 0) { mcAbsorbProp[INDEX_CURRENT_COMP]++; ABSORB; }; \
    mcPROP_DT(mc_dt); \
    mcnly = 0; \
    DISALLOW_BACKPROP; \
  } while(0)

/*moved from mccode-r.h*/
void mcsetstate(double x, double y, double z, double vx, double vy, double vz,
                double t, double sx, double sy, double sz, double p);

#ifdef DEBUG

#define mcDEBUG_STATE(x,y,z,vx,vy,vz,t,sx,sy,sz,p) if(!mcdotrace); else \
  printf("STATE: %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n", \
         x,y,z,vx,vy,vz,t,sx,sy,sz,p);
#define mcDEBUG_SCATTER(x,y,z,vx,vy,vz,t,sx,sy,sz,p) if(!mcdotrace); else \
  printf("SCATTER: %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n", \
         x,y,z,vx,vy,vz,t,sx,sy,sz,p);

#else

#define mcDEBUG_STATE(x,y,z,vx,vy,vz,t,sx,sy,sz,p)
#define mcDEBUG_SCATTER(x,y,z,vx,vy,vz,t,sx,sy,sz,p)

#endif

#endif /* !MCCODE_H */

#endif /* MCSTAS_R_H */
/* End of file "mcstas-r.h". */

#line 945 "./reverse_test.c"

#line 1 "mccode-r.c"
/*******************************************************************************
*
* McCode, neutron/xray ray-tracing package
*         Copyright (C) 1997-2009, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Runtime: share/mccode-r.c
*
* %Identification
* Written by: KN
* Date:    Aug 29, 1997
* Release: McStas X.Y/McXtrace X.Y
* Version: $Revision$
*
* Runtime system for McStas and McXtrace.
* Embedded within instrument in runtime mode.
* Contains SECTIONS:
*   MPI handling (sum, send, recv)
*   format definitions
*   I/O
*   mcdisplay support
*   random numbers
*   coordinates handling
*   vectors math (solve 2nd order, normals, randvec...)
*   parameter handling
*   signal and main handlers
*
* Usage: Automatically embbeded in the c code whenever required.
*
* $Id$
*
*******************************************************************************/

/*******************************************************************************
* The I/O format definitions and functions
*******************************************************************************/


/** Include header files to avoid implicit declarations (not allowed on LLVM) */
#include <ctype.h>
#include <sys/types.h>

// UNIX specific headers (non-Windows)
#if defined(__unix__) || defined(__APPLE__)
#include <unistd.h>
#endif

#include <sys/stat.h>

#ifdef _WIN32 
#include <direct.h>
# define  mkdir( D, M )   _mkdir( D ) 
#endif 

#ifndef DANSE
#ifdef MC_ANCIENT_COMPATIBILITY
int mctraceenabled = 0;
int mcdefaultmain  = 0;
#endif
/* else defined directly in the McCode generated C code */

static   long mcseed                 = 0; /* seed for random generator */
static   long mcstartdate            = 0; /* start simulation time */
static   int  mcdisable_output_files = 0; /* --no-output-files */
mcstatic int  mcgravitation          = 0; /* use gravitation flag, for PROP macros */
int      mcMagnet                    = 0; /* magnet stack flag */
mcstatic int  mcdotrace              = 0; /* flag for --trace and messages for DISPLAY */
int      mcallowbackprop             = 0;         /* flag to enable negative/backprop */

/* Number of particle histories to simulate. */
#ifdef NEUTRONICS
mcstatic unsigned long long int mcncount             = 1;
mcstatic unsigned long long int mcrun_num            = 0;
#else
mcstatic unsigned long long int mcncount             = 1000000;
mcstatic unsigned long long int mcrun_num            = 0;
#endif /* NEUTRONICS */

#else
#include "mcstas-globals.h"
#endif /* !DANSE */

/* SECTION: MPI handling ==================================================== */

#ifdef USE_MPI
/* MPI rank */
static int mpi_node_rank;
static int mpi_node_root = 0;


/*******************************************************************************
* mc_MPI_Reduce: Gathers arrays from MPI nodes using Reduce function.
*******************************************************************************/
int mc_MPI_Sum(double *sbuf, long count)
{
  if (!sbuf || count <= 0) return(MPI_SUCCESS); /* nothing to reduce */
  else {
    /* we must cut the buffer into blocks not exceeding the MPI max buffer size of 32000 */
    long   offset=0;
    double *rbuf=NULL;
    int    length=MPI_REDUCE_BLOCKSIZE; /* defined in mccode-r.h */
    int    i=0;
    rbuf = calloc(count, sizeof(double));
    if (!rbuf)
      exit(-fprintf(stderr, "Error: Out of memory %li (mc_MPI_Sum)\n", count*sizeof(double)));
    while (offset < count) {
      if (!length || offset+length > count-1) length=count-offset;
      else length=MPI_REDUCE_BLOCKSIZE;
      if (MPI_Reduce((double*)(sbuf+offset), (double*)(rbuf+offset),
              length, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
        return MPI_ERR_COUNT;
      offset += length;
    }

    for (i=0; i<count; i++) sbuf[i] = rbuf[i];
    free(rbuf);
  }
  return MPI_SUCCESS;
} /* mc_MPI_Sum */

/*******************************************************************************
* mc_MPI_Send: Send array to MPI node by blocks to avoid buffer limit
*******************************************************************************/
int mc_MPI_Send(void *sbuf,
                  long count, MPI_Datatype dtype,
                  int dest)
{
  int dsize;
  long offset=0;
  int  tag=1;
  int  length=MPI_REDUCE_BLOCKSIZE; /* defined in mccode-r.h */

  if (!sbuf || count <= 0) return(MPI_SUCCESS); /* nothing to send */
  MPI_Type_size(dtype, &dsize);

  while (offset < count) {
    if (offset+length > count-1) length=count-offset;
    else length=MPI_REDUCE_BLOCKSIZE;
    if (MPI_Send((void*)(sbuf+offset*dsize), length, dtype, dest, tag++, MPI_COMM_WORLD) != MPI_SUCCESS)
      return MPI_ERR_COUNT;
    offset += length;
  }

  return MPI_SUCCESS;
} /* mc_MPI_Send */

/*******************************************************************************
* mc_MPI_Recv: Receives arrays from MPI nodes by blocks to avoid buffer limit
*             the buffer must have been allocated previously.
*******************************************************************************/
int mc_MPI_Recv(void *sbuf,
                  long count, MPI_Datatype dtype,
                  int source)
{
  int dsize;
  long offset=0;
  int  tag=1;
  int  length=MPI_REDUCE_BLOCKSIZE; /* defined in mccode-r.h */

  if (!sbuf || count <= 0) return(MPI_SUCCESS); /* nothing to recv */
  MPI_Type_size(dtype, &dsize);

  while (offset < count) {
    if (offset+length > count-1) length=count-offset;
    else length=MPI_REDUCE_BLOCKSIZE;
    if (MPI_Recv((void*)(sbuf+offset*dsize), length, dtype, source, tag++,
            MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS)
      return MPI_ERR_COUNT;
    offset += length;
  }

  return MPI_SUCCESS;
} /* mc_MPI_Recv */

#endif /* USE_MPI */

/* SECTION: parameters handling ============================================= */

/* Instrument input parameter type handling. */
/*******************************************************************************
* mcparm_double: extract double value from 's' into 'vptr'
*******************************************************************************/
static int
mcparm_double(char *s, void *vptr)
{
  char *p;
  double *v = (double *)vptr;

  if (!s) { *v = 0; return(1); }
  *v = strtod(s, &p);
  if(*s == '\0' || (p != NULL && *p != '\0') || errno == ERANGE)
    return 0;                        /* Failed */
  else
    return 1;                        /* Success */
}

/*******************************************************************************
* mcparminfo_double: display parameter type double
*******************************************************************************/
static char *
mcparminfo_double(char *parmname)
{
  return "double";
}

/*******************************************************************************
* mcparmerror_double: display error message when failed extract double
*******************************************************************************/
static void
mcparmerror_double(char *parm, char *val)
{
  fprintf(stderr, "Error: Invalid value '%s' for floating point parameter %s (mcparmerror_double)\n",
          val, parm);
}

/*******************************************************************************
* mcparmprinter_double: convert double to string
*******************************************************************************/
static void
mcparmprinter_double(char *f, void *vptr)
{
  double *v = (double *)vptr;
  sprintf(f, "%g", *v);
}

/*******************************************************************************
* mcparm_int: extract int value from 's' into 'vptr'
*******************************************************************************/
static int
mcparm_int(char *s, void *vptr)
{
  char *p;
  int *v = (int *)vptr;
  long x;

  if (!s) { *v = 0; return(1); }
  *v = 0;
  x = strtol(s, &p, 10);
  if(x < INT_MIN || x > INT_MAX)
    return 0;                        /* Under/overflow */
  *v = x;
  if(*s == '\0' || (p != NULL && *p != '\0') || errno == ERANGE)
    return 0;                        /* Failed */
  else
    return 1;                        /* Success */
}

/*******************************************************************************
* mcparminfo_int: display parameter type int
*******************************************************************************/
static char *
mcparminfo_int(char *parmname)
{
  return "int";
}

/*******************************************************************************
* mcparmerror_int: display error message when failed extract int
*******************************************************************************/
static void
mcparmerror_int(char *parm, char *val)
{
  fprintf(stderr, "Error: Invalid value '%s' for integer parameter %s (mcparmerror_int)\n",
          val, parm);
}

/*******************************************************************************
* mcparmprinter_int: convert int to string
*******************************************************************************/
static void
mcparmprinter_int(char *f, void *vptr)
{
  int *v = (int *)vptr;
  sprintf(f, "%d", *v);
}

/*******************************************************************************
* mcparm_string: extract char* value from 's' into 'vptr' (copy)
*******************************************************************************/
static int
mcparm_string(char *s, void *vptr)
{
  char **v = (char **)vptr;
  if (!s) { *v = NULL; return(1); }
  *v = (char *)malloc(strlen(s) + 1);
  if(*v == NULL)
  {
    exit(-fprintf(stderr, "Error: Out of memory %li (mcparm_string).\n", (long)strlen(s) + 1));
  }
  strcpy(*v, s);
  return 1;                        /* Success */
}

/*******************************************************************************
* mcparminfo_string: display parameter type string
*******************************************************************************/
static char *
mcparminfo_string(char *parmname)
{
  return "string";
}

/*******************************************************************************
* mcparmerror_string: display error message when failed extract string
*******************************************************************************/
static void
mcparmerror_string(char *parm, char *val)
{
  fprintf(stderr, "Error: Invalid value '%s' for string parameter %s (mcparmerror_string)\n",
          val, parm);
}

/*******************************************************************************
* mcparmprinter_string: convert string to string (including esc chars)
*******************************************************************************/
static void
mcparmprinter_string(char *f, void *vptr)
{
  char **v = (char **)vptr;
  char *p;

  if (!*v) { *f='\0'; return; }
  strcpy(f, "");
  for(p = *v; *p != '\0'; p++)
  {
    switch(*p)
    {
      case '\n':
        strcat(f, "\\n");
        break;
      case '\r':
        strcat(f, "\\r");
        break;
      case '"':
        strcat(f, "\\\"");
        break;
      case '\\':
        strcat(f, "\\\\");
        break;
      default:
        strncat(f, p, 1);
    }
  }
  /* strcat(f, "\""); */
} /* mcparmprinter_string */

/* now we may define the parameter structure, using previous functions */
static struct
  {
    int (*getparm)(char *, void *);
    char * (*parminfo)(char *);
    void (*error)(char *, char *);
    void (*printer)(char *, void *);
} mcinputtypes[] = {
  {
    mcparm_double, mcparminfo_double, mcparmerror_double,
    mcparmprinter_double
  }, {
    mcparm_int, mcparminfo_int, mcparmerror_int,
    mcparmprinter_int
  }, {
    mcparm_string, mcparminfo_string, mcparmerror_string,
    mcparmprinter_string
  }
};

/*******************************************************************************
* mcestimate_error: compute sigma from N,p,p2 in Gaussian large numbers approx
*******************************************************************************/
double mcestimate_error(double N, double p1, double p2)
{
  double pmean, n1;
  if(N <= 1)
    return p1;
  pmean = p1 / N;
  n1 = N - 1;
  /* Note: underflow may cause p2 to become zero; the fabs() below guards
     against this. */
  return sqrt((N/n1)*fabs(p2 - pmean*pmean));
}

double (*mcestimate_error_p)
  (double V2, double psum, double p2sum)=mcestimate_error;

/* ========================================================================== */

/*                               MCCODE_R_IO_C                                */

/* ========================================================================== */

#ifndef MCCODE_R_IO_C
#define MCCODE_R_IO_C "$Revision$"

/* SECTION: file i/o handling ================================================ */

#ifndef HAVE_STRCASESTR
// from msysgit: https://code.google.com/p/msysgit/source/browse/compat/strcasestr.c
char *strcasestr(const char *haystack, const char *needle)
{
  int nlen = strlen(needle);
  int hlen = strlen(haystack) - nlen + 1;
  int i;

  for (i = 0; i < hlen; i++) {
    int j;
    for (j = 0; j < nlen; j++) {
            unsigned char c1 = haystack[i+j];
            unsigned char c2 = needle[j];
            if (toupper(c1) != toupper(c2))
                    goto next;
    }
    return (char *) haystack + i;
  next:
    ;
  }
  return NULL;
}


#endif
#ifndef HAVE_STRCASECMP
int strcasecmp( const char *s1, const char *s2 )
{
  int c1, c2;
  do {
    c1 = tolower( (unsigned char) *s1++ );
    c2 = tolower( (unsigned char) *s2++ );
  } while (c1 == c2 && c1 != 0);
  return c2 > c1 ? -1 : c1 > c2;
}
#endif

/*******************************************************************************
* mcfull_file: allocates a full file name=mcdirname+file. Catenate extension if missing.
*******************************************************************************/
char *mcfull_file(char *name, char *ext)
{
  int   dirlen=0;
  char *mem   =NULL;

  dirlen = mcdirname ? strlen(mcdirname) : 0;
  mem = (char*)malloc(dirlen + strlen(name) + CHAR_BUF_LENGTH);
  if(!mem) {
    exit(-fprintf(stderr, "Error: Out of memory %li (mcfull_file)\n", (long)(dirlen + strlen(name) + 256)));
  }
  strcpy(mem, "");

  /* prepend directory name to path if name does not contain a path */
  if (dirlen > 0 && !strchr(name, MC_PATHSEP_C)) {
    strcat(mem, mcdirname);
    strcat(mem, MC_PATHSEP_S);
  } /* dirlen */

  strcat(mem, name);
  if (!strchr(name, '.') && ext && strlen(ext))
  { /* add extension if not in file name already */
    strcat(mem, ".");
    strcat(mem, ext);
  }
  return(mem);
} /* mcfull_file */

/*******************************************************************************
* mcnew_file: opens a new file within mcdirname if non NULL
*             the file is opened in "a" (append, create if does not exist)
*             the extension 'ext' is added if the file name does not include one.
*             the last argument is set to 0 if file did not exist, else to 1.
*******************************************************************************/
FILE *mcnew_file(char *name, char *ext, int *exists)
{
  char *mem;
  FILE *file=NULL;

  if (!name || strlen(name) == 0 || mcdisable_output_files) return(NULL);
  
  mem  = mcfull_file(name, ext); /* create mcdirname/name.ext */
  
  /* check for existence */
  file = fopen(mem, "r"); /* for reading -> fails if does not exist */
  if (file) {
    fclose(file);
    *exists=1;
  } else
    *exists=0;
  
  /* open the file for writing/appending */
#ifdef USE_NEXUS
  if (mcformat && strcasestr(mcformat, "NeXus")) {
    /* NXhandle nxhandle is defined in the .h with USE_NEXUS */
    NXaccess mode = (*exists ? NXACC_CREATE5 | NXACC_RDWR : NXACC_CREATE5);
      
    if (NXopen(mem, mode, &nxhandle) != NX_OK)
      file = NULL;
    else
      file = (FILE*)&nxhandle; /* to make it non NULL */
  } else
#endif
    file = fopen(mem, "a+"); 
    
  if(!file)
    fprintf(stderr, "Warning: could not open output file '%s' for %s (mcnew_file)\n", 
      mem, *exists ? "append" : "create");
  free(mem);

  return file;
} /* mcnew_file */

/*******************************************************************************
* mcdetector_statistics: compute detector statistics, error bars, [x I I_err N] 1D
* RETURN:            updated detector structure
* Used by: mcdetector_import
*******************************************************************************/
MCDETECTOR mcdetector_statistics(
  MCDETECTOR detector)
{

  if (!detector.p1 || !detector.m)
    return(detector);
  
  /* compute statistics and update MCDETECTOR structure ===================== */
  double sum_z  = 0, min_z  = 0, max_z  = 0;
  double fmon_x =0,  smon_x = 0, fmon_y =0, smon_y=0, mean_z=0;
  double Nsum=0, P2sum=0;

  double sum_xz = 0, sum_yz = 0, sum_x = 0, sum_y = 0, sum_x2z = 0, sum_y2z = 0;
  int    i,j;
  char   hasnan=0, hasinf=0;
  char   israw = ((char*)strcasestr(detector.format,"raw") != NULL);
  double *this_p1=NULL; /* new 1D McCode array [x I E N]. Freed after writing data */

  /* if McCode/PGPLOT and rank==1 we create a new m*4 data block=[x I E N] */
  if (detector.rank == 1 && strcasestr(detector.format,"McCode")) {
    this_p1 = (double *)calloc(detector.m*detector.n*detector.p*4, sizeof(double));
    if (!this_p1)
      exit(-fprintf(stderr, "Error: Out of memory creating %li 1D " MCCODE_STRING " data set for file '%s' (mcdetector_import)\n",
        detector.m*detector.n*detector.p*4*sizeof(double*), detector.filename));
  }

  max_z = min_z = detector.p1[0];
  
  /* compute sum and moments (not for lists) */
  if (!strcasestr(detector.format,"list") && detector.m)
  for(j = 0; j < detector.n*detector.p; j++)
  {
    for(i = 0; i < detector.m; i++)
    {
      double x,y,z;
      double N, E;
      long   index= !detector.istransposed ? i*detector.n*detector.p + j : i+j*detector.m;
      char   hasnaninf=0;

      if (detector.m) 
        x = detector.xmin + (i + 0.5)/detector.m*(detector.xmax - detector.xmin); 
      else x = 0;
      if (detector.n && detector.p) 
        y = detector.ymin + (j + 0.5)/detector.n/detector.p*(detector.ymax - detector.ymin); 
      else y = 0;
      z = detector.p1[index];
      N = detector.p0 ? detector.p0[index] : 1;
      E = detector.p2 ? detector.p2[index] : 0;
      if (detector.p2 && !israw) 
        detector.p2[index] = (*mcestimate_error_p)(detector.p0[index],detector.p1[index],detector.p2[index]); /* set sigma */
      
      if (detector.rank == 1 && this_p1 && strcasestr(detector.format,"McCode")) {
        /* fill-in 1D McCode array [x I E N] */
        this_p1[index*4]   = x;
        this_p1[index*4+1] = z;
        this_p1[index*4+2] = detector.p2 ? detector.p2[index] : 0;
        this_p1[index*4+3] = N;
      }
      
      if (isnan(z) || isnan(E) || isnan(N)) hasnaninf=hasnan=1;
      if (isinf(z) || isinf(E) || isinf(N)) hasnaninf=hasinf=1;

      /* compute stats integrals */
      if (!hasnaninf) {
        sum_xz += x*z;
        sum_yz += y*z;
        sum_x  += x;
        sum_y  += y;
        sum_z  += z;
        sum_x2z += x*x*z;
        sum_y2z += y*y*z;
        if (z > max_z) max_z = z;
        if (z < min_z) min_z = z;

        Nsum += N;
        P2sum += E;
      }

    }
  } /* for j */

  /* compute 1st and 2nd moments. For lists, sum_z=0 so this is skipped. */
  if (sum_z && detector.n*detector.m*detector.p)
  {
    fmon_x = sum_xz/sum_z;
    fmon_y = sum_yz/sum_z;
    smon_x = sum_x2z/sum_z-fmon_x*fmon_x; smon_x = smon_x > 0 ? sqrt(smon_x) : 0;
    smon_y = sum_y2z/sum_z-fmon_y*fmon_y; smon_y = smon_y > 0 ? sqrt(smon_y) : 0;
    mean_z = sum_z/detector.n/detector.m/detector.p;
  }
  /* store statistics into detector */
  detector.intensity = sum_z;
  detector.error     = Nsum ? (*mcestimate_error_p)(Nsum, sum_z, P2sum) : 0;
  detector.events    = Nsum;
  detector.min       = min_z;
  detector.max       = max_z;
  detector.mean      = mean_z;
  detector.centerX   = fmon_x;
  detector.halfwidthX= smon_x;
  detector.centerY   = fmon_y;
  detector.halfwidthY= smon_y;

  /* if McCode/PGPLOT and rank==1 replace p1 with new m*4 1D McCode and clear others */
  if (detector.rank == 1 && this_p1 && strcasestr(detector.format,"McCode")) {
    
    detector.p1 = this_p1;
    detector.n  = detector.m; detector.m  = 4;
    detector.p0 = detector.p2 = NULL;
    detector.istransposed = 1;
  }

  if (detector.n*detector.m*detector.p > 1)
    snprintf(detector.signal, CHAR_BUF_LENGTH, 
      "Min=%g; Max=%g; Mean=%g;", detector.min, detector.max, detector.mean);
  else
    strcpy(detector.signal, "None");
  snprintf(detector.values, CHAR_BUF_LENGTH,
    "%g %g %g", detector.intensity, detector.error, detector.events);

  switch (detector.rank) {
    case 1:  snprintf(detector.statistics, CHAR_BUF_LENGTH, "X0=%g; dX=%g;",
      detector.centerX, detector.halfwidthX); break;
    case 2:
    case 3:  snprintf(detector.statistics, CHAR_BUF_LENGTH, "X0=%g; dX=%g; Y0=%g; dY=%g;",
      detector.centerX, detector.halfwidthX, detector.centerY, detector.halfwidthY);
      break;
    default: strcpy(detector.statistics, "None");
  }
  
  if (hasnan)
    printf("WARNING: Nan detected in component/file %s %s\n", 
      detector.component, strlen(detector.filename) ? detector.filename : "");
  if (hasinf)
    printf("WARNING: Inf detected in component/file %s %s\n", 
      detector.component, strlen(detector.filename) ? detector.filename : "");
  
  return(detector);
  
} /* mcdetector_statistics */

/*******************************************************************************
* mcdetector_import: build detector structure, merge non-lists from MPI
*                    compute basic stat, write "Detector:" line
* RETURN:            detector structure. Invalid data if detector.p1 == NULL
*                    Invalid detector sets m=0 and filename=""
*                    Simulation data  sets m=0 and filename=mcsiminfo_name
* This function is equivalent to the old 'mcdetector_out', returning a structure
*******************************************************************************/
MCDETECTOR mcdetector_import(
  char *format,
  char *component, char *title,
  long m, long n,  long p,
  char *xlabel, char *ylabel, char *zlabel,
  char *xvar, char *yvar, char *zvar,
  double x1, double x2, double y1, double y2, double z1, double z2,
  char *filename,
  double *p0, double *p1, double *p2,
  Coords position)
{
  time_t t;       /* for detector.date */
  long   date_l;  /* date as a long number */
  char   istransposed=0;
  char   c[CHAR_BUF_LENGTH]; /* temp var for signal label */

  MCDETECTOR detector;

  /* build MCDETECTOR structure ============================================= */
  /* make sure we do not have NULL for char fields */

  /* these also apply to simfile */
  strncpy (detector.filename,  filename ? filename : "",        CHAR_BUF_LENGTH);
  strncpy (detector.format,    format   ? format   : "McCode" , CHAR_BUF_LENGTH);
  /* add extension if missing */
  if (strlen(detector.filename) && !strchr(detector.filename, '.'))
  { /* add extension if not in file name already */
    strcat(detector.filename, ".dat");
  }
  strncpy (detector.component, component ? component : MCCODE_STRING " component", CHAR_BUF_LENGTH);

  snprintf(detector.instrument, CHAR_BUF_LENGTH, "%s (%s)", mcinstrument_name, mcinstrument_source);
  snprintf(detector.user, CHAR_BUF_LENGTH,      "%s on %s",
        getenv("USER") ? getenv("USER") : MCCODE_NAME,
        getenv("HOST") ? getenv("HOST") : "localhost");
  time(&t);         /* get current write time */
  date_l = (long)t; /* same but as a long */
  snprintf(detector.date, CHAR_BUF_LENGTH, "%s", ctime(&t));
  if (strlen(detector.date))   detector.date[strlen(detector.date)-1] = '\0'; /* remove last \n in date */
  detector.date_l = date_l;

  if (!mcget_run_num() || mcget_run_num() >= mcget_ncount())
    snprintf(detector.ncount, CHAR_BUF_LENGTH, "%llu", mcget_ncount()
#ifdef USE_MPI
*mpi_node_count
#endif
  );
  else
    snprintf(detector.ncount, CHAR_BUF_LENGTH, "%g/%g", (double)mcget_run_num(), (double)mcget_ncount());

  detector.p0         = p0;
  detector.p1         = p1;
  detector.p2         = p2;

  /* handle transposition (not for NeXus) */
  if (!strcasestr(detector.format, "NeXus")) {
    if (m<0 || n<0 || p<0)             istransposed = !istransposed;
    if (strcasestr(detector.format, "transpose")) istransposed = !istransposed;
    if (istransposed) { /* do the swap once for all */
      long i=m; m=n; n=i;
    }
  }

  m=labs(m); n=labs(n); p=labs(p); /* make sure dimensions are positive */
  detector.istransposed = istransposed;

  /* determine detector rank (dimensionality) */
  if (!m || !n || !p || !p1) detector.rank = 4; /* invalid: exit with m=0 filename="" */
  else if (m*n*p == 1)       detector.rank = 0; /* 0D */
  else if (n == 1 || m == 1) detector.rank = 1; /* 1D */
  else if (p == 1)           detector.rank = 2; /* 2D */
  else                       detector.rank = 3; /* 3D */

  /* from rank, set type */
  switch (detector.rank) {
    case 0:  strcpy(detector.type,  "array_0d"); m=n=p=1; break;
    case 1:  snprintf(detector.type, CHAR_BUF_LENGTH, "array_1d(%ld)", m*n*p); m *= n*p; n=p=1; break;
    case 2:  snprintf(detector.type, CHAR_BUF_LENGTH, "array_2d(%ld, %ld)", m, n*p); n *= p; p=1; break;
    case 3:  snprintf(detector.type, CHAR_BUF_LENGTH, "array_3d(%ld, %ld, %ld)", m, n, p); break;
    default: m=0; strcpy(detector.type, ""); strcpy(detector.filename, "");/* invalid */
  }

  detector.m    = m;
  detector.n    = n;
  detector.p    = p;

  /* these only apply to detector files ===================================== */

  snprintf(detector.position, CHAR_BUF_LENGTH, "%g %g %g", position.x, position.y, position.z);
  /* may also store actual detector orientation in the future */

  strncpy(detector.title,      title && strlen(title) ? title : component,       CHAR_BUF_LENGTH);
  strncpy(detector.xlabel,     xlabel && strlen(xlabel) ? xlabel : "X", CHAR_BUF_LENGTH); /* axis labels */
  strncpy(detector.ylabel,     ylabel && strlen(ylabel) ? ylabel : "Y", CHAR_BUF_LENGTH);
  strncpy(detector.zlabel,     zlabel && strlen(zlabel) ? zlabel : "Z", CHAR_BUF_LENGTH);
  strncpy(detector.xvar,       xvar && strlen(xvar) ? xvar :       "x", CHAR_BUF_LENGTH); /* axis variables */
  strncpy(detector.yvar,       yvar && strlen(yvar) ? yvar :       detector.xvar, CHAR_BUF_LENGTH);
  strncpy(detector.zvar,       zvar && strlen(zvar) ? zvar :       detector.yvar, CHAR_BUF_LENGTH);

  /* set "variables" as e.g. "I I_err N" */
  strcpy(c, "I ");
  if (strlen(detector.zvar))      strncpy(c, detector.zvar,32);
  else if (strlen(detector.yvar)) strncpy(c, detector.yvar,32);
  else if (strlen(detector.xvar)) strncpy(c, detector.xvar,32);

  if (detector.rank == 1)
    snprintf(detector.variables, CHAR_BUF_LENGTH, "%s %s %s_err N", detector.xvar, c, c);
  else
    snprintf(detector.variables, CHAR_BUF_LENGTH, "%s %s_err N", c, c);

  /* limits */
  detector.xmin = x1;
  detector.xmax = x2;
  detector.ymin = y1;
  detector.ymax = y2;
  detector.zmin = z1;
  detector.zmax = z2;
  if (abs(detector.rank) == 1)
    snprintf(detector.limits, CHAR_BUF_LENGTH, "%g %g", x1, x2);
  else if (detector.rank == 2)
    snprintf(detector.limits, CHAR_BUF_LENGTH, "%g %g %g %g", x1, x2, y1, y2);
  else
    snprintf(detector.limits, CHAR_BUF_LENGTH, "%g %g %g %g %g %g", x1, x2, y1, y2, z1, z2);

  /* if MPI and nodes_nb > 1: reduce data sets when using MPI =============== */
#ifdef USE_MPI
  if (!strcasestr(detector.format,"list") && mpi_node_count > 1 && m) {
    /* we save additive data: reduce everything into mpi_node_root */
    if (p0) mc_MPI_Sum(p0, m*n*p);
    if (p1) mc_MPI_Sum(p1, m*n*p);
    if (p2) mc_MPI_Sum(p2, m*n*p);
    if (!p0) {  /* additive signal must be then divided by the number of nodes */
      int i;
      for (i=0; i<m*n*p; i++) {
        p1[i] /= mpi_node_count;
        if (p2) p2[i] /= mpi_node_count;
      }
    }
  }
#endif /* USE_MPI */

  /* compute statistics, Nsum, intensity, Error bars */
  detector = mcdetector_statistics(detector);

#ifdef USE_MPI
  /* slaves are done */
  if(mpi_node_rank != mpi_node_root) {
    return detector;
  }
#endif

  /* output "Detector:" line ================================================ */
  /* when this is a detector written by a component (not the SAVE from instrument),
     not an event lists */
  if (!m) return(detector);
  if (!strcasestr(detector.format,"list")) {
    if (!strcmp(detector.component, mcinstrument_name)) {
      if (strlen(detector.filename))  /* we name it from its filename, or from its title */
        strncpy(c, detector.filename, CHAR_BUF_LENGTH);
      else
        snprintf(c, CHAR_BUF_LENGTH, "%s", mcinstrument_name);
    } else
      strncpy(c, detector.component, CHAR_BUF_LENGTH);  /* usual detectors written by components */

    printf("Detector: %s_I=%g %s_ERR=%g %s_N=%g",
           c, detector.intensity,
           c, detector.error,
           c, detector.events);
    printf(" \"%s\"\n", strlen(detector.filename) ? detector.filename : detector.component);
  }
  

  return(detector);
} /* mcdetector_import */

/* end MCDETECTOR import section ============================================ */

















/* ========================================================================== */

/*                               ASCII output                                 */
/*     The SIM file is YAML based, the data files have '#' headers            */

/* ========================================================================== */


/*******************************************************************************
* mcinfo_out: output instrument tags/info (only in SIM)
* Used in: mcsiminfo_init (ascii), mcinfo(stdout)
*******************************************************************************/
static void mcinfo_out(char *pre, FILE *f)
{
  char Parameters[CHAR_BUF_LENGTH] = "";
  int  i;

  if (!f || mcdisable_output_files) return;

  /* create parameter string ================================================ */
  for(i = 0; i < mcnumipar; i++)
  {
    char ThisParam[CHAR_BUF_LENGTH];
    if (strlen(mcinputtable[i].name) > CHAR_BUF_LENGTH) break;
    snprintf(ThisParam, CHAR_BUF_LENGTH, " %s(%s)", mcinputtable[i].name,
            (*mcinputtypes[mcinputtable[i].type].parminfo)
                (mcinputtable[i].name));
    strcat(Parameters, ThisParam);
    if (strlen(Parameters) >= CHAR_BUF_LENGTH-64) break;
  }

  /* output data ============================================================ */
  if (f != stdout)
    fprintf(f, "%sFile: %s%c%s\n",    pre, mcdirname, MC_PATHSEP_C, mcsiminfo_name);
  else
    fprintf(f, "%sCreator: %s\n",     pre, MCCODE_STRING);

  fprintf(f, "%sSource: %s\n",   pre, mcinstrument_source);
  fprintf(f, "%sParameters: %s\n",    pre, Parameters);
  
  fprintf(f, "%sTrace_enabled: %s\n", pre, mctraceenabled ? "yes" : "no");
  fprintf(f, "%sDefault_main: %s\n",  pre, mcdefaultmain ?  "yes" : "no");
  fprintf(f, "%sEmbedded_runtime: %s\n", pre, 
#ifdef MC_EMBEDDED_RUNTIME
         "yes"
#else
         "no"
#endif
         );

  fflush(f);
} /* mcinfo_out */

/*******************************************************************************
* mcruninfo_out_backend: output simulation tags/info (both in SIM and data files)
* Used in: mcsiminfo_init (ascii case), mcdetector_out_xD_ascii, mcinfo(stdout)
*******************************************************************************/
static void mcruninfo_out_backend(char *pre, FILE *f, int info)
{
  int i;
  char Parameters[CHAR_BUF_LENGTH];

  if (!f || mcdisable_output_files) return;

  fprintf(f, "%sFormat: %s%s\n",      pre, 
    mcformat && strlen(mcformat) ? mcformat : MCCODE_NAME,
    mcformat && strcasestr(mcformat,"McCode") ? " with text headers" : "");
  fprintf(f, "%sURL: %s\n",         pre, "http://www.mccode.org");
  fprintf(f, "%sCreator: %s\n",     pre, MCCODE_STRING);
  fprintf(f, "%sInstrument: %s\n", pre, mcinstrument_source);
  fprintf(f, "%sNcount: %llu\n",        pre, mcget_ncount());
  fprintf(f, "%sTrace: %s\n",       pre, mcdotrace ? "yes" : "no");
  fprintf(f, "%sGravitation: %s\n", pre, mcgravitation ? "yes" : "no");
  snprintf(Parameters, CHAR_BUF_LENGTH, "%ld", mcseed);
  fprintf(f, "%sSeed: %s\n",        pre, Parameters);
  fprintf(f, "%sDirectory: %s\n",        pre, mcdirname ? mcdirname : ".");
#ifdef USE_MPI
  if (mpi_node_count > 1)
    fprintf(f, "%sNodes: %i\n",        pre, mpi_node_count);
#endif

  /* output parameter string ================================================ */
  for(i = 0; i < mcnumipar; i++) {
      if (!info){
          (*mcinputtypes[mcinputtable[i].type].printer)(Parameters, mcinputtable[i].par);
          fprintf(f, "%sParam: %s=%s\n", pre, mcinputtable[i].name, Parameters);
      }else{
        /*if an info run, some variables might not have values. Flag these by "NULL"*/
	if(mcinputtable[i].val && strlen(mcinputtable[i].val)){
            /* ... those with defautl values*/
            (*mcinputtypes[mcinputtable[i].type].printer)(Parameters, mcinputtable[i].par);
            fprintf(f, "%sParam: %s=%s\n", pre, mcinputtable[i].name, Parameters);
        }else{
            /* ... and those without */
            fprintf(f, "%sParam: %s=NULL\n", pre, mcinputtable[i].name);
	}
      }
  }
} /* mcruninfo_out_backend */

/************************
* wrapper function to mcruninfo_out_backend
*  Regular runs use this whereas the single call from mcinfo is directly to the backend
*************************/
static void mcruninfo_out(char *pre, FILE *f){
    mcruninfo_out_backend(pre,f,0);
}

/*******************************************************************************
* mcsiminfo_out:    wrapper to fprintf(mcsiminfo_file)
*******************************************************************************/
void mcsiminfo_out(char *format, ...)
{
  va_list ap;

  if(mcsiminfo_file && !mcdisable_output_files)
  {
    va_start(ap, format);
    vfprintf(mcsiminfo_file, format, ap);
    va_end(ap);
  }
} /* mcsiminfo_out */


/*******************************************************************************
* mcdatainfo_out: output detector header
*   mcdatainfo_out(prefix, file_handle, detector) writes info to data file
*******************************************************************************/
static void
mcdatainfo_out(char *pre, FILE *f, MCDETECTOR detector)
{
  if (!f || !detector.m || mcdisable_output_files) return;
  
  /* output data ============================================================ */
  fprintf(f, "%sDate: %s (%li)\n",       pre, detector.date, detector.date_l);
  fprintf(f, "%stype: %s\n",       pre, detector.type);
  fprintf(f, "%sSource: %s\n",     pre, detector.instrument);
  fprintf(f, "%scomponent: %s\n",  pre, detector.component);
  fprintf(f, "%sposition: %s\n",   pre, detector.position);

  fprintf(f, "%stitle: %s\n",      pre, detector.title);
  fprintf(f, !mcget_run_num() || mcget_run_num() >= mcget_ncount() ?
             "%sNcount: %s\n" : 
             "%sratio: %s\n",  pre, detector.ncount);

  if (strlen(detector.filename)) {
    fprintf(f, "%sfilename: %s\n", pre, detector.filename);
  }

  fprintf(f, "%sstatistics: %s\n", pre, detector.statistics);
  fprintf(f, "%ssignal: %s\n",     pre, detector.signal);
  fprintf(f, "%svalues: %s\n",     pre, detector.values);

  if (detector.rank >= 1)
  {
    fprintf(f, "%sxvar: %s\n",     pre, detector.xvar);
    fprintf(f, "%syvar: %s\n",     pre, detector.yvar);
    fprintf(f, "%sxlabel: %s\n",   pre, detector.xlabel);
    fprintf(f, "%sylabel: %s\n",   pre, detector.ylabel);
    if (detector.rank > 1) {
      fprintf(f, "%szvar: %s\n",   pre, detector.zvar);
      fprintf(f, "%szlabel: %s\n", pre, detector.zlabel);
    }
  }

  fprintf(f, 
    abs(detector.rank)==1 ?
             "%sxlimits: %s\n" : 
             "%sxylimits: %s\n", pre, detector.limits);
  fprintf(f, "%svariables: %s\n", pre, 
    strcasestr(detector.format, "list") ? detector.ylabel : detector.variables);
    
  fflush(f);

} /* mcdatainfo_out */

/* mcdetector_out_array_ascii: output a single array to a file
 *   m: columns
 *   n: rows
 *   p: array
 *   f: file handle (already opened)
 */
static void mcdetector_out_array_ascii(long m, long n, double *p, FILE *f, char istransposed)
{
  if(f)
  {
    int i,j;
    for(j = 0; j < n; j++)
    {
      for(i = 0; i < m; i++)
      {
          fprintf(f, "%.10g ", p[!istransposed ? i*n + j : j*m+i]);
      }
      fprintf(f,"\n");
    }
  }
} /* mcdetector_out_array_ascii */

/*******************************************************************************
* mcdetector_out_0D_ascii: called by mcdetector_out_0D for ascii output
*******************************************************************************/
MCDETECTOR mcdetector_out_0D_ascii(MCDETECTOR detector)
{
  int exists=0;
  FILE *outfile = NULL;
  
  /* Write data set information to simulation description file. */
  MPI_MASTER(
    mcsiminfo_out("\nbegin data\n"); // detector.component
    mcdatainfo_out("  ", mcsiminfo_file, detector);
    mcsiminfo_out("end data\n");
    /* Don't write if filename is NULL: mcnew_file handles this (return NULL) */
    outfile = mcnew_file(detector.component, "dat", &exists);
    if(outfile)
    {
      /* write data file header and entry in simulation description file */
      mcruninfo_out( "# ", outfile);
      mcdatainfo_out("# ", outfile, detector);
      /* write I I_err N */
      fprintf(outfile, "%g %g %g\n", 
        detector.intensity, detector.error, detector.events);
      fclose(outfile);
    }
  ); /* MPI_MASTER */
  return(detector);
} /* mcdetector_out_0D_ascii */

/*******************************************************************************
* mcdetector_out_1D_ascii: called by mcdetector_out_1D for ascii output
*******************************************************************************/
MCDETECTOR mcdetector_out_1D_ascii(MCDETECTOR detector)
{
  int exists=0;
  FILE *outfile = NULL;

  MPI_MASTER(
    /* Write data set information to simulation description file. */
    mcsiminfo_out("\nbegin data\n"); // detector.filename
    mcdatainfo_out("  ", mcsiminfo_file, detector);
    mcsiminfo_out("end data\n");
    /* Loop over array elements, writing to file. */
    /* Don't write if filename is NULL: mcnew_file handles this (return NULL) */
    outfile = mcnew_file(detector.filename, "dat", &exists);
    if(outfile)
    {
      /* write data file header and entry in simulation description file */
      mcruninfo_out( "# ", outfile);
      mcdatainfo_out("# ", outfile, detector);
      /* output the 1D array columns */
      mcdetector_out_array_ascii(detector.m, detector.n, detector.p1, outfile, detector.istransposed);
      
      fclose(outfile);
    }
  ); /* MPI_MASTER */
  return(detector);
  
}  /* mcdetector_out_1D_ascii */

/*******************************************************************************
* mcdetector_out_2D_ascii: called by mcdetector_out_2D for ascii output
*******************************************************************************/
MCDETECTOR mcdetector_out_2D_ascii(MCDETECTOR detector)
{
  int exists=0;
  FILE *outfile = NULL;
  
  MPI_MASTER(
    /* Loop over array elements, writing to file. */
    /* Don't write if filename is NULL: mcnew_file handles this (return NULL) */
    outfile = mcnew_file(detector.filename, "dat", &exists);
    if(outfile)
    {
      /* write header only if file has just been created (not appending) */
      if (!exists) {
        /* Write data set information to simulation description file. */
        mcsiminfo_out("\nbegin data\n"); // detector.filename
        mcdatainfo_out("  ", mcsiminfo_file, detector);
        mcsiminfo_out("end data\n");
      
        mcruninfo_out( "# ", outfile);
        mcdatainfo_out("# ", outfile,   detector);
      }
      fprintf(outfile, "# Data [%s/%s] %s:\n", detector.component, detector.filename, detector.zvar);
      mcdetector_out_array_ascii(detector.m, detector.n*detector.p, detector.p1, 
        outfile, detector.istransposed);
      if (detector.p2) {
        fprintf(outfile, "# Errors [%s/%s] %s_err:\n", detector.component, detector.filename, detector.zvar);
        mcdetector_out_array_ascii(detector.m, detector.n*detector.p, detector.p2, 
          outfile, detector.istransposed);
      }
      if (detector.p0) {
        fprintf(outfile, "# Events [%s/%s] N:\n", detector.component, detector.filename);
        mcdetector_out_array_ascii(detector.m, detector.n*detector.p, detector.p0, 
          outfile, detector.istransposed);
      }
      fclose(outfile);
      
      if (!exists) {
        if (strcasestr(detector.format, "list"))
          printf("Events:   \"%s\"\n",  
            strlen(detector.filename) ? detector.filename : detector.component);
      }
    } /* if outfile */
  ); /* MPI_MASTER */
#ifdef USE_MPI
  if (strcasestr(detector.format, "list") && mpi_node_count > 1) {
    int node_i=0;
    /* loop along MPI nodes to write sequentially */
    for(node_i=0; node_i<mpi_node_count; node_i++) {
      /* MPI: slaves wait for the master to write its block, then append theirs */
      MPI_Barrier(MPI_COMM_WORLD);
      if (node_i != mpi_node_root && node_i == mpi_node_rank) {
        if(strlen(detector.filename) && !mcdisable_output_files)	/* Don't write if filename is NULL */
          outfile = mcnew_file(detector.filename, "dat", &exists);
        if (!exists)
          fprintf(stderr, "Warning: [MPI node %i] file '%s' does not exist yet, "
                          "MASTER should have opened it before.\n",
            mpi_node_rank, detector.filename);
        if(outfile) {
          mcdetector_out_array_ascii(detector.m, detector.n*detector.p, detector.p1, 
            outfile, detector.istransposed);
          fclose(outfile);
        }
      }
    }
  } /* if strcasestr list */
#endif
  return(detector);
} /* mcdetector_out_2D_ascii */

/*******************************************************************************
* strcpy_valid: makes a valid string for variable names.
*   copy 'original' into 'valid', replacing invalid characters by '_'
*   char arrays must be pre-allocated
*******************************************************************************/
static char *strcpy_valid(char *valid, char *original)
{
  long i;
  int  n=32; /* max length of valid names */

  if (original == NULL || !strlen(original)) return(NULL);

  if (n > strlen(original)) n = strlen(original);
  else original += strlen(original)-n;
  strncpy(valid, original, n);

  for (i=0; i < n; i++)
  {
    if ( (valid[i] > 122)
      || (valid[i] < 32)
      || (strchr("!\"#$%&'()*+,-.:;<=>?@[\\]^`/ \n\r\t", valid[i]) != NULL) )
    {
      if (i) valid[i] = '_'; else valid[i] = 'm';
    }
  }
  valid[i] = '\0';

  return(valid);
} /* strcpy_valid */

/* end ascii output section ================================================= */







#ifdef USE_NEXUS

/* ========================================================================== */

/*                               NeXus output                                 */

/* ========================================================================== */

#define nxprintf(...)    nxstr('d', __VA_ARGS__)
#define nxprintattr(...) nxstr('a', __VA_ARGS__)

/*******************************************************************************
* nxstr: output a tag=value data set (char) in NeXus/current group
*   when 'format' is larger that 1024 chars it is used as value for the 'tag'
*   else the value is assembled with format and following arguments.
*   type='d' -> data set
*        'a' -> attribute for current data set
*******************************************************************************/
static int nxstr(char type, NXhandle *f, char *tag, char *format, ...)
{
  va_list ap;
  char value[CHAR_BUF_LENGTH];
  int  i;
  int  ret=NX_OK;
  
  if (!tag || !format || !strlen(tag) || !strlen(format)) return(NX_OK);
  
  /* assemble the value string */
  if (strlen(format) < CHAR_BUF_LENGTH) {
    va_start(ap, format);
    ret = vsnprintf(value, CHAR_BUF_LENGTH, format, ap);
    va_end(ap);
  
    i = strlen(value);
  } else {
    i = strlen(format);
  }

  if (type == 'd') {
    /* open/put/close data set */
    if (NXmakedata (f, tag, NX_CHAR, 1, &i) != NX_OK) return(NX_ERROR);
    NXopendata (f, tag);
    if (strlen(format) < CHAR_BUF_LENGTH)
      ret = NXputdata  (f, value);
    else
      ret = NXputdata  (f, format);
    NXclosedata(f);
  } else {
    if (strlen(format) < CHAR_BUF_LENGTH)
      ret = NXputattr  (f, tag, value, strlen(value), NX_CHAR);
    else
      ret = NXputattr  (f, tag, format, strlen(format), NX_CHAR);
  }
  
  return(ret);
  
} /* nxstr */

/*******************************************************************************
* mcinfo_readfile: read a full file into a string buffer which is allocated
*   Think to free the buffer after use.
* Used in: mcinfo_out_nexus (nexus)
*******************************************************************************/
char *mcinfo_readfile(char *filename)
{
  FILE *f = fopen(filename, "rb");
  if (!f) return(NULL);
  fseek(f, 0, SEEK_END);
  long fsize = ftell(f);
  rewind(f);
  char *string = malloc(fsize + 1);
  if (string) {
    int n = fread(string, fsize, 1, f);
    fclose(f);

    string[fsize] = 0;
  }
  return(string);
}

/*******************************************************************************
* mcinfo_out: output instrument/simulation groups in NeXus file
* Used in: mcsiminfo_init (nexus)
*******************************************************************************/
static void mcinfo_out_nexus(NXhandle f)
{
  FILE  *fid;     /* for intrument source code/C/IDF */
  char  *buffer=NULL;
  time_t t     =time(NULL); /* for date */
  char   entry0[CHAR_BUF_LENGTH];
  int    count=0;
  char   name[CHAR_BUF_LENGTH];
  char   class[CHAR_BUF_LENGTH];
  
  if (!f || mcdisable_output_files) return;
  
  /* write NeXus NXroot attributes */
  /* automatically added: file_name, HDF5_Version, file_time, NeXus_version */ 
  nxprintattr(f, "creator",   "%s generated with " MCCODE_STRING, mcinstrument_name);
  
  /* count the number of existing NXentry and create the next one */
  NXgetgroupinfo(f, &count, name, class);
  sprintf(entry0, "entry%i", count+1);

  /* create the main NXentry (mandatory in NeXus) */
  if (NXmakegroup(f, entry0, "NXentry") == NX_OK) 
  if (NXopengroup(f, entry0, "NXentry") == NX_OK) {
    
    nxprintf(nxhandle, "program_name", MCCODE_STRING);
    nxprintf(f, "start_time", ctime(&t));
    nxprintf(f, "title", "%s%s%s simulation generated by instrument %s", 
      mcdirname && strlen(mcdirname) ? mcdirname : ".", MC_PATHSEP_S, mcsiminfo_name,
      mcinstrument_name);
    nxprintattr(f, "program_name", MCCODE_STRING);
    nxprintattr(f, "instrument",   mcinstrument_name);
    nxprintattr(f, "simulation",   "%s%s%s",
        mcdirname && strlen(mcdirname) ? mcdirname : ".", MC_PATHSEP_S, mcsiminfo_name);

    /* write NeXus instrument group */
    if (NXmakegroup(f, "instrument", "NXinstrument") == NX_OK)
    if (NXopengroup(f, "instrument", "NXinstrument") == NX_OK) {
      int   i;
      char *string=NULL;

      /* write NeXus parameters(types) data =================================== */
      string = (char*)malloc(CHAR_BUF_LENGTH);
      if (string) {
        strcpy(string, "");
        for(i = 0; i < mcnumipar; i++)
        {
          char ThisParam[CHAR_BUF_LENGTH];
          snprintf(ThisParam, CHAR_BUF_LENGTH, " %s(%s)", mcinputtable[i].name,
                  (*mcinputtypes[mcinputtable[i].type].parminfo)
                      (mcinputtable[i].name));
          if (strlen(string) + strlen(ThisParam) < CHAR_BUF_LENGTH)
            strcat(string, ThisParam);
        }
        nxprintattr(f, "Parameters",    string);
        free(string);
      }
        
      nxprintattr(f, "name",          mcinstrument_name);
      nxprintf   (f, "name",          mcinstrument_name);
      nxprintattr(f, "Source",        mcinstrument_source);
      
      nxprintattr(f, "Trace_enabled", mctraceenabled ? "yes" : "no");
      nxprintattr(f, "Default_main",  mcdefaultmain ?  "yes" : "no");
      nxprintattr(f, "Embedded_runtime",  
  #ifdef MC_EMBEDDED_RUNTIME
           "yes"
  #else
           "no"
  #endif
           );
           
      /* add instrument source code when available */
      buffer = mcinfo_readfile(mcinstrument_source);
      if (buffer && strlen(buffer)) {
        long length=strlen(buffer);
        nxprintf (f, "description", buffer);
        NXopendata(f,"description");
        nxprintattr(f, "file_name", mcinstrument_source);
        nxprintattr(f, "file_size", "%li", length);
        nxprintattr(f, "MCCODE_STRING", MCCODE_STRING);
        NXclosedata(f);
        nxprintf (f,"instrument_source", "%s " MCCODE_NAME " " MCCODE_PARTICLE " Monte Carlo simulation", mcinstrument_name);
        free(buffer);
      } else
        nxprintf (f, "description", "File %s not found (instrument description %s is missing)", 
          mcinstrument_source, mcinstrument_name);
      
      /* add Mantid/IDF.xml when available */
      char *IDFfile=NULL;
      IDFfile = (char*)malloc(CHAR_BUF_LENGTH);
      sprintf(IDFfile,"%s%s",mcinstrument_source,".xml");
      buffer = mcinfo_readfile(IDFfile);
      if (buffer && strlen(buffer)) {
        NXmakegroup (nxhandle, "instrument_xml", "NXnote");
        NXopengroup (nxhandle, "instrument_xml", "NXnote");
        nxprintf(f, "data", buffer);
        nxprintf(f, "description", "IDF.xml file found with instrument %s", mcinstrument_source);
        nxprintf(f, "type", "text/xml");
        NXclosegroup(f); /* instrument_xml */
        free(buffer);
      }
      free(IDFfile);
      NXclosegroup(f); /* instrument */
    } /* NXinstrument */

    /* write NeXus simulation group */
    if (NXmakegroup(f, "simulation", "NXnote") == NX_OK)
    if (NXopengroup(f, "simulation", "NXnote") == NX_OK) {

      nxprintattr(f, "name",   "%s%s%s",
        mcdirname && strlen(mcdirname) ? mcdirname : ".", MC_PATHSEP_S, mcsiminfo_name);
      
      nxprintf   (f, "name",      "%s",     mcsiminfo_name);
      nxprintattr(f, "Format",    mcformat && strlen(mcformat) ? mcformat : MCCODE_NAME);
      nxprintattr(f, "URL",       "http://www.mccode.org");
      nxprintattr(f, "program",   MCCODE_STRING);
      nxprintattr(f, "Instrument",mcinstrument_source);
      nxprintattr(f, "Trace",     mcdotrace ?     "yes" : "no");
      nxprintattr(f, "Gravitation",mcgravitation ? "yes" : "no");
      nxprintattr(f, "Seed",      "%li", mcseed);
      nxprintattr(f, "Directory", mcdirname);
    #ifdef USE_MPI
      if (mpi_node_count > 1)
        nxprintf(f, "Nodes", "%i",        mpi_node_count);
    #endif
    
      /* output parameter string ================================================ */
      if (NXmakegroup(f, "Param", "NXparameters") == NX_OK)
      if (NXopengroup(f, "Param", "NXparameters") == NX_OK) {
        int i;
        char string[CHAR_BUF_LENGTH];
        for(i = 0; i < mcnumipar; i++) {
          if (mcget_run_num() || (mcinputtable[i].val && strlen(mcinputtable[i].val))) {
            if (mcinputtable[i].par == NULL)
              strncpy(string, (mcinputtable[i].val ? mcinputtable[i].val : ""), CHAR_BUF_LENGTH);
            else
              (*mcinputtypes[mcinputtable[i].type].printer)(string, mcinputtable[i].par);

            nxprintf(f,  mcinputtable[i].name, "%s", string);
            nxprintattr(f, mcinputtable[i].name, string);
          }
        }
        NXclosegroup(f); /* Param */
      } /* NXparameters */
      
      NXclosegroup(f); /* simulation */
    } /* NXsimulation */
    
    /* create a group to hold all monitors */
    NXmakegroup(f, "data", "NXdetector");

    /* leave the NXentry opened (closed at exit) */
  } /* NXentry */
} /* mcinfo_out_nexus */

/*******************************************************************************
* mcdatainfo_out_nexus: output detector header
*   mcdatainfo_out_nexus(detector) create group and write info to NeXus data file
*   open data:NXdetector then filename:NXdata and write headers/attributes
*   requires: NXentry to be opened
*******************************************************************************/
static void
mcdatainfo_out_nexus(NXhandle f, MCDETECTOR detector)
{
  char data_name[32];
  if (!f || !detector.m || mcdisable_output_files) return;
  
  strcpy_valid(data_name, 
    detector.filename && strlen(detector.filename) ? 
      detector.filename : detector.component);

  /* the NXdetector group has been created in mcinfo_out_nexus (mcsiminfo_init) */
  if (NXopengroup(f, "data", "NXdetector") == NX_OK) {

    /* create and open the data group */
    /* this may fail when appending to list -> ignore/skip */
    NXMDisableErrorReporting(); /* unactivate NeXus error messages, as creation may fail */
    
    if (NXmakegroup(f, data_name, "NXdata") == NX_OK)
    if (NXopengroup(f, data_name, "NXdata") == NX_OK) {
    
      /* output metadata (as attributes) ======================================== */
      nxprintattr(f, "Date",       detector.date);
      nxprintattr(f, "type",       detector.type);
      nxprintattr(f, "Source",     detector.instrument);
      nxprintattr(f, "component",  detector.component);
      nxprintattr(f, "position",   detector.position);

      nxprintattr(f, "title",      detector.title);
      nxprintattr(f, !mcget_run_num() || mcget_run_num() >= mcget_ncount() ?
                 "Ncount" : 
                 "ratio",  detector.ncount);

      if (strlen(detector.filename)) {
        nxprintattr(f, "filename", detector.filename);
      }

      nxprintattr(f, "statistics", detector.statistics);
      nxprintattr(f, "signal",     detector.signal);
      nxprintattr(f, "values",     detector.values);

      if (detector.rank >= 1)
      {
        nxprintattr(f, "xvar",     detector.xvar);
        nxprintattr(f, "yvar",     detector.yvar);
        nxprintattr(f, "xlabel",   detector.xlabel);
        nxprintattr(f, "ylabel",   detector.ylabel);
        if (detector.rank > 1) {
          nxprintattr(f, "zvar",   detector.zvar);
          nxprintattr(f, "zlabel", detector.zlabel);
        }
      }

      nxprintattr(f, abs(detector.rank)==1 ?
                 "xlimits" : 
                 "xylimits", detector.limits);
      nxprintattr(f, "variables", 
        strcasestr(detector.format, "list") ? detector.ylabel : detector.variables);
      nxprintf(f, "distance", detector.position);
      nxprintf(f, "acquisition_mode",
        strcasestr(detector.format, "list") ? "event" : "summed");
        
      NXclosegroup(f);
    } /* NXdata (filename) */
    NXMEnableErrorReporting();  /* re-enable NeXus error messages */
    NXclosegroup(f);
  } /* NXdetector (data) */
  
} /* mcdatainfo_out_nexus */

/*******************************************************************************
* mcdetector_out_axis_nexus: write detector axis into current NXdata
*   requires: NXdata to be opened
*******************************************************************************/
int mcdetector_out_axis_nexus(NXhandle f, char *label, char *var, int rank, long length, double min, double max)
{
  if (!f || length <= 1 || mcdisable_output_files || max == min) return(NX_OK);
  else {
    double axis[length];
    char valid[32];
    int dim=(int)length;
    int i;
    int nprimary=1;
    /* create an axis from [min:max] */
    for(i = 0; i < length; i++)
      axis[i] = min+(max-min)*(i+0.5)/length;
    /* create the data set */
    strcpy_valid(valid, label);
    NXcompmakedata(f, valid, NX_FLOAT64, 1, &dim, NX_COMP_LZW, &dim);
    /* open it */
    if (NXopendata(f, valid) != NX_OK) {
      fprintf(stderr, "Warning: could not open axis rank %i '%s' (NeXus)\n",
        rank, valid);
      return(NX_ERROR);
    }
    /* put the axis and its attributes */
    NXputdata  (f, axis);
    nxprintattr(f, "long_name",  label);
    nxprintattr(f, "short_name", var);
    NXputattr  (f, "axis",       &rank,     1, NX_INT32);
    nxprintattr(f, "units",      var);
    NXputattr  (f, "primary",    &nprimary, 1, NX_INT32);
    NXclosedata(f);
    
    return(NX_OK);
  }
} /* mcdetector_out_axis_nexus */

/*******************************************************************************
* mcdetector_out_array_nexus: write detector array into current NXdata (1D,2D)
*   requires: NXdata to be opened
*******************************************************************************/
int mcdetector_out_array_nexus(NXhandle f, char *part, double *data, MCDETECTOR detector)
{
  
  int dims[3]={detector.m,detector.n,detector.p};  /* number of elements to write */
  int signal=1;
  int exists=0;
  int current_dims[3]={0,0,0};
  int ret=NX_OK;
  
  if (!f || !data || !detector.m || mcdisable_output_files) return(NX_OK);
  
  /* when this is a list, we set 1st dimension to NX_UNLIMITED for creation */
  if (strcasestr(detector.format, "list")) dims[0] = NX_UNLIMITED;
  
  /* create the data set in NXdata group */
  NXMDisableErrorReporting(); /* unactivate NeXus error messages, as creation may fail */
  /* NXcompmakedata fails with NX_UNLIMITED */
  if (strcasestr(detector.format, "list"))
    ret = NXmakedata(    f, part, NX_FLOAT64, detector.rank, dims);
  else
    ret = NXcompmakedata(f, part, NX_FLOAT64, detector.rank, dims, NX_COMP_LZW, dims);
  if (ret != NX_OK) {
    /* failed: data set already exists */
    int datatype=0;
    int rank=0;
    exists=1;
    /* inquire current size of data set (nb of events stored) */
    NXopendata(f, part);
    NXgetinfo(f, &rank, current_dims, &datatype);
    NXclosedata(f);
  }
  NXMEnableErrorReporting();  /* re-enable NeXus error messages */
  dims[0] = detector.m; /* restore actual dimension from data writing */
  
  /* open the data set */
  if (NXopendata(f, part) == NX_ERROR) {
    fprintf(stderr, "Warning: could not open DataSet %s '%s' (NeXus)\n",
      part, detector.title);
    return(NX_ERROR);
  }
  if (strcasestr(detector.format, "list")) {
    current_dims[1] = current_dims[2] = 0; /* set starting location for writing slab */
    NXputslab(f, data, current_dims, dims);
    if (!exists)
      printf("Events:   \"%s\"\n",  
        strlen(detector.filename) ? detector.filename : detector.component);
  } else {
    NXputdata (f, data);
  }
  
  if (strstr(part,"data") || strstr(part, "events")) {
    NXputattr(f, "signal", &signal, 1, NX_INT32);
    nxprintattr(f, "short_name", detector.filename && strlen(detector.filename) ? 
      detector.filename : detector.component);
  }
  nxprintattr(f, "long_name", "%s '%s'", part, detector.title);
  NXclosedata(f);
  
  return(NX_OK);
} /* mcdetector_out_array_nexus */

/*******************************************************************************
* mcdetector_out_data_nexus: write detector axes+data into current NXdata
*   The data:NXdetector is opened, then filename:NXdata
*   requires: NXentry to be opened
*******************************************************************************/
int mcdetector_out_data_nexus(NXhandle f, MCDETECTOR detector)
{
  char data_name[32];
  
  if (!f || !detector.m || mcdisable_output_files) return(NX_OK);
  
  strcpy_valid(data_name, 
    detector.filename && strlen(detector.filename) ? 
      detector.filename : detector.component);

  /* the NXdetector group has been created in mcinfo_out_nexus (mcsiminfo_init) */
  if (NXopengroup(f, "data", "NXdetector") == NX_OK) {

    /* the NXdata group has been created in mcdatainfo_out_nexus */
    if (NXopengroup(f, data_name, "NXdata") == NX_OK) {
  
      /* write axes, for histogram data sets, not for lists */
      if (!strcasestr(detector.format, "list")) {
        mcdetector_out_axis_nexus(f, detector.xlabel, detector.xvar, 
          1, detector.m, detector.xmin, detector.xmax);
          
        mcdetector_out_axis_nexus(f, detector.ylabel, detector.yvar, 
          2, detector.n, detector.ymin, detector.ymax);
          
        mcdetector_out_axis_nexus(f, detector.zlabel, detector.zvar, 
          3, detector.p, detector.zmin, detector.zmax);

      } /* !list */
      
      /* write the actual data (appended if already exists) */
      if (!strcasestr(detector.format, "list")) {
        mcdetector_out_array_nexus(f, "data", detector.p1, detector);
        mcdetector_out_array_nexus(f, "errors", detector.p2, detector);
        mcdetector_out_array_nexus(f, "ncount", detector.p0, detector);
      } else
        mcdetector_out_array_nexus(  f, "events", detector.p1, detector);
      
      NXclosegroup(f);
    } /* NXdata */
    NXclosegroup(f);
  } /* NXdetector */
  
  return(NX_OK);
} /* mcdetector_out_array_nexus */

#ifdef USE_MPI
/*******************************************************************************
* mcdetector_out_list_slaves: slaves send their list data to master which writes
*   requires: NXentry to be opened
* WARNING: this method has a flaw: it requires all nodes to flush the lists
*   the same number of times. In case one node is just below the buffer size
*   when finishing (e.g. monitor_nd), it may not trigger save but others may. 
*   Then the number of recv/send is not constant along nodes, and simulation stalls.  
*******************************************************************************/
MCDETECTOR mcdetector_out_list_slaves(MCDETECTOR detector)
{
  int     node_i=0;
  MPI_MASTER(
	     printf("\n** MPI master gathering slave node list data ** \n");
  );
  
  if (mpi_node_rank != mpi_node_root) {
    /* MPI slave: slaves send their data to master: 2 MPI_Send calls */
    /* m, n, p must be sent first, since all slaves do not have the same number of events */
    int mnp[3]={detector.m,detector.n,detector.p};

    if (mc_MPI_Send(mnp, 3, MPI_INT, mpi_node_root)!= MPI_SUCCESS)
      fprintf(stderr, "Warning: proc %i to master: MPI_Send mnp list error (mcdetector_out_list_slaves)\n", mpi_node_rank);
    if (!detector.p1
     || mc_MPI_Send(detector.p1, mnp[0]*mnp[1]*mnp[2], MPI_DOUBLE, mpi_node_root) != MPI_SUCCESS)
      fprintf(stderr, "Warning: proc %i to master: MPI_Send p1 list error: mnp=%i (mcdetector_out_list_slaves)\n", mpi_node_rank, abs(mnp[0]*mnp[1]*mnp[2]));
    /* slaves are done: sent mnp and p1 */
    return (detector);
  } /* end slaves */

  /* MPI master: receive data from slaves sequentially: 2 MPI_Recv calls */

  if (mpi_node_rank == mpi_node_root) {
    for(node_i=0; node_i<mpi_node_count; node_i++) {
      double *this_p1=NULL;                               /* buffer to hold the list from slaves */
      int     mnp[3]={0,0,0};  /* size of this buffer */
      if (node_i != mpi_node_root) { /* get data from slaves */
	if (mc_MPI_Recv(mnp, 3, MPI_INT, node_i) != MPI_SUCCESS)
	  fprintf(stderr, "Warning: master from proc %i: "
		  "MPI_Recv mnp list error (mcdetector_write_data)\n", node_i);
	if (mnp[0]*mnp[1]*mnp[2]) {
	  this_p1 = (double *)calloc(mnp[0]*mnp[1]*mnp[2], sizeof(double));
	  if (!this_p1 || mc_MPI_Recv(this_p1, abs(mnp[0]*mnp[1]*mnp[2]), MPI_DOUBLE, node_i)!= MPI_SUCCESS)
	    fprintf(stderr, "Warning: master from proc %i: "
		    "MPI_Recv p1 list error: mnp=%i (mcdetector_write_data)\n", node_i, mnp[0]*mnp[1]*mnp[2]);
	  else {
	    printf(". MPI master writing data for slave node %i\n",node_i);
	    detector.p1 = this_p1;
	    detector.m  = mnp[0]; detector.n  = mnp[1]; detector.p  = mnp[2];
	    
	    mcdetector_out_data_nexus(nxhandle, detector);
	  }
	}
      } /* if not master */
    } /* for */
  MPI_MASTER(
	     printf("\n** Done ** \n");
  );   
  }
}
#endif

MCDETECTOR mcdetector_out_0D_nexus(MCDETECTOR detector)
{
  /* Write data set information to NeXus file. */
  MPI_MASTER(
    mcdatainfo_out_nexus(nxhandle, detector);
  );
  
  return(detector);
} /* mcdetector_out_0D_ascii */

MCDETECTOR mcdetector_out_1D_nexus(MCDETECTOR detector_inc)
{
  MCDETECTOR detector = detector_inc;
  MPI_MASTER(
  mcdatainfo_out_nexus(nxhandle, detector);
  mcdetector_out_data_nexus(nxhandle, detector);
  );
  return(detector);
} /* mcdetector_out_1D_ascii */

MCDETECTOR mcdetector_out_2D_nexus(MCDETECTOR detector_inc)
{
  MCDETECTOR detector = detector_inc;
  MPI_MASTER(
  mcdatainfo_out_nexus(nxhandle, detector);
  mcdetector_out_data_nexus(nxhandle, detector);
  );
  
#ifdef USE_MPI // and USE_NEXUS
  /* NeXus: slave nodes have master write their lists */
  if (strcasestr(detector.format, "list") && mpi_node_count > 1) {
    mcdetector_out_list_slaves(detector);
  }
#endif /* USE_MPI */

  return(detector);
} /* mcdetector_out_2D_nexus */

#endif /* USE_NEXUS*/








/* ========================================================================== */

/*                            Main input functions                            */
/*            DETECTOR_OUT_xD function calls -> ascii or NeXus                */

/* ========================================================================== */

/*******************************************************************************
* mcsiminfo_init:   open SIM and write header
*******************************************************************************/
FILE *mcsiminfo_init(FILE *f)
{
  int exists=0;
  int index;
  
  /* check format */      
  if (!mcformat || !strlen(mcformat) 
   || !strcasecmp(mcformat, "MCSTAS") || !strcasecmp(mcformat, "MCXTRACE") 
   || !strcasecmp(mcformat, "PGPLOT") || !strcasecmp(mcformat, "GNUPLOT") || !strcasecmp(mcformat, "MCCODE")
   || !strcasecmp(mcformat, "MATLAB")) {
    mcformat="McCode";
#ifdef USE_NEXUS
  } else if (strcasestr(mcformat, "NeXus")) {
    /* Do nothing */
#endif
  } else {
    fprintf(stderr,
	    "Warning: You have requested the output format %s which is unsupported by this binary. Resetting to standard %s format.\n",mcformat ,"McCode");
    mcformat="McCode";
  }
  
  /* open the SIM file if not defined yet */
  if (mcsiminfo_file || mcdisable_output_files) 
    return (mcsiminfo_file);
    
#ifdef USE_NEXUS
  /* only master writes NeXus header: calls NXopen(nxhandle) */
  if (mcformat && strcasestr(mcformat, "NeXus")) {
	  MPI_MASTER(
	  mcsiminfo_file = mcnew_file(mcsiminfo_name, "h5", &exists);
    if(!mcsiminfo_file)
      fprintf(stderr,
	      "Warning: could not open simulation description file '%s'\n",
	      mcsiminfo_name);
	  else
	    mcinfo_out_nexus(nxhandle);
	  );
    return(mcsiminfo_file); /* points to nxhandle */
  }
#endif
  
  /* write main description file (only MASTER) */
  MPI_MASTER(

  mcsiminfo_file = mcnew_file(mcsiminfo_name, "sim", &exists);
  if(!mcsiminfo_file)
    fprintf(stderr,
	    "Warning: could not open simulation description file '%s'\n",
	    mcsiminfo_name);
  else
  {
    /* write SIM header */
    time_t t=time(NULL);
    mcsiminfo_out("%s simulation description file for %s.\n", 
      MCCODE_NAME, mcinstrument_name);
    mcsiminfo_out("Date:    %s", ctime(&t)); /* includes \n */
    mcsiminfo_out("Program: %s\n\n", MCCODE_STRING);
    
    mcsiminfo_out("begin instrument: %s\n", mcinstrument_name);
    mcinfo_out(   "  ", mcsiminfo_file);
    mcsiminfo_out("end instrument\n");

    mcsiminfo_out("\nbegin simulation: %s\n", mcdirname);
    mcruninfo_out("  ", mcsiminfo_file);
    mcsiminfo_out("end simulation\n");

  }
  return (mcsiminfo_file);
  
  ); /* MPI_MASTER */
  
} /* mcsiminfo_init */

/*******************************************************************************
*   mcsiminfo_close:  close SIM
*******************************************************************************/
void mcsiminfo_close()
{
  MPI_MASTER(
  if(mcsiminfo_file && !mcdisable_output_files) {
#ifdef USE_NEXUS
    if (mcformat && strcasestr(mcformat, "NeXus")) {
      time_t t=time(NULL);
      nxprintf(nxhandle, "end_time", ctime(&t));
      nxprintf(nxhandle, "duration", "%li", (long)t-mcstartdate);
      NXclosegroup(nxhandle); /* NXentry */
      NXclose(&nxhandle);
    } else
#endif
      fclose(mcsiminfo_file);
    );
    mcsiminfo_file = NULL;
  }
} /* mcsiminfo_close */

/*******************************************************************************
* mcdetector_out_0D: wrapper for 0D (single value).
*   Output single detector/monitor data (p0, p1, p2).
*   Title is t, component name is c.
*******************************************************************************/
MCDETECTOR mcdetector_out_0D(char *t, double p0, double p1, double p2,
                         char *c, Coords posa)
{
  /* import and perform basic detector analysis (and handle MPI reduce) */
  MCDETECTOR detector = mcdetector_import(mcformat,
    c, (t ? t : MCCODE_STRING " data"),
    1, 1, 1,
    "I", "", "",
    "I", "", "",
    0, 0, 0, 0, 0, 0, "",
    &p0, &p1, &p2, posa); /* write Detector: line */

#ifdef USE_NEXUS
  if (strcasestr(detector.format, "NeXus"))
    return(mcdetector_out_0D_nexus(detector));
  else
#endif
    return(mcdetector_out_0D_ascii(detector));
    
} /* mcdetector_out_0D */



/*******************************************************************************
* mcdetector_out_1D: wrapper for 1D.
*   Output 1d detector data (p0, p1, p2) for n bins linearly
*   distributed across the range x1..x2 (x1 is lower limit of first
*   bin, x2 is upper limit of last bin). Title is t, axis labels are xl
*   and yl. File name is f, component name is c.
*******************************************************************************/
MCDETECTOR mcdetector_out_1D(char *t, char *xl, char *yl,
        char *xvar, double x1, double x2,
        long n,
        double *p0, double *p1, double *p2, char *f,
        char *c, Coords posa)
{
  /* import and perform basic detector analysis (and handle MPI_Reduce) */
  MCDETECTOR detector = mcdetector_import(mcformat,
    c, (t ? t : MCCODE_STRING " 1D data"),
    n, 1, 1,
    xl, yl, (n > 1 ? "Signal per bin" : " Signal"),
    xvar, "(I,I_err)", "I",
    x1, x2, 0, 0, 0, 0, f,
    p0, p1, p2, posa); /* write Detector: line */
  if (!detector.p1 || !detector.m) return(detector);

#ifdef USE_NEXUS
  if (strcasestr(detector.format, "NeXus"))
    return(mcdetector_out_1D_nexus(detector));
  else
#endif
    return(mcdetector_out_1D_ascii(detector));
  
} /* mcdetector_out_1D */

/*******************************************************************************
* mcdetector_out_2D: wrapper for 2D.
*   special case for list: master creates file first, then slaves append their blocks without header
*******************************************************************************/
MCDETECTOR mcdetector_out_2D(char *t, char *xl, char *yl,
                  double x1, double x2, double y1, double y2,
                  long m, long n,
                  double *p0, double *p1, double *p2, char *f,
                  char *c, Coords posa)
{
  char xvar[CHAR_BUF_LENGTH];
  char yvar[CHAR_BUF_LENGTH];
  
  /* create short axes labels */
  if (xl && strlen(xl)) { strncpy(xvar, xl, CHAR_BUF_LENGTH); xvar[strcspn(xvar,"\n\r ")]='\0'; }
  else strcpy(xvar, "x");
  if (yl && strlen(yl)) { strncpy(yvar, yl, CHAR_BUF_LENGTH); yvar[strcspn(yvar,"\n\r ")]='\0'; }
  else strcpy(yvar, "y");

  MCDETECTOR detector;

  /* import and perform basic detector analysis (and handle MPI_Reduce) */
  if (labs(m) == 1) {/* n>1 on Y, m==1 on X: 1D, no X axis*/
    detector = mcdetector_import(mcformat,
      c, (t ? t : MCCODE_STRING " 1D data"),
      n, 1, 1,
      yl, "", "Signal per bin",
      yvar, "(I,Ierr)", "I",
      y1, y2, x1, x2, 0, 0, f,
      p0, p1, p2, posa); /* write Detector: line */
  } else if (labs(n)==1) {/* m>1 on X, n==1 on Y: 1D, no Y axis*/
    detector = mcdetector_import(mcformat,
      c, (t ? t : MCCODE_STRING " 1D data"),
      m, 1, 1,
      xl, "", "Signal per bin",
      xvar, "(I,Ierr)", "I",
      x1, x2, y1, y2, 0, 0, f,
      p0, p1, p2, posa); /* write Detector: line */
  }else {
    detector = mcdetector_import(mcformat,
      c, (t ? t : MCCODE_STRING " 2D data"),
      m, n, 1,
      xl, yl, "Signal per bin",
      xvar, yvar, "I",
      x1, x2, y1, y2, 0, 0, f,
      p0, p1, p2, posa); /* write Detector: line */
  }

  if (!detector.p1 || !detector.m) return(detector);

#ifdef USE_NEXUS
  if (strcasestr(detector.format, "NeXus"))
    return(mcdetector_out_2D_nexus(detector));
  else
#endif
    return(mcdetector_out_2D_ascii(detector));
  
} /* mcdetector_out_2D */

/*******************************************************************************
* mcdetector_out_list: wrapper for list output (calls out_2D with mcformat+"list").
*   m=number of events, n=size of each event
*******************************************************************************/
MCDETECTOR mcdetector_out_list(char *t, char *xl, char *yl,
                  long m, long n,
                  double *p1, char *f,
                  char *c, Coords posa)
{
  char       format_new[CHAR_BUF_LENGTH];
  char      *format_org;
  MCDETECTOR detector;
  
  format_org = mcformat;
  strcpy(format_new, mcformat);
  strcat(format_new, " list");
  mcformat = format_new;

  detector = mcdetector_out_2D(t, xl, yl,
                  1,labs(m),1,labs(n),
                  m,n,
                  NULL, p1, NULL, f,
                  c, posa);
  
  mcformat = format_org;
  return(detector);
}

/*******************************************************************************
 * mcuse_dir: set data/sim storage directory and create it,
 * or exit with error if exists
 ******************************************************************************/
static void
mcuse_dir(char *dir)
{
  if (!dir || !strlen(dir)) return;
#ifdef MC_PORTABLE
  fprintf(stderr, "Error: "
          "Directory output cannot be used with portable simulation (mcuse_dir)\n");
  exit(1);
#else  /* !MC_PORTABLE */
  /* handle file://directory URL type */
  if (strncmp(dir, "file://", strlen("file://")))
    mcdirname = dir;
  else
    mcdirname = dir+strlen("file://");
  
  
  
  MPI_MASTER(
    if(mkdir(mcdirname, 0777)) {
#ifndef DANSE
      fprintf(stderr, "Error: unable to create directory '%s' (mcuse_dir)\n", dir);
      fprintf(stderr, "(Maybe the directory already exists?)\n");
#endif
#ifdef USE_MPI
    MPI_Abort(MPI_COMM_WORLD, -1);
#endif
    exit(-1);
    }
  ); /* MPI_MASTER */
  
  /* remove trailing PATHSEP (if any) */
  while (strlen(mcdirname) && mcdirname[strlen(mcdirname) - 1] == MC_PATHSEP_C)
    mcdirname[strlen(mcdirname) - 1]='\0';
#endif /* !MC_PORTABLE */
} /* mcuse_dir */

/*******************************************************************************
* mcinfo: display instrument simulation info to stdout and exit
*******************************************************************************/
static void
mcinfo(void)
{
  fprintf(stdout, "begin instrument: %s\n", mcinstrument_name);
  mcinfo_out("  ", stdout);
  fprintf(stdout, "end instrument\n");
  fprintf(stdout, "begin simulation: %s\n", mcdirname ? mcdirname : ".");
  mcruninfo_out_backend("  ", stdout,1);
  fprintf(stdout, "end simulation\n");
  exit(0); /* includes MPI_Finalize in MPI mode */
} /* mcinfo */

#endif /* ndef MCCODE_R_IO_C */

/* end of the I/O section =================================================== */







/*******************************************************************************
* mcset_ncount: set total number of rays to generate
*******************************************************************************/
void mcset_ncount(unsigned long long int count)
{
  mcncount = count;
}

/* mcget_ncount: get total number of rays to generate */
unsigned long long int mcget_ncount(void)
{
  return mcncount;
}

/* mcget_run_num: get curent number of rays in TRACE */
unsigned long long int mcget_run_num(void)
{
  return mcrun_num;
}

/* mcsetn_arg: get ncount from a string argument */
static void
mcsetn_arg(char *arg)
{
  mcset_ncount((long long int) strtod(arg, NULL));
}

/* mcsetseed: set the random generator seed from a string argument */
static void
mcsetseed(char *arg)
{
  mcseed = atol(arg);
  if(mcseed) {
    srandom(mcseed);
  } else {
    fprintf(stderr, "Error: seed must not be zero (mcsetseed)\n");
    exit(1);
  }
}

/* Following part is only embedded when not redundent with mccode-r.h ========= */

#ifndef MCCODE_H

/* SECTION: MCDISPLAY support. =============================================== */

/*******************************************************************************
* Just output MCDISPLAY keywords to be caught by an external plotter client.
*******************************************************************************/

void mcdis_magnify(char *what){
  // Do nothing here, better use interactive zoom from the tools
}

void mcdis_line(double x1, double y1, double z1,
                double x2, double y2, double z2){
  printf("MCDISPLAY: multiline(2,%g,%g,%g,%g,%g,%g)\n",
         x1,y1,z1,x2,y2,z2);
}

void mcdis_dashed_line(double x1, double y1, double z1,
		       double x2, double y2, double z2, int n){
  int i;
  const double dx = (x2-x1)/(2*n+1);
  const double dy = (y2-y1)/(2*n+1);
  const double dz = (z2-z1)/(2*n+1);

  for(i = 0; i < n+1; i++)
    mcdis_line(x1 + 2*i*dx,     y1 + 2*i*dy,     z1 + 2*i*dz,
	       x1 + (2*i+1)*dx, y1 + (2*i+1)*dy, z1 + (2*i+1)*dz);
}

void mcdis_multiline(int count, ...){
  va_list ap;
  double x,y,z;

  printf("MCDISPLAY: multiline(%d", count);
  va_start(ap, count);
  while(count--)
    {
    x = va_arg(ap, double);
    y = va_arg(ap, double);
    z = va_arg(ap, double);
    printf(",%g,%g,%g", x, y, z);
    }
  va_end(ap);
  printf(")\n");
}

void mcdis_rectangle(char* plane, double x, double y, double z,
		     double width, double height){
  /* draws a rectangle in the plane           */
  /* x is ALWAYS width and y is ALWAYS height */
  if (strcmp("xy", plane)==0) {
    mcdis_multiline(5,
		    x - width/2, y - height/2, z,
		    x + width/2, y - height/2, z,
		    x + width/2, y + height/2, z,
		    x - width/2, y + height/2, z,
		    x - width/2, y - height/2, z);
  } else if (strcmp("xz", plane)==0) {
    mcdis_multiline(5,
		    x - width/2, y, z - height/2,
		    x + width/2, y, z - height/2,
		    x + width/2, y, z + height/2,
		    x - width/2, y, z + height/2,
		    x - width/2, y, z - height/2);
  } else if (strcmp("yz", plane)==0) {
    mcdis_multiline(5,
		    x, y - height/2, z - width/2,
		    x, y - height/2, z + width/2,
		    x, y + height/2, z + width/2,
		    x, y + height/2, z - width/2,
		    x, y - height/2, z - width/2);
  } else {

    fprintf(stderr, "Error: Definition of plane %s unknown\n", plane);
    exit(1);
  }
}

/*  draws a box with center at (x, y, z) and
    width (deltax), height (deltay), length (deltaz) */
void mcdis_box(double x, double y, double z,
	       double width, double height, double length){

  mcdis_rectangle("xy", x, y, z-length/2, width, height);
  mcdis_rectangle("xy", x, y, z+length/2, width, height);
  mcdis_line(x-width/2, y-height/2, z-length/2,
	     x-width/2, y-height/2, z+length/2);
  mcdis_line(x-width/2, y+height/2, z-length/2,
	     x-width/2, y+height/2, z+length/2);
  mcdis_line(x+width/2, y-height/2, z-length/2,
	     x+width/2, y-height/2, z+length/2);
  mcdis_line(x+width/2, y+height/2, z-length/2,
	     x+width/2, y+height/2, z+length/2);
}

void mcdis_circle(char *plane, double x, double y, double z, double r){
  printf("MCDISPLAY: circle('%s',%g,%g,%g,%g)\n", plane, x, y, z, r);
}

/* Draws a circle with center (x,y,z), radius (r), and in the plane
 * with normal (nx,ny,nz)*/
void mcdis_Circle(double x, double y, double z, double r, double nx, double ny, double nz){
    int i;
    if(nx==0 && ny && nz==0){
        for (i=0;i<24; i++){
            mcdis_line(x+r*sin(i*2*M_PI/24),y,z+r*cos(i*2*M_PI/24),
                    x+r*sin((i+1)*2*M_PI/24),y,z+r*cos((i+1)*2*M_PI/24));
        }
    }else{
        double mx,my,mz;
        /*generate perpendicular vector using (nx,ny,nz) and (0,1,0)*/
        vec_prod(mx,my,mz, 0,1,0, nx,ny,nz);
        NORM(mx,my,mz);
        /*draw circle*/
        for (i=0;i<24; i++){
            double ux,uy,uz;
            double wx,wy,wz;
            rotate(ux,uy,uz, mx,my,mz, i*2*M_PI/24, nx,ny,nz);
            rotate(wx,wy,wz, mx,my,mz, (i+1)*2*M_PI/24, nx,ny,nz);
            mcdis_line(x+ux*r,y+uy*r,z+uz*r,
                    x+wx*r,y+wy*r,z+wz*r);
        }
    }
}

/* Draws a cylinder with center at (x,y,z) with extent (r,height).
 * The cylinder axis is along the vector nx,ny,nz.
 * N determines how many vertical lines are drawn.*/
void mcdis_cylinder( double x, double y, double z,
        double r, double height, int N, double nx, double ny, double nz){
    int i;
    /*no lines make little sense - so trigger the default*/
    if(N<=0) N=5;

    NORM(nx,ny,nz);
    double h_2=height/2.0;
    mcdis_Circle(x+nx*h_2,y+ny*h_2,z+nz*h_2,r,nx,ny,nz);
    mcdis_Circle(x-nx*h_2,y-ny*h_2,z-nz*h_2,r,nx,ny,nz);

    double mx,my,mz;
    /*generate perpendicular vector using (nx,ny,nz) and (0,1,0)*/
    if(nx==0 && ny && nz==0){
        mx=my=0;mz=1;
    }else{
        vec_prod(mx,my,mz, 0,1,0, nx,ny,nz);
        NORM(mx,my,mz);
    }
    /*draw circle*/
    for (i=0; i<24; i++){
        double ux,uy,uz;
        rotate(ux,uy,uz, mx,my,mz, i*2*M_PI/24, nx,ny,nz);
        mcdis_line(x+nx*h_2+ux*r, y+ny*h_2+uy*r, z+nz*h_2+uz*r,
                 x-nx*h_2+ux*r, y-ny*h_2+uy*r, z-nz*h_2+uz*r);
    }
}

/* draws a sphere with center at (x,y,z) with extent (r)
 * The sphere is drawn using N longitudes and N latitudes.*/
void mcdis_sphere(double x, double y, double z, double r, int N){
    double nx,ny,nz;
    int i;
    /*no lines make little sense - so trigger the default*/
    if(N<=0) N=5;

    nx=0;ny=0;nz=1;
    mcdis_Circle(x,y,z,r,nx,ny,nz);
    for (i=1;i<N;i++){
        rotate(nx,ny,nz, nx,ny,nz, M_PI/N, 0,1,0);
        mcdis_Circle(x,y,z,r,nx,ny,nz);
    }
    /*lastly draw a great circle perpendicular to all N circles*/
    //mcdis_Circle(x,y,z,radius,1,0,0);

    for (i=1;i<=N;i++){
        double yy=-r+ 2*r*((double)i/(N+1));
        mcdis_Circle(x,y+yy ,z,  sqrt(r*r-yy*yy) ,0,1,0);
    }
}

/* SECTION: coordinates handling ============================================ */

/*******************************************************************************
* Since we use a lot of geometric calculations using Cartesian coordinates,
* we collect some useful routines here. However, it is also permissible to
* work directly on the underlying struct coords whenever that is most
* convenient (that is, the type Coords is not abstract).
*
* Coordinates are also used to store rotation angles around x/y/z axis.
*
* Since coordinates are used much like a basic type (such as double), the
* structure itself is passed and returned, rather than a pointer.
*
* At compile-time, the values of the coordinates may be unknown (for example
* a motor position). Hence coordinates are general expressions and not simple
* numbers. For this we used the type Coords_exp which has three CExp
* fields. For runtime (or calculations possible at compile time), we use
* Coords which contains three double fields.
*******************************************************************************/

/* coords_set: Assign coordinates. */
Coords
coords_set(MCNUM x, MCNUM y, MCNUM z)
{
  Coords a;

  a.x = x;
  a.y = y;
  a.z = z;
  return a;
}

/* coords_get: get coordinates. Required when 'x','y','z' are #defined as ray pars */
Coords
coords_get(Coords a, MCNUM *x, MCNUM *y, MCNUM *z)
{
  *x = a.x;
  *y = a.y;
  *z = a.z;
  return a;
}

/* coords_add: Add two coordinates. */
Coords
coords_add(Coords a, Coords b)
{
  Coords c;

  c.x = a.x + b.x;
  c.y = a.y + b.y;
  c.z = a.z + b.z;
  if (fabs(c.z) < 1e-14) c.z=0.0;
  return c;
}

/* coords_sub: Subtract two coordinates. */
Coords
coords_sub(Coords a, Coords b)
{
  Coords c;

  c.x = a.x - b.x;
  c.y = a.y - b.y;
  c.z = a.z - b.z;
  if (fabs(c.z) < 1e-14) c.z=0.0;
  return c;
}

/* coords_neg: Negate coordinates. */
Coords
coords_neg(Coords a)
{
  Coords b;

  b.x = -a.x;
  b.y = -a.y;
  b.z = -a.z;
  return b;
}

/* coords_scale: Scale a vector. */
Coords coords_scale(Coords b, double scale) {
  Coords a;

  a.x = b.x*scale;
  a.y = b.y*scale;
  a.z = b.z*scale;
  return a;
}

/* coords_sp: Scalar product: a . b */
double coords_sp(Coords a, Coords b) {
  double value;

  value = a.x*b.x + a.y*b.y + a.z*b.z;
  return value;
}

/* coords_xp: Cross product: a = b x c. */
Coords coords_xp(Coords b, Coords c) {
  Coords a;

  a.x = b.y*c.z - c.y*b.z;
  a.y = b.z*c.x - c.z*b.x;
  a.z = b.x*c.y - c.x*b.y;
  return a;
}

/* coords_len: Gives length of coords set. */
double coords_len(Coords a) {
  return sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
}

/* coords_mirror: Mirror a in plane (through the origin) defined by normal n*/
Coords coords_mirror(Coords a, Coords n) {
  double t = scalar_prod(n.x, n.y, n.z, n.x, n.y, n.z);
  Coords b;
  if (t!=1) {
    t = sqrt(t);
    n.x /= t;
    n.y /= t;
    n.z /= t;
  }
  t=scalar_prod(a.x, a.y, a.z, n.x, n.y, n.z);
  b.x = a.x-2*t*n.x;
  b.y = a.y-2*t*n.y;
  b.z = a.z-2*t*n.z;
  return b;
}

/* coords_print: Print out vector values. */
void coords_print(Coords a) {

  fprintf(stdout, "(%f, %f, %f)\n", a.x, a.y, a.z);
  return;
}

mcstatic void coords_norm(Coords* c) {
	double temp = coords_sp(*c,*c);

	// Skip if we will end dividing by zero
	if (temp == 0) return;

	temp = sqrt(temp);

	c->x /= temp;
	c->y /= temp;
	c->z /= temp;
}

/*******************************************************************************
* The Rotation type implements a rotation transformation of a coordinate
* system in the form of a double[3][3] matrix.
*
* Contrary to the Coords type in coords.c, rotations are passed by
* reference. Functions that yield new rotations do so by writing to an
* explicit result parameter; rotations are not returned from functions. The
* reason for this is that arrays cannot by returned from functions (though
* structures can; thus an alternative would have been to wrap the
* double[3][3] array up in a struct). Such are the ways of C programming.
*
* A rotation represents the tranformation of the coordinates of a vector when
* changing between coordinate systems that are rotated with respect to each
* other. For example, suppose that coordinate system Q is rotated 45 degrees
* around the Z axis with respect to coordinate system P. Let T be the
* rotation transformation representing a 45 degree rotation around Z. Then to
* get the coordinates of a vector r in system Q, apply T to the coordinates
* of r in P. If r=(1,0,0) in P, it will be (sqrt(1/2),-sqrt(1/2),0) in
* Q. Thus we should be careful when interpreting the sign of rotation angles:
* they represent the rotation of the coordinate systems, not of the
* coordinates (which has opposite sign).
*******************************************************************************/

/*******************************************************************************
* rot_set_rotation: Get transformation for rotation first phx around x axis,
* then phy around y, then phz around z.
*******************************************************************************/
void
rot_set_rotation(Rotation t, double phx, double phy, double phz)
{
  if ((phx == 0) && (phy == 0) && (phz == 0)) {
    t[0][0] = 1.0;
    t[0][1] = 0.0;
    t[0][2] = 0.0;
    t[1][0] = 0.0;
    t[1][1] = 1.0;
    t[1][2] = 0.0;
    t[2][0] = 0.0;
    t[2][1] = 0.0;
    t[2][2] = 1.0;
  } else {
    double cx = cos(phx);
    double sx = sin(phx);
    double cy = cos(phy);
    double sy = sin(phy);
    double cz = cos(phz);
    double sz = sin(phz);

    t[0][0] = cy*cz;
    t[0][1] = sx*sy*cz + cx*sz;
    t[0][2] = sx*sz - cx*sy*cz;
    t[1][0] = -cy*sz;
    t[1][1] = cx*cz - sx*sy*sz;
    t[1][2] = sx*cz + cx*sy*sz;
    t[2][0] = sy;
    t[2][1] = -sx*cy;
    t[2][2] = cx*cy;
  }
}

/*******************************************************************************
* rot_test_identity: Test if rotation is identity
*******************************************************************************/
int
rot_test_identity(Rotation t)
{
  return (t[0][0] + t[1][1] + t[2][2] == 3);
}

/*******************************************************************************
* rot_mul: Matrix multiplication of transformations (this corresponds to
* combining transformations). After rot_mul(T1, T2, T3), doing T3 is
* equal to doing first T2, then T1.
* Note that T3 must not alias (use the same array as) T1 or T2.
*******************************************************************************/
void
rot_mul(Rotation t1, Rotation t2, Rotation t3)
{
  if (rot_test_identity(t1)) {
    rot_copy(t3, t2);
  } else if (rot_test_identity(t2)) {
    rot_copy(t3, t1);
  } else {
    int i,j;
    for(i = 0; i < 3; i++)
      for(j = 0; j < 3; j++)
	t3[i][j] = t1[i][0]*t2[0][j] + t1[i][1]*t2[1][j] + t1[i][2]*t2[2][j];
  }
}

/*******************************************************************************
* rot_copy: Copy a rotation transformation (arrays cannot be assigned in C).
*******************************************************************************/
void
rot_copy(Rotation dest, Rotation src)
{
  int i,j;
  for(i = 0; i < 3; i++)
    for(j = 0; j < 3; j++)
      dest[i][j] = src[i][j];
}

/*******************************************************************************
* rot_transpose: Matrix transposition, which is inversion for Rotation matrices
*******************************************************************************/
void
rot_transpose(Rotation src, Rotation dst)
{
  dst[0][0] = src[0][0];
  dst[0][1] = src[1][0];
  dst[0][2] = src[2][0];
  dst[1][0] = src[0][1];
  dst[1][1] = src[1][1];
  dst[1][2] = src[2][1];
  dst[2][0] = src[0][2];
  dst[2][1] = src[1][2];
  dst[2][2] = src[2][2];
}

/*******************************************************************************
* rot_apply: returns t*a
*******************************************************************************/
Coords
rot_apply(Rotation t, Coords a)
{
  Coords b;
  if (rot_test_identity(t)) {
    return a;
  } else {
    b.x = t[0][0]*a.x + t[0][1]*a.y + t[0][2]*a.z;
    b.y = t[1][0]*a.x + t[1][1]*a.y + t[1][2]*a.z;
    b.z = t[2][0]*a.x + t[2][1]*a.y + t[2][2]*a.z;
    return b;
  }
}

/**
 * Pretty-printing of rotation matrices.
 */
void rot_print(Rotation rot) {
	printf("[ %4.2f %4.2f %4.2f ]\n",
			rot[0][0], rot[0][1], rot[0][2]);
	printf("[ %4.2f %4.2f %4.2f ]\n",
			rot[1][0], rot[1][1], rot[1][2]);
	printf("[ %4.2f %4.2f %4.2f ]\n\n",
			rot[2][0], rot[2][1], rot[2][2]);
}

/**
 * Vector product: used by vec_prod (mccode-r.h). Use coords_xp for Coords.
 */
mcstatic void vec_prod_func(double *x, double *y, double *z,
		double x1, double y1, double z1,
		double x2, double y2, double z2) {
    *x = (y1)*(z2) - (y2)*(z1);
    *y = (z1)*(x2) - (z2)*(x1);
    *z = (x1)*(y2) - (x2)*(y1);
}

/**
 * Scalar product: use coords_sp for Coords.
 */
mcstatic double scalar_prod(
		double x1, double y1, double z1,
		double x2, double y2, double z2) {
	return ((x1 * x2) + (y1 * y2) + (z1 * z2));
}

/*******************************************************************************
* mccoordschange: applies rotation to (x y z) and (vx vy vz) and Spin (sx,sy,sz)
*******************************************************************************/
void
mccoordschange(Coords a, Rotation t, double *x, double *y, double *z,
               double *vx, double *vy, double *vz, double *sx, double *sy, double *sz)
{
  Coords b, c;

  b.x = *x;
  b.y = *y;
  b.z = *z;
  c = rot_apply(t, b);
  b = coords_add(c, a);
  *x = b.x;
  *y = b.y;
  *z = b.z;

  if ( (vz && vy  && vx) && (*vz != 0.0 || *vx != 0.0 || *vy != 0.0) ) mccoordschange_polarisation(t, vx, vy, vz);

  if ( (sz && sy  && sx) && (*sz != 0.0 || *sx != 0.0 || *sy != 0.0) ) mccoordschange_polarisation(t, sx, sy, sz);

}

/*******************************************************************************
* mccoordschange_polarisation: applies rotation to vector (sx sy sz)
*******************************************************************************/
void
mccoordschange_polarisation(Rotation t, double *sx, double *sy, double *sz)
{
  Coords b, c;

  b.x = *sx;
  b.y = *sy;
  b.z = *sz;
  c = rot_apply(t, b);
  *sx = c.x;
  *sy = c.y;
  *sz = c.z;
}

/* SECTION: vector math  ==================================================== */

/* normal_vec_func: Compute normal vector to (x,y,z). */
mcstatic void normal_vec_func(double *nx, double *ny, double *nz,
                double x, double y, double z)
{
  double ax = fabs(x);
  double ay = fabs(y);
  double az = fabs(z);
  double l;
  if(x == 0 && y == 0 && z == 0)
  {
    *nx = 0;
    *ny = 0;
    *nz = 0;
    return;
  }
  if(ax < ay)
  {
    if(ax < az)
    {                           /* Use X axis */
      l = sqrt(z*z + y*y);
      *nx = 0;
      *ny = z/l;
      *nz = -y/l;
      return;
    }
  }
  else
  {
    if(ay < az)
    {                           /* Use Y axis */
      l = sqrt(z*z + x*x);
      *nx = z/l;
      *ny = 0;
      *nz = -x/l;
      return;
    }
  }
  /* Use Z axis */
  l = sqrt(y*y + x*x);
  *nx = y/l;
  *ny = -x/l;
  *nz = 0;
} /* normal_vec */

/*******************************************************************************
 * solve_2nd_order: second order equation solve: A*t^2 + B*t + C = 0
 * solve_2nd_order(&t1, NULL, A,B,C)
 *   returns 0 if no solution was found, or set 't1' to the smallest positive
 *   solution.
 * solve_2nd_order(&t1, &t2, A,B,C)
 *   same as with &t2=NULL, but also returns the second solution.
 * EXAMPLE usage for intersection of a trajectory with a plane in gravitation
 * field (gx,gy,gz):
 * The neutron starts at point r=(x,y,z) with velocityv=(vx vy vz). The plane
 * has a normal vector n=(nx,ny,nz) and contains the point W=(wx,wy,wz).
 * The problem consists in solving the 2nd order equation:
 *      1/2.n.g.t^2 + n.v.t + n.(r-W) = 0
 * so that A = 0.5 n.g; B = n.v; C = n.(r-W);
 * Without acceleration, t=-n.(r-W)/n.v
 ******************************************************************************/
int solve_2nd_order(double *t1, double *t2,
                  double A,  double B,  double C)
{
  int ret=0;

  if (!t1) return 0;
  *t1 = 0;
  if (t2) *t2=0;

  if (fabs(A) < 1E-10) /* approximate to linear equation: A ~ 0 */
  {
    if (B) {  *t1 = -C/B; ret=1; if (t2) *t2=*t1; }
    /* else no intersection: A=B=0 ret=0 */
  }
  else
  {
    double D;
    D = B*B - 4*A*C;
    if (D >= 0) /* Delta > 0: two solutions */
    {
      double sD, dt1, dt2;
      sD = sqrt(D);
      dt1 = (-B + sD)/2/A;
      dt2 = (-B - sD)/2/A;
      /* we identify very small values with zero */
      if (fabs(dt1) < 1e-10) dt1=0.0;
      if (fabs(dt2) < 1e-10) dt2=0.0;

      /* now we choose the smallest positive solution */
      if      (dt1<=0.0 && dt2>0.0) ret=2; /* dt2 positive */
      else if (dt2<=0.0 && dt1>0.0) ret=1; /* dt1 positive */
      else if (dt1> 0.0 && dt2>0.0)
      {  if (dt1 < dt2) ret=1; else ret=2; } /* all positive: min(dt1,dt2) */
      /* else two solutions are negative. ret=-1 */
      if (ret==1) { *t1 = dt1;  if (t2) *t2=dt2; }
      else        { *t1 = dt2;  if (t2) *t2=dt1; }
      ret=2;  /* found 2 solutions and t1 is the positive one */
    } /* else Delta <0: no intersection. ret=0 */
  }
  return(ret);
} /* solve_2nd_order */

/*******************************************************************************
 * randvec_target_circle: Choose random direction towards target at (x,y,z)
 * with given radius.
 * If radius is zero, choose random direction in full 4PI, no target.
 ******************************************************************************/
void
randvec_target_circle(double *xo, double *yo, double *zo, double *solid_angle,
               double xi, double yi, double zi, double radius)
{
  double l2, phi, theta, nx, ny, nz, xt, yt, zt, xu, yu, zu;

  if(radius == 0.0)
  {
    /* No target, choose uniformly a direction in full 4PI solid angle. */
    theta = acos (1 - rand0max(2));
    phi = rand0max(2 * PI);
    if(solid_angle)
      *solid_angle = 4*PI;
    nx = 1;
    ny = 0;
    nz = 0;
    yi = sqrt(xi*xi+yi*yi+zi*zi);
    zi = 0;
    xi = 0;
  }
  else
  {
    double costheta0;
    l2 = xi*xi + yi*yi + zi*zi; /* sqr Distance to target. */
    costheta0 = sqrt(l2/(radius*radius+l2));
    if (radius < 0) costheta0 *= -1;
    if(solid_angle)
    {
      /* Compute solid angle of target as seen from origin. */
        *solid_angle = 2*PI*(1 - costheta0);
    }

    /* Now choose point uniformly on circle surface within angle theta0 */
    theta = acos (1 - rand0max(1 - costheta0)); /* radius on circle */
    phi = rand0max(2 * PI); /* rotation on circle at given radius */
    /* Now, to obtain the desired vector rotate (xi,yi,zi) angle theta around a
       perpendicular axis u=i x n and then angle phi around i. */
    if(xi == 0 && zi == 0)
    {
      nx = 1;
      ny = 0;
      nz = 0;
    }
    else
    {
      nx = -zi;
      nz = xi;
      ny = 0;
    }
  }

  /* [xyz]u = [xyz]i x n[xyz] (usually vertical) */
  vec_prod(xu,  yu,  zu, xi, yi, zi,        nx, ny, nz);
  /* [xyz]t = [xyz]i rotated theta around [xyz]u */
  rotate  (xt,  yt,  zt, xi, yi, zi, theta, xu, yu, zu);
  /* [xyz]o = [xyz]t rotated phi around n[xyz] */
  rotate (*xo, *yo, *zo, xt, yt, zt, phi, xi, yi, zi);
} /* randvec_target_circle */

/*******************************************************************************
 * randvec_target_rect_angular: Choose random direction towards target at
 * (xi,yi,zi) with given ANGULAR dimension height x width. height=phi_x=[0,PI],
 * width=phi_y=[0,2*PI] (radians)
 * If height or width is zero, choose random direction in full 4PI, no target.
 *******************************************************************************/
void
randvec_target_rect_angular(double *xo, double *yo, double *zo, double *solid_angle,
               double xi, double yi, double zi, double width, double height, Rotation A)
{
  double theta, phi, nx, ny, nz, xt, yt, zt, xu, yu, zu;
  Coords tmp;
  Rotation Ainverse;

  rot_transpose(A, Ainverse);

  if(height == 0.0 || width == 0.0)
  {
    randvec_target_circle(xo, yo, zo, solid_angle,
               xi, yi, zi, 0);
    return;
  }
  else
  {
    if(solid_angle)
    {
      /* Compute solid angle of target as seen from origin. */
      *solid_angle = 2*fabs(width*sin(height/2));
    }

    /* Go to global coordinate system */

    tmp = coords_set(xi, yi, zi);
    tmp = rot_apply(Ainverse, tmp);
    coords_get(tmp, &xi, &yi, &zi);

    /* Now choose point uniformly on the unit sphere segment with angle theta/phi */
    phi   = width*randpm1()/2.0;
    theta = asin(randpm1()*sin(height/2.0));
    /* Now, to obtain the desired vector rotate (xi,yi,zi) angle theta around
       n, and then phi around u. */
    if(xi == 0 && zi == 0)
    {
      nx = 1;
      ny = 0;
      nz = 0;
    }
    else
    {
      nx = -zi;
      nz = xi;
      ny = 0;
    }
  }

  /* [xyz]u = [xyz]i x n[xyz] (usually vertical) */
  vec_prod(xu,  yu,  zu, xi, yi, zi,        nx, ny, nz);
  /* [xyz]t = [xyz]i rotated theta around [xyz]u */
  rotate  (xt,  yt,  zt, xi, yi, zi, theta, nx, ny, nz);
  /* [xyz]o = [xyz]t rotated phi around n[xyz] */
  rotate (*xo, *yo, *zo, xt, yt, zt, phi, xu,  yu,  zu);

  /* Go back to local coordinate system */
  tmp = coords_set(*xo, *yo, *zo);
  tmp = rot_apply(A, tmp);
  coords_get(tmp, &*xo, &*yo, &*zo);

} /* randvec_target_rect_angular */

/*******************************************************************************
 * randvec_target_rect_real: Choose random direction towards target at (xi,yi,zi)
 * with given dimension height x width (in meters !).
 *
 * Local emission coordinate is taken into account and corrected for 'order' times.
 * (See remarks posted to mcstas-users by George Apostolopoulus <gapost@ipta.demokritos.gr>)
 *
 * If height or width is zero, choose random direction in full 4PI, no target.
 *
 * Traditionally, this routine had the name randvec_target_rect - this is now a
 * a define (see mcstas-r.h) pointing here. If you use the old rouine, you are NOT
 * taking the local emmission coordinate into account.
*******************************************************************************/

void
randvec_target_rect_real(double *xo, double *yo, double *zo, double *solid_angle,
               double xi, double yi, double zi,
               double width, double height, Rotation A,
               double lx, double ly, double lz, int order)
{
  double dx, dy, dist, dist_p, nx, ny, nz, mx, my, mz, n_norm, m_norm;
  double cos_theta;
  Coords tmp;
  Rotation Ainverse;

  rot_transpose(A, Ainverse);

  if(height == 0.0 || width == 0.0)
  {
    randvec_target_circle(xo, yo, zo, solid_angle,
               xi, yi, zi, 0);
    return;
  }
  else
  {

    /* Now choose point uniformly on rectangle within width x height */
    dx = width*randpm1()/2.0;
    dy = height*randpm1()/2.0;

    /* Determine distance to target plane*/
    dist = sqrt(xi*xi + yi*yi + zi*zi);
    /* Go to global coordinate system */

    tmp = coords_set(xi, yi, zi);
    tmp = rot_apply(Ainverse, tmp);
    coords_get(tmp, &xi, &yi, &zi);

    /* Determine vector normal to trajectory axis (z) and gravity [0 1 0] */
    vec_prod(nx, ny, nz, xi, yi, zi, 0, 1, 0);

    /* This now defines the x-axis, normalize: */
    n_norm=sqrt(nx*nx + ny*ny + nz*nz);
    nx = nx/n_norm;
    ny = ny/n_norm;
    nz = nz/n_norm;

    /* Now, determine our y-axis (vertical in many cases...) */
    vec_prod(mx, my, mz, xi, yi, zi, nx, ny, nz);
    m_norm=sqrt(mx*mx + my*my + mz*mz);
    mx = mx/m_norm;
    my = my/m_norm;
    mz = mz/m_norm;

    /* Our output, random vector can now be defined by linear combination: */

    *xo = xi + dx * nx + dy * mx;
    *yo = yi + dx * ny + dy * my;
    *zo = zi + dx * nz + dy * mz;

    /* Go back to local coordinate system */
    tmp = coords_set(*xo, *yo, *zo);
    tmp = rot_apply(A, tmp);
    coords_get(tmp, &*xo, &*yo, &*zo);

    /* Go back to local coordinate system */
    tmp = coords_set(xi, yi, zi);
    tmp = rot_apply(A, tmp);
    coords_get(tmp, &xi, &yi, &zi);

    if (solid_angle) {
      /* Calculate vector from local point to remote random point */
      lx = *xo - lx;
      ly = *yo - ly;
      lz = *zo - lz;
      dist_p = sqrt(lx*lx + ly*ly + lz*lz);

      /* Adjust the 'solid angle' */
      /* 1/r^2 to the chosen point times cos(\theta) between the normal */
      /* vector of the target rectangle and direction vector of the chosen point. */
      cos_theta = (xi * lx + yi * ly + zi * lz) / (dist * dist_p);
      *solid_angle = width * height / (dist_p * dist_p);
      int counter;
      for (counter = 0; counter < order; counter++) {
	*solid_angle = *solid_angle * cos_theta;
      }
    }
  }
} /* randvec_target_rect_real */

/* SECTION: random numbers ================================================== */

/*
 * Copyright (c) 1983 Regents of the University of California.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms are permitted
 * provided that the above copyright notice and this paragraph are
 * duplicated in all such forms and that any documentation,
 * advertising materials, and other materials related to such
 * distribution and use acknowledge that the software was developed
 * by the University of California, Berkeley.  The name of the
 * University may not be used to endorse or promote products derived
 * from this software without specific prior written permission.
 * THIS SOFTWARE IS PROVIDED ``AS IS'' AND WITHOUT ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, WITHOUT LIMITATION, THE IMPLIED
 * WARRANTIES OF MERCHANTIBILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 */

/*
 * This is derived from the Berkeley source:
 *        @(#)random.c        5.5 (Berkeley) 7/6/88
 * It was reworked for the GNU C Library by Roland McGrath.
 * Rewritten to use reentrant functions by Ulrich Drepper, 1995.
 */

/*******************************************************************************
* Modified for McStas from glibc 2.0.7pre1 stdlib/random.c and
* stdlib/random_r.c.
*
* This way random() is more than four times faster compared to calling
* standard glibc random() on ix86 Linux, probably due to multithread support,
* ELF shared library overhead, etc. It also makes McStas generated
* simulations more portable (more likely to behave identically across
* platforms, important for parrallel computations).
*******************************************************************************/


#define        TYPE_3                3
#define        BREAK_3                128
#define        DEG_3                31
#define        SEP_3                3

static mc_int32_t randtbl[DEG_3 + 1] =
  {
    TYPE_3,

    -1726662223, 379960547, 1735697613, 1040273694, 1313901226,
    1627687941, -179304937, -2073333483, 1780058412, -1989503057,
    -615974602, 344556628, 939512070, -1249116260, 1507946756,
    -812545463, 154635395, 1388815473, -1926676823, 525320961,
    -1009028674, 968117788, -123449607, 1284210865, 435012392,
    -2017506339, -911064859, -370259173, 1132637927, 1398500161,
    -205601318,
  };

static mc_int32_t *fptr = &randtbl[SEP_3 + 1];
static mc_int32_t *rptr = &randtbl[1];
static mc_int32_t *state = &randtbl[1];
#define rand_deg DEG_3
#define rand_sep SEP_3
static mc_int32_t *end_ptr = &randtbl[sizeof (randtbl) / sizeof (randtbl[0])];

mc_int32_t
mc_random (void)
{
  mc_int32_t result;

  *fptr += *rptr;
  /* Chucking least random bit.  */
  result = (*fptr >> 1) & 0x7fffffff;
  ++fptr;
  if (fptr >= end_ptr)
  {
    fptr = state;
    ++rptr;
  }
  else
  {
    ++rptr;
    if (rptr >= end_ptr)
      rptr = state;
  }
  return result;
}

void
mc_srandom (unsigned int x)
{
  /* We must make sure the seed is not 0.  Take arbitrarily 1 in this case.  */
  state[0] = x ? x : 1;
  {
    long int i;
    for (i = 1; i < rand_deg; ++i)
    {
      /* This does:
         state[i] = (16807 * state[i - 1]) % 2147483647;
         but avoids overflowing 31 bits.  */
      long int hi = state[i - 1] / 127773;
      long int lo = state[i - 1] % 127773;
      long int test = 16807 * lo - 2836 * hi;
      state[i] = test + (test < 0 ? 2147483647 : 0);
    }
    fptr = &state[rand_sep];
    rptr = &state[0];
    for (i = 0; i < 10 * rand_deg; ++i)
      random ();
  }
}

/* "Mersenne Twister", by Makoto Matsumoto and Takuji Nishimura. */
/* See http://www.math.keio.ac.jp/~matumoto/emt.html for original source. */


/*
   A C-program for MT19937, with initialization improved 2002/1/26.
   Coded by Takuji Nishimura and Makoto Matsumoto.

   Before using, initialize the state by using mt_srandom(seed)
   or init_by_array(init_key, key_length).

   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote
        products derived from this software without specific prior written
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


   Any feedback is very welcome.
   http://www.math.keio.ac.jp/matumoto/emt.html
   email: matumoto@math.keio.ac.jp
*/

#include <stdio.h>

/* Period parameters */
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

static unsigned long mt[N]; /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

/* initializes mt[N] with a seed */
void mt_srandom(unsigned long s)
{
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<N; mti++) {
        mt[mti] =
            (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
void init_by_array(unsigned long init_key[], unsigned long key_length)
{
    int i, j, k;
    mt_srandom(19650218UL);
    i=1; j=0;
    k = (N>key_length ? N : key_length);
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL))
          + init_key[j] + j; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=N-1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL))
          - i; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
    }

    mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long mt_random(void)
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) { /* generate N words at one time */
        int kk;

        if (mti == N+1)   /* if mt_srandom() has not been called, */
            mt_srandom(5489UL); /* a default initial seed is used */

        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }

    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

#undef N
#undef M
#undef MATRIX_A
#undef UPPER_MASK
#undef LOWER_MASK

/* End of "Mersenne Twister". */

/* End of McCode random number routine. */

/* randnorm: generate a random number from normal law */
double
randnorm(void)
{
  static double v1, v2, s;
  static int phase = 0;
  double X, u1, u2;

  if(phase == 0)
  {
    do
    {
      u1 = rand01();
      u2 = rand01();
      v1 = 2*u1 - 1;
      v2 = 2*u2 - 1;
      s = v1*v1 + v2*v2;
    } while(s >= 1 || s == 0);

    X = v1*sqrt(-2*log(s)/s);
  }
  else
  {
    X = v2*sqrt(-2*log(s)/s);
  }

  phase = 1 - phase;
  return X;
}

/**
 * Generate a random number from -1 to 1 with triangle distribution
 */
double randtriangle(void) {
	double randnum = rand01();
	if (randnum>0.5) return(1-sqrt(2*(randnum-0.5)));
	else return(sqrt(2*randnum)-1);
}

/**
 * Random number between 0.0 and 1.0 (including?)
 */
double rand01() {
	double randnum;
	randnum = (double) random();
	randnum /= (double) MC_RAND_MAX + 1;
	return randnum;
}

/**
 * Return a random number between 1 and -1
 */
double randpm1() {
	double randnum;
	randnum = (double) random();
	randnum /= ((double) MC_RAND_MAX + 1) / 2;
	randnum -= 1;
	return randnum;
}

/**
 * Return a random number between 0 and max.
 */
double rand0max(double max) {
	double randnum;
	randnum = (double) random();
	randnum /= ((double) MC_RAND_MAX + 1) / max;
	return randnum;
}

/**
 * Return a random number between min and max.
 */
double randminmax(double min, double max) {
	return rand0max(max - min) + max;
}

/* SECTION: main and signal handlers ======================================== */

/*******************************************************************************
* mchelp: displays instrument executable help with possible options
*******************************************************************************/
static void
mchelp(char *pgmname)
{
  int i;

  fprintf(stderr, "%s (%s) instrument simulation, generated with " MCCODE_STRING " (" MCCODE_DATE ")\n", mcinstrument_name, mcinstrument_source);
  fprintf(stderr, "Usage: %s [options] [parm=value ...]\n", pgmname);
  fprintf(stderr,
"Options are:\n"
"  -s SEED   --seed=SEED      Set random seed (must be != 0)\n"
"  -n COUNT  --ncount=COUNT   Set number of " MCCODE_PARTICLE "s to simulate.\n"
"  -d DIR    --dir=DIR        Put all data files in directory DIR.\n"
"  -t        --trace          Enable trace of " MCCODE_PARTICLE "s through instrument.\n"
"  -g        --gravitation    Enable gravitation for all trajectories.\n"
"  --no-output-files          Do not write any data files.\n"
"  -h        --help           Show this help message.\n"
"  -i        --info           Detailed instrument information.\n"
"  --format=FORMAT            Output data files using FORMAT="
   FLAVOR_UPPER
#ifdef USE_NEXUS
   " NEXUS"
#endif
"\n\n"
);
#ifdef USE_MPI
  fprintf(stderr,
  "This instrument has been compiled with MPI support.\n  Use 'mpirun %s [options] [parm=value ...]'.\n", pgmname);
#endif
  if(mcnumipar > 0)
  {
    fprintf(stderr, "Instrument parameters are:\n");
    for(i = 0; i < mcnumipar; i++)
      if (mcinputtable[i].val && strlen(mcinputtable[i].val))
        fprintf(stderr, "  %-16s(%s) [default='%s']\n", mcinputtable[i].name,
        (*mcinputtypes[mcinputtable[i].type].parminfo)(mcinputtable[i].name),
        mcinputtable[i].val);
      else
        fprintf(stderr, "  %-16s(%s)\n", mcinputtable[i].name,
        (*mcinputtypes[mcinputtable[i].type].parminfo)(mcinputtable[i].name));
  }

#ifndef NOSIGNALS
  fprintf(stderr, "Known signals are: "
#ifdef SIGUSR1
  "USR1 (status) "
#endif
#ifdef SIGUSR2
  "USR2 (save) "
#endif
#ifdef SIGBREAK
  "BREAK (save) "
#endif
#ifdef SIGTERM
  "TERM (save and exit)"
#endif
  "\n");
#endif /* !NOSIGNALS */
} /* mchelp */


/* mcshowhelp: show help and exit with 0 */
static void
mcshowhelp(char *pgmname)
{
  mchelp(pgmname);
  exit(0);
}

/* mcusage: display usage when error in input arguments and exit with 1 */
static void
mcusage(char *pgmname)
{
  fprintf(stderr, "Error: incorrect command line arguments\n");
  mchelp(pgmname);
  exit(1);
}

/* mcenabletrace: enable trace/mcdisplay or error if requires recompile */
static void
mcenabletrace(void)
{
 if(mctraceenabled)
  mcdotrace = 1;
 else
 {
   fprintf(stderr,
           "Error: trace not enabled (mcenabletrace)\n"
           "Please re-run the " MCCODE_NAME " compiler "
                   "with the --trace option, or rerun the\n"
           "C compiler with the MC_TRACE_ENABLED macro defined.\n");
   exit(1);
 }
}

/*******************************************************************************
* mcreadparams: request parameters from the prompt (or use default)
*******************************************************************************/
void
mcreadparams(void)
{
  int i,j,status;
  static char buf[CHAR_BUF_LENGTH];
  char *p;
  int len;

  MPI_MASTER(printf("Instrument parameters for %s (%s)\n",
                    mcinstrument_name, mcinstrument_source));

  for(i = 0; mcinputtable[i].name != 0; i++)
  {
    do
    {
      MPI_MASTER(
                 if (mcinputtable[i].val && strlen(mcinputtable[i].val))
                   printf("Set value of instrument parameter %s (%s) [default='%s']:\n",
                          mcinputtable[i].name,
                          (*mcinputtypes[mcinputtable[i].type].parminfo)
                          (mcinputtable[i].name), mcinputtable[i].val);
                 else
                   printf("Set value of instrument parameter %s (%s):\n",
                          mcinputtable[i].name,
                          (*mcinputtypes[mcinputtable[i].type].parminfo)
                          (mcinputtable[i].name));
                 fflush(stdout);
                 );
#ifdef USE_MPI
      if(mpi_node_rank == mpi_node_root)
        {
          p = fgets(buf, CHAR_BUF_LENGTH, stdin);
          if(p == NULL)
            {
              fprintf(stderr, "Error: empty input for paramater %s (mcreadparams)\n", mcinputtable[i].name);
              exit(1);
            }
        }
      else
        p = buf;
      MPI_Bcast(buf, CHAR_BUF_LENGTH, MPI_CHAR, mpi_node_root, MPI_COMM_WORLD);
#else /* !USE_MPI */
      p = fgets(buf, CHAR_BUF_LENGTH, stdin);
      if(p == NULL)
        {
          fprintf(stderr, "Error: empty input for paramater %s (mcreadparams)\n", mcinputtable[i].name);
          exit(1);
        }
#endif /* USE_MPI */
      len = strlen(buf);
      if (!len || (len == 1 && (buf[0] == '\n' || buf[0] == '\r')))
      {
        if (mcinputtable[i].val && strlen(mcinputtable[i].val)) {
          strncpy(buf, mcinputtable[i].val, CHAR_BUF_LENGTH);  /* use default value */
          len = strlen(buf);
        }
      }
      for(j = 0; j < 2; j++)
      {
        if(len > 0 && (buf[len - 1] == '\n' || buf[len - 1] == '\r'))
        {
          len--;
          buf[len] = '\0';
        }
      }

      status = (*mcinputtypes[mcinputtable[i].type].getparm)
                   (buf, mcinputtable[i].par);
      if(!status)
      {
        (*mcinputtypes[mcinputtable[i].type].error)(mcinputtable[i].name, buf);
        if (!mcinputtable[i].val || strlen(mcinputtable[i].val)) {
          fprintf(stderr, "       Change %s default value in instrument definition.\n", mcinputtable[i].name);
          exit(1);
        }
      }
    } while(!status);
  }
} /* mcreadparams */

/*******************************************************************************
* mcparseoptions: parse command line arguments (options, parameters)
*******************************************************************************/
void
mcparseoptions(int argc, char *argv[])
{
  int i, j;
  char *p;
  int paramset = 0, *paramsetarray;
  char *usedir=NULL;

  /* Add one to mcnumipar to avoid allocating zero size memory block. */
  paramsetarray = (int*)malloc((mcnumipar + 1)*sizeof(*paramsetarray));
  if(paramsetarray == NULL)
  {
    fprintf(stderr, "Error: insufficient memory (mcparseoptions)\n");
    exit(1);
  }
  for(j = 0; j < mcnumipar; j++)
    {
      paramsetarray[j] = 0;
      if (mcinputtable[j].val != NULL && strlen(mcinputtable[j].val))
      {
        int  status;
        char buf[CHAR_BUF_LENGTH];
        strncpy(buf, mcinputtable[j].val, CHAR_BUF_LENGTH);
        status = (*mcinputtypes[mcinputtable[j].type].getparm)
                   (buf, mcinputtable[j].par);
        if(!status) fprintf(stderr, "Invalid '%s' default value %s in instrument definition (mcparseoptions)\n", mcinputtable[j].name, buf);
        else paramsetarray[j] = 1;
      } else {
        (*mcinputtypes[mcinputtable[j].type].getparm)
          (NULL, mcinputtable[j].par);
        paramsetarray[j] = 0;
      }
    }
  for(i = 1; i < argc; i++)
  {
    if(!strcmp("-s", argv[i]) && (i + 1) < argc)
      mcsetseed(argv[++i]);
    else if(!strncmp("-s", argv[i], 2))
      mcsetseed(&argv[i][2]);
    else if(!strcmp("--seed", argv[i]) && (i + 1) < argc)
      mcsetseed(argv[++i]);
    else if(!strncmp("--seed=", argv[i], 7))
      mcsetseed(&argv[i][7]);
    else if(!strcmp("-n", argv[i]) && (i + 1) < argc)
      mcsetn_arg(argv[++i]);
    else if(!strncmp("-n", argv[i], 2))
      mcsetn_arg(&argv[i][2]);
    else if(!strcmp("--ncount", argv[i]) && (i + 1) < argc)
      mcsetn_arg(argv[++i]);
    else if(!strncmp("--ncount=", argv[i], 9))
      mcsetn_arg(&argv[i][9]);
    else if(!strcmp("-d", argv[i]) && (i + 1) < argc)
      usedir=argv[++i];  /* will create directory after parsing all arguments (end of this function) */
    else if(!strncmp("-d", argv[i], 2))
      usedir=&argv[i][2];
    else if(!strcmp("--dir", argv[i]) && (i + 1) < argc)
      usedir=argv[++i];
    else if(!strncmp("--dir=", argv[i], 6))
      usedir=&argv[i][6];
    else if(!strcmp("-h", argv[i]))
      mcshowhelp(argv[0]);
    else if(!strcmp("--help", argv[i]))
      mcshowhelp(argv[0]);
    else if(!strcmp("-i", argv[i])) {
      mcformat=FLAVOR_UPPER;
      mcinfo();
    }
    else if(!strcmp("--info", argv[i]))
      mcinfo();
    else if(!strcmp("-t", argv[i]))
      mcenabletrace();
    else if(!strcmp("--trace", argv[i]))
      mcenabletrace();
    else if(!strcmp("--gravitation", argv[i]))
      mcgravitation = 1;
    else if(!strcmp("-g", argv[i]))
      mcgravitation = 1;
    else if(!strncmp("--format=", argv[i], 9)) {
      mcformat=&argv[i][9];
    }
    else if(!strcmp("--format", argv[i]) && (i + 1) < argc) {
      mcformat=argv[++i];
    }
    else if(!strcmp("--no-output-files", argv[i]))
      mcdisable_output_files = 1;
    else if(argv[i][0] != '-' && (p = strchr(argv[i], '=')) != NULL)
    {
      *p++ = '\0';

      for(j = 0; j < mcnumipar; j++)
        if(!strcmp(mcinputtable[j].name, argv[i]))
        {
          int status;
          status = (*mcinputtypes[mcinputtable[j].type].getparm)(p,
                        mcinputtable[j].par);
          if(!status || !strlen(p))
          {
            (*mcinputtypes[mcinputtable[j].type].error)
              (mcinputtable[j].name, p);
            exit(1);
          }
          paramsetarray[j] = 1;
          paramset = 1;
          break;
        }
      if(j == mcnumipar)
      {                                /* Unrecognized parameter name */
        fprintf(stderr, "Error: unrecognized parameter %s (mcparseoptions)\n", argv[i]);
        exit(1);
      }
    }
    else if(argv[i][0] == '-') {
      fprintf(stderr, "Error: unrecognized option argument %s (mcparseoptions). Ignored.\n", argv[i++]);
    }
    else {
      fprintf(stderr, "Error: unrecognized argument %s (mcparseoptions). Aborting.\n", argv[i]);
      mcusage(argv[0]);
    }
  }
  if(!paramset)
    mcreadparams();                /* Prompt for parameters if not specified. */
  else
  {
    for(j = 0; j < mcnumipar; j++)
      if(!paramsetarray[j])
      {
        fprintf(stderr, "Error: Instrument parameter %s left unset (mcparseoptions)\n",
                mcinputtable[j].name);
        exit(1);
      }
  }
  free(paramsetarray);
#ifdef USE_MPI
  if (mcdotrace) mpi_node_count=1; /* disable threading when in trace mode */
#endif
  if (usedir && strlen(usedir) && !mcdisable_output_files) mcuse_dir(usedir);
} /* mcparseoptions */

#ifndef NOSIGNALS
mcstatic char  mcsig_message[256];


/*******************************************************************************
* sighandler: signal handler that makes simulation stop, and save results
*******************************************************************************/
void sighandler(int sig)
{
  /* MOD: E. Farhi, Sep 20th 2001: give more info */
  time_t t1, t0;
#define SIG_SAVE 0
#define SIG_TERM 1
#define SIG_STAT 2
#define SIG_ABRT 3

  printf("\n# " MCCODE_STRING ": [pid %i] Signal %i detected", getpid(), sig);
#ifdef USE_MPI
  printf(" [proc %i]", mpi_node_rank);
#endif
#if defined(SIGUSR1) && defined(SIGUSR2) && defined(SIGKILL)
  if (!strcmp(mcsig_message, "sighandler") && (sig != SIGUSR1) && (sig != SIGUSR2))
  {
    printf("\n# Fatal : unrecoverable loop ! Suicide (naughty boy).\n");
    kill(0, SIGKILL); /* kill myself if error occurs within sighandler: loops */
  }
#endif
  switch (sig) {
#ifdef SIGINT
    case SIGINT : printf(" SIGINT (interrupt from terminal, Ctrl-C)"); sig = SIG_TERM; break;
#endif
#ifdef SIGILL
    case SIGILL  : printf(" SIGILL (Illegal instruction)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGFPE
    case SIGFPE  : printf(" SIGFPE (Math Error)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGSEGV
    case SIGSEGV : printf(" SIGSEGV (Mem Error)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGTERM
    case SIGTERM : printf(" SIGTERM (Termination)"); sig = SIG_TERM; break;
#endif
#ifdef SIGABRT
    case SIGABRT : printf(" SIGABRT (Abort)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGQUIT
    case SIGQUIT : printf(" SIGQUIT (Quit from terminal)"); sig = SIG_TERM; break;
#endif
#ifdef SIGTRAP
    case SIGTRAP : printf(" SIGTRAP (Trace trap)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGPIPE
    case SIGPIPE : printf(" SIGPIPE (Broken pipe)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGUSR1
    case SIGUSR1 : printf(" SIGUSR1 (Display info)"); sig = SIG_STAT; break;
#endif
#ifdef SIGUSR2
    case SIGUSR2 : printf(" SIGUSR2 (Save simulation)"); sig = SIG_SAVE; break;
#endif
#ifdef SIGHUP
    case SIGHUP  : printf(" SIGHUP (Hangup/update)"); sig = SIG_SAVE; break;
#endif
#ifdef SIGBUS
    case SIGBUS  : printf(" SIGBUS (Bus error)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGURG
    case SIGURG  : printf(" SIGURG (Urgent socket condition)"); sig = SIG_ABRT; break;
#endif
#ifdef SIGBREAK
    case SIGBREAK: printf(" SIGBREAK (Break signal, Ctrl-Break)"); sig = SIG_SAVE; break;
#endif
    default : printf(" (look at signal list for signification)"); sig = SIG_ABRT; break;
  }
  printf("\n");
  printf("# Simulation: %s (%s) \n", mcinstrument_name, mcinstrument_source);
  printf("# Breakpoint: %s ", mcsig_message);
  if (strstr(mcsig_message, "Save") && (sig == SIG_SAVE))
    sig = SIG_STAT;
  SIG_MESSAGE("sighandler");
  if (mcget_ncount() == 0)
    printf("(0 %%)\n" );
  else
  {
    printf("%.2f %% (%10.1f/%10.1f)\n", 100.0*mcget_run_num()/mcget_ncount(), 1.0*mcget_run_num(), 1.0*mcget_ncount());
  }
  t0 = (time_t)mcstartdate;
  t1 = time(NULL);
  printf("# Date:      %s", ctime(&t1));
  printf("# Started:   %s", ctime(&t0));

  if (sig == SIG_STAT)
  {
    printf("# " MCCODE_STRING ": Resuming simulation (continue)\n");
    fflush(stdout);
    return;
  }
  else
  if (sig == SIG_SAVE)
  {
    printf("# " MCCODE_STRING ": Saving data and resume simulation (continue)\n");
    mcsave(NULL);
    fflush(stdout);
    return;
  }
  else
  if (sig == SIG_TERM)
  {
    printf("# " MCCODE_STRING ": Finishing simulation (save results and exit)\n");
    mcfinally();
    exit(0);
  }
  else
  {
    fflush(stdout);
    perror("# Last I/O Error");
    printf("# " MCCODE_STRING ": Simulation stop (abort).\n");
// This portion of the signal handling only works on UNIX
#if defined(__unix__) || defined(__APPLE__)
    signal(sig, SIG_DFL); /* force to use default sighandler now */
    kill(getpid(), sig);  /* and trigger it with the current signal */
#endif
    exit(-1);
  }
#undef SIG_SAVE
#undef SIG_TERM
#undef SIG_STAT
#undef SIG_ABRT

} /* sighandler */
#endif /* !NOSIGNALS */

/*******************************************************************************
* mccode_main: McCode main() function.
*******************************************************************************/
int mccode_main(int argc, char *argv[])
{
/*  double run_num = 0; */
  time_t  t;
#ifdef USE_MPI
  char mpi_node_name[MPI_MAX_PROCESSOR_NAME];
  int  mpi_node_name_len;
#endif /* USE_MPI */

#ifdef MAC
  argc = ccommand(&argv);
#endif

#ifdef USE_MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_node_count); /* get number of nodes */
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_node_rank);
  MPI_Comm_set_name(MPI_COMM_WORLD, mcinstrument_name);
  MPI_Get_processor_name(mpi_node_name, &mpi_node_name_len);
#endif /* USE_MPI */

t = time(NULL);
mcseed = (long)t+(long)getpid();

#ifdef USE_MPI
/* *** print number of nodes *********************************************** */
  if (mpi_node_count > 1) {
    MPI_MASTER(
    printf("Simulation '%s' (%s): running on %i nodes (master is '%s', MPI version %i.%i).\n",
      mcinstrument_name, mcinstrument_source, mpi_node_count, mpi_node_name, MPI_VERSION, MPI_SUBVERSION);
    );
  }
#endif /* USE_MPI */
  
  mcstartdate = (long)t;  /* set start date before parsing options and creating sim file */

/* *** parse options ******************************************************* */
  SIG_MESSAGE("main (Start)");
  mcformat=getenv(FLAVOR_UPPER "_FORMAT") ?
           getenv(FLAVOR_UPPER "_FORMAT") : FLAVOR_UPPER;
  mcinstrument_exe = argv[0]; /* store the executable path */
  /* read simulation parameters and options */
  mcparseoptions(argc, argv); /* sets output dir and format */
  
#ifdef USE_MPI
  if (mpi_node_count > 1) {
    /* share the same seed, then adapt random seed for each node */
    MPI_Bcast(&mcseed, 1, MPI_LONG, 0, MPI_COMM_WORLD); /* root sends its seed to slaves */
    mcseed += mpi_node_rank; /* make sure we use different seeds per node */
  }
#endif
  srandom(mcseed);

/* *** install sig handler, but only once !! after parameters parsing ******* */
#ifndef NOSIGNALS
#ifdef SIGQUIT
  if (signal( SIGQUIT ,sighandler) == SIG_IGN)
    signal( SIGQUIT,SIG_IGN);   /* quit (ASCII FS) */
#endif
#ifdef SIGABRT
  if (signal( SIGABRT ,sighandler) == SIG_IGN)
    signal( SIGABRT,SIG_IGN);   /* used by abort, replace SIGIOT in the future */
#endif
#ifdef SIGTERM
  if (signal( SIGTERM ,sighandler) == SIG_IGN)
    signal( SIGTERM,SIG_IGN);   /* software termination signal from kill */
#endif
#ifdef SIGUSR1
  if (signal( SIGUSR1 ,sighandler) == SIG_IGN)
    signal( SIGUSR1,SIG_IGN);   /* display simulation status */
#endif
#ifdef SIGUSR2
  if (signal( SIGUSR2 ,sighandler) == SIG_IGN)
    signal( SIGUSR2,SIG_IGN);
#endif
#ifdef SIGHUP
  if (signal( SIGHUP ,sighandler) == SIG_IGN)
    signal( SIGHUP,SIG_IGN);
#endif
#ifdef SIGILL
  if (signal( SIGILL ,sighandler) == SIG_IGN)
    signal( SIGILL,SIG_IGN);    /* illegal instruction (not reset when caught) */
#endif
#ifdef SIGFPE
  if (signal( SIGFPE ,sighandler) == SIG_IGN)
    signal( SIGSEGV,SIG_IGN);    /* floating point exception */
#endif
#ifdef SIGBUS
  if (signal( SIGBUS ,sighandler) == SIG_IGN)
    signal( SIGSEGV,SIG_IGN);    /* bus error */
#endif
#ifdef SIGSEGV
  if (signal( SIGSEGV ,sighandler) == SIG_IGN)
    signal( SIGSEGV,SIG_IGN);   /* segmentation violation */
#endif
#endif /* !NOSIGNALS */
  mcsiminfo_init(NULL); /* open SIM */
  SIG_MESSAGE("main (Init)");
  mcinit();
#ifndef NOSIGNALS
#ifdef SIGINT
  if (signal( SIGINT ,sighandler) == SIG_IGN)
    signal( SIGINT,SIG_IGN);    /* interrupt (rubout) only after INIT */
#endif
#endif /* !NOSIGNALS */

/* ================ main particle generation/propagation loop ================ */
#if defined (USE_MPI)
  /* sliced Ncount on each MPI node */
  mcncount = mpi_node_count > 1 ?
    floor(mcncount / mpi_node_count) :
    mcncount; /* number of rays per node */
#endif

/* main particle event loop */
while(mcrun_num < mcncount || mcrun_num < mcget_ncount())
  {
#ifndef NEUTRONICS
    mcgenstate();
#endif
    /* old init: mcsetstate(0, 0, 0, 0, 0, 1, 0, sx=0, sy=1, sz=0, 1); */
    mcraytrace();
    mcrun_num++;
  }

#ifdef USE_MPI
 /* merge run_num from MPI nodes */
  if (mpi_node_count > 1) {
  double mcrun_num_double = (double)mcrun_num;
  mc_MPI_Sum(&mcrun_num_double, 1);
  mcrun_num = (unsigned long long)mcrun_num_double;
  }
#endif

/* save/finally executed by master node/thread */
  mcfinally();

#ifdef USE_MPI
  MPI_Finalize();
#endif /* USE_MPI */

  return 0;
} /* mccode_main */

#ifdef NEUTRONICS
/*Main neutronics function steers the McStas calls, initializes parameters etc */
/* Only called in case NEUTRONICS = TRUE */
void neutronics_main_(float *inx, float *iny, float *inz, float *invx, float *invy, float *invz, float *intime, float *insx, float *insy, float *insz, float *inw, float *outx, float *outy, float *outz, float *outvx, float *outvy, float *outvz, float *outtime, float *outsx, float *outsy, float *outsz, float *outwgt)
{

  extern double mcnx, mcny, mcnz, mcnvx, mcnvy, mcnvz;
  extern double mcnt, mcnsx, mcnsy, mcnsz, mcnp;

  /* External code governs iteration - McStas is iterated once per call to neutronics_main. I.e. below counter must be initiancated for each call to neutronics_main*/
  mcrun_num=0;

  time_t t;
  t = (time_t)mcstartdate;
  mcstartdate = t;  /* set start date before parsing options and creating sim file */
  mcinit();

  /* *** parse options *** */
  SIG_MESSAGE("main (Start)");
  mcformat=getenv(FLAVOR_UPPER "_FORMAT") ?
           getenv(FLAVOR_UPPER "_FORMAT") : FLAVOR_UPPER;

  /* Set neutron state based on input from neutronics code */
  mcsetstate(*inx,*iny,*inz,*invx,*invy,*invz,*intime,*insx,*insy,*insz,*inw);

  /* main neutron event loop - runs only one iteration */

  //mcstas_raytrace(&mcncount); /* prior to McStas 1.12 */

  mcallowbackprop = 1; //avoid absorbtion from negative dt
  int argc=1;
  char *argv[0];
  int dummy = mccode_main(argc, argv);

  *outx =  mcnx;
  *outy =  mcny;
  *outz =  mcnz;
  *outvx =  mcnvx;
  *outvy =  mcnvy;
  *outvz =  mcnvz;
  *outtime =  mcnt;
  *outsx =  mcnsx;
  *outsy =  mcnsy;
  *outsz =  mcnsz;
  *outwgt =  mcnp;

  return;
} /* neutronics_main */

#endif /*NEUTRONICS*/

#endif /* !MCCODE_H */
/* End of file "mccode-r.c". */
/* End of file "mccode-r.c". */

#line 4977 "./reverse_test.c"

#line 1 "mcstas-r.c"
/*******************************************************************************
*
* McStas, neutron ray-tracing package
*         Copyright (C) 1997-2009, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Runtime: share/mcstas-r.c
*
* %Identification
* Written by: KN
* Date:    Aug 29, 1997
* Release: McStas X.Y
* Version: $Revision$
*
* Runtime system for McStas.
* Embedded within instrument in runtime mode.
*
* Usage: Automatically embbeded in the c code whenever required.
*
* $Id$
*
*******************************************************************************/

#ifndef MCSTAS_R_H
#include "mcstas-r.h"
#endif
#ifdef DANSE
#include "mcstas-globals.h"
#endif

/*******************************************************************************
* The I/O format definitions and functions
*******************************************************************************/

/*the magnet stack*/
#ifdef MC_POL_COMPAT
void (*mcMagnetPrecession) (double, double, double, double, double, double,
    double, double*, double*, double*, double, Coords, Rotation)=NULL;
Coords   mcMagnetPos;
Rotation mcMagnetRot;
double*  mcMagnetData                = NULL;
/* mcMagneticField(x, y, z, t, Bx, By, Bz) */
int (*mcMagneticField) (double, double, double, double,
    double*, double*, double*, void *) = NULL;
#endif

#ifndef MCSTAS_H

/*******************************************************************************
* mcstore_neutron: stores neutron coodinates into global array (per component)
*******************************************************************************/
void
mcstore_neutron(MCNUM *s, int index, double x, double y, double z,
               double vx, double vy, double vz, double t,
               double sx, double sy, double sz, double p)
{
    double *dptr = &s[11*index];
    *dptr++  = x;
    *dptr++  = y ;
    *dptr++  = z ;
    *dptr++  = vx;
    *dptr++  = vy;
    *dptr++  = vz;
    *dptr++  = t ;
    *dptr++  = sx;
    *dptr++  = sy;
    *dptr++  = sz;
    *dptr    = p ;
} /* mcstore_neutron */

/*******************************************************************************
* mcrestore_neutron: restores neutron coodinates from global array
*******************************************************************************/
void
mcrestore_neutron(MCNUM *s, int index, double *x, double *y, double *z,
               double *vx, double *vy, double *vz, double *t,
               double *sx, double *sy, double *sz, double *p)
{
    double *dptr = &s[11*index];
    *x  =  *dptr++;
    *y  =  *dptr++;
    *z  =  *dptr++;
    *vx =  *dptr++;
    *vy =  *dptr++;
    *vz =  *dptr++;
    *t  =  *dptr++;
    *sx =  *dptr++;
    *sy =  *dptr++;
    *sz =  *dptr++;
    *p  =  *dptr;
} /* mcrestore_neutron */

/*******************************************************************************
* mcsetstate: transfer parameters into global McStas variables 
*******************************************************************************/
void
mcsetstate(double x, double y, double z, double vx, double vy, double vz,
           double t, double sx, double sy, double sz, double p)
{
  extern double mcnx, mcny, mcnz, mcnvx, mcnvy, mcnvz;
  extern double mcnt, mcnsx, mcnsy, mcnsz, mcnp;

  mcnx = x;
  mcny = y;
  mcnz = z;
  mcnvx = vx;
  mcnvy = vy;
  mcnvz = vz;
  mcnt = t;
  mcnsx = sx;
  mcnsy = sy;
  mcnsz = sz;
  mcnp = p;
} /* mcsetstate */

/*******************************************************************************
* mcgenstate: set default neutron parameters 
*******************************************************************************/
void
mcgenstate(void)
{
  mcsetstate(0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1);
  /* old initialisation: mcsetstate(0, 0, 0, 0, 0, 1, 0, sx=0, sy=1, sz=0, 1); */
}

/* intersection routines ==================================================== */

/*******************************************************************************
* inside_rectangle: Check if (x,y) is inside rectangle (xwidth, yheight) 
* return 0 if outside and 1 if inside 
*******************************************************************************/
int inside_rectangle(double x, double y, double xwidth, double yheight)
{
  if (x>-xwidth/2 && x<xwidth/2 && y>-yheight/2 && y<yheight/2)
    return 1;
  else
    return 0;
}

/*******************************************************************************
 * box_intersect: compute time intersection with a box
 * returns 0 when no intersection is found
 *      or 1 in case of intersection with resulting times dt_in and dt_out
 * This function written by Stine Nyborg, 1999. 
 *******************************************************************************/
int box_intersect(double *dt_in, double *dt_out,
                  double x, double y, double z,
                  double vx, double vy, double vz,
                  double dx, double dy, double dz)
{
  double x_in, y_in, z_in, tt, t[6], a, b;
  int i, count, s;

      /* Calculate intersection time for each of the six box surface planes
       *  If the box surface plane is not hit, the result is zero.*/

  if(vx != 0)
   {
    tt = -(dx/2 + x)/vx;
    y_in = y + tt*vy;
    z_in = z + tt*vz;
    if( y_in > -dy/2 && y_in < dy/2 && z_in > -dz/2 && z_in < dz/2)
      t[0] = tt;
    else
      t[0] = 0;

    tt = (dx/2 - x)/vx;
    y_in = y + tt*vy;
    z_in = z + tt*vz;
    if( y_in > -dy/2 && y_in < dy/2 && z_in > -dz/2 && z_in < dz/2)
      t[1] = tt;
    else
      t[1] = 0;
   }
  else
    t[0] = t[1] = 0;

  if(vy != 0)
   {
    tt = -(dy/2 + y)/vy;
    x_in = x + tt*vx;
    z_in = z + tt*vz;
    if( x_in > -dx/2 && x_in < dx/2 && z_in > -dz/2 && z_in < dz/2)
      t[2] = tt;
    else
      t[2] = 0;

    tt = (dy/2 - y)/vy;
    x_in = x + tt*vx;
    z_in = z + tt*vz;
    if( x_in > -dx/2 && x_in < dx/2 && z_in > -dz/2 && z_in < dz/2)
      t[3] = tt;
    else
      t[3] = 0;
   }
  else
    t[2] = t[3] = 0;

  if(vz != 0)
   {
    tt = -(dz/2 + z)/vz;
    x_in = x + tt*vx;
    y_in = y + tt*vy;
    if( x_in > -dx/2 && x_in < dx/2 && y_in > -dy/2 && y_in < dy/2)
      t[4] = tt;
    else
      t[4] = 0;

    tt = (dz/2 - z)/vz;
    x_in = x + tt*vx;
    y_in = y + tt*vy;
    if( x_in > -dx/2 && x_in < dx/2 && y_in > -dy/2 && y_in < dy/2)
      t[5] = tt;
    else
      t[5] = 0;
   }
  else
    t[4] = t[5] = 0;

  /* The intersection is evaluated and *dt_in and *dt_out are assigned */

  a = b = s = 0;
  count = 0;

  for( i = 0; i < 6; i = i + 1 )
    if( t[i] == 0 )
      s = s+1;
    else if( count == 0 )
    {
      a = t[i];
      count = 1;
    }
    else
    {
      b = t[i];
      count = 2;
    }

  if ( a == 0 && b == 0 )
    return 0;
  else if( a < b )
  {
    *dt_in = a;
    *dt_out = b;
    return 1;
  }
  else
  {
    *dt_in = b;
    *dt_out = a;
    return 1;
  }

} /* box_intersect */

/*******************************************************************************
 * cylinder_intersect: compute intersection with a cylinder
 * returns 0 when no intersection is found
 *      or 2/4/8/16 bits depending on intersection,
 *     and resulting times t0 and t1
 * Written by: EM,NB,ABA 4.2.98 
  *******************************************************************************/
int
cylinder_intersect(double *t0, double *t1, double x, double y, double z,
                   double vx, double vy, double vz, double r, double h)
{
  double D, t_in, t_out, y_in, y_out;
  int ret=1;

  D = (2*vx*x + 2*vz*z)*(2*vx*x + 2*vz*z)
    - 4*(vx*vx + vz*vz)*(x*x + z*z - r*r);

  if (D>=0)
  {
    if (vz*vz + vx*vx) {
      t_in  = (-(2*vz*z + 2*vx*x) - sqrt(D))/(2*(vz*vz + vx*vx));
      t_out = (-(2*vz*z + 2*vx*x) + sqrt(D))/(2*(vz*vz + vx*vx));
    } else if (vy) { /* trajectory parallel to cylinder axis */
      t_in = (-h/2-y)/vy;
      t_out = (h/2-y)/vy;
      if (t_in>t_out){
        double tmp=t_in;
        t_in=t_out;t_out=tmp;
      }
    } else return 0;
    y_in = vy*t_in + y;
    y_out =vy*t_out + y;

    if ( (y_in > h/2 && y_out > h/2) || (y_in < -h/2 && y_out < -h/2) )
      return 0;
    else
    {
      if (y_in > h/2)
        { t_in = ((h/2)-y)/vy; ret += 2; }
      else if (y_in < -h/2)
        { t_in = ((-h/2)-y)/vy; ret += 4; }
      if (y_out > h/2)
        { t_out = ((h/2)-y)/vy; ret += 8; }
      else if (y_out < -h/2)
        { t_out = ((-h/2)-y)/vy; ret += 16; }
    }
    *t0 = t_in;
    *t1 = t_out;
    return ret;
  }
  else
  {
    *t0 = *t1 = 0;
    return 0;
  }
} /* cylinder_intersect */


/*******************************************************************************
 * sphere_intersect: Calculate intersection between a line and a sphere.
 * returns 0 when no intersection is found
 *      or 1 in case of intersection with resulting times t0 and t1 
 *******************************************************************************/
int
sphere_intersect(double *t0, double *t1, double x, double y, double z,
                 double vx, double vy, double vz, double r)
{
  double A, B, C, D, v;

  v = sqrt(vx*vx + vy*vy + vz*vz);
  A = v*v;
  B = 2*(x*vx + y*vy + z*vz);
  C = x*x + y*y + z*z - r*r;
  D = B*B - 4*A*C;
  if(D < 0)
    return 0;
  D = sqrt(D);
  *t0 = (-B - D) / (2*A);
  *t1 = (-B + D) / (2*A);
  return 1;
} /* sphere_intersect */

/*******************************************************************************
 * plane_intersect: Calculate intersection between a plane and a line.
 * returns 0 when no intersection is found (i.e. line is parallel to the plane)
 * returns 1 or -1 when intersection time is positive and negative respectively
 *******************************************************************************/
int
plane_intersect(double *t, double x, double y, double z,
                 double vx, double vy, double vz, double nx, double ny, double nz, double wx, double wy, double wz)
{
  double s;
  if (fabs(s=scalar_prod(nx,ny,nz,vx,vy,vz))<FLT_EPSILON) return 0;
  *t = - scalar_prod(nx,ny,nz,x-wx,y-wy,z-wz)/s;
  if (*t<0) return -1;
  else return 1;
} /* plane_intersect */

#endif /* !MCSTAS_H */
/* End of file "mcstas-r.c". */

#line 5337 "./reverse_test.c"
#ifdef MC_TRACE_ENABLED
int mctraceenabled = 1;
#else
int mctraceenabled = 0;
#endif
#define MCSTAS "/usr/share/mcstas/2.7/tools/Python/mcrun/../mccodelib/../../../"
int mcdefaultmain = 1;
char mcinstrument_name[] = "logspir_test";
char mcinstrument_source[] = "reverse_test.instr";
char *mcinstrument_exe=NULL; /* will be set to argv[0] in main */
int main(int argc, char *argv[]){return mccode_main(argc, argv);}
void mcinit(void);
void mcraytrace(void);
void mcsave(FILE *);
void mcfinally(void);
void mcdisplay(void);

/* Shared user declarations for all components 'PSD_monitor'. */
#line 57 "/usr/share/mcstas/2.7/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"

#ifndef ARRAYS_H
#define ARRAYS_H
typedef double* DArray1d;
DArray1d create_darr1d(int n);
void destroy_darr1d(DArray1d a);

typedef double** DArray2d;
DArray2d create_darr2d(int nx, int ny);
void destroy_darr2d(DArray2d a);

typedef double*** DArray3d;
DArray3d create_darr3d(int nx, int ny, int nz);
void destroy_darr3d(DArray3d a);
#endif
#ifndef ARRAYS_C
#define ARRAYS_C
#include <stdlib.h>

DArray1d create_darr1d(int n){
  DArray1d arr2d;
  arr2d = calloc(n, sizeof(double));
  return arr2d;
}
void destroy_darr1d(DArray1d a){
  free(a);
}

DArray2d create_darr2d(int nx, int ny){
  DArray2d arr2d;
  arr2d = calloc(nx, sizeof(double *));

  double *p1;
  p1 = calloc(nx*ny, sizeof(double));

  int i;
  for (i=0; i<nx; i++){
    arr2d[i] = &(p1[i*ny]);
  }
  return arr2d;
}
void destroy_darr2d(DArray2d a){
  free(a[0]);
  free(a);
}

DArray3d create_darr3d(int nx, int ny, int nz){
  DArray3d arr3d;
  int i, j;

  // 1d
  arr3d = calloc(nx, sizeof(double **));

  // d2
  double **p1;
  p1 = calloc(nx*ny, sizeof(double *));

  for (i=0; i<nx; i++){
    arr3d[i] = &(p1[i*ny]);
  }

  // 3d
  double *p2;
  p2 = calloc(nx*ny*nz, sizeof(double));
  for (i=0; i<nx; i++){
    for (j=0; j<ny; j++){
      arr3d[i][j] = &(p2[(i*ny+j)*nz]);
    }
  }
  return arr3d;
}
void destroy_darr3d(DArray3d a){
  free(a[0][0]);
  free(a[0]);
  free(a);
}
#endif

#line 5435 "./reverse_test.c"

/* Shared user declarations for all components 'FlatEllipse_finite_mirror'. */
#line 54 "FlatEllipse_finite_mirror.comp"
/**
\mainpage
Simple Meta-Conic Neutron Raytracer is a framework for raytracing geometries of the form: @f$ r^2=k_1 + k_2 z + k_3 z^2 @f$.

<h3>General Notes</h3>
To use the software you must make a Scene element using the function makeScene(). You must then add items to this scene element using the various add function (addDisk(), addParaboloid(), etc...). Next you must call the function traceSingleNeutron() for every neutron you would like to trace through the geometry. The maximum number of each geometry you can place in a scene is defined by the MAX_CONICSURF, MAX_DISK and MAX_DETECTOR definitions in the conic.h file.

<h3>TODO</h3>

@todo
       Name variable for each component <br/>
       Normalize Detector Events by weight of neutron <br/>

<h3>Known Bugs</h3>
@bug  HPPH works incorrectly for M != 1 <br/>
      Neutrons t=1e-11 away from a surface will pass through <br/>

<h3>Note on Pointers</h3>
This framework uses pointers extensivly. It is important to be familiar with how to use pointers.
<br/>Here is a quick overview.

@code
//Making an Element
ConicSurf c = makeParaboloid(...);
int k = 10;

//Making a Pointer from an Element
ConicSurf * c_pointer = &c;
int * k_pointer = &k;

//Making an Element from a Pointer
ConicSurf c2 = *c_pointer;
int ten = *k_pointer;

//Getting Item in Element
double k1 = c.k1;

//Getting Item in Element from a Pointer
double k1 = c_pointer->k1;

//Functions that have pointers as parameters can modify the element being pointed to.
//Functions that have elements as parameters can not modify the value of the element.
@endcode

<h3>Stand Alone Example Code</h3>
This framework can be used entirely by itself as the following example demonstrates.
@code
/////////////////////////////////////////////////////////////////////////////////
// Giacomo Resta <gresta@mit.edu>
//
// Basic standalone exmaple of a raytracer. It is advised to modify
// the getRandom() function to use a better random number generator (i.e. glib).
// The systems random generator was used to preserve clarity and conciseness.
/////////////////////////////////////////////////////////////////////////////////

#include <math.h>
#include <stdio.h>

#include "conic.h"
#include "w1_general.h"

#define NUM_NEUTRON 1000000
*/
//Function to get random number
double getRandom() {
    return (double)lrand48()/RAND_MAX;
}
/*
//Function to get new particle from source
Particle generateParticleFromSource(double radius, double div, double Lambda0,
        double num_neutrons) {

    double chi = 2*M_PI*(getRandom());
    double r = sqrt(getRandom()) * radius;

    double v = 3956.036 / Lambda0;

    double theta = sqrt(getRandom())*div;
    double phi = getRandom()*2*M_PI;
    double tan_r = tan(theta);

    double vz = v/sqrt(1+tan_r*tan_r);
    double vr = tan_r * vz;

    return makeParticle(r*cos(chi),r*sin(chi),0,cos(phi)*vr,sin(phi)*vr,
                vz,0,0.0,0.0,0.0,1.0/num_neutrons);
}

//Function to add items to scene
void addItems(Scene* s, double instr_len,double r,double f,double M,
        double max_mir_len,double m, double mirr_thick) {

    //Change code here to make different Scenes
    PP p = addPPShell(0.0, instr_len, r, f, M, max_mir_len, m, m, mirr_thick, s);
    addDisk(p.p0->zs, 0.0, rConic(p.p0->ze, *p.p0)*p.p0->zs/p.p0->ze, s);
    addDisk(p.p1->zs, rConic(p.p1->zs, *p.p1), 10000,s);
    addDetector(10.0, -0.01, 0.01, -0.01, 0.01, 600, 600, NUM_NEUTRON, "test.txt", s);
}

//Main Function
int main() {
    //seed random generator (starts the generator)
    srand48((long)time(0));

    //Make a new Scene
    Scene s = makeScene();

    //Add Items and initialize Scene
    addItems(&s,10.0,0.068,4.2,1,0.7,3,0.001);
    initSimulation(&s);

    //Raytrace all particles through the Scene
    double i;
    for (i = 0; i < NUM_NEUTRON; i++) {
        Particle p = generateParticleFromSource(0.005, 0.02422, 4, NUM_NEUTRON);
        traceSingleNeutron(&p,s);
    }

    //Finish Simulation of the Scene
    finishSimulation(&s);

    return 0;
}
@endcode

*/

/**
    @file conic.h
    \brief General Framework for generating and raytracing geometries of the
    form @f$ r = k_1+k_2 z+k_3 z^2 @f$

    @author Giacomo Resta <gresta@mit.edu>
    @version 0.2

    @section LICENSE
    Permission is hereby granted, free of charge, to any person obtaining a
    copy of this software and associated documentation files (the "Software"),
    to deal in the Software without restriction, including without limitation
    the rights to use, copy, modify, merge, publish, distribute, sublicense,
    and/or sell copies of the Software, and to permit persons to whom the Software
    is furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
    INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
    PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE
    FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
    OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
    DEALINGS IN THE SOFTWARE.

    @section DESCRIPTION
     General Framework for generating and raytracing geometries of the form
     @f$ r = k_1+k_2 z+k_3 z^2 @f$
*/

/////////////////////////////////////
// Simulation
/////////////////////////////////////

#ifndef MIT_CONICS
#define MIT_CONICS

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/** @defgroup simgroup Simulator Internals
    Contains items general to the simulation
    @{
*/
//! Max number of ConicSurf allowed in a Scene
#define MAX_FLATSURF 200

//! Max number of ConicSurf allowed in a Scene
#define MAX_CONICSURF 100

//! Max number of Disks allowed in a Scene
#define MAX_DISK 100

//! Max number of Detectors allowed in a Scene
#define MAX_DETECTOR 10

//! If "1" simulator will record z location where neutron with greatest grazing angle reflected for each ConicSurf
/*! The information is stored in the max_ga and max_ga_z0 members of each ConicSurf, which are only present if
    this flag is 1. See source code for clarification.

@note You must use default traceNeutronConic function or implement saving routine yourself */
#define REC_MAX_GA 0

#define V2Q_conic 1.58825361e-3
#define Q2V_conic 629.622368
#define POT_V 46.839498800356665//theta_crit = 4 pi^2 kappa^2/m_n^2
#define m_Si 0.478 // m-value of pure silicon
//! Stucture to represent a point
typedef struct {
    double x; //!< x-axis position of point
    double y; //!< y-axis position of point
    double z; //!< z-axis position of point
} Point;

//! Structure to represent a vector
typedef struct {
    double x; //!< x-axis length of vector
    double y; //!< y-axis length of vector
    double z; //!< z-axis length of vector
} Vec;

//! Structure to represent a particle
typedef struct {
    double _x; //!< x axis position of particle
    double _y; //!< y axis position of particle
    double _z; //!< z axis position of particle
    double _vx; //!< x axis components of velocity
    double _vy; //!< y axis components of velocity
    double _vz; //!< z axis components of velocity
    double _sx; //!< x spin vector components
    double _sy; //!< y spin vector components
    double _sz; //!< z spin vector components
    double w; //!< Weight of particle
    int silicon; //!< +1 if Particle is in silicon -1 if Particle is in air
    int absorb; //!< Absorb flag (0 is not absorbed)
    double _t; //!< Time of flight of particle
} Particle;

/*! \brief Function to make a point

@param x x-axis position
@param y y-axis position
@param z z-axis position
*/
Point makePoint(double x, double y, double z) {
    Point p;
    p.x = x;
    p.y = y;
    p.z = z;

    return p;
}

/*! \brief Function to make a vector

@param x x-axis length
@param y y-axis length
@param z z-axis length
*/
Vec makeVec(double x, double y, double z) {
    Vec p;
    p.x = x;
    p.y = y;
    p.z = z;

    return p;
}

//! Function to compute length of a vector
double getMagVec(Vec v) {
    return sqrt(v.x*v.x+v.y*v.y+v.z*v.z);
}

//! Function to compute dot product of two vectors
double dotVec(Vec v1, Vec v2) {
    return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}

//! Function to compute the sum of two vectors
Vec difVec(Vec v1, Vec v2){
    return makeVec(v1.x-v2.x, v1.y-v2.y, v1.z-v2.z);
}

//! Function to compute the sum of two vectors
Vec sumVec(Vec v1, Vec v2){
    return makeVec(v1.x+v2.x, v1.y+v2.y, v1.z+v2.z);
}

//! Function to compute the sum of two vectors
Vec skalarVec(Vec v1, double k){
    return makeVec(v1.x*k, v1.y*k, v1.z*k);
}
/*! \brief Function to make a particle

@param x x-axis position
@param y y-axis position
@param z z-axis position
@param vx x-axis velocity
@param vy y-axis velocity
@param vz z-axis velocity
@param t time of flight of neutron
@param sx x-axis component of spin vector
@param sy y-axis component of spin vector
@param sz z-axis component of spin vector
@param w weight of particle
*/
Particle makeParticle(double x, double y, double z,
    double vx, double vy, double vz, double t,
    double sx, double sy, double sz, int silicon, double w) {

    Particle pa;

    pa._x = x;
    pa._y = y;
    pa._z = z;

    pa._vx = vx;
    pa._vy = vy;
    pa._vz = vz;

    pa._sx = sx;
    pa._sy = sy;
    pa._sz = sz;
    pa.silicon = silicon;
    pa.w = w;

    pa.absorb = 0;
    pa._t = t;

    return pa;
}

//! Function to get position of particle
Point getParticlePos(Particle p) {
    return makePoint(p._x,p._y,p._z);
}

//! Function to get velocity vector of particle
Vec getParticleVel(Particle p) {
    return makeVec(p._vx, p._vy, p._vz);
}

/*! \brief Function to move particle a specific time step.

Will not move particle if t < 0. Does not simulate
gravity.

@param t time step to move particle
@param p pointer to particle
*/
void moveParticleT(double t, Particle* p) {
    if (t < 0)
        return;
    if (p->silicon>0){
        p->_x = p->_x+p->_vx*t;
        p->_y = p->_y+p->_vy*t;
        p->_z = p->_z+p->_vz*t;
        p->_t = p->_t+t;
        p->w *= exp(-t*98900/52.338);//ugly hard coding in the penetration depth of silicon at 989 m/s
        //absorbParticle(p);
    }
    else{
        p->_x = p->_x+p->_vx*t;
        p->_y = p->_y+p->_vy*t;
        p->_z = p->_z+p->_vz*t;
        p->_t = p->_t+t;
    }
}

/*! \brief Function to move particle to position z.

Will not move particle if moving particle to z
position would require negative time.
Does not simulate gravity.

@param z z-position to move particle
@param p pointer to particle
*/
void moveParticleZ(double z, Particle* p) {
    double t = (z-p->_z)/p->_vz;
    moveParticleT(t, p);
}

/*! \brief Function to compute new position of particle
without modifying the position of the actual particle.

Will not move particle if t < 0. Does not simulate gravity.

@param t timestep to move particle
@param p particle
*/
Particle copyMoveParticleT(double t, Particle p) {
    Particle p2 = p;
    moveParticleT(t,&p2);
    return p2;
}

/*! \brief Function to move particle to position z
without modifying the position of the actual particle.

Will not move particle if moving particle to z
position would require negative time.
Does not simulate gravity.

@param z z-position to move particle
@param p pointer to particle
*/
Particle copyMoveParticleZ(double z, Particle p) {
    Particle p2 = p;
    moveParticleZ(z,&p2);
    return p2;
}

/*! \brief Mathematical Aid for Snell's Law for reflection.

Does not take into account grazing angle constrictions.
Only computes mathematical reflection.

@param n Normal vector
@param p Pointer to particle
*/
void reflectParticle(Vec n, Particle* p) {
    double vn = dotVec(getParticleVel(*p),n);

    p->_vx = p->_vx-2*vn*n.x;
    p->_vy = p->_vy-2*vn*n.y;
    p->_vz = p->_vz-2*vn*n.z;
}

/*! \brief Function to mark particle as absorbed

@param p Pointer to particle to be absorbed */
void absorbParticle(Particle* p)  {
    p->_vx = 0;
    p->_vy = 0;
    p->_vz = 0;
    p->w = 0;
    p->absorb = 1;
}

/*! \brief Function to set weight of particle.

Will set the weight of the particle to w.  */
void setWeightParticle(double w, Particle* pa) {
    pa->w = w;
}

/*! \brief Function to solve quadratic equations for smallest positive result.

If no positive result returns -1. Parameters are coefficents such that
@f$ 0=A z^2 + B z + C @f$

@return Will return smallest positive value or -1 if no smallest positive value
*/
double solveQuad(double A, double B, double C) {
    if (fabs(A) < 1e-11 && B != 0)           //FIXME: 1e-11 cutoff may cause problems
        return -C/B;
    else {
        double det = B*B - 4*A*C;
        if (det < 0)
            return -1;
        else {
            double sdet = sqrt(det);
            double s1 = (-B+sdet)/(2*A);
            double s2 = (-B-sdet)/(2*A);

            if (fabs(s1) < 1e-11) s1=0.0;     //FIXME: 1e-11 cutoff may cause problems
            if (fabs(s2) < 1e-11) s2=0.0;     //FIXME: 1e-11 cutoff may cause problems

            if (s1 > 0.0) {
                if (s2 > 0.0) {
                    if (s1 > s2)
                        return s2;
                    return s1;
                }
                else
                    return s1;
           }
           if (s2 > 0.0)
               return s2;
        }
    }
    return -1;
}

//! Returns sign of x, either -1 or 1
int sign(double x) {
    if (x < 0)
        return -1;
    return 1;
}

/** @} */ //end of simgroup

/*! @ingroup detectorgroup
\brief Structure for representing inline detectors

@warning Do not directly modify this structure*/
typedef struct {
    double z0; //!< z-axis position of detector
    double xmin; //!< Smallest x value to detect
    double xmax; //!< Largest x value to detect
    double xstep; //!< x size of subsampling pixel
    double ymin; //!< Smallest y value to detect
    double ymax; //!< Largest y value to detect
    double ystep; //!< y size of subsampling pixel
    int nrows;    //!< Number of pixels along y axis
    int ncols;    //!< Number of pixels along x axis
    double num_particles; //!< Number of particles being emitted from source
    double *num_count; //!< Pointer to the number of particles that hit detector
    double **data; //!< Pointer to 2d data, pixel_x_y = data[x][y]
    char* filename; //!< Name of output file of detector (should end in .txt)
} Detector;

/*! @ingroup diskgroup
\brief Structure for representing Disk geometry

 Creates a doughnut with inner radius r0 and outer radus r1 at position z0.
Neutrons between r0 and r1 are absorbed. */
typedef struct {
    double r0; //!< Inner radius of doughnut
    double r1; //!< Outer radius of doughnut
    double z0; //!< z-axis position of Disk
} Disk;

/*! @ingroup conicgroup */
enum ConicType {
    PARA,
    HYPER,
    ELLIP
};


/*! @ingroup conicgroup
\brief Structure to contain z-axis symetric conic sections

Contains any geometry that can be expressed as
@f$ r^2=k_1 + k_2 z + k_3 z^2 @f$

@warning Do not directly modify values in this structure directly */
typedef struct {
    double k1; //!< @f$ k_1 @f$ in equation below
    double k2; //!< @f$ k_2 @f$ in equation below
    double k3; //!< @f$ k_3 @f$ in equation below
    double zs; //!< z-axis position of start of mirror
    double ze; //!< z-axis position of end of mirror
    double m;  //!< m value for mirror (1.0 for Nickel)
    int doubleReflections; //!< 0 if reflections from the back of the surface cannot happen 1 otherwise
    //Only for reference
    double f1; //!< z-axis position of first focus
    double f2; //!< z-axis position of second focus, for paraboloid this is unassigned
    double a;  //!< Value of a, specific to geometry type
    double c;  //!< Value of c, for paraboloid this is unassigned
    enum ConicType type; //!< Type of mirror geometry

    #if REC_MAX_GA
    double max_ga; //!< Max Grazing Angle of Reflected Neutron (Exists only if REC_MAX_GA)
    double max_ga_z0; //!< Collision point of Max Grazing Neutron (Exists only if REC_MAX_GA)
    #endif

} ConicSurf;


/*! @ingroup flatgroup
\brief Structure to contain z-axis symetric flat sections

Contains flat geometries which can be expressed as
@f$ x^2 = k_1 + k_2 z + k_3 z^2 @f$ or
@f$ y^2 = k_1 + k_2 z + k_3 z^2 @f$

@warning Do not directly modify values in this structure directly */
typedef struct {
    double k1; //!< @f$ k_1 @f$ in equation below
    double k2; //!< @f$ k_2 @f$ in equation below
    double k3; //!< @f$ k_3 @f$ in equation below
    double zs; //!< z-axis position of start of mirror
    double ze; //!< z-axis position of end of mirror
    double ll; //!< left/lower limit of mirror along translational symmetry
    double rl; //!< right/upper limit of mirror along translational symmetry
    double m;  //!< m value for mirror (1.0 for Nickel)

    //Only for reference
    double f1; //!< z-axis position of first focus
    double f2; //!< z-axis position of second focus, for paraboloid this is unassigned
    double a;  //!< Value of a, specific to geometry type
    double c;  //!< Value of c, for paraboloid this is unassigned
    //enum FlatType type; //!< Type of mirror geometry
    int doubleReflections; // will determine whether the geometry allows double reflections
    #if REC_MAX_GA
    double max_ga; //!< Max Grazing Angle of Reflected Neutron (Exists only if REC_MAX_GA)
    double max_ga_z0; //!< Collision point of Max Grazing Neutron (Exists only if REC_MAX_GA)
    #endif

} FlatSurf;
/*! @ingroup simgroup
\brief Structure to hold all scene geometry

The number of possible ConicSurf, Disk and Detector in the Scene are
determined by MAX_CONICSURF, MAX_DISK and MAX_DETECTOR.
*/
typedef struct {
    FlatSurf f[MAX_FLATSURF]; //!< Array of all ConicSurf in Scene
    int num_f;                  //!< Number of ConicSurf in Scene

    ConicSurf c[MAX_CONICSURF]; //!< Array of all ConicSurf in Scene
    int num_c;                  //!< Number of ConicSurf in Scene

    Disk di[MAX_DISK];          //!< Array of all Disk in Scene
    int num_di;                 //!< Number of Disk in Scene

    Detector d[MAX_DETECTOR];  //!< Array of all Detector in Scene
    int num_d;                 //!< Number of Detector in Scene

    //! Function called to handle Neutron-Flat Interaction
    void (*traceNeutronFlat)(Particle*,FlatSurf);
    //! Function called to handle Neutron-Conic Interaction
    void (*traceNeutronConic)(Particle*,ConicSurf);
    //! Function called to handle Neutron-Disk Interaction
    void (*traceNeutronDisk)(Particle*,Disk);
    //! Function called to handle Neutron-Detector Interaction
    void (*traceNeutronDetector)(Particle*,Detector);

} Scene;

/////////////////////////////////////
// Inline Detector
/////////////////////////////////////

/** @defgroup detectorgroup Detector
    Contains code related to the inline detectors
    @{
*/

/*! \brief Function to make Detector

@param z0 z-axis position of detector
@param xmin Smallest x value to detect
@param xmax Largest x value to detect
@param ymin Smallest y value to detect
@param ymax Largest y value to detect
@param xres Number of pixels along x axis
@param yres Number of pixels along y axis
@param num_particles Total number of particles being emitted
@param filename Name of output file of detector (should end in .txt)
*/
Detector makeDetector(double z0,double xmin, double xmax, double ymin, double ymax, int xres,
    int yres, double num_particles, char* filename) {

    Detector d;
    d.z0 = z0;
    d.xmin = xmin;
    d.xmax = xmax;
    d.xstep = (xmax-xmin)/xres;

    d.ymin = ymin;
    d.ymax = ymax;
    d.ystep = (ymax-ymin)/yres;

    d.ncols = xres;
    d.nrows = yres;
    d.filename = filename;
    d.num_particles = num_particles;

    d.num_count = (double*)malloc(sizeof(double));
    if (d.num_count == NULL) {
        fprintf(stderr, "MEMORY ALLOCATION PROBLEM\n");
        exit(-1);
    }
    (*d.num_count) = 0;

    d.data = (double**)malloc(d.ncols*sizeof(double *));
    if (d.data == NULL) {
        fprintf(stderr, "MEMORY ALLOCATION PROBLEM\n");
        exit(-1);
    }
    int x;
    for(x = 0; x  < d.ncols; x++) {
        d.data[x] = (double*)malloc(d.ncols*sizeof(double));
        if (d.data[x] == NULL) {
            fprintf(stderr, "MEMORY ALLOCATION PROBLEM\n");
            exit(-1);
        }
        (*d.data[x]) = 0;
    }

    return d;
}

/*! \brief Function to make and add Detector

@param z0 z-axis position of detector
@param xmin Smallest x value to detect
@param xmax Largest x value to detect
@param ymin Smallest y value to detect
@param ymax Largest y value to detect
@param xres Number of pixels along x axis
@param yres Number of pixels along y axis
@param num_particles Total number of particles being emitted
@param filename Name of output file of detector (should end in .txt)
@param s Scene to add Detector to
*/
Detector* addDetector(double z0, double xmin, double xmax, double ymin, double ymax, double xres,
    double yres, double num_particles, char* filename, Scene* s) {
    if (s->num_d >= MAX_DETECTOR-1) {
        fprintf(stderr,"TOO MANY DETECTORS IN SCENE");
        exit(-1);
    }
    s->d[s->num_d] = makeDetector(z0,xmin,xmax,ymin,ymax,xres,yres,num_particles,filename);
    s->num_d++;
    return &s->d[s->num_d-1];
}

/*! \brief Function to compute time of first collision for a Detector.

@param p Particle to consider
@param d Detector to consider

@return Time until the propogation or -1 if particle will not hit detector
*/
double getTimeOfFirstCollisionDetector(Particle p, Detector d) {
    double t = (d.z0-p._z)/p._vz;
    if (t <= 0)
        return -1;
    Particle p2 = copyMoveParticleT(t,p);
    if (p2._x > d.xmax || p2._x < d.xmin || p2._y > d.ymax || p2._y < d.ymin)
        return -1;
    return t;
}

/*! \brief Function to raytrace Detector

@param p Pointer to particle to be traced
@param d Detector to be traced
*/
void traceNeutronDetector(Particle* p, Detector d) {
    double t = getTimeOfFirstCollisionDetector(*p, d);
    if (t < 0)
        return;
    moveParticleT(t,p);
    d.data[(int)floor((p->_x-d.xmin)/d.xstep)][(int)floor((p->_y-d.ymin)/d.ystep)] += p->w;
    (*d.num_count) += p->w;
}

/*! \brief Function to finalize detector

Will write data and free data array.

@param d Detector to finalize
*/
void finishDetector(Detector d) {
    int x,y;
    if (d.filename != "") {
        FILE *file;
        file = fopen(d.filename,"w");

        double intensity = (*d.num_count);
        fprintf(file, "#I=%e I_ERR=%e xmin=%f xmax=%f ymin=%f ymax=%f ncols=%i nrows=%i\n",
            intensity, sqrt(intensity/d.num_particles), d.xmin, d.xmax, d.ymin, d.ymax, d.ncols, d.nrows); //FIXME: check I_ERR sqrt(I/num_particles)

        //Write data
        for (x=0; x < d.ncols; x++) {
            for (y=0; y < d.nrows; y++)
                fprintf(file, "%e ", d.data[x][y]);
            fprintf(file, "\n");
        }
        fclose(file);
    }
    for (x=0; x < d.ncols; x++)
        free(d.data[x]);
    free(d.data);
    free(d.num_count);
}

/** @} */ //end of detectorgroup

/////////////////////////////////////
// Geometry Types
/////////////////////////////////////

/////////////////////////////////////
// Disks
/////////////////////////////////////

/** @defgroup diskgroup Disk
    Contains code related to Disks
    @{
*/

/*! \brief Function for creating a Disk structure

@param z0 z-axis position of Disk
@param r0 Inner radius of doughnut
@param r1 Outer radius of doughnut

@see Disk
*/
Disk makeDisk(double z0, double r0, double r1) {
    Disk d;

    d.r0 = r0;
    d.z0 = z0;
    d.r1 = r1;

    return d;
}

/*! \brief Function for making and adding Disk to Scene

@param z0 z-axis position of Disk
@param r0 Inner radius of doughnut
@param r1 Outer radius of doughnut
@param s Scene to add Disk to

@see Disk
*/
Disk* addDisk(double z0, double r0, double r1, Scene* s) {
    if (s->num_di >= MAX_DISK-1) {
        fprintf(stderr,"TOO MANY DISKS IN SCENE");
        exit(-1);
    }
    s->di[s->num_di] = makeDisk(z0, r0, r1);
    s->num_di++;
    return &s->di[s->num_di -1];
}

/*! \brief Function to compute time of first collision for a disk

@param p Particle to consider
@param d Disk to consider
@return Time until the propogation or -1 if particle will not hit disk
*/
double getTimeOfFirstCollisionDisk(Particle p, Disk d) {
    double tz = (d.z0-p._z)/p._vz;
    if (tz <= 0)
        return -1;
    Particle p2 = copyMoveParticleT(tz, p);
    double rp = sqrt(p2._x*p2._x+p2._y*p2._y);
    if (rp > d.r0 && rp < d.r1 && fabs(p2._z-d.z0) < 1e-11)
        return (d.z0-p._z)/p._vz;
    return -1;
}

/*! \brief Function to raytrace Disks

@param p Pointer to particle to be traced
@param d Disk to be traced
*/
void traceNeutronDisk(Particle* p, Disk d) {
    double t = getTimeOfFirstCollisionDisk(*p, d);

    if (t <= 0)
        return;

    moveParticleT(t, p);
    //absorbParticle(p); //Disk will only be used to propagate neutrons somewhere
}

/** @} */ //end of diskgroup

/////////////////////////////////////
// Z-Axis Symetric Conic Sections
/////////////////////////////////////

/** @defgroup conicgroup ConicSurf
    Contains code related to ConicSurfs
    @{
*/

/*! \brief Function to return radius of ConicSurf at a z-axis position.

Will return radius even if z is outside the bounds of zs and ze
for the particular ConicSurf.

@param z z-axis position to compute radius
@param s ConicSurf to compute radius of
*/
double rConic(double z, ConicSurf s) {
    return sqrt(s.k1+s.k2*z+s.k3*z*z);
}

/*! \brief Function for generating Hyperboloid ConicSurf.

@param f1 z position of focus closest to actual mirror surface
@param f2 z position of focus furthest from actual mirror surface
@param p A Point on the actual surface of the mirror
@param zstart z position of start of mirror surface
@param zend z position of end of mirror surface
@param m m value for reflectivity of the surface

@see ConicSurf
*/
ConicSurf makeHyperboloid(double f1, double f2, Point p,
   double zstart, double zend, double m) {
    ConicSurf s;
    s.zs = zstart;
    s.ze = zend;

    double r2 = p.x*p.x+p.y*p.y;
    double c = (f1-f2)/2;

    double u = p.z+c-f1;
    double a = sqrt(((u*u+c*c+r2)-sqrt(pow(u*u+c*c+r2,2)-4*c*c*u*u))/2);

    s.k3 = c*c/(a*a)-1;
    s.k2 = 2*s.k3*(c-f1);
    s.k1 = (s.k3)*(c-f1)*(c-f1)-c*c+a*a;

    s.m = m;
    s.f1 = f1;
    s.f2 = f2;
    s.a = a;
    s.c = c;

    s.type = HYPER;

    #if REC_MAX_GA
    s.max_ga = -1;
    s.max_ga_z0 = -1;
    #endif

    return s;
}

/*! \brief Function for generating Ellipsoid ConicSurf.

@param f1 z position of focus closest to actual mirror surface
@param f2 z position of focus furthest from actual mirror surface
@param p A Point on the actual surface of the mirror
@param zstart z position of start of mirror surface
@param zend z position of end of mirror surface
@param m m value for reflectivity of the surface

@see ConicSurf
*/
ConicSurf makeEllipsoid(double f1, double f2, Point p,
    double zstart, double zend, double m, int doubleReflections) {
    ConicSurf s;
    s.zs = zstart;
    s.ze = zend;
    s.doubleReflections = doubleReflections;
    double r2 = p.x*p.x+p.y*p.y;
    double c = (f1-f2)/2;

    double u = p.z+c-f1;
    double a = sqrt(((u*u+c*c+r2)+sqrt(pow(u*u+c*c+r2,2)-4*c*c*u*u))/2);

    s.k3 = c*c/(a*a)-1;
    s.k2 = 2*s.k3*(c-f1);
    s.k1 = (s.k3)*(c-f1)*(c-f1)-c*c+a*a;

    s.m = m;
    s.f1 = f1;
    s.f2 = f2;
    s.a = a;
    s.c = c;

    s.type = ELLIP;

    #if REC_MAX_GA
    s.max_ga = -1;
    s.max_ga_z0 = -1;
    #endif

    return s;
}


/*! \brief Function for generating Flat Ellipse with symmetry along the vertical y.

@param f1 z position of focus closest to actual mirror surface
@param f2 z position of focus furthest from actual mirror surface
@param b the short half axis of the ellipse, positive for translational symmetry along y, negative for translational symmetry along x
@param p A Point on the actual surface of the mirror
@param zstart z position of start of mirror surface
@param zend z position of end of mirror surface
@param ll the left/lower limit of the ellipse surface (in either y or x for b > 0 or b < 0
@param rl the right/upper limit of the ellipse surface (in either y or x for b > 0 or b < 0
@param m m value for reflectivity of the surface


@see ConicSurf
*/
FlatSurf makeFlatEllipse(
    double f1,
    double f2,
    Point p,
    double zstart,
    double zend,
    double ll,
    double rl,
    double m,
    int doubleReflections) {
        FlatSurf s;
        s.zs = zstart;
        s.ze = zend;
        s.doubleReflections = doubleReflections;
        s.ll = ll;
        s.rl = rl;
        double r2 = p.x*p.x + p.y*p.y;
        double c = (f1-f2)/2;
        double u = p.z+c-f1;
        double a = sqrt(((u*u+c*c+r2)+sqrt(pow(u*u+c*c+r2,2)-4*c*c*u*u))/2);

        s.k3 = c*c/(a*a)-1;
        s.k2 = 2*s.k3*(c-f1);
        s.k1 = (s.k3)*(c-f1)*(c-f1)-c*c+a*a;

        s.m = m;
        s.f1 = f1;
        s.f2 = f2;
        s.a = a;
        s.c = c;

        //s.type = ELLIP;

        #if REC_MAX_GA
        s.max_ga = -1;
        s.max_ga_z0 = -1;
        #endif

        return s;
}

/*! \brief Function for generating Paraboloid ConicSurf.

@param f z position of focus closest to actual mirror surface
@param p A Point on the actual surface of the mirror
@param zstart z position of start of mirror surface
@param zend z position of end of mirror surface
@param m m value for reflectivity of the surface

@see ConicSurf
*/
ConicSurf makeParaboloid(double f, Point p, double zstart,
    double zend, double m, int doubleReflections) {

    ConicSurf s;
    s.zs = zstart;
    s.ze = zend;
    s.doubleReflections = doubleReflections;
    double r2 = p.x*p.x+p.y*p.y;
    double a = (-(p.z-f)+sign(p.z-f)*sqrt((p.z-f)*(p.z-f)+r2))/2;

    s.k3 = 0.0;
    s.k2 = 4*a;
    s.k1 = s.k2*(a-f);

    s.m = m;
    s.f1 = f;
    s.a = a;

    s.type = PARA;

    #if REC_MAX_GA
    s.max_ga = -1;
    s.max_ga_z0 = -1;
    #endif

    return s;
}

/*! \brief Function for generating Flat Parabola for FlatSurf.

@param f z position of focus closest to actual mirror surface
@param p A Point on the actual surface of the mirror, putting one of x or y to 0 results in the surface being parallel to said coordinate
@param zstart z position of start of mirror surface
@param zend z position of end of mirror surface
@param ll the left/lower limit of the ellipse surface (in either y or x for b > 0 or b < 0
@param rl the right/upper limit of the ellipse surface (in either y or x for b > 0 or b < 0
@param m m value for reflectivity of the surface
@param doubleReflections wether double reflections are allowed

@see FlatSurf
*/
FlatSurf makeFlatparbola(
    double f,
    Point p,
    double zstart,
    double zend,
    double ll,
    double rl,
    double m,
    int doubleReflections) {

    FlatSurf s;
    s.zs = zstart;
    s.ze = zend;
    s.doubleReflections = doubleReflections;
    double r2 = p.x*p.x+p.y*p.y;
    double a = (-(p.z-f)+sign(p.z-f)*sqrt((p.z-f)*(p.z-f)+r2))/2;

    s.k3 = 0.0;
    s.k2 = 4*a;
    s.k1 = s.k2*(a-f);

    s.m = m;
    s.f1 = f;
    s.a = a;
    s.ll = ll;
    s.rl = rl;
    //s.type = PARA;

    #if REC_MAX_GA
    s.max_ga = -1;
    s.max_ga_z0 = -1;
    #endif

    return s;
}

/*! \brief Function for generating and adding Paraboloid ConicSurf.

@param f1 z position of focus closest to actual mirror surface
@param p A Point on the actual surface of the mirror
@param zstart z position of start of mirror surface
@param zend z position of end of mirror surface
@param m m value for reflectivity of the surface
@param s Scene to add Paraboloid to

@see ConicSurf
*/
ConicSurf* addParaboloid(double f1, Point p, double zstart, double zend,
    double m, int doubleReflections, Scene* s) {
    if (s->num_c >= MAX_CONICSURF-1) {
        fprintf(stderr,"TOO MANY CONICSURF IN SCENE");
        exit(-1);
    }
    s->c[s->num_c] = makeParaboloid(f1,p,zstart,zend,m, doubleReflections);
    s->num_c++;
    return &s->c[s->num_c-1];
}

/*! \brief Function for generating and adding a flat Parabolic FlatSurf.

@param f1 z position of focus closest to actual mirror surface
@param p A Point on the actual surface of the mirror
@param zstart z position of start of mirror surface
@param zend z position of end of mirror surface
@param ll the lower bound of the mirror
@param rl the upper bound of the mirror
@param m m value for reflectivity of the surface
@param s Scene to add Ellipsoid to
@param doubleReflections wether double reflections can occur
@see FlatSurf
*/
FlatSurf* addFlatParabola(
    double f,
    Point p,
    double zstart,
    double zend,
    double ll,
    double rl,
    double m,
    Scene* s,
    int doubleReflections) {
    if (s->num_f >= MAX_FLATSURF-1) {
        fprintf(stderr,"TOO MANY FLATSURF IN SCENE");
        exit(-1);
    }
    s->f[s->num_f] = makeFlatparbola(f,p,zstart,zend,ll,rl,m,doubleReflections);
    s->num_f++;
    return &s->f[s->num_f-1];
}

/*! \brief Function for generating and adding Hyperboloid ConicSurf.

@param f1 z position of focus closest to actual mirror surface
@param f2 z position of focus furthest from actual mirror surface
@param p A Point on the actual surface of the mirror
@param zstart z position of start of mirror surface
@param zend z position of end of mirror surface
@param m m value for reflectivity of the surface
@param s Scene to add Hyperboloid to

@see ConicSurf
*/
ConicSurf* addHyperboloid(double f1, double f2, Point p, double zstart,
    double zend, double m, Scene* s) {
    if (s->num_c >= MAX_CONICSURF-1) {
        fprintf(stderr,"TOO MANY CONICSURF IN SCENE");
        exit(-1);
    }
    s->c[s->num_c] = makeHyperboloid(f1,f2,p,zstart,zend,m);
    s->num_c++;
    return &s->c[s->num_c-1];
}

/*! \brief Function for generating and adding Ellipsoid ConicSurf.

@param f1 z position of focus closest to actual mirror surface
@param f2 z position of focus furthest from actual mirror surface
@param p A Point on the actual surface of the mirror
@param zstart z position of start of mirror surface
@param zend z position of end of mirror surface
@param m m value for reflectivity of the surface
@param s Scene to add Ellipsoid to

@see ConicSurf
*/
ConicSurf* addEllipsoid(double f1, double f2, Point p, double zstart,
    double zend, double m, int doubleReflections, Scene* s) {
    if (s->num_c >= MAX_CONICSURF-1) {
        fprintf(stderr,"TOO MANY CONICSURF IN SCENE");
        exit(-1);
    }
    s->c[s->num_c] = makeEllipsoid(f1,f2,p,zstart,zend,m,doubleReflections);
    s->num_c++;
    return &s->c[s->num_c-1];
}

/*! \brief Function for generating and adding a flat Ellipse FlatSurf.

@param f1 z position of focus closest to actual mirror surface
@param f2 z position of focus furthest from actual mirror surface
@param p A Point on the actual surface of the mirror
@param zstart z position of start of mirror surface
@param zend z position of end of mirror surface
@param m m value for reflectivity of the surface
@param s Scene to add Ellipsoid to

@see ConicSurf
*/
FlatSurf* addFlatEllipse(
    double f1,
    double f2,
    Point p,
    double zstart,
    double zend,
    double ll,
    double rl,
    double m,
    int doubleReflections,
    Scene* s) {
    if (s->num_f >= MAX_FLATSURF-1) {
        fprintf(stderr,"TOO MANY FLATSURF IN SCENE");
        exit(-1);
    }
    s->f[s->num_f] = makeFlatEllipse(f1,f2,p,zstart,zend,ll,rl,m,doubleReflections);
    s->num_f++;
    return &s->f[s->num_f-1];
}
//!TODO
double getGrazeAngleConic(Particle p, ConicSurf s) {
    /*
    double v = sqrt(dotVec(getParticleVel(p),getParticleVel(p)));
    double vn = dotVec(getParticleVel(p),n);
    return fabs(acos(vn/v)) - M_PI/2;
    */
}

/*! \brief Function for returning normal vector of ConicSurf at Point p

Will compute vector even if p is not on surface.
MAKE SURE p IS ON SURFACE

@param p Point to compute normal vector
@param s ConicSurf to compute normal vector of
*/
Vec getNormConic(Point p, ConicSurf s) {
    double det = s.k2*s.k2+4*s.k3*(p.x*p.x+p.y*p.y-s.k1);
    if (det <= 0.){

        return makeVec(-p.x/sqrt(p.x*p.x + p.y*p.y),-p.y/(p.x*p.x + p.y*p.y),0);
    }
    double den = sqrt(det);
    double nx = -2*p.x/den;
    double ny = -2*p.y/den;
    double nz = sign(2*s.k3*p.z+s.k2);
    double n = sqrt(nx*nx+ny*ny+nz*nz);
    //printf("%f,%f,%f \n", nx/n,ny/n,nz/n);
    return makeVec(nx/n,ny/n,nz/n);
}

/*! \brief Function for returning normal vector of FlatSurf at Point p

Will compute vector even if p is not on surface.
MAKE SURE p IS ON SURFACE

@param p Point to compute normal vector
@param s FlatSurf to compute normal vector of; for s.b > 0 surface posseses translation symmetry along y direction
*/
Vec getNormFlat(Point p, FlatSurf s) {
    double r;
    //if(s.b > 0){
    r = p.x;
    //else{
    //    r = p.y;
    //};
    double den;
    double det = s.k2*s.k2+4*s.k3*(r*r-s.k1);
    if (det > 0){
        den = sqrt(det);
    }
    else{
    return makeVec(-1,0,0); // if the neutron hits the apex of the ellipse we run into a divide by zero problem
    }
    double nx = -2*p.x/den;
    double ny = 0;
    double nz = sign(2*s.k3*p.z+s.k2);
    double n = sqrt(nx*nx+ny*ny+nz*nz);
    return makeVec(nx/n,ny/n,nz/n);
}

/*! \brief Function to compute time of first collision for a ConicSurf

@param p Particle to consider
@param s ConicSurf to consider
@return Time until the propogation or -1 if particle will not hit surface
*/
double getTimeOfFirstCollisionConic(Particle p, ConicSurf s) {
    double tz = (s.zs-p._z)/p._vz;
    if (tz < 0) {
       tz = 0;
       if (p._z > s.ze)
            return -1;
    }

    Particle p2 = copyMoveParticleT(tz,p);

    double A = p2._vx*p2._vx+p2._vy*p2._vy-s.k3*p2._vz*p2._vz;
    double B = 2*(p2._vx*p2._x+p2._vy*p2._y-s.k3*p2._vz*p2._z)-s.k2*p2._vz;
    double C = p2._x*p2._x+p2._y*p2._y-s.k3*p2._z*p2._z-s.k2*p2._z-s.k1;

    double t = solveQuad(A,B,C);

    if (t <= 0 || p2._vz*t+p2._z > s.ze || p2._vz*t+p2._z < s.zs)
        return -1;
    return t+tz;
}

/*! \brief Function to compute time of first collision for a FlatSurf

@param p Particle to consider
@param s FlatSurf to consider
@return Time until the propogation or -1 if particle will not hit surface
*/
//TODO
double getTimeOfFirstCollisionFlat(Particle p, FlatSurf s) {
    double tz = (s.zs-p._z)/p._vz;
    if (tz < 0) {
       tz = 0;
       if (p._z > s.ze)
            return -1;
    }

    Particle p2 = copyMoveParticleT(tz,p);
    double vs = 0;//the vector important for calculating the intersection with the ellipse
    double s0 = 0;
    double vt = 0;//the other component only important for testing whether the mirror is hit
    double t0 = 0;
    //if(s.b > 0){
    vs = p2._vx;
    s0 = p2._x;
    vt = p2._vy;
    t0 = p2._y;

    //}
    /*else{
    vs = p2._vy;
    s0 = p2._y;
    vt = p2._vx;
    t0 = p2._x;
    };
    */
    double A = vs*vs-s.k3*p2._vz*p2._vz;
    double B = 2*(vs*s0-s.k3*p2._vz*p2._z)-s.k2*p2._vz;
    double C = s0*s0-s.k3*p2._z*p2._z-s.k2*p2._z-s.k1;

    double t = solveQuad(A,B,C);

    if (t <= 0 || p2._vz*t+p2._z > s.ze || p2._vz*t+p2._z < s.zs||vt*t+t0 < s.ll||vt*t +t0 > s.rl)
        return -1;
    return t+tz;
}

/*! \brief Function to handle supermirror reflectivity copied from mcstas.
@note Uses only m-value for calculating reflectivity curve TODO more sophisticated formulae in the future

@param q k_i - k_f momentum transfer of the neutron at the super mirror surface
@param m supermirror m-value
@param R_0 low angle reflectivity
@param Q_c critical momentum transfer of the super mirror

@return p weight reduction of the neutron for further simulation
*/
double calcSupermirrorReflectivity(double q, double m, double R_0, double Q_c){
    double arg;
    double beta = 0;
    double alpha = 2.5;
    double W = 0.004;
    double weight = 1.0; //neutron weight to be transformed
    q = fabs(q);
    if (m >= 10){
        weight = 1.0;
        return weight;
    }
    if (W==0 && alpha==0) {
      m=m*0.9853+0.1978;
      W=-0.0002*m+0.0022;
      alpha=0.2304*m+5.0944;
      beta=-7.6251*m+68.1137;
      if (m<=3) {
	    alpha=m;
	    beta=0;
        }
    }
    arg = W > 0 ? (q - m*Q_c)/W : 11;
    if (arg > 10 || m <= 0 || Q_c <=0 || R_0 <= 0) {
      weight = 0.0;
      return weight;
    }

    if (m < 1) { Q_c *= m; m=1; }

    if(q <= Q_c) {
      weight = R_0;
      return weight;
    }


    weight = R_0*0.5*(1 - tanh(arg))*(1 - alpha*(q - Q_c) + beta*(q - Q_c)*(q - Q_c));
    return weight;
}

/*! \brief Function to handle reflection of neutron for a ConicSurf.

@note Uses step function for reflectivity

@warning Make sure particle has been moved to surface of mirror
before computing reflection

@param p Pointer of particle to reflect
@param s ConicSurf to use

@return Value of critical angle of the neutron or -1 if neutron is absorbed

@see traceNeutronConic()
*/
double reflectNeutronConic(Particle* p, ConicSurf s) {//TODO add super mirror reflectivity an passing neutrons
    Vec n = getNormConic(getParticlePos(*p),s);
    Vec pv = getParticleVel(*p);
    double weight = 0;
    

    double v = getMagVec(pv);
    double vn = dotVec(pv,n);
    
    weight = calcSupermirrorReflectivity(V2Q_conic*vn*2, s.m, 1.0, 0.0218);

    //Hitting shell from outside

    if (vn > 0 && !s.doubleReflections) {
        absorbParticle(p);
        return -1;
    }



    double ga = fabs(acos(vn/v)) - M_PI/2;
    double gc = 6.84459399932*s.m/v;

    if (weight <= 0) {
        printf("weight <0");
        absorbParticle(p);
        return -1;
    }
    else {
        p->_vx = p->_vx-2*vn*n.x;
        p->_vy = p->_vy-2*vn*n.y;
        p->_vz = p->_vz-2*vn*n.z;
        p->w *= weight;
    }
    return ga;
}

/*! \brief Function to handle refraction of neutron.

@warning Make sure particle has been moved to surface of mirror
before computing reflection

@param p Pointer of particle to reflect
@param n Normal vector
@param m1 m-value of the original material we come from
@param m2 m-value of the goal material we go to

@return Value of critical angle of the neutron or -1 if neutron is absorbed

@see traceNeutronConic()
*/
void refractNeutronFlat(Particle* p, Vec n, double m1, double m2) {
    //Vec n = getNormFlat(getParticlePos(*p),s);
    Vec pv = getParticleVel(*p);
    //printf("incoming %.14g %.14g %.14g\n", pv.x, pv.y, pv.z);
    //printf("normal %.14g %.14g %.14g\n", n.x, n.y, n.z);
    double v = getMagVec(pv);
    double vn = dotVec(pv, n);
    //printf("vn %.9g", vn);
    Vec v_p = difVec(pv, skalarVec(n, vn));

    double k2_perp = POT_V*(m1*m1-m2*m2)+vn*vn;
    //printf("the magnitude %f\n", k2_perp);
    if (k2_perp>0){//refraction
        
        k2_perp = sqrt(k2_perp)*sign(vn);
        p->_vx = v_p.x + n.x*k2_perp;
        p->_vy = v_p.y + n.y*k2_perp;
        p->_vz = v_p.z + n.z*k2_perp;
        p-> silicon *= -1;// from silicon to air or vice versa
        //printf("resulting vector %f %f %f\n", p->_vx, p->_vy, p->_vz);
    }else{//total reflection, no change in material
        p->_vx = p->_vx-2*vn*n.x;
        p->_vy = p->_vy-2*vn*n.y;
        p->_vz = p->_vz-2*vn*n.z;
        //no change in material
    }

    
}

/*! \brief Function to handle reflection of neutron for a FlatSurf.

@note Uses step function for reflectivity

@warning Make sure particle has been moved to surface of mirror
before computing reflection

@param p Pointer of particle to reflect
@param s FlatSurf to use

@return Value of critical angle of the neutron or -1 if neutron is absorbed

@see traceNeutronConic()
*/
double reflectNeutronFlat(Particle* p, FlatSurf s) {
    Vec n = getNormFlat(getParticlePos(*p),s);
    Vec pv = getParticleVel(*p);
    //printf("nothing");
    double v = getMagVec(pv);
    double vn = dotVec(pv,n);
    //printf("before %f \n", p->w);
    //Hitting shell from outside For FlatSurface this has to be checked
    // make it able to reflect from the outside
    double ga = fabs(acos(vn/v)) - M_PI/2;
    double gc = 6.84459399932*s.m/v;


    double weight = 0;
    weight = calcSupermirrorReflectivity(V2Q_conic*2*vn, s.m, 0.995, 0.0218);

    if (vn > 0 && !s.doubleReflections) {
        absorbParticle(p);
        return -1;
    }


    if (weight < 0) {
        printf("this happens?");
        absorbParticle(p);
        return -1;
    } else {//here we need to implement the refraction
        if (getRandom() <= weight){//to be updated to use the mcstas random function or quasoi deterministic model
            //printf("oh a reflections\n");
            //printf("this total reflection?");
            p->_vx = p->_vx-2*vn*n.x;
            p->_vy = p->_vy-2*vn*n.y;
            p->_vz = p->_vz-2*vn*n.z;
            return ga;
        }
        else{//
            //printf("enter the refraction process");
            double m1 = 0;
            double m2 = 0;
            //if no mirrorwidth is specified no refraction has to be calc
            if (p->silicon == 0){
                return ga;
            }
            if (p->silicon == 1){
                m1 = m_Si;
                m2 = 0;
            }
            if (p->silicon == -1){
                m1 = 0;
                m2 = m_Si;                
            }
            //printf("k1 %f k2 %f k3 %f x %f y %f z %f\n", s.k1, s.k2, s.k3, p->_x, p->_y, p->_z);
            refractNeutronFlat(p, n, m1, m2);// this can still lead to total reflection, we miss the supermirror, but are still reflected by the silicon takes care of change of material for refraction
            return ga;
        }
    }
    return ga;
}


/*! \brief Function to handle raytracing of neutron for a ConicSurf.

@param p Pointer of particle to reflect
@param c ConicSurf to use

*/
void traceNeutronConic(Particle* p, ConicSurf c) {
    double t = getTimeOfFirstCollisionConic(*p, c);
    if (t < 0)
        return;
    else {
        moveParticleT(t, p);
        double ga = reflectNeutronConic(p, c);
#if REC_MAX_GA
        if (ga > c.max_ga) {
            c.max_ga = ga;
            c.max_ga_z0 = p->_z;
        }
#endif
    }
}

/*! \brief Function to handle raytracing of neutron for a FlatSurf.

@param p Pointer of particle to reflect
@param f FlatSurf to use

*/
void traceNeutronFlat(Particle* p, FlatSurf f) {
    double t = getTimeOfFirstCollisionFlat(*p, f);
    if (t < 0)
        return;
    else {

        moveParticleT(t, p);

        //printf("weight before reflect %f", p->w);
        double ga = reflectNeutronFlat(p, f);
        //printf("weight after reflect %f\n", p->w);
#if REC_MAX_GA
        if (ga > f.max_ga) {
            f.max_ga = ga;
            f.max_ga_z0 = p->_z;
        }
#endif
    }
}
/** @} */ //end of conicgroup
/////////////////////////////////////
// Scene Functions
/////////////////////////////////////
/** @ingroup simgroup
    @{
*/
enum GEO {
    NONE,
    DETECTOR,
    DISK,
    CONIC,
    FLAT
};

//! Function to generate an empty Scene
Scene makeScene() {
    Scene s;
    s.num_f = 0;
    s.num_c = 0;
    s.num_di = 0;
    s.num_d = 0;

    s.traceNeutronFlat = traceNeutronFlat;
    s.traceNeutronConic = traceNeutronConic;
    s.traceNeutronDisk = traceNeutronDisk;
    s.traceNeutronDetector = traceNeutronDetector;

    return s;
}

//! Function to init simulation items
/*! Should be called after all items
have been added to scene but before
neutrons are traced.

@param s Pointer of Scene to init
*/
void initSimulation(Scene* s) {
    //
}

/*! \brief Function to raytrace single neutron through geometries specified by d, di and c.

@param p Pointer of particle to trace
@param s Scene to trace
*/
void traceSingleNeutron(Particle* p, Scene s) {

    int contact = 1;
    do {
        double t;
        enum  GEO type = NONE;
        int index = -1;
        int i;

        for (i = 0; i < s.num_c; i++) {
            double t2 = getTimeOfFirstCollisionConic(*p,s.c[i]);

            if (t2 <= 0)
                continue;
            if (index == -1 || t2 < t) {
                type = CONIC;
                index = i;
                t = t2;
            }
        }

        for (i = 0; i < s.num_f; i++) {
            double t2 = getTimeOfFirstCollisionFlat(*p,s.f[i]);

            if (t2 <= 0)
                continue;
            if (index == -1 || t2 < t) {
                type = FLAT;
                index = i;
                t = t2;
            }
        }

        for (i = 0; i < s.num_di; i++)  {
            double t2 = getTimeOfFirstCollisionDisk(*p,s.di[i]);

            if (t2 <= 0)
                continue;
            else if (index == -1 || t2 < t) {
                type = DISK;
                index = i;
                t = t2;
            }
        }

        for (i = 0; i < s.num_d; i++) {
            double t2 = getTimeOfFirstCollisionDetector(*p,s.d[i]);

            if (t2 <= 0)
                continue;
            else if (index == -1 || t2 < t) {
                type = DETECTOR;
                index = i;
                t = t2;
            }
        }

        switch (type) {
            case DETECTOR:
                s.traceNeutronDetector(p, s.d[index]);
                break;
            case FLAT:
                s.traceNeutronFlat(p, s.f[index]);
                break;
            case DISK:
                s.traceNeutronDisk(p, s.di[index]);
                break;
            case CONIC:
                s.traceNeutronConic(p, s.c[index]);
                break;
            default:
                contact = 0;
                break;
        }
    } while (contact && !p->absorb);

}

//!Finishes tracing the scene
/*! This function should be called after all of the
particles have been raytraced.

@param s Pointer of Scene to finish tracing
*/
void finishSimulation(Scene* s) {
    int i;

    //Finish Detectors
    for (i=0; i < s->num_d; i++)
        finishDetector(s->d[i]);
}

/** @} */ //end of ingroup simgroup

#endif

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
double * get_r_at_z0(int number, double z_0, double r_0, double z_extract, double LStart, double LEnd, double lStart, double lEnd) {
    /*
    number: how many mirrors,
    z_0: z-position of the initial point on mirror
    r_0: r-distance of the initial point on mirror
    z_extract: z-position of which the distances are to be extracted
    LStart: Position of the first focal point
    LEnd: Postion of the second focal point
    lStart: Beginning of the mirror
    lEnd: End of the mirror
    */
    int n = number;
    double *r_zExtracts = malloc(n*sizeof(double_t)); /* n is an array of 10 integers */
	r_zExtracts[0] = r_0;
    //helper variables as in conic_finite_mirror.h and explained in swissneutronics_überlegungen
    double k1;
    double k2;
    double k3;
    double c;
    double u;
    double a;
    double r_lEnd;
    double r_lStart;
    //initial mirror is calculated from the initial point z0, r0
    c = (LEnd - LStart)/2;
    u = (z_0 + c - LEnd);
    a = sqrt((u*u+c*c+r_0*r_0+sqrt(pow(u*u+c*c+r_0*r_0, 2)-4*c*c*u*u))/2);
    k3 = c*c/(a*a)-1;
    k2 = 2*k3*(c-LEnd);
    k1 = k3*(c-LEnd)*(c-LEnd)-c*c+a*a;
    printf("k1 %f k2 %f k3 %f\n", k1, k2, k3);
	//next mirror will be calculated with the point on the surface being lStart, r_lStart
	for( int k = 0; k < number;++k){
        r_zExtracts[k] = sqrt(k1 + k2*z_extract + k3*z_extract*z_extract); 
        r_lEnd = sqrt(k1+ k2*lEnd + k3*lEnd*lEnd);//calculate the radius at the end
        r_lStart = r_lEnd*(lStart-LStart)/(lEnd-LStart);//

        c = (LEnd - LStart)/2;
        u = (lStart + c - LEnd);
        a = sqrt((u*u+c*c+r_lStart*r_lStart+sqrt(pow(u*u+c*c+r_lStart*r_lStart, 2)-4*c*c*u*u))/2);
        k3 = c*c/(a*a)-1;
        k2 = 2*k3*(c-LEnd);
        k1 = k3*(c-LEnd)*(c-LEnd)-c*c+a*a;
        printf("k1 %f k2 %f k3 %f\n", k1, k2, k3);
        //r_lEnd = sqrt(k1+ k2*lEnd + k3*lEnd*lEnd);
        //r_lStart = r_lEnd*(lStart-LStart)/(lEnd-LStart);
	};
   return r_zExtracts;
}




    //%include "w1_general.h"
    //%include "read_table-lib"
#line 7269 "./reverse_test.c"

/* Shared user declarations for all components 'LogSpiral'. */
#line 58 "LogSpiral.comp"
/*****************************************************************************
*
* McStas, neutron ray-tracing package
*         Copyright 1997-2006, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Library: share/ref-lib.h
*
* %Identification
* Written by: Peter Christiansen
* Date: August, 2006
* Origin: RISOE
* Release: McStas 1.10
* Version: $Revision$
*
* Commonly used reflection functions are declared in this file which
* are used by some guide and mirror components.
*
* Depends on read_table-lib
*
* Usage: within SHARE
* %include "ref-lib"
*
****************************************************************************/
/*******************************************************************************
*
* McStas, neutron ray-tracing package
*         Copyright 1997-2002, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Library: share/read_table-lib.h
*
* %Identification
* Written by: EF
* Date: Aug 28, 2002
* Origin: ILL
* Release: McStas 1.6
* Version: $Revision$
*
* This file is to be imported by components that may read data from table files
* It handles some shared functions.
*
* This library may be used directly as an external library. It has no dependency
*
* Usage: within SHARE
* %include "read_table-lib"
*
*******************************************************************************/

#ifndef READ_TABLE_LIB_H
#define READ_TABLE_LIB_H "$Revision$"

#define READ_TABLE_STEPTOL  0.04 /* tolerancy for constant step approx */

#ifndef MC_PATHSEP_C
#ifdef WIN32
#define MC_PATHSEP_C '\\'
#define MC_PATHSEP_S "\\"
#else  /* !WIN32 */
#ifdef MAC
#define MC_PATHSEP_C ':'
#define MC_PATHSEP_S ":"
#else  /* !MAC */
#define MC_PATHSEP_C '/'
#define MC_PATHSEP_S "/"
#endif /* !MAC */
#endif /* !WIN32 */
#endif /* !MC_PATHSEP_C */

#ifndef MCSTAS
#ifdef WIN32
#define MCSTAS "C:\\mcstas\\lib"
#else  /* !WIN32 */
#ifdef MAC
#define MCSTAS ":mcstas:lib" /* ToDo: What to put here? */
#else  /* !MAC */
#define MCSTAS "/usr/local/lib/mcstas"
#endif /* !MAC */
#endif /* !WIN32 */
#endif /* !MCSTAS */

#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>

  typedef struct struct_table
  {
    char    filename[1024];
    long    filesize;
    char   *header;  /* text header, e.g. comments */
    double *data;    /* vector { x[0], y[0], ... x[n-1], y[n-1]... } */
    double  min_x;   /* min value of first column */
    double  max_x;   /* max value of first column */
    double  step_x;  /* minimal step value of first column */
    long    rows;    /* number of rows in matrix block */
    long    columns; /* number of columns in matrix block */

    long    begin;   /* start fseek index of block */
    long    end;     /* stop  fseek index of block */
    long    block_number;  /* block index. 0 is catenation of all */
    long    array_length;  /* number of elements in the t_Table array */
    char    monotonic;     /* true when 1st column/vector data is monotonic */
    char    constantstep;  /* true when 1st column/vector data has constant step */
    char    method[32];    /* interpolation method: nearest, linear */
  } t_Table;

/*maximum number of rows to rebin a table = 1M*/
enum { mcread_table_rebin_maxsize = 1000000 };

typedef struct t_Read_table_file_item {
    int ref_count;
    t_Table *table_ref;
} t_Read_table_file_item;

typedef enum enum_Read_table_file_actions {STORE,FIND,GC}  t_Read_table_file_actions;

/* read_table-lib function prototypes */
/* ========================================================================= */

/* 'public' functions */
long     Table_Read              (t_Table *Table, char *File, long block_number);
long     Table_Read_Offset       (t_Table *Table, char *File, long block_number,
                                  long *offset, long max_lines);
long     Table_Read_Offset_Binary(t_Table *Table, char *File, char *Type,
                                  long *Offset, long Rows, long Columns);
long     Table_Rebin(t_Table *Table); /* rebin table with regular 1st column and interpolate all columns 2:end */
long     Table_Info (t_Table Table);
double   Table_Index(t_Table Table,   long i, long j); /* get indexed value */
double   Table_Value(t_Table Table, double X, long j); /* search X in 1st column and return interpolated value in j-column */
t_Table *Table_Read_Array(char *File, long *blocks);
void     Table_Free_Array(t_Table *Table);
long     Table_Info_Array(t_Table *Table);
int      Table_SetElement(t_Table *Table, long i, long j, double value);
long     Table_Init(t_Table *Table, long rows, long columns); /* create a Table */
double   Table_Value2d(t_Table Table, double X, double Y);    /* same as Table_Index with non-integer indices and 2d interpolation */
MCDETECTOR Table_Write(t_Table Table, char*file, char*xl, char*yl, 
           double x1, double x2, double y1, double y2); /* write Table to disk */
void * Table_File_List_Handler(t_Read_table_file_actions action, void *item, void *item_modifier);
t_Table *Table_File_List_find(char *name, int block, int offset);
int Table_File_List_gc(t_Table *tab);
void *Table_File_List_store(t_Table *tab);

#define Table_ParseHeader(header, ...) \
  Table_ParseHeader_backend(header,__VA_ARGS__,NULL);

char **Table_ParseHeader_backend(char *header, ...);

/* private functions */
void Table_Free(t_Table *Table);
long Table_Read_Handle(t_Table *Table, FILE *fid, long block_number, long max_lines, char *name);
static void Table_Stat(t_Table *Table);
double Table_Interp1d(double x, double x1, double y1, double x2, double y2);
double Table_Interp1d_nearest(double x, double x1, double y1, double x2, double y2);
double Table_Interp2d(double x, double y, double x1, double y1, double x2, double y2,
double z11, double z12, double z21, double z22);

#endif

/* end of read_table-lib.h */
/*******************************************************************************
*
* McStas, neutron ray-tracing package
*         Copyright (C) 1997-2009, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Library: share/read_table-lib.c
*
* %Identification
* Written by: EF
* Date: Aug 28, 2002
* Origin: ILL
* Release: McStas CVS_090504
* Version: $Revision: 5052 $
*
* This file is to be imported by components that may read data from table files
* It handles some shared functions. Embedded within instrument in runtime mode.
*
* Usage: within SHARE
* %include "read_table-lib"
*
*******************************************************************************/

#ifndef READ_TABLE_LIB_H
#include "read_table-lib.h"
#endif


/*******************************************************************************
 * void *Table_File_List_Handler(action, item, item_modifier)
 *   ACTION: handle file entries in the read_table-lib file list. If a file is read - it is supposed to be
 *   stored in a list such that we can avoid reading the same file many times.
 *   input  action: FIND, STORE, GC. check if file exists in the list, store an item in the list, or check if it can be garbage collected.
 *   input item: depends on the action.
 *    FIND)  item is a filename, and item_modifier is the block number
 *    STORE) item is the Table to store - item_modifier is ignored
 *    GC)    item is the Table to check. If it has a ref_count >1 then this is simply decremented.
 *   return  depends on the action
 *    FIND)  return a reference to a table+ref_count item if found - NULL otherwise. I.e. NULL means the file has not been read before and must be read again.
 *    STORE) return NULL always
 *    GC)    return NULL if no garbage collection is needed, return an adress to the t_Table which should be garbage collected. 0x1 is returned if
 *           the item is not found in the list
*******************************************************************************/
void * Table_File_List_Handler(t_Read_table_file_actions action, void *item, void *item_modifier){

    /* logic here is Read_Table should include a call to FIND. If found the return value should just be used as
     * if the table had been read from disk. If not found then read the table and STORE.
     * Table_Free should include a call to GC. If this returns non-NULL then we should proceed with freeing the memory
     * associated with the table item - otherwise only decrement the reference counter since there are more references
     * that may need it.*/

    static t_Read_table_file_item read_table_file_list[1024];  
    static int read_table_file_count=0;

    t_Read_table_file_item *tr;
    switch(action){
        case FIND:
            /*interpret data item as a filename, if it is found return a pointer to the table and increment refcount.
             * if not found return the item itself*/
            tr=read_table_file_list;
            while ( tr->table_ref!=NULL ){
                int i=*((int*) item_modifier);
                int j=*( ((int*) item_modifier)+1);
                if ( !strcmp(tr->table_ref->filename,(char *) item) &&
                        tr->table_ref->block_number==i && tr->table_ref->begin==j ){
                    tr->ref_count++;
                    return (void *) tr;
                }
                tr++;
            }
            return NULL;
        case STORE:
            /*find an available slot and store references to table there*/
            tr=&(read_table_file_list[read_table_file_count++]);
            tr->table_ref = ((t_Table *) item);
            tr->ref_count++;
            return NULL;
        case GC:
            /* Should this item be garbage collected (freed) - if so scratch the entry and return the address of the item - 
             * else decrement ref_count and return NULL.
             * A non-NULL return expects the item to actually be freed afterwards.*/
            tr=read_table_file_list;
            while ( tr->table_ref!=NULL ){
                if ( tr->table_ref->data ==((t_Table *)item)->data && 
                        tr->table_ref->block_number == ((t_Table *)item)->block_number){
                    /*matching item found*/
                    if (tr->ref_count>1){
                        /*the item is found and no garbage collection needed*/
                        tr->ref_count--;
                        return NULL;
                    }else{
                        /* The item is found and the reference counter is 1.
                         * This means we should garbage collect. Move remaining list items up one slot,
                         * and return the table for garbage collection by caller*/
                        while (tr->table_ref!=NULL){
                            *tr=*(tr+1);
                            tr++;
                        }
                        read_table_file_count--;
                        return (t_Table *) item;
                    }
                }
                tr++;
            }
            /* item not found, and so should be garbage collected. This could be the case if freeing a
             * Table that has been constructed from code - not read from file. Return 0x1 to flag it for
             * collection.*/
            return (void *) 0x1 ;
    }
}

/* Access functions to the handler*/

/********************************************
 * t_Table *Table_File_List_find(char *name, int block, int offset)
 * input name: filename to search for in the file list
 * input block: data block in the file as each file may contain more than 1 data block.
 * return a ref. to a table if it is found (you may use this pointer and skip reading the file), NULL otherwise (i.e. go ahead and read the file)
*********************************************/
t_Table *Table_File_List_find(char *name, int block, int offset){
    int vars[2]={block,offset};
    t_Read_table_file_item *item = Table_File_List_Handler(FIND,name, vars);
    if (item == NULL){
        return NULL;
    }else{
        return item->table_ref;
    }
}
/********************************************
 * int Table_File_List_gc(t_Table *tab)
 * input tab: the table to check for references.
 * return 0: no garbage collection needed
 *        1: Table's data and header (at least) should be freed.
*********************************************/
int Table_File_List_gc(t_Table *tab){
    void *rval=Table_File_List_Handler(GC,tab,0);
    if (rval==NULL) return 0;
    else return 1;
}


/*****************************************************************************
 * void *Table_File_List_store(t_Table *tab)
 * input tab: pointer to table to store.
 * return None. 
*******************************************************************************/
void *Table_File_List_store(t_Table *tab){
    return Table_File_List_Handler(STORE,tab,0);
}


/*******************************************************************************
* FILE *Open_File(char *name, char *Mode, char *path)
*   ACTION: search for a file and open it. Optionally return the opened path.
*   input   name:  file name from which table should be extracted
*           mode: "r", "w", "a" or any valid fopen mode
*           path:  NULL or a pointer to at least 1024 allocated chars
*   return  initialized file handle or NULL in case of error
*******************************************************************************/

  FILE *Open_File(char *File, const char *Mode, char *Path)
  {
    char path[1024];
    FILE *hfile = NULL;
    
    if (!File || File[0]=='\0')                     return(NULL);
    if (!strcmp(File,"NULL") || !strcmp(File,"0"))  return(NULL);
    
    /* search in current or full path */
    strncpy(path, File, 1024);
    hfile = fopen(path, Mode);
    if(!hfile)
    {
      char dir[1024];

      if (!hfile && mcinstrument_source[0] != '\0' && strlen(mcinstrument_source)) /* search in instrument source location */
      {
        char *path_pos   = NULL;
        /* extract path: searches for last file separator */
        path_pos    = strrchr(mcinstrument_source, MC_PATHSEP_C);  /* last PATHSEP */
        if (path_pos) {
          long path_length = path_pos +1 - mcinstrument_source;  /* from start to path+sep */
          if (path_length) {
            strncpy(dir, mcinstrument_source, path_length);
            dir[path_length] = '\0';
            snprintf(path, 1024, "%s%c%s", dir, MC_PATHSEP_C, File);
            hfile = fopen(path, Mode);
          }
        }
      }
      if (!hfile && mcinstrument_exe[0] != '\0' && strlen(mcinstrument_exe)) /* search in PWD instrument executable location */
      {
        char *path_pos   = NULL;
        /* extract path: searches for last file separator */
        path_pos    = strrchr(mcinstrument_exe, MC_PATHSEP_C);  /* last PATHSEP */
        if (path_pos) {
          long path_length = path_pos +1 - mcinstrument_exe;  /* from start to path+sep */
          if (path_length) {
            strncpy(dir, mcinstrument_exe, path_length);
            dir[path_length] = '\0';
            snprintf(path, 1024, "%s%c%s", dir, MC_PATHSEP_C, File);
            hfile = fopen(path, Mode);
          }
        }
      }
      if (!hfile) /* search in HOME or . */
      {
        strcpy(dir, getenv("HOME") ? getenv("HOME") : ".");
        snprintf(path, 1024, "%s%c%s", dir, MC_PATHSEP_C, File);
        hfile = fopen(path, Mode);
      }
      if (!hfile) /* search in MCSTAS/data */
      {
        strcpy(dir, getenv(FLAVOR_UPPER) ? getenv(FLAVOR_UPPER) : MCSTAS);
        snprintf(path, 1024, "%s%c%s%c%s", dir, MC_PATHSEP_C, "data", MC_PATHSEP_C, File);
        hfile = fopen(path, Mode);
      }
      if (!hfile) /* search in MVCSTAS/contrib */
      {
        strcpy(dir, getenv(FLAVOR_UPPER) ? getenv(FLAVOR_UPPER) : MCSTAS);
        snprintf(path, 1024, "%s%c%s%c%s", dir, MC_PATHSEP_C, "contrib", MC_PATHSEP_C, File);
        hfile = fopen(path, Mode);
      }
      if(!hfile)
      {
        fprintf(stderr, "Error: Could not open input file '%s' (Open_File)\n", File);
        return (NULL);
      }
    }
    if (Path) strncpy(Path, path, 1024);
    return(hfile);
  } /* end Open_File */

/*******************************************************************************
* long Read_Table(t_Table *Table, char *name, int block_number)
*   ACTION: read a single Table from a text file
*   input   Table: pointer to a t_Table structure
*           name:  file name from which table should be extracted
*           block_number: if the file does contain more than one
*                 data block, then indicates which one to get (from index 1)
*                 a 0 value means append/catenate all
*   return  initialized single Table t_Table structure containing data, header, ...
*           number of read elements (-1: error, 0:header only)
* The routine stores any line starting with '#', '%' and ';' into the header
* File is opened, read and closed
* Other lines are interpreted as numerical data, and stored.
* Data block should be a rectangular matrix or vector.
* Data block may be rebinned with Table_Rebin (also sort in ascending order)
*******************************************************************************/
  long Table_Read(t_Table *Table, char *File, long block_number)
  { /* reads all or a single data block from 'file' and returns a Table structure  */
    return(Table_Read_Offset(Table, File, block_number, NULL, 0));
  } /* end Table_Read */

/*******************************************************************************
* long Table_Read_Offset(t_Table *Table, char *name, int block_number, long *offset
*                        long max_rows)
*   ACTION: read a single Table from a text file, starting at offset
*     Same as Table_Read(..) except:
*   input   offset:    pointer to an offset (*offset should be 0 at start)
*           max_rows: max number of data rows to read from file (0 means all)
*   return  initialized single Table t_Table structure containing data, header, ...
*           number of read elements (-1: error, 0:header only)
*           updated *offset position (where end of reading occured)
*******************************************************************************/
  long Table_Read_Offset(t_Table *Table, char *File,
                         long block_number, long *offset,
                         long max_rows)
  { /* reads all/a data block in 'file' and returns a Table structure  */
    FILE *hfile;
    long  nelements=0;
    long  begin=0;
    long  filesize=0;
    char  name[1024];
    char  path[1024];
    struct stat stfile;

    /*Need to be able to store the pointer*/
    if (!Table) return(-1);
    
    //if (offset && *offset) snprintf(name, 1024, "%s@%li", File, *offset);
    //else                   
    strncpy(name, File, 1024);
    if(offset && *offset){
        begin=*offset;
    }
    /* Check if the table has already been read from file.
     * If so just reuse the table, if not (this is flagged by returning NULL
     * set up a new table and read the data into it */
    t_Table *tab_p= Table_File_List_find(name,block_number,begin);
    if ( tab_p!=NULL ){
        /*table was found in the Table_File_List*/
        // printf("Reusing input file '%s' (Table_Read_Offset)\n", name);
        *Table=*tab_p;
        return Table->rows*Table->columns;
    }

    /* open the file */
    hfile = Open_File(File, "r", path);
    if (!hfile) return(-1);
    else {
      MPI_MASTER(
      printf("Opening input file '%s' (Table_Read_Offset)\n", path);
      );
    }
    
    /* read file state */
    stat(path,&stfile); filesize = stfile.st_size;
    if (offset && *offset) fseek(hfile, *offset, SEEK_SET);
    begin     = ftell(hfile);
    
    Table_Init(Table, 0, 0);

    /* read file content and set the Table */
    nelements = Table_Read_Handle(Table, hfile, block_number, max_rows, name);
    Table->begin = begin;
    Table->end   = ftell(hfile);
    Table->filesize = (filesize>0 ? filesize : 0);
    Table_Stat(Table);
    
    Table_File_List_store(Table);

    if (offset) *offset=Table->end;
    fclose(hfile);
    return(nelements);

  } /* end Table_Read_Offset */

/*******************************************************************************
* long Table_Read_Offset_Binary(t_Table *Table, char *File, char *type,
*                               long *offset, long rows, long columns)
*   ACTION: read a single Table from a binary file, starting at offset
*     Same as Table_Read_Offset(..) except that it handles binary files.
*   input   type: may be "float"/NULL or "double"
*           offset: pointer to an offset (*offset should be 0 at start)
*           rows   : number of rows (0 means read all)
*           columns: number of columns
*   return  initialized single Table t_Table structure containing data, header, ...
*           number of read elements (-1: error, 0:header only)
*           updated *offset position (where end of reading occured)
*******************************************************************************/
  long Table_Read_Offset_Binary(t_Table *Table, char *File, char *type,
                                long *offset, long rows, long columns)
  { /* reads all/a data block in binary 'file' and returns a Table structure  */
    long    nelements, sizeofelement;
    long    filesize;
    FILE   *hfile;
    char    path[1024];
    struct stat stfile;
    double *data;
    long    i;
    long    begin;

    if (!Table) return(-1);

    Table_Init(Table, 0, 0);
    
    /* open the file */
    hfile = Open_File(File, "r", path);
    if (!hfile) return(-1);
    else {
      MPI_MASTER(
      printf("Opening input file '%s' (Table_Read, Binary)\n", path);
      );
    }
    
    /* read file state */
    stat(File,&stfile);
    filesize = stfile.st_size;
    Table->filesize=filesize;
    
    /* read file content */
    if (type && !strcmp(type,"double")) sizeofelement = sizeof(double);
    else  sizeofelement = sizeof(float);
    if (offset && *offset) fseek(hfile, *offset, SEEK_SET);
    begin     = ftell(hfile);
    if (rows && filesize > sizeofelement*columns*rows)
      nelements = columns*rows;
    else nelements = (long)(filesize/sizeofelement);
    if (!nelements || filesize <= *offset) return(0);
    data    = (double*)malloc(nelements*sizeofelement);
    if (!data) {
      fprintf(stderr,"Error: allocating %ld elements for %s file '%s'. Too big (Table_Read_Offset_Binary).\n", nelements, type, File);
      exit(-1);
    }
    nelements = fread(data, sizeofelement, nelements, hfile);

    if (!data || !nelements)
    {
      fprintf(stderr,"Error: reading %ld elements from %s file '%s' (Table_Read_Offset_Binary)\n", nelements, type, File);
      exit(-1);
    }
    Table->begin   = begin;
    Table->end     = ftell(hfile);
    if (offset) *offset=Table->end;
    fclose(hfile);
    data = (double*)realloc(data, (double)nelements*sizeofelement);
    /* copy file data into Table */
    if (type && !strcmp(type,"double")) Table->data = data;
    else {
      float  *s;
      double *dataf;
      s     = (float*)data;
      dataf = (double*)malloc(sizeof(double)*nelements);
      for (i=0; i<nelements; i++)
        dataf[i]=s[i];
      free(data);
      Table->data = dataf;
    }
    strncpy(Table->filename, File, 1024);
    Table->rows    = nelements/columns;
    Table->columns = columns;
    Table->array_length = 1;
    Table->block_number = 1;

    Table_Stat(Table);

    return(nelements);
  } /* end Table_Read_Offset_Binary */

/*******************************************************************************
* long Table_Read_Handle(t_Table *Table, FILE *fid, int block_number, long max_rows, char *name)
*   ACTION: read a single Table from a text file handle (private)
*   input   Table:pointer to a t_Table structure
*           fid:  pointer to FILE handle
*           block_number: if the file does contain more than one
*                 data block, then indicates which one to get (from index 1)
*                 a 0 value means append/catenate all
*           max_rows: if non 0, only reads that number of lines
*   return  initialized single Table t_Table structure containing data, header, ...
*           modified Table t_Table structure containing data, header, ...
*           number of read elements (-1: error, 0:header only)
* The routine stores any line starting with '#', '%' and ';' into the header
* Other lines are interpreted as numerical data, and stored.
* Data block should be a rectangular matrix or vector.
* Data block may be rebined with Table_Rebin (also sort in ascending order)
*******************************************************************************/
  long Table_Read_Handle(t_Table *Table, FILE *hfile,
                         long block_number, long max_rows, char *name)
  { /* reads all/a data block from 'file' handle and returns a Table structure  */
    double *Data;
    char *Header              = NULL;
    long  malloc_size         = CHAR_BUF_LENGTH;
    long  malloc_size_h       = 4096;
    long  Rows = 0,   Columns = 0;
    long  count_in_array      = 0;
    long  count_in_header     = 0;
    long  block_Current_index = 0;
    char  flag_End_row_loop   = 0;

    if (!Table) return(-1);
    Table_Init(Table, 0, 0);
    if (name && name[0]!='\0') strncpy(Table->filename, name, 1024);

    if(!hfile) {
       fprintf(stderr, "Error: File handle is NULL (Table_Read_Handle).\n");
       return (-1);
    }
    Header = (char*)  calloc(malloc_size_h, sizeof(char));
    Data   = (double*)calloc(malloc_size,   sizeof(double));
    if ((Header == NULL) || (Data == NULL)) {
       fprintf(stderr, "Error: Could not allocate Table and Header (Table_Read_Handle).\n");
       return (-1);
    }

    int flag_In_array = 0;
    do { /* while (!flag_End_row_loop) */
      char  line[1024*CHAR_BUF_LENGTH];
      long  back_pos=0;   /* ftell start of line */

      back_pos = ftell(hfile);
      if (fgets(line, 1024*CHAR_BUF_LENGTH, hfile) != NULL) { /* analyse line */
        /* first skip blank and tabulation characters */
        int i = strspn(line, " \t");

        /* handle comments: stored in header */
        if (NULL != strchr("#%;/", line[i]))
        { /* line is a comment */
          count_in_header += strlen(line);
          if (count_in_header >= malloc_size_h) {
            /* if succeed and in array : add (and realloc if necessary) */
            malloc_size_h = count_in_header+4096;
            Header        = (char*)realloc(Header, malloc_size_h*sizeof(char));
          }
          strncat(Header, line, 4096);
          flag_In_array=0;
          /* exit line and file if passed desired block */
          if (block_number > 0 && block_number == block_Current_index) {
            flag_End_row_loop = 1;
          }

          /* Continue with next line */
          continue;
        }

        /* get the number of columns splitting line with strtok */
        char  *lexeme;
        char  flag_End_Line = 0;
        long  block_Num_Columns = 0;
        const char seps[] = " ,;\t\n\r";

        lexeme = strtok(line, seps);
        while (!flag_End_Line) {
          if ((lexeme != NULL) && (lexeme[0] != '\0')) {
            /* reading line: the token is not empty */
            double X;
            int    count=1;
            /* test if we have 'NaN','Inf' */
            if (!strncasecmp(lexeme,"NaN",3))
              X = 0;
            else if (!strncasecmp(lexeme,"Inf",3) || !strncasecmp(lexeme,"+Inf",4))
              X = FLT_MAX;
            else if (!strncasecmp(lexeme,"-Inf",4))
              X = -FLT_MAX;
            else
              count = sscanf(lexeme,"%lg",&X);
            if (count == 1) {
              /* reading line: the token is a number in the line */
              if (!flag_In_array) {
                /* reading num: not already in a block: starts a new data block */
                block_Current_index++;
                flag_In_array    = 1;
                block_Num_Columns= 0;
                if (block_number > 0) {
                  /* initialise a new data block */
                  Rows = 0;
                  count_in_array = 0;
                } /* else append */
              }
              /* reading num: all blocks or selected block */
              if (flag_In_array && (block_number == 0 ||
                  block_number == block_Current_index)) {
                /* starting block: already the desired number of rows ? */
                if (block_Num_Columns == 0 &&
                    max_rows > 0 && Rows >= max_rows) {
                  flag_End_Line      = 1;
                  flag_End_row_loop  = 1;
                  flag_In_array      = 0;
                  /* reposition to begining of line (ignore line) */
                  fseek(hfile, back_pos, SEEK_SET);
                } else { /* store into data array */
                  if (count_in_array >= malloc_size) {
                    /* realloc data buffer if necessary */
                    malloc_size = count_in_array*1.5;
                    Data = (double*) realloc(Data, malloc_size*sizeof(double));
                    if (Data == NULL) {
                      fprintf(stderr, "Error: Can not re-allocate memory %li (Table_Read_Handle).\n",
                              malloc_size*sizeof(double));
                      return (-1);
                    }
                  }
                  if (0 == block_Num_Columns) Rows++;
                  Data[count_in_array] = X;
                  count_in_array++;
                  block_Num_Columns++;
                }
              } /* reading num: end if flag_In_array */
            } /* end reading num: end if sscanf lexeme -> numerical */
            else {
              /* reading line: the token is not numerical in that line. end block */
              if (block_Current_index == block_number) {
                flag_End_Line = 1;
                flag_End_row_loop = 1;
              } else {
                flag_In_array = 0;
                flag_End_Line = 1;
              }
            }
          }
          else {
            /* no more tokens in line */
            flag_End_Line = 1;
            if (block_Num_Columns > 0) Columns = block_Num_Columns;
          }

          // parse next token
          lexeme = strtok(NULL, seps);

        } /* while (!flag_End_Line) */
      } /* end: if fgets */
      else flag_End_row_loop = 1; /* else fgets : end of file */

    } while (!flag_End_row_loop); /* end while flag_End_row_loop */

    Table->block_number = block_number;
    Table->array_length = 1;

    // shrink header to actual size (plus terminating 0-byte)
    if (count_in_header) {
      Header = (char*)realloc(Header, count_in_header*sizeof(char) + 1);
    }
    Table->header = Header;

    if (count_in_array*Rows*Columns == 0)
    {
      Table->rows         = 0;
      Table->columns      = 0;
      free(Data);
      return (0);
    }
    if (Rows * Columns != count_in_array)
    {
      fprintf(stderr, "Warning: Read_Table :%s %s Data has %li values that should be %li x %li\n",
        (Table->filename[0] != '\0' ? Table->filename : ""),
        (!block_number ? " catenated" : ""),
        count_in_array, Rows, Columns);
      Columns = count_in_array; Rows = 1;
    }
    Data     = (double*)realloc(Data, count_in_array*sizeof(double));
    Table->data         = Data;
    Table->rows         = Rows;
    Table->columns      = Columns;

    return (count_in_array);

  } /* end Table_Read_Handle */

/*******************************************************************************
* long Table_Rebin(t_Table *Table)
*   ACTION: rebin a single Table, sorting 1st column in ascending order
*   input   Table: single table containing data.
*                  The data block is reallocated in this process
*   return  updated Table with increasing, evenly spaced first column (index 0)
*           number of data elements (-1: error, 0:empty data)
*******************************************************************************/
  long Table_Rebin(t_Table *Table)
  {
    double new_step=0;
    long   i;
    /* performs linear interpolation on X axis (0-th column) */

    if (!Table) return(-1);
    if (!Table->data 
    || Table->rows*Table->columns == 0 || !Table->step_x)
      return(0);
    Table_Stat(Table); /* recompute statitstics and minimal step */
    new_step = Table->step_x; /* minimal step in 1st column */

    if (!(Table->constantstep)) /* not already evenly spaced */
    {
      long Length_Table;
      double *New_Table;

      Length_Table = ceil(fabs(Table->max_x - Table->min_x)/new_step)+1;
      /*return early if the rebinned table will become too large*/
      if (Length_Table > mcread_table_rebin_maxsize){
        fprintf(stderr,"WARNING: (Table_Rebin): Rebinning table from %s would exceed 1M rows. Skipping.\n", Table->filename); 
        return(Table->rows*Table->columns);
      }
      New_Table    = (double*)malloc(Length_Table*Table->columns*sizeof(double));

      for (i=0; i < Length_Table; i++)
      {
        long   j;
        double X;
        X = Table->min_x + i*new_step;
        New_Table[i*Table->columns] = X;
        for (j=1; j < Table->columns; j++)
          New_Table[i*Table->columns+j]
                = Table_Value(*Table, X, j);
      } /* end for i */

      Table->rows = Length_Table;
      Table->step_x = new_step;
      Table->max_x = Table->min_x + (Length_Table-1)*new_step; 
      /*max might not be the same anymore
       * Use Length_Table -1 since the first and laset rows are the limits of the defined interval.*/
      free(Table->data);
      Table->data = New_Table;
      Table->constantstep=1;
    } /* end else (!constantstep) */
    return (Table->rows*Table->columns);
  } /* end Table_Rebin */

/*******************************************************************************
* double Table_Index(t_Table Table, long i, long j)
*   ACTION: read an element [i,j] of a single Table
*   input   Table: table containing data
*           i : index of row      (0:Rows-1)
*           j : index of column   (0:Columns-1)
*   return  Value = data[i][j]
* Returns Value from the i-th row, j-th column of Table
* Tests are performed on indexes i,j to avoid errors
*******************************************************************************/

#ifndef MIN
#define MIN(a, b)  (((a) < (b)) ? (a) : (b))
#endif
#ifndef MAX
#define MAX(a, b)  (((a) > (b)) ? (a) : (b))
#endif

double Table_Index(t_Table Table, long i, long j)
{
  long AbsIndex;

  if (Table.rows == 1 || Table.columns == 1) {
    /* vector */
    j = MIN(MAX(0, i+j), Table.columns*Table.rows - 1);
    i = 0;
  } else {
    /* matrix */
    i = MIN(MAX(0, i), Table.rows - 1);
    j = MIN(MAX(0, j), Table.columns - 1);
  }

  /* handle vectors specifically */
  AbsIndex = i*(Table.columns)+j;

  if (Table.data != NULL)
    return (Table.data[AbsIndex]);
  else
    return 0;
} /* end Table_Index */

/*******************************************************************************
* void Table_SetElement(t_Table *Table, long i, long j, double value)
*   ACTION: set an element [i,j] of a single Table
*   input   Table: table containing data
*           i : index of row      (0:Rows-1)
*           j : index of column   (0:Columns-1)
*           value = data[i][j]
* Returns 0 in case of error
* Tests are performed on indexes i,j to avoid errors
*******************************************************************************/
int Table_SetElement(t_Table *Table, long i, long j,
                     double value)
{
  long AbsIndex;

  if (Table->rows == 1 || Table->columns == 1) {
    /* vector */
    j = MIN(MAX(0, i+j), Table->columns*Table->rows - 1); i=0;
  } else {
    /* matrix */
    i = MIN(MAX(0, i), Table->rows - 1);
    j = MIN(MAX(0, j), Table->columns - 1);
  }

  AbsIndex = i*(Table->columns)+j;
  if (Table->data != NULL) {
    Table->data[AbsIndex] = value;
    return 1;
  }

  return 0;
} /* end Table_SetElement */

/*******************************************************************************
* double Table_Value(t_Table Table, double X, long j)
*   ACTION: read column [j] of a single Table at row which 1st column is X
*   input   Table: table containing data.
*           X : data value in the first column (index 0)
*           j : index of column from which is extracted the Value (0:Columns-1)
*   return  Value = data[index for X][j] with linear interpolation
* Returns Value from the j-th column of Table corresponding to the
* X value for the 1st column (index 0)
* Tests are performed (within Table_Index) on indexes i,j to avoid errors
* NOTE: data should rather be monotonic, and evenly sampled.
*******************************************************************************/
double Table_Value(t_Table Table, double X, long j)
{
  long   Index = -1;
  double X1=0, Y1=0, X2=0, Y2=0;
  double ret=0;

  if (X > Table.max_x) return Table_Index(Table,Table.rows-1  ,j);
  if (X < Table.min_x) return Table_Index(Table,0  ,j);

  // Use constant-time lookup when possible
  if(Table.constantstep) {
    Index = (long)floor(
              (X - Table.min_x) / (Table.max_x - Table.min_x) * (Table.rows-1));
    X1 = Table_Index(Table,Index  ,0);
    X2 = Table_Index(Table,Index+1,0);
  }
  // Use binary search on large, monotonic tables
  else if(Table.monotonic && Table.rows > 100) {
    long left = Table.min_x;
    long right = Table.max_x;

    while (!((X1 <= X) && (X < X2)) && (right - left > 1)) {
      Index = (left + right) / 2;

      X1 = Table_Index(Table, Index-1, 0);
      X2 = Table_Index(Table, Index,   0);

      if (X < X1) {
        right = Index;
      } else {
        left  = Index;
      }
    }
  }

  // Fall back to linear search, if no-one else has set X1, X2 correctly
  if (!((X1 <= X) && (X < X2))) {
    /* look for index surrounding X in the table -> Index */
    for (Index=1; Index <= Table.rows-1; Index++) {
        X1 = Table_Index(Table, Index-1,0);
        X2 = Table_Index(Table, Index  ,0);
        if ((X1 <= X) && (X < X2)) break;
      } /* end for Index */
  }

  Y1 = Table_Index(Table,Index, j);
  Y2 = Table_Index(Table,Index+1, j);

  if (!strcmp(Table.method,"linear")) {
    ret = Table_Interp1d(X, X1,Y1, X2,Y2);
  }
  else if (!strcmp(Table.method,"nearest")) {
    ret = Table_Interp1d_nearest(X, X1,Y1, X2,Y2);
  }

  return ret;
} /* end Table_Value */

/*******************************************************************************
* double Table_Value2d(t_Table Table, double X, double Y)
*   ACTION: read element [X,Y] of a matrix Table
*   input   Table: table containing data.
*           X : row index, may be non integer
*           Y : column index, may be non integer
*   return  Value = data[index X][index Y] with bi-linear interpolation
* Returns Value for the indices [X,Y]
* Tests are performed (within Table_Index) on indexes i,j to avoid errors
* NOTE: data should rather be monotonic, and evenly sampled.
*******************************************************************************/
  double Table_Value2d(t_Table Table, double X, double Y)
  {
    long   x1,x2,y1,y2;
    double z11,z12,z21,z22;
    double ret=0;

    x1 = (long)floor(X);
    y1 = (long)floor(Y);

    if (x1 > Table.rows-1 || x1 < 0) {
      x2 = x1;
    } else {
      x2 = x1 + 1;
    }

    if (y1 > Table.columns-1 || y1 < 0) {
      y2 = y1;
    } else {
      y2 = y1 + 1;
    }

    z11 = Table_Index(Table, x1, y1);

    if (y2 != y1) z12=Table_Index(Table, x1, y2); else z12 = z11;
    if (x2 != x1) z21=Table_Index(Table, x2, y1); else z21 = z11;
    if (y2 != y1) z22=Table_Index(Table, x2, y2); else z22 = z21;

    if (!strcmp(Table.method,"linear"))
      ret = Table_Interp2d(X,Y, x1,y1,x2,y2, z11,z12,z21,z22);
    else {
      if (fabs(X-x1) < fabs(X-x2)) {
        if (fabs(Y-y1) < fabs(Y-y2)) ret = z11; else ret = z12;
      } else {
        if (fabs(Y-y1) < fabs(Y-y2)) ret = z21; else ret = z22;
      }
    }
    return ret;
  } /* end Table_Value2d */


/*******************************************************************************
* void Table_Free(t_Table *Table)
*   ACTION: free a single Table. First Call Table_File_list_gc. If this returns
*   non-zero it means there are more refernces to the table, and so the table
*   should not bee freed.
*   return: empty Table
*******************************************************************************/
  void Table_Free(t_Table *Table)
  {
    if( !Table_File_List_gc(Table) ){
       return;
    } 
    if (!Table) return;
    if (Table->data   != NULL) free(Table->data);
    if (Table->header != NULL) free(Table->header);
    Table->data   = NULL;
    Table->header = NULL;
  } /* end Table_Free */

/******************************************************************************
* void Table_Info(t_Table Table)
*    ACTION: print informations about a single Table
*******************************************************************************/
  long Table_Info(t_Table Table)
  {
    char buffer[256];
    long ret=0;

    if (!Table.block_number) strcpy(buffer, "catenated");
    else sprintf(buffer, "block %li", Table.block_number);
    printf("Table from file '%s' (%s)",
      Table.filename[0] != '\0' ? Table.filename : "", buffer);
    if ((Table.data != NULL) && (Table.rows*Table.columns))
    {
      printf(" is %li x %li ", Table.rows, Table.columns);
      if (Table.rows*Table.columns > 1)
        printf("(x=%g:%g)", Table.min_x, Table.max_x);
      else printf("(x=%g) ", Table.min_x);
      ret = Table.rows*Table.columns;
      if (Table.monotonic)    printf(", monotonic");
      if (Table.constantstep) printf(", constant step");
      printf(". interpolation: %s\n", Table.method);
    }
    else printf(" is empty.\n");

    if (Table.header && strlen(Table.header)) {
      char *header;
      int  i;
      header = malloc(80);
      if (!header) return(ret);
      for (i=0; i<80; header[i++]=0);
      strncpy(header, Table.header, 75);
      if (strlen(Table.header) > 75) {
        strcat( header, " ...");
      }
      for (i=0; i<strlen(header); i++)
        if (header[i] == '\n' || header[i] == '\r') header[i] = ';';
      printf("  '%s'\n", header);
      free(header);
    }

    return(ret);
  } /* end Table_Info */

/******************************************************************************
* long Table_Init(t_Table *Table, m, n)
*   ACTION: initialise a Table to empty m by n table
*   return: empty Table
******************************************************************************/
long Table_Init(t_Table *Table, long rows, long columns)
{
  double *data=NULL;
  long   i;

  if (!Table) return(0);

  Table->header  = NULL;
  Table->filename[0]= '\0';
  Table->filesize= 0;
  Table->min_x   = 0;
  Table->max_x   = 0;
  Table->step_x  = 0;
  Table->block_number = 0;
  Table->array_length = 0;
  Table->monotonic    = 0;
  Table->constantstep = 0;
  Table->begin   = 0;
  Table->end     = 0;
  strcpy(Table->method,"linear");

  if (rows*columns >= 1) {
    data    = (double*)malloc(rows*columns*sizeof(double));
    if (data) for (i=0; i < rows*columns; data[i++]=0);
    else {
      fprintf(stderr,"Error: allocating %ld double elements."
                     "Too big (Table_Init).\n", rows*columns);
      rows = columns = 0;
    }
  }
  Table->rows    = (rows >= 1 ? rows : 0);
  Table->columns = (columns >= 1 ? columns : 0);
  Table->data    = data;
  return(Table->rows*Table->columns);
} /* end Table_Init */

/******************************************************************************
* long Table_Write(t_Table Table, char *file, x1,x2, y1,y2)
*   ACTION: write a Table to disk (ascii).
*     when x1=x2=0 or y1=y2=0, the table default limits are used.
*   return: 0=all is fine, non-0: error
*******************************************************************************/
MCDETECTOR Table_Write(t_Table Table, char *file, char *xl, char *yl, 
  double x1, double x2, double y1, double y2)
{
  long    i =0;
  MCDETECTOR detector;

  if ((Table.data == NULL) && (Table.rows*Table.columns)) {
    detector.m = 0;
    return(detector); /* Table is empty - nothing to do */
  }
  if (!x1 && !x2) {
    x1 = Table.min_x;
    x2 = Table.max_x;
  }
  if (!y1 && !y2) {
    y1 = 1;
    y2 = Table.columns;
  }

  /* transfer content of the Table into a 2D detector */
  Coords coords = { 0, 0, 0};

  if (Table.rows == 1 || Table.columns == 1) {
    detector = mcdetector_out_1D(Table.filename,
                      xl ? xl : "", yl ? yl : "",
                      "x", x1, x2,
                      Table.rows * Table.columns,
                      NULL, Table.data, NULL,
                      file, file, coords);
  } else {
    detector = mcdetector_out_2D(Table.filename,
                      xl ? xl : "", yl ? yl : "",
                      x1, x2, y1, y2,
                      Table.rows, Table.columns,
                      NULL, Table.data, NULL,
                      file, file, coords);
  }
  return(detector);
}

/******************************************************************************
* void Table_Stat(t_Table *Table)
*   ACTION: computes min/max/mean step of 1st column for a single table (private)
*   return: updated Table
*******************************************************************************/
  static void Table_Stat(t_Table *Table)
  {
    long   i;
    double max_x, min_x;
    double row=1;
    char   monotonic=1;
    char   constantstep=1;
    double step=0;
    long n;

    if (!Table) return;
    if (!Table->rows || !Table->columns) return;
    if (Table->rows == 1) row=0; // single row
    max_x = -FLT_MAX;
    min_x =  FLT_MAX;
    n     = (row ? Table->rows : Table->columns);
    /* get min and max of first column/vector */
    for (i=0; i < n; i++)
    {
      double X;
      X = (row ? Table_Index(*Table,i  ,0)
                               : Table_Index(*Table,0, i));
      if (X < min_x) min_x = X;
      if (X > max_x) max_x = X;
    } /* for */
    
    /* test for monotonicity and constant step if the table is an XY or single vector */
    if (n > 1) {
      /* mean step */
      step = (max_x - min_x)/(n-1);
      /* now test if table is monotonic on first column, and get minimal step size */
      for (i=0; i < n-1; i++) {
        double X, diff;;
        X    = (row ? Table_Index(*Table,i  ,0)
                    : Table_Index(*Table,0,  i));
        diff = (row ? Table_Index(*Table,i+1,0)
                    : Table_Index(*Table,0,  i+1)) - X;
        if (diff && fabs(diff) < fabs(step)) step = diff;
        /* change sign ? */
        if ((max_x - min_x)*diff < 0 && monotonic)
          monotonic = 0;
      } /* end for */
      
      /* now test if steps are constant within READ_TABLE_STEPTOL */
      if(!step){
        /*means there's a disconitnuity -> not constantstep*/
        constantstep=0;
      }else if (monotonic) {
        for (i=0; i < n-1; i++) {
          double X, diff;
          X    = (row ? Table_Index(*Table,i  ,0)
              : Table_Index(*Table,0,  i));
          diff = (row ? Table_Index(*Table,i+1,0)
              : Table_Index(*Table,0,  i+1)) - X;
          if ( fabs(step)*(1+READ_TABLE_STEPTOL) < fabs(diff) ||
                fabs(diff) < fabs(step)*(1-READ_TABLE_STEPTOL) )
          { constantstep = 0; break; }
        }
      }

    }
    Table->step_x= step;
    Table->max_x = max_x;
    Table->min_x = min_x;
    Table->monotonic = monotonic;
    Table->constantstep = constantstep;
  } /* end Table_Stat */

/******************************************************************************
* t_Table *Table_Read_Array(char *File, long *blocks)
*   ACTION: read as many data blocks as available, iteratively from file
*   return: initialized t_Table array, last element is an empty Table.
*           the number of extracted blocks in non NULL pointer *blocks
*******************************************************************************/
  t_Table *Table_Read_Array(char *File, long *blocks)
  {
    t_Table *Table_Array=NULL;
    long offset=0;
    long block_number=0;
    long allocated=256;
    long nelements=1;

    /* first allocate an initial empty t_Table array */
    Table_Array = (t_Table *)malloc(allocated*sizeof(t_Table));
    if (!Table_Array) {
      fprintf(stderr, "Error: Can not allocate memory %li (Table_Read_Array).\n",
         allocated*sizeof(t_Table));
      *blocks = 0;
      return (NULL);
    }

    while (nelements > 0)
    {
      t_Table Table;

      /* if ok, set t_Table block number else exit loop */
      block_number++;
      Table.block_number = block_number;
      
      /* access file at offset and get following block. Block number is from the set offset
       * hence the hardcoded 1 - i.e. the next block counted from offset.*/
      nelements = Table_Read_Offset(&Table, File, 1, &offset,0);
      /*if the block is empty - don't store it*/
      if (nelements>0){
          /* if t_Table array is not long enough, expand and realocate */
          if (block_number >= allocated-1) {
              allocated += 256;
              Table_Array = (t_Table *)realloc(Table_Array,
                      allocated*sizeof(t_Table));
              if (!Table_Array) {
                  fprintf(stderr, "Error: Can not re-allocate memory %li (Table_Read_Array).\n",
                          allocated*sizeof(t_Table));
                  *blocks = 0;
                  return (NULL);
              }
          }
          /* store it into t_Table array */
          //snprintf(Table.filename, 1024, "%s#%li", File, block_number-1);
          Table_Array[block_number-1] = Table;
      }
      /* continues until we find an empty block */
    }
    /* send back number of extracted blocks */
    if (blocks) *blocks = block_number-1;

    /* now store total number of elements in Table array */
    for (offset=0; offset < block_number;
      Table_Array[offset++].array_length = block_number-1);

    return(Table_Array);
  } /* end Table_Read_Array */
/*******************************************************************************
* void Table_Free_Array(t_Table *Table)
*   ACTION: free a Table array
*******************************************************************************/
  void Table_Free_Array(t_Table *Table)
  {
    long index;
    if (!Table) return;
    for (index=0;index < Table[0].array_length; index++){
            Table_Free(&Table[index]);
    }
    free(Table);
  } /* end Table_Free_Array */

/******************************************************************************
* long Table_Info_Array(t_Table *Table)
*    ACTION: print informations about a Table array
*    return: number of elements in the Table array
*******************************************************************************/
  long Table_Info_Array(t_Table *Table)
  {
    long index=0;

    if (!Table) return(-1);
    while (index < Table[index].array_length
       && (Table[index].data || Table[index].header)
       && (Table[index].rows*Table[index].columns) ) {
      Table_Info(Table[index]);
      index++;
    }
    printf("This Table array contains %li elements\n", index);
    return(index);
  } /* end Table_Info_Array */

/******************************************************************************
* char **Table_ParseHeader(char *header, symbol1, symbol2, ..., NULL)
*    ACTION: search for char* symbols in header and return their value or NULL
*            the search is not case sensitive.
*            Last argument MUST be NULL
*    return: array of char* with line following each symbol, or NULL if not found
*******************************************************************************/
#ifndef MyNL_ARGMAX
#define MyNL_ARGMAX 50
#endif

char **Table_ParseHeader_backend(char *header, ...){
  va_list ap;
  char exit_flag=0;
  int counter   =0;
  char **ret    =NULL;
  if (!header || header[0]=='\0') return(NULL);

  ret = (char**)calloc(MyNL_ARGMAX, sizeof(char*));
  if (!ret) {
    printf("Table_ParseHeader: Cannot allocate %i values array for Parser (Table_ParseHeader).\n",
      MyNL_ARGMAX);
    return(NULL);
  }
  for (counter=0; counter < MyNL_ARGMAX; ret[counter++] = NULL);
  counter=0;

  va_start(ap, header);
  while(!exit_flag && counter < MyNL_ARGMAX-1)
  {
    char *arg_char=NULL;
    char *pos     =NULL;
    /* get variable argument value as a char */
    arg_char = va_arg(ap, char *);
    if (!arg_char || arg_char[0]=='\0'){
      exit_flag = 1; break;
    }
    /* search for the symbol in the header */
    pos = (char*)strcasestr(header, arg_char);
    if (pos) {
      char *eol_pos;
      eol_pos = strchr(pos+strlen(arg_char), '\n');
      if (!eol_pos)
        eol_pos = strchr(pos+strlen(arg_char), '\r');
      if (!eol_pos)
        eol_pos = pos+strlen(pos)-1;
      ret[counter] = (char*)malloc(eol_pos - pos);
      if (!ret[counter]) {
        printf("Table_ParseHeader: Cannot allocate value[%i] array for Parser searching for %s (Table_ParseHeader).\n",
          counter, arg_char);
        exit_flag = 1; break;
      }
      strncpy(ret[counter], pos+strlen(arg_char), eol_pos - pos - strlen(arg_char));
      ret[counter][eol_pos - pos - strlen(arg_char)]='\0';
    }
    counter++;
  }
  va_end(ap);
  return(ret);
} /* Table_ParseHeader */

/******************************************************************************
* double Table_Interp1d(x, x1, y1, x2, y2)
*    ACTION: interpolates linearly at x between y1=f(x1) and y2=f(x2)
*    return: y=f(x) value
*******************************************************************************/
double Table_Interp1d(double x,
  double x1, double y1,
  double x2, double y2)
{
  double slope;
  if (x2 == x1) return (y1+y2)/2;
  if (y1 == y2) return  y1;
  slope = (y2 - y1)/(x2 - x1);
  return y1+slope*(x - x1);
} /* Table_Interp1d */

/******************************************************************************
* double Table_Interp1d_nearest(x, x1, y1, x2, y2)
*    ACTION: table lookup with nearest method at x between y1=f(x1) and y2=f(x2)
*    return: y=f(x) value
*******************************************************************************/
double Table_Interp1d_nearest(double x,
  double x1, double y1,
  double x2, double y2)
{
  if (fabs(x-x1) < fabs(x-x2)) return (y1);
  else return(y2);
} /* Table_Interp1d_nearest */

/******************************************************************************
* double Table_Interp2d(x,y, x1,y1, x2,y2, z11,z12,z21,z22)
*    ACTION: interpolates bi-linearly at (x,y) between z1=f(x1,y1) and z2=f(x2,y2)
*    return: z=f(x,y) value
*    x,y |   x1   x2
*    ----------------
*     y1 |   z11  z21
*     y2 |   z12  z22
*******************************************************************************/
double Table_Interp2d(double x, double y,
  double x1, double y1,
  double x2, double y2,
  double z11, double z12, double z21, double z22)
{
  double ratio_x, ratio_y;
  if (x2 == x1) return Table_Interp1d(y, y1,z11, y2,z12);
  if (y1 == y2) return Table_Interp1d(x, x1,z11, x2,z21);

  ratio_y = (y - y1)/(y2 - y1);
  ratio_x = (x - x1)/(x2 - x1);
  return (1-ratio_x)*(1-ratio_y)*z11 + ratio_x*(1-ratio_y)*z21
    + ratio_x*ratio_y*z22         + (1-ratio_x)*ratio_y*z12;
} /* Table_Interp2d */

/* end of read_table-lib.c */


#ifndef REF_LIB_H
#define REF_LIB_H "$Revision$"

void StdReflecFunc(double, double*, double*);
void TableReflecFunc(double, t_Table*, double*);

#endif

/* end of ref-lib.h */
/****************************************************************************
*
* McStas, neutron ray-tracing package
*         Copyright 1997-2006, All rights reserved
*         Risoe National Laboratory, Roskilde, Denmark
*         Institut Laue Langevin, Grenoble, France
*
* Library: share/ref-lib.c
*
* %Identification
* Written by: Peter Christiansen
* Date: August, 2006
* Origin: RISOE
* Release: McStas 1.10
* Version: $Revision$
*
* Commonly used reflection functions are declared in this file which
* are used by some guide and mirror components.
*
* Variable names have prefix 'mc_ref_' for 'McStas Reflection' 
* to avoid conflicts
*
* Usage: within SHARE
* %include "ref-lib"
*
****************************************************************************/

#ifndef REF_LIB_H
#include "ref-lib.h"
#endif

#ifndef READ_TABLE_LIB_H
#include "read_table-lib.h"
#include "read_table-lib.c"
#endif

/****************************************************************************
* void StdReflecFunc(double q, double *par, double *r)
* 
* The McStas standard analytic parametrization of the reflectivity.
* The parameters are:
* R0:      [1]    Low-angle reflectivity
* Qc:      [AA-1] Critical scattering vector
* alpha:   [AA]   Slope of reflectivity
* m:       [1]    m-value of material. Zero means completely absorbing.
* W:       [AA-1] Width of supermirror cut-off
*****************************************************************************/
void StdReflecFunc(double mc_pol_q, double *mc_pol_par, double *mc_pol_r) {
    double R0    = mc_pol_par[0];
    double Qc    = mc_pol_par[1];
    double alpha = mc_pol_par[2];
    double m     = mc_pol_par[3];
    double W     = mc_pol_par[4];
    double beta  = 0;
    mc_pol_q     = fabs(mc_pol_q);
    double arg;
        
    /* Simpler parametrization from Henrik Jacobsen uses these values that depend on m only.
       double m_value=m*0.9853+0.1978;
       double W=-0.0002*m_value+0.0022;
       double alpha=0.2304*m_value+5.0944;
       double beta=-7.6251*m_value+68.1137; 
       If W and alpha are set to 0, use Henrik's approach for estimating these parameters
       and apply the formulation:
       arg = R0*0.5*(1-tanh(arg))*(1-alpha*(q-Qc)+beta*(q-Qc)*(q-Qc));
    */  
    if (W==0 && alpha==0) {
      m=m*0.9853+0.1978;
      W=-0.0002*m+0.0022;
      alpha=0.2304*m+5.0944;
      beta=-7.6251*m+68.1137;
      if (m<=3) {
	alpha=m;
	beta=0;
      }
    }
    
    arg = W > 0 ? (mc_pol_q - m*Qc)/W : 11;

    if (arg > 10 || m <= 0 || Qc <=0 || R0 <= 0) {
      *mc_pol_r = 0;
      return;
    }
    
    if (m < 1) { Qc *= m; m=1; }
    
    if(mc_pol_q <= Qc) {      
      *mc_pol_r = R0;
      return;
    }
    
    
    *mc_pol_r = R0*0.5*(1 - tanh(arg))*(1 - alpha*(mc_pol_q - Qc) + beta*(mc_pol_q - Qc)*(mc_pol_q - Qc));
    
    return;
  }

/****************************************************************************
* void TableReflecFunc(double q, t_Table *par, double *r) {
* 
* Looks up the reflectivity in a table using the routines in read_table-lib.
*****************************************************************************/
void TableReflecFunc(double mc_pol_q, t_Table *mc_pol_par, double *mc_pol_r) {
    
  *mc_pol_r = Table_Value(*mc_pol_par, mc_pol_q, 1);
  if(*mc_pol_r>1)
    *mc_pol_r = 1;
  return;
}

/* end of ref-lib.c */


#define NUMBERMIRROR 20
#define V2Q_conic 1.58825361e-3
#define Q2V_conic 629.622368
#define DEG2RAD 3.1415927 / 180

///////////////////////////////////////////////////////////////////////////
/////////////// Some Structures
///////////////////////////////////////////////////////////////////////////

		//! Stucture to represent one Branch of the spiral
typedef struct _LogSpir{
	double zmin;
	double zmax;
	double ymin;
	double ymax;
	double k;			//=cotan(psi_rad);
	double precision;	// which precision should be reached by the various newton methods
	double theta_end;	// angle (rad) under which the endpoint of the spiral is seen
	double mValue;		// m-value of the mirror coatings
	double mindistance; // minimum distace between two mirrors, useful for the minimum time
	double phi_rot;		// the rotation of the spiral branch with respect to the optical axis
	int mirrored;		// whether the spiral is mirrored the optical axis
	} LogSpir;

	//! Structure to collect all the LogSpiral branches and count them
typedef struct _SceneLog{
	int number_spirals;
	LogSpir all_branches[NUMBERMIRROR];
} SceneLog;

	//! Structure resembling a neutron with position direction and spin (sometimes misused as a vector)
typedef struct _part_log{
	double customz;
	double customy;
	double customx;
	double customvz;
	double customvy;
	double customvx;
	double customsx;
	double customsy;
	double customsz;
} part_log; // Structure for a 2D part_log

typedef struct _BranchTime{
	LogSpir logspir;  // the logspiral branch the neutron interacts with
	double t;		  // time of interaction
	double theta_int; // angle under which the interaction takes place
} BranchTime;		  // Structure for a 2D part_log

	double getRandomLog();
	part_log Neutron2Dinit(part_log * neutron, double z, double y, double x, double vz, double vy, double vx);
	void propagate_neutron(part_log * incoming, double dt);
	double calcSupermirrorReflectivityLog(double q, double m);
	void rotate_vector(part_log * neutron, part_log initneut, double theta_rot);
	void mirror_vector(part_log * neutron, part_log initneut, double invert);
	void reflect_neutron(part_log * neutron, part_log normal, double mValue);
	double newton_theta_end(LogSpir logspir, double theta);
	double newton_theta_end_derivative(LogSpir logspir, double theta);
	double return_r(double theta, double k, double zmin);
	part_log return_cart_coordsspiral(LogSpir logspir, double theta);
	double return_approx_theta_end(LogSpir logspir);
	double return_precise_theta_end(LogSpir logspir, int max_iterations);
	part_log return_normal_vec(LogSpir logspir, double theta);
	void LogSpirinit(LogSpir * logspir, double zmin, double zmax, double ymin, double ymax,
					 double psi, double phi_rot, double mValue, double precision, double max_iterations, int mirrored);
	void initialize_scene(SceneLog * s);
	void add_logspir(SceneLog * s, LogSpir logspir);
	double return_approx_theta_int(LogSpir logspir, part_log neutron);
	double return_precise_theta_int(LogSpir logspir, part_log n, double max_iterations);
	float return_intersection_time(LogSpir logspir, part_log neutron, double theta_int);
	BranchTime return_first_interaction(SceneLog s, part_log init_neut);
	int evaluate_first_interaction(SceneLog s, part_log * neutron);
	void run_scene(SceneLog s, part_log * n, int max_interactions);
	void populate_scene(SceneLog * s, int branches, int doublesided,
						float zmin, float zmax, float ymin, float ymax, float psi, float phi_rot, float m);

	double getRandomLog()
	{
		return (double)lrand48() / RAND_MAX;
	}

	///////////////////////////////////////////////////////////////////////////
	/////////////// Some auxiliary functions
	///////////////////////////////////////////////////////////////////////////

	/*! \brief Function to create a Neutron part_log
	@param *neutron pointer to neutron, which is initialized
	@param z z-cooridinate of neutron (m)
	@param y
	@param x
	@param vz speed in z-direction (m/s)
	@param vy
	@param vx
	@return part_log Neutron
	*/
	part_log Neutron2Dinit(part_log * neutron, double z, double y, double x, double vz, double vy, double vx){ // creating the neutron TODO polarization
		neutron->customz = z;
		neutron->customx = x;
		neutron->customvz = vz;
		neutron->customvx = vx;
		neutron->customy = y;
		neutron->customvy = vy;
		return *neutron;
	}

	/*! \brief Function to propagate the neutron in a straight line; McStas Builtin?

	@param *incoming pointer to incoming neutron, which is propagated
	@param dt (s) time increment, in which the neutron is propagated
	*/
	void propagate_neutron(part_log * incoming, double dt){
		incoming->customx += incoming->customvx * dt;
		incoming->customz += incoming->customvz * dt;
		incoming->customy += incoming->customvy * dt;
	}

	/*! \brief Function returns the supermirrorreflectivity for a given perpendicular impulse q and m-value of the mirror coating

	@param q perp impulse component
	@param m Value of the mirror
	@return Probability of reflection from 0 to 1
	*/
	double calcSupermirrorReflectivityLog(double q, double m){
		double arg;
		double R_0 = 0.995;
		double Q_c = 0.0218;
		double beta = 0;
		double alpha = 2.5;
		double W = 0.004;
		double weight = 1.0; // neutron weight to be transformed
		q = fabs(q);
		if (m >= 10){ // above m=10 we just assume perfect reflectivity for the mirror
			return weight;
		}
		if (W == 0 && alpha == 0)
		{ // approximation for relfectivity curve
			m = m * 0.9853 + 0.1978;
			W = -0.0002 * m + 0.0022;
			alpha = 0.2304 * m + 5.0944;
			beta = -7.6251 * m + 68.1137;
			if (m <= 3){
				alpha = m;
				beta = 0;
			}
		}
		arg = W > 0 ? (q - m * Q_c) / W : 11;
		if (arg > 10 || m <= 0 || Q_c <= 0 || R_0 <= 0){
			weight = 0.0;
			return weight;
		}

		if (m < 1){
			Q_c *= m;
			m = 1;
		}

		if (q <= Q_c){
			weight = R_0;
			return weight;
		}
		weight = R_0 * 0.5 * (1 - tanh(arg)) * (1 - alpha * (q - Q_c) + beta * (q - Q_c) * (q - Q_c));
		return weight;
	}
	/*! \brief Function rotates initneut by theta_rot around (0, 1, 0) and stores the result in *neutron
	@param *neutron pointer to neutron in which the result is stored
	@param initneut neutron which is rotated
	@param theta_rot (rad) angle by which to rotate initneut
	*/
	void rotate_vector(part_log * neutron, part_log initneut, double theta_rot){ // rotating the initneut by thetarot and saving the result in neutron
		double sint = sin(theta_rot);
		double cost = cos(theta_rot);
		double new_z = cost * initneut.customz - sint * initneut.customx;
		double new_x = sint * initneut.customz + cost * initneut.customx;
		double new_vz = cost * initneut.customvz - sint * initneut.customvx;
		double new_vx = sint * initneut.customvz + cost * initneut.customvx;
		neutron->customz = new_z;
		neutron->customx = new_x;
		neutron->customvz = new_vz;
		neutron->customvx = new_vx;
	}

	/*! \brief Function rotates initneut by 180 deg around (0, 0, 1)/(mirror zy-plane) and stores the result in *neutron
	@param *neutron pointer to neutron in which the result is stored
	@param initneut neutron which is rotated/mirrored
	@param invert; if +1 does nothing to the vector, if -1 inverts x component of vector
	*/
	void mirror_vector(part_log * neutron, part_log initneut, double invert){
		neutron->customz = initneut.customz;
		neutron->customx = initneut.customx * invert;
		neutron->customvz = initneut.customvz;
		neutron->customvx = initneut.customvx * invert;
		neutron->customvy = initneut.customvy;
		neutron->customy = initneut.customy;
	}

	/*! \brief Function determines probability of reflection at surface and changes velocities accordingly
	@param *neutron the neutron which to rotate
	@param normal normalvector of surface (only vx, vz and vy are used as directions)
	@param mValue m-value of the surface
	*/
	void reflect_neutron(part_log * neutron, part_log normal, double mValue){
		double vz = neutron->customvz;
		double vx = neutron->customvx;
		double vdotn = vz * normal.customvz + vx * normal.customvx;
		double weight = calcSupermirrorReflectivityLog(V2Q_conic * 2 * vdotn, mValue);
		if (rand01() <= weight){
			neutron->customvx = vx - 2 * vdotn * normal.customvx;
			neutron->customvz = vz - 2 * vdotn * normal.customvz;
		}
		else{
			; // if no reflection takes place at the mirror we dont have to change the direction
		}
	}

	/*! \brief Function for Newton Raphson determination of angle theta at the end of the spiral
		zmax
	@param logspir Logarithmic spiral
	@param theta (rad) angle theta
	@return function of which to find a root for
	*/
	double newton_theta_end(LogSpir logspir, double theta){
		return cos(theta) * logspir.zmin * exp(logspir.k * theta) - logspir.zmax;
	}
	/*! \brief Function for Newton Raphson
	@param logspir Logarithmic spiral
	@param theta (rad) angle theta
	@return derivative of function
	*/
	double newton_theta_end_derivative(LogSpir logspir, double theta){
		return logspir.zmin * exp(logspir.k * theta) * (cos(theta) * logspir.k - sin(theta));
	}

	/*! \brief Function returns r = zmin*exp(theta*cotan(psi))

	@param logspir logarithmic spiral struct
	@param theta angle under which the point is seen
	@return r coordinate
	*/
	double return_r(double theta, double k, double zmin){
		return zmin * exp(k * theta);
	};

	/*! \brief Function returns carthesian Coordinates of a point on the logspiral
	@param logspir Logarithmic spiral
	@param theta angle under which the point is seen
	@return part_log with z and x corrdinate corresponding to z and x coordinates of the spiral
	*/
	part_log return_cart_coordsspiral(LogSpir logspir, double theta){
		part_log n;
		double r = return_r(theta, logspir.k, logspir.zmin);
		n.customz = cos(theta) * r;
		n.customx = sin(theta) * r;
		return n;
	}

	/*! \brief Function returns an approximation of the angle of the end of the spiral to be refined by Newton-Raphson
	@param logspir Logarithmic spiral
	@return approximation of theta_end
	*/
	double return_approx_theta_end(LogSpir logspir){
		return log(logspir.zmax / logspir.zmin) / logspir.k;
	};

	/*! \brief Function returns the precise angle of the end of the spiral by Newton-Raphson, needed for the alignment of spirals
	@param logspir Logarithmic spiral
	@param max_iterations number of iterations after which the algorithms gives up
	@return precise angle if convergence is reached, -10 else
	*/
	double return_precise_theta_end(LogSpir logspir, int max_iterations){
		double theta_0 = return_approx_theta_end(logspir);
		double theta_n;
		for (int ii; ii < max_iterations; ii++){
			// printf("theta0=%f",theta_0);
			theta_n = theta_0 - newton_theta_end(logspir, theta_0) / newton_theta_end_derivative(logspir, theta_0);
			if (fabs(theta_0 - theta_n) < logspir.precision){
				return theta_n;
			}
			theta_0 = theta_n;
		}
		return -10.0;
	}

	/*! \brief Function returns normal vector to the logspir
	@param logspir Logarithmic spiral
	@param theta (rad) angle of the point on the spiral at which the normal vector is to be determined
	*/
	part_log return_normal_vec(LogSpir logspir, double theta){ // returns the normalized normal vector to the surface
		part_log n2d;
		double prefac = 1 / sqrt(1 + logspir.k * logspir.k);
		n2d.customvz = (cos(theta) + logspir.k * sin(theta)) * prefac;
		n2d.customvx = (sin(theta) - logspir.k * cos(theta)) * prefac;
		return n2d;
	}

	/*! \brief Function initializes a Logarithmic Spiral with all information necessary to calculate reflection
	@param logspir pointer of LogSpir object which to initialize
	@param zmin z-coordinate (m) at which the spiral starts
	@param zmax z-coordinate (m) at which the spiral ends
	@param ymin minimum y-coordinate (m) at which the mirror is realized
	@param ymax maximum y-coordinate (m) at which the mirror is realized
	@param psi angle (deg) under which the spiral is hit by neutrons originating from the origin
	@param phi_rot (rad) angle under which the individual spiral branches are rotated by
	@param mValue m-value of the reflecting surfaces
	@param precision (1) precision of the many Newton-Raphson implementations
	@param max_iterations (1) number of the maximum iterations for Newton-Raphson
	@param branches number of spiral branches
	*/
	void LogSpirinit(LogSpir * logspir, double zmin, double zmax, double ymin, double ymax,
		double psi, double phi_rot, double mValue, double precision, double max_iterations, int mirrored){
		logspir->zmin = zmin;
		logspir->zmax = zmax;
		logspir->ymin = ymin;
		logspir->ymax = ymax;
		logspir->mValue = mValue;
		logspir->k = 1 / tan(psi * DEG2RAD);
		logspir->precision = precision;
		logspir->theta_end = return_precise_theta_end(*logspir, max_iterations);
		logspir->phi_rot = phi_rot;
		logspir->mirrored = mirrored;
		logspir->mindistance = 2 * sin(logspir->theta_end / 2) * zmin;
	}
	/*! \brief Function initializes a scene to hold all the LogSpiral objects
	@param s pointer to the to initialize scene
	*/
	void initialize_scene(SceneLog * s){
		s->number_spirals = 0;
	}

	/*! \brief Function adding a single LogSpiral object to the scene and incrementing the counter
	@param s pointer to the to initialize scene
	*/
	void add_logspir(SceneLog * s, LogSpir logspir){
		s->all_branches[s->number_spirals] = logspir;
		s->number_spirals += 1;
	}

	///////////////////////////////////////////////////////////////////////////
	/////////////// Functions assisting in the reflection
	///////////////////////////////////////////////////////////////////////////

	/*! \brief Function returns the approximate angle of intersection spiral to be refined by Newton-Raphson
	@param logspir Logarithmic spiral
	*/
	double return_approx_theta_int(LogSpir logspir, part_log neutron){ // approximating by intersection with line seems comp expensive + if the intersection is far away
		part_log *mirror_approx;
		double denom;
		double lam;
		double z_int;
		double x_int;
		double m = logspir.zmax * tan(logspir.theta_end) / (logspir.zmax - logspir.zmin);
		denom = neutron.customvx - neutron.customvz * m;
		if (denom == 0){
			return logspir.theta_end / 2; // neutron parallel to the spiral arm
		}
		else{
			lam = ((-logspir.zmin + neutron.customz) * m - neutron.customx) / denom; // approximate mirror as a line calculate the intersection and determine the angle
			z_int = neutron.customz + neutron.customvz * lam;
			// printf("z intersection %f x intersection %f lam %f theta_end %f\n", z_int, x_int, lam, logspir.theta_end);
			if (z_int >= logspir.zmin && z_int <= logspir.zmax){
				x_int = neutron.customx + neutron.customvx * lam;
				// printf("init theta %f theta_end/2 %f\n", atan((x_int)/(z_int)), logspir.theta_end/2);
				return atan((x_int) / (z_int));
			}
		}
		return -10; // if the intersection is not on the line return the False Value
	}

	/*! \brief Function returns the precise angle of intersection by Newton-Raphson
	@param logspir Logarithmic spiral
	@param n neutron impinging on the mirror branch
	@param max_iterations maximum numbers of iterations after which Newton gives up
	@return angle of intersection if convergence is reached -10 else
	*/
	double return_precise_theta_int(LogSpir logspir, part_log n, double max_iterations){ //
		double theta_0, theta_n;
		double m = n.customvx / n.customvz;
		double x0 = n.customx - n.customz / n.customvz * n.customvx;
		double f(double theta){ // function which to minimize via Newton-Raphson
			return logspir.zmin * exp(logspir.k * theta) * (sin(theta) - m * cos(theta)) - x0;
		}
		double f_derivative(double theta){ // derivative of the function for Newton-Raphson
			return logspir.zmin * exp(logspir.k * theta) * (cos(theta) * (1 - logspir.k * m) + sin(theta) * (logspir.k + m));
		}
		theta_0 = return_approx_theta_int(logspir, n);
		if (theta_0 < -9){
			return -10;
		}
		for (int ii; ii < max_iterations; ii++){
			theta_n = theta_0 - f(theta_0) / f_derivative(theta_0);
			if (ii > 1 && (theta_n < 0 || theta_n > logspir.theta_end)){
				return -10;
			}
			if (fabs(theta_n - theta_0) < logspir.precision){
				if ((0 < theta_n) && (theta_n < logspir.theta_end)){ // if the intersection theta is outside the bounds of the spiral--> no interaction
					return theta_n;
				}
				// printf("no htis\n");
				return -10;
			}
			theta_0 = theta_n;
		}
		return -10;
	}

	/*! \brief Function returns the time until interaction, only if interaction is possible, i.e., neutron hits mirror (ycoord)
	@param logspir Logarithmic spiral
	@param neutron intersecting neutron
	@param theta_int angle of intersection
	@return time (s) if the intersection takes place, -10 else
	*/
	float return_intersection_time(LogSpir logspir, part_log neutron, double theta_int){
		part_log int_coord;
		double t;
		double y;
		if (theta_int > 0){
			int_coord = return_cart_coordsspiral(logspir, theta_int); // cart coords of the interaction
			////printf("z of spiral %f z0 = %f\n", int_coord.z, z0);
			t = (int_coord.customz - neutron.customz) / neutron.customvz;
			y = neutron.customy + neutron.customvy * t;
			if (t * t > 0.25 * logspir.mindistance * logspir.mindistance /
							(neutron.customvz * neutron.customvz + neutron.customvx * neutron.customvx) &&
				y >= logspir.ymin && y <= logspir.ymax){ // minimum distance between two mirrors (sensible distance, smaller makes no sense)
				return t;
			}
		}
		return -10.0;
	}

	/*! \brief Function returns the BranchTime object of the first Interaction
	@param logspir Logarithmic spiral
	@param part_log incoming neutron
	*/
	BranchTime return_first_interaction(SceneLog s, part_log init_neut){
		double t;
		double theta_rot;
		double theta_int;
		int mirrored;
		BranchTime brt; // stores the first interaction and saves the according branch
		part_log rotneut;
		brt.t = -1.0;
		for (int kk = 0; kk < s.number_spirals; kk++){ // go through all the orientations
			theta_rot = -s.all_branches[kk].phi_rot;
			mirrored = s.all_branches[kk].mirrored;
			// rotate and mirror the neutron according to the branches orientation
			mirror_vector(&rotneut, init_neut, mirrored);
			rotate_vector(&rotneut, rotneut, theta_rot);
			theta_int = return_precise_theta_int(s.all_branches[kk], rotneut, 10);
			// printf("thetaint %f vz berfore %f\n", theta_int, init_neut.vz);
			t = return_intersection_time(s.all_branches[kk], rotneut, theta_int);
			if (t > 0){
				if (t < brt.t || brt.t < 0){
					brt.t = t;
					brt.logspir = s.all_branches[kk];
					brt.theta_int = theta_int;
				}
			}
		}
		return brt;
	}

	int evaluate_first_interaction(SceneLog s, part_log * neutron){
		BranchTime brt;
		brt = return_first_interaction(s, *neutron);
		double t_prop;
		double theta_int;
		double phi_rot;
		int mirrored;
		part_log n;
		part_log n_rot;
		t_prop = brt.t;
		theta_int = brt.theta_int;		 // the angle (unrotated) in which the interaction takes place
		phi_rot = brt.logspir.phi_rot;	 // the angle between the the logspir and the McStas Co-ordinate system
		mirrored = brt.logspir.mirrored; // is the spiral mirrored?
		// printf("proptime %f", t_prop);
		if (t_prop < 0){
			return 0;
		}
		else{
			// if the time is valid we propagate the neutron to the surface
			propagate_neutron(neutron, t_prop); // propagate neutron
			n = return_normal_vec(brt.logspir, theta_int);
			rotate_vector(&n_rot, n, phi_rot); //
			mirror_vector(&n_rot, n_rot, mirrored);
			reflect_neutron(neutron, n_rot, brt.logspir.mValue);
			return 1;
		}
	}

	void run_scene(SceneLog s, part_log * n, int max_interactions){
		int interaction = 1;
		for (int kk; kk < max_interactions; kk++){
			interaction = evaluate_first_interaction(s, n);
			if (interaction == 0){
				break;
			}
		}
	}

	/*! \brief Function to populate the scene with mirrors, this can be changed if one wants to use a non standard geometry
	@param *s pointer to Scene to populate
	@param branches number of branches per side
	@param doublesided spirals on both sides of the optical axis
	@param zmin zmin of all spirals
	@param zmax zmax of all spirals
	@param ymin
	@param ymax
	@param psi (deg) angle of attack for all mirrors
	@param phi_rot (rad) angle between the mirrors
	@param m m-value of coating of all mirrors
	*/
	void populate_scene(SceneLog * s, int branches, int doublesided, float zmin, float zmax,
		float ymin, float ymax, float psi, float phi_rot, float m){
		LogSpir logspir;
		LogSpirinit(&logspir, zmin, zmax, ymin, ymax, psi, 0, m, 1e-7, 10, 1);
		double phi_rot_base = phi_rot == 0 ? logspir.theta_end : phi_rot;
		// LogSpirinit(logspir, zmin, zmax, ymin, ymax, psi, phi, m,  1e-7, 10, mirrored);
		initialize_scene(s);
		for (int kk = 0; kk < ((doublesided + 1) * branches); kk++){ // TODO reading orientations from a file?
			logspir.phi_rot = kk < branches ? kk * phi_rot_base : (kk - branches) * phi_rot_base;
			logspir.mirrored = kk < branches ? 1 : -1;
			add_logspir(s, logspir);
			printf("phi %f, inv %d total_spir %d\n", s->all_branches[kk].phi_rot, s->all_branches[kk].mirrored, s->number_spirals);
		}
	}

	///////////////////////////////////////////////////////////////////////////
	/////////////// End of functions
	///////////////////////////////////////////////////////////////////////////
	// part_log n;
#line 9393 "./reverse_test.c"

/* Instrument parameters. */
MCNUM mcipsource_width;
MCNUM mcipsource_divergence;
MCNUM mcipL_source;
MCNUM mcipdL;
MCNUM mcipflux;
MCNUM mcipzmin;
MCNUM mcipzmax;
MCNUM mcipymin;
MCNUM mcipymax;
MCNUM mcippsi;
MCNUM mcipphi_rot;
MCNUM mcipprecision;
MCNUM mcipmax_iterations;
MCNUM mcipmValue;
MCNUM mcipbranches;
MCNUM mcipdoublesided;
MCNUM mcipplaceholder;

#define mcNUMIPAR 17
int mcnumipar = 17;
struct mcinputtable_struct mcinputtable[mcNUMIPAR+1] = {
  "source_width", &mcipsource_width, instr_type_double, "0.000001", 
  "source_divergence", &mcipsource_divergence, instr_type_double, "3.7116398665896515", 
  "L_source", &mcipL_source, instr_type_double, "5", 
  "dL", &mcipdL, instr_type_double, "2", 
  "flux", &mcipflux, instr_type_double, "1", 
  "zmin", &mcipzmin, instr_type_double, "0.15", 
  "zmax", &mcipzmax, instr_type_double, "0.3", 
  "ymin", &mcipymin, instr_type_double, "-1", 
  "ymax", &mcipymax, instr_type_double, "1", 
  "psi", &mcippsi, instr_type_double, "1", 
  "phi_rot", &mcipphi_rot, instr_type_double, "0", 
  "precision", &mcipprecision, instr_type_double, "1e-7", 
  "max_iterations", &mcipmax_iterations, instr_type_double, "10", 
  "mValue", &mcipmValue, instr_type_double, "6.2", 
  "branches", &mcipbranches, instr_type_double, "5", 
  "doublesided", &mcipdoublesided, instr_type_double, "1", 
  "placeholder", &mcipplaceholder, instr_type_double, "0", 
  NULL, NULL, instr_type_double, ""
};

/* User declarations from instrument definition. */
#define mccompcurname  logspir_test
#define mccompcurtype  INSTRUMENT
#define mccompcurindex 0
#define mcposalogspir_test coords_set(0,0,0)
#define source_width mcipsource_width
#define source_divergence mcipsource_divergence
#define L_source mcipL_source
#define dL mcipdL
#define flux mcipflux
#define zmin mcipzmin
#define zmax mcipzmax
#define ymin mcipymin
#define ymax mcipymax
#define psi mcippsi
#define phi_rot mcipphi_rot
#define precision mcipprecision
#define max_iterations mcipmax_iterations
#define mValue mcipmValue
#define branches mcipbranches
#define doublesided mcipdoublesided
#define placeholder mcipplaceholder
#undef placeholder
#undef doublesided
#undef branches
#undef mValue
#undef max_iterations
#undef precision
#undef phi_rot
#undef psi
#undef ymax
#undef ymin
#undef zmax
#undef zmin
#undef flux
#undef dL
#undef L_source
#undef source_divergence
#undef source_width
#undef mcposalogspir_test
#undef mccompcurindex
#undef mccompcurtype
#undef mccompcurname

/* neutron state table at each component input (local coords) */
/* [x, y, z, vx, vy, vz, t, sx, sy, sz, p] */
MCNUM mccomp_storein[11*12];
/* Components position table (absolute and relative coords) */
Coords mccomp_posa[12];
Coords mccomp_posr[12];
/* Counter for each comp to check for inactive ones */
MCNUM  mcNCounter[12];
MCNUM  mcPCounter[12];
MCNUM  mcP2Counter[12];
#define mcNUMCOMP 11 /* number of components */
/* Counter for PROP ABSORB */
MCNUM  mcAbsorbProp[12];
/* Flag true when previous component acted on the neutron (SCATTER) */
MCNUM mcScattered=0;
/* Flag true when neutron should be restored (RESTORE) */
MCNUM mcRestore=0;
/* Declarations of component definition and setting parameters. */

/* Setting parameters for component 'origin' [1]. */
char mccorigin_profile[16384];
MCNUM mccorigin_percent;
MCNUM mccorigin_flag_save;
MCNUM mccorigin_minutes;

/* Setting parameters for component 'source_div' [3]. */
MCNUM mccsource_div_xwidth;
MCNUM mccsource_div_yheight;
MCNUM mccsource_div_focus_aw;
MCNUM mccsource_div_focus_ah;
MCNUM mccsource_div_E0;
MCNUM mccsource_div_dE;
MCNUM mccsource_div_lambda0;
MCNUM mccsource_div_dlambda;
MCNUM mccsource_div_gauss;
MCNUM mccsource_div_flux;

/* Setting parameters for component 'slit' [4]. */
MCNUM mccslit_xmin;
MCNUM mccslit_xmax;
MCNUM mccslit_ymin;
MCNUM mccslit_ymax;
MCNUM mccslit_radius;
MCNUM mccslit_xwidth;
MCNUM mccslit_yheight;

/* Setting parameters for component 'psd_before_optic' [5]. */
int mccpsd_before_optic_nx;
int mccpsd_before_optic_ny;
char mccpsd_before_optic_filename[16384];
MCNUM mccpsd_before_optic_xmin;
MCNUM mccpsd_before_optic_xmax;
MCNUM mccpsd_before_optic_ymin;
MCNUM mccpsd_before_optic_ymax;
MCNUM mccpsd_before_optic_xwidth;
MCNUM mccpsd_before_optic_yheight;
MCNUM mccpsd_before_optic_restore_neutron;

/* Definition parameters for component 'flat_ellipse_horizontal' [6]. */
#define mccflat_ellipse_horizontal_reflect "supermirror_m3.rfl" /* declared as a string. May produce warnings at compile */
/* Setting parameters for component 'flat_ellipse_horizontal' [6]. */
MCNUM mccflat_ellipse_horizontal_sourceDist;
MCNUM mccflat_ellipse_horizontal_LStart;
MCNUM mccflat_ellipse_horizontal_LEnd;
MCNUM mccflat_ellipse_horizontal_lStart;
MCNUM mccflat_ellipse_horizontal_lEnd;
MCNUM mccflat_ellipse_horizontal_r_0;
MCNUM mccflat_ellipse_horizontal_nummirror;
MCNUM mccflat_ellipse_horizontal_mf;
MCNUM mccflat_ellipse_horizontal_mb;
MCNUM mccflat_ellipse_horizontal_mirror_width;
MCNUM mccflat_ellipse_horizontal_mirror_sidelength;
MCNUM mccflat_ellipse_horizontal_doubleReflections;

/* Setting parameters for component 'slit1' [7]. */
MCNUM mccslit1_xmin;
MCNUM mccslit1_xmax;
MCNUM mccslit1_ymin;
MCNUM mccslit1_ymax;
MCNUM mccslit1_radius;
MCNUM mccslit1_xwidth;
MCNUM mccslit1_yheight;

/* Setting parameters for component 'psd_monitor_beforelog' [8]. */
int mccpsd_monitor_beforelog_nx;
int mccpsd_monitor_beforelog_ny;
char mccpsd_monitor_beforelog_filename[16384];
MCNUM mccpsd_monitor_beforelog_xmin;
MCNUM mccpsd_monitor_beforelog_xmax;
MCNUM mccpsd_monitor_beforelog_ymin;
MCNUM mccpsd_monitor_beforelog_ymax;
MCNUM mccpsd_monitor_beforelog_xwidth;
MCNUM mccpsd_monitor_beforelog_yheight;
MCNUM mccpsd_monitor_beforelog_restore_neutron;

/* Setting parameters for component 'logspir' [9]. */
MCNUM mcclogspir_zmin;
MCNUM mcclogspir_zmax;
MCNUM mcclogspir_ymin;
MCNUM mcclogspir_ymax;
MCNUM mcclogspir_psi;
MCNUM mcclogspir_phi_rot;
MCNUM mcclogspir_precision;
MCNUM mcclogspir_max_iterations;
MCNUM mcclogspir_mValue;
MCNUM mcclogspir_branches;
MCNUM mcclogspir_doublesided;
MCNUM mcclogspir_placeholder;

/* Setting parameters for component 'psd_monitor' [10]. */
int mccpsd_monitor_nx;
int mccpsd_monitor_ny;
char mccpsd_monitor_filename[16384];
MCNUM mccpsd_monitor_xmin;
MCNUM mccpsd_monitor_xmax;
MCNUM mccpsd_monitor_ymin;
MCNUM mccpsd_monitor_ymax;
MCNUM mccpsd_monitor_xwidth;
MCNUM mccpsd_monitor_yheight;
MCNUM mccpsd_monitor_restore_neutron;

/* User component declarations. */

/* User declarations for component 'origin' [1]. */
#define mccompcurname  origin
#define mccompcurtype  Progress_bar
#define mccompcurindex 1
#define IntermediateCnts mccorigin_IntermediateCnts
#define StartTime mccorigin_StartTime
#define EndTime mccorigin_EndTime
#define CurrentTime mccorigin_CurrentTime
#define profile mccorigin_profile
#define percent mccorigin_percent
#define flag_save mccorigin_flag_save
#define minutes mccorigin_minutes
#line 44 "/usr/share/mcstas/2.7/tools/Python/mcrun/../mccodelib/../../../misc/Progress_bar.comp"
#ifndef PROGRESS_BAR
#define PROGRESS_BAR
#else
#error Only one Progress_bar component may be used in an instrument definition.
#endif

double IntermediateCnts;
time_t StartTime;
time_t EndTime;
time_t CurrentTime;
#line 9627 "./reverse_test.c"
#undef minutes
#undef flag_save
#undef percent
#undef profile
#undef CurrentTime
#undef EndTime
#undef StartTime
#undef IntermediateCnts
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'source' [2]. */
#define mccompcurname  source
#define mccompcurtype  Arm
#define mccompcurindex 2
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'source_div' [3]. */
#define mccompcurname  source_div
#define mccompcurtype  Source_div
#define mccompcurindex 3
#define thetah mccsource_div_thetah
#define thetav mccsource_div_thetav
#define sigmah mccsource_div_sigmah
#define sigmav mccsource_div_sigmav
#define tan_h mccsource_div_tan_h
#define tan_v mccsource_div_tan_v
#define p_init mccsource_div_p_init
#define dist mccsource_div_dist
#define focus_xw mccsource_div_focus_xw
#define focus_yh mccsource_div_focus_yh
#define xwidth mccsource_div_xwidth
#define yheight mccsource_div_yheight
#define focus_aw mccsource_div_focus_aw
#define focus_ah mccsource_div_focus_ah
#define E0 mccsource_div_E0
#define dE mccsource_div_dE
#define lambda0 mccsource_div_lambda0
#define dlambda mccsource_div_dlambda
#define gauss mccsource_div_gauss
#define flux mccsource_div_flux
#line 69 "/usr/share/mcstas/2.7/tools/Python/mcrun/../mccodelib/../../../sources/Source_div.comp"
double thetah, thetav, sigmah, sigmav, tan_h, tan_v, p_init, dist, focus_xw, focus_yh;
#line 9674 "./reverse_test.c"
#undef flux
#undef gauss
#undef dlambda
#undef lambda0
#undef dE
#undef E0
#undef focus_ah
#undef focus_aw
#undef yheight
#undef xwidth
#undef focus_yh
#undef focus_xw
#undef dist
#undef p_init
#undef tan_v
#undef tan_h
#undef sigmav
#undef sigmah
#undef thetav
#undef thetah
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'slit' [4]. */
#define mccompcurname  slit
#define mccompcurtype  Slit
#define mccompcurindex 4
#define xmin mccslit_xmin
#define xmax mccslit_xmax
#define ymin mccslit_ymin
#define ymax mccslit_ymax
#define radius mccslit_radius
#define xwidth mccslit_xwidth
#define yheight mccslit_yheight
#undef yheight
#undef xwidth
#undef radius
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'psd_before_optic' [5]. */
#define mccompcurname  psd_before_optic
#define mccompcurtype  PSD_monitor
#define mccompcurindex 5
#define PSD_N mccpsd_before_optic_PSD_N
#define PSD_p mccpsd_before_optic_PSD_p
#define PSD_p2 mccpsd_before_optic_PSD_p2
#define nx mccpsd_before_optic_nx
#define ny mccpsd_before_optic_ny
#define filename mccpsd_before_optic_filename
#define xmin mccpsd_before_optic_xmin
#define xmax mccpsd_before_optic_xmax
#define ymin mccpsd_before_optic_ymin
#define ymax mccpsd_before_optic_ymax
#define xwidth mccpsd_before_optic_xwidth
#define yheight mccpsd_before_optic_yheight
#define restore_neutron mccpsd_before_optic_restore_neutron
#line 62 "/usr/share/mcstas/2.7/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
  DArray2d PSD_N;
  DArray2d PSD_p;
  DArray2d PSD_p2;
#line 9742 "./reverse_test.c"
#undef restore_neutron
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef filename
#undef ny
#undef nx
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'flat_ellipse_horizontal' [6]. */
#define mccompcurname  flat_ellipse_horizontal
#define mccompcurtype  FlatEllipse_finite_mirror
#define mccompcurindex 6
#define reflect mccflat_ellipse_horizontal_reflect
#define s mccflat_ellipse_horizontal_s
#define pTable mccflat_ellipse_horizontal_pTable
#define R0 mccflat_ellipse_horizontal_R0
#define Qc mccflat_ellipse_horizontal_Qc
#define W mccflat_ellipse_horizontal_W
#define alpha mccflat_ellipse_horizontal_alpha
#define transmit mccflat_ellipse_horizontal_transmit
#define sourceDist mccflat_ellipse_horizontal_sourceDist
#define LStart mccflat_ellipse_horizontal_LStart
#define LEnd mccflat_ellipse_horizontal_LEnd
#define lStart mccflat_ellipse_horizontal_lStart
#define lEnd mccflat_ellipse_horizontal_lEnd
#define r_0 mccflat_ellipse_horizontal_r_0
#define nummirror mccflat_ellipse_horizontal_nummirror
#define mf mccflat_ellipse_horizontal_mf
#define mb mccflat_ellipse_horizontal_mb
#define mirror_width mccflat_ellipse_horizontal_mirror_width
#define mirror_sidelength mccflat_ellipse_horizontal_mirror_sidelength
#define doubleReflections mccflat_ellipse_horizontal_doubleReflections
#line 62 "FlatEllipse_finite_mirror.comp"
    //Scene where all geometry is added to
    Scene s;
    //Variables for Reflectivity from McStas Tables
    //t_Table pTable;
    //double R0 = 0.99;
    //double Qc = 0.021;
    //double W = 0.003;
    //double alpha = 6.07;
    //double transmit = 0;
    double *pointer_lStart;// lStart has to be used in Trace later, this requires a pointer
    //Function to handle Conic-Neutron collisions with reflectivity from McStas Tables
    void traceNeutronConicWithTables(Particle* p, ConicSurf c);
    double *rfront_inner;
    
    double dt;
    int silicon; // +1: neutron in silicon, -1: neutron in air, 0: mirrorwidth is 0; neutron cannot be in silicon
#line 9801 "./reverse_test.c"
#undef doubleReflections
#undef mirror_sidelength
#undef mirror_width
#undef mb
#undef mf
#undef nummirror
#undef r_0
#undef lEnd
#undef lStart
#undef LEnd
#undef LStart
#undef sourceDist
#undef transmit
#undef alpha
#undef W
#undef Qc
#undef R0
#undef pTable
#undef s
#undef reflect
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'slit1' [7]. */
#define mccompcurname  slit1
#define mccompcurtype  Slit
#define mccompcurindex 7
#define xmin mccslit1_xmin
#define xmax mccslit1_xmax
#define ymin mccslit1_ymin
#define ymax mccslit1_ymax
#define radius mccslit1_radius
#define xwidth mccslit1_xwidth
#define yheight mccslit1_yheight
#undef yheight
#undef xwidth
#undef radius
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'psd_monitor_beforelog' [8]. */
#define mccompcurname  psd_monitor_beforelog
#define mccompcurtype  PSD_monitor
#define mccompcurindex 8
#define PSD_N mccpsd_monitor_beforelog_PSD_N
#define PSD_p mccpsd_monitor_beforelog_PSD_p
#define PSD_p2 mccpsd_monitor_beforelog_PSD_p2
#define nx mccpsd_monitor_beforelog_nx
#define ny mccpsd_monitor_beforelog_ny
#define filename mccpsd_monitor_beforelog_filename
#define xmin mccpsd_monitor_beforelog_xmin
#define xmax mccpsd_monitor_beforelog_xmax
#define ymin mccpsd_monitor_beforelog_ymin
#define ymax mccpsd_monitor_beforelog_ymax
#define xwidth mccpsd_monitor_beforelog_xwidth
#define yheight mccpsd_monitor_beforelog_yheight
#define restore_neutron mccpsd_monitor_beforelog_restore_neutron
#line 62 "/usr/share/mcstas/2.7/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
  DArray2d PSD_N;
  DArray2d PSD_p;
  DArray2d PSD_p2;
#line 9869 "./reverse_test.c"
#undef restore_neutron
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef filename
#undef ny
#undef nx
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'logspir' [9]. */
#define mccompcurname  logspir
#define mccompcurtype  LogSpiral
#define mccompcurindex 9
#define zmin mcclogspir_zmin
#define zmax mcclogspir_zmax
#define ymin mcclogspir_ymin
#define ymax mcclogspir_ymax
#define psi mcclogspir_psi
#define phi_rot mcclogspir_phi_rot
#define precision mcclogspir_precision
#define max_iterations mcclogspir_max_iterations
#define mValue mcclogspir_mValue
#define branches mcclogspir_branches
#define doublesided mcclogspir_doublesided
#define placeholder mcclogspir_placeholder
#line 593 "LogSpiral.comp"
	SceneLog s;
#line 9905 "./reverse_test.c"
#undef placeholder
#undef doublesided
#undef branches
#undef mValue
#undef max_iterations
#undef precision
#undef phi_rot
#undef psi
#undef ymax
#undef ymin
#undef zmax
#undef zmin
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

/* User declarations for component 'psd_monitor' [10]. */
#define mccompcurname  psd_monitor
#define mccompcurtype  PSD_monitor
#define mccompcurindex 10
#define PSD_N mccpsd_monitor_PSD_N
#define PSD_p mccpsd_monitor_PSD_p
#define PSD_p2 mccpsd_monitor_PSD_p2
#define nx mccpsd_monitor_nx
#define ny mccpsd_monitor_ny
#define filename mccpsd_monitor_filename
#define xmin mccpsd_monitor_xmin
#define xmax mccpsd_monitor_xmax
#define ymin mccpsd_monitor_ymin
#define ymax mccpsd_monitor_ymax
#define xwidth mccpsd_monitor_xwidth
#define yheight mccpsd_monitor_yheight
#define restore_neutron mccpsd_monitor_restore_neutron
#line 62 "/usr/share/mcstas/2.7/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
  DArray2d PSD_N;
  DArray2d PSD_p;
  DArray2d PSD_p2;
#line 9943 "./reverse_test.c"
#undef restore_neutron
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef filename
#undef ny
#undef nx
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

Coords mcposaorigin, mcposrorigin;
Rotation mcrotaorigin, mcrotrorigin;
Coords mcposasource, mcposrsource;
Rotation mcrotasource, mcrotrsource;
Coords mcposasource_div, mcposrsource_div;
Rotation mcrotasource_div, mcrotrsource_div;
Coords mcposaslit, mcposrslit;
Rotation mcrotaslit, mcrotrslit;
Coords mcposapsd_before_optic, mcposrpsd_before_optic;
Rotation mcrotapsd_before_optic, mcrotrpsd_before_optic;
Coords mcposaflat_ellipse_horizontal, mcposrflat_ellipse_horizontal;
Rotation mcrotaflat_ellipse_horizontal, mcrotrflat_ellipse_horizontal;
Coords mcposaslit1, mcposrslit1;
Rotation mcrotaslit1, mcrotrslit1;
Coords mcposapsd_monitor_beforelog, mcposrpsd_monitor_beforelog;
Rotation mcrotapsd_monitor_beforelog, mcrotrpsd_monitor_beforelog;
Coords mcposalogspir, mcposrlogspir;
Rotation mcrotalogspir, mcrotrlogspir;
Coords mcposapsd_monitor, mcposrpsd_monitor;
Rotation mcrotapsd_monitor, mcrotrpsd_monitor;

MCNUM mcnx, mcny, mcnz, mcnvx, mcnvy, mcnvz, mcnt, mcnsx, mcnsy, mcnsz, mcnp;

/* end declare */

void mcinit(void) {
#define mccompcurname  logspir_test
#define mccompcurtype  INSTRUMENT
#define mccompcurindex 0
#define mcposalogspir_test coords_set(0,0,0)
#define source_width mcipsource_width
#define source_divergence mcipsource_divergence
#define L_source mcipL_source
#define dL mcipdL
#define flux mcipflux
#define zmin mcipzmin
#define zmax mcipzmax
#define ymin mcipymin
#define ymax mcipymax
#define psi mcippsi
#define phi_rot mcipphi_rot
#define precision mcipprecision
#define max_iterations mcipmax_iterations
#define mValue mcipmValue
#define branches mcipbranches
#define doublesided mcipdoublesided
#define placeholder mcipplaceholder
#undef placeholder
#undef doublesided
#undef branches
#undef mValue
#undef max_iterations
#undef precision
#undef phi_rot
#undef psi
#undef ymax
#undef ymin
#undef zmax
#undef zmin
#undef flux
#undef dL
#undef L_source
#undef source_divergence
#undef source_width
#undef mcposalogspir_test
#undef mccompcurindex
#undef mccompcurtype
#undef mccompcurname
  /* Computation of coordinate transformations. */
  {
    Coords mctc1, mctc2, mcLastComp;
    Rotation mctr1;
    double mcAccumulatedILength = 0;
    /* Initialize "last" component origin as (0,0,0) */
    mcLastComp = coords_set(0,0,0);

    mcDEBUG_INSTR()
  /* Component initializations. */
    /* Component origin. */
  /* Setting parameters for component origin. */
  SIG_MESSAGE("origin (Init:SetPar)");
#line 39 "reverse_test.instr"
  if("NULL") strncpy(mccorigin_profile, "NULL" ? "NULL" : "", 16384); else mccorigin_profile[0]='\0';
#line 39 "reverse_test.instr"
  mccorigin_percent = 10;
#line 39 "reverse_test.instr"
  mccorigin_flag_save = 0;
#line 39 "reverse_test.instr"
  mccorigin_minutes = 0;
#line 10050 "./reverse_test.c"

  SIG_MESSAGE("origin (Init:Place/Rotate)");
  rot_set_rotation(mcrotaorigin,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 10057 "./reverse_test.c"
  rot_copy(mcrotrorigin, mcrotaorigin);
  mcposaorigin = coords_set(
#line 57 "reverse_test.instr"
    0,
#line 57 "reverse_test.instr"
    0,
#line 57 "reverse_test.instr"
    0);
#line 10066 "./reverse_test.c"
  mctc1 = coords_neg(mcposaorigin);
  mcposrorigin = rot_apply(mcrotaorigin, mctc1);
  mcDEBUG_COMPONENT("origin", mcposaorigin, mcrotaorigin)
  mccomp_posa[1] = mcposaorigin;
  mccomp_posr[1] = mcposrorigin;
  mcNCounter[1]  = mcPCounter[1] = mcP2Counter[1] = 0;
  mcAbsorbProp[1]= 0;
    /* Component source. */
  /* Setting parameters for component source. */
  SIG_MESSAGE("source (Init:SetPar)");

  SIG_MESSAGE("source (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 10083 "./reverse_test.c"
  rot_mul(mctr1, mcrotaorigin, mcrotasource);
  rot_transpose(mcrotaorigin, mctr1);
  rot_mul(mcrotasource, mctr1, mcrotrsource);
  mctc1 = coords_set(
#line 61 "reverse_test.instr"
    0,
#line 61 "reverse_test.instr"
    0,
#line 61 "reverse_test.instr"
    0);
#line 10094 "./reverse_test.c"
  rot_transpose(mcrotaorigin, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposasource = coords_add(mcposaorigin, mctc2);
  mctc1 = coords_sub(mcposaorigin, mcposasource);
  mcposrsource = rot_apply(mcrotasource, mctc1);
  mcDEBUG_COMPONENT("source", mcposasource, mcrotasource)
  mccomp_posa[2] = mcposasource;
  mccomp_posr[2] = mcposrsource;
  mcNCounter[2]  = mcPCounter[2] = mcP2Counter[2] = 0;
  mcAbsorbProp[2]= 0;
    /* Component source_div. */
  /* Setting parameters for component source_div. */
  SIG_MESSAGE("source_div (Init:SetPar)");
#line 65 "reverse_test.instr"
  mccsource_div_xwidth = mcipsource_width;
#line 64 "reverse_test.instr"
  mccsource_div_yheight = 0.024;
#line 66 "reverse_test.instr"
  mccsource_div_focus_aw = mcipsource_divergence;
#line 67 "reverse_test.instr"
  mccsource_div_focus_ah = 0.000001;
#line 64 "reverse_test.instr"
  mccsource_div_E0 = 0.0;
#line 64 "reverse_test.instr"
  mccsource_div_dE = 0.0;
#line 68 "reverse_test.instr"
  mccsource_div_lambda0 = mcipL_source;
#line 70 "reverse_test.instr"
  mccsource_div_dlambda = 0;
#line 64 "reverse_test.instr"
  mccsource_div_gauss = 0;
#line 69 "reverse_test.instr"
  mccsource_div_flux = 1e11;
#line 10128 "./reverse_test.c"

  SIG_MESSAGE("source_div (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 72 "reverse_test.instr"
    (0)*DEG2RAD,
#line 72 "reverse_test.instr"
    (0)*DEG2RAD,
#line 72 "reverse_test.instr"
    (0)*DEG2RAD);
#line 10138 "./reverse_test.c"
  rot_mul(mctr1, mcrotasource, mcrotasource_div);
  rot_transpose(mcrotasource, mctr1);
  rot_mul(mcrotasource_div, mctr1, mcrotrsource_div);
  mctc1 = coords_set(
#line 71 "reverse_test.instr"
    0,
#line 71 "reverse_test.instr"
    0,
#line 71 "reverse_test.instr"
    0);
#line 10149 "./reverse_test.c"
  rot_transpose(mcrotasource, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposasource_div = coords_add(mcposasource, mctc2);
  mctc1 = coords_sub(mcposasource, mcposasource_div);
  mcposrsource_div = rot_apply(mcrotasource_div, mctc1);
  mcDEBUG_COMPONENT("source_div", mcposasource_div, mcrotasource_div)
  mccomp_posa[3] = mcposasource_div;
  mccomp_posr[3] = mcposrsource_div;
  mcNCounter[3]  = mcPCounter[3] = mcP2Counter[3] = 0;
  mcAbsorbProp[3]= 0;
    /* Component slit. */
  /* Setting parameters for component slit. */
  SIG_MESSAGE("slit (Init:SetPar)");
#line 75 "reverse_test.instr"
  mccslit_xmin = -0.0198;
#line 76 "reverse_test.instr"
  mccslit_xmax = 0.0198;
#line 77 "reverse_test.instr"
  mccslit_ymin = -1;
#line 78 "reverse_test.instr"
  mccslit_ymax = 1;
#line 46 "reverse_test.instr"
  mccslit_radius = 0;
#line 46 "reverse_test.instr"
  mccslit_xwidth = 0;
#line 46 "reverse_test.instr"
  mccslit_yheight = 0;
#line 10177 "./reverse_test.c"

  SIG_MESSAGE("slit (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 82 "reverse_test.instr"
    (0)*DEG2RAD,
#line 82 "reverse_test.instr"
    (0)*DEG2RAD,
#line 82 "reverse_test.instr"
    (0)*DEG2RAD);
#line 10187 "./reverse_test.c"
  rot_mul(mctr1, mcrotasource, mcrotaslit);
  rot_transpose(mcrotasource_div, mctr1);
  rot_mul(mcrotaslit, mctr1, mcrotrslit);
  mctc1 = coords_set(
#line 81 "reverse_test.instr"
    0,
#line 81 "reverse_test.instr"
    0,
#line 81 "reverse_test.instr"
    0.675 -0.060001);
#line 10198 "./reverse_test.c"
  rot_transpose(mcrotasource, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaslit = coords_add(mcposasource, mctc2);
  mctc1 = coords_sub(mcposasource_div, mcposaslit);
  mcposrslit = rot_apply(mcrotaslit, mctc1);
  mcDEBUG_COMPONENT("slit", mcposaslit, mcrotaslit)
  mccomp_posa[4] = mcposaslit;
  mccomp_posr[4] = mcposrslit;
  mcNCounter[4]  = mcPCounter[4] = mcP2Counter[4] = 0;
  mcAbsorbProp[4]= 0;
    /* Component psd_before_optic. */
  /* Setting parameters for component psd_before_optic. */
  SIG_MESSAGE("psd_before_optic (Init:SetPar)");
#line 85 "reverse_test.instr"
  mccpsd_before_optic_nx = 500;
#line 86 "reverse_test.instr"
  mccpsd_before_optic_ny = 500;
#line 87 "reverse_test.instr"
  if("beforeoptic.dat") strncpy(mccpsd_before_optic_filename, "beforeoptic.dat" ? "beforeoptic.dat" : "", 16384); else mccpsd_before_optic_filename[0]='\0';
#line 50 "reverse_test.instr"
  mccpsd_before_optic_xmin = -0.05;
#line 50 "reverse_test.instr"
  mccpsd_before_optic_xmax = 0.05;
#line 50 "reverse_test.instr"
  mccpsd_before_optic_ymin = -0.05;
#line 50 "reverse_test.instr"
  mccpsd_before_optic_ymax = 0.05;
#line 88 "reverse_test.instr"
  mccpsd_before_optic_xwidth = 1;
#line 89 "reverse_test.instr"
  mccpsd_before_optic_yheight = 1;
#line 90 "reverse_test.instr"
  mccpsd_before_optic_restore_neutron = 1;
#line 10232 "./reverse_test.c"

  SIG_MESSAGE("psd_before_optic (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD,
    (0.0)*DEG2RAD);
#line 10239 "./reverse_test.c"
  rot_mul(mctr1, mcrotasource, mcrotapsd_before_optic);
  rot_transpose(mcrotaslit, mctr1);
  rot_mul(mcrotapsd_before_optic, mctr1, mcrotrpsd_before_optic);
  mctc1 = coords_set(
#line 91 "reverse_test.instr"
    0,
#line 91 "reverse_test.instr"
    0,
#line 91 "reverse_test.instr"
    0.675 -0.06);
#line 10250 "./reverse_test.c"
  rot_transpose(mcrotasource, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposapsd_before_optic = coords_add(mcposasource, mctc2);
  mctc1 = coords_sub(mcposaslit, mcposapsd_before_optic);
  mcposrpsd_before_optic = rot_apply(mcrotapsd_before_optic, mctc1);
  mcDEBUG_COMPONENT("psd_before_optic", mcposapsd_before_optic, mcrotapsd_before_optic)
  mccomp_posa[5] = mcposapsd_before_optic;
  mccomp_posr[5] = mcposrpsd_before_optic;
  mcNCounter[5]  = mcPCounter[5] = mcP2Counter[5] = 0;
  mcAbsorbProp[5]= 0;
    /* Component flat_ellipse_horizontal. */
  /* Setting parameters for component flat_ellipse_horizontal. */
  SIG_MESSAGE("flat_ellipse_horizontal (Init:SetPar)");
#line 94 "reverse_test.instr"
  mccflat_ellipse_horizontal_sourceDist = - ( 0.675 );
#line 95 "reverse_test.instr"
  mccflat_ellipse_horizontal_LStart = - ( 0.675 );
#line 96 "reverse_test.instr"
  mccflat_ellipse_horizontal_LEnd = 0.675;
#line 97 "reverse_test.instr"
  mccflat_ellipse_horizontal_lStart = -0.06;
#line 98 "reverse_test.instr"
  mccflat_ellipse_horizontal_lEnd = 0.06;
#line 99 "reverse_test.instr"
  mccflat_ellipse_horizontal_r_0 = 0.02;
#line 102 "reverse_test.instr"
  mccflat_ellipse_horizontal_nummirror = 30;
#line 104 "reverse_test.instr"
  mccflat_ellipse_horizontal_mf = 4;
#line 105 "reverse_test.instr"
  mccflat_ellipse_horizontal_mb = 0;
#line 100 "reverse_test.instr"
  mccflat_ellipse_horizontal_mirror_width = 0;
#line 101 "reverse_test.instr"
  mccflat_ellipse_horizontal_mirror_sidelength = 10;
#line 103 "reverse_test.instr"
  mccflat_ellipse_horizontal_doubleReflections = 1;
#line 10288 "./reverse_test.c"

  SIG_MESSAGE("flat_ellipse_horizontal (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 108 "reverse_test.instr"
    (0)*DEG2RAD,
#line 108 "reverse_test.instr"
    (0)*DEG2RAD,
#line 108 "reverse_test.instr"
    (0)*DEG2RAD);
#line 10298 "./reverse_test.c"
  rot_mul(mctr1, mcrotasource, mcrotaflat_ellipse_horizontal);
  rot_transpose(mcrotapsd_before_optic, mctr1);
  rot_mul(mcrotaflat_ellipse_horizontal, mctr1, mcrotrflat_ellipse_horizontal);
  mctc1 = coords_set(
#line 107 "reverse_test.instr"
    0,
#line 107 "reverse_test.instr"
    0,
#line 107 "reverse_test.instr"
    0.675);
#line 10309 "./reverse_test.c"
  rot_transpose(mcrotasource, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaflat_ellipse_horizontal = coords_add(mcposasource, mctc2);
  mctc1 = coords_sub(mcposapsd_before_optic, mcposaflat_ellipse_horizontal);
  mcposrflat_ellipse_horizontal = rot_apply(mcrotaflat_ellipse_horizontal, mctc1);
  mcDEBUG_COMPONENT("flat_ellipse_horizontal", mcposaflat_ellipse_horizontal, mcrotaflat_ellipse_horizontal)
  mccomp_posa[6] = mcposaflat_ellipse_horizontal;
  mccomp_posr[6] = mcposrflat_ellipse_horizontal;
  mcNCounter[6]  = mcPCounter[6] = mcP2Counter[6] = 0;
  mcAbsorbProp[6]= 0;
    /* Component slit1. */
  /* Setting parameters for component slit1. */
  SIG_MESSAGE("slit1 (Init:SetPar)");
#line 111 "reverse_test.instr"
  mccslit1_xmin = -0.0198;
#line 112 "reverse_test.instr"
  mccslit1_xmax = 0.0198;
#line 113 "reverse_test.instr"
  mccslit1_ymin = -1;
#line 114 "reverse_test.instr"
  mccslit1_ymax = 1;
#line 46 "reverse_test.instr"
  mccslit1_radius = 0;
#line 46 "reverse_test.instr"
  mccslit1_xwidth = 0;
#line 46 "reverse_test.instr"
  mccslit1_yheight = 0;
#line 10337 "./reverse_test.c"

  SIG_MESSAGE("slit1 (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 118 "reverse_test.instr"
    (0)*DEG2RAD,
#line 118 "reverse_test.instr"
    (0)*DEG2RAD,
#line 118 "reverse_test.instr"
    (0)*DEG2RAD);
#line 10347 "./reverse_test.c"
  rot_mul(mctr1, mcrotasource, mcrotaslit1);
  rot_transpose(mcrotaflat_ellipse_horizontal, mctr1);
  rot_mul(mcrotaslit1, mctr1, mcrotrslit1);
  mctc1 = coords_set(
#line 117 "reverse_test.instr"
    0,
#line 117 "reverse_test.instr"
    0,
#line 117 "reverse_test.instr"
    0.675 + 0.060001);
#line 10358 "./reverse_test.c"
  rot_transpose(mcrotasource, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposaslit1 = coords_add(mcposasource, mctc2);
  mctc1 = coords_sub(mcposaflat_ellipse_horizontal, mcposaslit1);
  mcposrslit1 = rot_apply(mcrotaslit1, mctc1);
  mcDEBUG_COMPONENT("slit1", mcposaslit1, mcrotaslit1)
  mccomp_posa[7] = mcposaslit1;
  mccomp_posr[7] = mcposrslit1;
  mcNCounter[7]  = mcPCounter[7] = mcP2Counter[7] = 0;
  mcAbsorbProp[7]= 0;
    /* Component psd_monitor_beforelog. */
  /* Setting parameters for component psd_monitor_beforelog. */
  SIG_MESSAGE("psd_monitor_beforelog (Init:SetPar)");
#line 121 "reverse_test.instr"
  mccpsd_monitor_beforelog_nx = 1001;
#line 122 "reverse_test.instr"
  mccpsd_monitor_beforelog_ny = 1001;
#line 123 "reverse_test.instr"
  if("psdbeforelog.dat") strncpy(mccpsd_monitor_beforelog_filename, "psdbeforelog.dat" ? "psdbeforelog.dat" : "", 16384); else mccpsd_monitor_beforelog_filename[0]='\0';
#line 50 "reverse_test.instr"
  mccpsd_monitor_beforelog_xmin = -0.05;
#line 50 "reverse_test.instr"
  mccpsd_monitor_beforelog_xmax = 0.05;
#line 50 "reverse_test.instr"
  mccpsd_monitor_beforelog_ymin = -0.05;
#line 50 "reverse_test.instr"
  mccpsd_monitor_beforelog_ymax = 0.05;
#line 124 "reverse_test.instr"
  mccpsd_monitor_beforelog_xwidth = 0.5;
#line 125 "reverse_test.instr"
  mccpsd_monitor_beforelog_yheight = 0.5;
#line 126 "reverse_test.instr"
  mccpsd_monitor_beforelog_restore_neutron = 1;
#line 10392 "./reverse_test.c"

  SIG_MESSAGE("psd_monitor_beforelog (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 128 "reverse_test.instr"
    (0)*DEG2RAD,
#line 128 "reverse_test.instr"
    (0)*DEG2RAD,
#line 128 "reverse_test.instr"
    (0)*DEG2RAD);
#line 10402 "./reverse_test.c"
  rot_mul(mctr1, mcrotasource, mcrotapsd_monitor_beforelog);
  rot_transpose(mcrotaslit1, mctr1);
  rot_mul(mcrotapsd_monitor_beforelog, mctr1, mcrotrpsd_monitor_beforelog);
  mctc1 = coords_set(
#line 127 "reverse_test.instr"
    0,
#line 127 "reverse_test.instr"
    0,
#line 127 "reverse_test.instr"
    0.675 + 0.0601);
#line 10413 "./reverse_test.c"
  rot_transpose(mcrotasource, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposapsd_monitor_beforelog = coords_add(mcposasource, mctc2);
  mctc1 = coords_sub(mcposaslit1, mcposapsd_monitor_beforelog);
  mcposrpsd_monitor_beforelog = rot_apply(mcrotapsd_monitor_beforelog, mctc1);
  mcDEBUG_COMPONENT("psd_monitor_beforelog", mcposapsd_monitor_beforelog, mcrotapsd_monitor_beforelog)
  mccomp_posa[8] = mcposapsd_monitor_beforelog;
  mccomp_posr[8] = mcposrpsd_monitor_beforelog;
  mcNCounter[8]  = mcPCounter[8] = mcP2Counter[8] = 0;
  mcAbsorbProp[8]= 0;
    /* Component logspir. */
  /* Setting parameters for component logspir. */
  SIG_MESSAGE("logspir (Init:SetPar)");
#line 132 "reverse_test.instr"
  mcclogspir_zmin = 0.15;
#line 133 "reverse_test.instr"
  mcclogspir_zmax = 0.3;
#line 134 "reverse_test.instr"
  mcclogspir_ymin = -0.05;
#line 135 "reverse_test.instr"
  mcclogspir_ymax = 0.05;
#line 136 "reverse_test.instr"
  mcclogspir_psi = mcippsi;
#line 138 "reverse_test.instr"
  mcclogspir_phi_rot = mcipphi_rot;
#line 44 "reverse_test.instr"
  mcclogspir_precision = 1e-7;
#line 45 "reverse_test.instr"
  mcclogspir_max_iterations = 10;
#line 139 "reverse_test.instr"
  mcclogspir_mValue = 6.2;
#line 137 "reverse_test.instr"
  mcclogspir_branches = mcipbranches;
#line 140 "reverse_test.instr"
  mcclogspir_doublesided = mcipdoublesided;
#line 49 "reverse_test.instr"
  mcclogspir_placeholder = 0;
#line 10451 "./reverse_test.c"

  SIG_MESSAGE("logspir (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 142 "reverse_test.instr"
    (0)*DEG2RAD,
#line 142 "reverse_test.instr"
    (180)*DEG2RAD,
#line 142 "reverse_test.instr"
    (0)*DEG2RAD);
#line 10461 "./reverse_test.c"
  rot_mul(mctr1, mcrotasource, mcrotalogspir);
  rot_transpose(mcrotapsd_monitor_beforelog, mctr1);
  rot_mul(mcrotalogspir, mctr1, mcrotrlogspir);
  mctc1 = coords_set(
#line 141 "reverse_test.instr"
    0,
#line 141 "reverse_test.instr"
    0,
#line 141 "reverse_test.instr"
    0.675 * 2);
#line 10472 "./reverse_test.c"
  rot_transpose(mcrotasource, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposalogspir = coords_add(mcposasource, mctc2);
  mctc1 = coords_sub(mcposapsd_monitor_beforelog, mcposalogspir);
  mcposrlogspir = rot_apply(mcrotalogspir, mctc1);
  mcDEBUG_COMPONENT("logspir", mcposalogspir, mcrotalogspir)
  mccomp_posa[9] = mcposalogspir;
  mccomp_posr[9] = mcposrlogspir;
  mcNCounter[9]  = mcPCounter[9] = mcP2Counter[9] = 0;
  mcAbsorbProp[9]= 0;
    /* Component psd_monitor. */
  /* Setting parameters for component psd_monitor. */
  SIG_MESSAGE("psd_monitor (Init:SetPar)");
#line 148 "reverse_test.instr"
  mccpsd_monitor_nx = 1001;
#line 149 "reverse_test.instr"
  mccpsd_monitor_ny = 1001;
#line 150 "reverse_test.instr"
  if("psdafterlog.dat") strncpy(mccpsd_monitor_filename, "psdafterlog.dat" ? "psdafterlog.dat" : "", 16384); else mccpsd_monitor_filename[0]='\0';
#line 50 "reverse_test.instr"
  mccpsd_monitor_xmin = -0.05;
#line 50 "reverse_test.instr"
  mccpsd_monitor_xmax = 0.05;
#line 50 "reverse_test.instr"
  mccpsd_monitor_ymin = -0.05;
#line 50 "reverse_test.instr"
  mccpsd_monitor_ymax = 0.05;
#line 151 "reverse_test.instr"
  mccpsd_monitor_xwidth = 0.026;
#line 152 "reverse_test.instr"
  mccpsd_monitor_yheight = 0.026;
#line 153 "reverse_test.instr"
  mccpsd_monitor_restore_neutron = 1;
#line 10506 "./reverse_test.c"

  SIG_MESSAGE("psd_monitor (Init:Place/Rotate)");
  rot_set_rotation(mctr1,
#line 155 "reverse_test.instr"
    (0)*DEG2RAD,
#line 155 "reverse_test.instr"
    (0)*DEG2RAD,
#line 155 "reverse_test.instr"
    (0)*DEG2RAD);
#line 10516 "./reverse_test.c"
  rot_mul(mctr1, mcrotasource, mcrotapsd_monitor);
  rot_transpose(mcrotalogspir, mctr1);
  rot_mul(mcrotapsd_monitor, mctr1, mcrotrpsd_monitor);
  mctc1 = coords_set(
#line 154 "reverse_test.instr"
    0,
#line 154 "reverse_test.instr"
    0,
#line 154 "reverse_test.instr"
    2 * 0.675);
#line 10527 "./reverse_test.c"
  rot_transpose(mcrotasource, mctr1);
  mctc2 = rot_apply(mctr1, mctc1);
  mcposapsd_monitor = coords_add(mcposasource, mctc2);
  mctc1 = coords_sub(mcposalogspir, mcposapsd_monitor);
  mcposrpsd_monitor = rot_apply(mcrotapsd_monitor, mctc1);
  mcDEBUG_COMPONENT("psd_monitor", mcposapsd_monitor, mcrotapsd_monitor)
  mccomp_posa[10] = mcposapsd_monitor;
  mccomp_posr[10] = mcposrpsd_monitor;
  mcNCounter[10]  = mcPCounter[10] = mcP2Counter[10] = 0;
  mcAbsorbProp[10]= 0;
  /* Component initializations. */
  /* Initializations for component origin. */
  SIG_MESSAGE("origin (Init)");
#define mccompcurname  origin
#define mccompcurtype  Progress_bar
#define mccompcurindex 1
#define IntermediateCnts mccorigin_IntermediateCnts
#define StartTime mccorigin_StartTime
#define EndTime mccorigin_EndTime
#define CurrentTime mccorigin_CurrentTime
#define profile mccorigin_profile
#define percent mccorigin_percent
#define flag_save mccorigin_flag_save
#define minutes mccorigin_minutes
#line 57 "/usr/share/mcstas/2.7/tools/Python/mcrun/../mccodelib/../../../misc/Progress_bar.comp"
{
IntermediateCnts=0;
StartTime=0;
EndTime=0;
CurrentTime=0;

fprintf(stdout, "[%s] Initialize\n", mcinstrument_name);
  if (percent*mcget_ncount()/100 < 1e5) {
    percent=1e5*100.0/mcget_ncount();
  }
}
#line 10564 "./reverse_test.c"
#undef minutes
#undef flag_save
#undef percent
#undef profile
#undef CurrentTime
#undef EndTime
#undef StartTime
#undef IntermediateCnts
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component source. */
  SIG_MESSAGE("source (Init)");

  /* Initializations for component source_div. */
  SIG_MESSAGE("source_div (Init)");
#define mccompcurname  source_div
#define mccompcurtype  Source_div
#define mccompcurindex 3
#define thetah mccsource_div_thetah
#define thetav mccsource_div_thetav
#define sigmah mccsource_div_sigmah
#define sigmav mccsource_div_sigmav
#define tan_h mccsource_div_tan_h
#define tan_v mccsource_div_tan_v
#define p_init mccsource_div_p_init
#define dist mccsource_div_dist
#define focus_xw mccsource_div_focus_xw
#define focus_yh mccsource_div_focus_yh
#define xwidth mccsource_div_xwidth
#define yheight mccsource_div_yheight
#define focus_aw mccsource_div_focus_aw
#define focus_ah mccsource_div_focus_ah
#define E0 mccsource_div_E0
#define dE mccsource_div_dE
#define lambda0 mccsource_div_lambda0
#define dlambda mccsource_div_dlambda
#define gauss mccsource_div_gauss
#define flux mccsource_div_flux
#line 72 "/usr/share/mcstas/2.7/tools/Python/mcrun/../mccodelib/../../../sources/Source_div.comp"
{
sigmah = DEG2RAD*focus_aw/(sqrt(8.0*log(2.0)));
  sigmav = DEG2RAD*focus_ah/(sqrt(8.0*log(2.0)));

  if (xwidth < 0 || yheight < 0 || focus_aw < 0 || focus_ah < 0) {
      printf("Source_div: %s: Error in input parameter values!\n"
             "ERROR       Exiting\n",
           NAME_CURRENT_COMP);
      exit(-1);
  }
  if ((!lambda0 && !E0 && !dE && !dlambda)) {
    printf("Source_div: %s: You must specify either a wavelength or energy range!\n ERROR - Exiting\n",
           NAME_CURRENT_COMP);
    exit(-1);
  }
  if ((!lambda0 && !dlambda && (E0 <= 0 || dE < 0 || E0-dE <= 0))
    || (!E0 && !dE && (lambda0 <= 0 || dlambda < 0 || lambda0-dlambda <= 0))) {
    printf("Source_div: %s: Unmeaningful definition of wavelength or energy range!\n ERROR - Exiting\n",
           NAME_CURRENT_COMP);
      exit(-1);
  }
  /* compute distance to next component */
  Coords ToTarget;
  double tx,ty,tz;
  ToTarget = coords_sub(POS_A_COMP_INDEX(INDEX_CURRENT_COMP+1),POS_A_CURRENT_COMP);
  ToTarget = rot_apply(ROT_A_CURRENT_COMP, ToTarget);
  coords_get(ToTarget, &tx, &ty, &tz);
  dist=sqrt(tx*tx+ty*ty+tz*tz);
  /* compute target area */
  if (dist) {
    focus_xw=dist*tan(focus_aw*DEG2RAD);
    focus_yh=dist*tan(focus_ah*DEG2RAD);
  }

  p_init  = flux*1e4*xwidth*yheight/mcget_ncount();
  if (!focus_aw || !focus_ah)
    exit(printf("Source_div: %s: Zero divergence defined. \n"
                "ERROR       Use non zero values for focus_aw and focus_ah.\n",
           NAME_CURRENT_COMP));
  p_init *= 2*fabs(DEG2RAD*focus_aw*sin(DEG2RAD*focus_ah/2));  /* solid angle */
  if (dlambda)
    p_init *= 2*dlambda;
  else if (dE)
    p_init *= 2*dE;
}
#line 10651 "./reverse_test.c"
#undef flux
#undef gauss
#undef dlambda
#undef lambda0
#undef dE
#undef E0
#undef focus_ah
#undef focus_aw
#undef yheight
#undef xwidth
#undef focus_yh
#undef focus_xw
#undef dist
#undef p_init
#undef tan_v
#undef tan_h
#undef sigmav
#undef sigmah
#undef thetav
#undef thetah
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component slit. */
  SIG_MESSAGE("slit (Init)");
#define mccompcurname  slit
#define mccompcurtype  Slit
#define mccompcurindex 4
#define xmin mccslit_xmin
#define xmax mccslit_xmax
#define ymin mccslit_ymin
#define ymax mccslit_ymax
#define radius mccslit_radius
#define xwidth mccslit_xwidth
#define yheight mccslit_yheight
#line 50 "/usr/share/mcstas/2.7/tools/Python/mcrun/../mccodelib/../../../optics/Slit.comp"
{
if (xwidth > 0)  { 
  if (!xmin && !xmax) {
    xmax=xwidth/2;  xmin=-xmax;
  } else {
    fprintf(stderr,"Slit: %s: Error: please specify EITHER xmin & xmax or xwidth\n", NAME_CURRENT_COMP); exit(-1);
  }
 }
 if (yheight > 0) { 
   if (!ymin && !ymax) {
     ymax=yheight/2; ymin=-ymax; 
   } else {
     fprintf(stderr,"Slit: %s: Error: please specify EITHER ymin & ymax or ywidth\n", NAME_CURRENT_COMP); exit(-1);
   }
 }
 if (xmin == 0 && xmax == 0 && ymin == 0 && ymax == 0 && radius == 0)
    { fprintf(stderr,"Slit: %s: Warning: Running with CLOSED slit - is this intentional?? \n", NAME_CURRENT_COMP); }

}
#line 10708 "./reverse_test.c"
#undef yheight
#undef xwidth
#undef radius
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component psd_before_optic. */
  SIG_MESSAGE("psd_before_optic (Init)");
#define mccompcurname  psd_before_optic
#define mccompcurtype  PSD_monitor
#define mccompcurindex 5
#define PSD_N mccpsd_before_optic_PSD_N
#define PSD_p mccpsd_before_optic_PSD_p
#define PSD_p2 mccpsd_before_optic_PSD_p2
#define nx mccpsd_before_optic_nx
#define ny mccpsd_before_optic_ny
#define filename mccpsd_before_optic_filename
#define xmin mccpsd_before_optic_xmin
#define xmax mccpsd_before_optic_xmax
#define ymin mccpsd_before_optic_ymin
#define ymax mccpsd_before_optic_ymax
#define xwidth mccpsd_before_optic_xwidth
#define yheight mccpsd_before_optic_yheight
#define restore_neutron mccpsd_before_optic_restore_neutron
#line 68 "/usr/share/mcstas/2.7/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
  if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }

  if ((xmin >= xmax) || (ymin >= ymax)){
    printf("PSD_monitor: %s: Null detection area !\n"
           "ERROR        (xwidth,yheight,xmin,xmax,ymin,ymax). Exiting",
    NAME_CURRENT_COMP);
    exit(-1);
  }

  PSD_N = create_darr2d(nx, ny);
  PSD_p = create_darr2d(nx, ny);
  PSD_p2 = create_darr2d(nx, ny);

  int i, j;
  for (i=0; i<nx; i++){
    for (j=0; j<ny; j++){
      PSD_N[i][j] = 0;
      PSD_p[i][j] = 0;
      PSD_p2[i][j] = 0;
    }
  }
}
#line 10763 "./reverse_test.c"
#undef restore_neutron
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef filename
#undef ny
#undef nx
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component flat_ellipse_horizontal. */
  SIG_MESSAGE("flat_ellipse_horizontal (Init)");
#define mccompcurname  flat_ellipse_horizontal
#define mccompcurtype  FlatEllipse_finite_mirror
#define mccompcurindex 6
#define reflect mccflat_ellipse_horizontal_reflect
#define s mccflat_ellipse_horizontal_s
#define pTable mccflat_ellipse_horizontal_pTable
#define R0 mccflat_ellipse_horizontal_R0
#define Qc mccflat_ellipse_horizontal_Qc
#define W mccflat_ellipse_horizontal_W
#define alpha mccflat_ellipse_horizontal_alpha
#define transmit mccflat_ellipse_horizontal_transmit
#define sourceDist mccflat_ellipse_horizontal_sourceDist
#define LStart mccflat_ellipse_horizontal_LStart
#define LEnd mccflat_ellipse_horizontal_LEnd
#define lStart mccflat_ellipse_horizontal_lStart
#define lEnd mccflat_ellipse_horizontal_lEnd
#define r_0 mccflat_ellipse_horizontal_r_0
#define nummirror mccflat_ellipse_horizontal_nummirror
#define mf mccflat_ellipse_horizontal_mf
#define mb mccflat_ellipse_horizontal_mb
#define mirror_width mccflat_ellipse_horizontal_mirror_width
#define mirror_sidelength mccflat_ellipse_horizontal_mirror_sidelength
#define doubleReflections mccflat_ellipse_horizontal_doubleReflections
#line 81 "FlatEllipse_finite_mirror.comp"
{
    if (sourceDist == 0){
        sourceDist = LStart;
    }
    pointer_lStart = &lStart;
    //Load Reflectivity Data File
    /*if (reflect && strlen(reflect)) {
        if (Table_Read(&pTable, reflect, 1) <= 0)
            exit(fprintf(stderr, "Can not read file: %s\n", reflect));
    }

    //Custom function for tracing neutrons using table data for reflectivity
    void traceNeutronConicWithTables(Particle* pa, ConicSurf c) {
        double tl = getTimeOfFirstCollisionConic(*pa, c);
        if (tl < 0)
            return;
        else {
            //Move Particle to Surface Edge
            moveParticleT(tl,pa);

            if (c.m==0) {
                absorbParticle(pa);
                return;
            }

            //Handle Reflectivity
            Vec n = getNormConic(getParticlePos(*pa),c);
            double vdotn = dotVec(n,getParticleVel(*pa));

            double q = fabs(2*vdotn*V2Q);
            double B;
            if (reflect && strlen(reflect))
                B=Table_Value(pTable, q, 1);
            else {
                B = R0;
                if (q > Qc) {
                    double arg = (q-c.m*Qc)/W;
                    if(arg < 10)
                        B *= .5*(1-tanh(arg))*(1-alpha*(q-Qc));
                    else
                        B=0;
                }
            }
            if (B < 0)
                B=0;
            else if (B > 1)
                B=1;
            if (!transmit) {
                if (!B) absorbParticle(pa);
                pa->w *= B;
                reflectParticle(n,pa);
            } else {
                if (B == 0 || rand01() >= B) { /*unreflected*/ /*}
                else { reflectParticle(n,pa); }
            }
        }
    }
    */
    //Make new scene
    silicon = (mirror_width==0) ? 0 : -1; //neutron starts in air by default
    s = makeScene();
    rfront_inner = get_r_at_z0(nummirror, 0, r_0, lStart, sourceDist, LEnd, lStart, lEnd);
    
    //Set Scene to use custom trace function for conic
    //s.traceNeutronConic = traceNeutronConicWithTables;

    //Add Geometry Here
    Point p1;
    for (int i = 0; i < nummirror; i++) {
		    p1 = makePoint(rfront_inner[i], 0, lStart);
            addFlatEllipse(LStart, LEnd, p1, lStart, lEnd, -mirror_sidelength/2, mirror_sidelength/2, mf, doubleReflections, &s); //inner side of the mirror
		    printf("b[%d] = %f\n", i, rfront_inner[i]);
    }
    if (mirror_width > 0){
        for (int i = 0; i < nummirror; i++){
            p1 = makePoint(rfront_inner[i]+mirror_width, 0, lStart);
            addFlatEllipse(LStart, LEnd, p1, lStart, lEnd, -mirror_sidelength/2, mirror_sidelength/2, mb, doubleReflections, &s); //backside of the above mirror
        }
    }
    addDisk(lEnd, 0.0, 2.0, &s); //neutrons will be propagated important if the are in silicon
	//addEllipsoid(-L, L,p1, -l,+l, 40,&s);
}
#line 10889 "./reverse_test.c"
#undef doubleReflections
#undef mirror_sidelength
#undef mirror_width
#undef mb
#undef mf
#undef nummirror
#undef r_0
#undef lEnd
#undef lStart
#undef LEnd
#undef LStart
#undef sourceDist
#undef transmit
#undef alpha
#undef W
#undef Qc
#undef R0
#undef pTable
#undef s
#undef reflect
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component slit1. */
  SIG_MESSAGE("slit1 (Init)");
#define mccompcurname  slit1
#define mccompcurtype  Slit
#define mccompcurindex 7
#define xmin mccslit1_xmin
#define xmax mccslit1_xmax
#define ymin mccslit1_ymin
#define ymax mccslit1_ymax
#define radius mccslit1_radius
#define xwidth mccslit1_xwidth
#define yheight mccslit1_yheight
#line 50 "/usr/share/mcstas/2.7/tools/Python/mcrun/../mccodelib/../../../optics/Slit.comp"
{
if (xwidth > 0)  { 
  if (!xmin && !xmax) {
    xmax=xwidth/2;  xmin=-xmax;
  } else {
    fprintf(stderr,"Slit: %s: Error: please specify EITHER xmin & xmax or xwidth\n", NAME_CURRENT_COMP); exit(-1);
  }
 }
 if (yheight > 0) { 
   if (!ymin && !ymax) {
     ymax=yheight/2; ymin=-ymax; 
   } else {
     fprintf(stderr,"Slit: %s: Error: please specify EITHER ymin & ymax or ywidth\n", NAME_CURRENT_COMP); exit(-1);
   }
 }
 if (xmin == 0 && xmax == 0 && ymin == 0 && ymax == 0 && radius == 0)
    { fprintf(stderr,"Slit: %s: Warning: Running with CLOSED slit - is this intentional?? \n", NAME_CURRENT_COMP); }

}
#line 10946 "./reverse_test.c"
#undef yheight
#undef xwidth
#undef radius
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component psd_monitor_beforelog. */
  SIG_MESSAGE("psd_monitor_beforelog (Init)");
#define mccompcurname  psd_monitor_beforelog
#define mccompcurtype  PSD_monitor
#define mccompcurindex 8
#define PSD_N mccpsd_monitor_beforelog_PSD_N
#define PSD_p mccpsd_monitor_beforelog_PSD_p
#define PSD_p2 mccpsd_monitor_beforelog_PSD_p2
#define nx mccpsd_monitor_beforelog_nx
#define ny mccpsd_monitor_beforelog_ny
#define filename mccpsd_monitor_beforelog_filename
#define xmin mccpsd_monitor_beforelog_xmin
#define xmax mccpsd_monitor_beforelog_xmax
#define ymin mccpsd_monitor_beforelog_ymin
#define ymax mccpsd_monitor_beforelog_ymax
#define xwidth mccpsd_monitor_beforelog_xwidth
#define yheight mccpsd_monitor_beforelog_yheight
#define restore_neutron mccpsd_monitor_beforelog_restore_neutron
#line 68 "/usr/share/mcstas/2.7/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
  if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }

  if ((xmin >= xmax) || (ymin >= ymax)){
    printf("PSD_monitor: %s: Null detection area !\n"
           "ERROR        (xwidth,yheight,xmin,xmax,ymin,ymax). Exiting",
    NAME_CURRENT_COMP);
    exit(-1);
  }

  PSD_N = create_darr2d(nx, ny);
  PSD_p = create_darr2d(nx, ny);
  PSD_p2 = create_darr2d(nx, ny);

  int i, j;
  for (i=0; i<nx; i++){
    for (j=0; j<ny; j++){
      PSD_N[i][j] = 0;
      PSD_p[i][j] = 0;
      PSD_p2[i][j] = 0;
    }
  }
}
#line 11001 "./reverse_test.c"
#undef restore_neutron
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef filename
#undef ny
#undef nx
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component logspir. */
  SIG_MESSAGE("logspir (Init)");
#define mccompcurname  logspir
#define mccompcurtype  LogSpiral
#define mccompcurindex 9
#define zmin mcclogspir_zmin
#define zmax mcclogspir_zmax
#define ymin mcclogspir_ymin
#define ymax mcclogspir_ymax
#define psi mcclogspir_psi
#define phi_rot mcclogspir_phi_rot
#define precision mcclogspir_precision
#define max_iterations mcclogspir_max_iterations
#define mValue mcclogspir_mValue
#define branches mcclogspir_branches
#define doublesided mcclogspir_doublesided
#define placeholder mcclogspir_placeholder
#line 597 "LogSpiral.comp"
{
	populate_scene(&s, branches, doublesided, zmin, zmax, ymin, ymax, psi, phi_rot, mValue);
	// TODO Test user input for illegal values
}
#line 11041 "./reverse_test.c"
#undef placeholder
#undef doublesided
#undef branches
#undef mValue
#undef max_iterations
#undef precision
#undef phi_rot
#undef psi
#undef ymax
#undef ymin
#undef zmax
#undef zmin
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* Initializations for component psd_monitor. */
  SIG_MESSAGE("psd_monitor (Init)");
#define mccompcurname  psd_monitor
#define mccompcurtype  PSD_monitor
#define mccompcurindex 10
#define PSD_N mccpsd_monitor_PSD_N
#define PSD_p mccpsd_monitor_PSD_p
#define PSD_p2 mccpsd_monitor_PSD_p2
#define nx mccpsd_monitor_nx
#define ny mccpsd_monitor_ny
#define filename mccpsd_monitor_filename
#define xmin mccpsd_monitor_xmin
#define xmax mccpsd_monitor_xmax
#define ymin mccpsd_monitor_ymin
#define ymax mccpsd_monitor_ymax
#define xwidth mccpsd_monitor_xwidth
#define yheight mccpsd_monitor_yheight
#define restore_neutron mccpsd_monitor_restore_neutron
#line 68 "/usr/share/mcstas/2.7/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  if (xwidth  > 0) { xmax = xwidth/2;  xmin = -xmax; }
  if (yheight > 0) { ymax = yheight/2; ymin = -ymax; }

  if ((xmin >= xmax) || (ymin >= ymax)){
    printf("PSD_monitor: %s: Null detection area !\n"
           "ERROR        (xwidth,yheight,xmin,xmax,ymin,ymax). Exiting",
    NAME_CURRENT_COMP);
    exit(-1);
  }

  PSD_N = create_darr2d(nx, ny);
  PSD_p = create_darr2d(nx, ny);
  PSD_p2 = create_darr2d(nx, ny);

  int i, j;
  for (i=0; i<nx; i++){
    for (j=0; j<ny; j++){
      PSD_N[i][j] = 0;
      PSD_p[i][j] = 0;
      PSD_p2[i][j] = 0;
    }
  }
}
#line 11101 "./reverse_test.c"
#undef restore_neutron
#undef yheight
#undef xwidth
#undef ymax
#undef ymin
#undef xmax
#undef xmin
#undef filename
#undef ny
#undef nx
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if(mcdotrace) mcdisplay();
    mcDEBUG_INSTR_END()
  }

} /* end init */

void mcraytrace(void) {
  /* Neutronics-specific defines */
#ifdef NEUTRONICS
extern double mcnx, mcny, mcnz, mcnvx, mcnvy, mcnvz;
extern double mcnt, mcnsx, mcnsy, mcnsz, mcnp;
#endif
  /* End of Neutronics-specific defines */
  /* Copy neutron state to local variables. */
  MCNUM mcnlx = mcnx;
  MCNUM mcnly = mcny;
  MCNUM mcnlz = mcnz;
  MCNUM mcnlvx = mcnvx;
  MCNUM mcnlvy = mcnvy;
  MCNUM mcnlvz = mcnvz;
  MCNUM mcnlt = mcnt;
  MCNUM mcnlsx = mcnsx;
  MCNUM mcnlsy = mcnsy;
  MCNUM mcnlsz = mcnsz;
  MCNUM mcnlp = mcnp;

  mcDEBUG_ENTER()
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define mcabsorb mcabsorbAll
  /* TRACE Component origin [1] */
  mccoordschange(mcposrorigin, mcrotrorigin,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component origin (without coords transformations) */
  mcJumpTrace_origin:
  SIG_MESSAGE("origin (Trace)");
  mcDEBUG_COMP("origin")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbComporigin
  STORE_NEUTRON(1,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[1]++;
  mcPCounter[1] += p;
  mcP2Counter[1] += p*p;
#define mccompcurname  origin
#define mccompcurtype  Progress_bar
#define mccompcurindex 1
#define IntermediateCnts mccorigin_IntermediateCnts
#define StartTime mccorigin_StartTime
#define EndTime mccorigin_EndTime
#define CurrentTime mccorigin_CurrentTime
{   /* Declarations of origin=Progress_bar() SETTING parameters. */
char* profile = mccorigin_profile;
MCNUM percent = mccorigin_percent;
MCNUM flag_save = mccorigin_flag_save;
MCNUM minutes = mccorigin_minutes;
#line 70 "/usr/share/mcstas/2.7/tools/Python/mcrun/../mccodelib/../../../misc/Progress_bar.comp"
{
  double ncount;
  ncount = mcget_run_num();
  if (!StartTime) {
    time(&StartTime); /* compute starting time */
    IntermediateCnts = 1e3;
  }
  time_t NowTime;
  time(&NowTime);
  /* compute initial estimate of computation duration */
  if (!EndTime && ncount >= IntermediateCnts) {
    CurrentTime = NowTime;
    if (difftime(NowTime,StartTime) > 10 && ncount) { /* wait 10 sec before writing ETA */
      EndTime = StartTime + (time_t)(difftime(NowTime,StartTime)
				     *(double)mcget_ncount()/ncount);
      IntermediateCnts = 0;
      fprintf(stdout, "\nTrace ETA ");
      if (difftime(EndTime,StartTime) < 60.0)
        fprintf(stdout, "%g [s] %% ", difftime(EndTime,StartTime));
      else if (difftime(EndTime,StartTime) > 3600.0)
        fprintf(stdout, "%g [h] %% ", difftime(EndTime,StartTime)/3600.0);
      else
        fprintf(stdout, "%g [min] %% ", difftime(EndTime,StartTime)/60.0);
    } else IntermediateCnts += 1e3;
    fflush(stdout);
  }

  /* display percentage when percent or minutes have reached step */
  if (EndTime && mcget_ncount() &&
    (    (minutes && difftime(NowTime,CurrentTime) > minutes*60)
      || (percent && !minutes && ncount >= IntermediateCnts))   )
  {
    fprintf(stdout, "%d ", (int)(ncount*100.0/mcget_ncount())); fflush(stdout);
    CurrentTime = NowTime;

    IntermediateCnts = ncount + percent*mcget_ncount()/100;
    /* check that next intermediate ncount check is a multiple of the desired percentage */
    IntermediateCnts = floor(IntermediateCnts*100/percent/mcget_ncount())*percent*mcget_ncount()/100;
    /* raise flag to indicate that we did something */
    SCATTER;
    if (flag_save) mcsave(NULL);
  }
}
#line 11272 "./reverse_test.c"
}   /* End of origin=Progress_bar() SETTING parameter declarations. */
#undef CurrentTime
#undef EndTime
#undef StartTime
#undef IntermediateCnts
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbComporigin:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(1,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component source [2] */
  mccoordschange(mcposrsource, mcrotrsource,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component source (without coords transformations) */
  mcJumpTrace_source:
  SIG_MESSAGE("source (Trace)");
  mcDEBUG_COMP("source")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompsource
  STORE_NEUTRON(2,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[2]++;
  mcPCounter[2] += p;
  mcP2Counter[2] += p*p;
#define mccompcurname  source
#define mccompcurtype  Arm
#define mccompcurindex 2
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompsource:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(2,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component source_div [3] */
  mccoordschange(mcposrsource_div, mcrotrsource_div,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component source_div (without coords transformations) */
  mcJumpTrace_source_div:
  SIG_MESSAGE("source_div (Trace)");
  mcDEBUG_COMP("source_div")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompsource_div
  STORE_NEUTRON(3,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[3]++;
  mcPCounter[3] += p;
  mcP2Counter[3] += p*p;
#define mccompcurname  source_div
#define mccompcurtype  Source_div
#define mccompcurindex 3
#define thetah mccsource_div_thetah
#define thetav mccsource_div_thetav
#define sigmah mccsource_div_sigmah
#define sigmav mccsource_div_sigmav
#define tan_h mccsource_div_tan_h
#define tan_v mccsource_div_tan_v
#define p_init mccsource_div_p_init
#define dist mccsource_div_dist
#define focus_xw mccsource_div_focus_xw
#define focus_yh mccsource_div_focus_yh
{   /* Declarations of source_div=Source_div() SETTING parameters. */
MCNUM xwidth = mccsource_div_xwidth;
MCNUM yheight = mccsource_div_yheight;
MCNUM focus_aw = mccsource_div_focus_aw;
MCNUM focus_ah = mccsource_div_focus_ah;
MCNUM E0 = mccsource_div_E0;
MCNUM dE = mccsource_div_dE;
MCNUM lambda0 = mccsource_div_lambda0;
MCNUM dlambda = mccsource_div_dlambda;
MCNUM gauss = mccsource_div_gauss;
MCNUM flux = mccsource_div_flux;
#line 118 "/usr/share/mcstas/2.7/tools/Python/mcrun/../mccodelib/../../../sources/Source_div.comp"
{
  double E,lambda,v;

  p=p_init;
  z=0;
  t=0;

  x=randpm1()*xwidth/2.0;
  y=randpm1()*yheight/2.0;
  if(lambda0==0) {
    if (!gauss) {
      E=E0+dE*randpm1();              /*  Choose from uniform distribution */
    } else {
      E=E0+randnorm()*dE;
    }
    v=sqrt(E)*SE2V;
  } else {
    if (!gauss) {
      lambda=lambda0+dlambda*randpm1();
    } else {
      lambda=lambda0+randnorm()*dlambda;
    }
    v = K2V*(2*PI/lambda);
  }

  if (gauss==1) {
      thetah = randnorm()*sigmah;
      thetav = randnorm()*sigmav;
  } else {
      /*find limits of uniform sampling scheme for vertical divergence.
        thetav should be acos(1-2*U) for U\in[0,1]. for theta measured from vertical axis
        we only use a sub-interval for U and measure from horizontal plane.*/
      double sample_lim1,u2;
      sample_lim1=(1-cos(M_PI_2 - focus_ah/2.0*DEG2RAD))*0.5;
      u2=randpm1()*(sample_lim1-0.5) + 0.5;
      thetav = acos(1-2*u2) - M_PI_2;
      thetah = randpm1()*focus_aw*DEG2RAD/2;
  }

  tan_h = tan(thetah);
  tan_v = tan(thetav);

  /* Perform the correct treatment - no small angle approx. here! */
  vz = v / sqrt(1 + tan_v*tan_v + tan_h*tan_h);
  vy = tan_v * vz;
  vx = tan_h * vz;
}
#line 11553 "./reverse_test.c"
}   /* End of source_div=Source_div() SETTING parameter declarations. */
#undef focus_yh
#undef focus_xw
#undef dist
#undef p_init
#undef tan_v
#undef tan_h
#undef sigmav
#undef sigmah
#undef thetav
#undef thetah
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompsource_div:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(3,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component slit [4] */
  mccoordschange(mcposrslit, mcrotrslit,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component slit (without coords transformations) */
  mcJumpTrace_slit:
  SIG_MESSAGE("slit (Trace)");
  mcDEBUG_COMP("slit")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompslit
  STORE_NEUTRON(4,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[4]++;
  mcPCounter[4] += p;
  mcP2Counter[4] += p*p;
#define mccompcurname  slit
#define mccompcurtype  Slit
#define mccompcurindex 4
{   /* Declarations of slit=Slit() SETTING parameters. */
MCNUM xmin = mccslit_xmin;
MCNUM xmax = mccslit_xmax;
MCNUM ymin = mccslit_ymin;
MCNUM ymax = mccslit_ymax;
MCNUM radius = mccslit_radius;
MCNUM xwidth = mccslit_xwidth;
MCNUM yheight = mccslit_yheight;
#line 71 "/usr/share/mcstas/2.7/tools/Python/mcrun/../mccodelib/../../../optics/Slit.comp"
{
    PROP_Z0;
    if (((radius == 0) && (x<xmin || x>xmax || y<ymin || y>ymax))
	|| ((radius != 0) && (x*x + y*y > radius*radius))) {
      RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
      ABSORB;
    }
    else
        SCATTER;
}
#line 11687 "./reverse_test.c"
}   /* End of slit=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompslit:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(4,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component psd_before_optic [5] */
  mccoordschange(mcposrpsd_before_optic, mcrotrpsd_before_optic,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component psd_before_optic (without coords transformations) */
  mcJumpTrace_psd_before_optic:
  SIG_MESSAGE("psd_before_optic (Trace)");
  mcDEBUG_COMP("psd_before_optic")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbComppsd_before_optic
  STORE_NEUTRON(5,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[5]++;
  mcPCounter[5] += p;
  mcP2Counter[5] += p*p;
#define mccompcurname  psd_before_optic
#define mccompcurtype  PSD_monitor
#define mccompcurindex 5
#define PSD_N mccpsd_before_optic_PSD_N
#define PSD_p mccpsd_before_optic_PSD_p
#define PSD_p2 mccpsd_before_optic_PSD_p2
{   /* Declarations of psd_before_optic=PSD_monitor() SETTING parameters. */
int nx = mccpsd_before_optic_nx;
int ny = mccpsd_before_optic_ny;
char* filename = mccpsd_before_optic_filename;
MCNUM xmin = mccpsd_before_optic_xmin;
MCNUM xmax = mccpsd_before_optic_xmax;
MCNUM ymin = mccpsd_before_optic_ymin;
MCNUM ymax = mccpsd_before_optic_ymax;
MCNUM xwidth = mccpsd_before_optic_xwidth;
MCNUM yheight = mccpsd_before_optic_yheight;
MCNUM restore_neutron = mccpsd_before_optic_restore_neutron;
#line 94 "/usr/share/mcstas/2.7/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  PROP_Z0;
  if (x>xmin && x<xmax && y>ymin && y<ymax){
    int i = floor((x - xmin)*nx/(xmax - xmin));
    int j = floor((y - ymin)*ny/(ymax - ymin));
    PSD_N[i][j]++;
    PSD_p[i][j] += p;
    PSD_p2[i][j] += p*p;
    SCATTER;
  }
  if (restore_neutron) {
    RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
  }
}
#line 11821 "./reverse_test.c"
}   /* End of psd_before_optic=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbComppsd_before_optic:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(5,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component flat_ellipse_horizontal [6] */
  mccoordschange(mcposrflat_ellipse_horizontal, mcrotrflat_ellipse_horizontal,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component flat_ellipse_horizontal (without coords transformations) */
  mcJumpTrace_flat_ellipse_horizontal:
  SIG_MESSAGE("flat_ellipse_horizontal (Trace)");
  mcDEBUG_COMP("flat_ellipse_horizontal")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompflat_ellipse_horizontal
  STORE_NEUTRON(6,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[6]++;
  mcPCounter[6] += p;
  mcP2Counter[6] += p*p;
#define mccompcurname  flat_ellipse_horizontal
#define mccompcurtype  FlatEllipse_finite_mirror
#define mccompcurindex 6
#define reflect mccflat_ellipse_horizontal_reflect
#define s mccflat_ellipse_horizontal_s
#define pTable mccflat_ellipse_horizontal_pTable
#define R0 mccflat_ellipse_horizontal_R0
#define Qc mccflat_ellipse_horizontal_Qc
#define W mccflat_ellipse_horizontal_W
#define alpha mccflat_ellipse_horizontal_alpha
#define transmit mccflat_ellipse_horizontal_transmit
{   /* Declarations of flat_ellipse_horizontal=FlatEllipse_finite_mirror() SETTING parameters. */
MCNUM sourceDist = mccflat_ellipse_horizontal_sourceDist;
MCNUM LStart = mccflat_ellipse_horizontal_LStart;
MCNUM LEnd = mccflat_ellipse_horizontal_LEnd;
MCNUM lStart = mccflat_ellipse_horizontal_lStart;
MCNUM lEnd = mccflat_ellipse_horizontal_lEnd;
MCNUM r_0 = mccflat_ellipse_horizontal_r_0;
MCNUM nummirror = mccflat_ellipse_horizontal_nummirror;
MCNUM mf = mccflat_ellipse_horizontal_mf;
MCNUM mb = mccflat_ellipse_horizontal_mb;
MCNUM mirror_width = mccflat_ellipse_horizontal_mirror_width;
MCNUM mirror_sidelength = mccflat_ellipse_horizontal_mirror_sidelength;
MCNUM doubleReflections = mccflat_ellipse_horizontal_doubleReflections;
#line 165 "FlatEllipse_finite_mirror.comp"
{
    dt = (-z + *pointer_lStart)/vz; // first propagate neutron to the entrance window
    if (dt < 0) {
        printf("negative time\n");
    }
    PROP_DT(dt); //propagate neutron to the entrance window of the NMO
    Particle pa = makeParticle(x, y, z, vx, vy, vz, t, sx, sy, sz, silicon, p); //Assume the particle is not arriving in silicon
    if (mirror_width>0){ // if the width of the mirrors is finite neutrons have to know whether they are in silicon or not
        for (int i = 0; i < nummirror; i++){//TODO check if the real frontside is hit
            dt = fabs(rfront_inner[i]); //make sure the mirror distance to check against is positive
            if (dt +mirror_width >= fabs(x)){ //backside of the mirror further out than neutron
                if (dt <= fabs(x)) { // mirror itself closer to the optical axis than the mirro, i.e., we arrive in silicon
                    //ABSORB;
                    pa = makeParticle(x,y,z,vx,vy,vz, t, sx, sy, sz, silicon, p); //create a particle knowing it is in silicon
                    //First we have to refract at the entrance
                    Vec nStart = makeVec(0, 0, 1); //surface normal is oriented in beam direction hopefully
                    Vec init_vec = getParticleVel(pa);
                

                    //printf("before vx = %f\n", pa._vx);
                    //printf("before vy = %f\n", pa._vy);
                    //printf("before vz = %f\n", pa._vz);
                    refractNeutronFlat(&pa, nStart, 0, 0.478);//m_{silicon} =  0.478 laut Peter
                    //printf("after vx = %f\n", pa._vx);
                    //printf("after vy = %f\n", pa._vy);
                    //printf("after vz = %f\n", pa._vz);
                    break;
                    }
                }
            else{   // we do not arrive in Silicon
                //printf("arrived in Air\n");
                    break;
                    }

        }
    }
    
    traceSingleNeutron(&pa, s);//trace the neutron through the mirror assembly
    Vec nEnd = makeVec(0, 0, 1);
    if (pa.silicon==1){
        refractNeutronFlat(&pa, nEnd, 0.478, 0);//TODO
    }


    //Communicate particle state to McStas
    x = pa._x;
    y = pa._y;
    z = pa._z;
    vx = pa._vx;
    vy = pa._vy;
    vz = pa._vz;
    t = pa._t;
    sx = pa._sx;
    sy = pa._sy;
    sz = pa._sz;
    p = pa.w;

    if (pa.absorb)
        ABSORB;

    SCATTER;
}
#line 12013 "./reverse_test.c"
}   /* End of flat_ellipse_horizontal=FlatEllipse_finite_mirror() SETTING parameter declarations. */
#undef transmit
#undef alpha
#undef W
#undef Qc
#undef R0
#undef pTable
#undef s
#undef reflect
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompflat_ellipse_horizontal:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(6,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component slit1 [7] */
  mccoordschange(mcposrslit1, mcrotrslit1,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component slit1 (without coords transformations) */
  mcJumpTrace_slit1:
  SIG_MESSAGE("slit1 (Trace)");
  mcDEBUG_COMP("slit1")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbCompslit1
  STORE_NEUTRON(7,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[7]++;
  mcPCounter[7] += p;
  mcP2Counter[7] += p*p;
#define mccompcurname  slit1
#define mccompcurtype  Slit
#define mccompcurindex 7
{   /* Declarations of slit1=Slit() SETTING parameters. */
MCNUM xmin = mccslit1_xmin;
MCNUM xmax = mccslit1_xmax;
MCNUM ymin = mccslit1_ymin;
MCNUM ymax = mccslit1_ymax;
MCNUM radius = mccslit1_radius;
MCNUM xwidth = mccslit1_xwidth;
MCNUM yheight = mccslit1_yheight;
#line 71 "/usr/share/mcstas/2.7/tools/Python/mcrun/../mccodelib/../../../optics/Slit.comp"
{
    PROP_Z0;
    if (((radius == 0) && (x<xmin || x>xmax || y<ymin || y>ymax))
	|| ((radius != 0) && (x*x + y*y > radius*radius))) {
      RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
      ABSORB;
    }
    else
        SCATTER;
}
#line 12145 "./reverse_test.c"
}   /* End of slit1=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbCompslit1:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(7,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component psd_monitor_beforelog [8] */
  mccoordschange(mcposrpsd_monitor_beforelog, mcrotrpsd_monitor_beforelog,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component psd_monitor_beforelog (without coords transformations) */
  mcJumpTrace_psd_monitor_beforelog:
  SIG_MESSAGE("psd_monitor_beforelog (Trace)");
  mcDEBUG_COMP("psd_monitor_beforelog")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbComppsd_monitor_beforelog
  STORE_NEUTRON(8,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[8]++;
  mcPCounter[8] += p;
  mcP2Counter[8] += p*p;
#define mccompcurname  psd_monitor_beforelog
#define mccompcurtype  PSD_monitor
#define mccompcurindex 8
#define PSD_N mccpsd_monitor_beforelog_PSD_N
#define PSD_p mccpsd_monitor_beforelog_PSD_p
#define PSD_p2 mccpsd_monitor_beforelog_PSD_p2
{   /* Declarations of psd_monitor_beforelog=PSD_monitor() SETTING parameters. */
int nx = mccpsd_monitor_beforelog_nx;
int ny = mccpsd_monitor_beforelog_ny;
char* filename = mccpsd_monitor_beforelog_filename;
MCNUM xmin = mccpsd_monitor_beforelog_xmin;
MCNUM xmax = mccpsd_monitor_beforelog_xmax;
MCNUM ymin = mccpsd_monitor_beforelog_ymin;
MCNUM ymax = mccpsd_monitor_beforelog_ymax;
MCNUM xwidth = mccpsd_monitor_beforelog_xwidth;
MCNUM yheight = mccpsd_monitor_beforelog_yheight;
MCNUM restore_neutron = mccpsd_monitor_beforelog_restore_neutron;
#line 94 "/usr/share/mcstas/2.7/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  PROP_Z0;
  if (x>xmin && x<xmax && y>ymin && y<ymax){
    int i = floor((x - xmin)*nx/(xmax - xmin));
    int j = floor((y - ymin)*ny/(ymax - ymin));
    PSD_N[i][j]++;
    PSD_p[i][j] += p;
    PSD_p2[i][j] += p*p;
    SCATTER;
  }
  if (restore_neutron) {
    RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
  }
}
#line 12279 "./reverse_test.c"
}   /* End of psd_monitor_beforelog=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbComppsd_monitor_beforelog:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(8,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component logspir [9] */
  mccoordschange(mcposrlogspir, mcrotrlogspir,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component logspir (without coords transformations) */
  mcJumpTrace_logspir:
  SIG_MESSAGE("logspir (Trace)");
  mcDEBUG_COMP("logspir")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbComplogspir
  STORE_NEUTRON(9,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[9]++;
  mcPCounter[9] += p;
  mcP2Counter[9] += p*p;
#define mccompcurname  logspir
#define mccompcurtype  LogSpiral
#define mccompcurindex 9
{   /* Declarations of logspir=LogSpiral() SETTING parameters. */
MCNUM zmin = mcclogspir_zmin;
MCNUM zmax = mcclogspir_zmax;
MCNUM ymin = mcclogspir_ymin;
MCNUM ymax = mcclogspir_ymax;
MCNUM psi = mcclogspir_psi;
MCNUM phi_rot = mcclogspir_phi_rot;
MCNUM precision = mcclogspir_precision;
MCNUM max_iterations = mcclogspir_max_iterations;
MCNUM mValue = mcclogspir_mValue;
MCNUM branches = mcclogspir_branches;
MCNUM doublesided = mcclogspir_doublesided;
MCNUM placeholder = mcclogspir_placeholder;
#line 603 "LogSpiral.comp"
{
	part_log n;
	// PROP_Z0;this must not be used to allow neutrons from the other side to hit the logspiral
	// printf("\n z=%f y=%f x=%f vz=%f vy=%f vx=%f \n", z, y, x, vz, vy, vx);
	n = Neutron2Dinit(&n, z, y, x, vz, vy, vx); // puts all of the neutron info into the pointer ptrn pointing to n
	n.customsx = sx;
	n.customsy = sy;
	n.customsz = sz;
	run_scene(s, &n, 10);
	z = n.customz;
	vz = n.customvz;
	y = n.customy;
	vy = n.customvy;
	x = n.customx;
	vx = n.customvx;
	sx = n.customsx;
	sy = n.customsy;
	sz = n.customsz;
}
#line 12420 "./reverse_test.c"
}   /* End of logspir=LogSpiral() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbComplogspir:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(9,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  /* TRACE Component psd_monitor [10] */
  mccoordschange(mcposrpsd_monitor, mcrotrpsd_monitor,
    &mcnlx,
    &mcnly,
    &mcnlz,
    &mcnlvx,
    &mcnlvy,
    &mcnlvz,
    &mcnlsx,
    &mcnlsy,
    &mcnlsz);
  /* define label inside component psd_monitor (without coords transformations) */
  mcJumpTrace_psd_monitor:
  SIG_MESSAGE("psd_monitor (Trace)");
  mcDEBUG_COMP("psd_monitor")
  mcDEBUG_STATE(
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp)
#define x mcnlx
#define y mcnly
#define z mcnlz
#define vx mcnlvx
#define vy mcnlvy
#define vz mcnlvz
#define t mcnlt
#define sx mcnlsx
#define sy mcnlsy
#define sz mcnlsz
#define p mcnlp

#define mcabsorbComp mcabsorbComppsd_monitor
  STORE_NEUTRON(10,
    mcnlx,
    mcnly,
    mcnlz,
    mcnlvx,
    mcnlvy,
    mcnlvz,
    mcnlt,
    mcnlsx,
    mcnlsy,
    mcnlsz,
    mcnlp);
  mcScattered=0;
  mcRestore=0;
  mcNCounter[10]++;
  mcPCounter[10] += p;
  mcP2Counter[10] += p*p;
#define mccompcurname  psd_monitor
#define mccompcurtype  PSD_monitor
#define mccompcurindex 10
#define PSD_N mccpsd_monitor_PSD_N
#define PSD_p mccpsd_monitor_PSD_p
#define PSD_p2 mccpsd_monitor_PSD_p2
{   /* Declarations of psd_monitor=PSD_monitor() SETTING parameters. */
int nx = mccpsd_monitor_nx;
int ny = mccpsd_monitor_ny;
char* filename = mccpsd_monitor_filename;
MCNUM xmin = mccpsd_monitor_xmin;
MCNUM xmax = mccpsd_monitor_xmax;
MCNUM ymin = mccpsd_monitor_ymin;
MCNUM ymax = mccpsd_monitor_ymax;
MCNUM xwidth = mccpsd_monitor_xwidth;
MCNUM yheight = mccpsd_monitor_yheight;
MCNUM restore_neutron = mccpsd_monitor_restore_neutron;
#line 94 "/usr/share/mcstas/2.7/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  PROP_Z0;
  if (x>xmin && x<xmax && y>ymin && y<ymax){
    int i = floor((x - xmin)*nx/(xmax - xmin));
    int j = floor((y - ymin)*ny/(ymax - ymin));
    PSD_N[i][j]++;
    PSD_p[i][j] += p;
    PSD_p2[i][j] += p*p;
    SCATTER;
  }
  if (restore_neutron) {
    RESTORE_NEUTRON(INDEX_CURRENT_COMP, x, y, z, vx, vy, vz, t, sx, sy, sz, p);
  }
}
#line 12554 "./reverse_test.c"
}   /* End of psd_monitor=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex
  /* Label for restoring  neutron */
  mcabsorbComppsd_monitor:
  if (RESTORE) /* restore if needed */
  { RESTORE_NEUTRON(10,
      mcnlx,
      mcnly,
      mcnlz,
      mcnlvx,
      mcnlvy,
      mcnlvz,
      mcnlt,
      mcnlsx,
      mcnlsy,
      mcnlsz,
      mcnlp); }
#undef mcabsorbComp
#undef p
#undef sz
#undef sy
#undef sx
#undef t
#undef vz
#undef vy
#undef vx
#undef z
#undef y
#undef x
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)

  mcabsorbAll:
  mcDEBUG_LEAVE()
  mcDEBUG_STATE(
mcnlx,
mcnly,
mcnlz,
mcnlvx,
mcnlvy,
mcnlvz,
mcnlt,
mcnlsx,
mcnlsy,
mcnlsz,
mcnlp)
  /* Copy neutron state to global variables. */
  mcnx = mcnlx;
  mcny = mcnly;
  mcnz = mcnlz;
  mcnvx = mcnlvx;
  mcnvy = mcnlvy;
  mcnvz = mcnlvz;
  mcnt = mcnlt;
  mcnsx = mcnlsx;
  mcnsy = mcnlsy;
  mcnsz = mcnlsz;
  mcnp = mcnlp;

} /* end trace */

void mcsave(FILE *handle) {
  if (!handle) mcsiminfo_init(NULL);
  /* User component SAVE code. */

  /* User SAVE code for component 'origin'. */
  SIG_MESSAGE("origin (Save)");
#define mccompcurname  origin
#define mccompcurtype  Progress_bar
#define mccompcurindex 1
#define IntermediateCnts mccorigin_IntermediateCnts
#define StartTime mccorigin_StartTime
#define EndTime mccorigin_EndTime
#define CurrentTime mccorigin_CurrentTime
{   /* Declarations of origin=Progress_bar() SETTING parameters. */
char* profile = mccorigin_profile;
MCNUM percent = mccorigin_percent;
MCNUM flag_save = mccorigin_flag_save;
MCNUM minutes = mccorigin_minutes;
#line 115 "/usr/share/mcstas/2.7/tools/Python/mcrun/../mccodelib/../../../misc/Progress_bar.comp"
{
  MPI_MASTER(fprintf(stdout, "\nSave [%s]\n", mcinstrument_name););
  if (profile && strlen(profile) && strcmp(profile,"NULL") && strcmp(profile,"0")) {
    char filename[256];
    if (!strlen(profile) || !strcmp(profile,"NULL") || !strcmp(profile,"0")) strcpy(filename, mcinstrument_name);
    else strcpy(filename, profile);
    DETECTOR_OUT_1D(
        "Intensity profiler",
        "Component index [1]",
        "Intensity",
        "prof", 1, mcNUMCOMP, mcNUMCOMP-1,
        &mcNCounter[1],&mcPCounter[1],&mcP2Counter[1],
        filename);

  }
}
#line 12666 "./reverse_test.c"
}   /* End of origin=Progress_bar() SETTING parameter declarations. */
#undef CurrentTime
#undef EndTime
#undef StartTime
#undef IntermediateCnts
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'psd_before_optic'. */
  SIG_MESSAGE("psd_before_optic (Save)");
#define mccompcurname  psd_before_optic
#define mccompcurtype  PSD_monitor
#define mccompcurindex 5
#define PSD_N mccpsd_before_optic_PSD_N
#define PSD_p mccpsd_before_optic_PSD_p
#define PSD_p2 mccpsd_before_optic_PSD_p2
{   /* Declarations of psd_before_optic=PSD_monitor() SETTING parameters. */
int nx = mccpsd_before_optic_nx;
int ny = mccpsd_before_optic_ny;
char* filename = mccpsd_before_optic_filename;
MCNUM xmin = mccpsd_before_optic_xmin;
MCNUM xmax = mccpsd_before_optic_xmax;
MCNUM ymin = mccpsd_before_optic_ymin;
MCNUM ymax = mccpsd_before_optic_ymax;
MCNUM xwidth = mccpsd_before_optic_xwidth;
MCNUM yheight = mccpsd_before_optic_yheight;
MCNUM restore_neutron = mccpsd_before_optic_restore_neutron;
#line 110 "/usr/share/mcstas/2.7/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  DETECTOR_OUT_2D(
    "PSD monitor",
    "X position [cm]",
    "Y position [cm]",
    xmin*100.0, xmax*100.0, ymin*100.0, ymax*100.0,
    nx, ny,
    &PSD_N[0][0],&PSD_p[0][0],&PSD_p2[0][0],
    filename);
}
#line 12706 "./reverse_test.c"
}   /* End of psd_before_optic=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'psd_monitor_beforelog'. */
  SIG_MESSAGE("psd_monitor_beforelog (Save)");
#define mccompcurname  psd_monitor_beforelog
#define mccompcurtype  PSD_monitor
#define mccompcurindex 8
#define PSD_N mccpsd_monitor_beforelog_PSD_N
#define PSD_p mccpsd_monitor_beforelog_PSD_p
#define PSD_p2 mccpsd_monitor_beforelog_PSD_p2
{   /* Declarations of psd_monitor_beforelog=PSD_monitor() SETTING parameters. */
int nx = mccpsd_monitor_beforelog_nx;
int ny = mccpsd_monitor_beforelog_ny;
char* filename = mccpsd_monitor_beforelog_filename;
MCNUM xmin = mccpsd_monitor_beforelog_xmin;
MCNUM xmax = mccpsd_monitor_beforelog_xmax;
MCNUM ymin = mccpsd_monitor_beforelog_ymin;
MCNUM ymax = mccpsd_monitor_beforelog_ymax;
MCNUM xwidth = mccpsd_monitor_beforelog_xwidth;
MCNUM yheight = mccpsd_monitor_beforelog_yheight;
MCNUM restore_neutron = mccpsd_monitor_beforelog_restore_neutron;
#line 110 "/usr/share/mcstas/2.7/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  DETECTOR_OUT_2D(
    "PSD monitor",
    "X position [cm]",
    "Y position [cm]",
    xmin*100.0, xmax*100.0, ymin*100.0, ymax*100.0,
    nx, ny,
    &PSD_N[0][0],&PSD_p[0][0],&PSD_p2[0][0],
    filename);
}
#line 12745 "./reverse_test.c"
}   /* End of psd_monitor_beforelog=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* User SAVE code for component 'psd_monitor'. */
  SIG_MESSAGE("psd_monitor (Save)");
#define mccompcurname  psd_monitor
#define mccompcurtype  PSD_monitor
#define mccompcurindex 10
#define PSD_N mccpsd_monitor_PSD_N
#define PSD_p mccpsd_monitor_PSD_p
#define PSD_p2 mccpsd_monitor_PSD_p2
{   /* Declarations of psd_monitor=PSD_monitor() SETTING parameters. */
int nx = mccpsd_monitor_nx;
int ny = mccpsd_monitor_ny;
char* filename = mccpsd_monitor_filename;
MCNUM xmin = mccpsd_monitor_xmin;
MCNUM xmax = mccpsd_monitor_xmax;
MCNUM ymin = mccpsd_monitor_ymin;
MCNUM ymax = mccpsd_monitor_ymax;
MCNUM xwidth = mccpsd_monitor_xwidth;
MCNUM yheight = mccpsd_monitor_yheight;
MCNUM restore_neutron = mccpsd_monitor_restore_neutron;
#line 110 "/usr/share/mcstas/2.7/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  DETECTOR_OUT_2D(
    "PSD monitor",
    "X position [cm]",
    "Y position [cm]",
    xmin*100.0, xmax*100.0, ymin*100.0, ymax*100.0,
    nx, ny,
    &PSD_N[0][0],&PSD_p[0][0],&PSD_p2[0][0],
    filename);
}
#line 12784 "./reverse_test.c"
}   /* End of psd_monitor=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  if (!handle) mcsiminfo_close(); 
} /* end save */
void mcfinally(void) {
  /* User component FINALLY code. */
  mcsiminfo_init(NULL);
  mcsave(mcsiminfo_file); /* save data when simulation ends */

  /* User FINALLY code for component 'origin'. */
  SIG_MESSAGE("origin (Finally)");
#define mccompcurname  origin
#define mccompcurtype  Progress_bar
#define mccompcurindex 1
#define IntermediateCnts mccorigin_IntermediateCnts
#define StartTime mccorigin_StartTime
#define EndTime mccorigin_EndTime
#define CurrentTime mccorigin_CurrentTime
{   /* Declarations of origin=Progress_bar() SETTING parameters. */
char* profile = mccorigin_profile;
MCNUM percent = mccorigin_percent;
MCNUM flag_save = mccorigin_flag_save;
MCNUM minutes = mccorigin_minutes;
#line 133 "/usr/share/mcstas/2.7/tools/Python/mcrun/../mccodelib/../../../misc/Progress_bar.comp"
{
  time_t NowTime;
  time(&NowTime);
  fprintf(stdout, "\nFinally [%s: %s]. Time: ", mcinstrument_name, mcdirname ? mcdirname : ".");
  if (difftime(NowTime,StartTime) < 60.0)
    fprintf(stdout, "%g [s] ", difftime(NowTime,StartTime));
  else if (difftime(NowTime,StartTime) > 3600.0)
    fprintf(stdout, "%g [h] ", difftime(NowTime,StartTime)/3660.0);
  else
    fprintf(stdout, "%g [min] ", difftime(NowTime,StartTime)/60.0);
  fprintf(stdout, "\n");
}
#line 12827 "./reverse_test.c"
}   /* End of origin=Progress_bar() SETTING parameter declarations. */
#undef CurrentTime
#undef EndTime
#undef StartTime
#undef IntermediateCnts
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[1]) fprintf(stderr, "Warning: No neutron could reach Component[1] origin\n");
    if (mcAbsorbProp[1]) fprintf(stderr, "Warning: %g events were removed in Component[1] origin=Progress_bar()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[1]);
    if (!mcNCounter[2]) fprintf(stderr, "Warning: No neutron could reach Component[2] source\n");
    if (mcAbsorbProp[2]) fprintf(stderr, "Warning: %g events were removed in Component[2] source=Arm()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[2]);
    if (!mcNCounter[3]) fprintf(stderr, "Warning: No neutron could reach Component[3] source_div\n");
    if (mcAbsorbProp[3]) fprintf(stderr, "Warning: %g events were removed in Component[3] source_div=Source_div()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[3]);
    if (!mcNCounter[4]) fprintf(stderr, "Warning: No neutron could reach Component[4] slit\n");
    if (mcAbsorbProp[4]) fprintf(stderr, "Warning: %g events were removed in Component[4] slit=Slit()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[4]);
  /* User FINALLY code for component 'psd_before_optic'. */
  SIG_MESSAGE("psd_before_optic (Finally)");
#define mccompcurname  psd_before_optic
#define mccompcurtype  PSD_monitor
#define mccompcurindex 5
#define PSD_N mccpsd_before_optic_PSD_N
#define PSD_p mccpsd_before_optic_PSD_p
#define PSD_p2 mccpsd_before_optic_PSD_p2
{   /* Declarations of psd_before_optic=PSD_monitor() SETTING parameters. */
int nx = mccpsd_before_optic_nx;
int ny = mccpsd_before_optic_ny;
char* filename = mccpsd_before_optic_filename;
MCNUM xmin = mccpsd_before_optic_xmin;
MCNUM xmax = mccpsd_before_optic_xmax;
MCNUM ymin = mccpsd_before_optic_ymin;
MCNUM ymax = mccpsd_before_optic_ymax;
MCNUM xwidth = mccpsd_before_optic_xwidth;
MCNUM yheight = mccpsd_before_optic_yheight;
MCNUM restore_neutron = mccpsd_before_optic_restore_neutron;
#line 122 "/usr/share/mcstas/2.7/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  destroy_darr2d(PSD_N);
  destroy_darr2d(PSD_p);
  destroy_darr2d(PSD_p2);
}
#line 12870 "./reverse_test.c"
}   /* End of psd_before_optic=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[5]) fprintf(stderr, "Warning: No neutron could reach Component[5] psd_before_optic\n");
    if (mcAbsorbProp[5]) fprintf(stderr, "Warning: %g events were removed in Component[5] psd_before_optic=PSD_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[5]);
  /* User FINALLY code for component 'flat_ellipse_horizontal'. */
  SIG_MESSAGE("flat_ellipse_horizontal (Finally)");
#define mccompcurname  flat_ellipse_horizontal
#define mccompcurtype  FlatEllipse_finite_mirror
#define mccompcurindex 6
#define reflect mccflat_ellipse_horizontal_reflect
#define s mccflat_ellipse_horizontal_s
#define pTable mccflat_ellipse_horizontal_pTable
#define R0 mccflat_ellipse_horizontal_R0
#define Qc mccflat_ellipse_horizontal_Qc
#define W mccflat_ellipse_horizontal_W
#define alpha mccflat_ellipse_horizontal_alpha
#define transmit mccflat_ellipse_horizontal_transmit
{   /* Declarations of flat_ellipse_horizontal=FlatEllipse_finite_mirror() SETTING parameters. */
MCNUM sourceDist = mccflat_ellipse_horizontal_sourceDist;
MCNUM LStart = mccflat_ellipse_horizontal_LStart;
MCNUM LEnd = mccflat_ellipse_horizontal_LEnd;
MCNUM lStart = mccflat_ellipse_horizontal_lStart;
MCNUM lEnd = mccflat_ellipse_horizontal_lEnd;
MCNUM r_0 = mccflat_ellipse_horizontal_r_0;
MCNUM nummirror = mccflat_ellipse_horizontal_nummirror;
MCNUM mf = mccflat_ellipse_horizontal_mf;
MCNUM mb = mccflat_ellipse_horizontal_mb;
MCNUM mirror_width = mccflat_ellipse_horizontal_mirror_width;
MCNUM mirror_sidelength = mccflat_ellipse_horizontal_mirror_sidelength;
MCNUM doubleReflections = mccflat_ellipse_horizontal_doubleReflections;
#line 228 "FlatEllipse_finite_mirror.comp"
{
    //Mainly Writes Inline Detector Data
    finishSimulation(&s);
}
#line 12912 "./reverse_test.c"
}   /* End of flat_ellipse_horizontal=FlatEllipse_finite_mirror() SETTING parameter declarations. */
#undef transmit
#undef alpha
#undef W
#undef Qc
#undef R0
#undef pTable
#undef s
#undef reflect
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[6]) fprintf(stderr, "Warning: No neutron could reach Component[6] flat_ellipse_horizontal\n");
    if (mcAbsorbProp[6]) fprintf(stderr, "Warning: %g events were removed in Component[6] flat_ellipse_horizontal=FlatEllipse_finite_mirror()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[6]);
    if (!mcNCounter[7]) fprintf(stderr, "Warning: No neutron could reach Component[7] slit1\n");
    if (mcAbsorbProp[7]) fprintf(stderr, "Warning: %g events were removed in Component[7] slit1=Slit()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[7]);
  /* User FINALLY code for component 'psd_monitor_beforelog'. */
  SIG_MESSAGE("psd_monitor_beforelog (Finally)");
#define mccompcurname  psd_monitor_beforelog
#define mccompcurtype  PSD_monitor
#define mccompcurindex 8
#define PSD_N mccpsd_monitor_beforelog_PSD_N
#define PSD_p mccpsd_monitor_beforelog_PSD_p
#define PSD_p2 mccpsd_monitor_beforelog_PSD_p2
{   /* Declarations of psd_monitor_beforelog=PSD_monitor() SETTING parameters. */
int nx = mccpsd_monitor_beforelog_nx;
int ny = mccpsd_monitor_beforelog_ny;
char* filename = mccpsd_monitor_beforelog_filename;
MCNUM xmin = mccpsd_monitor_beforelog_xmin;
MCNUM xmax = mccpsd_monitor_beforelog_xmax;
MCNUM ymin = mccpsd_monitor_beforelog_ymin;
MCNUM ymax = mccpsd_monitor_beforelog_ymax;
MCNUM xwidth = mccpsd_monitor_beforelog_xwidth;
MCNUM yheight = mccpsd_monitor_beforelog_yheight;
MCNUM restore_neutron = mccpsd_monitor_beforelog_restore_neutron;
#line 122 "/usr/share/mcstas/2.7/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  destroy_darr2d(PSD_N);
  destroy_darr2d(PSD_p);
  destroy_darr2d(PSD_p2);
}
#line 12955 "./reverse_test.c"
}   /* End of psd_monitor_beforelog=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[8]) fprintf(stderr, "Warning: No neutron could reach Component[8] psd_monitor_beforelog\n");
    if (mcAbsorbProp[8]) fprintf(stderr, "Warning: %g events were removed in Component[8] psd_monitor_beforelog=PSD_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[8]);
    if (!mcNCounter[9]) fprintf(stderr, "Warning: No neutron could reach Component[9] logspir\n");
    if (mcAbsorbProp[9]) fprintf(stderr, "Warning: %g events were removed in Component[9] logspir=LogSpiral()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[9]);
  /* User FINALLY code for component 'psd_monitor'. */
  SIG_MESSAGE("psd_monitor (Finally)");
#define mccompcurname  psd_monitor
#define mccompcurtype  PSD_monitor
#define mccompcurindex 10
#define PSD_N mccpsd_monitor_PSD_N
#define PSD_p mccpsd_monitor_PSD_p
#define PSD_p2 mccpsd_monitor_PSD_p2
{   /* Declarations of psd_monitor=PSD_monitor() SETTING parameters. */
int nx = mccpsd_monitor_nx;
int ny = mccpsd_monitor_ny;
char* filename = mccpsd_monitor_filename;
MCNUM xmin = mccpsd_monitor_xmin;
MCNUM xmax = mccpsd_monitor_xmax;
MCNUM ymin = mccpsd_monitor_ymin;
MCNUM ymax = mccpsd_monitor_ymax;
MCNUM xwidth = mccpsd_monitor_xwidth;
MCNUM yheight = mccpsd_monitor_yheight;
MCNUM restore_neutron = mccpsd_monitor_restore_neutron;
#line 122 "/usr/share/mcstas/2.7/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  destroy_darr2d(PSD_N);
  destroy_darr2d(PSD_p);
  destroy_darr2d(PSD_p2);
}
#line 12993 "./reverse_test.c"
}   /* End of psd_monitor=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

    if (!mcNCounter[10]) fprintf(stderr, "Warning: No neutron could reach Component[10] psd_monitor\n");
    if (mcAbsorbProp[10]) fprintf(stderr, "Warning: %g events were removed in Component[10] psd_monitor=PSD_monitor()\n"
"         (negative time, miss next components, rounding errors, Nan, Inf).\n", mcAbsorbProp[10]);
  mcsiminfo_close(); 
} /* end finally */
#define magnify mcdis_magnify
#define line mcdis_line
#define dashed_line mcdis_dashed_line
#define multiline mcdis_multiline
#define rectangle mcdis_rectangle
#define box mcdis_box
#define circle mcdis_circle
#define cylinder mcdis_cylinder
#define sphere mcdis_sphere
void mcdisplay(void) {
  printf("MCDISPLAY: start\n");
  /* Components MCDISPLAY code. */

  /* MCDISPLAY code for component 'origin'. */
  SIG_MESSAGE("origin (McDisplay)");
  printf("MCDISPLAY: component %s\n", "origin");
#define mccompcurname  origin
#define mccompcurtype  Progress_bar
#define mccompcurindex 1
#define IntermediateCnts mccorigin_IntermediateCnts
#define StartTime mccorigin_StartTime
#define EndTime mccorigin_EndTime
#define CurrentTime mccorigin_CurrentTime
{   /* Declarations of origin=Progress_bar() SETTING parameters. */
char* profile = mccorigin_profile;
MCNUM percent = mccorigin_percent;
MCNUM flag_save = mccorigin_flag_save;
MCNUM minutes = mccorigin_minutes;
#line 147 "/usr/share/mcstas/2.7/tools/Python/mcrun/../mccodelib/../../../misc/Progress_bar.comp"
{
  
}
#line 13038 "./reverse_test.c"
}   /* End of origin=Progress_bar() SETTING parameter declarations. */
#undef CurrentTime
#undef EndTime
#undef StartTime
#undef IntermediateCnts
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'source'. */
  SIG_MESSAGE("source (McDisplay)");
  printf("MCDISPLAY: component %s\n", "source");
#define mccompcurname  source
#define mccompcurtype  Arm
#define mccompcurindex 2
#line 40 "/usr/share/mcstas/2.7/tools/Python/mcrun/../mccodelib/../../../optics/Arm.comp"
{
  /* A bit ugly; hard-coded dimensions. */
  
  line(0,0,0,0.2,0,0);
  line(0,0,0,0,0.2,0);
  line(0,0,0,0,0,0.2);
}
#line 13062 "./reverse_test.c"
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'source_div'. */
  SIG_MESSAGE("source_div (McDisplay)");
  printf("MCDISPLAY: component %s\n", "source_div");
#define mccompcurname  source_div
#define mccompcurtype  Source_div
#define mccompcurindex 3
#define thetah mccsource_div_thetah
#define thetav mccsource_div_thetav
#define sigmah mccsource_div_sigmah
#define sigmav mccsource_div_sigmav
#define tan_h mccsource_div_tan_h
#define tan_v mccsource_div_tan_v
#define p_init mccsource_div_p_init
#define dist mccsource_div_dist
#define focus_xw mccsource_div_focus_xw
#define focus_yh mccsource_div_focus_yh
{   /* Declarations of source_div=Source_div() SETTING parameters. */
MCNUM xwidth = mccsource_div_xwidth;
MCNUM yheight = mccsource_div_yheight;
MCNUM focus_aw = mccsource_div_focus_aw;
MCNUM focus_ah = mccsource_div_focus_ah;
MCNUM E0 = mccsource_div_E0;
MCNUM dE = mccsource_div_dE;
MCNUM lambda0 = mccsource_div_lambda0;
MCNUM dlambda = mccsource_div_dlambda;
MCNUM gauss = mccsource_div_gauss;
MCNUM flux = mccsource_div_flux;
#line 167 "/usr/share/mcstas/2.7/tools/Python/mcrun/../mccodelib/../../../sources/Source_div.comp"
{
  
  multiline(5, -xwidth/2.0, -yheight/2.0, 0.0,
                xwidth/2.0, -yheight/2.0, 0.0,
                xwidth/2.0,  yheight/2.0, 0.0,
               -xwidth/2.0,  yheight/2.0, 0.0,
               -xwidth/2.0, -yheight/2.0, 0.0);
  if (dist) {
    dashed_line(0,0,0, -focus_xw/2,-focus_yh/2,dist, 4);
    dashed_line(0,0,0,  focus_xw/2,-focus_yh/2,dist, 4);
    dashed_line(0,0,0,  focus_xw/2, focus_yh/2,dist, 4);
    dashed_line(0,0,0, -focus_xw/2, focus_yh/2,dist, 4);
  }
}
#line 13109 "./reverse_test.c"
}   /* End of source_div=Source_div() SETTING parameter declarations. */
#undef focus_yh
#undef focus_xw
#undef dist
#undef p_init
#undef tan_v
#undef tan_h
#undef sigmav
#undef sigmah
#undef thetav
#undef thetah
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'slit'. */
  SIG_MESSAGE("slit (McDisplay)");
  printf("MCDISPLAY: component %s\n", "slit");
#define mccompcurname  slit
#define mccompcurtype  Slit
#define mccompcurindex 4
{   /* Declarations of slit=Slit() SETTING parameters. */
MCNUM xmin = mccslit_xmin;
MCNUM xmax = mccslit_xmax;
MCNUM ymin = mccslit_ymin;
MCNUM ymax = mccslit_ymax;
MCNUM radius = mccslit_radius;
MCNUM xwidth = mccslit_xwidth;
MCNUM yheight = mccslit_yheight;
#line 83 "/usr/share/mcstas/2.7/tools/Python/mcrun/../mccodelib/../../../optics/Slit.comp"
{
  
  if (radius == 0) {
    double xw, yh;
    xw = (xmax - xmin)/2.0;
    yh = (ymax - ymin)/2.0;
    multiline(3, xmin-xw, (double)ymax, 0.0,
              (double)xmin, (double)ymax, 0.0,
              (double)xmin, ymax+yh, 0.0);
    multiline(3, xmax+xw, (double)ymax, 0.0,
              (double)xmax, (double)ymax, 0.0,
              (double)xmax, ymax+yh, 0.0);
    multiline(3, xmin-xw, (double)ymin, 0.0,
              (double)xmin, (double)ymin, 0.0,
              (double)xmin, ymin-yh, 0.0);
    multiline(3, xmax+xw, (double)ymin, 0.0,
              (double)xmax, (double)ymin, 0.0,
              (double)xmax, ymin-yh, 0.0);
  } else {
    circle("xy",0,0,0,radius);
  }
}
#line 13162 "./reverse_test.c"
}   /* End of slit=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'psd_before_optic'. */
  SIG_MESSAGE("psd_before_optic (McDisplay)");
  printf("MCDISPLAY: component %s\n", "psd_before_optic");
#define mccompcurname  psd_before_optic
#define mccompcurtype  PSD_monitor
#define mccompcurindex 5
#define PSD_N mccpsd_before_optic_PSD_N
#define PSD_p mccpsd_before_optic_PSD_p
#define PSD_p2 mccpsd_before_optic_PSD_p2
{   /* Declarations of psd_before_optic=PSD_monitor() SETTING parameters. */
int nx = mccpsd_before_optic_nx;
int ny = mccpsd_before_optic_ny;
char* filename = mccpsd_before_optic_filename;
MCNUM xmin = mccpsd_before_optic_xmin;
MCNUM xmax = mccpsd_before_optic_xmax;
MCNUM ymin = mccpsd_before_optic_ymin;
MCNUM ymax = mccpsd_before_optic_ymax;
MCNUM xwidth = mccpsd_before_optic_xwidth;
MCNUM yheight = mccpsd_before_optic_yheight;
MCNUM restore_neutron = mccpsd_before_optic_restore_neutron;
#line 129 "/usr/share/mcstas/2.7/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  multiline(5,
    (double)xmin, (double)ymin, 0.0,
    (double)xmax, (double)ymin, 0.0,
    (double)xmax, (double)ymax, 0.0,
    (double)xmin, (double)ymax, 0.0,
    (double)xmin, (double)ymin, 0.0);
}
#line 13197 "./reverse_test.c"
}   /* End of psd_before_optic=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'flat_ellipse_horizontal'. */
  SIG_MESSAGE("flat_ellipse_horizontal (McDisplay)");
  printf("MCDISPLAY: component %s\n", "flat_ellipse_horizontal");
#define mccompcurname  flat_ellipse_horizontal
#define mccompcurtype  FlatEllipse_finite_mirror
#define mccompcurindex 6
#define reflect mccflat_ellipse_horizontal_reflect
#define s mccflat_ellipse_horizontal_s
#define pTable mccflat_ellipse_horizontal_pTable
#define R0 mccflat_ellipse_horizontal_R0
#define Qc mccflat_ellipse_horizontal_Qc
#define W mccflat_ellipse_horizontal_W
#define alpha mccflat_ellipse_horizontal_alpha
#define transmit mccflat_ellipse_horizontal_transmit
{   /* Declarations of flat_ellipse_horizontal=FlatEllipse_finite_mirror() SETTING parameters. */
MCNUM sourceDist = mccflat_ellipse_horizontal_sourceDist;
MCNUM LStart = mccflat_ellipse_horizontal_LStart;
MCNUM LEnd = mccflat_ellipse_horizontal_LEnd;
MCNUM lStart = mccflat_ellipse_horizontal_lStart;
MCNUM lEnd = mccflat_ellipse_horizontal_lEnd;
MCNUM r_0 = mccflat_ellipse_horizontal_r_0;
MCNUM nummirror = mccflat_ellipse_horizontal_nummirror;
MCNUM mf = mccflat_ellipse_horizontal_mf;
MCNUM mb = mccflat_ellipse_horizontal_mb;
MCNUM mirror_width = mccflat_ellipse_horizontal_mirror_width;
MCNUM mirror_sidelength = mccflat_ellipse_horizontal_mirror_sidelength;
MCNUM doubleReflections = mccflat_ellipse_horizontal_doubleReflections;
#line 234 "FlatEllipse_finite_mirror.comp"
{
    //Enlarge xy-plane when mcdisplay is ran with --zoom
	magnify("xy");

	//Draw xy-axis contour for Conic Surfaces
	int i;
    for (i = 0; i < s.num_c; i++) {
        double step = (s.c[i].ze-s.c[i].zs)/100;
        double cz;
	    for (cz = s.c[i].zs+step; cz <= s.c[i].ze; cz+= step) {
            double rp = rConic(cz-step,s.c[i]);
            double rc = rConic(cz, s.c[i]);

            line(0,rp,cz-step,0,rc,cz);
            line(0,-rp,cz-step,0,-rc,cz);

            line(rp,0,cz-step,rc,0,cz);
            line(-rp,0,cz-step,-rc,0,cz);
        }
    }

    //Draw xy-axis cross hairs for Disks
    for (i = 0; i < s.num_di; i++) {
        line(s.di[i].r0, 0, s.di[i].z0, s.di[i].r1, 0, s.di[i].z0);
        line(-s.di[i].r0, 0, s.di[i].z0, -s.di[i].r1, 0, s.di[i].z0);
        line(0, s.di[i].r0, s.di[i].z0, 0, s.di[i].r1,s.di[i].z0);
        line(0, -s.di[i].r0, s.di[i].z0, 0, -s.di[i].r1,s.di[i].z0);
    }

}
#line 13264 "./reverse_test.c"
}   /* End of flat_ellipse_horizontal=FlatEllipse_finite_mirror() SETTING parameter declarations. */
#undef transmit
#undef alpha
#undef W
#undef Qc
#undef R0
#undef pTable
#undef s
#undef reflect
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'slit1'. */
  SIG_MESSAGE("slit1 (McDisplay)");
  printf("MCDISPLAY: component %s\n", "slit1");
#define mccompcurname  slit1
#define mccompcurtype  Slit
#define mccompcurindex 7
{   /* Declarations of slit1=Slit() SETTING parameters. */
MCNUM xmin = mccslit1_xmin;
MCNUM xmax = mccslit1_xmax;
MCNUM ymin = mccslit1_ymin;
MCNUM ymax = mccslit1_ymax;
MCNUM radius = mccslit1_radius;
MCNUM xwidth = mccslit1_xwidth;
MCNUM yheight = mccslit1_yheight;
#line 83 "/usr/share/mcstas/2.7/tools/Python/mcrun/../mccodelib/../../../optics/Slit.comp"
{
  
  if (radius == 0) {
    double xw, yh;
    xw = (xmax - xmin)/2.0;
    yh = (ymax - ymin)/2.0;
    multiline(3, xmin-xw, (double)ymax, 0.0,
              (double)xmin, (double)ymax, 0.0,
              (double)xmin, ymax+yh, 0.0);
    multiline(3, xmax+xw, (double)ymax, 0.0,
              (double)xmax, (double)ymax, 0.0,
              (double)xmax, ymax+yh, 0.0);
    multiline(3, xmin-xw, (double)ymin, 0.0,
              (double)xmin, (double)ymin, 0.0,
              (double)xmin, ymin-yh, 0.0);
    multiline(3, xmax+xw, (double)ymin, 0.0,
              (double)xmax, (double)ymin, 0.0,
              (double)xmax, ymin-yh, 0.0);
  } else {
    circle("xy",0,0,0,radius);
  }
}
#line 13315 "./reverse_test.c"
}   /* End of slit1=Slit() SETTING parameter declarations. */
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'psd_monitor_beforelog'. */
  SIG_MESSAGE("psd_monitor_beforelog (McDisplay)");
  printf("MCDISPLAY: component %s\n", "psd_monitor_beforelog");
#define mccompcurname  psd_monitor_beforelog
#define mccompcurtype  PSD_monitor
#define mccompcurindex 8
#define PSD_N mccpsd_monitor_beforelog_PSD_N
#define PSD_p mccpsd_monitor_beforelog_PSD_p
#define PSD_p2 mccpsd_monitor_beforelog_PSD_p2
{   /* Declarations of psd_monitor_beforelog=PSD_monitor() SETTING parameters. */
int nx = mccpsd_monitor_beforelog_nx;
int ny = mccpsd_monitor_beforelog_ny;
char* filename = mccpsd_monitor_beforelog_filename;
MCNUM xmin = mccpsd_monitor_beforelog_xmin;
MCNUM xmax = mccpsd_monitor_beforelog_xmax;
MCNUM ymin = mccpsd_monitor_beforelog_ymin;
MCNUM ymax = mccpsd_monitor_beforelog_ymax;
MCNUM xwidth = mccpsd_monitor_beforelog_xwidth;
MCNUM yheight = mccpsd_monitor_beforelog_yheight;
MCNUM restore_neutron = mccpsd_monitor_beforelog_restore_neutron;
#line 129 "/usr/share/mcstas/2.7/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  multiline(5,
    (double)xmin, (double)ymin, 0.0,
    (double)xmax, (double)ymin, 0.0,
    (double)xmax, (double)ymax, 0.0,
    (double)xmin, (double)ymax, 0.0,
    (double)xmin, (double)ymin, 0.0);
}
#line 13350 "./reverse_test.c"
}   /* End of psd_monitor_beforelog=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  /* MCDISPLAY code for component 'psd_monitor'. */
  SIG_MESSAGE("psd_monitor (McDisplay)");
  printf("MCDISPLAY: component %s\n", "psd_monitor");
#define mccompcurname  psd_monitor
#define mccompcurtype  PSD_monitor
#define mccompcurindex 10
#define PSD_N mccpsd_monitor_PSD_N
#define PSD_p mccpsd_monitor_PSD_p
#define PSD_p2 mccpsd_monitor_PSD_p2
{   /* Declarations of psd_monitor=PSD_monitor() SETTING parameters. */
int nx = mccpsd_monitor_nx;
int ny = mccpsd_monitor_ny;
char* filename = mccpsd_monitor_filename;
MCNUM xmin = mccpsd_monitor_xmin;
MCNUM xmax = mccpsd_monitor_xmax;
MCNUM ymin = mccpsd_monitor_ymin;
MCNUM ymax = mccpsd_monitor_ymax;
MCNUM xwidth = mccpsd_monitor_xwidth;
MCNUM yheight = mccpsd_monitor_yheight;
MCNUM restore_neutron = mccpsd_monitor_restore_neutron;
#line 129 "/usr/share/mcstas/2.7/tools/Python/mcrun/../mccodelib/../../../monitors/PSD_monitor.comp"
{
  multiline(5,
    (double)xmin, (double)ymin, 0.0,
    (double)xmax, (double)ymin, 0.0,
    (double)xmax, (double)ymax, 0.0,
    (double)xmin, (double)ymax, 0.0,
    (double)xmin, (double)ymin, 0.0);
}
#line 13388 "./reverse_test.c"
}   /* End of psd_monitor=PSD_monitor() SETTING parameter declarations. */
#undef PSD_p2
#undef PSD_p
#undef PSD_N
#undef mccompcurname
#undef mccompcurtype
#undef mccompcurindex

  printf("MCDISPLAY: end\n");
} /* end display */
#undef magnify
#undef line
#undef dashed_line
#undef multiline
#undef rectangle
#undef box
#undef circle
#undef cylinder
#undef sphere
/* end of generated C code ./reverse_test.c */
