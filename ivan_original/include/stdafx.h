// stdafx.h : include file seldom to be modified.

#ifndef _STDAFX_H
#define _STDAFX_H

#define ENCOURS 0

#define PRINTBUG 0
#define PRINTMEM 0
#define PRINTBUGPARA 0


#include <unistd.h> // -for new gcc 87
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <stdexcept>
#include <cmath>
#include <cassert>
#include <cstring>
#include <iostream>
#include <fstream>
#include <istream>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <valarray>
#include <list>
#include <limits>
#include <utility>
#include <string>
#include<sys/time.h>
#include<sys/types.h>
#include <omp.h>

using namespace std;
typedef double* C;


//- mesh parameters container
struct mrp {
  unsigned  DIM;
  int       PRECOMPUTE_COORDS;
  int       MESH_TYPE;
  unsigned  *SIZES;
  double    *XMIN;
  double    *XMAX;
  int       *PERIODIC;
  int       *BORDERSIZE;        // [2016] bordersize[d] gives the number of boundary points in dir.
  int       OBSTACLE;
  double    (*G_OBSTACLE)       (double, const double*);
  int       OBSTACLE_TILDE;
  double    (*G_OBSTACLE_TILDE) (double, const double*);
  int       PRECOMPUTE_OBSTACLE;
  int       COMPUTE_IN_SUBDOMAIN;
  double    (*G_DOMAIN)         (const double*);
};

//- control parameters container
struct cp {
  int       DIM;
  int       *SIZES;
  int       NC;
  double    *UMIN;
  double    *UMAX;
};

//- parallel parameters container
struct pp {
  unsigned  DIM;
  unsigned  *SIZES;
  unsigned  MY_ID;
  unsigned  PROC_NUMS;
  unsigned  OMP_NUM_THREADS;
};

//- general parameters container
struct gp {
  int       METHOD;
  const char* DATA_NAME;
  int       EXTERNALV0;
  int       COMPUTE_MAIN_LOOP;
  int       OPTIM;
  int       COMPUTE_TOPT;
  int       TOPT_TYPE;               //- july 2015
  int       COMPUTE_VEX;
  int       SAVE_VF_ALL;
  int       SAVE_VF_ALL_STEP;
  int       SAVE_VF_FINAL;
  int       SAVE_VF_FINAL_ONSET;
  int       VALUE_PB;             //- savings for value pb
  int       SAVE_VALUE_ALL;       //- savings for value pb
  int       SAVE_VALUE_ALL_STEP;  //- savings for value pb
  int       SAVE_VALUE_FINAL;     //- savings for value pb
  int       FORMAT_FULLDATA;      //- 1 for full data savings / 0 for short
  int       CHECK_ERROR;
  int       CHECK_ERROR_STEP;
  double    C_THRESHOLD;

  double    T;
  double    DT;
  int       BOUNDARY;
  int       MAX_ITERATION;
  double    EPSILON;

  void      (*DYNAMICS)         (const double*, C, double, double*);
  void      (*DYNAMICS2)        (const double*, C, C, double, double*);
  double    (*DISTRIBUTED_COST) (const double*, C, double);
  double    (*DISTRIBUTED_COST2)(const double*, C, C, double);
  double    (*V0)               (const double*);
  double    (*VEX)              (double, const double*);
  double    (*G_BORDER)         (double, const double*);
  double    (*G_BORDERMIX)      (double, const double*, double);

  // string    FILE_PREFIX;              //- prefix for output filenames 
  char*     FILE_PREFIX;              //- prefix for output filenames 
  char*     VFile;
  char*     VFALLFile;
  char*     tminFile;
  char*     dataFile;

  //- ICI 2015
  //char      VFile[100];
  //char*     VALLFile;
  //char*     tminFile;
  //char*     dataFile;
};

//- FD parameters container
struct fdp {
  int       TYPE_SCHEME;      //space scheme : 1: LF ; 2: ENO2
  int       TYPE_RK;          //time scheme  : 1: RK1 ; 2: RK2 ; 3: RK3
  int       COMMANDS;
  double    CFL;
  double    (*HNUM)             (const double t, const double* x, const double vi, const double* Dv);
  double    (*DISCOUNT_FACTOR)  (const double*);
  void      (*COMPUTEAMAX)      (double*, double);
};

//- SL parameters container
struct slp {
  int       TYPE_SCHEME;      //schéma espace : 1: STA(TIONNAIRE) ; 2: EVO(LUTIVE)
  int       TYPE_RK;          // 1: RK1_EULER ; 2: RK2_HEUN ; 3: RK2_PM (Middle Point)
  int       ORDER;
  int       TYPE_STA_LOOP;    // 1:normal : 2: 2*dim mesh loop per mainloop's iteration
  int       INTERPOLATION;    // 1:BILINEAR ; 2:PRECOMPBL ; 3:DIRPERDIR
  int       P_INTERMEDIATE;
  double    PARAMP;           // number of diffusions for second order case
  double    (*FUNCR)            (const double*, C, double);
  void      (*FUNCY)            (const double*, int, double, C, double, double, double*);
  double    (*DISCOUNT_FACTOR)  (const double*);
};

//- UB parameters container
struct ubp {
  double    CFL;
  int       SUB_N;
  int       TARGET;
  int       DEB;
  int       DEC;
  double    DIST_THRES_COEFF;
  double    (*DISTOTARG)        (const double*);
  void      (*TARGET_BOUNDARIES)(double*, double*);
};

//- DG parameters container
struct dgp {
  int       TYPE_SCHEME;      //space scheme : 1: LF ; 2: ENO2
  int       TYPE_RK;          //time scheme  : 1: RK1 ; 2: RK2 ; 3: RK3
  int       COMMANDS;
  double    CFL;
  double    (*HNUM)             (const double*, const double*, const double);
  void      (*COMPUTEAMAX)      (double*, double);
};

//- post processing parameters container
struct ppp {
  int       COMPUTE_TRAJ_POINTS;
  double    *POINTS;
  int       TRAJ_METHOD;      // 2014
  double    time_TRAJ_START;  // 2014

  double    min_TRAJ_STOP;    // 2017
  double    max_TRAJ_STOP;    // 2017
  int       TARGET_STOP;      // 2017
  double    (*G_TARGET)         (double, const double*); // 2017

  int       ADVERSE_METHOD;   // 2017
  void      (*U2_ADVERSE)       (double, const double*, C&); // 2017
  //void      (*U2_ADVERSE)       (double, const double*, double*); // 2017

  int       PRINTTRAJ;        // 2017 

  int       *COUPE_DIMS;
  double    *COUPE_VALS;
};


//- UB chained list node structure
struct tris {
  double val_init;
  double val_its;
  double val_its2;
  double val;
};


//Discontinued optimization constants
const int       MINIMUM     = 1;
const int       MAXIMUM     = 2;
const int       MINMAX      = 3;
const int       MAXMIN      = 4;

//used method constants
const int       MFD         = 1;
const int       MSL         = 2;
const int       MUB         = 3;
const int       MDG         = 4;

//SL constants
const int       STA         = 1;
const int       EVO         = 2;

const int       RK1_EULER   = 1;
const int       RK2_HEUN    = 2;
const int       RK2_PM      = 3;

const int       NORMAL      = 1;
const int       SPECIAL     = 2;

//- BOKA COMMENT 2014
//const int       ONEVAL      = 1;
//const int       FUNC        = 2;
//const int       PERIOD      = 3;

const int       BILINEAR    = 1;
const int       PRECOMPBL   = 2;
const int       DIRPERDIR   = 3;

//FD constants
const int       LF          = 1;
const int       ENO2        = 2;
const int       RK1         = 1;
const int       RK2         = 2;
const int       RK3         = 3;


//Discontinued Galerkin constants
//const int       DGP        = 1;
//const int       DGQ        = 2;

const double    pi          = 3.14159265359;
const double    INF         = 1.e5;

// -------------------
// - output filenames
// -------------------
//- ICI 2015 : filenames 
//- complete filename will be {OUTPUTdir}{FILE_PREFIX}{FILENAME}.dat where FILENAME is one of :
//- [if numbering (for VF or traj) the filename may be ended with -n.dat instead of .dat]
const char      OUTPUTdir   []= "OUTPUT/";  // will be completed by ".dat" or "n.dat" extention 
const char      OUTPUTdirWin[]= "OUTPUT\\";  // will be completed by ".dat" or "n.dat" extention

const char      VFile0      []= "VF";       // will be completed by ".dat" or "n.dat" extention 
const char      XFile_user0 []= "X_user";   // will be completed by ".dat" or "n.dat" extention 
const char      VFile_user0 []= "XV_user";  // will be completed by ".dat" or "n.dat" extention 
const char      VEXFile0    []= "VEX";      // will be completed by ".dat" or "n.dat" extention 
const char      tminFile0   []= "topt";     // will be completed by ".dat" extention 
const char      dataFile0   []= "data";     // will be completed by ".dat" extention 
const char      ValueFile0  []= "Value";    // will be completed by ".dat" extention 
const char      coupeFile0  []= "coupe";    // will be completed by ".dat" extention 
const char      coupeexFile0[]= "coupeex";  // will be completed by ".dat" extention 
const char      trajFile0   []= "traj";     // will be completed by "-n.dat" extention for nth trajectory
const char      trajFileMore0[]="moretraj"; // will be completed by "-n.dat" extention for nth trajectory
const char      DtFile0     []= "Dt";       // will be completed by ".dat" extention 
const char      successTrajectoriesFile0 []= "successTrajectories";// will be completed by ".dat" extention 
//const char      VFALLFile0  []= "VFALL";  // not used

const char      toptcoupeFile0  []= "toptcoupe";    // [DEC 2016]

/*
// -------------------
// - output filenames
// -------------------
const char      VFile      []= "OUTPUT/VF.dat";
const char      VEXFile    []= "OUTPUT/VEX.dat";  
const char      VFALLFile  []= "OUTPUT/VFALL.dat";
const char      tminFile   []= "OUTPUT/tmin.dat";
const char      dataFile   []= "OUTPUT/data.dat";
const char      ValueFile  []= "OUTPUT/Value.dat";
const char      coupeFile  []= "OUTPUT/coupe.dat";
const char      coupeexFile[]= "OUTPUT/coupeex.dat";
const char      trajFile   []= "OUTPUT/traj";  // will be completed by "-n.dat" extention for nth trajectory
const char      DtFile     []= "OUTPUT/Dt.dat";  // will be completed by "-n.dat" extention for nth trajectory
const char      successTrajectoriesFile []= "OUTPUT/successTrajectories.dat";
*/

//const char      VTKFile[]   = "OUTPUT/res";

const char      fFORMAT[]   = "%15.10f"; //- Ex:  "%15.10f", ou "%10.4e"
const char      iFORMAT[]   = "%2i";

inline double   MAX(double a, double b) { return a>=b ?  a : b; }
inline double   MIN(double a, double b) { return a<=b ?  a : b; }
inline double   ABS(double a)           { return a<0. ? -a : a; }
inline double   SIGN(double a)          { return a<0. ? -1.: 1.; }


#endif

