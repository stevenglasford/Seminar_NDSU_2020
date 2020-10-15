//----------------------------------------
//- Eikonal Equation example, dimension 2
//----------------------------------------
//-    u_t + c(x) || \nabla u || = 0 or  u_t + max_a ( -f(x,a) nabla u ) = 0
//-    u(0,x)=u0(x)
//-    ex_1b : using numerical hamiltonian Hnum (COMMANDS=0)
//----------------------------------------
#include "data_default.h"	// do not modify for basic examples
//----------------------------------------

const char    NAME[150]         = "data_FD_2d_ex1_advanced.h (Eikonal equation example) Feb. 2013\nAuthors   : O. Bokanowski, H. Zidani";
const int     DIM               = 2;   //- space dimension

const int     COMMANDS          = 0;   //- 0: local Hnum function ; 1: Hnum defined with dynamics and distributed_cost functions
const int     OPTIM             = MAXIMUM;


//----------------------------
//- DF method parameters:
//----------------------------
const int     METHOD            = MFD;
const int     TYPE_SCHEME       = ENO2;        //- space scheme : 1 or LF ; 2 or ENO2
const int     TYPE_RK           = RK1;         //- time scheme (ENO2 case) : 1 or ENO2RK1, 2 or ENO2RK2,  3 or ENO2RK3

//----------------------------
//- stopping criteria parameters:
//----------------------------
const double  EPSILON           = 0.0;
const int     MAX_ITERATION     = 100000;
const double                   T= 0.5;         //- Terminal time


//----------------------------
//- discretization parameters:
//----------------------------
const double  DT                = 0.;          //- if DT=0 the program will compute a DT based on space steps.
const double  CFL               = 0.5;

const int     NN=100;
const int     ND[DIM]           = { NN , NN};
const double  XMIN[DIM]         = { -2. , -2.};
const double  XMAX[DIM]         = {  2. ,  2.};           //- bound of the domain

const int     PERIODIC[DIM]     = {  0  ,  0 };           //- periodic mesh (1:periodic, 0:otherwise)
const int     MESH              = 1;                      //- 0 : xi = center of cell ; 1 : xi contains the boundary

const int     cDIM              = 1;
const int     NCD[cDIM]         = { 1 };
const double  UMIN[cDIM]        = { 0.};
const double  UMAX[cDIM]        = { 2.*pi};

const int     cDIM2             = 1;
const int     NCD2[cDIM2]       = { 1 };
const double  UMIN2[cDIM2]      = { 0.};
const double  UMAX2[cDIM2]      = { 1.};

//-----------------------------
//- BOUNDARY : 0 (Void, FD only) or 1 (Dirichlet, V=g_border, for FD/SL) or 2 (Vx=g_bordermix, for FD), or 3 (Vxx=0, for FD/SL) 
//-----------------------------
const int BOUNDARY=1;

// Dirichlet boundary condition (case BOUNDARY=1)
const double  VBORD             = 0.4; 
double g_border(double t, const double* arg){
  return VBORD;
}

// Mixed Neumann bc :  ux=g(t,x,u) (case BOUNDARY=2)
double g_bordermix(double t, const double* x, double val){return 0.0;}


//-----------------------------
//- mainloop parameters
//-----------------------------
const int     COMPUTE_MAIN_LOOP = 1;
const int     COMPUTE_VEX       = 1;
const int     COMPUTE_TOPT      = 1;
const int     TOPT_TYPE    = 0;// 0= min time , 1= exit time

const int     SAVE_VF_ALL       = 0;       //- 1 to save "VFn.dat" every  SAVE_VF_ALL_STEP iterations
const int     SAVE_VF_ALL_STEP  = 10;      //- (a large number will save only VF0.dat and final: VF1.dat)
const int     SAVE_VF_FINAL     = 1;

const int     CHECK_ERROR       = 1;       //- 1 to compute errors every CHECK_ERROR_STEP iterations
const int     CHECK_ERROR_STEP  = 10;


//-------------------
//- "coupe" : used only for exemples with d>=3 (cut into some plane)
//-------------------
const int     COUPE_DIMS[DIM] = {1 ,1 };
const double  COUPE_VALS[DIM] = {0.,0.};



inline void dynamics(const double* x, C u, double t, double* res)
{
  res[0]=  cos(u[0]);
  res[1]=  sin(u[0]);
}


inline double distributed_cost(const double* arg, C u, double t)
{
  return 0.;
}


//- unused functions since only 1 control here (COMMANDS=1)
inline void dynamics2(const double* x, C u, C u2, double t, double* res) { return; }
inline double distributed_cost2(const double* arg, C u, C u2, double t) { return 0.; }


//---------------
//- initial data
//---------------
const  double   X0=0.0, Y0=0.0, RAYON=0.5;
inline double   q(double x)
{
  return MIN(VBORD,x-RAYON);
}

inline double   v0(const double * arg)
{
   double nx=sqrt( (arg[0]-X0)*(arg[0]-X0) + (arg[1]-Y0)*(arg[1]-Y0) );
   return q(nx);
}


//---------------------------
//- Exact solution (if known)
//---------------------------
inline double Vex(double t, const double* arg)
{
    double nx= sqrt(arg[0] *arg[0] + arg[1] *arg[1]);
    if      (OPTIM==MAXIMUM) return q( MAX(0, nx-t) );
    else if (OPTIM==MINIMUM) return q( MAX(0, nx+t) );
    else {printf("this OPTIM not programmed!\n"); return q( MAX(0, nx-t) );}
}


//--------------------
//- Hamiltonian (to be used with the LF numerical hamiltonian here)
//--------------------
inline double H(const double* x, const double* p)
{
  double z;
  double c=1.;  //- speed
  z=c*sqrt(p[0]*p[0]+ p[1]*p[1]);
  return z;
}

//--------------------
//- Stability constants to be used for the numerical hamiltonian.
//--------------------
inline void compute_Hconst(double* aMAX, double t)
{
    aMAX[0]=1.;
    aMAX[1]=1.;
}

//--------------------
//- Numerical Hamiltonian (Lax Friedriech) in the case H is defined
//--------------------
inline double Hnum(const double t, const double* x, const double vi, const double* Dv)
{
  double z=0.;
  double p[DIM];
  double amax[DIM];
  int i;
  for(i=0;i<DIM;i++) {
    amax[i]=1.;
    p[i]=(Dv[2*i] + Dv[2*i+1])/2.;
    z += amax[i]*(Dv[2*i+1] - Dv[2*i])/2.;
  }
  return H(x,p)-z;
}


//-----------------------------------------------
//- PARAMETERS : compute trajectory 
//-----------------------------------------------
// method of trajectory reconstruction and starting time
int           TRAJ_METHOD       = 0; 
double        time_TRAJ_START   = 0.00;
// initial point
const int     TRAJPT            = 0;  
const double  initialpoint[TRAJPT*DIM] = {}; //- Initial point for HJB
// stopping criteria
double        min_TRAJ_STOP     = 0.00;    //- to stop traj reconstruction when val(x) <=min (val=topt here)
double        max_TRAJ_STOP     = 100.00;  //- to stop traj reconstruction when val(x) >=max (val=topt here)
int           TARGET_STOP       = 1;       //- to stop traj reconstruction when g_target(x)<=0
inline double g_target(double t, const double* x){ return v0(x);}  
// adverse control case (only for COMMANDS=2 & TRAJ_METHOD=0)
int           ADVERSE_METHOD    = 0;
void u2_adverse(double t, const double* x, C& u){}


//-----------------------------------------------
//- PARAMETERS FOR ADVANCED USERS
//- not necessary to modify these for basic usage
//------------------------------------------------

//--------------------------------
//- PARAMETERS FOR ADVANCED USERS
//--------------------------------
const int     EXTERNALV0        = 0;           //- if 1 then starts computation from data in VF.dat and not from the v0 function
const int     PRECOMPUTE_COORDS = 1;           //- to precompute mesh coordinates (faster but needs more memory).

//--------------------
//- OBSTACLE g 
//--------------------
const int OBSTACLE=0;
inline double g_obstacle(double t, const double* arg){ return 0.;}

//--------------------------------
//- border: number of ghost cells in each direction (default is {2,2,...} for FD)
//--------------------------------
int BORDERSIZE[DIM] = {2,2};

//--------------------
//- Special parmeter for method SL
//--------------------
const int     P_INTERMEDIATE    = 1;    //- number of discretisation steps for approximating each trajectory in method MSL 

