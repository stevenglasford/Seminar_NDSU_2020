//----------------------------------------
//- 2D example 
//----------------------------------------
//-    v_t + max_a < - f(x,a) . nabla v >  = 0;    f(x,a) = DYN 
//-    v(0,x)=v0(x)
//-    2d advection test 
//----------------------------------------
#include "data_FD_2d_advection_default.h"	// specific to this example
//----------------------------------------
//- problem parameter f(x,a)=dyn
const double dyn0=1.0; //- dynamics x
const double dyn1=1.5; //- dynamics y

const char    NAME[150]         = "data_FD_2d_advection.h, Mar 2018, for 2d test\nAuthor : Boka";
const int     DIM               = 2;   //- space dimension

const int     COMMANDS          = 1;   //- 0: local Hnum function ; 1: Hnum defined with dynamics and distributed_cost functions
const int     OPTIM             = MAXIMUM;

//----------------------------
//- FD method parameters:
//----------------------------
const int     METHOD            = MFD;
const int     TYPE_SCHEME       = ENO2;        //- space scheme : 1 or LF ; 2 or ENO2 ; 3 or ENO3
const int     TYPE_RK           = RK1;         //- time scheme (ENO2 case) : 1 or ENO2RK1, 2 or ENO2RK2,  3 or ENO2RK3

//----------------------------
//- stopping criteria parameters:
//----------------------------
const double  EPSILON           = 0.0;
const int     MAX_ITERATION     = 100000;

const double                   T= 0.50;        //- Terminal time

//----------------------------
//- discretization parameters:
//----------------------------
const double  DT                = 0.0;         //- if DT=0 the program will compute a DT based on space steps.
const double  CFL               = 0.1;

const int     NN=40;
const int     ND[DIM]           = { NN , NN};
const double  XMIN[DIM]         = { -2. , -2.};
const double  XMAX[DIM]         = {  2. ,  2.};           //- bound of the domain

const int     PERIODIC[DIM]     = {  1  ,  1 };           //- periodic mesh (1:periodic, 0:otherwise)
const int     MESH              = 0;                      //- 0 : xi = center of cell ; 1 : xi contains the boundary

const int     cDIM              = 1;
const int     NCD[cDIM]         = { 1 };
const double  UMIN[cDIM]        = { 0.0 };
const double  UMAX[cDIM]        = { 0.0 };

const int     cDIM2             = 1;
const int     NCD2[cDIM2]       = { 1 };
const double  UMIN2[cDIM2]      = { 0.};
const double  UMAX2[cDIM2]      = { 1.};

//-----------------------------
//- BOUNDARY : 0 (Void, FD only) or 1 (Dirichlet, V=g_border, for FD/SL) or 2 (Vx=g_bordermix, for FD), or 3 (Vxx=0, for FD/SL) 
//-----------------------------
const int BOUNDARY=0;

// Dirichlet boundary condition (case BOUNDARY=1)
const double  VBORD             = 0.5; 
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
const int     COMPUTE_TOPT      = 0;
const int     TOPT_TYPE         = 0;       //- 0= min time pb., 1= exit time pb.

const int     SAVE_VF_ALL       = 0;       //- 1 to save "VFn.dat" every  SAVE_VF_ALL_STEP iterations
const int     SAVE_VF_ALL_STEP  = 10;      //- (a large number will save only VF0.dat and final: VF1.dat)
const int     SAVE_VF_FINAL     = 1;

const int     CHECK_ERROR       = 1;       //- 1 to compute errors every CHECK_ERROR_STEP iterations
const int     CHECK_ERROR_STEP  = 100;
const int     CHECK_NEG_POINTS  = 0;       //- 1/0 : to print number of negative point in domain and domain + boundary (default : 0)

//-------------------
//- "coupe" : (cut into some plane ==> results in files "coupe{}{ex}{topt}.dat"
//-------------------
const int     COUPE_DIMS[DIM] = {1 ,1 }; //- put 1 to keep the dimensions in the final savings.
const double  COUPE_VALS[DIM] = {0.,0.}; //- if COUPE_DIM[d]=0 then make precise the value of the cut in direction x_d.

//---------------
//- initial data
//---------------
/*
const  double   X0=0.0, Y0=0.0, RAYON=0.5;
inline double   q(double x)
{
  return MIN(VBORD,x-RAYON);
}
*/

inline double   v0(const double *x)
{
  //double nx=sqrt( (arg[0]-X0)*(arg[0]-X0) + (arg[1]-Y0)*(arg[1]-Y0) );
  double y[2];
  y[0]=x[0]-floor((x[0]-XMIN[0])/(XMAX[0]-XMIN[0]))*(XMAX[0]-XMIN[0]);
  y[1]=x[1]-floor((x[1]-XMIN[1])/(XMAX[1]-XMIN[1]))*(XMAX[1]-XMIN[1]);
  double X0=0.0;
  double xx=y[0]+2*y[1]-X0;
  double res=sin(pi*xx/2.0);
  //double res=max(1.0-abs(xx/1.0),0.0);
  return res;
}

//------------------------------------------------------------------------------------------------------
//- dynamics and distributed cost functions used in the case of COMMANDS=0 or 1 (assumes only 1 control)
//------------------------------------------------------------------------------------------------------
inline void dynamics(const double* x, C u, double t, double* res)
{
  res[0]=  -dyn0;
  res[1]=  -dyn1;
}

inline double distributed_cost(const double* x, C u, double t){
  return 0.0;
}

inline void feedback(double t, const double* x, const double* p, C& u){  //- only needed for NBEE-OC
  //- input  : t,x,p
  //- output : u, that maximizes  max_u (- p. dynamics(t,x,u) - dist_cost(t,x,u)) 
  //- Example 
  //u[0]=1.0;
}

//---------------------------------------------------------------
//- if COMMANDS=2 (2 player games): dynamics and distributed cost
//---------------------------------------------------------------
inline void dynamics2(const double* x, C u, C u2, double t, double* res) { return; }
inline double distributed_cost2(const double* x, C u, C u2, double t) { return 0.; }


//---------------------------
//- Exact solution (if known)
//---------------------------
/*
inline double Vex(double t, const double* x)
{
    double nx= sqrt(arg[0] *arg[0] + arg[1] *arg[1]);
    if      (OPTIM==MAXIMUM) return q( MAX(0, nx-t) );
    else if (OPTIM==MINIMUM) return q( MAX(0, nx+t) );
    else {printf("this OPTIM not programmed!\n"); return q( MAX(0, nx-t) );}
}
*/

inline double Vex(double t, const double* x){
  double y[2]={x[0]-dyn0*t,x[1]-dyn1*t};
  return v0(y);
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
double        min_TRAJ_STOP     =   0.00;  //- to stop traj reconstruction when val(x) <=min (val=topt here)
double        max_TRAJ_STOP     = 100.00;  //- to stop traj reconstruction when val(x) >=max (val=topt here)
int           TARGET_STOP       = 1;       //- to stop traj reconstruction when g_target(x)<=0
inline double g_target(double t, const double* x){ return v0(x);}  
// adverse control case (only for COMMANDS=2 & TRAJ_METHOD=0)
int           ADVERSE_METHOD    = 0;
inline void u2_adverse(double t, const double* x, C& u){}

//-----------------------------------------------
//- PARAMETERS FOR ADVANCED USERS
//- not necessary to modify these for basic usage
//------------------------------------------------

//--------------------
//- Numerical Hamiltonian and stability constants (case COMMANDS=0);
//--------------------
inline void compute_Hconst(double* aMAX, double t){}
inline double Hnum(const double t, const double* x, const double vi, const double* v){return 0.0;}

//--------------------
//- OBSTACLE g
//--------------------
const int OBSTACLE=0;
inline double g_obstacle(double t, const double* arg){ return 0.;}

//--------------------------------
//- border: number of ghost cells in each direction (default is {2,2,...} for FD)
//--------------------------------
int BORDERSIZE[DIM] = {3,3};

//--------------------
//- Special parameter for method SL
//--------------------
const int     P_INTERMEDIATE    = 1;    //- number of discretisation steps for approximating each trajectory in method MSL 

