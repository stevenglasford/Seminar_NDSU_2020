//----------------------------------------
//- Rotation example , dimension 2
//----------------------------------------
//-    v_t + < f, nabla v >  = 0;    f= -(-y,x);
//-    v(0,x)=v0(x)
//----------------------------------------
#include "data_default.h"	// do not modify for basic examples
//----------------------------------------

const char    NAME[150]         = "data_SL_2d_rot.h, Advection-Rotation Example, 2014 (O. Bokanowski, A. Desilles, H. Zidani)";
const int     DIM               = 2;   //- space dimension

const int     COMMANDS          = 0;   //- 0: local Hnum function ; 1: Hnum defined with dynamics and distributed_cost functions
const int     OPTIM             = MAXIMUM;



//----------------------------
//- SL method parameters:
//----------------------------
const int     METHOD            = 3; //MSL;         //- Method : Semi-Lagrangian
const int     TYPE_SCHEME       = EVO;         //- space scheme  : 1 or STA (static); 2 or EVO
      int     TYPE_RK           = RK2_HEUN;    //- time scheme   : 1 or RK1_EULER, 2 or RK2_HEUN,  3 or RK2_PM

//--------------------
//- Special parmeter for method SL
//--------------------
      int     P_INTERMEDIATE    = 1;           //- number of discretisation steps for approximating each trajectory in method MSL 

//----------------------------
//- stopping criteria parameters:
//----------------------------
const double  EPSILON           = 0.0;
const int     MAX_ITERATION     = 100000;
double                         T= 0.25;        //- Terminal time


//----------------------------
//- discretization parameters:
//----------------------------
// (for MFD, the program will compute a DT based on space steps).
const int nn=4;
      double                  DT= 1/float(10*nn)/2;  //- if DT=0, MSL will not work

const int                     NN= 25*nn;
      int     ND[DIM]           = { NN , NN};
      double  XMIN[DIM]         = { -2. , -2.};      //- bounds of the domain
      double  XMAX[DIM]         = {  2. ,  2.};      

const int     PERIODIC[DIM]     = {  0  ,  0 };           //- periodic mesh (1:periodic, 0:otherwise)
const int     MESH              = 1;                      //- 0 : xi = center of cell ; 1 : xi contains the boundary

const int     cDIM              = 1;
const int     NCD[cDIM]         = { 1 };
const double  UMIN[cDIM]        = { 0.   };
const double  UMAX[cDIM]        = { 2.*pi};

const int     cDIM2             = 1;
const int     NCD2[cDIM2]       = { 1 };
const double  UMIN2[cDIM2]      = { 0.};
const double  UMAX2[cDIM2]      = { 1.};

//-----------------------------
//- BOUNDARY : 0 (Void, FD only) or 1 (Dirichlet, V=g_border, for FD/SL) or 2 (Vx=g_bordermix, for FD), or 3 (Vxx=0, for FD/SL) 
//-----------------------------
const int BOUNDARY=3;

// Dirichlet boundary condition (case BOUNDARY=1)
const double  VBORD             = 0.3; 
double g_border(double t, const double* arg){
  return VBORD;
}

// Mixed Neumann bc :  ux=g(t,x,u) (case BOUNDARY=2)
double g_bordermix(double t, const double* x, double val){return 1.0;}

// xterm printings of trajectory (default should be 0)
//int           PRINTTRAJ         = 0;        


//-----------------------------
//- mainloop parameters
//-----------------------------
const int     COMPUTE_MAIN_LOOP = 1;
const int     COMPUTE_VEX       = 1;
const int     COMPUTE_TOPT      = 0;
const int     TOPT_TYPE         = 0; // 0= min time , 1= exit time

const int     SAVE_VF_ALL       = 1;        //- 1 to save "VFn.dat" every  SAVE_VFALL_STEP iterations
const int     SAVE_VF_ALL_STEP  = 1000;    //- (a large number will save only VF0.dat and final: VF1.dat)
const int     SAVE_VF_FINAL     = 1;

const int     CHECK_ERROR       = 1;        //- 1 to compute errors every CHECK_ERROR_STEP iterations
const int     CHECK_ERROR_STEP  = 10;

//-------------------
//- "coupe" : used only for exemples with d>=3 (cut into some plane)
//-------------------
const int     COUPE_DIMS[DIM] = {1 ,1 };
const double  COUPE_VALS[DIM] = {0.,0.};


//---------------
//- initial data
//---------------
inline double   v0(const double * x) {
  double X0=1.0, Y0=0.0, RAYON=0.5;
  double nx=sqrt( (x[0]-X0)*(x[0]-X0) +(x[1]-Y0)*(x[1]-Y0))-RAYON; //- ball target
  //double nx=max( abs(arg[0]-X0) -RAYON ,abs(arg[1]-Y0) -RAYON); //- square target
  return MIN(VBORD,nx);
}


inline void dynamics(const double* x, C u, double t, double* res)
{
  double X=x[0], Y=x[1];
  double f0,f1;
  f0=-(-2*pi*Y);
  f1=-(+2*pi*X);
  //f0=-(-1.0)*2;
  //f1=-(-2.0)*2;
  res[0]=f0, res[1]=f1;
}

inline double distributed_cost(const double* arg, C u, double t)
{
  return 0.;
}


//---------------------------------------------------------------
//- unused functions since only 1 control here (COMMANDS=1)
//---------------------------------------------------------------
inline void dynamics2(const double* x, C u, C u2, double t, double* res) { return; }
inline double distributed_cost2(const double* arg, C u, C u2, double t) { return 0.; }


//---------------------------
//- Exact solution (if known)
//---------------------------
inline double Vex(double t, const double* x)
{
  double y[2];
  double t1=-t;
  y[0]=cos(2*pi*t1)*x[0] - sin(2*pi*t1)*x[1];
  y[1]=sin(2*pi*t1)*x[0] + cos(2*pi*t1)*x[1];
  return v0(y);
}


//-------------------
//- compute trajectory parameters
//-------------------
const int     TRAJPT            = 0;
const double  initialpoint[TRAJPT*DIM] = {};      //- Initial point for HJB
// method of trajectory reconstruction and starting time
int           TRAJ_METHOD       = 0; 
double        time_TRAJ_START   = 0.0;
// stopping criteria
double        min_TRAJ_STOP     = 0.00;    //- to stop traj reconstruction when val(x) <=min (val=topt here)
double        max_TRAJ_STOP     = 100.00;  //- to stop traj reconstruction when val(x) >=max (val=topt here)
int           TARGET_STOP       = 1;       //- to stop traj reconstruction when g_target(x)<=0
inline double g_target(double t, const double* x){ return v0(x);}  
// adverse control case (only for COMMANDS=2 & TRAJ_METHOD=0)
int           ADVERSE_METHOD    = 0;
inline void u2_adverse(double t, const double* x, C& u){}


//---------------------------------------------------
//- These parameters are not used for the SL method -
//---------------------------------------------------
inline void compute_Hconst(double* aMAX, double t){} 
inline double Hnum(const double t, const double* x, const double vi, const double* Dv){return 0.0;}
double  CFL               = 0.8;		//- not used for SL method


//--------------------------------
//- PARAMETERS FOR ADVANCED USERS
//--------------------------------
const int     EXTERNALV0        = 0;           //- if 1 then starts computation from data in VF.dat and not from the v0 function
const int     PRECOMPUTE_COORDS = 1;           //- to precompute mesh coordinates (faster but needs more memory).

//--------------------
//- OBSTACLE g
//--------------------
const int OBSTACLE              = 0;
inline double g_obstacle(double t, const double* arg){ return 0.;}

//--------------------------------
//- border: number of ghost cells in each direction (default is {2,2,...} for FD)
//--------------------------------
int BORDERSIZE[DIM] = {0,0};


