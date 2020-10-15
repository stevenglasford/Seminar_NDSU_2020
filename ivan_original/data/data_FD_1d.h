//----------------------------------------
//- 1D example 
//----------------------------------------
//-    v_t + max_a / min_a < - f(x,a), nabla v >  = 0;    f(x_1,a)= a;  a in {0,1}
//-    v(0,x)=v0(x)
//-    "max_a"  or  "min_a"  1d tests problems
//----------------------------------------
#include "data_default.h"	// do not modify for basic examples
//----------------------------------------

const char    NAME[150]         = "data_FD_1d.h, April. 2013, for 1d test\nAuthors : O. Bokanowski, H. Zidani";
const int     DIM               = 1;   //- space dimension

const int     COMMANDS          = 1;   //- 0: local Hnum function ; 1/2: Hnum defined with dynamics and distributed_cost functions
const int     OPTIM             = MAXIMUM;
//const int     OPTIM             = MINIMUM;


//----------------------------
//- DF method parameters:
//----------------------------
const int     METHOD            = MFD;         //- Method : Finite differences 
const int     TYPE_SCHEME       = ENO2;        //- space scheme : 1 or LF ; 2 or ENO2
int           TYPE_RK           = RK2;         //- time scheme (ENO2 case) : 1 or ENO2RK1, 2 or ENO2RK2,  3 or ENO2RK3

//----------------------------
//- stopping criteria parameters:
//----------------------------
const double  EPSILON           = 0.0;
const int     MAX_ITERATION     = 100000;
double                         T= 1.0;         //- Terminal time


//----------------------------
//- discretization parameters:
//----------------------------
const double  DT                = 0.;          //- if DT=0 the program will compute a DT based on space steps.
      double  CFL               = 0.5;

const int     NN=80;
      int     ND[DIM]           = { NN };
      double  XMIN[DIM]         = { -2.};
      double  XMAX[DIM]         = {  2.};                 //- bound of the domain

const int     PERIODIC[DIM]     = {  0 };                 //- periodic mesh (1:periodic, 0:otherwise)
const int     MESH              = 1;                      //- 0 : xi = center of cell ; 1 : xi contains the boundary

const int     cDIM              = 1;
const int     NCD[cDIM]         = { 2 };
const double  UMIN[cDIM]        = { -1.0 };
const double  UMAX[cDIM]        = {  1.0 };

const int     cDIM2             = 1;
const int     NCD2[cDIM2]       = { 1 };
const double  UMIN2[cDIM2]      = { 0.};
const double  UMAX2[cDIM2]      = { 1.};


//-----------------------------
//- BOUNDARY : 0 (Void, FD only) or 1 (Dirichlet, V=g_border, for FD/SL) or 2 (Vx=g_bordermix, for FD), or 3 (Vxx=0, for FD/SL) 
//-----------------------------
const int BOUNDARY=2;

// Dirichlet boundary condition (case BOUNDARY=1)
const double  VBORD             = 1.0; 
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

const int     SAVE_VF_ALL       = 1;       //- 1 to save "VFn.dat" every  SAVE_VFALL_STEP iterations
const int     SAVE_VF_ALL_STEP  = 5;       //- (a large number will save only VF0.dat and final: VF1.dat)
const int     SAVE_VF_FINAL     = 1;

const int     CHECK_ERROR       = 1;       //- 1 to compute errors every CHECK_ERROR_STEP iterations
const int     CHECK_ERROR_STEP  = 10;


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
inline double g_target(double t, const double* x){ return 0.0;}  
// adverse control case (only for COMMANDS=2 & TRAJ_METHOD=0)
int           ADVERSE_METHOD    = 0;
void u2_adverse(double t, const double* x, C& u){}


//-------------------
//- "coupe" : used only for exemples with d>=3 (cut into some plane ==> results in files "coupe.dat"/"coupeex.dat")
//-------------------
const int     COUPE_DIMS[DIM] = {1 };
const double  COUPE_VALS[DIM] = {0.};

//---------------
//- initial data
//---------------

inline double   v0(const double * x) {
  double X0=0.0;
  //double nx=abs(x[0]-X0);
  double nx=sqrt((x[0]-X0)*(x[0]-X0));
  return min(nx*nx,2*VBORD);
}

//------------------------------------------------------------------------------------------------------
//- dynamics and distributed cost functions used in the case of COMMANDS=0 or 1 (assumes only 1 control)
//------------------------------------------------------------------------------------------------------
inline void dynamics(const double* x, C u, double t, double* res){
  //double  X=x[0];
  double f0;
  f0=-u[0];
  res[0]=f0;
}

inline double distributed_cost(const double* arg, C u, double t){
  return 0.;
}


//---------------------------------------------------------------
//- if COMMANDS=2 (2 player games): dynamics and distributed cost
//---------------------------------------------------------------
inline void dynamics2(const double* x, C u, C u2, double t, double* res) { return; }
inline double distributed_cost2(const double* arg, C u, C u2, double t) { return 0.; }


//---------------------------
//- Exact solution (if known)
//---------------------------
inline double Vex(double t, const double* x)
{
  double y[1]={x[0]};
  double res=0.0;
  //double nx=max(0.,abs(x)-t);
  double res1,res2;
  //if (OPTIM==MAXIMUM){
    double xx=x[0];
    if (xx> t) {res=(xx-t)*(xx-t);}
    if (xx<-t) {res=(xx+t)*(xx+t);}
  //}
  if (OPTIM==MINIMUM){
     y[0] = x[0]-t;
     res1=v0(y);
     y[0] = x[0]+t;
     res2=v0(y);
     res=max(res1,res2);
  }
  //else{
  //  printf("This OPTIM=%i not programmed. Abort.\n",OPTIM);
  //}
  return res;
}

//-----------------------------------------------
//- MORE ADVANCED PARAMETERS
//-----------------------------------------------

//- Case COMMANDS=0: numerical Hamiltonian function Hnum.
inline void compute_Hconst(double* aMAX, double t){};
inline double Hnum(const double t, const double* x, const double vi, const double* v){return 0.0;}

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
int BORDERSIZE[DIM] = {2};

//--------------------
//- Special parameter for method SL
//--------------------
const int     P_INTERMEDIATE    = 1;    //- number of discretisation steps for approximating each trajectory in method MSL 

