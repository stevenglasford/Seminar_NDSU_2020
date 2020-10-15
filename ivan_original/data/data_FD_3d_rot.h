//----------------------------------------
//- Advection Equation, dimension 3 (2d advection in x,y + nothing in z - for testing purposes)
//----------------------------------------
//-    u_t + < f, nabla u >  = 0;    f= ( -(-y,x), 0);
//-    u(0,x)=u0(x)
//----------------------------------------
#include "data_default.h"	// do not modify for basic examples
//----------------------------------------

const char    NAME[200]         = "data_FD_3d_rot.h: 3D example (rotation, bound-test, coup-test), Feb/Aug 2016 [OB]";
const int     DIM               = 3;   //- space dimension   

const int     COMMANDS          = 0;   //- 0: user Hnum;  1: Hnum from dynamics()/distributed_cost(); (2: for min/max pbs.)
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
const double                   T= 0.10;        //- Terminal time


//----------------------------
//- discretization parameters:
//----------------------------
const double  DT                = 0.000;       //- if DT=0 the program will compute a DT based on space steps.  
const double  CFL               = 0.5;

const int     NN=50;
const int     ND[DIM]           = { NN , NN , NN };
const double  XMIN[DIM]         = { -2., -2., -2.};
const double  XMAX[DIM]         = {  2.,  2.,  2.};           //- bound of the domain
const int     PERIODIC[DIM]     = {  0 ,  0 ,  0 };           //- periodic mesh (1:periodic, 0:otherwise)

//- for a DIM=4 test:
//const int     ND[DIM]           = { NN , NN , NN ,  3};
//const double  XMIN[DIM]         = { -2., -2., -2., -2};
//const double  XMAX[DIM]         = {  2.,  2.,  2.,  2};       //- bound of the domain
//const int     PERIODIC[DIM]     = {  0 ,  0 ,  0 ,  0};       //- periodic mesh (1:periodic, 0:otherwise)

const int     MESH              = 1;                          //- 0 : xi = center of cell ; 1 : xi contains the boundary

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
const int BOUNDARY=3;

//- VBORD value 
const double  VBORD             = 0.3; 

//- Dirichlet boundary condition:
double g_border(double t, const double* arg){
  return VBORD;
}

//- Mixed Neumann bc :  ux=g(t,x,u) (case BOUNDARY=2)
double g_bordermix(double t, const double* x, double val){return 0.0;}


//-----------------------------
//- mainloop parameters
//-----------------------------
const int     COMPUTE_MAIN_LOOP = 1;
const int     COMPUTE_VEX       = 1;
const int     COMPUTE_TOPT      = 1;
const int     TOPT_TYPE         = 0;       //- 0= min time , 1= exit time


const int     SAVE_VF_ALL       = 0;       //- 1 to save "VFn.dat" every  SAVE_VF_ALL_STEP iterations
const int     SAVE_VF_ALL_STEP  = 1;       //- (a large number will save only VF0.dat and final: VF1.dat)
const int     SAVE_VF_FINAL     = 1;

const int     CHECK_ERROR       = 1;       //- 1 to compute errors every CHECK_ERROR_STEP iterations
const int     CHECK_ERROR_STEP  = 10;


//-------------------
//- "coupe" : useful for exemples with DIM>=3 (cut into some plane)
//-------------------
const int     COUPE_DIMS[DIM] = {  1,      1,     0};
const double  COUPE_VALS[DIM] = {0.12, 0.0234, 0.222};
//const int     COUPE_DIMS[DIM] = {  1,   0,   1};
//const double  COUPE_VALS[DIM] = {0.0, 0.15, 0.0};
//const int     COUPE_DIMS[DIM] = {  1,      1,     0,    0};
//const double  COUPE_VALS[DIM] = {0.0, 0.0234, 0.122,  0.0};

//------------------------------------------------------------------------------------------------------
//- dynamics and distributed cost functions used in the case of COMMANDS=0 or 1 (assumes only 1 control)
//------------------------------------------------------------------------------------------------------
inline void dynamics(const double* x, C u, double t, double* res)
{
  double X=x[0], Y=x[1];
  double f0,f1,f2;
  f0=-(-2*pi*Y);
  f1=-(+2*pi*X);
  f2=       0.0;
  res[0]=f0, res[1]=f1, res[2]=f2;
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
inline double v0(const double * xx)
{

  //- uncomment desired example:
  //- rotating circle:
  double x=xx[0], y=xx[1], z=xx[2];
  double X0=1.0, Y0=0.0, Z0=0.1, RAYON=0.5;
  double nx=sqrt( (x-X0)*(x-X0) + (y-Y0)*(y-Y0) + (z-Z0)*(z-Z0) );
  return min(VBORD,nx-RAYON);

  //- rotating plane:
  //double x=xx[0],y=xx[1],z=xx[2];
  //return z-x+2*y;
}


//---------------------------
//- Exact solution (if known)
//---------------------------
inline double Vex(double t, const double* x) {
  double y[3];
  double t1=-t;
  y[0]=cos(2*pi*t1)*x[0] - sin(2*pi*t1)*x[1];
  y[1]=sin(2*pi*t1)*x[0] + cos(2*pi*t1)*x[1];
  y[2]=x[2];
  return v0(y);
}


//--------------------
//- Stability constants to be used for the numerical hamiltonian.
//--------------------
inline void compute_Hconst(double* aMAX, double t)
{
    double xm=max(ABS(XMIN[0]),abs(XMAX[0]));
    double ym=max(ABS(XMIN[1]),abs(XMAX[1]));
    aMAX[0]=2*pi*ym;
    aMAX[1]=2*pi*xm;
    aMAX[2]=0.;
}

//--------------------
//- Numerical Hamiltonian :
//--------------------
inline double Hnum(const double t, const double* x, const double vi, const double* Dv)
{
  // DIRECT METHOD
  double res;
  double X=x[0], Y=x[1]; 
  double f0=-(-2*pi*Y);
  double f1=-(+2*pi*X);
  //double Z=x[2];
  //double f0=-(+2*pi*Z);
  //double f2=-(+2*pi*X);
  res=  max(0.,-f0)*Dv[0] + min(0.,-f0)*Dv[1] 
      + max(0.,-f1)*Dv[2] + min(0.,-f1)*Dv[3];
      //+ max(0.,-f2)*Dv[4] + min(0.,-f2)*Dv[5];  
  return res;
}


//-------------------
//- compute trajectory parameters
//-------------------
// method of trajectory reconstruction and starting time
int           TRAJ_METHOD       = 0; 
double        time_TRAJ_START   = 0.0;
const int     TRAJPT            = 0;
const double  initialpoint[TRAJPT*DIM] = {};      //- Initial point for HJB

// stopping criteria
double        min_TRAJ_STOP     =   0.00;  //- to stop traj reconstruction when val(x) <=min (val=topt here)
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
//- BOUNDARY : size of boundary in each direction (default is {2,2,...} for FD)
//--------------------------------
int BORDERSIZE[DIM] = {2,2,2};

//--------------------
//- Special parmeter for method SL
//--------------------
const int     P_INTERMEDIATE    = 1;    //- number of discretisation steps for approximating each trajectory in method MSL 

