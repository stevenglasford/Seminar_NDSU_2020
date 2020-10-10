//----------------------------------------
//- Basic max(advection,0) example , dimension 2
//----------------------------------------
//-    v_t + max_a < - f(x,a), nabla v >  = 0;    f(x_1,x_2,a)= a (-DYN0,-DYN1);  a in {0,1}
//-    v(0,x)=v0(x)
//----------------------------------------
//- Note: here periodic case and SL approximation scheme.
//-    ==> For SL, better to set MESH=1 (mesh includes boundary values) 
//-        in order to avoid interpolation issues at the boundary
//-        (in particular for trajectory reconstruction)
//-        Then, because of the periodic boundary conditions, the BOUNDARY value has no effect
//-    ==> For SL, if still want to use MESH=0 (mesh at center of cells), better to set BOUNDARY = 3 (Vxx=0 , or linear extention)
//----------------------------------------
#include "data_default.h"                      //- do not modify for basic examples
//----------------------------------------

const char    NAME[150]         = "data_SL_2d_periodic.h, Feb 2017 (TEST ADV/ROT+BOUNDARY TEST - Author : Olivier)";
const int     DIM               = 2;   //- space dimension

const int     COMMANDS          = 0;   //- 0: local Hnum function ; 1/2: Hnum defined with dynamics and distributed_cost functions
const int     OPTIM             = MAXIMUM;



//----------------------------
//- SL method parameters:
//----------------------------
const int     METHOD            = MSL;         //- Method : Semi-Lagrangian
const int     TYPE_SCHEME       = EVO;         //- (SL) : STA (1) or EVO (2)
int           TYPE_RK           = RK2_HEUN;    //- (SL) : type of RK approximation

double        DT                = 0.05;        //- if DT=0.0, MSL will not work
double        CFL               = 0.00;        //- CFL not used for SL method

/*
//----------------------------
//- DF method parameters:
//----------------------------
const int     METHOD            = MFD;         //- Method : Finite differences 
const int     TYPE_SCHEME       = ENO2;        //- space scheme : 1 or LF ; 2 or ENO2
int           TYPE_RK           = RK1;         //- time scheme (ENO2 case) : 1 or ENO2RK1, 2 or ENO2RK2,  3 or ENO2RK3
const double  DT                = 0.00;        //- if DT=0 the program will compute a DT based on space steps.
      double  CFL               = 0.4;
*/

//----------------------------
//- stopping criteria parameters:
//----------------------------
const double  EPSILON           = 0.0;
const int     MAX_ITERATION     = 100000;
double                         T= 3.0;         //- Terminal time


//----------------------------
//- discretization parameters:
//----------------------------
//const double  DT                = 0.00;          //- if DT=0 the program will compute a DT based on space steps.

const int     NN=80;
      int     ND[DIM]           = { NN , NN};
      double  XMIN[DIM]         = { -2. , -2.};
      double  XMAX[DIM]         = {  2. ,  2.};           //- bound of the domain

const int     PERIODIC[DIM]     = {  1  ,  1 };           //- periodic mesh (1:periodic, 0:otherwise)
const int     MESH              = 1;                      //- 0 : xi = center of cell ; 1 : xi contains the boundary

const int     cDIM              = 1;
const int     NCD[cDIM]         = { 2 };
const double  UMIN[cDIM]        = { 0.};
const double  UMAX[cDIM]        = { 1.};

const int     cDIM2             = 1;
const int     NCD2[cDIM2]       = { 1 };
const double  UMIN2[cDIM2]      = { 0.};
const double  UMAX2[cDIM2]      = { 1.};


//-----------------------------
//- mainloop parameters
//-----------------------------
const int     COMPUTE_MAIN_LOOP = 1;
const int     COMPUTE_VEX       = 0;
const int     COMPUTE_TOPT      = 1;
const int     TOPT_TYPE         = 0;       //- 0= min time pb , 1= exit time pb

const int     SAVE_VF_ALL       = 1;       //- 1 to save "VFn.dat" every  SAVE_VF_ALL_STEP iterations
const int     SAVE_VF_ALL_STEP  = 5;       //- (a large number will save only VF0.dat and final: VF1.dat)
const int     SAVE_VF_FINAL     = 1;

const int     CHECK_ERROR       = 1;       //- 1 to compute errors every CHECK_ERROR_STEP iterations
const int     CHECK_ERROR_STEP  = 10;


//-------------------
//- "coupe" : used only for exemples with d>=3 (cut into some plane ==> results in files "coupe.dat"/"coupeex.dat")
//-------------------
const int     COUPE_DIMS[DIM] = {1 ,1 };
const double  COUPE_VALS[DIM] = {0.,0.};


//---------------
//- initial data
//---------------
inline double   v0(const double* x) {
  double VBORD=0.3;
  double X0=1.0, Y0=1.0, RAYON=0.50;
  double y[2]={x[0],x[1]};
  for (int i=0;i<=1;i++) 
    if (PERIODIC[i]==1) y[i]=x[i] - floor((x[i]-XMIN[i])/(XMAX[i]-XMIN[i]))*(XMAX[i]-XMIN[i]);
  //double nx=(y[0]-X0)*(y[0]-X0) +(y[1]-Y0)*(y[1]-Y0)-RAYON*RAYON; 	//- ball   target |x|^2-R^2
  double nx=sqrt( (y[0]-X0)*(y[0]-X0) +(y[1]-Y0)*(y[1]-Y0))-RAYON; 	//- ball   target |x|-R
  //double nx=sqrt( (y[0]-X0)*(y[0]-X0) +(y[1]-Y0)*(y[1]-Y0)+0.1)-sqrt(RAYON*RAYON+0.1); //- ball   target sqrt(|x|^2+e)-sqrt(R^2+e)
  //double nx=max( abs(y[0]-X0) -RAYON ,abs(y[1]-Y0) -RAYON);		//- square target |x|_inf-R
  return min(VBORD,nx);
}

//- dynamics specific to this data file:
const double DYN0=  2.0;
const double DYN1= -1.0;

//-----------------------------
//- BOUNDARY : 0 (Void, FD only) or 1 (Dirichlet, V=g_border, for FD/SL) or 2 (Neumann Vx=g_bordermix, for FD), or 3 (Vxx=0, for FD/SL) 
//-----------------------------
const int BOUNDARY=4;

// Dirichlet boundary condition (case BOUNDARY=1)
double g_border(double t, const double *x){
  double y[2];
  y[0]=x[0]+t*DYN0;
  y[1]=x[1]+t*DYN1;
  for (int i=0;i<=1;i++)
    if (PERIODIC[i]==1) y[i]=x[i] - floor((x[i]-XMIN[i])/(XMAX[i]-XMIN[i]))*(XMAX[i]-XMIN[i]);
  return v0(y);
}

// Mixed Neumann bc :  ux=g(t,x,u) (case BOUNDARY=2)
double g_bordermix(double t, const double* x, double val){return 0.0;}


//------------------------------------------------------------------------------------------------------
//- dynamics and distributed cost functions used in the case of COMMANDS=0 or 1 (assumes only 1 control)
//------------------------------------------------------------------------------------------------------
inline void dynamics(const double* x, C u, double t, double* res)
{
  //double X=x[0], Y=x[1];
  double f0,f1;
  //f0=-(-2*pi*Y);
  //f1=-(+2*pi*X);
  //double eps=0.5;
  //f0=u[1]*(-(-2*pi*Y)*u[0] + eps*(1-u[0]));
  //f1=u[1]*(-(+2*pi*X)*u[0]);// + eps*2*pi*Y*(1-u[0]);
  f0=DYN0;
  f1=DYN1;
  //res[0]=f0, res[1]=f1;
  res[0]=f0*u[0], res[1]=f1*u[0];
}

inline double distributed_cost(const double* arg, C u, double t)
{
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
  double y[2];
  //double t1=-t;
  //y[0]=cos(2*pi*t1)*x[0] - sin(2*pi*t1)*x[1];
  //y[1]=sin(2*pi*t1)*x[0] + cos(2*pi*t1)*x[1];
  y[0]=x[0]+DYN0*t;
  y[1]=x[1]+DYN1*t;
  //- periodization in x and y: in v0
  int i;
  i=0; if (PERIODIC[i]==1) y[i]=y[i] - floor((y[i]-XMIN[i])/(XMAX[i]-XMIN[i]))*(XMAX[i]-XMIN[i]);
  i=1; if (PERIODIC[i]==1) y[i]=y[i] - floor((y[i]-XMIN[i])/(XMAX[i]-XMIN[i]))*(XMAX[i]-XMIN[i]);
  return v0(y);
}


//-------------------
//- compute trajectory parameters
//-------------------
// method of trajectory reconstruction and starting time
int           TRAJ_METHOD       = 0; 
double        time_TRAJ_START   = 0.0;
// initial point 
const int TRAJPT = 1; const double  initialpoint[TRAJPT*DIM] = {-1.0,0.0}; 
//const double  initialpoint[TRAJPT*DIM] = {-2.0+1e-2,-1.0};      //- Initial point for HJB
//const double  initialpoint[TRAJPT*DIM] = {-1.0,1.0};      //- Initial point for HJB
//const double  initialpoint[TRAJPT*DIM] = {-1.0,2.0-5e-2};      //- Initial point for HJB
//const int TRAJPT = 0; const double  initialpoint[TRAJPT*DIM] = {};      //- Initial point for HJB
// stopping criteria
double        min_TRAJ_STOP     = 0.00;    //- to stop traj reconstruction when val(x) <=min (val=topt here)
double        max_TRAJ_STOP     = 100.00;  //- to stop traj reconstruction when val(x) >=max (val=topt here)
int           TARGET_STOP       = 1;       //- to stop traj reconstruction when g_target(x)<=0
inline double g_target(double t, const double* x){ return v0(x);}
// adverse control case (only for COMMANDS=2 & TRAJ_METHOD=0)
int           ADVERSE_METHOD    = 0;
inline void u2_adverse(double t, const double* x, C& u){}


//-----------------------------------------------
//- MORE ADVANCED PARAMETERS
//-----------------------------------------------

//- Case COMMANDS=0: numerical Hamiltonian function Hnum.
//inline void compute_Hconst(double* aMAX, double t){};
//inline double Hnum(const double t, const double* x, const double vi, const double* Dv){return 0.0};

//--------------------
//- Hamiltonian (LF numerical hamiltonian here)
//--------------------
//- case COMMANDS=0; 

inline double H(const double* x, const double* p)
{
  //double X=x[0], Y=x[1];
  double res,f0,f1;
  //f0=-(-2*pi*Y);
  //f1=-(+2*pi*X);
  //f0=(1.0);
  //f1=(0.0);
  f0=DYN0;
  f1=DYN1;
  res=-f0*p[0]-f1*p[1];
  //return res;
  return max(0.,res);
}

//--------------------
//- Stability constants to be used for the numerical hamiltonian.
//--------------------
inline void compute_Hconst(double* aMAX, double t)
{
  aMAX[0]=abs(DYN0);
  aMAX[1]=abs(DYN1);
}

//--------------------
//- Numerical Hamiltonian (Lax Friedriech) in the case H is defined
//--------------------
inline double Hnum(const double t, const double* x, const double vi, const double* Dv)
{
  //- METHOD : Hamiltonian H(.) function + LF stabilisation
  //- for each direction x_i, Dv[2*i]   is the left  centered approximation of dv/dx_i,
  //-                         Dv[2*i+1] is the right centered approximation of dv/dx_i,
  double z=0.;
  double p[DIM],
  aMax[DIM];
  int i;
  for(i=0;i<DIM;i++){
    aMax[i]=1;//2.0*pi*max(abs(XMIN[i]),abs(XMAX[i]));
    p[i]=(Dv[2*i] + Dv[2*i+1])/2.;
    z += aMax[i]*(Dv[2*i+1] - Dv[2*i])/2.;
  }
  return H(x,p)-z;
}

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
int BORDERSIZE[DIM] = {2,2};


//--------------------
//- Special parameter for method SL
//--------------------
const int     P_INTERMEDIATE    = 1;    //- number of discretisation steps for approximating each trajectory in method MSL 

