//----------------------------------------
//- Rotation example , dimension 2
//----------------------------------------
//-    v_t + max_a < - f(x,a), nabla v >  = 0;    f(x_1,x_2,a)= a (-x_2,x_1);  a in {0,1}
//-    v(0,x)=v0(x)
//----------------------------------------
#include "data_default.h"	// do not modify for basic examples
//----------------------------------------

const char    NAME[]            = "data_advancedmodel.h, Jun. 2014 (O. Bokanowski, H. Zidani)";
const int     DIM               = 2;   //- space dimension

const int     COMMANDS          = 0;   //- 0: local Hnum function ; 1/2: Hnum defined with dynamics and distributed_cost functions
const int     OPTIM             = MAXIMUM;

//----------------------------
//- DF method parameters:
//----------------------------
const int     METHOD            = MFD;         //- Method : Finite differences(MFD)/Semi-Lagrangian(MSL)
const int     TYPE_SCHEME       = ENO2;        //- space scheme : 1 or LF ; 2 or ENO2
const int     TYPE_RK           = RK1;         //- time scheme (ENO2 case) : 1 or ENO2RK1, 2 or ENO2RK2,  3 or ENO2RK3

//----------------------------
//- stopping criteria parameters:
//----------------------------
const double  EPSILON           = 0.0;
const int     MAX_ITERATION     = 100000;

double                         T= 0.5;  //0.5;         //- Terminal time

//----------------------------
//- discretization parameters:
//----------------------------
const double  DT                = 0.0;         //- if DT=0 the program will compute a DT based on space steps.
      double  CFL               = 0.5;

const int     NN=80;
      int     ND[DIM]           = { NN , NN};
      double  XMIN[DIM]         = { -2. , -2.};
      double  XMAX[DIM]         = {  2. ,  2.};           //- bound of the domain

const int     PERIODIC[DIM]     = {  0  ,  0 };           //- periodic mesh (1:periodic, 0:otherwise)
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
//- BOUNDARY : 0 (Void, FD only) or 1 (Dirichlet, V=g_border, for FD/SL) or 2 (Vx=g_bordermix, for FD), or 3 (Vxx=0, for FD/SL) 
//-----------------------------
const int BOUNDARY=1;

const double  VBORD             = 0.2;

//- For Dirichlet boundary condition: u=g(t,x) (case BOUNDARY=1)
double g_border(double t, const double* x){
  return VBORD;
}

//- For mixed Neumann b.c. :  ux=g(t,x,u) (case BOUNDARY=2)
double g_bordermix(double t, const double* x, double val){return 0.0;}

//-----------------------------
//- mainloop parameters
//-----------------------------
const int     COMPUTE_MAIN_LOOP = 1;
const int     COMPUTE_VEX       = 0;
const int     COMPUTE_TOPT      = 1;
const int     TOPT_TYPE         = 0;       //- 0= min time , 1= exit time

const int     SAVE_VF_ALL       = 1;       //- 1 to save "VFn.dat" every  SAVE_VF_ALL_STEP iterations (useful for movie)
const int     SAVE_VF_ALL_STEP  = 100;     //- (a large number will save only VF0.dat and final: VF1.dat)
const int     SAVE_VF_FINAL     = 1;       //- 1 to save last VF file as "VF.dat" 

const int     CHECK_ERROR       = 1;       //- 1 to compute errors every CHECK_ERROR_STEP iterations
const int     CHECK_ERROR_STEP  = 100;

//-------------------
//- "coupe" : used to make some cut into some plane ==> results in files "coupe.dat"/"coupeex.dat" (if COUPE_DIM not equal to {0,...,0}) 
//-------------------
const int     COUPE_DIMS[DIM] = {1 ,1 }; //- put 1 to keep the dimensions in the final savings.
const double  COUPE_VALS[DIM] = {0.,0.}; //- if COUPE_DIM[d]=0 then make precise the value of the cut in direction x_d.

//---------------
//- initial data
//---------------

inline double   v0(const double* x) {
  double X0=1.0, Y0=0.0, RAYON=0.25;
  //double nx=sqrt( (x[0]-X0)*(x[0]-X0) +(x[1]-Y0)*(x[1]-Y0))-RAYON; 	//- ball target
  //double nx=max( abs(x[0]-X0) -RAYON ,abs(x[1]-Y0) -RAYON);		//- square target
  //return min(VBORD,nx);
  double r02=RAYON*RAYON, nx2=(x[0]-X0)*(x[0]-X0) +(x[1]-Y0)*(x[1]-Y0); //- smoother ball target
  return VBORD*(1-pow(max(0.0, (1-nx2)/(1-r02)),4));
}


//------------------------------------------------------------------------------------------------------
//- dynamics and distributed cost functions used in the case of COMMANDS=0 or 1 (assumes only 1 control)
//------------------------------------------------------------------------------------------------------
inline void dynamics(const double* x, C u, double t, double* res){
  double  X=x[0], Y=x[1];
  double f0,f1;
  f0=-(-2*pi*Y)*u[0];
  f1=-(+2*pi*X)*u[0];
  res[0]=f0, res[1]=f1;
}

inline double distributed_cost(const double* x, C u, double t){
  return 0.;
}


//---------------------------------------------------------------
//- if COMMANDS=2 (2 player games): dynamics and distributed cost
//---------------------------------------------------------------
inline void dynamics2(const double* x, C u, C u2, double t, double* res) { return; }
inline double distributed_cost2(const double* x, C u, C u2, double t) { return 0.; }


//---------------------------
//- Exact solution (if known)
//---------------------------
inline double Vex(double t, const double* x)
{
  double x2[2];
  double t1=-t;
  x2[0]=cos(2*pi*t1)*x[0] - sin(2*pi*t1)*x[1];
  x2[1]=sin(2*pi*t1)*x[0] + cos(2*pi*t1)*x[1];
  return v0(x2);
}

//-----------------------------------------------
//- PARAMETERS : compute trajectory 
//-----------------------------------------------
// method of trajectory reconstruction and starting time
int           TRAJ_METHOD       = 0; 
double        time_TRAJ_START   = 0.00;
//const int     TRAJPT     = 0;  const double  initialpoint[TRAJPT*DIM] = {}; //- Initial point for HJB
const int     TRAJPT     = 1;  const double  initialpoint[TRAJPT*DIM] = {-1.,0.0};

// stopping criteria
double        min_TRAJ_STOP     = 0.00;    //- to stop traj reconstruction when val(x) <=min (val=topt here)
double        max_TRAJ_STOP     = 100.00;  //- to stop traj reconstruction when val(x) >=max (val=topt here)
int           TARGET_STOP       = 1;       //- to stop traj reconstruction when g_target(x)<=0
inline double g_target(double t, const double* x){ return v0(x);}  

// adverse control case (only for COMMANDS=2 & TRAJ_METHOD=0): if ADVERSE_METHOD=1, allows to define some u2 adverse control
int           ADVERSE_METHOD    = 0;
void u2_adverse(double t, const double* x, C& u){}


//-----------------------------------------------
//- MORE ADVANCED PARAMETERS: user-defined numerical Hamiltonian
//-----------------------------------------------

//--------------------
//- Stability constants to be used for the numerical hamiltonian (bounds for dH/dp_i, H=H(x,p))
//--------------------
inline void compute_Hconst(double* aMAX, double t)
{
  aMAX[0]=2.0*pi*max(abs(XMIN[1]),abs(XMAX[1]));
  aMAX[1]=2.0*pi*max(abs(XMIN[0]),abs(XMAX[0]));
}

/*
//--------------------
//- Hamiltonian (LF numerical hamiltonian here)
//--------------------
//- case COMMANDS=0; 

inline double H(const double* x, const double* p)
{
  double X=x[0], Y=x[1];
  double res,f0,f1;
  f0=-(-2*pi*Y);
  f1=-(+2*pi*X);
  res=-f0*p[0]-f1*p[1];
  return max(0.,res);
  //return res;
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
    aMax[i]=2.0*pi*max(abs(XMIN[i]),abs(XMAX[i]));
    p[i]=(Dv[2*i] + Dv[2*i+1])/2.;
    z += aMax[i]*(Dv[2*i+1] - Dv[2*i])/2.;
  }
  return H(x,p)-z;
}
*/


//--------------------
//- Numerical Hamiltonian : direct approximation (better than Lax-Driedriech)
//--------------------
inline double Hnum(const double t, const double* x, const double vi, const double* Dv)
{
  double X=x[0], Y=x[1];
  double f0,f1,res;
  f0=-(-2*pi*Y);
  f1=-(+2*pi*X);
  res=  max(-f0,0.0)*Dv[0]+min(-f0,0.0)*Dv[1]
      + max(-f1,0.0)*Dv[2]+min(-f1,0.0)*Dv[3];
  return max(0.0,res);
}


//--------------------------------
//- PARAMETERS FOR ADVANCED USERS
//--------------------------------
const int     EXTERNALV0        = 0;           //- if 1 then starts computation from data in VF.dat and not from the v0 function
const int     PRECOMPUTE_COORDS = 1;           //- to precompute mesh coordinates (faster but needs more memory).

//--------------------
//- OBSTACLE g
//--------------------
const int OBSTACLE=0;
inline double g_obstacle(double t, const double* x){
  double X0=0.0, Y0=0.75, RAYON=0.25;
  double nx=max( abs(x[0]-X0) -RAYON ,abs(x[1]-Y0) -RAYON);		//- small square obstacle
  return MIN(VBORD,-nx);
}

//--------------------------------
//- BORDERSIZE : size of boundary in each direction (default is {2,2,...} for FD)
//--------------------------------
int BORDERSIZE[DIM] = {2,2};

//--------------------
//- Special parameter for method SL
//--------------------
const int     P_INTERMEDIATE    = 1;    //- number of discretisation steps for approximating each trajectory in method MSL

