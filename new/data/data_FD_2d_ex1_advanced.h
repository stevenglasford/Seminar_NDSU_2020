//----------------------------------------
//- Eikonal Equation example, dimension 2
//----------------------------------------
//-    u_t + c(x) || \nabla u || = 0 or  u_t + max_a ( -f(x,a) nabla u ) = 0
//-    u(0,x)=u0(x)
//-    More advanced approach :  using numerical hamiltonian Hnum (COMMANDS=0)
//----------------------------------------
#include "data_default.h"	  // for basic examples
//----------------------------------------

const char    NAME[150]         = "data_FD_2d_ex1_advanced.h (Eikonal equation example) Feb. 2013\nAuthors   : O. Bokanowski, H. Zidani";
const int     DIM               = 2;   //- space dimension

const int     COMMANDS          = 0;   //- 0: local Hnum function ; 1: Hnum defined with dynamics and distributed_cost functions
const int     OPTIM             = MAXIMUM;

//----------------------------
//- DF method parameters:
//----------------------------
const int     METHOD            = MFD;         //- Method : Finite differences(MFD)/Semi-Lagrangian(MSL)
const int     TYPE_SCHEME       = ENO2;        //- space scheme : 1 or LF ; 2 or ENO2 ; 3 or ENO3
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
const double  DT                = 0.0;         //- if DT=0 the program will compute a DT based on space steps.
      double  CFL               = 0.5;

const int     NN=50;
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
const int     TOPT_TYPE         = 0;       //- 0= min time , 1= exit time

const int     SAVE_VF_ALL       = 0;       //- 1 : to save "VFn.dat" every  SAVE_VF_ALL_STEP iterations (useful for movie)
const int     SAVE_VF_ALL_STEP  = 100;     //- (set to a large number to save only VF0.dat and ev. a final: VF1.dat)
const int     SAVE_VF_FINAL     = 1;       //- 1 : to save furthermore the last VF file as "VF.dat" 

const int     CHECK_ERROR       = 1;       //- 1 to compute errors every CHECK_ERROR_STEP iterations (if COMPUTE_VEX=1)
const int     CHECK_ERROR_STEP  = 10;
const int     CHECK_NEG_POINTS  = 0;       //- 1/0 : to print number of negative point in domain and domain + boundary (default : 0)

//-------------------
//- "coupe" : used to make some cut into some plane ==> results in files "coupe.dat"/"coupeex.dat" (if COUPE_DIM not equal to {0,...,0}) 
//-------------------
const int     COUPE_DIMS[DIM] = {1 ,1 }; //- put 1 to keep the dimensions in the final savings.
const double  COUPE_VALS[DIM] = {0.,0.}; //- if COUPE_DIM[d]=0 then make precise the value of the cut in direction x_d.

//---------------
//- initial data
//---------------
const  double   X0=0.0, Y0=0.0, RAYON=0.5;
inline double   q(double x){return max(-0.4,min(VBORD,x-RAYON));}

inline double   v0(const double * arg){
   double nx=sqrt( (arg[0]-X0)*(arg[0]-X0) + (arg[1]-Y0)*(arg[1]-Y0) );
   return q(nx);
}

//------------------------------------------------------------------------------------------------------
//- dynamics and distributed cost functions used in the case of COMMANDS=0 or 1 (assumes only 1 control)
//------------------------------------------------------------------------------------------------------
inline void dynamics(const double* x, C u, double t, double* res)
{
  res[0]=  cos(u[0]);
  res[1]=  sin(u[0]);
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
inline double Vex(double t, const double* x)
{
    double nx= sqrt(x[0]*x[0] + x[1]*x[1]);
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
//- MORE ADVANCED PARAMETERS
//-----------------------------------------------

//--------------------
//- Numerical Hamiltonian and stability constants (case COMMANDS=0);
//--------------------
//- Case COMMANDS=0: numerical Hamiltonian function Hnum.
//inline void compute_Hconst(double* aMAX, double t){};
//inline double Hnum(const double t, const double* x, const double vi, const double* dv){return 0.0;}

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
//- border: number of ghost cells in each direction (default is {2,2,...} for FD)
//--------------------------------
int BORDERSIZE[DIM] = {3,3};

//--------------------
//- Special parameter for method SL
//--------------------
const int     P_INTERMEDIATE    = 1;    //- number of discretisation steps for approximating each trajectory in method MSL

// -------------------
// - For parameter intialization/termination
// -------------------
void post_data(){}
void init_data(){
  SAVE_COUPE_ALL_STEP= 10; //- can be different from SAVE_VF_ALL_STEP; used for all "coupe" files
  SAVE_COUPE_ALL     = 1;
  SAVE_COUPE_FINAL   = 1;
  SAVE_COUPEEX_ALL   = 1;
  SAVE_COUPEEX_FINAL = 1;
}

