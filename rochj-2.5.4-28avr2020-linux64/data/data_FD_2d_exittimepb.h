//----------------------------------------
//- Exit time pb example , dimension 2
//---------------------------------------- //-  We look for viable points, such that y^u(.)\subset K
//-  dynamics is f(x,u)= (1-u)*(1,0) + u* (0,-1)  for u in [0,1]
//-  Value pb is v(t,x)= inf_u max_s g(y^u_x(s))
//-  PDE solved is 
//-    min( v_t + max_u <  - f(x,u), nabla v >,  v- g(x))  = 0;    
//-    v(0,x)=g(x) 
//-  Tmax is computed  (==> we set TOPT_TYPE=1)
//----------------------------------------
#include "data_default.h"	// do not modify for basic examples
//----------------------------------------

const char    NAME[]            = "data_FD_2d_exittimepb.h, Jul. 2015 (O. Bokanowski, A Desilles)";
const int     DIM               = 2;   //- space dimension

const int     COMMANDS          = 1;   //- 0: local Hnum function ; 1/2: Hnum defined with dynamics and distributed_cost functions
const int     OPTIM             = MAXIMUM;

//----------------------------
//- DF method parameters:
//----------------------------
const int     METHOD            = MFD;         //- Method : Finite differences
const int     TYPE_SCHEME       = ENO2;        //- space scheme : 1 or LF ; 2 or ENO2
const int     TYPE_RK           = RK1;         //- time scheme (ENO2 case) : 1 or ENO2RK1, 2 or ENO2RK2,  3 or ENO2RK3

//----------------------------
//- stopping criteria parameters:
//----------------------------
const double  EPSILON           = 0.0;
const int     MAX_ITERATION     = 100000;

      double                   T= 1.2;         //- Terminal time

//----------------------------
//- discretization parameters:
//----------------------------
const double  DT                = 0.0;         //- if DT=0 the program will compute a DT based on space steps.
      double  CFL               = 0.5;

const int     NN=80;
      int     ND[DIM]           = { NN , NN};
      double  XMIN[DIM]         = { -2.5 , -2.5};
      double  XMAX[DIM]         = {  2.5 ,  2.5};         //- bound of the domain

const int     PERIODIC[DIM]     = {  0  ,  0 };           //- periodic mesh (1:periodic, 0:otherwise)
const int     MESH              = 1;                      //- 0 : xi = center of cell ; 1 : xi contains the boundary

/*
const int     cDIM              = 1;
const int     NCD[cDIM]         = { 3 };
const double  UMIN[cDIM]        = { 0.};
const double  UMAX[cDIM]        = { 1.};
*/
const int     cDIM              = 2;
const int     NCD[cDIM]         = { 3  ,   10  };
const double  UMIN[cDIM]        = { 0.0,  0.0  };
const double  UMAX[cDIM]        = { 0.4, 2*pi  };


const int     cDIM2             = 1;
const int     NCD2[cDIM2]       = { 1 };
const double  UMIN2[cDIM2]      = { 0.};
const double  UMAX2[cDIM2]      = { 1.};

//-----------------------------
//- BOUNDARY : 0 (Void, FD only) or 1 (Dirichlet, V=g_border, for FD/SL) or 2 (Vx=g_bordermix, for FD), or 3 (Vxx=0, for FD/SL) 
//-----------------------------
const int BOUNDARY=1;

const double  VBORD             = 0.5;

//-  Dirichlet boundary condition:
double g_border(double t, const double* x){
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
const int     TOPT_TYPE         = 1;       // 0= min time , 1= max time (=exit time pb)

const int     SAVE_VF_ALL       = 1;       //- 1 to save "VFn.dat" every  SAVE_VF_ALL_STEP iterations (useful for movie)
const int     SAVE_VF_ALL_STEP  = 2;       //- (a large number will save only VF0.dat and final: VF1.dat)
const int     SAVE_VF_FINAL     = 1;

const int     CHECK_ERROR       = 1;       //- 1 to compute errors every CHECK_ERROR_STEP iterations
const int     CHECK_ERROR_STEP  = 10;
const int     CHECK_NEG_POINTS  = 0;       //- 1/0 : to print number of negative point in domain and domain + boundary (default : 0)

//-------------------
//- "coupe" : used to make some cut into some plane ==> results in files "coupe.dat"/"coupeex.dat" (if COUPE_DIM not equal to {0,...,0}) 
//-------------------
const int     COUPE_DIMS[DIM] = {1 ,1 };
const double  COUPE_VALS[DIM] = {0.,0.};

//---------------
//- initial data
//---------------

inline double   v0(const double * x) {
  double X0=0.0, Y0=0.0, R0=2.0;
  double nx0=max( abs(x[0]-X0) -R0,abs(x[1]-Y0) -R0);		//- square target Omega0:=[-2,2]^2: nx<0 inside Omega0

  double X1=-0.5, Y1=-0.5, R1=0.5;
  double nx1=max( abs(x[0]-X1) -R1,abs(x[1]-Y1) -R1);		//- square obstacle Omega1:=[-1,0]^2: -nx1<0 OUTside Omega1

  return MIN(VBORD, max(nx0,-nx1));
}

//------------------------------------------------------------------------------------------------------
//- dynamics and distributed cost functions used in the case of COMMANDS=0 or 1 (assumes only 1 control)
//------------------------------------------------------------------------------------------------------
inline void dynamics(const double* x, C u, double t, double* res){
  //double  X=x[0], Y=x[1];
  //res[0]=(1-u[0])*(-1.0);
  //res[1]=   u[0] *(-1.0);
  /*
  res[0]=( -1.0);
  res[1]=( -1.0);
  */
  res[0]=-1.0 +  u[0]*cos(u[1]);
  res[1]=-1.0 +  u[0]*sin(u[1]);
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
  return v0(x);
}


//-------------------
//- compute trajectory parameters
//-------------------
// method of trajectory reconstruction and starting time
int           TRAJ_METHOD       = 0;
double        time_TRAJ_START   = 0.00;

//const int     TRAJPT     = 0;  const double  initialpoint[TRAJPT*DIM] = {}; //- Initial point for HJB
//const int     TRAJPT     = 1;  const double  initialpoint[TRAJPT*DIM] = {0.0,1.0}; //- Initial point for HJB
//const int     TRAJPT     = 1;  const double  initialpoint[TRAJPT*DIM] = {-0.5,1.0};
//const int     TRAJPT     = 3;  const double  initialpoint[TRAJPT*DIM] = {-0.5,1.0,  -1.0,-1.5,  1.0, 0.5};
const int     TRAJPT     = 4;  const double  initialpoint[TRAJPT*DIM] = {-0.5,1.0,  -1.0,-1.5,  1.0, 0.5,  0.4, 0.6};

// stopping criteria
double        min_TRAJ_STOP     = -10.00;   //- to stop traj reconstruction when val(x) <=min_TRAJ_STOP
double        max_TRAJ_STOP     =  10.00;   //- to stop traj reconstruction when val(x) >=max_TRAJ_STOP
int           TARGET_STOP       =      0;   //- to stop traj reconstruction when g_target(x)<=0
inline double g_target(double t, const double* x){ return -v0(x);}   //! CAREFULL : MINUS SIGN IN FRONT OF v0(x)
                                                                     //! since v0(x)<=0 is the viability domain here

// adverse control case (only for COMMANDS=2 & TRAJ_METHOD=0)
int           ADVERSE_METHOD    = 0;
inline void u2_adverse(double t, const double* x, C& u){}


//-----------------------------------------------
//- MORE ADVANCED PARAMETERS
//-----------------------------------------------
//- Case COMMANDS=0: numerical Hamiltonian function Hnum.
inline void compute_Hconst(double* aMAX, double t){};
inline double Hnum(const double t, const double* x, const double vi, const double* dv){ return 0.; };

//--------------------
//- OBSTACLE g
//--------------------
const int OBSTACLE=1;
inline double g_obstacle(double t, const double* x){
  //double X0=0.0, Y0=0.75, RAYON=0.25;
  //double nx=max( abs(x[0]-X0) -RAYON ,abs(x[1]-Y0) -RAYON);		//- small square obstacle
  //return MIN(VBORD,-nx);
  return v0(x);
}

//--------------------------------
//- BORDERSIZE : size of boundary in each direction (default is {2,2,...} for FD)
//--------------------------------
int BORDERSIZE[DIM] = {3,3};

//--------------------
//- Special parameter for method SL
//--------------------
const int     P_INTERMEDIATE    = 1;    //- number of discretisation steps for approximating each trajectory in method MSL

void post_data(){}
void init_data(){}

