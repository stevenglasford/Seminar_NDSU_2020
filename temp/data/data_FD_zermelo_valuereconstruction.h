//----------------------------------------
//- Zermelo Example (WITH RECONSTRUCTION USING THE VALUE)
//----------------------------------------
//-    u_t + max_a f(x,a) . nabla u(x) >  = 0;    f= zermelo dynamics
//-    u(0,x)=u0(x)
//----------------------------------------
//#include "data_default.h"	//- do not modify for basic examples
#include "data_FD_zermelo_default.h"	//- do not modify for basic examples
//----------------------------------------

const char    NAME[200]         = "data_FD_zermelo_valuereconstruction.h;  Date: Apr. 2014\nO. Bokanowski, A. Desilles, H. Zidani";
const int     DIM               = 2;   //- space dimension

const int     COMMANDS          = 0;   //- 0: user Hnum;  1: Hnum from dynamics()/distributed_cost(); (2: for min/max pbs.)
const int     OPTIM             = MAXIMUM;


//----------------------------
//- DF method parameters:
//----------------------------
const int     METHOD            = MFD;         //- Method : Finite differences(MFD)/Semi-Lagrangian(MSL)
const int     TYPE_SCHEME       = ENO2;        //- scheme : 1 or LF ; 2 or ENO2
const int     TYPE_RK           = RK1;         //- time scheme (ENO2 case) : 1 or ENO2RK1, 2 or ENO2RK2,  3 or ENO2RK3

//----------------------------
//- stopping criteria parameters:
//----------------------------
const double  EPSILON           = 0.0;
const int     MAX_ITERATION     = 100000;

      double                   T= 4.0;         //- Terminal time

//----------------------------
//- discretization parameters:
//----------------------------
const double  DT                = 0.00;        //- if DT=0 the program will compute a DT based on space steps.  
const double  CFL               = 0.5;

int     NN=80;
int     ND[DIM]           = {  NN, NN};
double  XMIN[DIM]         = { -5.0, -2.2 };
double  XMAX[DIM]         = {  2.0,  2.2 };           //- bound of the domain

const int     PERIODIC[DIM]     = {  0 ,  0   };           //- periodic mesh (1:periodic, 0:otherwise)
const int     MESH              = 1;                       //- 0 : xi = center of cell ; 1 : xi contains the boundary
const int     OBSTACLE          = 1;


//- controls:
const int     cDIM         = 2;
const int     NC1          = 1;     //- number of controls 1
const int     NC2          = 10;    //- number of controls 2
const int     NC           = NC1*NC2;
      int     NCD[cDIM]    = { NC1 , NC2};
      double  UMIN[cDIM]   = {1.0,  0.};
      double  UMAX[cDIM]   = {1.0,  2*pi};

const int     cDIM2             = 1;
const int     NCD2[cDIM2]       = { 1 };
const double  UMIN2[cDIM2]      = { 0.};
const double  UMAX2[cDIM2]      = { 1.};

//-----------------------------
//- BOUNDARY : 0 (Void, FD only) or 1 (Dirichlet, V=g_border, for FD/SL) or 2 (Vx=g_bordermix, for FD), or 3 (Vxx=0, for FD/SL) 
//-----------------------------
const int BOUNDARY=1;

//- Dirichlet boundary condition : u=g(t,x) (case BOUNDARY=1) 
const double  VBORD             = 0.2;

double g_border(double t, const double* x){
  return VBORD;
}

//- Mixed Neumann bc :  vx=g(t,x,v) (case BOUNDARY=2)
double g_bordermix(double t, const double* x, double val){return 0.0;}


//-----------------------------
//- mainloop parameters
//-----------------------------
      int     COMPUTE_MAIN_LOOP = 1;
const int     COMPUTE_VEX       = 0;
const int     COMPUTE_TOPT      = 0;
const int     TOPT_TYPE         = 0;         //- 0= min time , 1= exit time

const int     SAVE_VF_ALL       = 1;         //- 1 to save "VFn.dat" every  SAVE_VF_ALL_STEP iterations
const int     SAVE_VF_ALL_STEP  = 5;         //-   (a large number will save only VF0.dat and final: VF1.dat)
const int     SAVE_VF_FINAL     = 1;         //- 1 to save last VF file as "VF.dat" 

const int     CHECK_ERROR       = 1;         //- 1 to compute errors every CHECK_ERROR_STEP iterations
const int     CHECK_ERROR_STEP  = 100;
const int     CHECK_NEG_POINTS  = 0;         //- 1/0 : to print number of negative point in domain and domain + boundary (default : 0)


//-------------------
//- "coupe" : in general used for exemples with dimension d>=3 (cut into some plane)
//-------------------
const int     COUPE_DIMS[DIM] = {1 , 1  };
const double  COUPE_VALS[DIM] = {0., 0. };


//- Modify this for obstacles.
int nbobs = 1; double obstacle[1*4]={-1.7, 1.0, 0.5, 0.5};


double c=2.0, a=0.5;
inline void dynamics(const double* x, C u, double t, double* res)
{
  res[0]= u[0]*cos(u[1]) + c-a*x[1]*x[1];  //- x
  res[1]= u[0]*sin(u[1]);                  //- y
}


inline double distributed_cost(const double* x, C u, double t)
{
  double res=0.0;
  //- Here modification of the distributed_cost that can be used for reconstruction by tmin with a slight cost in the dynamics
  /* 
  double Fnorm=0.0;
  double y[DIM];
  dynamics(x,u,t,y);
  for(int d=0; d<DIM; d++)
    Fnorm=max(Fnorm, abs(y[d]));
  res=0.1 * Fnorm;
  */
  return res;
}

inline void feedback(double t, const double* x, const double* p, C& u){  //- only needed for NBEE-OC
  //- input  : t,x,p
  //- output : u, that maximizes  max_u (- p. dynamics(t,x,u) - dist_cost(t,x,u)) 
  //- Example
  u[0]=0.0;
}


//-----------------------
//- fonction obstacle
//-----------------------
double ymin=-2.0, ymax=2.0, deltay=VBORD;
double xmin=-4.5, xmax=1.5;

inline double g_water(const double* X){
  // fonction obstacle pour la contrainte y dans [ymin ymax]:
  //  g(y)<= dans [ymin ymax]
  //return  ABS(y-(ymin+ymax)/2.) - (ymax-ymin)/2.
  //double x=X[0];
  double y=X[1];
  double res=max(y-ymax,ymin-y);
  //res=max(res,xmin-x);
  //res=max(res,-VBORD); //cutoff part - not necessary
  return min(res,VBORD);
}

inline double g_obstacle(double t, const double* arg)
{
  double x1=arg[0], x2=arg[1];
  double res=-INF, current;
  for(int i=0;i<nbobs;i++){
    current = -MAX(ABS(x1-(obstacle[i*4]))-obstacle[i*4+2], ABS(x2-(obstacle[i*4+1]))-obstacle[i*4+3]);
    res=MAX(res,current);
  }
  res= MAX(res,g_water(arg));
  return MIN(res, VBORD);
}

//- unused functions since only 1 control here (COMMANDS=1)
inline void dynamics2(const double* x, C u, C u2, double t, double* res) { return; }
inline double distributed_cost2(const double* arg, C u, C u2, double t) { return 0.; }


//---------------
//- initial data
//---------------
double X0=0.0, Y0=0.0, RAYONx=1.0,RAYONy=1.0, RAYON=0.5;
int targetShape=0; // 0<=> disc, 1 <=> rectangle

inline double   v0(const double * arg)
{
  double nx=0.0;
  switch(targetShape)
  {
    case 0:
      nx=sqrt( (arg[0]-X0)*(arg[0]-X0) +(arg[1]-Y0)*(arg[1]-Y0))-RAYON;
      break;

    case 1:
      nx=max( abs(arg[0]-X0) -RAYONx ,abs(arg[1]-Y0) -RAYONy);
      break;
  }
  return MIN(VBORD,nx);
}


//---------------------------
//- Exact solution (if known)
//---------------------------
inline double Vex(double t, const double* x){
  return v0(x);
}

//--------------------
//- Stability constants to be used for the numerical hamiltonian.
//--------------------
inline void compute_Hconst(double* aMAX, double t)
{
  double tempMax=max(abs(c-a*XMAX[1]*XMAX[1]),abs(c-a*XMIN[1]*XMIN[1]));
  aMAX[0]=UMAX[0]+max(1.0,tempMax);
  aMAX[1]=UMAX[0];
}


//--------------------
//- Numerical Hamiltonian (Lax Friedriech) in the case H is defined
//--------------------
inline double Hnum(const double t, const double* x, const double vi, const double* Dv)
{
  double cc[DIM];
  double p[DIM];
  double z=0;
  double umax=max(UMAX[0], UMIN[0]);
  compute_Hconst(cc,t);
  int i;
  for(i=0;i<DIM;i++) {
    p[i]=(Dv[2*i] + Dv[2*i+1])/2.0;
    z +=cc[i]* (Dv[2*i+1] - Dv[2*i])/2.0;
  }
  double res=max(0.0,-(c-a*x[1]*x[1]))*Dv[0]+min(0.0,-(c-a*x[1]*x[1]))*Dv[1]+umax*(sqrt(p[0]*p[0]+p[1]*p[1])-z) ;
  return res;
  //return MAX(0.0,res);
}

//-----------------------------------------------
//- PARAMETERS : compute trajectory 
//-----------------------------------------------
//- number of initial points and list of coordinates.
//int TRAJPT = 0;     double initialpoint[0]={};     //- Initial point for HJB
//int TRAJPT = 1;     double initialpoint[2]={-4.0, 2.0};     //- Initial point for HJB
//int TRAJPT = 1;     double initialpoint[2]={-3.0, 2.0};     //- Initial point for HJB
//int TRAJPT = 2;     double initialpoint[2*DIM]={-3.0, 1.8, 0.0, 1.5};     //- Initial point for HJB
//int TRAJPT = 3;     double initialpoint[3*DIM]={-3.0, 1.8, 0.0, 1.5, -0.5, 1.25};     //- Initial point for HJB
//int TRAJPT = 1;     double initialpoint[1*DIM]={0.0, 1.5};     //- Initial point for HJB
int TRAJPT = 3;     double initialpoint[3*DIM]={-3.0, 1.8, 0.0, 1.5, -3.0, -1.50};     //- Initial point for HJB
//int TRAJPT = 1;     double initialpoint[1*DIM]={-3.0, 1.0};     //- Initial point for HJB

// method of trajectory reconstruction and starting time
//-   TRAJ_METHOD=0 <==> based on tmin, used for minimal time as well as exit time problems
//-   TRAJ_METHOD=1 <==> based on value (in that case will use time_TRAJ_START)
//-     : If time_TRAJ_START<0 then look for first time such that V(t,x)>0 to start trajectory reconstruction
int TRAJ_METHOD=1; 
double  time_TRAJ_START= -1.0;

// stopping criteria
double        min_TRAJ_STOP      =   -1.00; //- to stop traj reconstruction when val(x)<=min 
double        max_TRAJ_STOP      =  100.00; //- to stop traj reconstruction when val(x)>=max 
int           TARGET_STOP        =  0;      //- to stop traj reconstruction if g_target(x)<=0
inline double g_target(double t, const double* x){ return v0(x);}  

// adverse control case (only for COMMANDS=2 & TRAJ_METHOD=0): if ADVERSE_METHOD=1, allows to define some u2 adverse control
int           ADVERSE_METHOD    = 0;
void u2_adverse(double t, const double* x, C& u){}


//-----------------------------------------------
//- PARAMETERS FOR ADVANCED USERS
//- not necessary to modify these for basic usage
//------------------------------------------------

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
  SAVE_COUPE_ALL     = 1;
  SAVE_COUPE_ALL_STEP= 10; //- can be different from SAVE_VF_ALL_STEP; used for all "coupe" files
  SAVE_COUPE_FINAL   = 1;
  SAVE_COUPEEX_ALL   = 0;
  SAVE_COUPEEX_FINAL = 1;
  SAVE_COUPETOPT_ALL   = 0;
  SAVE_COUPETOPT_FINAL = 1;
  BINARY             = 1; //- for BINARY(1) or TEXT/ASCII(0) save & load 
  PRINTTRAJ          = 0; //- xterm printings of trajectory (default should be 0)
}

