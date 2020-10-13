//----------------------------------------
//-  steady equation  L u + H(x,nabla u) =0, d=1/2 
//----------------------------------------
//- for solving  L*u + H(x,nabla u) = 0
//- where H(x,nabla u) = max(-f(x,a) nabla u - l(x,a))
//- here L=1;  f(x,a)=a in [-1,1]  l(x,a)=CC;
//-    L u  + max_a (-f(x,a) . nabla u - l(x,a))  = 0
//-    u(1) = u(-1) = 0;
//----------------------------------------
//#include "data_default.h"	// do not modify for basic examples
//----------------------------------------


const char    NAME[150]         = "data_FD_steady.h, May. 2013\nAuthors  : O. Bokanowski, H. Zidani";
//const int     DIM               = 1;   //- space dimension: see below

const int     COMMANDS          = 1;   //- 0: local Hnum function ; 1: Hnum defined with dynamics and distributed_cost functions
//const int     OPTIM             = MINIMUM; 
const int     OPTIM             = MAXIMUM; 


//----------------------------
//- FD method parameters:
//----------------------------
const int     METHOD            = MFD;         //- MFD = Method Finite Differences
const int     TYPE_SCHEME       = LF;          //- space scheme : 1 or LF ; 2 or ENO2
const int     TYPE_RK           = RK1;         //- time scheme (ENO2 case) : 1 or ENO2RK1, 2 or ENO2RK2,  3 or ENO2RK3


//----------------------------
//- stopping criteria parameters:
//----------------------------
const double  EPSILON           = 1.e-4;//1.e-8; //1e-10; //- 0.0 not OK ???!
const int     MAX_ITERATION     = 1e4;
const double                   T= 5.00;        //- Terminal time


//----------------------------
//- discretization parameters:
//----------------------------
const double                  DT= 0.00;  //- if DT=0, MSL : will not work.
                                         //-          MFD : the program will compute a DT based on space steps.


//----------------------
//- CHOOSING THE EXAMPLE 
//----------------------

//- EXAMPLE WITH DIM=2:
/*
const int     DIM               = 2;            //- space dimension
const int                     NN= 100-2;
      int     ND[DIM]           = { NN , NN};
const double  XMINout[DIM]      = { -1. , -1.}; //- bounds of the domain
const double  XMAXout[DIM]      = {  1. ,  1.};      
const double  dx=(XMAXout[0]-XMINout[0])/NN;
//const double  dx=0.1;
      double  XMIN[DIM]         = { XMINout[0]+dx, XMINout[1]};
      double  XMAX[DIM]         = { XMAXout[0]-dx, XMAXout[1]};
const int     COUPE_DIMS[DIM]   = { 1 , 1 };    //- "coupe"
const double  COUPE_VALS[DIM]   = { 0., 0.};
const int     PERIODIC[DIM]     = { 0 , 0 };    //- periodic mesh (1:periodic, 0:otherwise)
//--------------------------------
//- border: number of ghost cells in each direction (default is {2,2,...} for FD)
//--------------------------------
int BORDERSIZE[DIM] = {2,2};
*/

//- EXAMPLE WITH DIM=1:

const int     DIM               = 1;         //- space dimension
const int                     NN= 100;       //- total of NN points. If MESH=1. Put NN=20.==> 20 intervals, 21 points. 
const int     ND[DIM]           = { NN-2}; 	   
const double  XMINout[DIM]      = { -1.};       
const double  XMAXout[DIM]      = {  1.};      
const double  dx=(XMAXout[0]-XMINout[0])/NN;
const double  XMIN[DIM]         = { XMINout[0]+dx};      
const double  XMAX[DIM]         = { XMAXout[0]-dx};      

const int     COUPE_DIMS[DIM]   = { 0  };    //- "coupe" : used only for exemples with d>=2 (cut into some plane)
const double  COUPE_VALS[DIM]   = { 0. }; 
const int     PERIODIC[DIM]     = { 0  };    //- periodic mesh (1:periodic, 0:otherwise)
//const int     ND[DIM]           = { NN }; 	 //- will have NN+1 points, boundary included.
//const double  XMIN[DIM]         = { -1.0};     //- bounds of the domain
//const double  XMAX[DIM]         = {  1.0};      
//--------------------------------
//- border: number of ghost cells in each direction (default is {2,2,...} for FD)
//--------------------------------
int BORDERSIZE[DIM] = {2};


const int     MESH              = 1;      //- 0 : nodes = center of cell ; 1 : nodes contain the boundary

const int     cDIM              = 1;
const int     NCD[cDIM]         = { 2 };
const double  UMIN[cDIM]        = { -1.0};
const double  UMAX[cDIM]        = { +1.0};

const int     cDIM2             = 1;
const int     NCD2[cDIM2]       = { 1 };
const double  UMIN2[cDIM2]      = { 0.};
const double  UMAX2[cDIM2]      = { 1.};


//-----------------------------
//- mainloop parameters
//-----------------------------
const int     COMPUTE_MAIN_LOOP = 1;
const int     COMPUTE_VEX       = 1;
const int     COMPUTE_TOPT      = 0;
const int     TOPT_TYPE    = 0;// 0= min time , 1= exit time

const int     SAVE_VF_ALL       = 1;       //- 1 to save "VFn.dat" every  SAVE_VFALL_STEP iterations
const int     SAVE_VF_ALL_STEP  = 10000;   //- (a large number will save only VF0.dat and final: VF1.dat)
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



//---------------
//- initial data
//---------------

inline double   v0(const double * arg) {
  //double X0=0.0, Y0=0.0, RAYON=0.5;
  //double nx=sqrt( (arg[0]-X0)*(arg[0]-X0) +(arg[1]-Y0)*(arg[1]-Y0))-RAYON; //- ball target
  //double nx=max( abs(arg[0]-X0) -RAYON ,abs(arg[1]-Y0) -RAYON); //- square target
  //double nx=abs(arg[0]-X0)/XMAX[0]; //- square target
  //double res=MIN(VBORD,nx);
  //return nx;
  //return 1-abs(arg[0]);
  return 0.0;

}


//- for solving  LAMBDA*u + H(x,nabla u) = 0
//- where H(x,nabla u) = max(-f(x,a) nabla u - l(x,a))
//- here L=1;  f(x,a)=a;  l(x,a)=0;

const double xmin=XMINout[0];
const double xmax=XMAXout[0];
const double LAMBDA=1.0;
const double CC=2.0;

//- HERE PUT LAMBDA:
inline double discount_factor(const double* x){
  return LAMBDA;
}

inline void dynamics(const double* x, C u, double t, double* res) {
  //double  X=x[0], Y=x[1];
  //f0=-(-2*pi*Y);
  //f1=-(+2*pi*X);
  res[0]=u[0]*1.0;
  if(DIM==2) res[1]=0.0;
  //res[1]=0.0;
}

inline double distributed_cost(const double* arg, C u, double t) {
  return CC;
}



//- unused functions since only 1 control here (COMMANDS=1)
inline void dynamics2(const double* x, C u, C u2, double t, double* res) { return; }
inline double distributed_cost2(const double* arg, C u, C u2, double t) { return 0.; }


//---------------------------
//- Exact solution (if known)
//---------------------------
inline double Vex(double t, const double* arg)
{
  //double arg2[2];
  //double t1=-t;
  //arg2[0]=cos(2*pi*t1)*arg[0] - sin(2*pi*t1)*arg[1];
  //arg2[1]=sin(2*pi*t1)*arg[0] + cos(2*pi*t1)*arg[1];
  //return v0(arg2);
  double x=arg[0];
  double L1=CC/exp(-xmin);
  double L2=CC/exp( xmax);
  double res =min(CC-L1*exp(-x), CC-L2*exp(x));
  return res;
}


//-----------------------------
//- BOUNDARY : 0 (Void, FD only) or 1 (Dirichlet, V=g_border, for FD/SL) or 2 (Vx=g_bordermix, for FD), or 3 (Vxx=0, for FD/SL) 
//-----------------------------
const int BOUNDARY=1;


// Dirichlet boundary condition (case BOUNDARY=1)
const double  VBORD             = 0.0;
double g_border(double t, const double* arg){
  //return Vex(0.,arg);
  return 0.0;
  /*
  double eps=1e-5;
  if (DIM==1){ 
    return 0.0; 
  }
  if (DIM==2){
    if (abs(arg[0])>1-eps)
      return 0.0;
    else
      return Vex(0.,arg);
  }
  */
}


// Mixed Neumann bc :  ux=g(t,x,u) (case BOUNDARY=2)
double g_bordermix(double t, const double* x, double val){return 0.0;}

// xterm printings of trajectory (default should be 0)
int           PRINTTRAJ         = 0;        



//--------------------------------------------------
//- These parameters are not used for the SL method
//--------------------------------------------------
inline void compute_Hconst(double* aMAX, double t){} 
inline double Hnum(const double t, const double* x, const double v, const double* Dv){return 0.0;}
double  CFL               = 0.8;		//- not used for SL method


//-----------------------------------------------
//- PARAMETERS FOR ADVANCED USERS
//- not necessary to modify these for basic usage
//------------------------------------------------
const int     EXTERNALV0        = 0;

//--------------------
//- OBSTACLE g
//--------------------
const int     OBSTACLE          = 0;
inline double g_obstacle(double t, const double* arg){
  return 0.;
}

// -------------------
// - BEGIN TEMPORARY (FOR DISTRIB VERSION)
// -------------------
//const double  VBORD=0.0; 	//- obsolete

// special parameters for SL METHOD
const int     TYPE_STA_LOOP     = NORMAL; 
const int     INTERPOLATION     = BILINEAR;
const int     ORDER             = 1;
// special parameters for SL METHOD, SECOND ORDER (case ORDER=2)
const int     PARAMP            = 0;
inline double funcR             (const double* x, C u, double t) {  return 0.; }
inline void   funcY             (const double* x, int k, double eps, C u, double t, double h, double* res){}
// ---------------------------
// - OTHER MAINLOOP PARAMETERS
// ---------------------------
const int     VALUE_PB           = 0; //- to compute value v(t,x) s.t. 0 = w(t,x,z)=v(t,x)-z (assumes DIM=d+1, x in R^d, z in R^1)
const int     SAVE_VALUE_ALL     = 0; //- to save at each iteration the value v(t,x) (assumes VALUE_PB=1)
const int     SAVE_VALUE_ALL_STEP= 0; //- to save final value v(t,x) (case VALUE_PB=1)
const int     SAVE_VALUE_FINAL   = 0; //- to save final value v(t,x) (case VALUE_PB=1)
const int     SAVE_VF_FINAL_ONSET= 0; //- to save final value on some user defined points (X_user.txt)
// -------------------
// - format for saving .dat files (VF.dat,VFxx.dat)  1=default(coord+val);  0=only values
// -------------------
const int     FORMAT_FULLDATA         = 1;
// --------------------
// - OBSTACLE g tilde
// --------------------
const int     OBSTACLE_TILDE          = 0;
inline double g_obstacle_tilde(double t, const double* arg){ return 0.;}
// --------------------
// - PRECOMPUTE_OBSTACLE (to precompute obstacle terms)
// --------------------
const int PRECOMPUTE_OBSTACLE=0;

// ----------------------------------
// - RESTRICTED COMPUTATIONAL DOMAIN:  
// - if RESTRICTED_DOMAIN==1 will compute only for x s.t. g_domain(x)<=0
// ----------------------------------
const int     COMPUTE_IN_SUBDOMAIN    = 0;
inline double g_domain(const double *x){return 0.0;}
// -------------------
// - END TEMPORARY 
// -------------------

// -------------------
// - For error computations
// -------------------
const double  C_THRESHOLD       = 2.0;     //- Error threshold for error computations

// -------------------
// - For parameter intialization
// -------------------
void post_data(){};
void init_data(){
  /*
  const char   file[]= "param.dat";
  FILE * pFile;
  pFile = fopen (file,"w");
  if (pFile==NULL){
    cout << "Impossible d'ouvrir le fichier en écriture !" << endl;
    exit(1);
  }
  else {
    double val=0.87877;
    fprintf(pFile,"%8.5f\n",  val);
  }
  fclose (pFile);
  */
};

// -------------------
// - prefix in front of the default output names 
// -------------------
char          FILE_PREFIX[]     = "";  //- prefix for the name of all data files ("" = no prefix is used.)

//--------------------
//- Special parmeter for method SL
//--------------------
const int     P_INTERMEDIATE    = 1;    //- number of discretisation steps for approximating each trajectory in method MSL 

//--------------------
//- Other parameters
//--------------------
const int     PRECOMPUTE_COORDS = 1;                      //- to precompute mesh coordinates (faster but needs more memory).


