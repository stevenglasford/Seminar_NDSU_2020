//----------------------------------------
//- advection + diffusion Equation, dimension 2
//----------------------------------------
//-    u_t - 0.5 Tr([ms][ms]^t D2u)  - b . Du + r u - ell = 0
//-    u(0,x)=u0(x)
//- Here with r=ell=0, ms : 2-2 constant matrix and b=constant vector
//- Notice that if PERIODIC={1,1} then no boundary function is used.
//----------------------------------------
//#include "data_default.h"	// do not modify for basic examples
//----------------------------------------

//- initial data and Vex here defined in order to have Vex for the (dirichlet) boundary conditions for testing purposes (in that case put PERIODIC={0,0}).
//---------------
//- initial data
//---------------
const double v0x1=-1*1.0;	//- advection speeds in x1 and x2 to define drift b=(-v0x1,-v0x2)^t.
const double v0x2=-1*2.0;

const double a1=2.0, a2= -1.0; //- Diffusion matrix: composantes of vect s1=sigma1
const double b1=1.0, b2=  0.0; //- Diffusion matrix: composantes of vect s2=sigma2

const double ma1= b2/(a1*b2-a2*b1), ma2=-a2/(a1*b2-a2*b1); //- inverse M of matrix sigma=[s1,s2]
const double mb1=-b1/(a1*b2-a2*b1), mb2= a1/(a1*b2-a2*b1);

const double w1=2*pi;
const double w2=4*pi;

const double coef11=1./(1.+1.), coef12= 1./(1.+2.);
const double coef21=1./(2.+1.), coef22=-1./(2.+2.);

inline double   v01(double r)
{
  return coef11*cos(w1*r)+coef12*cos(w2*r);
}
inline double   v02(double r)
{
  return coef21*cos(w1*r)+coef22*cos(w2*r);
}
inline double   v0(const double *x)
{
  double y1,y2;  //- y=Mx where M = (sigma)^-1
  y1=ma1*x[0]+mb1*x[1]; 
  y2=ma2*x[0]+mb2*x[1]; 
  return v01(y1) + v02(y2); 
}


//---------------------------
//- Exact solution (if known)
//---------------------------

const double K1=0.5;
const double K2=0.5;
inline double   vex1(double t, double r)
{
  return coef11*cos(w1*r)*exp(-K1*w1*w1*t) + coef12*cos(w2*r)*exp(-K2*w2*w2*t);
}
inline double   vex2(double t, double r)
{
  return coef21*cos(w1*r)*exp(-K1*w1*w1*t) + coef22*cos(w2*r)*exp(-K2*w2*w2*t);
}

inline double Vex(double t, const double* x)
{
  //-  vex1(t,a1 x1 + a2 x2)  is the solution of v_t - 1/2 Tr(ms1 ms1^T D^2 v) = 0 where  K=1/2 (a1^2+a2^2)^2
  double x0=x[0]-t*v0x1;
  double x1=x[1]-t*v0x2;
  double y1,y2;
  y1=ma1*x0+mb1*x1; 
  y2=ma2*x0+mb2*x1; 
  return vex1(t,y1) + vex2(t,y2);
}



const char    NAME[150]         = "data_SL_order2_2d_advectiondiffusion.h";  //- File name 
const int     DIM               = 2;                      //- space dimension

const int     COMMANDS          = 0; //- unused for SL
const int     OPTIM             = MAXIMUM;
//const int     OPTIM             = MINIMUM;


//----------------------------
//- SL method parameters:
//----------------------------
const int     METHOD            = MSL;
const int     TYPE_SCHEME       = EVO;                 //- space scheme  : 1 or STA ; 2 or EVO
const int     TYPE_RK           = RK1_EULER;           //- time scheme   : 1 or RK1_EULER, 2 or RK2_HEUN,  3 or RK2_PM

//--------------------
//- Special parmeter for method SL
//--------------------
      int     P_INTERMEDIATE    = 1;    //- number of discretisation steps for approximating each trajectory in method MSL 

//----------------------------
//- stopping criteria parameters:
//----------------------------
const double  EPSILON           = 0.0;
const int     MAX_ITERATION     = 100000;
const double                   T= 0.10;                  //- Terminal time


//----------------------------
//- discretization parameters:
//----------------------------

const int nn=1*4;
const int     NN                = nn*20;
const int     ND[DIM]           = {NN,NN};

const int     N                 = nn*10;
const double  DT                = T/float(N);                   //- if DT=0 the program will compute a DT based on space steps (if METHOD=MFD)
const double  XMIN[DIM]         = {  0. , 0.};
const double  XMAX[DIM]         = {  1. , 1.};            //- bound of the domain
const int     PERIODIC[DIM]     = {  1 ,  1 };            //- periodic mesh (1:periodic, 0:otherwise)
const int     MESH              = 1;                      //- 0 : xi = center of cell ; 1 : xi contains the boundary
const int     OBSTACLE          = 0;

const int     cDIM              = 2;
const int     NCD[cDIM]         = { 1 ,  1 };
const double  UMIN[cDIM]        = {  0., 1.};
const double  UMAX[cDIM]        = {2*pi, 1.};

const int     cDIM2             = 1;
const int     NCD2[cDIM2]       = { 1 };
const double  UMIN2[cDIM2]      = { 0.};
const double  UMAX2[cDIM2]      = { 1.};

//-----------------------------
//- BOUNDARY : 0 (Void, FD only) or 1 (Dirichlet, V=g_border, for FD/SL) or 2 (Vx=g_bordermix, for FD), or 3 (Vxx=0, for FD/SL) 
//-----------------------------
const int BOUNDARY=1;

// Dirichlet boundary condition (case BOUNDARY=1)
double g_border(double t, const double* x){
  return Vex(t,x);  // putting Vex for testing purposes
}

// Mixed Neumann bc :  ux=g(t,x,u) (case BOUNDARY=2)
double g_bordermix(double t, const double* x, double val){return 0.0;}

// xterm printings of trajectory (default should be 0)
int           PRINTTRAJ         = 0;        



//-----------------------------
//- mainloop parameters
//-----------------------------


const int     COMPUTE_MAIN_LOOP = 1;
const int     COMPUTE_VEX       = 1;
const int     COMPUTE_TOPT      = 0;
const int     TOPT_TYPE         = 0;        //- 0= min time pb, 1= exit time pb


const int     SAVE_VF_ALL        = 0;       //- 1 to save  V every SAVE_VFALL_STEP iterations in files "VFxx.dat"
const int     SAVE_VF_ALL_STEP   = 10;
const int     SAVE_VF_FINAL      = 1;

const int     CHECK_ERROR       = 1;        //- 1 to compute errors every CHECK_ERROR_STEP iterations
const int     CHECK_ERROR_STEP  = 10;


//-------------------
//- compute trajectory parameters
//-------------------
const int     TRAJPT            = 0;
const double  initialpoint[TRAJPT*DIM] = {};//- Initial point for HJB
// method of trajectory reconstruction and starting time
int           TRAJ_METHOD       = 0; 
double        time_TRAJ_START   = 0.0;
// stopping criteria
double        min_TRAJ_STOP     = 0.00;     //- to stop traj reconstruction when val(x) <=min (val=topt here)
double        max_TRAJ_STOP     = 100.00;   //- to stop traj reconstruction when val(x) >=max (val=topt here)
int           TARGET_STOP       = 1;        //- to stop traj reconstruction when g_target(x)<=0
inline double g_target(double t, const double* x){return 0.0;}
// adverse control case (only for COMMANDS=2 & TRAJ_METHOD=0)
int           ADVERSE_METHOD    = 0;
inline void u2_adverse(double t, const double* x, C& u){}

//-------------------
//- "coupe" : used only for exemples with d>=3 (cut into some plane)
//-------------------
const int     COUPE_DIMS[DIM] = {1 ,1 };
const double  COUPE_VALS[DIM] = {0.,0.};



//- FOR SL_EVO FUNCTION (1ST ORDER)
inline void dynamics(const double* x, C u, double t, double* res)
{
  res[0] = -0.0;
  res[1] = -0.0;
}


inline double distributed_cost(const double* x, C u, double t)
{
  return 0.;
}

//---------------------------------
//- two controls case (only for FD)
//---------------------------------
inline void dynamics2(const double* x, C u, C u2, double t, double* res)
{
  return;
}

inline double distributed_cost2(const double* arg, C u, C u2, double t)
{
  return 0.;
}

//---------------
//- obstacle function  (not used if OBSTACLE==0)
//---------------
inline double g_obstacle(double t, const double* arg)
{
  return 0.;
}


//-----------------------------------------------------------------
//- second order scheme parameters => FUNCTION secondorder_itSL_evo
//-----------------------------------------------------------------
const int ORDER = 2;
const int PARAMP= 2; //DIM;

//const double v0x1=-1.0;
//const double v0x2=-3.0;

inline void Sigma(const double* x, int k, C u, double t, double* res)
{
  double coef=sqrt(PARAMP);
  if(k==0){
    res[0] =  a1;
    res[1] =  a2;
  }
  else if(k==1){
    res[0] =  b1;
    res[1] =  b2;
  }
  res[0]=res[0]*coef;
  res[1]=res[1]*coef;
}

inline void Drift(const double* x, int k, C u, double t, double* res)
{
  res[0] = -v0x1;
  res[1] = -v0x2;
}
inline double funcR(const double* x, C u, double t)
{
  return 0.;
}
//- two point formula (basic example)
inline void funcY(const double* x, int k, double eps, C u, double t, double h, double* res)
{
  //- first order
  
  double resDrift[DIM], resSigma[DIM];
  Drift(x, k, u, t, resDrift);
  Sigma(x, k, u, t, resSigma);
  for(int d=0;d<DIM;d++) {
    res[d] = x[d] + resDrift[d]*h + eps*resSigma[d]*sqrt(abs(h));
  }
  
  //- second order
  /*
  double resDrift[DIM], resSigma[DIM];
  Drift(x, k, u, t, resDrift);
  Sigma(x, k, u, t, resSigma);
  double x1[DIM];
  for(int d=0;d<DIM;d++) 
    x1[d] = x[d] + resDrift[d]*h;
  double resDrift1[DIM];
  Drift(x1, k, u, t, resDrift1);
  for(int d=0;d<DIM;d++) {
    //res[d] = x[d] + resDrift[d]*h + eps*resSigma[d]*sqrt(abs(h)); //- first order 
    res[d] = x[d] + 0.5*(resDrift[d]+resDrift1[d])*h + eps*resSigma[d]*sqrt(abs(h)); //- second order
  }
  */

}

/*
//- 3 point formula (Kloeden & Platen 1999, Debrabant)
inline void Sigma(const double* x, int k, C u, double t, double* res)
{
  if(k==0){
    res[0] =  2.;
    res[1] = -1.;
  }
  else if(k==1){
    res[0] =  1.;
    res[1] =  0.;
  }
}
inline void Drift(const double* x, int k, C u, double t, double* res)
{
  res[0] = 0.;
  res[1] = 0.;
}
inline void funcY(const double* x, int k, double eps, C u, double t, double h, double* res)
{
  double resDrift[DIM], resSigma[DIM];
  Drift(x, k, u, t, resDrift);
  Sigma(x, k, u, t, resSigma);
  for(int d=0;d<DIM;d++) {
    res[d] = x[d] + resDrift[d]*h + resSigma[d]*(eps* sqrt(3.0*abs(h)));
  }
}
*/

/*
//---------------
//- initial data
//---------------
const double coef11=1./(1.+1.), coef12=1./(1.+2.);
const double coef21=1./(2.+1.), coef22=1./(2.+2.);
inline double   u01(double r)
{
  return coef11*cos(2.*pi*r)+coef12*cos(2.*pi*2.*r);
}
inline double   u02(double r)
{
  return coef21*cos(2.*pi*r)+coef22*cos(2.*pi*2.*r);
}
inline double   v0(const double * arg)
{
  return u01(arg[0]+2.*arg[1])+u02(-arg[1]);
}


//---------------------------
//- Exact solution (if known)
//---------------------------

inline double   vex1(double t, double r)
{
  return coef11*cos(2.*pi*r)*exp(-(2.*pi*1.)*(2.*pi*1.)*t/2.) +coef12*cos(2.*pi*2.*r)*exp(-(2.*pi*2.)*(2.*pi*2.)*t/2.);
}
inline double   vex2(double t, double r)
{
  return coef21*cos(2.*pi*r)*exp(-(2.*pi*1.)*(2.*pi*1.)*t/2.) +coef22*cos(2.*pi*2.*r)*exp(-(2.*pi*2.)*(2.*pi*2.)*t/2.);
}

inline double Vex(double t, const double* x)
{
  double x1=x[0],x2=x[1];
  x1=x1-t*v0x1;
  x2=x2-t*v0x2;
  return vex1(t,x1+2.*x2)+vex2(t,-x2);
}
*/

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
//- UNUSED
//--------------------
//inline double velocity(const double* arg, double t) {  return 0.;}    //FMM
//const int     SUB_N             = 9;                                  //UB
//const int     TARGET            = 1;                                  //UB
//double distance_to_target(const double *x){ return 0; }               //UB

// -------------------
// - BEGIN TEMPORARY (FOR DISTRIB VERSION)
// - [data_default:]
// -------------------
// special parameters for SL METHOD
const int     TYPE_STA_LOOP     = NORMAL; 
const int     INTERPOLATION     = BILINEAR;
//const int     ORDER             = 1;
//- special parameters for SL METHOD, SECOND ORDER (case ORDER=2)
//const int     PARAMP            = 0;
//inline double funcR             (const double* x, C u, double t) {return 0.;}
//inline void   funcY             (const double* x, int k, double eps, C u, double t, double h, double* res) {;}
inline double discount_factor             (const double* x){return 0.0;}
// ---------------------------
// - OTHER MAINLOOP PARAMETERS
// ---------------------------
const int     VALUE_PB           = 0; //- to compute value v(t,x) s.t. 0 = w(t,x,z)=v(t,x)-z (assumes DIM=d+1, x in R^d, z in R^1)
const int     SAVE_VALUE_ALL     = 0; //- to save at each iteration the value v(t,x) (assumes VALUE_PB=1)
const int     SAVE_VALUE_ALL_STEP= 0; //- to save final value v(t,x) (case VALUE_PB=1)
const int     SAVE_VALUE_FINAL   = 0; //- to save final value v(t,x) (case VALUE_PB=1)
const int     SAVE_VF_FINAL_ONSET= 0;
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
// - Restricted computational domain:
// - if COMPUTE_IN_SUBDOMAIN==1 will compute only for x s.t. g_domain(x)<0
// ----------------------------------
const int     COMPUTE_IN_SUBDOMAIN    = 0;
inline double g_domain(const double *x){return 0.0;}
// -------------------
// - END TEMPORARY 
// -------------------

// -------------------
// - For error computations
// -------------------
const double  C_THRESHOLD       = 1.00;     //- Error threshold for error computations

// -------------------
// - For trajectory reconstruction
// -------------------
//- method of reconstruction : 0 (based on tmin) or 1 (based on value)
//- starting time : t_TRAJ_START (starting time) may be used, in [0,T[,  if TRAJ_METHOD==1 
//int TRAJ_METHOD=0; 
//double time_TRAJ_START=0.0;


// -------------------
// - For parameter intialization/termination
// -------------------
void init_data(){};
void post_data(){};

// -------------------
// - prefix in front of the default output names 
// -------------------
char          FILE_PREFIX[]     = "";  //- prefix for the name of all data files ("" = no prefix is used.)

// --------------------------------
//- border: number of ghost cells in each direction (default is {2,2,...} for FD)
// --------------------------------
int BORDERSIZE[DIM] = {0,0};

