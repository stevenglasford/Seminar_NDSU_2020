//----------------------------------------
//- 2D 2controls 
//----------------------------------------
//-    v_t + max_{u} min_{u2}< f, nabla v >  = 0;    x in R^2
//-    v(0,x)=v0(x)
//----------------------------------------


//- Here programmation for two controls and Max(Min(..)) Hamiltonian:
//- set COMMANDS=2 and OPTIM=MAXMIN : dynamics2 and distributed_cost2 are used 
//- (for advanced users, put COMMANDS=0, Hnum function is used,
//- options 0/1/2/ are proposed in the end of this file)

const char    NAME[150]         = "data_DF_2d_2c.h (HJ Issacs equation, Mar. 2013 (OB, AD, HZ)";
const int     DIM               = 2;                      //- space dimension

const int     COMMANDS          = 2;           //- 0: local Hnum function ; 1,2: Hnum defined with dynamics and distributed_cost functions
const int     OPTIM             = MAXMIN;


//----------------------------
//- FD method parameters:
//----------------------------
const int     METHOD            = MFD;
const int     TYPE_SCHEME       = ENO2;        //- space scheme : 1 or LF ; 2 or ENO2; 3 or ENO3; 
const int     TYPE_RK           = RK1;         //- time scheme (ENO2 case) : 1 or ENO2RK1, 2 or ENO2RK2,  3 or ENO2RK3

//----------------------------
//- stopping criteria parameters:
//----------------------------
const double  EPSILON           = 0.0;
const int     MAX_ITERATION     = 100000;
const double                   T= 0.50;        //- Terminal time


//----------------------------
//- discretization parameters:
//----------------------------
const double  DT                = 0.0;         //- if DT=0 the program will compute a DT based on space steps.
const double  CFL               = 0.5;

const int     NN=80;
const int     ND[DIM]           = { NN  , NN };
const double  XMIN[DIM]         = { -2. , -2.};
const double  XMAX[DIM]         = {  2. ,  2.};           //- bound of the domain
const int     PERIODIC[DIM]     = {  1  ,  1 };           //- periodic mesh (1:periodic, 0:otherwise)
const int     MESH              = 1;                      //- 0 : xi = center of cell ; 1 : xi contains the boundary
const int     OBSTACLE          = 0;

const int     cDIM              = 1;
const int     NCD[cDIM]         = { 2 };
const double  UMIN[cDIM]        = {-1.};
const double  UMAX[cDIM]        = { 1.};


const int     cDIM2             = 1;
const int     NCD2[cDIM2]       = { 2 };
const double  UMIN2[cDIM2]      = {-1.};
const double  UMAX2[cDIM2]      = { 1.};

//-----------------------------
//- BOUNDARY : 0 (Void, FD only) or 1 (Dirichlet, V=g_border, for FD/SL) or 2 (Vx=g_bordermix, for FD), or 3 (Vxx=0, for FD/SL) 
//-----------------------------
const int BOUNDARY=1;

//-  Dirichlet boundary condition:
const double  VBORD             = 0.5; 
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
const int     SAVE_VF_ALL_STEP  = 10000;   //- (a large number will save only VF0.dat and final: VF1.dat)
const int     SAVE_VF_FINAL     = 1;

const int     CHECK_ERROR       = 1;       //- 1 to compute errors every CHECK_ERROR_STEP iterations
const int     CHECK_ERROR_STEP  = 10;
const int     CHECK_NEG_POINTS  = 0;       //- 1/0 : to print number of negative point in domain and domain + boundary (default : 0)

//-------------------
//- "coupe" : (cut into some plane ==> results in files "coupe{}{ex}{topt}.dat"
//-------------------
const int     COUPE_DIMS[DIM] = {1 ,1 };
const double  COUPE_VALS[DIM] = {0.,0.};


inline void dynamics(const double* x, C u, double t, double* res)
{
  res[0]=0.0;
  res[1]=0.0;
}
inline double distributed_cost(const double* arg, C u, double t)
{
  return 0.;
}


//----------------------------------------
//- these are used in the case COMMANDS=2
//----------------------------------------
inline void dynamics2(const double* x, C u, C u2, double t, double* res)
{
  //res[0]=-1.0; //+u[0];
  //res[1]= 0.0; //+u2[0];
  res[0]=u[0];
  res[1]=u2[0];
  return;
}

inline double distributed_cost2(const double* arg, C u, C u2, double t)
{
  return 0.;
}

inline void feedback(double t, const double* x, const double* p, C& u){  //- only needed for NBEE-OC
  //- input  : t,x,p
  //- output : u, that maximizes  max_u (- p. dynamics(t,x,u) - dist_cost(t,x,u)) 
  //- Example 
  //u[0]=1.0;
}


//---------------
//- obstacle function
//---------------
inline double g_obstacle(double t, const double* arg)
{
  return 0.;
}


//---------------
//- initial data
//---------------
const  double   X0=0.0, Y0=0.0, RAYON=1.;
inline double   q(double x) { return MIN(VBORD,x-RAYON); }
inline double   v0(const double * arg) {
   double nx=MAX(abs(arg[0]-X0), abs(arg[1]-Y0));
   return q(nx);
}

//---------------------------
//- Exact solution (if known)
//---------------------------
inline double Vex(double t, const double* arg) {
    double arg2[2];
    //double v = abs(arg[0]+arg[1])/sqrt(2.);
    //double u = abs(arg[0]-arg[1])/sqrt(2.);
    //arg2[0]  = max(0., abs(v)-t*sqrt(2.));
    //arg2[1]  = u;
    double u = arg[0];
    double v = arg[1];
    arg2[0]  = max(0., abs(u)-t);
    arg2[1]  = max(0., abs(v)+t);
    //double nx = sqrt(arg2[0]*arg2[0]+arg2[1]*arg2[1]);
    double nx = MAX(abs(arg2[0]),abs(arg2[1]));
    return q(nx);
}


//--------------------
//- Hamiltonian (to be used with the LF numerical hamiltonian here)
//--------------------
inline double H(const double* x, const double* p)
{
  double z = abs(p[0]) - abs(p[1]);
  //double z = p[1];
  return z;
}

//--------------------
//- Stability constants to be used for the numerical hamiltonian.
//--------------------
inline void compute_Hconst(double* aMAX, double t)
{
    aMAX[0] = 1.;
    aMAX[1] = 1.;
}

//--------------------
//- Hnum : option 0 : Numerical Hamiltonian (Lax Friedriech) in the case H is defined
//--------------------
inline double Hnum(const double t, const double* x, const double vi,  const double* Dv)
{
  double z=0.;
  double p[DIM];
  double amax[DIM];
  int i;
  for(i=0;i<DIM;i++){
    amax[i]=1.;
    p[i]=(Dv[2*i] + Dv[2*i+1])/2.;
    z += amax[i]*(Dv[2*i+1] - Dv[2*i])/2.;
  }
  return H(x,p)-z;
}

//--------------------
//- Hnum : option 1 :  Numerical Hamiltonian (Lax Friedriech) in the case H is defined - TWO VARIANTS FOR EXPERIMENTED USERS
//--------------------
/*
const int ncall1=2;
const double u1[ncall1][1]={{-1.},{1.}};
const int ncall2=10;
const double u2[ncall2][1]={{-1.},{-0.5},{0.0},{0.0},{0.0},{0.0},{0.0},{0.0},{0.5},{1.}};
//- local (for Hnum 2 control)
inline void local_dynamics2(const double* x, const double *u, const double *u2, double t, double* res)
{
  res[0]=u[0];
  res[1]=u2[0];
  return;
}
inline double local_distributed_cost2(const double* arg, const double *u, const double *u2, double t)
{
  return 0.;
}
inline double Hnum(const double* x, const double* dvn, double t)
//- for MAX(MIN()) approach: similar as:
//- double HJB_FD::Hnum_2C_MaxMin(const double* x, const double* dvn, const double t)
{
  int c1,c2,d;
  double z, maxA=-INF, minA, dvectdouble[DIM];
  for(c1=0;c1<ncall1;c1++){
    minA = INF;
    for(c2=0;c2<ncall2;c2++){
      local_dynamics2(x,u1[c1],u2[c2],t,dvectdouble);
      z=0.;
      for(d=0;d<DIM;d++)
        z += MAX(0.,dvectdouble[d])*dvn[2*d] + MIN(0.,dvectdouble[d])*dvn[2*d+1] + local_distributed_cost2(x,u1[c1],u2[c2],t);
      minA = MIN(minA,z);
    }
    maxA = MAX(maxA, minA);
  }
  return maxA;
}
*/

//-----------------------------------------------------------
//-  Hnum : option 2 : 2 controls, MAXMIN, with ANALYTIC FORM for the MIN  part:
//--------------------
/*
const int ncall1=2;
const double u1[ncall1][1]={{-1.},{1.}};
const int ncall2=10;
const double u2[ncall2][1]={{-1.},{-0.5},{0.0},{0.0},{0.0},{0.0},{0.0},{0.0},{0.5},{1.}};

inline void local_dynamics2(const double* x, const double *u, const double *u2, double t, double* res)
{
  res[0]=u[0];
  res[1]=u2[0];
  return;
}
inline double local_distributed_cost2(const double* arg, const double *u, const double *u2, double t)
{
  return 0.;
}
inline double Hnum(const double* x, const double* dvn, double t)
//- for MAX(MIN()) approach with ANALYTIC FORM for the MIN part:
//- double HJB_FD::Hnum_2C_MaxMin(const double* x, const double* dvn, const double t)
{
  double z, maxA=-INF;
  double f1;
  for(int c1=0;c1<ncall1;c1++){
    // //- minimization loop:
    // z=+INF;
    // double dvectdouble[DIM];
    // int d;
    // for(int c2=0;c2<ncall2;c2++){
    //   local_dynamics2(x,u1[c1],u2[c2],t,dvectdouble);
    //   double zloc=0.;
    //   for(d=0;d<DIM;d++)
    //     zloc += MAX(0.,dvectdouble[d])*dvn[2*d] + MIN(0.,dvectdouble[d])*dvn[2*d+1];
    //     //+ local_distributed_cost2(x,u1[c1],u2[c2],t);
    //  z= MIN(z,zloc);
    //}
    // 
    
    //- analytic minimization:
    f1=u1[c1][0];
    z = MAX(0.,f1)*dvn[0] + MIN(0.,f1)*dvn[1]  + 1.* MIN(dvn[2],-dvn[3]);
   
    
    maxA = MAX(maxA, z);
  }
  return maxA;
}
*/


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
double        min_TRAJ_STOP     =   0.00;  //- to stop traj reconstruction when val(x) <=min (val=topt here)
double        max_TRAJ_STOP     = 100.00;  //- to stop traj reconstruction when val(x) >=max (val=topt here)
int           TARGET_STOP       = 1;       //- to stop traj reconstruction when g_target(x)<=0
inline double g_target(double t, const double* x){ return 0.0;}

// adverse control case (only for COMMANDS=2): if ADVERSE_METHOD=1, allows to define some u2 adverse control
int           ADVERSE_METHOD    = 0;
void u2_adverse(double t, const double* x, C& u){}  // to define u as some adverse control (of type u2)


// xterm printings of trajectory (default should be 0)
int           PRINTTRAJ         = 0;

//--------------------
//- UNUSED
//--------------------
//void init_data(){};
//const int     TYPE_STA_LOOP     = 0;
//const int     INTERPOLATION     = 0;
//const int     P_INTERMEDIATE    = 0;
//inline double velocity(const double* arg, double t) {  return 0.;}
//const int     SUB_N             = 9;
//const int     TARGET            = 1;
//const int     DEB               = 0;
//const int     DEC               = 1;
//const double  DIST_THRES_COEFF  = 2.;
//const int     ORDER             = 2;
//const int     PARAMP            = DIM;
//inline double funcR             (const double* x, C u, double t) {  return 0.; }
//inline void   funcY             (const double* x, int k, double eps, C u, double t, double h, double* res) {  return; }
//inline double distance_to_target(const double *x) { return 0.; }
//inline void   set_target_boundaries(double* target_xmin, double* target_xmax) {}



// -------------------
// - BEGIN TEMPORARY (FOR DISTRIB VERSION)
// -------------------
// special parameters for SL METHOD
const int     TYPE_STA_LOOP     = NORMAL;
const int     INTERPOLATION     = BILINEAR;
const int     ORDER             = 1;
// special parameters for SL METHOD, SECOND ORDER (case ORDER=2)
const int     PARAMP            = 0;
inline double funcR             (const double* x, C u, double t) {return 0.;}
inline void   funcY             (const double* x, int k, double eps, C u, double t, double h, double* res){}
inline double discount_factor   (const double* x){return 0.0;}
const int     P_INTERMEDIATE    = 1;    //- number of discretisation steps for approximating each trajectory in method MSL 

// ---------------------------
// - OTHER MAINLOOP PARAMETERS
// ---------------------------

const int     VALUE_PB           = 0; //- to compute value v(t,x) s.t. 0 = w(t,x,z)=v(t,x)-z (assumes DIM=d+1, x in R^d, z in R^1)
const int     SAVE_VALUE_ALL     = 0; //- to save at each iteration the value v(t,x) (assumes VALUE_PB=1)
const int     SAVE_VALUE_ALL_STEP= 0; //- to save final value v(t,x) (case VALUE_PB=1)
const int     SAVE_VALUE_FINAL   = 0; //- to save final value v(t,x) (case VALUE_PB=1)

//- note : {VF,Vex,topt}.dat savings use the same step as VF savings : SAVE_VF_ALL_STEP
//- note : {coupe,coupeex,coupetopt}.dat savings use the same step as coupe savings: SAVE_COUPE_ALL_STEP
const int     SAVE_VEX_ALL       = 0;
const int     SAVE_VEX_FINAL     = 1;
const int     SAVE_TOPT_ALL      = 0;
const int     SAVE_TOPT_FINAL    = 1;
const int     SAVE_COUPE_ALL     = 0;
const int     SAVE_COUPE_ALL_STEP= 1; //- can be different from SAVE_VF_ALL_STEP
const int     SAVE_COUPE_FINAL   = 1;
const int     SAVE_COUPEEX_ALL   = 0;
const int     SAVE_COUPEEX_FINAL = 1;
const int     SAVE_COUPETOPT_ALL = 0;
const int     SAVE_COUPETOPT_FINAL=1;
const int     SAVE_VF_ONSET_ALL    = 0; //- the iteration step for savings is defined by "SAVE_VF_ALL_STEP" 
const int     SAVE_VF_ONSET_FINAL  = 0; //- to save final value on some user defined points (X_user.txt)
const int     SAVE_VEX_ONSET_ALL   = 0;
const int     SAVE_VEX_ONSET_FINAL = 0;
const int     SAVE_TOPT_ONSET_ALL  = 0;
const int     SAVE_TOPT_ONSET_FINAL= 0;
//char          XUSER_DATA_FILE[]= "Xuser.txt";
//const char      XFile_user0   []= "X_user";   // to be completed by ".txt" extention (to be distinguished from .dat files)
//const char      VFile_user0   []= "XV_user";  // will be completed by ".dat" or "n.dat" extention 
char          XuserFile0[]      = "Xuser";     // to be completed by ".txt" extention (to be distinguished from .dat files)
char          XuserVFile0[]     = "XuserV";    // will be completed by ".dat" or "n.dat" extention 
char          XuserVexFile0[]   = "XuserVex";  // will be completed by ".dat" or "n.dat" extention 
char          XuserToptFile0[]  = "XuserVtopt";// will be completed by ".dat" or "n.dat" extention 

const int     SAVE_COUPE_INTERPOLATION=0; //- if 1 then uses interpolation. Use 0 for testing purposes (compare VFxx.dat and coupexx.dat files)


      int     BINARY             = 0; //- for BINARY(1) or TEXT/ASCII(0) save & load 

// -------------------
// - format for saving .dat files (VF.dat,VFxx.dat)  1=default(coord+val);  0=only values
// -------------------
const int     FORMAT_FULLDATA         = 1; //- furthermore add (at left) the list of indexes corresponding to the current point {files VF.dat, etc.}
const int     SHOW_COORDINATES        = 0; //- will furthermore add (at right) the list of coordinates on each line of the {VF.dat, etc.} files 
// -------------------
// - END TEMPORARY (FOR DISTRIB VERSION)
// -------------------

// -------------------
//  OBSTACLE g tilde
// -------------------
const int     OBSTACLE_TILDE          = 0;
inline double g_obstacle_tilde(double t, const double* arg){ return 0.;}
// --------------------
// - PRECOMPUTE_OBSTACLE (to precompute obstacle terms : g_obstacle() and g_obstacle_tilde(), on the grid)
// --------------------
const int     PRECOMPUTE_OBSTACLE     = 0;
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
const double  C_THRESHOLD       = 0.1;     //- Error threshold for error computations

// -------------------
// - For parameter intialization/termination
// -------------------
void init_data(){};
void post_data(){};

//--------------------------------
//- border: number of ghost cells in each direction (default is {2,2,...} for FD)
//--------------------------------
int BORDERSIZE[DIM] = {2,2};

// -------------------
// - prefix in front of the default output names 
// -------------------
char          FILE_PREFIX[]     = "";      //- prefix for the name of all data files ("" = no prefix is used.)
char          EXTERNAL_FILE_PREFIX[]= "";  //- prefix for loading files (when using EXTERNALV0 and EXTENAL_TOPT)

// --------------------------------
// - MORE ADVANCED PARAMETERS: 
// --------------------------------
const int     EXTERNALV0        = 0;           //- if 1 then starts computation from data in VF.dat and not from the v0 function
const int     EXTERNAL_TOPT     = 0;           //- if 0 (default):  topt is set 0 when min(max(v0,g_obstacle),g_obstacle_tilde)<=0), if 1 load topt.dat
                                               //- (this parameter is ignored if EXTERNALV0=0)
const int     PRECOMPUTE_COORDS = 1;           //- to precompute mesh coordinates (faster but needs more memory).

// -------------------
// - parameter intialization/termination & others <2018>
// -------------------
void init_data_general(){}
void post_data_general(){}
int step_HJB;					//- general iterator for successive HJB computations (do not touch for standard HJB)
int beginStep_HJB=0;				//- starting at beginStep_HJB
int endStep_HJB  =0; 				//- ending   at endStep_HJB 

//- PMP link
void dynamicsGrad(const double* x, C u, double s, double** res){}

