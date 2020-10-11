//----------------------------------------
//- Rotation example , dimension 2
//----------------------------------------
//-    v_t + max_a < - f(x,a), nabla v >  = 0;    f(x_1,x_2,a)= a (-x_2,x_1);  a in {0,1}
//-    v(0,x)=v0(x)
//----------------------------------------
#include "data_default.h"	  // for basic examples
//----------------------------------------

const char    NAME[]            = "data_basicmodel.h, Jun. 2014 (O. Bokanowski, H. Zidani)";
const int     DIM               = 2;   //- space dimension

const int     COMMANDS          = 1;   //- 0: local Hnum function ; 1/2: Hnum defined with dynamics and distributed_cost functions
const int     OPTIM             = MAXIMUM;

//----------------------------
//- FD method parameters:
//----------------------------
const int     METHOD            = MFD;         //- Method : Finite differences(MFD)/Semi-Lagrangian(MSL)
const int     TYPE_SCHEME       = ENO2;        //- space scheme : 1 or LF ; 2 or ENO2 ; 3 or ENO3
const int     TYPE_RK           = RK1;         //- time scheme (ENO2 case) : 1 or ENO2RK1, 2 or ENO2RK2,  3 or ENO2RK3
//const int     TYPE_SCHEME       = NBEE;        //- space scheme : 1 or LF ; 2 or ENO2 ; 3 or ENO3

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
      double  CFL               = 0.4;

const int     NN=3*25;
      int     ND[DIM]           = { NN , NN};
      double  XMIN[DIM]         = { -2. , -2.};
      double  XMAX[DIM]         = {  2. ,  2.};           //- bound of the domain

const int     PERIODIC[DIM]     = {  0  ,  0 };           //- periodic mesh (1:periodic, 0:otherwise)
      int     MESH              = 0;                      //- 0 : xi = center of cell ; 1 : xi contains the boundary

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

//- Dirichlet boundary condition (case BOUNDARY=1)
const double  VBORD             = 0.2;

double g_border(double t, const double* x){
  return VBORD;
}

//- Mixed Neumann bc :  vx=g(t,x,v) (case BOUNDARY=2)
double g_bordermix(double t, const double* x, double val){return 0.0;}


//-----------------------------
//- mainloop parameters
//-----------------------------
const int     COMPUTE_MAIN_LOOP = 1;
const int     COMPUTE_VEX       = 1;
const int     COMPUTE_TOPT      = 1;
const int     TOPT_TYPE         = 0;       //- 0= min time , 1= exit time

const int     SAVE_VF_ALL       = 1;       //- 1 : to save "VFn.dat" every  SAVE_VF_ALL_STEP iterations (useful for movie)
const int     SAVE_VF_ALL_STEP  = 100;     //- (set to a large number to save only VF0.dat and ev. a final: VF1.dat)
const int     SAVE_VF_FINAL     = 1;       //- 1 : to save furthermore the last VF file as "VF.dat" 

const int     CHECK_ERROR       = 1;       //- 1 to compute errors every CHECK_ERROR_STEP iterations (if COMPUTE_VEX=1)
const int     CHECK_ERROR_STEP  = 100;
const int     CHECK_NEG_POINTS  = 0;       //- 2/1/0 : to print number of negative point in domain and domain + boundary (default : 0)

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
  //double X0=1.0, Y0=0.0, RAYON=0.05; //double r020=RAYON*RAYON; double resmin=1-pow((1-0)/(1-r020),4); //- test
  //double nx=sqrt( (x[0]-X0)*(x[0]-X0) +(x[1]-Y0)*(x[1]-Y0))-RAYON; 	//- ball target
  //double nx=max( abs(x[0]-X0) -RAYON ,abs(x[1]-Y0) -RAYON);		//- square target
  //return min(VBORD,nx);
  double r02=RAYON*RAYON, nx2=(x[0]-X0)*(x[0]-X0) +(x[1]-Y0)*(x[1]-Y0); //- smoother ball target
  double res=VBORD*(1-pow(max(0.0, (1-nx2)/(1-r02)),4));
  return res;
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
  return 0.0;
}

inline void feedback(double t, const double* x, const double* p, C& u){  //- only needed for NBEE-OC
  //- input  : t,x,p
  //- output : u, that maximizes  max_u (- p. dynamics(t,x,u) - dist_cost(t,x,u)) 
  //- Example rotation (no control)
  u[0]=1.0;
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

  //- rotation pure:
  /*
  double x2[2];  // here we compute v0(R_{- 2pit}x) (backward rotation)
  double t1=-t;
  //- Rotation backward de x=(x[0],x[1]):
  x2[0]=cos(2*pi*t1)*x[0] - sin(2*pi*t1)*x[1];
  x2[1]=sin(2*pi*t1)*x[0] + cos(2*pi*t1)*x[1];
  return v0(x2);
  */

  //- Euclidean distance to curved segment [0,th=2*pi*t] starting from X0,Y0
  double X0=1.0, Y0=0.0, RAYON=0.25;
  double y[2];
  //- computation of projection of x on curved segment:
  double th=atan2(x[1],x[0]);
  if (x[0]<=0 && x[1]<=0){
    th=th+2*pi;
  }
  th=min(max(th,0.0),2*pi*t);  //- projection of th on [0, 2*pi*t]
  y[0]=cos(th)*X0 - sin(th)*Y0;
  y[1]=sin(th)*X0 + cos(th)*Y0;
  double r=sqrt(pow(x[0]-y[0],2.0)+pow(x[1]-y[1],2.0));
  double r02=RAYON*RAYON, nx2=r*r; //- linked to v0 : smoother ball target 
  double res=VBORD*(1-pow(max(0.0, (1-nx2)/(1-r02)),4)); //- linked to v0 
  return res;

}


//-----------------------------------------------
//- PARAMETERS : compute trajectory 
//-----------------------------------------------
// method of trajectory reconstruction (0: based on topt; 1: based on "value"; 2: misc. methods, based on "topt")
int           TRAJ_METHOD       = 0; 

// list of starting points
//const int     TRAJPT     = 0;  const double  initialpoint[TRAJPT*DIM] = {}; //- no initial point 
const int     TRAJPT     = 1;  const double  initialpoint[TRAJPT*DIM] = {-1.,0.0}; //- 1 initial point 

// stopping criteria
double        min_TRAJ_STOP     =   0.00;  //- to stop traj reconstruction when val(x) <=min_TRAJ_STOP
double        max_TRAJ_STOP     = 100.00;  //- to stop traj reconstruction when val(x) >=max_TRAJ_STOP
int           TARGET_STOP       = 1;       //- to stop traj reconstruction when g_target(x)<=0
inline double g_target(double t, const double* x){ return v0(x);}  

//- parameters for TRAJ_METHOD=0 : adverse control case (only for COMMANDS=2 & TRAJ_METHOD=0)
int           ADVERSE_METHOD    = 0;
inline void u2_adverse(double t, const double* x, C& u){}

//- parameters for TRAJ_METHOD=1 : 
double        time_TRAJ_START   = 0.00;    //- Starting time, only for TRAJ_METHOD=1 (reconstruction by value)



//-----------------------------------------------
//- MORE ADVANCED PARAMETERS
//-----------------------------------------------

//--------------------
//- Stability constants to be used for the numerical hamiltonian (bounds for dH/dp_i, H=H(x,p))
//--------------------
//- Case COMMANDS=0: numerical Hamiltonian function Hnum.
inline void compute_Hconst(double* aMAX, double t){};
inline double Hnum(const double t, const double* x, const double vi, const double* dv){return 0.0;}

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
  SAVE_COUPE_ALL_STEP = 50; //- can be different from SAVE_VF_ALL_STEP; used for all "coupe" files
  SAVE_COUPE_ALL      = 1;
  SAVE_COUPE_FINAL    = 1;
  SAVE_COUPEEX_ALL    = 1;
  SAVE_COUPEEX_FINAL  = 1;
  SAVE_COUPETOPT_ALL  = 1;
  SAVE_COUPETOPT_FINAL= 1;
  SAVE_TOPT_FINAL     = 1;

  BINARY              = 1;  //- for BINARY(1) or TEXT/ASCII(0) save & load 
  FORMAT_FULLDATA     = 1;  //- furthermore add (at left) the list of indexes corresponding to the current point {files VF.dat, etc.}
  SHOW_COORDINATES    = 0;  //- will furthermore add (at right) the list of coordinates on each line of the {VF.dat, etc.} files 
  EXTERNALV0          = 0;  //- if 1 then starts computation from data in VF.dat and not from the v0 function

  SAVE_VF_ONSET_ALL    = 0; //- savings iteration step for {VF,Vex,topt (& onset-files) is defined by "SAVE_VF_ALL_STEP" 
  SAVE_VF_ONSET_FINAL  = 1; //- to save final value on some user defined points (such as Xuser.txt)
  SAVE_VEX_ONSET_ALL   = 0;
  SAVE_VEX_ONSET_FINAL = 1;
  SAVE_TOPT_ONSET_ALL  = 0;
  SAVE_TOPT_ONSET_FINAL= 1;

  NBPOINTS_Xuser=3; //100*100;
}

