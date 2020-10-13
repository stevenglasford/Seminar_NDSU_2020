//----------------------------------------
//- Advection Equation, dimension 3 (2d advection in x,y + 1 in z - for testing purposes)
//----------------------------------------
//-   u_t + < f, nabla u >  = 0;    f= ( -(-y,x), 1);
//-   u(0,x)=u0(x)
//- 
//- Abstract:
//-   2 rotations of 360 deg. (in x,y). 
//-   Good test to compare schemas ENO2, ENO3, NBEE, NBEE_OCOPT : observe numerical diffusion (or anti-diffusion...)
//-   (NBEE refers to a particular unpublished scheme; it is related to NBEE Bokanowski-Zidani JCP 2007) 
//-
//-   IMPORTANT: for ENO2 scheme use COMMANDS=0 or 1 (=> 0:H=Hnum definition)
//-   IMPORTANT: for NBEE scheme use COMMANDS=1 only (=> 1:H=max(-f(x,a).nabla u) definition)
//----------------------------------------
#include "data_FD_3d_rot_default.h"
//----------------------------------------

const char    NAME[200]         = "data_FD_3d_rot.h: 3D example (rotation, bound-test, coup-test), Feb/Aug 2016;2018 [OB]";
const int     DIM               = 3;   //- space dimension   

const int     COMMANDS          = 1;   //- 0: user Hnum;  1: Hnum from dynamics()/distributed_cost(); (2: for min/max pbs.)
const int     OPTIM             = MAXIMUM;


//----------------------------
//- DF method parameters:
//----------------------------
const int	METHOD	    = MFD;              //- Method : Finite differences(MFD)/Semi-Lagrangian(MSL)
//const int	TYPE_SCHEME = ENO2;       double CFL=0.7; //- space scheme : 1 or LF ; 2 or ENO2; 3 or ENO3;  21 or NBEE
//const int	TYPE_SCHEME = ENO3;       double CFL=0.7; //- space scheme : 1 or LF ; 2 or ENO2; 3 or ENO3;  21 or NBEE
const int	TYPE_SCHEME = NBEE;       double CFL=1.0; //- space scheme : 1 or LF ; 2 or ENO2; 3 or ENO3;  21 or NBEE (CFL=1 then)
//const int	TYPE_SCHEME = NBEE_OCOPT; double CFL=1.0/2;
const int     TYPE_RK= RK1;                     //- time scheme (ENO2 case) : 1 or ENO2RK1, 2 or ENO2RK2,  3 or ENO2RK3
//int BORDERSIZE[DIM] = {2,2,2};
int BORDERSIZE[DIM] = {3,3,3};

//----------------------------
//- stopping criteria parameters:
//----------------------------
const double  EPSILON           = 0.0;
const int     MAX_ITERATION     = 100000;

const double                   T= 2.00;         //- Terminal time

//----------------------------
//- discretization parameters:
//----------------------------
const double  DT                = 0.0;         //- if DT=0 the program will compute a DT based on space steps.
      //double  CFL               = 0.7;         //- for ENO
      //double  CFL               = 0.7;         //- for NBEE


//- for a DIM=3 test:
const int nn=1;
const int     NN=2*20;
//const int     ND[DIM]           = { NN , NN , NN };
const int     ND[DIM]           = { NN , NN , NN };
const double  XMIN[DIM]         = { -2.0, -2.0, -2.0};
const double  XMAX[DIM]         = {  2.0,  2.0,  2.0};           //- bound of the domain
const int     PERIODIC[DIM]     = {  0 ,  0 ,  0 };           //- periodic mesh (1:periodic, 0:otherwise)

//- for a DIM=4 test:
//const int     ND[DIM]           = { NN , NN , NN ,  3};
//const double  XMIN[DIM]         = { -2., -2., -2., -2};
//const double  XMAX[DIM]         = {  2.,  2.,  2.,  2};       //- bound of the domain
//const int     PERIODIC[DIM]     = {  0 ,  0 ,  0 ,  0};       //- periodic mesh (1:periodic, 0:otherwise)

const int     MESH              = 1;                          //- 0 : xi = center of cell ; 1 : xi contains the boundary

const int     cDIM              = 2;
const int     NCD[cDIM]         = { 2    , 2};
const double  UMIN[cDIM]        = { 0.0  , 0.0};
const double  UMAX[cDIM]        = { 1.0  , 1.0};

const int     cDIM2             = 1;
const int     NCD2[cDIM2]       = { 1 };
const double  UMIN2[cDIM2]      = { 0.};
const double  UMAX2[cDIM2]      = { 1.};

//-----------------------------
//- BOUNDARY : 0 (Void, FD only) or 1 (Dirichlet, V=g_border, for FD/SL) or 2 (Vx=g_bordermix, for FD), or 3 (Vxx=0, for FD/SL) 
//-----------------------------
const int BOUNDARY=0;

// Dirichlet boundary condition (case BOUNDARY=1)
//const double  VBORD             = 0.2;
//double g_border(double t, const double* x){ return VBORD; }



// Mixed Neumann bc :  vx=g(t,x,v) (case BOUNDARY=2)
double g_bordermix(double t, const double* x, double val){return 0.0;}


//-----------------------------
//- mainloop parameters
//-----------------------------
      int     COMPUTE_MAIN_LOOP = 1;
const int     COMPUTE_VEX       = 1;
const int     COMPUTE_TOPT      = 0;
const int     TOPT_TYPE         = 0;       //- 0= min time , 1= exit time

      int     SAVE_VF_ALL       = 0;       //- 1 to save "VFn.dat" every  SAVE_VF_ALL_STEP iterations
      int     SAVE_VF_ALL_STEP  = 10;      //- SHOULD BE >=1. Rem: can be set to a large number to save only VF0.dat and ev. final: VF1.dat)
      int     SAVE_VF_FINAL     = 1;       //- 1 to save last VF file as "VF.dat" 

const int     CHECK_ERROR       = 1;       //- 1 to compute errors every CHECK_ERROR_STEP iterations
const int     CHECK_ERROR_STEP  = 20;
const int     CHECK_NEG_POINTS  = 0;       //- 1/0 : to print number of negative point in domain and domain + boundary (default : 0)


//-------------------
//- "coupe" : useful for exemples with DIM>=3 (cut into some plane)
//-------------------
//const int     COUPE_DIMS[DIM] = {  1,      1,     0};
//const double  COUPE_VALS[DIM] = {0.12, 0.0234, 0.222};
//const int     COUPE_DIMS[DIM] = {  1,   0,   1};
//const double  COUPE_VALS[DIM] = {0.0, 0.15, 0.1};
//
const int     COUPE_DIMS[DIM] = {  1,   1,   0}; //- coupe 2d
const double  COUPE_VALS[DIM] = {0.0, 0.0, -0.0};
//const int     COUPE_DIMS[DIM] = {  1,   1,   1};
//const double  COUPE_VALS[DIM] = {0.0, 0.0, 0.0};
//const int     COUPE_DIMS[DIM] = {  1,      1,     0,    0};
//const double  COUPE_VALS[DIM] = {0.0, 0.0234, 0.122,  0.0};

//------------------------------------------------------------------------------------------------------
int i0=0,i1=1,i2=2;	//- direction choice for rotation pb : plane x-y
//int i0=1,i1=2,i2=0;	//- direction choice for rotation pb : plane y-z

//------------------------------------------------------------------------------------------------------
//- dynamics and distributed cost functions used in the case of COMMANDS=0 or 1 (assumes only 1 control)
//------------------------------------------------------------------------------------------------------
inline void dynamics(const double* x, C u, double t, double* res)
{
  double f0,f1,f2;
  //f0=-(-2*pi*x[1]);
  //f1=-(+2*pi*x[0]);
  //f2=       0.0;
  
  //double f0,f1,f2;
  f0=(+1.0       )*u[0];
  f1=(-2*pi*x[i2])*u[0];
  f2=(+2*pi*x[i1])*u[0];
  //f0= (-0.0)*u[0];            //- test 
  //f1= (-0.0)*u[0];
  //f2= (-1.0)*u[0];
  //f0= ( 1.0             )*u[1]; //- test
  //f1= ( 0.0 +(-1+2*u[0]))*u[1];
  //f2= (-0.0)             *u[1];
  res[i0]=f0, res[i1]=f1, res[i2]=f2;
}

inline double distributed_cost(const double* x, C u, double t){
  return 0.0;
}

//- feedback function for particular schemes like NBEE_OCOPT 
inline void feedback(double t, const double* x, const double* p, C& u){  //- only needed for NBEE-OC
  //- input  : t,x,p
  //- output : u, that maximizes  max_u (- p. dynamics(t,x,u) - dist_cost(t,x,u)) 
  //- Example
  //u[0]=0.0;
  double f0,f1,f2;
  //f1=(-2*pi*x[2]);
  //f2=(+2*pi*x[1]);
  //f0=(+1.0      );
  //f1=(-2*pi*x[i1]); //i0
  //f2=(+2*pi*x[i0]); //i1
  //f0=(+1.0       );
  //f0=(-2*pi*x[1]); //i0
  //f1=(+2*pi*x[0]); //i1
  //f2=(+1.0      );
  f0=(+1.0       );
  f1=(-2*pi*x[i2]);
  f2=(+2*pi*x[i1]);
  //f0= 0.0*u[0];
  //f1= 0.0*u[0];
  //f2= 1.0*u[0];
  u[0]=(-(f0*p[i0]+f1*p[i1]+f2*p[i2])>=0);
}


//---------------------------------------------------------------
//- if COMMANDS=2 (2 player games): dynamics and distributed cost
//---------------------------------------------------------------
inline void dynamics2(const double* x, C u, C u2, double t, double* res) { return; }
inline double distributed_cost2(const double* x, C u, C u2, double t) { return 0.; }


// Dirichlet boundary condition (case BOUNDARY=1)
//const double  VBORD             = 0.2;
//double g_border(double t, const double* x){ return VBORD; }

//---------------
//- initial data
//---------------

//double RAYON=0.5;
double RAYON=1*0.20/1;

double dx=(XMAX[0]-XMIN[0])/NN;
double resmax= RAYON*dx;
//double resmin=-RAYON*dx;
double resmin=-0.002;

const double  VBORD             = resmax;

double g_border(double t, const double* x){return VBORD; }

inline double v0(const double * xx)
{
  //- uncomment desired example:
  //- rotating circle:
  double x=xx[i0], y=xx[i1], z=xx[i2];
  double X0=1.0, Y0=0.0, Z0=-1.0; 
  double nx=sqrt( (x-X0)*(x-X0) + (y-Y0)*(y-Y0) + (z-Z0)*(z-Z0) ) - RAYON;  //- ball target
  //double nx=max(max(abs(x-X0),abs(y-Y0)),abs(z-Z0)) - RAYON; //- square target
  double factor=2*dx/10;
  return min(resmax,max(resmin,factor*nx));
  //return min(VBORD,nx-RAYON);
  //- rotating plane:
  //double x=xx[0],y=xx[1],z=xx[2];
  //return z-x+2*y;
}

//--------------------
//- OBSTACLE g
//--------------------
const int OBSTACLE=0;
inline double g_obstacle(double t, const double* x){
  //double X0=0.0, Y0=0.75, RAYON=0.25;
  //double nx=max( abs(x[0]-X0) -RAYON ,abs(x[1]-Y0) -RAYON);		//- small square obstacle
  //return MIN(VBORD,-nx);
  return resmin;
}

//---------------------------
//- Exact solution (if known)
//---------------------------

inline double Vex(double t, const double* x) {
  double y[3];
  double t1=-t;
  //y[0]=cos(2*pi*t1)*x[0] - sin(2*pi*t1)*x[1];
  //y[1]=sin(2*pi*t1)*x[0] + cos(2*pi*t1)*x[1];
  //y[2]=x[2];
  double f0=1.0;
  y[i0]=cos(2*pi*t1)*x[i0] - sin(2*pi*t1)*x[i1];
  y[i1]=sin(2*pi*t1)*x[i0] + cos(2*pi*t1)*x[i1];
  y[i2]=x[i2]-f0*t;
  return v0(y);
}


//--------------------
//- Stability constants to be used for the numerical hamiltonian.
//--------------------
inline void compute_Hconst(double* aMAX, double t)
{
  double xm=max(ABS(XMIN[0]),abs(XMAX[0]));
  double ym=max(ABS(XMIN[1]),abs(XMAX[1]));
  aMAX[i0]=2*pi*ym;
  aMAX[i1]=2*pi*xm;
  aMAX[i2]=1.0;
}

//--------------------
//- Numerical Hamiltonian :
//--------------------
inline double Hnum(const double t, const double* x, const double vi, const double* Dv)
{
  // DIRECT METHOD
  double res;
  //f0=-(-2*pi*x[1]);
  //f1=-(+2*pi*x[0]);
  //f2=       0.0;
  double f0,f1,f2;
  f1=(-2*pi*x[2]);
  f2=(+2*pi*x[1]);
  f0=(+1.0);
  //f0=0.0;
  //f1=0.0;
  //f2=1.0;
  res=  max(0.,-f0)*Dv[0] + min(0.,-f0)*Dv[1] 
      + max(0.,-f1)*Dv[2] + min(0.,-f1)*Dv[3]
      + max(0.,-f2)*Dv[4] + min(0.,-f2)*Dv[5];  
  res=max(res,0.0);
  return res;
}


//-------------------
//- compute trajectory parameters
//-------------------
// method of trajectory reconstruction and starting time
int           TRAJ_METHOD       = 0; 
const int     TRAJPT            = 0;
const double  initialpoint[TRAJPT*DIM] = {};      //- Initial point for HJB

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

//--------------------------------
//- border: number of ghost cells in each direction (default is {2,2,...} for FD)
//--------------------------------
//int BORDERSIZE[DIM] = {2,2,2};

//--------------------
//- Special parameter for method SL
//--------------------
const int     P_INTERMEDIATE    = 1;    //- number of discretisation steps for approximating each trajectory in method MSL

void init_data_general(){
  beginStep_HJB=0;     //- <Only in init_data_general(), which is before main HJB loop>  for 1 comput HJB: use beginStep_HJB=endStep_HJB=0
  endStep_HJB=0;          
}

void init_data(){

  //----------------------------------------
  //- some OUTPUTS parameter initialization:
  //----------------------------------------
  SAVE_COUPE_ALL_STEP= 1; //- SHOULD ALWAYS BE >= 1.  Can be different from SAVE_VF_ALL_STEP; used for all "coupe" files
  SAVE_COUPE_ALL     = 0;
  SAVE_COUPE_FINAL   = 1;
  //SAVE_COUPEEX_ALL   = 0;
  //SAVE_COUPEEX_FINAL = 1;

  SAVE_COUPE_INTERPOLATION=1; //- if 1 then uses interpolation. Use 0 for testing purposes (to compare VFxx.dat and coupexx.dat files)

  //---------------------------------
  //- 1ier calcul : calcul de VF.dat
  //---------------------------------

  COMPUTE_MAIN_LOOP   = 1;
  SAVE_VF_ONSET_FINAL = 0;   //- to save final value on some user defined points (such as Xuser.txt or .dat)
  //beginStep_HJB=0;         //- <Only in init_data_general(), which is before main HJB loop>  for 1 comput HJB: use beginStep_HJB=endStep_HJB=0
  //endStep_HJB=0;           //-

  //---------------------------------
  //- 2ie calcul  : calcul de XuserVREFj.dat pour j=1,2...
  //---------------------------------
/*
  COMPUTE_MAIN_LOOP   = 0;
  //EXTERNALV0          = 1;   //- Precise de charger dans v le contenu de VF.dat 
  SAVE_VF_ONSET_FINAL = 1;   //- to save final value on some user defined points (such as Xuser.txt or .dat)
  beginStep_HJB=0;
  endStep_HJB=0;             //- cas classique pour 1 calcul HJB: endStep_HJB=1
*/

  SAVE_COUPE_ALL_STEP= 20; //- SHOULD ALWAYS BE >= 1.  Can be different from SAVE_VF_ALL_STEP; used for all "coupe" files
  SAVE_COUPE_ALL     = 0;
  SAVE_COUPE_FINAL   = 1;
  SAVE_VF_ALL        = 0; //- 1 to save "VFn.dat" every  SAVE_VF_ALL_STEP iterations
  SAVE_VF_ALL_STEP   = 1; //- SHOULD BE >=1. Rem: can be set to a large number to save only VF0.dat and ev. final: VF1.dat)
  SAVE_VF_FINAL      = 1; //- 1 to save last VF file as "VF.dat" 
  BINARY      	     = 0; //- for BINARY(1) or TEXT/ASCII(0) save & load 
  FORMAT_FULLDATA    = 1; //- furthermore add (at left) the list of indexes corresponding to the current point {files VF.dat, etc.}
  SHOW_COORDINATES   = 0; //- will furthermore add (at right) the list of coordinates on each line of the {VF.dat, etc.} files 
}

