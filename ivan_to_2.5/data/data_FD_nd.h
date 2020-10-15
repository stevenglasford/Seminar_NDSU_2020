//----------------------------------------
//- Advection Equation, dimension d
//----------------------------------------
//-    u_t + < f, nabla u >  = 0;    f= 6D or more dynamics like ( -(-y,x), 0) for TESTING PURPOSES>
//-    u(0,x)=u0(x)
//----------------------------------------
#include "data_FD_nd_default.h"	// do not modify for basic examples
//----------------------------------------


const char    NAME[150]         = "data_FD_nd.h: high d example, Feb. 2014 (OB)";
//const int     DIM               = 4;   //- space dimension   

const int     COMMANDS          = 0;   //- 0: user Hnum;  1: Hnum from dynamics()/distributed_cost(); (2: for min/max pbs.)
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

const double                   T= 0.5;         //- Terminal time

//----------------------------
//- discretization parameters:
//----------------------------
const double  DT                = 0.000;       //- if DT=0 the program will compute a DT based on space steps.  
      double  CFL               = 0.5;



/*
const int     DIM=2;   //- space dimension   
const int     NN=20;
const int     ND[DIM]           = {NN,NN};
const double  XMIN[DIM]         = { -2., -2.};
const double  XMAX[DIM]         = {  2.,  2.};
const int     PERIODIC[DIM]     = {0,0};           //- periodic mesh (1:periodic, 0:otherwise)
const int     MESH              = 1;               //- 0 : xi = center of cell ; 1 : xi contains the boundary
const int     COUPE_DIMS[DIM]   = {  1,   1};
const double  COUPE_VALS[DIM]   = {0.0, 0.0};
      int     BORDERSIZE[DIM]   = {2,2};
*/


const int     DIM=3;   //- space dimension   
const int     NN=20;
const int     ND[DIM]           = {NN,NN,NN};
const double  XMIN[DIM]         = { -2., -2., -2.};
const double  XMAX[DIM]         = {  2.,  2.,  2.};
const int     PERIODIC[DIM]     = {0,0,0};           //- periodic mesh (1:periodic, 0:otherwise)
const int     MESH              = 1;                 //- 0 : xi = center of cell ; 1 : xi contains the boundary
const int     COUPE_DIMS[DIM]   = {  1,   1,   0};
const double  COUPE_VALS[DIM]   = {0.0, 0.0, 0.6};
      int     BORDERSIZE[DIM]   = {2,2,2};


/*
const int     DIM=4;   //- space dimension   
const int     NN=20;
const int     ND[DIM]           = {NN,NN,NN,NN};
const double  XMIN[DIM]         = { -2., -2., -2., -2.};
const double  XMAX[DIM]         = {  2.,  2.,  2.,  2.};
const int     PERIODIC[DIM]     = {0,0,0,0};           //- periodic mesh (1:periodic, 0:otherwise)
const int     MESH              = 1;                          //- 0 : xi = center of cell ; 1 : xi contains the boundary
const int     COUPE_DIMS[DIM]   = {  1,   1,   0,   0};
const double  COUPE_VALS[DIM]   = {0.0, 0.0, 0.0, 0.0};
      int     BORDERSIZE[DIM]   = {2,2,2,2};
*/


/*
const int     DIM=5;   //- space dimension   
const int     NN=20;
const int     ND[DIM]           = {NN,NN,NN,NN,NN};
const double  XMIN[DIM]         = { -2., -2., -2., -2., -2.};
const double  XMAX[DIM]         = {  2.,  2.,  2.,  2.,  2.};
const int     PERIODIC[DIM]     = {0,0,0,0,0};           //- periodic mesh (1:periodic, 0:otherwise)
const int     MESH              = 1;                          //- 0 : xi = center of cell ; 1 : xi contains the boundary
const int     COUPE_DIMS[DIM]   = {  1,   1,   0,   0,   0};
const double  COUPE_VALS[DIM]   = {0.0, 0.0, 0.0, 0.0, 0.0};
      int     BORDERSIZE[DIM]   = {2,2,2,2,2};
*/

/*
const int     DIM=6;   //- space dimension   
const int     NN=10;
const int     ND[DIM]           = {NN,NN,NN,NN,NN,NN};
const double  XMIN[DIM]         = { -2., -2., -2., -2., -2., -2.};
const double  XMAX[DIM]         = {  2.,  2.,  2.,  2.,  2.,  2.};
const int     PERIODIC[DIM]     = {0,0,0,0,0,0};          //- periodic mesh (1:periodic, 0:otherwise)
const int     MESH              = 1;                      //- 0 : xi = center of cell ; 1 : xi contains the boundary
const int     COUPE_DIMS[DIM]   = {  1,   1,   0,   0,   0,   0};
const double  COUPE_VALS[DIM]   = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
      int     BORDERSIZE[DIM]   = {2,2,2,2,2,2};
*/



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

double g_border(double t, const double* x){
  return VBORD;
}

// Mixed Neumann bc :  vx=g(t,x,v) (case BOUNDARY=2)
double g_bordermix(double t, const double* x, double val){return 0.0;}


//--------------------------------
//- border: number of ghost cells in each direction (default is {2,2,...} for FD)
//--------------------------------
//int BORDERSIZE[DIM] = {2,2};

//-----------------------------
//- mainloop parameters
//-----------------------------
const int     COMPUTE_MAIN_LOOP = 1;
const int     COMPUTE_VEX       = 1;
const int     COMPUTE_TOPT      = 0;
const int     TOPT_TYPE         = 0;       //- 0= min time pb , 1= exit time pb

const int     SAVE_VF_ALL       = 0;       //- 1 to save "VFn.dat" every  SAVE_VF_ALL_STEP iterations (useful for movie)
const int     SAVE_VF_ALL_STEP  = 100;     //- (a large number will save only VF0.dat and final: VF1.dat)
const int     SAVE_VF_FINAL     = 1;

const int     CHECK_ERROR       = 1;                      //- 1 to compute errors every CHECK_ERROR_STEP iterations
const int     CHECK_ERROR_STEP  = 10;
const int     CHECK_NEG_POINTS  = 0;       //- 2/1/0 : to print number of negative point in domain and domain + boundary (default : 0)


//-------------------
//- compute trajectory parameters
//-------------------
// method of trajectory reconstruction and starting time
int           TRAJ_METHOD       = 0; 
// initial points 
const int     TRAJPT            = 0;
const double  initialpoint[TRAJPT*DIM] = {};      //- Initial point for HJB
// stopping criteria
double        min_TRAJ_STOP     = 0.00;    //- to stop traj reconstruction when val(x) <=min (val=topt here)
double        max_TRAJ_STOP     = 100.00;  //- to stop traj reconstruction when val(x) >=max (val=topt here)
int           TARGET_STOP       = 1;       //- to stop traj reconstruction when g_target(x)<=0
inline double g_target(double t, const double* x){ return 0.0;}
// adverse control case (only for COMMANDS=2 & TRAJ_METHOD=0)
int           ADVERSE_METHOD    = 0;
inline void u2_adverse(double t, const double* x, C& u){}
// parameters for TRAJ_METHOD=1 : 
double        time_TRAJ_START   = 0.00;    //- Starting time, only for TRAJ_METHOD=1 (reconstruction by value)


//-------------------
//- "coupe" : (cut into some plane)
//-------------------
//const int     COUPE_DIMS[DIM] = {  1,   1,   0,   0};
//const double  COUPE_VALS[DIM] = {0.0, 0.1, 0.1, 0.0};

//------------------------------------------------------------------------------------------------------
//- dynamics and distributed cost functions used in the case of COMMANDS=0 or 1 (assumes only 1 control)
//------------------------------------------------------------------------------------------------------
inline void dynamics(const double* x, C u, double t, double* res)
{
  //- rotational dynamics here
  //double X=x[0], Y=x[1];
  //double f0,f1;
  //f0=-(-2*pi*Y);
  //f1=-(+2*pi*X);
  //res[0]=f0, res[1]=f1;
  //for (int i=2;i<DIM;i++) res[i]=0.0;
  res[0]=0.0;
}

inline double distributed_cost(const double* x, C u, double t)
{
  return 0.;
}

inline void feedback(double t, const double* x, const double* p, C& u){ 
  //- only needed for NBEE-OC
  //- input  : t,x,p
  //- output : u, that maximizes  max_u (- p. dynamics(t,x,u) - dist_cost(t,x,u)) 
  //- Example rotation (no control)
  u[0]=1.0;
}
 
//---------------------------------------------------------------
//- if COMMANDS=2 (2 player games): dynamics and distributed cost
//- unused functions if only 1 control here (COMMANDS<=1)
//---------------------------------------------------------------
inline void dynamics2(const double* x, C u, C u2, double t, double* res) { return; }
inline double distributed_cost2(const double* arg, C u, C u2, double t) { return 0.; }


//---------------
//- initial data
//---------------
double RAYON=0.5, X0=0.0, Y0=0.0;

inline double   q(double x){
  return MIN(VBORD,x-RAYON);
}

inline double v0(const double * x){
  //double X0=1.0, Y0=0.0, RAYON=0.5;
  //double nx=sqrt( (x[0]-X0)*(x[0]-X0) + (x[1]-Y0)*(x[1]-Y0) );
  //double X0=1.0, RAYON=0.5;
  //return MIN(VBORD,nx-RAYON);
  //- smooth ball target
  double r02=RAYON*RAYON; 
  double nx2=0.0; for (int d=0;d<DIM;d++) nx2=nx2+x[d]*x[d];
  return VBORD*(1-pow(max(0.0, (1-nx2)/(1-r02)),4));
}


//---------------------------
//- Exact solution (if known)
//---------------------------
/*
//- Rotation case (2d)
inline double Vex(double t, const double* x) {
  double y[3];
  double t1=-t;
  y[0]=cos(2*pi*t1)*x[0] - sin(2*pi*t1)*x[1];
  y[1]=sin(2*pi*t1)*x[0] + cos(2*pi*t1)*x[1];
  y[2]=0.0;
  return v0(y);
}*/
//- Eikonal case
inline double Vex(double t, const double* x)
{
  /*
  double nx2=0.0; for (int d=0;d<DIM;d++) nx2=nx2+x[d]*x[d];
  double nx=sqrt(nx2);
  if      (OPTIM==MAXIMUM) return q( MAX(0, nx-t) );
  else if (OPTIM==MINIMUM) return q( MAX(0, nx+t) );
  else {printf("this OPTIM not programmed!\n"); return q( MAX(0, nx-t) );}
  */
  //- smooth ball target
  double r02=RAYON*RAYON; 
  double nx2=0.0; for (int d=0;d<DIM;d++) nx2=nx2+x[d]*x[d]; 
  double nx=sqrt(nx2);
  nx=max(nx-t,0.0); nx2=nx*nx;
  return VBORD*(1-pow(max(0.0, (1-nx2)/(1-r02)),4));
  //return VBORD*(1-pow(max(0.0, (1-max(nx2-t,0.0))/(1-r02)),4));
}


//--------------------
//- Stability constants to be used for the numerical hamiltonian.
//--------------------
inline void compute_Hconst(double* aMAX, double t)
{
  //double xm=max(ABS(XMIN[0]),abs(XMAX[0]));
  //double ym=max(ABS(XMIN[1]),abs(XMAX[1]));
  //aMAX[0]=2*pi*ym;
  //aMAX[1]=2*pi*xm;
  //for (int i=2;i<DIM;i++) aMAX[i]=0.01;
  for (int i=0;i<DIM;i++) aMAX[i]=1.0;
}

//--------------------
//- Numerical Hamiltonian :
//--------------------
inline double Hnum(const double t, const double* x, const double vi, const double* Dv)
{
  //- Rotation
  //double res;
  //double X=x[0], Y=x[1]; 
  //double f0,f1;
  //f0=-(-2*pi*Y);
  //f1=-(+2*pi*X);
  //res=  max(0.,-f0)*v[0] + min(0.,-f0)*v[1] 
  //    + max(0.,-f1)*v[2] + min(0.,-f1)*v[3];
  //- Eikonale 
  double res=0.0;
  double p;
  for(int d=0;d<DIM;d++) {
    p=max(0.0,max(Dv[2*d],-Dv[2*d+1]));
    res += p*p;
  }
  res=sqrt(res);
  return res;
}



//-----------------------------------------------
//- MORE ADVANCED PARAMETERS
//-----------------------------------------------

//--------------------
//- OBSTACLE g
//--------------------
const int OBSTACLE=0;
inline double g_obstacle(double t, const double* x){ return 0.;}

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
  SAVE_COUPE_ALL_STEP= 20; //- can be different from SAVE_VF_ALL_STEP; used for all "coupe" files
  SAVE_COUPE_FINAL   = 1;
  SAVE_COUPEEX_ALL   = 0;
  SAVE_COUPEEX_FINAL = 1;
}

