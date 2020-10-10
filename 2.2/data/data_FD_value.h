//-------------------------------------------
//- Exemple  of value computation - SEPT 2011 / updated JUNE 2012 / NOV 2015
//-------------------------------------------
//- (I) CONSTRAINTS = K , minimise the value  = distributive cost + terminal cost:
//- w(t,x,z)= inf_y  max( int t^T ell(y(s),a(s))ds + phi(y(T)) - z) , max_{s in [t,T]} g(y(s)) )
//- gives the pde
//-
//-    - w_t + max_a (  (-f(x,a), ell(x,a)) . (\nabla_x u, du/dz) = 0
//-    v0(x)==u(T,x)=max ( phi(x)-z , g(x) )
//- 
//-   Reversing time : 
//-    + w_t + max_a (  (-f(x,a), ell(x,a)) . (\nabla_x u, du/dz) = 0
//-    v0(x)==u(0,x)=max ( phi(x)-z , g(x) )
//- 
//- (II) CONSTRAINTS = {stay in K} AND {reach target C}, minimise value = distributive cost :
//-      HERE terminal cost = 0, 
//-      HERE target C coded by phi(x)<=0, hence, the following approach:
//- w(t,x,z)= inf_y  MAX [  (int_t^T ell(y(s),a(s))ds - z), max_{s in [t,T]} g(y(s)), phi(y(T)) ] 
//- so that w <= 0 is equivalent to 
//-    | int_t^T ell(y(s),a(s))ds  <= z 
//-    | y(s) in K for s in [t,T]
//-    | y(T) in C
//-
//- gives the pde (since the term v0 deasappears in the formulation of the DPP)
//-
//-    min(  -w_t + max_a (  (-f(x,a), ell(x,a)) . (\nabla_x w, dw/dz) , w - g(x))  = 0
//-    w(T,x,z)=max ( -z, g(x), phi(x) )
//- 
//- Reversing time  t <--> T-t :
//-    min( + w_t + max_a (  (-f(x,a), ell(x,a)) . (\nabla_x w, dw/dz) , w - g(x))  = 0
//-    w(0,x,z)=max ( -z, g(x), phi(x) )
//- 
//- 
//- This file = approach (II)
//-
//- WE CONSIDER THE FOLLOWING CASE
//- EXEMPLE FOR PAPER BOKANOWSKI-ZIDANI-ALTAROVICI : Zermelo type pb for minimal fuel consumption
//-  Dynamics f(x,u)= (v cos(u) + c- a x_2^2, v sin(u))
//-  ell(x,a)= v
//-  phi = ball
//-  K = (square obstacle)
//-  Controls :  (u,v) in [0, 2pi] * [0,Vmax]
//-  
//----------------------------------------
//#include "data_default.h"

const char    NAME[] = "data_FD_value.h (3d) = minimal consumption pb with Zermelo dynamics";
const int     DIM               = 3;           //- space dimension

const int     COMMANDS          = 0;           //- 0: local Hnum function ; 1: Hnum defined with dynamics and distributed_cost functions
const int     OPTIM             = MAXIMUM;

//----------------------------
//- initialisation 
//----------------------------
const int     EXTERNALV0        = 0;           //- 1: for initialisation from VF.dat / tmin.dat;  0: from function v0

//----------------------------
//- DF method parameters:
//----------------------------
const int     METHOD            = MFD;         //- Method : Finite differences
const int     TYPE_SCHEME       = ENO2;        //- space scheme : 1 or LF ; 2 or ENO2
const int     TYPE_RK           = RK1;         //- time scheme (ENO2 case) : 1 or ENO2RK1, 2 or ENO2RK2,  3 or ENO2RK3

//----------------------------
//- discretization parameters:
//----------------------------
const double  DT                = 0.0;          //- if DT=0 the program will compute a DT based on space steps.
const double  CFL               = 0.7;

const int     NN                = 0;
const int     N1                = 80; 
const int     N2                = 80; 
const int     N3                = 20; 
const int     ND[DIM]           = {N1,N2,N3};

//----------------------------
//- stopping criteria parameters:
//----------------------------
const double                   T= 3.0;         //- Terminal time
const double  EPSILON           = 0.0;
const int     MAX_ITERATION     = 100000;

//----------------------------
//- domain : [Zmin, Zmax] interval should contain the expected values of V
//----------------------------
const double  Zmin    =-1.0;
const double  Zmax    = T+0.0;
const double  XMIN[DIM]         = { -3.0 , -2.0 ,  Zmin};
const double  XMAX[DIM]         = {  2.0 ,  2.0 ,  Zmax};  //- bound of the domain - last value = z value


const int     PERIODIC[DIM]     = { 0 , 0 , 0 };       //- periodic mesh (1:periodic, 0:otherwise)
const int     MESH              = 1;                   //- 0 : xi = center of cell ; 1 : xi contains the boundary
//const int     OBSTACLE          = 0;

//---------------
//- CONTROLS
//---------------
const int     cDIM              = 2;
const int     NC1               = 50;    //- number of controls 1
const int     NC2               = 50;    //- number of controls 2
const int     NC                = NC1*NC2;
const int     NCD[cDIM]         = {NC1,NC2};
const double  UMIN[cDIM]        = {0.,    0.};
const double  UMAX[cDIM]        = {2.*pi, 1.};

const int     cDIM2             = 1;
const int     NCD2[cDIM2]       = { 1 };
const double  UMIN2[cDIM2]      = { 0.};
const double  UMAX2[cDIM2]      = { 1.};

//-----------------------------
//- BOUNDARY : 0 (Void, FD only) or 1 (Dirichlet, V=g_border, for FD/SL) or 2 (Vx=g_bordermix, for FD), or 3 (Vxx=0, for FD/SL) 
//-----------------------------
const int BOUNDARY=1;

const double  VBORD             = 0.2;
const double  cutoff            = VBORD; // used in the definition of V0

// Dirichlet boundary condition (case BOUNDARY=1)
double g_border(double t, const double* x){
  return VBORD;
}

// Mixed Neumann bc :  ux=g(t,x,u) (case BOUNDARY=2)
double g_bordermix(double t, const double* x, double val){return 0.0;}


//-----------------------------
//- mainloop parameters
//-----------------------------
const int     COMPUTE_VEX       = 0;
const int     COMPUTE_TOPT      = 1;
const int     TOPT_TYPE         = 0;     // 0= min time pb, 1= exit time pb



//- For VALUE=1 (value problem) we have to set dist.cost=0. in order to compute correctly value and optimal control 
inline double distributed_cost(const double* arg, C u, double t){
  return 0.0;
}

const double Vmax=1.00;

inline void dynamics(const double* x, C u, double t, double* res){
  //double Z=x[2];
  //double f2=2.-0.5*x[1]*x[1]; 
  double f2=2.0;
  double dist_cost=1.0;
  res[0]= Vmax*u[1]*cos(u[0]) + f2;
  res[1]= Vmax*u[1]*sin(u[0]);
  res[2]= -dist_cost;
}


//- ENCOURS : dynamics2 and distributed_cost2 are unused function if only 1 control (COMMANDS=1)
inline void dynamics2(const double* x, C u, C u2, double t, double* res) {}
inline double distributed_cost2(const double* arg, C u, C u2, double t) {return 0.0;}


//const int target = 0;

inline double Hnum(const double t, const double* x, const double vi, const double* Dv) {
  double f2=2.0;
  double dist_cost=1.0;
  double z, z1, z2; 
  //double f1 = - (1.);
  //double f2 = - (1.);
  //double f3 = distributed_cost(x,NULL,t);   
  //C u;
  //double f3 = distributed_cost(x,u,t);   
  //double f3 = distributed_cost(x,(C)u,t);   
  //double f1 = - (1.);
  //double f2 = - (0.);
  //double f3 = 1.;
  //z=   MAX(0.,f1)*Dv[0] + MIN(0.,f1)*Dv[1] 
  //   + MAX(0.,f2)*Dv[2] + MIN(0.,f2)*Dv[3]  
  //   + MAX(0.,f3)*Dv[4] + MIN(0.,f3)*Dv[5];
  //z2=MAX(0.,f3)*Dv[4] + MIN(0.,f3)*Dv[5];
  //double Vmax=1.0;
  //z1=sqrt(pow(Dv[0]+Dv[1],2.) + pow(Dv[2]+Dv[3],2.)) - 0.5*( (Dv[1]-Dv[0])  +  (Dv[3]-Dv[2]) );
  z1=sqrt( pow(MAX(0.,MAX(Dv[0],-Dv[1])),2.) + pow(MAX(0.,MAX(Dv[2],-Dv[3])),2.));
  z2=Dv[4]; //- for term coming from dist. cost ell=1 (ell>=0)
  //f2=2.-0.5*x[1]*x[1];  //- f2 is positive because x[1]=x_2 is in [-2,2];
  z= Vmax*z1 + dist_cost*z2  - f2*Dv[1];
  return MAX(0.,z); 
}


const double  Hconst[DIM]={3.,1.,1.};

inline void compute_Hconst(double* cc, double t) {
  cc[0]=Hconst[0];
  cc[1]=Hconst[1];
  cc[2]=Hconst[2];
}


//- NEW : from R_general_HJ/DF_OLIVIER/datav.h

//----------------------------------------------
//- function phi corresponding to terminal cost:
//----------------------------------------------
inline double phi(const double* arg) {
    //const double  XC1=1.0, XC2=0.0, r0=0.5;
    //double x1 = arg[0], x2 = arg[1];
    //double r=sqrt((x1-XC1)*(x1-XC1) + (x2-XC2)*(x2-XC2));
    //return 3.0 - 1.0*(r-r0);
    return 0.0;
}

//----------------------------------------------
//- function psi for coding terminal constraint C (target)
//-   psi(x)<=0  <==>  x in  C
//----------------------------------------------
inline double psi(const double* arg) {
    //const double  XC1=1.0, XC2=0.0, r0=0.20;
    double XC1,XC2,r0;
    double x1=arg[0],x2=arg[1],r;
    XC1=1.5; XC2=0.0; r0=0.25;//- parametres a recopier dans complement.m
    r=sqrt((x1-XC1)*(x1-XC1)+(x2-XC2)*(x2-XC2));
    //double r=MAX(ABS(x1-XC1),ABS(x2-XC2));
    //return MIN(VBORD, r-r0);
    double C1=1.;//20.;
    return C1*(r-r0); //- truncation by VBORD will be done afterwards in v0
}

//- v0 to be defined after obstacle function
/*
inline double v0(const double* arg) {
  //double x1 = arg[0], x2 = arg[1]; 
  double x[2]; x[0] = arg[0]; x[1] = arg[1];//- parametres a recopier dans complement.m 
  double z=arg[2];
  return MIN(-z,phi(x),);
}
*/

inline double Vex(double t, const double* arg) {
  return 0;
}


//-------------------
//- OBSTACLE
//-------------------
const int OBSTACLE=1;
//inline double g_obstacle (double t , const double *x) {return 0.0;}

//- SQUARE obstacle CASE:
inline double g0(double t, double x, double y) {
  //double r=MAX(ABS(x-0.5),ABS(y-0.0));
  //double r=sqrt(pow(x-0.5,2.0)+pow(y-0.0,2.0));
  //return MAX(-VBORD, -(r-0.3));
  //return -1.;
  
  //double ra=MAX(ABS(x-( 0.0)),      ABS(y-( 0.8))); ra=0.4-ra;
  //double rb=MAX(ABS(x-(-1.0)),1./3.*ABS(y-(-0.8))); rb=0.2-rb;
  double C1=1.0;
  double r1,x1,y1,ax1,ay1;
  double r2,x2,y2,ax2,ay2;
  r1=0.4; x1=-0.5; y1= 0.5; ax1=1.0; ay1=1.0; //- parametres a recopier dans complement.m
  r2=0.2; x2=-1.0; y2=-1.5; ax2=1.0; ay2=5.0; //- parametres a recopier dans complement.m
  double ra=C1*( r1 - MAX(ABS(x-x1)/ax1 ,ABS(y-y1)/ay1) );
  double rb=C1*( r2 - MAX(ABS(x-x2)/ax2 ,ABS(y-y2)/ay2) );
  //double ra=C1* ( 0.4 - MAX( ABS(x-(-0.5))   ,  ABS(y-( 0.5))/1. ) );
  //double rb=C1* ( 0.2 - MAX( ABS(x-(-1.0))   ,  ABS(y-(-1.5))/5. ) );

  //double gmin=-0.2;
  //double res=MAX(gmin, C1*MAX(ra,rb));
  double res=MAX(ra,rb);
  return MAX(-cutoff,MIN(cutoff,res));
}

inline double g_obstacle(double t, const double* arg) {
  return g0(t,arg[0],arg[1]);
}

/*
//- LABYRINTH CASE:
inline double g0(double t, double x, double y) {
  double res=-INF;
  double dist;
  //if (max(abs_val(x-1.0),abs_val(y-0.0))<=0.2)
  //double x0=0.0, y0=0.0, r0=0.5;
  //return(float(pow(x-x0,2)+pow(y-y0,2)<=pow(r0,2)));
  //return(float(max(abs_val(x-x0),abs_val(y-y0))<=r0));
  //dx=(B0x-A0x)/Px;
  //dy=(B0y-A0y)/Py;
  //- Retirer les deux lignes suivantes si le bord est dans l'obstacle!
  //xmin_obs=x00-c_rayonobs-2*dx; xmax_obs=x00+c_rayonobs+2*dx;
  //ymin_obs=y00-c_rayonobs-2*dy; ymax_obs=y00+c_rayonobs+2*dy;

  double nr=sqrt( x*x+y*y );
  double r0=0.5, rdiff0=0.20, rouvert0=0.2;
  double r1=1.6, rdiff1=0.20, rouvert1=0.4;
  //function z=Vobstacle(x,y)
  //X=x*ones(y');
  //Y=ones(x)*(y');
  //i=find(nr>r0 & nr<r0+rdiff0 & (abs(Y)>rouvert0 |X<0)); z(i)=1;
  //i=find(nr>r1 & nr<r1+rdiff1 & (abs(X)>rouvert1 |Y>0)); z(i)=1;

  //if (nr>r0 && nr<r0+rdiff0 && (abs_val(y)>rouvert0 || x<0)) {res=1.0; return(res);}
  //if (nr>r1 && nr<r1+rdiff1 && (abs_val(x)>rouvert1 || y>0)) {res=1.0; return(res);}
  double ouv;

  //dist= -( ABS(nr-r0)-rdiff0 ); ouv=max(abs(y),abs(x-( r0)))-rouvert0;  dist=MIN( dist,  ouv); res=dist;
  //-  2 cercles ouverts que par 1 endroit:
  //dist= -( ABS(nr-r0)-rdiff0 ); ouv=max(abs(x),abs(y-(+r0)))-rouvert0;  dist=MIN( dist,  ouv); res=dist;
  //dist= -( ABS(nr-r1)-rdiff1 ); ouv=max(abs(x),abs(y-(-r1)))-rouvert1;  dist=MIN( dist,  ouv); res=MAX(res,dist);

  //-  2 cercles ouverts par 2 endroits:
  //dist= -( ABS(nr-r0)-rdiff0 ); ouv=max(abs(x),abs(y-(+r0)))-rouvert0;  dist=MIN( dist,  ouv); res=dist;
  dist= -( ABS(nr-r0)-rdiff0 ); res=dist;
  dist= -( ABS(nr-r1)-rdiff1 ); ouv=abs(x)-rouvert1;  dist=MIN( dist,  ouv); res=MAX(res,dist);

  res=MAX(-1.0,MIN(1.0,res));
  return res;
}


//- fonction obstacle
inline double g(double x1, double x2, double x3) {
    double res;
    //res=-INF; //- peut etre mettre -cutoff en general comme valeur inf.
    double t=0.;
    res=g0(t,x1,x2);
    res=MAX(res,MAX(g0(t,x1+Long*cos(x3)-Larg*sin(x3),x2+Long*sin(x3)+Larg*cos(x3)),g0(t,x1+Long*cos(x3)-Larg*sin(x3),x2+Long*sin(x3)+Larg*cos(x3))));
    res=MAX(res,MAX(g0(t,x1-Long*cos(x3)-Larg*sin(x3),x2-Long*sin(x3)+Larg*cos(x3)),g0(t,x1-Long*cos(x3)-Larg*sin(x3),x2-Long*sin(x3)+Larg*cos(x3))));
    return(res);
}

inline double g_obstacle(double t, const double* arg) {
    return g(arg[0],arg[1],arg[2]);
}
*/


inline double v0(const double* arg) {
  //double x1 = arg[0], x2 = arg[1]; 
  double x[2]; x[0] = arg[0]; x[1] = arg[1]; 
  double z=arg[2];
  //double res=MAX(-z,MAX(phi(x),g0(0.,x[0],x[1])));
  double res=MAX(MAX(phi(x)-z,psi(x)),g_obstacle(0,x)); //- psi is for adding target constraint to the value pb
  //return res;
  return MAX(-VBORD, MIN(VBORD,res)); //- truncation
  //return MAX(phi(x)-z,g0(0.,x[0],x[1]));
}


//--------------------
//- OUTPUT PARAMETERS:
//--------------------
const int     COMPUTE_MAIN_LOOP = 1;

const int     PRECOMPUTE_COORDS = 1; //- to precompute mesh coordinates (faster but needs more memory).

const int     SAVE_VF_ALL       = 1; //- 1 ==> save VFi.dat, i=0,1,2,... every SAVE_VFALL_STEP iterations 
const int     SAVE_VF_ALL_STEP  = 1;
const int     SAVE_VF_FINAL     = 1;
const int     VALUE_PB          = 1; //- 1 to compute value v / 0 = w(t,x,z)=v(t,x)-z (assumes DIM=d+1, x in R^d, z in R^1)
const int     SAVE_VALUE_ALL    = 1; //- 1 to save value v(t,x) at  all  time iterations (assumes VALUE_PB=1)
const int     SAVE_VALUE_ALL_STEP=2;
const int     SAVE_VALUE_FINAL  = 1; //- 1 to save value v(t,x) at final time (assumes VALUE_PB=1)
const int     SAVE_VF_FINAL_ONSET=0; //- to save final value on some user defined points (X_user.txt)

const int     CHECK_ERROR       = 0; //- 1 to compute errors every CHECK_ERROR_STEP iterations
const int     CHECK_ERROR_STEP  = 10;

//- "coupe" : used only for exemples with d>=3 (cut into some plane)
//const int     COUPE_DIMS[DIM] = {1 ,1 };
//const double  COUPE_VALS[DIM] = {0.,0.};
const int     COUPE_DIMS[DIM] = { 1 , 1 , 0 };
const double  COUPE_VALS[DIM] = { 0., 0., 0.};


//-------------------
//- miscellaneous parameters
//-------------------


// -------------------
// - BEGIN TEMPORARY (FOR DISTRIB VERSION)
// -------------------
//const double  VBORD=0.0; 	//- obsolete

// special parameters for SL METHOD
const int     TYPE_STA_LOOP     = NORMAL; 
const int     INTERPOLATION     = BILINEAR;
const int     ORDER             = 1;
const int     P_INTERMEDIATE    = 1;    //- number of discretisation steps for approximating each trajectory in method MSL
// special parameters for SL METHOD, SECOND ORDER (case ORDER=2)
const int     PARAMP            = 0;
inline double funcR             (const double* x, C u, double t) {return 0.;}
inline void   funcY             (const double* x, int k, double eps, C u, double t, double h, double* res) {;}
inline double discount_factor             (const double* x){return 0.0;}
// ---------------------------
// - OTHER MAINLOOP PARAMETERS
// ---------------------------
//const int     VALUE_PB                = 0; //- to compute value v(t,x) s.t. 0 = w(t,x,z)=v(t,x)-z (assumes DIM=d+1, x in R^d, z in R^1)
//const int     SAVE_FOR_VALUE_PB       = 0; //- to save at each iteration the value v(t,x) (assumes VALUE_PB=1)
//const int     SAVE_FOR_VALUE_PB_FINAL = 0; //- to save final value v(t,x) (assumes VALUE_PB=1)
// -------------------
// - format for saving .dat files (VF.dat,VFxx.dat)  1=default(coord+val);  0=only values
// -------------------
const int     FORMAT_FULLDATA         = 1;
//--------------------
//- OBSTACLE g tilde
//--------------------
const int OBSTACLE_TILDE=0;
inline double g_obstacle_tilde(double t, const double* arg)
{
  //double X1=1.0, Y1=1.0;
  //double nx=sqrt( (arg[0]-X1)*(arg[0]-X1) + (arg[1]-Y1)*(arg[1]-Y1) );
  //return MIN(VBORD,1.0-nx);
  return VBORD;
}

// --------------------
// - PRECOMPUTE_OBSTACLE (to precompute obstacle terms)
// --------------------
const int PRECOMPUTE_OBSTACLE=0;

//----------------------------------
//- Restricted computational domain:
//- if COMPUTE_IN_SUBDOMAIN==1 will compute only for x s.t. g_domain(x)<0
//----------------------------------
const int     COMPUTE_IN_SUBDOMAIN  = 0;

inline double  g_domain(const double *x){
  return g_obstacle(0.,x) - g_obstacle_tilde(0.,x);
}

// -------------------
// - END TEMPORARY 
// -------------------

// -------------------
// - For error computations
// -------------------
const double  C_THRESHOLD       = 1.0e5;     //- Error threshold for error computations

//-------------------
//- compute trajectory parameters
//-------------------
//- method of reconstruction : 0 (based on tmin) or 1 (based on value) 
//- starting time : t_TRAJ_START (starting time) may be used, in [0,T[,  if TRAJ_METHOD==1 

/*
int           TRAJ_METHOD       = 0; 
double        time_TRAJ_START   = 0.0;
// stopping criteria
double        min_TRAJ_STOP     = 0.0;     //- to stop traj reconstruction when val(x) <=min
double        max_TRAJ_STOP     = 1000.0;  //- to stop traj reconstruction when val(x) >=max
int           TARGET_STOP       = 1;       //- to stop traj reconstruction when g_target(x)<=0
inline double g_target(double t, const double* x){return psi(x);}
*/

int           TRAJ_METHOD       = 1; 
double        time_TRAJ_START   = 0.0;
// stopping criteria
double        min_TRAJ_STOP     =    0.05; //- to stop traj reconstruction when val(x) <=min
double        max_TRAJ_STOP     = 1000.00; //- to stop traj reconstruction when val(x) >=max
int           TARGET_STOP       = 0;       //- to stop traj reconstruction when g_target(x)<=0
inline double g_target(double t, const double* x){return psi(x);}


// xterm printings of trajectory (default should be 0)
int           PRINTTRAJ         = 1;


//- Initial points 
//const int TRAJPT  = 0;  const double initialpoint[TRAJPT*DIM] = {};
//const int TRAJPT  = 1;  const double initialpoint[TRAJPT*DIM] = {-1.5,1.0,3.0 };    // Initial point for HJB
//const int TRAJPT  = 1;  const double initialpoint[TRAJPT*DIM] = {-1.0,0.0,Zmax};    // Initial point for HJB
//const int TRAJPT  = 2;  const double initialpoint[TRAJPT*DIM] = {-1.0,1.0,Zmax, -2.0,-0.5,Zmax };  
const int TRAJPT  = 2;  const double initialpoint[TRAJPT*DIM] = {-1.0,1.0,Zmax, -2.0,-1.5,Zmax };  
//const int TRAJPT  = 1;  const double initialpoint[TRAJPT*DIM] = {-1.0,1.0,Zmax};



// adverse control case (only for COMMANDS=2 & TRAJ_METHOD=0)
int           ADVERSE_METHOD    = 0;
inline void u2_adverse(double t, const double* x, C& u){}


// -------------------
// - For parameter intialization/termination
// -------------------
void init_data(){};
void post_data(){};

// -------------------
// - prefix in front of the default output names 
// -------------------
char            FILE_PREFIX[]    = "";  //- prefix for the name of all data files ("" = no prefix is used.)

//--------------------------------
//- border: number of ghost cells in each direction (default is {2,2,...} for FD)
//--------------------------------
int BORDERSIZE[DIM] = {2,2,2};


