// -------------------
// - ADVANCED PARAMETERS
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

// ---------------------------
// - OTHER MAINLOOP PARAMETERS
// ---------------------------
const int     VALUE_PB           = 0; //- to compute value v(t,x) s.t. 0 = w(t,x,z)=v(t,x)-z (assumes DIM=d+1, x in R^d, z in R^1)
const int     SAVE_VALUE_ALL     = 0; //- to save at each iteration the value v(t,x) (assumes VALUE_PB=1)
const int     SAVE_VALUE_ALL_STEP= 0; //- to save final value v(t,x) (case VALUE_PB=1)
const int     SAVE_VALUE_FINAL   = 0; //- to save final value v(t,x) (case VALUE_PB=1)

//- note : {VF,Vex,topt}.dat savings use the same step as VF savings : SAVE_VF_ALL_STEP
//- note : {coupe,coupeex,coupetopt}.dat savings use the same step as coupe savings: SAVE_COUPE_ALL_STEP
      int     SAVE_VEX_ALL       = 0;
      int     SAVE_VEX_FINAL     = 1;
const int     SAVE_TOPT_ALL      = 0;
const int     SAVE_TOPT_FINAL    = 1;
//const int     SAVE_COUPE_ALL     = 0;
//const int     SAVE_COUPE_ALL_STEP= 50; //- can be different from SAVE_VF_ALL_STEP; used for all "coupe" files
      int     SAVE_COUPE_ALL_STEP= 1;   //- can be different from SAVE_VF_ALL_STEP; used for all "coupe" files
      int     SAVE_COUPE_ALL     = 1;
      int     SAVE_COUPE_FINAL   = 1;
      int     SAVE_COUPEEX_ALL   = 1;
      int     SAVE_COUPEEX_FINAL = 1;
      int     SAVE_COUPETOPT_ALL = 0;
      int     SAVE_COUPETOPT_FINAL=1;
const int     SAVE_VF_ONSET_ALL    = 0; //- savings iteration step for {VF,Vex,topt (& onset-files) is defined by "SAVE_VF_ALL_STEP" 
const int     SAVE_VF_ONSET_FINAL  = 0; //- to save final value on some user defined points (such as Xuser.txt)
const int     SAVE_VEX_ONSET_ALL   = 0;
const int     SAVE_VEX_ONSET_FINAL = 0;
const int     SAVE_TOPT_ONSET_ALL  = 0;
const int     SAVE_TOPT_ONSET_FINAL= 0;
//char          XUSER_DATA_FILE[] = "Xuser.txt"; // obsolete
//const char      XFile_user0[]   = "X_user";   // to be completed by ".txt" extention (to be distinguished from .dat files)
//const char      VFile_user0[]   = "XV_user";  // will be completed by ".dat" or "n.dat" extention 
      int     NBPOINTS_Xuser    = 0;
char          XuserFile0[]      = "Xuser";     // to be completed by ".txt" extention (to be distinguished from .dat files)
char          XuserVFile0[]     = "XuserV";    // will be completed by ".dat" or "n.dat" extention 
char          XuserVexFile0[]   = "XuserVex";  // will be completed by ".dat" or "n.dat" extention 
char          XuserToptFile0[]  = "XuserVtopt";// will be completed by ".dat" or "n.dat" extention 

const int     SAVE_COUPE_INTERPOLATION=1; //- if 1 then uses interpolation. Use 0 for testing purposes (compare VFxx.dat and coupexx.dat files)

const int     SAVE_MPI_FULL      = 1; //- save VFxx.dat (and toptyy.dat) -full matrices files- in the MPI mode, as in SEQ mode (names: VFxx_PROCyy.dat files)
                                      //- will otherwise use the same parameters as SAVE_VF** 
 
      int     BINARY             = 1; //- for BINARY(1) or TEXT/ASCII(0) save & load 

// -------------------
// - format for saving .dat files (VF.dat,VFxx.dat)  1=default(coord+val);  0=only values; see also "stdafx.h" for format
// -------------------
      int     FORMAT_FULLDATA         = 1; //- furthermore add (at left) the list of indexes corresponding to the current point {files VF.dat, etc.}
const int     SHOW_COORDINATES        = 1; //- will furthermore add (at right) the list of coordinates on each line of the {VF.dat, etc.} files 

// --------------------
// - OBSTACLE g tilde
// --------------------
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
// - For error computations
// -------------------
const double  C_THRESHOLD       = 1.0;//10.0;    //- Error threshold for error computations

// -------------------
// - For trajectory reconstruction
// -------------------
int           PRINTTRAJ         = 1;       //- xterm printings of trajectory (default should be 0)

// -------------------
// - For parameter intialization/termination
// -------------------
//void post_data(){}
//void init_data(){}

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

