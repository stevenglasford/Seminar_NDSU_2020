#ifndef _HJB_H
#define _HJB_H

#include "stdafx.h"

#include "RegularMesh.h"

#include "Commands.h"


//- Root class which contains mesh and commands parameters, interpolation and output functions, supports each single method (DF,SL)
class HJB {

protected:

  // -------------------
  // - output filenames
  // -------------------
  char      VFile      [100];
  char      XuserFile    [100];
  char      XuserVFile   [100];
  char      XuserVexFile [100];
  char      XuserToptFile[100];
  char      VEXFile    [100];
  char      VFALLFile  [100];
  char      tminFile   [100];
  char      dataFile   [100];
  char      ValueFile  [100];
  char      coupeFile  [100];
  char      coupeexFile[100];
  char      trajFile   [100];
  char      trajFileMore[100];
  char      DtFile     [100];
  char      extern_VFile[100];
  char      extern_tminFile[100];
  char      successTrajectoriesFile [100];
  char*     filePrefix; 
  char*     filePrefixExt; 

  char      toptcoupeFile[100]; // [DEC 2016]

  //const char      VTKFile[]   = "OUTPUT/res";
  char      dataName[100];


  int       type_scheme;        /*! Used method (between DF, SL, FMM and UB) */
  int       dim;                /*! problem dimension */
  int       OPTIM_PARAM;        /*! MAX or MIN or MINMAX or MAXMIN */
  int       epsOPTIM;		//- +1 if MAX(or MAXMIN) or -1 if MIN(or MINMAX) for OPTIM_PARAM
  Commands  *u;                 /*! Commmands */
  Commands  *u2;                /*! Commmands */

  int       ncall;              /*! Commmands number */
  int       ncall2; 		/*! in case of second set of commands */


  RegularMesh   *mesh;          /*! Mesh on which we are working */

  double    *divdx;             /*! Computational optimization for division 1/Dx[d] */
  double    *lb,*hb;
  double    spacesteps;

  //double    TYPE_VBORD;         /*! Border values' type (VBORD vlue, VBORD function or periodic mesh) */
  int       BOUNDARY;           //- type of boundary treatment

  double    *v;                 /*! Value function [SIZE] */
  double    *vold;              /*! Used to stocked the previous iteration's value function [SIZE] */
  double    *topt;              /*! Minimal time value function [SIZE] */
  double    *v_obstacle;        //- may2014: g_obstacle       precomputation
  double    *v_obstacle_tilde;  //- may2014: g_obstacle_tilde precomputation

  int       *rank;              /*! Inner mesh points' ranks in augmentated mesh [SIZE] */
                                /*! iterator of computerank */
  int       rkfirst;            /*! Mesh point's rank with integer coordinates at [0,0,..,0] in augmentated mesh */
  int       ranksize;           /*! size of rank (default is mesh->inn_nbPoints) */

  int       *dvectint;          /*! Local int     pointer */
  double    *dvectdouble;       /*! Local double  pointer */
  int       *iCoord;            /*! Local int     pointer for function 'mesh::setranks' output */
  double    *rCoord;            /*! Local double  pointer for function 'mesh::setrcoords' output */

  int       iteration;          /*! Main loop function iteration */
  double    t;                  /*! Current time */
  double    Dt;                 /*! Time step */
  int       OBSTACLE;           /*! Test if obstacle is included in the problem */
  int       OBSTACLE_TILDE;     /*! Test if obstacle(superior) is included in the problem */
  int       PRECOMPUTE_OBSTACLE;/*!  may2014: Precompute v_obstacle and v_obstacle_tilde */
  int       COMPUTE_IN_SUBDOMAIN;/*! 0 or 1 for subdomain computations */
  int       compute_main_loop;  /*! Test if algorithm mainloop is computed */
  int       compute_vex;        /*! Test if exact solution of the current problem is computed */
  int       compute_topt;       /*! Test if minimal/maximal time problem is studied (AD:) */
  int       TOPT_TYPE;          /*! Type of tmin/tmax time problem (0= tmin, 1= tmax of exit time)*/
  int       save_vfall;         /*! savings for v        - during the mainloop computation (data: SAVE_VF_ALL) */
  int       savevfallstep;      /*! savings for v : step - iteration step (data: SAVE_VF_ALL_STEP) */
  int       savevffinal;        /*! savings for v final  - (data: SAVE_VF_FINAL) */
  // 2018 
  int       savevffinalonset;   /*! save final value */
  int       VALUE_PB;           //- savings for value pb
  int       SAVE_VALUE_ALL;     //- savings for value pb
  int       SAVE_VALUE_ALL_STEP;//- savings for value pb
  int       SAVE_VALUE_FINAL;   //- savings for value pb
  int       format_fulldata;
  int       check_error;        /*! Test if, during the mainloop computation, errors evaluation needs to be called */
  int       checkerrorstep;     /*! Iteration step for errors evaluation */
  double    c_threshold;        /*! threshold for errors */
  int       computetraj;        /*! Test if optimal trajectories need to be computed in the post processing */
  int       TRAJ_METHOD;        //  reconstruction with topt (0) or value (1)
  double    time_TRAJ_START;    //  starting time for trajectory reconstruction
  double    min_TRAJ_STOP;      //  stopping threshold for trajectory reconstruction: to stop if val(x)<=min (rather than of val(x)<=0)
  double    max_TRAJ_STOP;      //  stopping threshold for trajectory reconstruction: to stop if val(x)>=max (rather than of val(x)>=INF)
  int       TARGET_STOP;        //  stopping on target (0/1) 
  int       ADVERSE_METHOD;     //  type of adverse control (0/1)
  int       PRINTTRAJ;          //  for further trajectory output printings 

  int       SAVE_VEX_ALL;        /*! savings for vex */  
  int       SAVE_VEX_FINAL;      /*! savings for vex final */
  int       SAVE_TOPT_ALL;       /*! savings for topt */
  int       SAVE_TOPT_FINAL;     /*! savings for topt final */
  int       SAVE_COUPE_ALL;      /*! savings for coupe */
  int       SAVE_COUPE_ALL_STEP; /*! savings for coupe : step */
  int       SAVE_COUPE_FINAL;    /*! savings for coupe final */
  int       SAVE_COUPEEX_ALL;    /*! savings for coupeex */
  int       SAVE_COUPEEX_FINAL;  /*! savings for coupeex final */
  int       SAVE_COUPETOPT_ALL;  /*! savings for coupetopt */
  int       SAVE_COUPETOPT_FINAL;/*! savings for coupetopt final */
  int       SAVE_VF_ONSET_ALL;    /*! savings for v - mesh user defined  */  
  int       SAVE_VF_ONSET_FINAL;  /*! savings for v final - mesh user defined */
  int       SAVE_VEX_ONSET_ALL;   /*! savings for vex - mesh user defined  */  
  int       SAVE_VEX_ONSET_FINAL; /*! savings for vex final - mesh user defined */
  int       SAVE_TOPT_ONSET_ALL;  /*! savings for topt - mesh user defined  */  
  int       SAVE_TOPT_ONSET_FINAL;/*! savings for topt final - mesh user defined */
  int       NBPOINTS_Xuser;
  char*     XuserFile0;
  char*     XuserVFile0;
  char*     XuserVexFile0;
  char*     XuserToptFile0;

  int       SAVE_COUPE_INTERPOLATION; /*! savings for coupe/coupeex : 1 ==> interpolates on grid */
  int       SHOW_COORDINATES;   /*! to show coordinates at right of the VF.dat (etc.) files */
  int       BINARY;             /*! for BINARY(1) or TEXT(0) save & load */

  //int       savevffinalonset;   /*! save final value on some user defined grid set */
  int       VALUE_PB;           //- savings for value pb
  int       SAVE_VALUE_ALL;     //- savings for value pb
  int       SAVE_VALUE_ALL_STEP;//- savings for value pb
  int       SAVE_VALUE_FINAL;   //- savings for value pb
  int       format_fulldata;
  int       check_error;        /*! Test if, during the mainloop computation, errors evaluation needs to be called */
  int       checkerrorstep;     /*! Iteration step for errors evaluation */
  int       CHECK_NEG_POINTS;   /*! Iteration step for errors evaluation */
  double    c_threshold;        /*! threshold for errors */
  int       computetraj;        /*! Test if optimal trajectories need to be computed in the post processing */
  int       TRAJ_METHOD;        //  reconstruction with topt (0) or value (1)
  double    time_TRAJ_START;    //  starting time for trajectory reconstruction
  double    min_TRAJ_STOP;      //  stopping threshold for trajectory reconstruction: to stop if val(x)<=min (rather than of val(x)<=0)
  double    max_TRAJ_STOP;      //  stopping threshold for trajectory reconstruction: to stop if val(x)>=max (rather than of val(x)>=INF)
  int       TARGET_STOP;        //  stopping on target (0/1) 
  int       ADVERSE_METHOD;     //  type of adverse control (0/1)
  int       ODE_SCHEME;
  int       CONTROL_METHOD;
  int       PRINTTRAJ;          //  for further trajectory output printings 


  //double    diffvmax;           /*! Stop criterion for value function's variation between 2 iterations */
  double    epsilon;
  double    T;                  /*! Stop criterion for time */
  int       MAX_ITERATION;      /*! Stop criterion for iteration number */
  double    *initialpoint;
  int       *coupe_dims;        /*! dimensions of the cut */
  double    *coupe_vals;        /*! values of the cut */
  double    norm1, normi;       /*! error values : norm one & infinite norm */

  int       *dvint;             /*! for SL only : local int* pointer */

  INT       *dvectint;          /*! MPI : local int* pointer */
  INT       OMP_ENABLE;         /*! Test if OpenMP parallelization is used */
  INT       OMP_NUM_THREADS;    /*! Thread number used in the MPI */
  INT       MPI_ENABLE;         /*! set to 0 when no MPI, or MPI on but only one process */
  INT       my_rank;            /*! MPI rank of current CPU in the MPI communicator intracomm */

  int       periodic_mesh;       /*! 1:if at least one of the dimension is set periodic, 0:otherwise */
  int       *periodsize;         /*! number of points for periodic mesh for each dimension */
  int       **period;            /*! indexes of the shared points for each dimension in the periodic mesh */
                                 /*! not used for MPI (only for MPI_ON==1 and MPI_ENABLE==0 for only 1 process)*/

  int       til;                 /*!  for bilinear interpolation [LONG] */
  double    *tabInter;
  int       *ddd;
  double    *P;
  int       *icorner, *cc;
  int       *p;
  int       iter;
  int       nbSaves;             /*! number of saved "VF" files  */
  int       nbSavesCoupe;        /*! number of saved "coupe" files  */
  int       **bordersData;

  vector<double> tictoc_stack; // for tic, toc, and gettimeofday use
  struct timeval s_time;


private:
  /*!
   *  \brief copy Constructor
   */
  HJB(const HJB &);

  /*!
   *  \brief copy assignment
   */
  void operator=(const HJB &);


protected:
  /*! constructeurs et destructeurs
   */

  /*!
   *  \brief default Constructor
   */
  HJB(){};


  /*!
   *  \brief alternate Constructor
   *
   *  \param RegularMesh_params : includes all regular mesh parameters
   *  \param commande_params    : includes all commands/controls parameters
   *  \param commande_params2   : includes all second commands/controls parameters
   *  \param parallel_params    : includes all parallelization parameters
   *  \param general_params     : includes all general parameters for mainloop
   *  \param postprocess_params : includes all post processing parameters
   */
  HJB(struct mrp RegularMesh_params, struct cp commande_params, struct cp commande_params2, struct pp parallel_params, struct gp general_params, struct ppp postprocess_params);

  /*!
   *  \brief default Destructor
   */
  virtual ~HJB();


public:

  /*!
   *  \brief execute the main algorithm of the chosen method
   */
  virtual void mainloop                 ()  = 0;

  /*!
   *  \brief execute init/post processing after computing the final value function or time minimal problem results
   */
  void      initprocess                 (); // TODO 2018
  void      postprocess                 ();
  //void      postprocess_user            (); // 2018

  /*!
   * \brief  implemented in user_{init,post}process(_default).h and is void by default; user can complete it with personal code
   */
  void      user_initprocess            (); // TODO 2018
  void      user_postprocess            ();
  void      user_stepprocess            (double t, double dt, double *vin, double *vout); 


  /*!
   *  \brief matlab tic/toc function to simply display on the console how long time took the program to do its processing.
   */
  void      tic                         ();
  double    toc                         (bool);

  /*!
   *  \brief MIN or MAX for optimisation mainloop
   *  \param a : first  double (to be compared with parameter b)
   *  \param b : second double (to be compared with parameter a)
   */
  double    (*OPTIM)                    (double a, double b);


  /*!
   *  \brief load in table v the value function saved in the string parameter
   *  \param file : name of the file to be readed
   */
  void      loadV0                      (const char* file, int PRINT);  // slow(old)
  void      loadV1                      (const char* file, int PRINT);  // fast(old)
  void      loadV                       (const char* file, int PRINT);  // old
  void      loadVlocal                  (const char* file, int PRINT);  // old

  void      loadVp                      (const char* filePrefix, const char* file0, double* vtab, int PRINT); // 2018
  void      loadVplocal                 (const char* filePrefix, double* vtab, int PRINT); // 2018
  void      (HJB::*loadVp)              (const char* filePrefix, const char* file0, double* vtab, int PRINT); // ENCOURS 2018
  void      loadVpTEX                   (const char* filePrefix, const char* file0, double* vtab, int PRINT); // ENCOURS 2018
  void      loadVpBIN                   (const char* filePrefix, const char* file0, double* vtab, int PRINT); // ENCOURS 2018
  void      loadVplocalTEX              (const char* filePrefix, double* vtab, int PRINT);                    // ENCOURS 2018
  void      loadVplocalBIN              (const char* filePrefix, double* vtab, int PRINT);                    // ENCOURS 2018

  /*!
   *  \brief save table v the value function in the string parameter
   *  \param file : name of the file to be written
   */
  void      saveV                       (const char* file);	       // to be removed
  void      saveVp                      (const char* file, double* vtab, int PRINT); 
  void      (HJB::*saveVp)              (const char* file, double* vtab, int PRINT); // 2018 ENCOURS
  void      saveVpTEX                   (const char* file, double* vtab, int PRINT); // 2018 ENCOURS
  void      saveVpBIN                   (const char* file, double* vtab, int PRINT); // 2018 ENCOURS

  void      saveVpOnSet                 (const char* file, double* vtab, int PRINT); // v or topt 
  void      (HJB::*saveVpOnSet)         (const char* file, double t, double* vtab, int PRINT); // 2018 ENCOURS
  void      saveVpOnSetTEX              (const char* file, double t, double* vtab, int PRINT); // 2018 ENCOURS
  void      saveVpOnSetBIN              (const char* file, double t, double* vtab, int PRINT); // 2018 ENCOURS

  //void      saveVexpOnSet               (const char* file, double t, int PRINT);     // vex 
  void      (HJB::*saveVexpOnSet)       (const char* file, double t, int PRINT); // 2018 ENCOURS
  void      saveVexpOnSetTEX            (const char* file, double t, int PRINT); // 2018 ENCOURS
  void      saveVexpOnSetBIN            (const char* file, double t, int PRINT); // 2018 ENCOURS

  /*!
   *  \brief load in table tmin the minimal time solution saved in the string parameter
   *  \param file : name of the file to be read
   */
  void      loadTmin                    (const char* file);            // old
  void      loadTminp0                  (const char* file, int PRINT); // old
  void      loadTminp                   (const char* file, int PRINT); // new - obsolete 2018
  void      loadTminlocal               (const char* file, int PRINT); //     - obsolete 2018

  /*!
   *  \brief save table tmin the minimal time solution in the string parameter
   *  \param file : name of the file to be written
   */
  void      saveTmin                    (const char* file);
  void      saveTminp                   (const char* file, int PRINT);


  /*!
   *  \brief save values {VF,Vex,topt} on hard disk - a priori linked to savetabs 
   *  \brief 2018 obsolete
   */
  void      (HJB::*savevfall)           (double tloc, int it, int PRINT); 
  //void      (HJB::*savevfallcoupe)      (double tloc, int it, int PRINT); 

  /*!
   *  \brief saving of {VF,Vex,topt}.dat on hard disk (SEQ/MPI)
   *  \brief saving of full matrices under MPI (if SAVE_MPI_FULL=1)
   *  \param  tloc : time (used for Vex savings)
   *  \param    it : current iteration in the mainloop algorithm
   *  \param PRINT : if 1 then prints a one-line comment for each saved file
   */
  void      savetabs                    (double tloc, int it, int PRINT);

/*!
   *  \brief in the MPI parallelization case (domain decomposition)
   *  \brief save table value function (ex:VF.dat) or minimal time (ex:tmin.dat) or Vex (ex:VEX.dat) in the file named after string parameter
   *  \brief every MPI procs sends its data to the master except itself
   *  \param iteration : current iteration in the mainloop algorithm
   */
  void      savetabsMPI0                (int it); //- memory limited
  void      savetabsMPI1                (int it); //- slow
  void      savetabsMPI                 (int it); //- last
  void      savetabsMPIIO               (int it); 	  //- 2014
  void      savetabsMPIIOlocal          (char*, double*); //- 2014
  void      savetabsMPIunordered        (int it); 	  //- 2014

  /*!
   *  \brief in the MPI parallelization case (domain decomposition)
   *  \brief savefulltabs sub function for saving table value function or minimal time in the file named after string parameter
   *  \param file : name of the file to be written
   *  \param tab  : table data to be saved
   *  \param tref : test if it is the table value function or minimal time which is to be saved
   */
  void      saveTab                     (const char* file, const double* tab , int tref);
 

  /*! \brief savings of coupe{it}.dat, coupeex{it}.dat, coupetopt{it}.dat on hard disk
   */
  void      savetabscoupe               (double tloc, int it, int PRINT);

  /*!
   *  \brief in the MPI parallelization case (domain decomposition)
   *  \brief load table value function or minimal time from the file named after string parameter
   *  \param file : name of the file to be written
   *  \param tab  : table data to be saved
   *  \param tref : test if it is the table value function or minimal time which is to be saved
   */
  void      loaDtab                     (const char* file, double* tab , int tref);


  /*!
   *  \brief Counting number of v(i)<=0 inside / total 
   */
  void      comptage_basic              (int &n1, int &n2, double* vin, int PRINT, int PRINTminmax);

  void      comptage                    (int &n1, int &n2, double* vin, double xmin[], double xmax[], double& vmini, int COMPUTEminmax, int PRINT, int PRINTminmax);

  /*!  \brief norms |v-vold| : Li,L1, and max(v)  */
  void      vdiff                       (double *v, double *vold, double &nLi, double &nL1, double &maxi);


  /*!  \brief get size of file in bites */
  long      getFileSize                 (const char* filename);
 
  /*!  \brief convert time (sec) into years, months, ans days */
  void convert_years_months_days
    (double t_seconds, int &ty, int &tm, int &td, int PRINT, string mess_beg, string mess_end);


protected:

  void      checkDiffMax                (const double*);

  /*!
   *  \brief copy the first argument table (usually the value function table itself) to the second argument table
   *
   *  \param vin  : table to be copied (usually the value function)
   *  \param vout : table to backup the fist argument table
   */
  void      (HJB::*copyvf)              (double *vin, double *vout);
  void      COPY_VF                     (double *vin, double *vout);  // is omp
  //  void      COPY_VF_omp                 (double *vin, double *vout);


  /*!
   *  \brief compute minimal time problem
   *  \param t : current time in the mainloop algorithm
   */
  void      (HJB::*compTOPT)            (double t); // tmax or tmin
  void      computeTmin                 (double t); // tmin time, omp
  //void      computeTmin_omp             (double t); 
  void      computeTmax                 (double t); // tmax time, omp  (=exit time)
  //void      computeTmax_omp             (double t); // exit time

  /*!
   *  \brief compute obstacle for the value function
   *  \param t : current time in the mainloop algorithm
   */
  void      (HJB::*compOBSTACLE)        (double t);
  void      computeObstacle             (double t);
  void      computeObstacle_omp         (double t);



  /*!
   *  \brief compute Hausdorff distance between two fronts
   *  \param front1 : first front's points coordinate table
   *  \param lfront1: points' number of the first front
   *  \param front2 : second front's points coordinate table
   *  \param lfront2: points' number of the second front
   *  \param d1     : result distance between first and second front
   *  \param d2     : result distance between second and first front
   *  \param z      : max(d1,d2)
   */
  void      Hausdorff                   (const double** front1, int lfront1, const double** front2, int lfront2, double& d1, double& d2, double& z);


  /*!
   *  \brief compute norm L1 and norm L infinite
   *  \param t : current time in the mainloop algorithm
   */
  void      (HJB::*compError)           (double t);
  void      computeError0               (double t);
  void      computeError                (double t);
  void      computeErrorMPI             (double t);


  /*!
   *  \brief stopping criterion, test if max_i (abs(V^(n)[i]-V^(n-1)[i])) is smaller than the epsilon threshold
   *  \param t : current time in the mainloop algorithm
   */
  bool      (HJB::*testStoppingCriteriaEps) (double t);
  bool      testStoppingCriteriaEpsilon (double t);
  bool      testStoppingCriteriaEpsilonMPI(double t);

  /*!
   *  \brief save mesh parameters of the problem
   *  \param file  : name of the file to be written
   */
  void      saveData                    (const char* file) const;
  /*!
   *  \brief save table Vex with exact value function in the string parameter
   *  \param file  : name of the file to be written
   *  \param t     : current time in the mainloop algorithm
   */
  void      saveVEX                     (const char* file, double t);
  //void      saveVexp                    (const char* file, double t, int PRINT);
  void      (HJB::*saveVexp)            (const char* file, double t, int PRINT); // 2018 ENCOURS
  void      saveVexpTEX                 (const char* file, double t, int PRINT); // 2018 ENCOURS
  void      saveVexpBIN                 (const char* file, double t, int PRINT); // 2018 ENCOURS


  /*!
   *  \brief save table value function or minimal time in the VTK format in the file named after string parameter
   *  \param file  : name of the file to be written
   *  \param tab   : data table to be saved
   */
  void      saveVTK                     (const char* file, int it, const double* tab);

  /*!
   *  \brief save in VTK format the 2D environment with satarting area / target / moving obstacles (depending of the time t)
   *  \param file  : name of the file to be written
   *  \param it    : current iteration in the mainloop algorithm
   *  \param t     : current time in the mainloop algorithm
   */
  void      saveVTK2D                   (const char* file, int it, double);

  //- ENCOURS
  void      saveVTK6                    (const char* file, int it, double t);

  /*! #VERSION2
   *  \brief save in VTK format the 2D vehicle (modelized with a rectangle) trajectory
   *  \param file  : name of the file to be written
   *  \param it    : VTK rank to start the vehicle trajectory animation
   *  \param rk    : computed points' number in the discretization of the optimal trajectory
   *  \param steps : steps' number between two trajectory's trace
   *  \param ot    : computed points in the discretization of the optimal trajectory
   */
  void      saveTrajCAR                 (const char* file, int it, int rk, int steps, double** ot);

  /*!
   *  \brief 2-dimensional savings for the value (datav.h)
   *  \param file  : name of the file to be written
   */
  void      saveV3d                     (const char* file);

  /*!
   *  \brief d-dimensional saving for the value (datav.h)
   *  \param file  : name of the file to be written
   */
  //- ENCOURS
  //void      saveVnd                     (const char* file);
  void      saveVnd                     (const char* file, double* vtab, int PRINT);
  int       computeVnd                  (const double* x,  double &y, const double* vtab);

  /*! #VERSION2
   *  \brief save a 2d projection of the value function
   *  \param file  : name of the file to be written
   *  \param coupe : axes on which the projection is made
   *  \param val   : discretization value of cutted axes to make the projection
   */
  //void      saveCoupeM                  (const char* file, const int* coupe, const double* val, double* vtab, int PRINT);
  //void      saveCoupeMex                (const char* file, const int* coupe, const double* val, int PRINT);
  void      (HJB::*saveCoupeM)          (const char* file, const int* coupe, const double* val, double* vtab, int PRINT);
  void      saveCoupeMTEX               (const char* file, const int* coupe, const double* val, double* vtab, int PRINT);
  void      saveCoupeMBIN               (const char* file, const int* coupe, const double* val, double* vtab, int PRINT);
  void      (HJB::*saveCoupeMex)        (const char* file, const int* coupe, const double* val, int PRINT);
  void      saveCoupeMexTEX             (const char* file, const int* coupe, const double* val, int PRINT);
  void      saveCoupeMexBIN             (const char* file, const int* coupe, const double* val, int PRINT);

  /*!
   *  \brief compute borders points in the augmented mesh
   */
  /*
  void      (HJB::*ComputeBorders)      (double t);
  void      ComputeBordersVoid          (double t); //- do nothing - default value for ComputeBorder
                                                    //-    : this will keep the initial value
  void      ComputeBordersDirichlet     (double t); //- based on the gborder(.) function
  void      ComputeBordersVxx           (double t); //- v_xi_xi=0 condition
  */
  //void      ComputeBordersSEQ         ();
  //void      ComputeBordersMPI         ();
  //void      ComputeBorders            ();
  //- FOR FD
  void      (HJB::*ComputeBorders)      (double t, double*);
  void      ComputeBordersVoid          (double t, double*); //- do nothing - will keep the initial value
  void      ComputeBordersDirichlet     (double t, double*); //- based on the gborder(.) function
  void      ComputeBordersVxx           (double t, double*); //- v_xi_xi=0 boundary condition
  void      ComputeBordersMix           (double t, double*); //- v_xi=F(t,x,v) mixed boundary condition

  double    (*g_border)                 (double, const double*);  //- boundary condition function for FD and/or SL
  double    (*g_bordermix)              (double, const double*, double);  //- b.c. for FD and/or SL

  //- FOR FD/SL : 
  double    (HJB::*TbordCompute)        (double t, const double* x, const double* vtab); 
  double    TbordCompute_tmin           (double t, const double* x, const double* vtab); 

  //- FOR SL [OR FD FOR SAVE (OnSet), etc.]
  double    (HJB::*VbordCompute)        (double t, const double* x, const double* vtab);  //- for SL only
  double    VbordComputeVoid            (double t, const double* x, const double* vtab);  //- for SL - this fonction should not be called.
  double    VbordComputeDirichlet       (double t, const double* x, const double* vtab);  //- for SL
  double    VbordComputeVx              (double t, const double* x, const double* vtab);  //- for SL
  double    VbordComputeVxx             (double t, const double* x, const double* vtab);  //- for SL

  void      ComputePeriodic             (const double*, double*);

  void      (*dynamics)                 (const double* x, C u, double t, double* res);
  void      (*dynamics2)                (const double* x, C u, C u2, double t, double* res);
  void      (*feedback)                 (double t, const double* x, const double* p, C& u);
  double    (*distributed_cost)         (const double*, C,    double);
  double    (*distributed_cost2)        (const double*, C, C, double);
  double    (*v0)                       (const double*);          //- initial data
  double    (*Vex)                      (double, const double*);
  double    (*g_obstacle)               (double, const double*);  //-
  double    (*g_obstacle_tilde)         (double, const double*);  //-
  double    (*g_domain)                 (const double*);          //- g_domain(x)<0 where we compute
  double    (*g_target)                 (double, const double*);  //- g_target: for traj. reconstruction, to stop if g_target(x)<=0.
  void      (*u2_adverse)               (double, const double*, C&);    //- adverse control - user defined
  //void      (*u2_adverse)               (double, const double*, double*);    //- adverse control - user defined


  /*!
   *   PMP link: adjoint dynamics
   */ 
  //void      (*dynamicsGrad)             (const double*, C, double, double**);
  void      (*dynamicsGrad)             (const double* x, C u, double, double** res);

  /*!
   *  \brief function pointer where one interpolation function is linked with
   */
  double    (HJB::*interpolation)       (const double*, const double*);


  /*!
   *  \brief funtion is independant of the mesh dimension
   *  \brief compute the bilinear interpolation for a given point with designed results table (value function or minimal time)
   *  \param PP  : current mesh point coordinates
   *  \param tab : table data for computation
   */
  double    bilinearInterpolation       (const double* x, const double* tab);

  /*! 
   *  \brief funtion is implemented only for 2d meshes
   *  \brief Works in critical cases where some mesh size are 1 in some directions
   */
  double    bilinearInterpolation_test  (const double* PP, const double* tab);

  /*!
   *  \brief funtion is implemented only for 2d meshes
   *  \brief compute the bilinear interpolation for a given point with designed results table (value function or minimal time)
   *  \param PP  : current mesh point coordinates
   *  \param tab : table data for computation
   */


  double    bilinearInterpolation1d   (const double* PP,const double* tab);
  double    bilinearInterpolation2d   (const double* PP,const double* tab);
  /*!
   *  \brief funtion is implemented only for 3d meshes
   *  \brief compute the bilinear interpolation for a given point with designed results table (value function or minimal time)
   *  \param PP  : current mesh point coordinates
   *  \param tab : table data for computation
   */
  double    bilinearInterpolation3d   (const double* PP,const double* tab);
  /*!
   *  \brief funtion is implemented only for 4d meshes
   *  \brief compute the bilinear interpolation for a given point with designed results table (value function or minimal time)
   *  \param PP  : current mesh point coordinates
   *  \param tab : table data for computation
   */
  double    bilinearInterpolation4d   (const double* PP,const double* tab);
  /*!
   *  \brief funtion is implemented only for 5d meshes
   *  \brief compute the bilinear interpolation for a given point with designed results table (value function or minimal time)
   *  \param PP  : current mesh point coordinates
   *  \param tab : table data for computation
   */
  double    bilinearInterpolation5d   (const double* PP,const double* tab);

  /*! TODO
  double    interpolation3              (const double* PP,const double* tab);
  double    interpolationENO            (const double* PP,const double* tab);
  double    interpolationWENO           (const double* PP,const double* tab);
  */

  /*!
   *  \brief funtion is independant of the mesh dimension
   *  \brief prepare the boundary value in a periodic mesh
   *  \param tab  : table data to be made periodic copies
   */
  void      setPeriodic                 (double* tab);



  /*! #VERSION2
   *  \brief compute the optimal control for a given point in optimal trajectory computing
   *  \param y   : current mesh point coordinates
   *  \param h   : discretized space step
   *  \param t   : current time in the mainloop algorithm
   *  \param val : computed value with the optimal control
   *  \return optimal control's id
   */
  //int       find_optimal_control        (const double* y, double h, double t, double& val);
  int       find_optimal_control_val    (const double* y, double h, double t, double& val, double* vtab);
  int       find_optimal_control_val_mpi(const double* y, double h, double t, double& val, double* vtab, int handler);


  // AD: jan 2017: adaptive reconstruction method
  /*
   * optimal trajectory reconstruction by value function with gradient computation
   */



  /*!
   *
   * \brief Numerical integration of dynamics equations with Heun method
   *
   * @param[in]  yInit   : initial state
   * @param[in]  h       : time step
   * @param[in]  tau     : initial time
   * @param[out] inTheBounds : logic indicator of exit / not exit from the computation domain
   * @param[in]  nbSteps : number of time steps for the intgration
   * @param[in]  c       : number of the control value to apply for all steps
   * @param[out] res     : final state
   *
   */
  void      heunIntegrator               (const double* yInit,int nbSteps, int c, double tau, double h, double *res, bool *inTheBounds);


  // ENCOURS [NOV 2016] : case 2 controls
  int       find_optimal_control_val2    (const double* y, double h, double t, double& val, double* vtab, int & c2opt);
  // -> [TO BECOME  find_optimal_control2]
  int       compute_optimal_trajectory2  (const double* initialpoint, const char* filename, int number, double* terminalpoint, double &terminaltime);

  /*! #VERSION2
   *  \brief compute the optimal trajectory for a given initial point and create the associated VTK files (cf saveVTK6 and saveTrajCAR)
   *  \param initialpoint : start point to compute the trajectory
   *  \param filename     : root name of files to be written
   *  \param start        : VTK rank to start the trajectory animation
   */
  int       compute_optimal_trajectory    (const double* initialpoint, const char* filename, int number, double* terminalpoint, double &terminaltime);
  int       compute_optimal_trajectory_val(const double* initialpoint, const char* filename, int number, double time_start, double* terminalpoint, double &terminaltime);




  /*! #VERSION2
   *  \brief Applies periodic conditions  to a given vector
   *  \param vect : given vector to correct
   */
  void      periodizePoint              (double * x);


  /*!
   *  \brief computes euclidean norm 
   *  \param point coordinates
   *  \return computed norm
   */
  double   norm                         (const double*)  const;


  /*!
   *  \brief test if the two double argument are equalled
   *  \param x   : point coordinates
   *  \param y   : point coordinates
   *  \return equal test between two doubles
   */
  bool    doublecomp                    (double x, double y)  const;




};


inline double HJB::norm(const double* point) const {
  double res=0;
  for(int d=0;d<dim;d++)
    res+=point[d]*point[d];
  return sqrt(res);
};

inline bool HJB::doublecomp(double x, double y)  const {
  return (abs(x-y)<1e-15) ? true : false;
}


inline double tk(double t)   {   return(1-exp(-t));  }
inline double atk(double v)  {   if (v>=1.) return(INF); else return(-log(1.-v));    }


#endif /* HJB_H */

