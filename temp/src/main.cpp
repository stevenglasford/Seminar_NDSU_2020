#include "stdafx.h"
#include "HJB_FD.h"
#include "HJB_SL.h"

//-----------------------------
//- INCLUDE of default file
//-----------------------------
//- include here the "data_xx.h" file
#include "../data/data_simulation.h"




///////////////////////////////////////////////
/////////////////////  MAIN ///////////////////
///////////////////////////////////////////////

int main(int argc, char *argv[], char *env[]){


  //- SEQ VERSION
  int num_procs=0;
  int my_id=0;

  //- xxFiles for main.cpp (xxFiles for HJB.cpp are defined separately)
  char      VFile    [100];
  char      tminFile [100];
  char*     filePrefix=FILE_PREFIX;
  sprintf(VFile,    "%s%s%s.dat",OUTPUTdir,filePrefix,VFile0);    // <== pour main.cpp
  sprintf(tminFile, "%s%s%s.dat",OUTPUTdir,filePrefix,tminFile0); // <== pour main.cpp

#ifdef _WIN32
  OUTPUTdir=OUTPUTdirWin; 
#endif


  if(my_id==0){

    printf("//---------------------------------------------------------------\n");
    printf("//- ROC-HJ Solver V2.5.4 (2019) beta                             \n");
    printf("//- Olivier Bokanowski, Hasnaa Zidani, Anya Desilles             \n");
    printf("//---------------------------------------------------------------\n");

    printf("HJB SOLVER : STABLE VERSION\n");

    if (METHOD==MFD) 
      printf("METHOD     : FD (Finite Differences)\n");
    else if (METHOD==MSL) 
      printf("METHOD     : SL (Semi Lagrangian)\n");
    else{
      printf("METHOD     : %i NOT PROGRAMMED. Abort.\n",METHOD);
      exit(1);
    }

  }


  //- general data init-treatment
  init_data_general();

  //----------------------------------------------------
  //- step_HJB: general index that can be used in the data to make several HJB computations
  //- beginStep_HJB: in data (default is 0)
  //- endStep_HJB:   in data (default is 0)
  //----------------------------------------------------

  //int beginStep_HJB=0; //- in case not defined
  for(step_HJB=beginStep_HJB; step_HJB<=endStep_HJB; step_HJB++)
  {

    //---------------------------------------------------
    //- declarations specifiques - dans le init_data() ou ici:
    //---------------------------------------------------

    //- for some general user data definitions
    init_data();

    int numthreads=1, griddim=2;
    unsigned taille_maillage[DIM];
    for(int i=0;i<DIM;i++){
      taille_maillage[i]=ND[i];
      //printf("taille_maillage[%i]=%i\n",i,taille_maillage[i]);
    }
    int NCDuser[cDIM], NCuser=1;
    for(int i=0;i<cDIM;i++){
      NCDuser[i]=NCD[i];
      NCuser*=NCDuser[i];
    }

    int NC2=1;
    for(int i=0;i<cDIM2;i++)
      NC2*=NCD2[i];

    if(argc>=3){
      for(int i=1 ; i<argc ; i+=2){
        if(strcmp(argv[i], "-nn") == 0){
          int nn = atoi(argv[i+1]);
          for(int i=0;i<DIM;i++)
            taille_maillage[i]=nn;
        }
        else if(strcmp(argv[i], "-nc") == 0){
          int nc = atoi(argv[i+1]);
          NCuser = 1;
          for(int i=0;i<cDIM;i++){
            NCDuser[i]=nc;
            NCuser*=nc;
          }
        }
        else if(strcmp(argv[i], "-nt") == 0){
          numthreads = atoi(argv[i+1]);
        }
        else if(strcmp(argv[i], "-nd") == 0){
          griddim = atoi(argv[i+1]);
        }
      }
    }

    //- SEQVERSION
    unsigned mpidatas[griddim];

    struct mrp  maillage_regulier_params =
    {  DIM, PRECOMPUTE_COORDS, MESH,
       (unsigned*) taille_maillage, (double*) XMIN, (double*) XMAX, (int*) PERIODIC,
       BORDERSIZE,
       OBSTACLE, g_obstacle, OBSTACLE_TILDE, g_obstacle_tilde,
       PRECOMPUTE_OBSTACLE,
       COMPUTE_IN_SUBDOMAIN, g_domain
    };

    struct cp   commande_params =
    {  cDIM, (int*) NCDuser, NCuser, (double*)  UMIN, (double*)  UMAX  };

    struct cp   commande_params2 =
    {  cDIM2, (int*) NCD2, NC2, (double*)  UMIN2, (double*)  UMAX2  };

    struct pp   parallel_params =
    {  (unsigned) griddim, mpidatas, (unsigned) my_id, (unsigned) num_procs, (unsigned) numthreads  };

    struct gp   general_params =
    {  METHOD, NAME, EXTERNALV0, EXTERNAL_TOPT,
       COMPUTE_MAIN_LOOP, 
       OPTIM, COMPUTE_TOPT, TOPT_TYPE, COMPUTE_VEX,
       SAVE_VF_ALL,    SAVE_VF_ALL_STEP,    SAVE_VF_FINAL, 
       SAVE_VEX_ALL,   SAVE_VEX_FINAL, 
       SAVE_TOPT_ALL,  SAVE_TOPT_FINAL, 
       SAVE_COUPE_ALL,     SAVE_COUPE_ALL_STEP,  SAVE_COUPE_FINAL,
       SAVE_COUPEEX_ALL,   SAVE_COUPEEX_FINAL, 
       SAVE_COUPETOPT_ALL, SAVE_COUPETOPT_FINAL, 
       SAVE_VF_ONSET_ALL,    SAVE_VF_ONSET_FINAL,
       SAVE_VEX_ONSET_ALL,   SAVE_VEX_ONSET_FINAL,
       SAVE_TOPT_ONSET_ALL,  SAVE_TOPT_ONSET_FINAL,
       NBPOINTS_Xuser, 
       XuserFile0, XuserVFile0, XuserVexFile0, XuserToptFile0,
       SAVE_COUPE_INTERPOLATION,
       SHOW_COORDINATES,
       BINARY,
       //SAVE_VF_FINAL_ONSET,
       //VALUE_PB, SAVE_VALUE_ALL, SAVE_VALUE_ALL_STEP, SAVE_VALUE_FINAL,
       VALUE_PB, SAVE_VALUE_ALL, SAVE_VALUE_FINAL,
       FORMAT_FULLDATA,
       CHECK_ERROR, CHECK_ERROR_STEP, CHECK_NEG_POINTS, C_THRESHOLD,
       T, DT, BOUNDARY, MAX_ITERATION, EPSILON,
       *dynamics, *dynamics2, 
       *feedback,
       *dynamicsGrad,
       *distributed_cost, *distributed_cost2, 
       *v0, *Vex,
       *g_border, *g_bordermix,
       FILE_PREFIX,EXTERNAL_FILE_PREFIX
    };

    struct ppp   postprocess_params =
    {  TRAJPT, (double*) initialpoint, TRAJ_METHOD, time_TRAJ_START,
       min_TRAJ_STOP, max_TRAJ_STOP, TARGET_STOP, *g_target, ADVERSE_METHOD, *u2_adverse,
       PRINTTRAJ,
       (int*) COUPE_DIMS, (double*) COUPE_VALS 
    };


    if(general_params.METHOD == MFD){
      struct fdp  FD_params =
      { TYPE_SCHEME, TYPE_RK, COMMANDS, CFL, *Hnum, *discount_factor,
        *compute_Hconst 
      };

      HJB_FD * pb = new HJB_FD(
         maillage_regulier_params, commande_params, commande_params2, parallel_params, general_params, postprocess_params, FD_params);

      {
        pb->initprocess();

        if(COMPUTE_MAIN_LOOP)
          pb->mainloop();

        pb->postprocess();
      }

      delete pb;

    }

    if(general_params.METHOD == MSL){

      if (num_procs>=2){ printf("MPI not programmed for SL method. Abort.\n"); exit(1);}

      struct slp  SL_params =
      {  TYPE_SCHEME, TYPE_RK, ORDER, TYPE_STA_LOOP, INTERPOLATION, P_INTERMEDIATE,
         PARAMP, *funcR, *funcY, *discount_factor };

      HJB_SL * pb = new HJB_SL(maillage_regulier_params, commande_params, commande_params2, parallel_params, general_params, postprocess_params, SL_params);

      pb->initprocess();

      if(COMPUTE_MAIN_LOOP)
        pb->mainloop();
      
      pb->postprocess();

      delete pb;

    }


    post_data();
    

    if( (general_params.METHOD != MFD) &&  (general_params.METHOD != MSL) )
    { 
      printf("ERROR: this METHOD (=%i) not programmed. Abort.\n",general_params.METHOD);
      exit(1);
    }

  }//- end of step_HJB loop

  post_data_general(); //- general data post-treatment

  return(0);
}

