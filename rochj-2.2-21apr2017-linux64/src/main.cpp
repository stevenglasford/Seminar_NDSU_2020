#include "stdafx.h"
#include "HJB_FD.h"
#include "HJB_SL.h"

//- HERE PUT THE DESIRED "data_xx.h" FILE
#include "../data/data_simulation.h"


///////////////////////////////////////////////
/////////////////////  MAIN ///////////////////
///////////////////////////////////////////////

int main(int argc, char *argv[], char *env[]){


    //- SEQVERSION
    int num_procs=0;
    int my_id=0;

    //- xxFiles for main.cpp (xxFiles for HJB.cpp are defined separately)
    char      VFile    [100];
    char      tminFile [100];
    char*     filePrefix=FILE_PREFIX;
    sprintf(VFile,    "%s%s%s.dat",OUTPUTdir,filePrefix,VFile0);    // <== pour main.cpp
    sprintf(tminFile, "%s%s%s.dat",OUTPUTdir,filePrefix,tminFile0); // <== pour main.cpp


    //MPI_NUM_PROCS       = parallel_params.PROC_NUMS;
    if (my_id==0){
      printf("//---------------------------------------------------------------\n");
      printf("//- HJB-Solver V2.2  (March, 2017)                               \n");
      printf("//- Olivier Bokanowski, Hasnaa Zidani, Junyi Zhao, Anya Desilles \n");
      printf("//---------------------------------------------------------------\n");
    }

    if (my_id==0)
      printf("HJB SOLVER : STABLE VERSION\n");

    if (my_id==0){
      if (METHOD==MFD) printf("METHOD : FD (Finite Differences)\n");
      if (METHOD==MSL) printf("METHOD : SL (Semi Lagrangian)\n");
    }

    //for some user data definitions
    init_data(); 

    int i, numthreads=1, griddim=2;
    int taille_maillage[DIM];
    for(i=0;i<DIM;i++){
      taille_maillage[i]=ND[i];
      //printf("taille_maillage[%i]=%i\n",i,taille_maillage[i]);
    }
    int NCDuser[cDIM], NCuser=1;
    for(i=0;i<cDIM;i++){
      NCDuser[i]=NCD[i];
      NCuser*=NCDuser[i];
    }

    if(argc>=3){
      for(i=1 ; i<argc ; i+=2){
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

    //- set the size of the mesh border
    //- 2016: BORDERSIZE[dim] initialized in data files 

    //- SEQVERSION
    unsigned mpidatas[griddim];

    static struct mrp  maillage_regulier_params =
    {  DIM, PRECOMPUTE_COORDS, MESH,
       (unsigned*) taille_maillage, (double*) XMIN, (double*) XMAX, (int*) PERIODIC,
       BORDERSIZE,
       OBSTACLE, g_obstacle, OBSTACLE_TILDE, g_obstacle_tilde,
       PRECOMPUTE_OBSTACLE,
       COMPUTE_IN_SUBDOMAIN, g_domain
    };

    static struct cp   commande_params =
    {  cDIM, (int*) NCDuser, NCuser, (double*)  UMIN, (double*)  UMAX  };

    int NC2=1;
    for(i=0;i<cDIM2;i++)
      NC2*=NCD2[i];
    static struct cp   commande_params2 =
    {  cDIM2, (int*) NCD2, NC2, (double*)  UMIN2, (double*)  UMAX2  };

    static struct pp   parallel_params =
    {  (unsigned) griddim, mpidatas, (unsigned) my_id, (unsigned) num_procs, (unsigned) numthreads  };

    static struct gp   general_params =
    {  METHOD, NAME, EXTERNALV0,
       COMPUTE_MAIN_LOOP, OPTIM, COMPUTE_TOPT, TOPT_TYPE, COMPUTE_VEX,
       SAVE_VF_ALL, SAVE_VF_ALL_STEP, SAVE_VF_FINAL, SAVE_VF_FINAL_ONSET,
       VALUE_PB, SAVE_VALUE_ALL,  SAVE_VALUE_ALL_STEP, SAVE_VALUE_FINAL,
       FORMAT_FULLDATA,
       CHECK_ERROR, CHECK_ERROR_STEP, C_THRESHOLD,
       T, DT, BOUNDARY, MAX_ITERATION, EPSILON,
       *dynamics, *dynamics2, *distributed_cost, *distributed_cost2, *v0, *Vex,
       *g_border, *g_bordermix,
       FILE_PREFIX
    };

    static struct ppp   postprocess_params =
    {  TRAJPT, (double*) initialpoint, TRAJ_METHOD, time_TRAJ_START,
       min_TRAJ_STOP, max_TRAJ_STOP, TARGET_STOP, *g_target, ADVERSE_METHOD, *u2_adverse,
       PRINTTRAJ,
       (int*) COUPE_DIMS, (double*) COUPE_VALS };

    /*
    static struct fdp  FD_params =
    {  TYPE_SCHEME, TYPE_RK, COMMANDS, CFL, *Hnum, *compute_Hconst };

    HJB_FD * pb = new HJB_FD(maillage_regulier_params, commande_params, commande_params2, parallel_params, general_params, postprocess_params, FD_params);

    if(COMPUTE_MAIN_LOOP){
      //pb->initFiles();
      pb->mainloop();
    }
    else {
      pb->loadTmin(tminFile);
      pb->loadV(VFile);
    }

    pb->postprocess();
    delete pb;
    */

    if(general_params.METHOD == MFD){
      static struct fdp  FD_params =
      {  TYPE_SCHEME, TYPE_RK, COMMANDS, CFL, *Hnum, *discount_factor,
         *compute_Hconst };
      HJB_FD * pb = new HJB_FD(maillage_regulier_params, commande_params, commande_params2, parallel_params, general_params, postprocess_params, FD_params);
      if(COMPUTE_MAIN_LOOP){
        //pb->initFiles();
        pb->mainloop();
      }
      pb->postprocess();
      delete pb;
    }
    else if(general_params.METHOD == MSL){

      if (num_procs>=2){ printf("MPI not programmed for SL method. Abort.\n"); exit(1);}
      
      static struct slp  SL_params =
      {  TYPE_SCHEME, TYPE_RK, ORDER, TYPE_STA_LOOP, INTERPOLATION, P_INTERMEDIATE,
         PARAMP, *funcR, *funcY, *discount_factor };
      HJB_SL * pb = new HJB_SL(maillage_regulier_params, commande_params, commande_params2, parallel_params, general_params, postprocess_params, SL_params);
      if(COMPUTE_MAIN_LOOP){
        //pb->initFiles();
        pb->mainloop();
      }
      pb->postprocess();
      delete pb;
    }
    else{ 
      printf("ERROR: this METHOD (=%i) not programmed. Abort.\n",general_params.METHOD);
      exit(1);
    }


    //for some user data post-treatment
    post_data();

    return(0);
}

