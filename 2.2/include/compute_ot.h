//---------------------------------
//- optimal trajectory for a given point and output in a  (.dat or VTK) format 
//- BASED ON topt
//- [MAR 2017] 1 control case
//---------------------------------


int HJB::compute_optimal_trajectory
  (const double* initialpoint, const char* filename, int number, double* terminalpoint, double &terminaltime)
{

  //- NON ADAPTATIVE VERSION 2014; 2017;
  //- outputs: 
  //-   terminalpoint[dim]
  //- 	terminaltime of the trajectory (if succcessful)
  //-   success=0/1 
  //-   "number" is not used.

  //PRINTTRAJ=1;       //- [now in data] use PRINTTRAJ=1 to more verbose
  int success=0;       //- success=0/1 for output
  double *time, **ot;  //- optimal trajectory savings;
  int n, *u_star;
  int handler=0;       //- sequential code: handler=0, my_rank=0 always

  time   = new double [MAX_ITERATION];
  ot     = new double*[MAX_ITERATION];
  u_star = new int    [MAX_ITERATION];


  for(n=0;n<MAX_ITERATION;n++)
    ot[n] = new double[dim];

  double position[dim],ff[dim],ff2[dim];
  int test, d, k, end=0;
  double hx=0., val;
  double step;

  if(PRINTTRAJ && my_rank==handler)
    printf("\n");

  if (PRINTBUG)
    printf("Compute optimal trajectory:\n");

  n=0;
  for(d=0; d<dim; d++){
    position[d]=initialpoint[d];
    ot[n][d]=position[d];
  }

  //- time[0]
  time[0]=0.0;

  //- time t (traj)
  t=(*this.*interpolation)(ot[n],topt); //t=T;

  if(PRINTTRAJ && my_rank==handler)
    printf("..Initialisation t=tmin=%8.3f (step n=0)\n",t);


  //-----------
  //- special treatment if VALUE_PB and z=INF
  //-----------
  int special_value=0;
  if (VALUE_PB && position[dim-1]>=INF){
    end=1;
    success=0;
    u_star[0]=0;
    special_value=1; // special case when last component is inf. 
  }


  n=0;
  while(!end && n<MAX_ITERATION-1){

    //int handler=0; //- sequential code: handler=0, my_rank=0 always

    for(d=0; d<dim; d++)
      position[d]=ot[n][d];

    //- computing max local speed in order to determine time step
    if(my_rank==handler) //- this test is always true in this (sequential) code.
    {
      //- 1 control
      double maximum=0.0;
      for(int c=0;c<ncall;c++){
        (*dynamics)(ot[n],(*u)[c],t,ff);
        for(d=0;d<dim;d++){
          maximum=max(maximum,abs(ff[d]));
        }
      }
      //- time adaptation 
      if (maximum<=1e-14) break; // stopping in case of small dynamics
      double h=mesh->Dx_min/(maximum);
      //double h=mesh->Dx_min/(1e-10+maximum);
      step=0.5;
      hx=h*step;
    }

    //- time[n+1]=time[n]+hx;
    //- t=t-hx;

    // Mesh adapt can be introduced here
   
    u_star[n]=find_optimal_control_val(ot[n],hx,t,val,topt);

    if (u_star[n]==-1) 
      break;

    //- TRAJECTORY PRINTINGS
    if(PRINTTRAJ && my_rank==handler){
      for(d=0;d<dim;d++)
        printf("x%i=%6.3f, ",d+1,ot[n][d]);
      printf("iter n=%3i, time[n]=%6.3f, hx=%6.3f; ",n,time[n],hx);
      for(int i=0;i<(*u).dimC;i++)                       // - 1st control
        printf("u%i=%6.3f, ",i+1,(*u)[u_star[n]][i]);
      printf("tmin=%5.3f, ", (*this.*interpolation)(ot[n],topt));
      printf("t=%6.3f (traj)", t);
      printf("\n");
    }

    if(my_rank==handler){
      //- 1 control:
      //-  compute next step
      //-  using optimal control found and RK2 scheme
      (*dynamics)(ot[n],(*u)[u_star[n]],t,ff);
      for(d=0;d<dim;d++)
        ot[n+1][d]=ot[n][d]+hx*ff[d];

      if(periodic_mesh)
        periodizePoint(ot[n+1]);

      // - Careful : remaining time so t+hx becomes t-hx
      (*dynamics)(ot[n+1],(*u)[u_star[n]],t-hx,ff2);

      for(d=0;d<dim;d++)
        ot[n+1][d]=ot[n][d]+0.5*hx*(ff[d]+ff2[d]);
      if(periodic_mesh)
        periodizePoint(ot[n+1]);

    }

    //--------------------------------
    //- stoppting iest: out of domain
    //--------------------------------
    //test = !mesh->isInDomain(ot[n+1]);
    test=0;
    for(d=0;d<dim;d++){
      if((ot[n+1][d] < mesh->GLOBAL_lowBounds[d] || ot[n+1][d] > mesh->GLOBAL_highBounds[d]) && !mesh->periodic[d])
        test=1;
    }

    if(test){
      if(my_rank==handler){
        printf("Getting out of the domain: n=%2i ",n); 
        for (d=0;d<dim;d++) 
          printf("x%i=%5.2f; ",d+1,ot[n][d]); 
        printf("\n");
      }
      success=0;
      end=1;
    }

    // 2017 
    //----------------------------------------------------------------------------
    //- recall that epsOPTIM=1 if OPTIM=MAXIMUM, and epsOPTIM=-1 if OPTIM=MINIMUM
    //----------------------------------------------------------------------------


    // stopping test : if val<=min_TRAJ_STOP (for instance for tmin pb, epsOPTIM= 1)
    //              or if val>=max_TRAJ_STOP (for instance for tmax pb, epsOPTIM=-1)
    if((val-min_TRAJ_STOP)*epsOPTIM<=0.0){
      if(my_rank==handler)
        printf("(val-min_TRAJ_STOP)*epsOPTIM<=0 ==> success! ");
      end=1;
      success=1;
    }
   

    if((val-max_TRAJ_STOP)*epsOPTIM>=0.0){
      if(my_rank==handler){
        printf("[(val-max_TRAJ_STOP)*epsOPTIM>=0 ==> stop; ");
        printf(" val=%5.1f, max_TRAJ_STOP=%5.1f, epsOPTIM=%i]; ",val,max_TRAJ_STOP,epsOPTIM);
      }
      end=1;
      success=0;
    }

    //- stopping test : target test (if TARGET_STOP=1)
    if (TARGET_STOP){
      double vv=(*g_target)(time[n],ot[n]);
      if(vv<=0.0){
        if(my_rank==handler)
          printf("[Target reached: success; Stopped with g_target(x)=%6.3f :<=0]; ",vv);
        end=1;
        success=1;
      }
    }

    //- stopping test (only for tmax problems) : if time[n]>=T
    //- stopping test (only for tmax problems) : or if target (v0) is used
    if(TOPT_TYPE==1){
      if(time[n]>=T){
        if(my_rank==handler)
          printf("[Maximal time reached: success.]; ");
        success=1;
        end=1;
      }
      else{
        double vv=(*v0)(ot[n]);
        if(vv>0.0){
          if(my_rank==handler)
            printf("[Target left: failure; Stopped with v0(ot[n])=%6.3f :<=0]; ",vv);
          success=0;
          end=1;
        }
      }
    }

    if (PRINTBUG && my_rank==handler)
      printf("end of step n=%4i ...\n",n);

    time[n+1]=time[n]+hx;
    t=t-hx;

    n++;

  }//- end n loop


  //----------------------------------------
  //- TRAJECTORY PRINTINGS - LAST POINT
  //----------------------------------------
  if(PRINTTRAJ && my_rank==handler){
      printf("\n");
      for(d=0;d<dim;d++)
        printf("x%i=%6.3f, ",d+1,ot[n][d]);
      printf("iter n=%3i, time[n]=%6.3f, hx=%6.3f; ",n,time[n],hx);
      for(int i=0;i<(*u).dimC;i++)                       // - 1st control
        printf("u%i=%6.3f, ",i+1,0.0);
      printf("tmin=%5.3f, ", (*this.*interpolation)(ot[n],topt));
      printf("t=%6.3f (traj)", t);
      printf(" [LAST POINT]\n");
  }

  
  //------------------BEGIN
  if (!special_value){

  if (u_star[n]==-1 && n==0){
    for(d=0; d<dim; d++)
      terminalpoint[d]=ot[n][d];
    terminaltime=time[n];
    u_star[n]=0;
  }
  else if (u_star[n]==-1 && n>=1){
    success=0;
    // no control found - tmin = +INF at terminal point.
    n--;
    // terminalpoint, terminal time:
    if (n==-1){
      for(d=0; d<dim; d++)
        terminalpoint[d]=ot[n][d];
    }
    else
      for(d=0; d<dim; d++)
        terminalpoint[d]=ot[n][d];
    terminaltime=time[n];
  }
  else{

    for(d=0; d<dim; d++)
      position[d]=ot[n][d];

    u_star[n]=u_star[n-1]; // this just to have final u_star[n] value to fprint
    if(PRINTTRAJ && my_rank==0)
      printf("..time[n]= %6.3f, n=%6i;\n",time[n],n);

    // terminalpoint, terminal time:
    for(d=0; d<dim; d++)
      terminalpoint[d]=ot[n][d];
    terminaltime=time[n];

  } 

  }
  else
  {
    for(d=0; d<dim; d++)
      position[d]=ot[n][d];

    u_star[n]=u_star[n-1]; // this just to have final u_star[n] value to fprint
    if(PRINTTRAJ && my_rank==0)
      printf("..time[n]= %6.3f, n=%6i;\n",time[n],n);

    // terminalpoint, terminal time:
    for(d=0; d<dim; d++)
      terminalpoint[d]=ot[n][d];
    terminaltime=time[n];
  }
  //------------------END
 


  if (my_rank==0){
    printf("(terminaltime=%8.5f, success=%1i)\n",terminaltime,success);
  }

  //- save trajectory
  if (PRINTTRAJ && my_rank==0)
    printf("Saving optimal trajectory...");

  if (my_rank==0){
    FILE * pFile;
    pFile = fopen (filename,"w");
    if(pFile==NULL){
      printf("19xx.Impossible to open the file %s for writting !\n",filename);
      exit(1);
    }
    else{
      for (k=0; k<=n; k++){
        for (d=0; d<dim; d++)
          fprintf(pFile,"%8.5f ", ot[k][d]);
        fprintf(pFile,"     ");
        for (int i=0; i<(*u).dimC; i++)
          fprintf(pFile,"%8.5f ", (*u)[u_star[k]][i]); //- 1st control value
        fprintf(pFile,"     ");
        fprintf(pFile,"%8.5f ", time[k]); //- time
        fprintf(pFile,"\n");
      }
    }
    fclose (pFile);

    if (PRINTTRAJ && my_rank==0)
      printf("done.\n");
  }


  for(n=0;n<MAX_ITERATION;n++)
    delete[] ot[n];
  delete[] ot;
  delete[] u_star;
  delete[] time;

  return success;
}

