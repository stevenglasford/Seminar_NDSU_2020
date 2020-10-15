//---------------------------------
//- optimal trajectory for a given point and output in a  (.dat or VTK) format 
//- BASED ON VALUE
//- [MAR 2017] 1 control case
//---------------------------------


int HJB::compute_optimal_trajectory_val
  (const double* initialpoint, const char* filename, int number, double time_start, double* terminalpoint, double &terminaltime)
{

  //- NON ADAPTATIVE VERSION 2014; 2017;
  //- outputs: 
  //-   terminalpoint[dim]
  //- 	terminaltime of the trajectory (if succcessful)
  //-   success=0/1 
  //-   "number" is not used.
  //- starting time "time_start" included : 
  //- if time_start<0        then starts traj. reconstruction at first time t>=0 such that Value(t,initialpoint)>0.
  //- if time_start in [0,T] then starts traj. reconstruction at time t=time_start.

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
  //double step;

  FILE * pFile;

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
  if (time_start<0)
    time[0]=0.0;
  else
    time[0]=time_start;

  //- time t (traj)
  t=T-time[0];

  if (PRINTTRAJ && my_rank==handler){
    printf("..Initialisation (T=%8.5f): time[0]=%8.5f, t=T-time[0]=%8.5f, time_start=%8.5f\n",T,time[0],t,time_start);
  }


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

  //-----------
  //- load Dt: 
  //-----------
  char  fileDt[99];
  sprintf(fileDt,"%s%s%s.dat",OUTPUTdir,filePrefix,DtFile0);
  printf("..Loading Dt (%s) : ",fileDt);
  pFile = fopen (fileDt,"r");
  if (pFile==NULL){ printf("19xx.Impossible to open the file %s for reading !\n",fileDt); exit(1);}
  int zz=fscanf(pFile,"%lf",&Dt); zz++;
  fclose(pFile);

  if (PRINTTRAJ && my_rank==handler)
    printf("Dt(HJB):=%8.5f; ",Dt);

  Dt=Dt*savevfallstep;

  if (PRINTTRAJ && my_rank==handler)
    printf("Dt<-Dt*SAVE_VFALL_STEP:=%8.5f;\n",Dt);

  //-----------
  //- using hx=Dt 
  //- (where hx=time step of trajectory construction
  //-  and Dt=time step used for value (HJB) computations)
  //-----------
  hx=Dt;

  int N=int(ceil(T/Dt));
  int ii=N; //- counter from N to 0

  if (PRINTTRAJ && my_rank==handler)
    printf("..Starting reading  VFN.dat file with N=%3i;\n",N);

  int start=0;

  while(!end && ii>=1 && time[n]<=T-Dt){

    int handler=0; //- sequential code: handler=0, my_rank=0 always

    ii=ii-1;

    //if (1){
    //  //- load v 
    //  char buffer [99];
    //  sprintf(buffer,"OUTPUT/VF%d.dat",ii);
    //  loadV(buffer,0); //- load in "v" the VFii.dat file ; 0=no printings
    //}
    
    //- Note:
    //- time_start<0:  will start only when v(.,x) becomes > 0
    //- time_start>=0: will start at approx at time t_n=time_start

    if (time_start<0.0){
      //- load v 
      char buffer [99];
  	sprintf(buffer,"%s%s%s%d.dat",OUTPUTdir,filePrefix,VFile0,ii);
      loadV(buffer,0); //- load in "v" the VFii.dat file ; 0=no printings
    }

    if (time_start>=0.0 && time[n]>=time_start){
      //- load v only if time[n]>time_start
      char buffer [99];
  	sprintf(buffer,"%s%s%s%d.dat",OUTPUTdir,filePrefix,VFile0,ii);
      loadV(buffer,0); //- load in "v" the VFii.dat file ; 0=no printings
    }


    //---------------------------------
    //- HERE FOR STARTING TRAJECTORY
    //---------------------------------
    //- will start trajectory reconstruction when v(ot)>0

    if (time_start<0.0){
    if (!start){
      if ((*this.*interpolation)(ot[n],v)>0.0){
	//printf("passage start=1\n");
        start=1;
	success=1;
	//time[n]=(N+1-ii)*Dt;
      }
      else{
        t=t-hx; //- go to next step (start=0 unchanged)
	time[n]=time[n]+Dt;
      }
    } 
    }


    // this is almost same as previous but to start after time_start
    if (time_start>=0.0){
    if (!start){
      if (time[n]>=time_start){
        start=1;
	success=1;
	printf("(start:=1 at time[n]=%f);\n",time[n]); 
      }
      else{
        t=t-hx; //- go to next step (start=0 unchanged)
        time[n]=time[n]+Dt;
      }
    } 
    }

    //----------------------------------
    //- WHEN TRAJECTORY HAS BEEN STARTED
    //----------------------------------
    if (start==1){

      u_star[n]=find_optimal_control_val(ot[n],hx,t,val,v);

      if (u_star[n]==-1){ // but this should never happend
        printf("no admissible controls. Abort traj. reconstruction\n");
	success=0;
	//exit(1);
      }

      for(d=0; d<dim; d++)
        position[d]=ot[n][d];
    
      //- computing local speed to determine a local adaptive time step (not used)
      // maximum=0.;
      // for(c=0;c<ncall;c++){
      //  (*dynamics)(ot[n],(*u)[c],t,ff);
      //   for(d=0;d<dim;d++){
      //     maximum=max(maximum,abs(ff[d]);
      //   }
      // }
      // h=mesh->Dx_min/maximum;
      // hx=0.5 * h;

      //- (using hx=Dt here)
      t=t-hx;
      time[n+1]=time[n]+Dt;

    
      //- TRAJECTORY PRINTINGS
      if(PRINTTRAJ && my_rank==handler){
        for(d=0;d<dim;d++)
          printf("x%i=%6.3f, ",d+1,ot[n][d]);
        printf("ii=%3i, ",ii);
        printf("iter n=%3i, time[n]=%6.3f, hx=%6.3f ",n,time[n],hx);
        for(int i=0;i<(*u).dimC;i++)
          printf("u%i=%5.3f, ",i+1,(*u)[u_star[n]][i]);
        printf("v=%6.3f, ", (*this.*interpolation)(ot[n],v));
        printf("t=%5.3f (traj)", t);
        printf("\n");
      }

  
      // compute next step by Heun scheme
      // using optimal control found and RK2 scheme
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


      //--------------------------------
      //- stoppting test: out of domain
      //--------------------------------
      test=0;
      for(d=0;d<dim;d++){
        if ((ot[n+1][d] < mesh->GLOBAL_lowBounds[d] || ot[n+1][d] > mesh->GLOBAL_highBounds[d]) && !mesh->periodic[d])
          test=1;
      }

      if(test){
        printf("Getting out of the domain: n=%2i ",n); 
        for (d=0;d<dim;d++) 
          printf("x%i=%5.2f; ",d+1,ot[n][d]); 
        printf("\n");
	success=0;
        end=1;
      }

      n++;
    }

    //------------------------------------------------------
    //- end of "if((*this.*interpolation)(ot[n],v)<0.0) ..."
    //------------------------------------------------------

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

    if (PRINTBUG && my_rank==handler)
      printf("end of step n=%4i ...\n",n);

  }//- end n loop


  //------------------BEGIN
  if (!special_value){

    if (u_star[n]==-1){
      success=0;
      // no control found - tmin = +INF at terminal point.
      n--;
      // terminalpoint, terminal time:
      for(d=0; d<dim; d++)
        terminalpoint[d]=ot[n][d];
      terminaltime=time[n];
    }
    else{
      for(d=0; d<dim; d++)
        position[d]=ot[n][d];
  
      u_star[n]=u_star[n-1]; // this just to have final u_star[n] value to fprint
      if(PRINTTRAJ && my_rank==0)
        printf("..time[n]= %6.3f, total steps n=%6i;\n",time[n],n+1);
  
      // terminalpoint, terminal time:
      for(d=0; d<dim; d++)
        terminalpoint[d]=ot[n][d];
      terminaltime=time[n];
    } 

  }
  

  if (my_rank==0){
    printf("(terminaltime=%8.5f, success=%1i)\n",terminaltime,success);
  }

  //- save trajectory
  if (PRINTTRAJ && my_rank==0)
    printf("Saving optimal trajectory...");

  if (my_rank==0){
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


  //- Load again VF.dat into "v"
  char buffer [99];
  sprintf(buffer,"%s%s%s.dat",OUTPUTdir,filePrefix,VFile0);
  loadV(buffer,0); //- load in "v" the VFii.dat file ; 0=no printings

  for(n=0;n<MAX_ITERATION;n++)
    delete[] ot[n];
  delete[] ot;
  delete[] u_star;
  delete[] time;

  return success;
}

