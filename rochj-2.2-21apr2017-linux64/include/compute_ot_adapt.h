//---------------------------------
//- optimal trajectory for a given point and output in a  (.dat or VTK) format 
//- based on : topt, adaptive case
//---------------------------------

/*!
 * \brief Optimal traj reconstruction based on topt, with adaptative time step
 *
 * @param initialpoint[in]  : initial state 
 * @param filename[in]      : file name to save trajectory
 * @param start
 * @param terminalpoint[out]: final state 
 * @param terminaltime[out] : final time of trajectory
 * @return succes           : success indicator if target reached or not
 */

int HJB::compute_optimal_trajectory_adapt
  (const double* initialpoint, const char* filename, int number, double* terminalpoint, double &terminaltime)
{

  //int PRINTTRAJ=0;     //- [now in data] PRINTTRAJ=1 to more verbose
  int PRINTloc=0;
  //int PRINTtarget=1;
  double *time, **ot, *position;
  int n, *u_star;

  time   = new double [MAX_ITERATION];
  ot     = new double*[MAX_ITERATION];
  u_star = new int    [MAX_ITERATION];


  for(n=0;n<MAX_ITERATION;n++)
    ot[n] = new double[dim];

  position = new double[dim];

  double dt0, dt;
  int success;
  success=0;


  if (PRINTloc)
    cout<< " compute traj debut\n *********************************\n";

  FILE * pFile;
  char  fileDt[99];
  sprintf(fileDt,"%s%s%s.dat",OUTPUTdir,filePrefix,DtFile0);
  if (PRINTloc) printf("Loading Dt   (%15s) : ",fileDt);
  pFile = fopen (fileDt,"r");
  if (pFile==NULL)
  {
    cout << "Impossible to open the file for reading!" << endl;
    exit(1);
  }
  int zz=fscanf(pFile,"%lf",&dt0); zz++;
  if (PRINTloc) 
    cout<<dt0<<endl;
  fclose(pFile);


  double *ff = new double [dim];
  double treste;

  int   c, d, k, end=0;                             //time
  double h=0., hx,   support, maximum, val;
  //- stopping threshold
  //double eps=0.00001;

  //double seuilTmin=0.001;//- stopping threshold

  ofstream outfile;

  if (PRINTloc) cout << endl << "compute optimal trajectory initial position = ";

  for(d=0; d<dim; d++)
  {
    position[d]=initialpoint[d];
    if (PRINTloc) cout<<  " "<< position[d];
  }
  if (PRINTloc) cout<<endl;

  n=0;
  for(d=0; d<dim; d++){
    ot[n][d]=position[d];
  }


  time[n]=0.;

  treste=T;//

  double tInterpol;


  //int special_value=0;
  if (VALUE_PB && position[dim-1]>=INF){
    end=1;
    u_star[0]=0;
    //special_value=1;
  }
  else
  {
    tInterpol=(*this.*interpolation)(position,topt);
    if (PRINTloc) cout<< " temps interpole "<<tInterpol<<endl;

    if(tInterpol>2.0*T)
    {
      if (PRINTloc) cout<< " out of capture basin abort \n ";
      end=1;
      u_star[n]=0;
    }
  }

  int uOpt;
  double tOpt=INF,hOpt;
  bool inTheBounds;


  int nbSteps,nbOpt=1, nbS;
  bool testTime ;
  int iter=0;

  if (PRINTloc) cout<< " dt0 = "<<dt0<<endl;
  while (!end && (time[n]<=T)&& (n<MAX_ITERATION-1))
  {
    int handler=0; //- sequential code: handler=0, my_rank=0 always

    /*!
     * calcul de la vitesse locale pour d�terminer le pas  de temps
     */
    maximum=0.;

    for(c=0;c<ncall;c++)
    {
      (*dynamics)(ot[n],(*u)[c],time[n],ff);
      for(d=0;d<dim;d++)
      {
        support=ABS(ff[d]);
        if(support>maximum)  maximum=support;
      }
    }

    /*!
     * On calcule  le pas  de temps  de façon adaptative
     *  en prenant en compte  la norme de la dynamique  au point  courent de
     *  la trajectoire
     */

    if(maximum>0)
    {
      h=min(T-time[n], sqrt(dt0)/3);
      if (h<dt0)
        dt=0.999*h;
      else
        dt=dt0;
      nbSteps=floor(h/dt);
    }
    else
    {
      h= min(T-time[n],sqrt(dt0));
      if (h<dt0)
        dt=0.999*h;
      else
        dt=dt0;
      nbSteps=floor(h/dt);
      end=1;
      break;
    }

    hx=dt*((double) nbSteps);

    if (PRINTloc) cout<< "  hx= "<<hx<< " nbsteps = "<<nbSteps<< " time = "<<time[n]<<endl;
    /*!
     * Une adaptation  du pas peut être nécessaire
     * si la trajectoire optimale passe près de la frontière de l'ensemble
     * atteignable ; Si  pour toutes  les valeurs possibles du contrôle
     * l'image se trouve  dans une maille  de frontière à l'extérieur
     * on diminue le pas jusqu'à ce  qu'on  arrive à obtenir au moins une image à l'intérieur
     * , de façon à avoir  une interpollation  de topt <T
     */
    testTime=false;
    iter=0;

    uOpt=find_optimal_control_adapt(ot[n],dt,time[n],val,nbSteps,topt);

    hOpt=hx;
    nbOpt=nbSteps;
    nbS=nbSteps;
    if(uOpt>=0)
    {
      heunIntegrator(ot[n], nbSteps,  uOpt, time[n], dt, ot[n+1],  &inTheBounds);
      /*!
       * Dans le cas où certaines variables sont périodiques
       * les composantes correspondantes du vecteur image  sont
       * corrigées par périodicité de façon à se trouver dans
       * l'intervalle  défini  comme période
       */
      if(inTheBounds)
        treste=(*this.*interpolation)(ot[n+1],topt);
      else
        treste=INF;
    }
    else
      treste=INF;

    if (PRINTloc)
      cout<<" iter= "<<iter<< " hx= "<<hx<< " topt dans le nouveau point "<<treste<<endl;

    if(   ((treste >2.0*T)))
    {
      if (PRINTloc) 
	cout<<  "treste= "<<treste<<  " hx= "<<hx<< " adapt\n";
      testTime=true;
    }

    while(testTime & (nbS>1))
    {
      /*!
       * On itere tant qu'il est nécessaire de diminuer le pas
       */
      testTime=false;
      nbSteps=nbS;
      hx=dt*((double) nbSteps);
      u_star[n]=find_optimal_control_adapt(ot[n],dt,time[n],val,nbSteps,topt);

      if(u_star[n]>=0)
      {
        heunIntegrator(ot[n], nbSteps,  u_star[n], time[n], dt, ot[n+1],  &inTheBounds);
        /*!
         * Dans le cas où certaines variables sont périodiques
         * les composantes correspondantes du vecteur image  sont
         * corrigées par périodicité de façon à se trouver dans
         * l'intervalle  défini  comme période
         */

        if(inTheBounds)
          treste=(*this.*interpolation)(ot[n+1],topt);
        else
          treste=INF;
      }
      else
        treste=INF;

      if(treste>1.*T)
      {
        testTime=true;
        if(treste<tOpt)
        {
          uOpt=u_star[n];
          hOpt=hx;
          tOpt=treste;
          nbOpt=nbSteps;
        }
        nbS-=1;
        iter++;
      }
      else
      {
        if (PRINTloc)
          cout<< " \n$$$$$$$$$$\ncorrection de pas avec succes !  nbSteps = "<<nbSteps<< " treste = " <<treste<<"\n$$$$$$$$$$\n";
        uOpt=u_star[n];
        hOpt=hx;
        tOpt=treste;
        nbOpt=nbSteps;
      }
    }

    if(uOpt>=0)
      u_star[n]=uOpt;
    else
      u_star[n]=0;

    hx=hOpt;
    nbSteps=nbOpt;
    time[n+1]=time[n]+hx;
    if (PRINTloc) 
      cout<< " u_star="<<u_star[n]<<endl;
    if(u_star[n]>=0)
    {
      heunIntegrator(ot[n], nbSteps,  u_star[n], time[n], dt, ot[n+1],  &inTheBounds);

      if(inTheBounds)
        treste=(*this.*interpolation)(ot[n+1],topt);
      else
        treste=INF;
    }
    else
      treste=INF;
    if(PRINTTRAJ)
    {
      printf(" calcul  de u: n=%3i,    time[n]=%5.2f, h=%6.3f, iu=%3i,  topt=%5.4f\n",
          n,  time[n],hx,u_star[n], treste);
      cout<< " x=";
      for(int lm=0;lm<dim;lm++)
        cout<<" "<< ot[n+1][lm];
      cout<< " \n  controle u= ";
      for(int lm=0;lm<2;lm++)
        cout<<" "<< (*u)[u_star[n]][lm];
      cout<<endl;
    }

    n++;
    tOpt=INF;

    if(!inTheBounds)
      end=1;
    if(val>=INF)
      end=1;


    // 2017 
    //----------------------------------------------------------------------------
    //- recall that epsOPTIM=1 if OPTIM=MAXIMUM, and epsOPTIM=-1 if OPTIM=MINIMUM
    //----------------------------------------------------------------------------

    // stopping test :
    if((val-max_TRAJ_STOP)*epsOPTIM>=0.0){
      if(my_rank==handler){
        printf("(val-max_TRAJ_STOP)*epsOPTIM>=0 ==> stop - not able to reach target. ");
        printf("(val=%10.5f, max_TRAJ_STOP=%10.5f, epsOPTIM=%2i)",val,max_TRAJ_STOP,epsOPTIM);
      }
      end=1;
      success=0;
    }

    // stopping test : if val<=min_TRAJ_STOP (for instance for tmin pb, epsOPTIM= 1)
    //              or if val>=max_TRAJ_STOP (for instance for tmax pb, epsOPTIM=-1)
    if(TOPT_TYPE==0){
      if((val-min_TRAJ_STOP)*epsOPTIM<=0.0){
        if(my_rank==handler)
          printf("(val-min_TRAJ_STOP)*epsOPTIM<=0 ==> success! ");
        end=1;
        success=1;
      }
    }

    // stopping test : target test (if TARGET_STOP=1)
    if (TARGET_STOP && TOPT_TYPE==0){
      double vv=(*g_target)(time[n],ot[n]);
      if(vv<=0.0){
        if(PRINTTRAJ && my_rank==handler)
          printf("(Target reached: success; [Stopped with g_target(ot[n])=%6.3f :<=0]; ",vv);
        end=1;
        success=1;
      }
    }


    // stopping test (only for tmax problems) : if time[n]>=T
    if(TOPT_TYPE==1){
      if(time[n]>=T){
        if(PRINTTRAJ && my_rank==handler)
          printf("Maximal time reached: success.\n");
        success=1;
        end=1;
      }
      // stopping test : if target (v0) is left  ICI
      double vv=(*v0)(ot[n]);
      if(vv>0.0){
        if(PRINTTRAJ && my_rank==handler)
          printf("(Target left: failure; [Stopped with v0(ot[n])=%6.3f :<=0]; ",vv);
        success=0;
        end=1;
      }
    }

    if (PRINTBUG && my_rank==handler)
      printf("end of step n=%4i .. end= %i  success = %i.\n",n, end, success);

  }

  /*
   *  on enregistre dans le tableau res  le dernier état  de la trajectoire
   *   et la durée pour un éventuel postraitement
   *   ces données  seront enregistrées dans un fichier par la fonction postprocess();
   */

  if (PRINTTRAJ)
    cout<< " final position  n= "<<n<< " = ";

  for(d=0; d<dim; d++){
    position[d]=ot[n][d];
    terminalpoint[d]=ot[n][d];
    if (PRINTTRAJ)
      cout<< "  "<<position[d];
  }
  terminaltime=time[n];
  if(n>0)
    u_star[n]=u_star[n-1];

  if (PRINTTRAJ)
    cout << "\n Terminal  time= " << time[n]<< "\n test de res  time= "  << n+1 << endl;


  //------------------END

  if (my_rank==0){
    printf("(terminaltime=%8.5f, success=%1i)\n",terminaltime,success);
  }

  //- save trajectory
  if (PRINTTRAJ && my_rank==0)
    printf("saving optimal trajectory...");

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

  for(n=0;n<MAX_ITERATION;n++)
    delete[] ot[n];
  delete[] ot;
  delete[] u_star;
  delete[] time;
  delete[] position;
  delete[] ff;

  return success;
}
