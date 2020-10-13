//---------------------------------
//- optimal trajectory for a given point and output in a  (.dat or VTK) format 
//- BASED ON tmin, adaptive case
//---------------------------------

/*!
 * \brief Reconstruction de trajectoire optimale basee sur une fonction  valeur , version avec pas adaptatif
 *
 * @param initialpoint[in] :etat initial pur la trajectoire
 * @param filename[in] : nom de fichier pour enregistrer la trajectoire
 * @param start
 * @param terminalpoint[out] : etat final de la trajectoire
 * @param terminaltime[out] : temps final de la trajectoire
 * @return succes: indicateur de succes pour l'atteinte ou non da la cible
 */

int HJB::compute_optimal_trajectory_texit(const double* initialpoint, const char* filename, int start, double * terminalpoint , double &terminaltime)
{

	cout<< "temps de sortie optim traj\n";

	double *time, **ot, *position;  // tableaux pour le stockage temporaire
	int n, *u_star;
	time   = new double [MAX_ITERATION];
	ot     = new double*[MAX_ITERATION];
	u_star = new int    [MAX_ITERATION];

	for(n=0;n<MAX_ITERATION;n++)
		ot[n] = new double[dim];

	int PRINTTRAJ=0;

	//  differents  pas de temps utilises
	double dt0, dt;
	//indicateur de succes d'atteinte de cible
	int succes=0;


	double *ff,*ff2,*maxConst;
	//double treste;
	ff = new double [dim];
	ff2 = new double [dim];
	maxConst = new double [dim];
	int   c, d, k, end=0;
	double h=0., hx,   support, maximum, val;


	//- stopping threshold
	//double eps=0.001; //- OB: commented because not really used

	FILE * pFile;


	position = new double[dim];

	for(d=0; d<dim-1; d++)
	{
		position[d]=initialpoint[d];
	}
	n=0;
	for(d=0; d<dim; d++)
	{
		ot[n][d]=position[d];
		maxConst[d]=0.0;                //starting point (t=0)
	}

	time[n]=0.;


/*
 * Cette partie est specifique aux problemes de type d+1
 * avec la dernire variable qui est auxiliaire
 * Il s'agit de  z au dessus de l'epigraphe de v(t,x)
 *
 * Ici on recherche la valeur optimale de z
 * pour laquelle la fonction  de temps
 * max est egale a T ( ou proche)
 */

	double z=mesh->lowBounds[dim-1], dz=0.05;
	bool initPointTest=false;
	while(!initPointTest & (z<mesh->highBounds[dim-1]))
	{
		position[dim-1]=z;
    	cout<< " z= "<< z << " t exit = " << (*this.*interpolation)(position,topt)<<endl;
		initPointTest=((*this.*interpolation)(position,topt)>=T);
		z+=dz;
	}

	if(!initPointTest)
	{
     /*
      * Si un tel z n'a pas ete trouve  il n'est pas possible de construire la trajectoire
      */
		cout<< " pas de z trouve pour la position initiale donc pas de trajectoire \n";
	}
	else
	{
		z-=dz;
		position[dim-1]=z;
		cout<< " z trouve "<<z<<endl;
		cout<< " temps interpole "<<(*this.*interpolation)(position,topt)<<endl;

		int uOpt;
		//double tOpt=INF; //- OB: not used 
		double hOpt;
		bool inTheBounds;
		double tReel=0.0;
		//double tEnCours;
		int nbSteps,nbOpt=1;
		//int nbS;
		//bool testTime; //- OB: commented because not really used
		//int iter=0;    //- OB: commented because not really used

		for(d=0; d<dim; d++)
			{
				ot[n][d]=position[d];
				maxConst[d]=0.0;                //starting point (t=0)
			}

		/*
		 * dt0 est le pas de temps de base pour l'integrateur
		 */
		/*!
		 * @todo :  a revoir cette conception du pas de temps de base
		 * variante 1, un peu testee: utiliser une determination de type CFL
		 * variante 2 : sortir cela comme parametre utilisateur
		 */
		dt0= 0.01;

		while (!end && (tReel<T)&& (n<MAX_ITERATION-1))
		{

			/*!
			 * calcul de la vitesse locale pour determiner le pas  de temps
			 */
			maximum=0.;

			for(c=0;c<ncall;c++)
			{
				(*dynamics)(ot[n],(*u)[c],tReel,ff);

				for(d=0;d<dim;d++)
				{
					support=ABS(ff[d]);
					if(support>maximum)  maximum=support;
				}
			}


			/*!
			 * On calcule  le pas  de temps
			 *  en prenant en compte  la norme de la dynamique  au point  courent de
			 *  la trajectoire;
			 *  tReel est le temps  de la trajectoire
			 *   h est le pas calcule a partir de l'estimation de la vitesse locale
			 *   et du pas de la grille spatiale
			 *   On devra calculer le successeur du point en cours sur la trajectoire
			 *   pour le temps t+h
			 *   Pour ce faire on determine avec quel pas d'integration dt on va integrer
			 *   les equations differentielles entre t et t+h
			 *   Si h < dt0 on prendra pour la pas d'avancement dt=0.5*h
			 *   Si h>dt0 on prendra comme pas dt=dt0
			 *   Le nombre de pas est determine de faÃ§on a realiser h avec les pas d'integration dt
			 */

			if(maximum>0)
			{
				h= min(T-tReel, (mesh->Dx_min/maximum)*0.5);

				if (h<dt0)
					dt=0.5*h;
				else
					dt=dt0;

				nbSteps=floor(h/dt);
			}
			else
			{
				h= min(T-tReel,10.0*dt0);
				if (h<dt0)
					dt=0.99*h;
				else
					dt=dt0;
				nbSteps=floor(h/dt);
				end=1;
				break;
			}

			// on calcule hx l'accroissement du temps de la trajectoire
			//en fonction du pas choisi et du nombre de pas a realiser
			hx=dt*((double) nbSteps);

			// le temps du successeur

			time[n+1]=time[n]+hx;

/*
 * On appelle la recherche de controle optimal qui
 * utilise juste l'integrateur de Heune avec le nombre precise de pas de base
 */
			uOpt=find_optimal_control_texit(ot[n],dt,tReel,val,nbSteps,topt);
			//- OB test:
			if (uOpt==-1) {
				printf("(! uopt=-1 => break)");
				printf("(=>break)");
				break;
			}
			hOpt=hx;
			nbOpt=nbSteps;
			//nbS=nbSteps; //- OB: not used

			u_star[n]=uOpt;
			hx=hOpt;
			nbSteps=nbOpt;
			time[n+1]=time[n]+hx;

			heunIntegrator(ot[n], nbSteps,  u_star[n], tReel, dt, ot[n+1],  &inTheBounds);

                        /*
			if(inTheBounds)
			  treste=(*this.*interpolation)(ot[n+1],topt);
			else
			  treste=INF;
                        */

			/*for(d=0;d<dim;d++)
				printf("x%i=%6.3f, ",d,ot[n][d]);
			//printf("iter n=%3i, time[n]=%5.2f, h=%6.3f, iu=%3i, u=%5.2f, tmin=%5.3f, val=%5.3f\n",
			//  n,time[n],hx,u_star[n],(*u)[u_star[n]][1],treste,val);
			printf("iter n=%3i, time[n]=%6.3f h=%6.3f, ",n,time[n],hx);

			for(int i=0;i<(*u).dimC;i++)
				printf("u%i=%6.3f, ",i,(*u)[u_star[n]][i]);

			printf("tmin=%5.3f, ", (*this.*interpolation)(ot[n],topt));

			printf("t=%5.3f (traj)", t);

			printf("\n");*/

			n++;

			tReel=min(tReel+hx,T);

			if(!inTheBounds)
			{
				end=1;
			}

			if(val>=INF)
			{
				end=1;
			}
			/*
			 *  test de cible
			 */

			if( (time[n]>=T))
			{
				//printf("Target reached."); //- removed since "success=1" will appear
				end=1;
				succes= 1;

			}


			if (PRINTBUG)
				printf("end of step n=%4i ...\n",n);
		}


	}

	//cout<< " final position =";

	for(d=0; d<dim; d++)
	{
		position[d]=ot[n][d];
		terminalpoint[d]=ot[n][d];
	}
	terminaltime=time[n];


	u_star[n]=u_star[n-1];

	printf("(terminaltime=%8.5f, success=%1i)\n",terminaltime,succes);


	//- save trajectory
	if (PRINTTRAJ && my_rank==0)
		printf("saving optimal trajectory...");

	if (my_rank==0){
		pFile = fopen (filename,"w");
		if(pFile==NULL){
			cout << "19xx.Impossible to write in the given file!" << endl;
			exit(1);
		}
		else{
			for (k=0; k<=n; k++)  {
				for (d=0; d<dim; d++)
					fprintf(pFile,"%8.5f ", ot[k][d]);
				fprintf(pFile,"     ");
				for (int i=0; i<(*u).dimC; i++)
					fprintf(pFile,"%8.5f ", (*u)[u_star[k]][i]); //- control value
				fprintf(pFile,"     ");
				fprintf(pFile,"%8.5f ", time[k]); //- time
				fprintf(pFile,"\n");
			}
		}
		fclose (pFile);

		if (PRINTTRAJ && my_rank==0)
			printf("done.\n");
	}

	//- save VTK files
	//cout << "saving VTK files...";
	//cout << "(VTK/tabxx)...";
	//for(k=0; k<=n; k++){
	//  saveVTK2D("VTK/tab",k,ot[k][dim-1]);
	//}
	//cout << "(VTK/trajxx)...";
	//for(k=0; k<=n; k++){
	//  saveTrajCAR(filename, k, k, 20, ot);
	//}
	//cout << "done." << endl;
	//- end

	for(n=0;n<MAX_ITERATION;n++)
		delete[] ot[n];
	delete[] ot;
	delete[] u_star;
	delete[] time;
	delete[] position;
	delete[] ff;
	delete[] ff2;

	return succes;
}
