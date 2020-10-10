/*
 * Fichier utilisateur associe  a la fonction de calcul de la trajectoire optimale
 * avec un pas adaptatif
 *
 * Ce fihcier contient deux fonctions :
 *
 *   -une implementation de la fontion  de recherche de contrele optimal, find_optimal_control_adapt()
 *   -une adaptation de l'integrateur de Heun pour plusieurs pas
 */
/*!
 *
 * \brief Cette fonction  recherche le contrele optimal
 *
 * @param[in] y : l'etat en cours
 * @param[in] h : pas de progression unitaire
 * @param[in] t : temps en cours
 * @param[out] val: la valeur optimale trouvee
 * @param[in] nbSteps : nombre de pas unitaires utilise pour l'integration a partir de y avec  un contrele donne
 * @param[in] vtab : la fonction valeur utilisee pour le calcul de la valeur optimale
 * @return le numero du contrele optimal trouve
 */
int HJB::find_optimal_control_texit(const double* y, double h, double t, double& val, int nbSteps, double *vtab)
{
	double val_star,*vector;
			bool inTheBounds;
			//double Fnorm;
			double fact1=1.0, fact3=0.0, rc;
			vector = new double[ncall];
			int c;//d; //test;
			int u_star;
			double t1; //rc_star=T;
			//double ff[dim], ff2[dim];
			double dvect[dim];

			val_star=INF; u_star=0;
			t1 =    (*this.*interpolation)(y,topt) ;

			//cout<< "   find control  t="<<t<<" t1= "<<t1<<endl;
			for(c=0;c<ncall;c++)
			{
				/*
				 * pour tous les conroles possibles
				 * calcul   de la positions suivante correspondante avec RK2
				 */


					heunIntegrator( y,nbSteps,  c,   t,   h, dvect, &inTheBounds);

					/*
					 * si le successeur est dans le domaine alors
					 *  on calcule la fonction tmin correspondante
					 *  par interpolation
					 */
					rc = (!inTheBounds) ? INF :  (*this.*interpolation)(dvect,topt) ;

					/*
					 *   si la valeur interpolée est nulle c'est qu'on est hors de
					 *   contraintes; on ne prend donc pas en compte ce contrôle
					 */
					if(rc>0)
					{

						 /*
						  *  Le contrôle optimal recherché maximise
						  *  la fonction de temps de sortie des successeurs;
						  *  il doit également correspondre  au principe de programmation dynamique qui
						  *  assure une sorte de continuité de la foncton de temps le long
						  *  de la trajectoire
						  *  c'est pourquoi  on peut minimiser ici
						  *  une fonction composite entre le max de t(x) ( min de (T-t(x)) et
						  *  la norme de l'erreur du principe de programmation dyn
						  */
					vector[c] = fact1*(T-rc) + fact3* (min(T,rc+nbSteps*h)-t1)* (min(T,rc+nbSteps*h)-t1);

						if(vector[c]<val_star){
							val_star=vector[c];
							u_star=c;
						}
					}

				}

			if(val_star>=INF){
				cout << "no admissible controls!!!" << endl;
			}
			val=val_star;
			return(u_star);
}
