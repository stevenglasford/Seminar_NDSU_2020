/*
 * Fichier utilisateur associe a la fonction de calcul de la trajectoire optimale avec un pas adaptatif
 *
 * Ce fichier contient deux fonctions :
 *   -une implementation de la fontion  de recherche de contrôle optimal, find_optimal_control_adapt()
 *   -une adaptation de l'integrateur de Heun pour plusieurs pas
 */

/*!
 *
 * \brief fonction de recherche du controle optimal
 *
 * @param[in] y : l'état en cours
 * @param[in] h : pas de progression unitaire
 * @param[in] t : temps en cours
 * @param[out] val: la valeur optimale trouvée
 * @param[in] nbSteps : nombre de pas unitaires utilisé pour l'intégration à partir de y avec  un contrôle donné
 * @param[in] vtab : la fonction valeur utilisée pour le calcul de la valeur optimale
 * @return le numéro du contrôle optimal trouvé
 */

int HJB::find_optimal_control_adapt(const double* y, double h, double t, double& val, int nbSteps, double *vtab)
{
  double val_star,vect[ncall];
  bool inTheBounds;

  int c;
  int u_star;
  double   rc;
  double dvect[dim];
  double currentVal;
  //int dimc=u->dimC;

  val_star=INF; u_star=-1;

  for(c=0;c<ncall;c++)
  {
    /*
     * pour tous les controles possibles
     * calcul de la positions suivante correspondante avec RK2
     */
    heunIntegrator(y, nbSteps, c, t, h, dvect, &inTheBounds);
    rc = (!inTheBounds) ? INF : (*this.*interpolation)(dvect,topt);
    if(rc<=2.0*T)
      vect[c]=rc;
    else
      vect[c]=INF;
  }

  for(c=0;c<ncall;c++)
  {
    if(vect[c]<INF)
    {
      currentVal=vect[c];

      if(currentVal<val_star)
      {
        val_star=currentVal;
        u_star=c;
      }
    }
  }

  val=val_star;
  return(u_star);
}

/*!
 *
 * \brief Cette fonction  réalise le calcul d'une intégration de numérique  de Heune sur un nombre donné de pas
 *
 * @param[in] yInit : l'état de début
 * @param[in] h : pas de progression unitaire
 * @param[in] tau : temps initial
 * @param[out] inTheBounds: indicateur booléen  de non dépassement des limites du domaine de calcul
 * @param[in] nbSteps : nombre de pas unitaires utilisé pour l'intégration à partir de yInit avec  un contrôle donné
 * @param[in] c : le numéro du contrôle à appliquer le long de tout le chemin
 * @param[out] res: l'état final (successeur de yInit pour le contrôle c )
 *
 */



void HJB::heunIntegrator(const double* yInit,int nbSteps, int c, double tau, double h, double *res, bool *inTheBounds)
{

  double ff2[dim], ff[dim],y[dim]; // tableaux pour le stockage temporaire des images intermediaires calculees
  // pour le schema de Heun
  int d;
  double tt=tau;                   // temps en cours pour le chemin a calculer
  int k=0;

  *inTheBounds=true;

  for(d=0;d<dim;d++)
    y[d]=yInit[d];

  while ((k<nbSteps )& (*inTheBounds))
  {

    (*dynamics)(y,(*u)[c],tt,ff2);
    for(d=0;d<dim;d++)
      res[d]=y[d]+h*ff2[d];
    (*dynamics)(res,(*u)[c],tt+h,ff);

    for(d=0;d<dim;d++)
    {
      res[d]=y[d]+0.5*h*(ff[d]+ff2[d]);
      y[d]=res[d];
    }

    if(periodic_mesh)
      periodizePoint(res);

    for(d=0; d<dim; d++)
    {
      if(res[d]<  mesh->lowBounds[d] || res[d]>mesh->highBounds[d])
      {
        (*inTheBounds)=false;
        break;
      }
    }

    k++;

    tt+=h;
  }

}


