#ifndef _HJB_SL_H
#define _HJB_SL_H

#include "stdafx.h"
#include "HJB.h"

class HJB_SL : public HJB {

  private:
  /*! \brief local HJB_SL variables */
  int     order;
  double  *tmpcoords;
  double  *ff;
  double  *ff2;
  double  *position;
  double  *hi;
  double  hmax;
  double  facteur1;
  double  *lb,*hb;
  int     interpolation_param;
  int     p_intermediate;
  double  beta_coef;
  double  P;


  /*! \brief  for the stationnary case problems, the 'compute' table defines the mesh points
   *          to be updated in process of the mainloop
   *
   *          for a given point i,
   *            compute[i]=0 if v0(i)<0
   *            compute[i]=1 elsewise
   */
  int     *compute;
  int     **rang;
  int     it_step;

  /*! \brief  'INTERPOLATION_PARAM' indicates if we use the basic bilinear interpolation
   *          or the direction per direction interolation
   *          if the bilinear interpolation is chosen, it tests if each point and each 
   *          command have to be precomputed and saved in 'tabPP' & 'tabRK'
   */
  double  ***tabPP;
  int     **tabRK;

  /*! \brief OpenMP pointer variables */
  double **rCoordomp;
  double **ffomp;
  double **ff2omp;
  double **dvectdoubleomp;


  public:
  /*!
  *  \brief default Constructor
  */
  HJB_SL(){};


  /*!
   *  \brief alternate Constructor
   *
   *  \param RegularMesh_params : includes all regular mesh parameters
   *  \param commande_params    : includes all commands/controls parameters
   *  \param commande_params2    : includes all second commands/controls parameters
   *  \param parallel_params    : includes all parallelization parameters
   *  \param general_params     : includes all general parameters for mainloop
   *  \param postprocess_params : includes all post processing parameters
   *  \param HJB_SL_params          : includes all Semi Lagrangian Method related parameters
   */
  HJB_SL(struct mrp maillage_regulier_params, struct cp commande_params, struct cp commande_params2, struct pp parallel_params, struct gp general_params, struct ppp postprocess_params, struct slp HJB_SL_params);

  /*!
   *  \brief default destructor
   */
  ~HJB_SL();

  private:

  /*!
   *  \brief copy Constructor
   */
  HJB_SL(const HJB_SL &);
  /*!
   *  \brief copy assignment
   */
  void operator=(const HJB_SL &);


  public:
  virtual void mainloop();

  private:

  /*!
   *  \brief MIN or MAX for optimisation (inverse of MAX / MIN from HJB::OPTIM()
   *  \param a : first  double (to be compared with parameter b)
   *  \param b : second double (to be compared with parameter a)
   */
  double    (*OPTIM_SL)                    (double a, double b);

  /*!
   *  \brief Semi Lagrangian iteration function for the stationnary case
   *
   *         Cadre u + H(x,Nabla u) =0 avec
   *         H(x,nabla u) == max_a (-f(x,a) . nabla u - L(x,a));
   *         Schema : iteration sur n de :
   *            u^{n+1}(x)  + max_a ( (u^{n+1}(x)- u^n(x+h f(x,a))/h  - L(x,a)) = 0
   *
   *         avec h=DT, on obtient
   *            u^{n+1}(x) =  min_a ( 1/(1+h) [u^n](x+h f(x,a)) + h L(x,a) )
   *
   *  \param t  : current time in the algorithm run
   *  \param dt : current time step in the algorithm run
   */
  void      itSL_stat         (double t, double dt);
  void      itSL_stat0        (double t, double dt);
  void      itSL_stat_omp     (double t, double dt);

  /*!
   *  \brief Semi Lagrangian 2xDIM iterations function for the stationnary case with 2 mesh loops
   *         in forward and backward for each dimension
   *
   *  \param t  : current time in the algorithm run
   *  \param dt : current time step in the algorithm run
   */
  void      itSL_stat2        (double t, double dt);
  void      itSL_stat2_omp    (double t, double dt);

  /*!
   *  \brief Semi Lagrangian iteration function for the dynamic case
   *
   *  \param t  : current time in the algorithm run
   *  \param dt : current time step in the algorithm run
   */
  void      itSL_evo          (double t, double dt);
  void      itSL_evo_omp      (double t, double dt);
  void      itSL_evo_omp2     (double t, double dt);


  /*!
   *  \brief Second order Semi Lagrangian iteration function for the dynamic case
   *
   *  \param t  : current time in the algorithm run
   *  \param dt : current time step in the algorithm run
   */
   void   secondorder_itSL_evo    (double t, double dt);
   void   secondorder_itSL_evo_omp(double t, double dt);


  /*!
   *  \brief Semi Lagrangian iteration function for the dynamic case with direction per direction interpolation
   *
   *  \param t  : current time in the algorithm run
   *  \param dt : current time step in the algorithm run
   */
  void      itSL_evo_dpd      (double t, double dt);

  /*!
   *  \brief Semi Lagrangian iteration function for the dynamic case where the bilinear interpolation for each
   *          point and for each command are precomputed and saved
   *
   *  \param t  : current time in the algorithm run
   *  \param dt : current time step in the algorithm run
   */
  void      itSL_evo_saved    (double t, double dt);
  void      itSL_evo_saved_omp(double t, double dt);

  /*!
   *  \brief function pointer where one the Semi Lagrangian iteration function is linked with
   *
   *  \param t  : current time in the algorithm run
   *  \param dt : current time step in the algorithm run
   */
  void      (HJB_SL::*PROC)   (double t, double dt);

  /*!
   *  \brief First order Rung-Kutta Euler method
   *
   *  \param h : space step
   *  \param t : current time in the algorithm run
   *  \param c : current used commmand ID
   */
  void      euler             (const double* coor, double* res, double h, double t, int c);

  /*!
   *  \brief Second order Rung-Kutta Heun method
   *
   *  \param h : space step
   *  \param t : current time in the algorithm run
   *  \param c : current used commmand ID
   */
  void      heun              (const double* coor, double* res, double h, double t, int c);

  /*!
   *  \brief Second order Rung-Kutta midpoint method
   *
   *  \param h : space step
   *  \param t : current time in the algorithm run
   *  \param c : current used commmand ID
   */
  void      midpoint          (const double* coor, double* res, double h, double t, int c);

  /*!
   *  \brief function pointer where one the Rung-Kutta method function is linked with
   */
  void      (HJB_SL::*METHODE)(const double*, double*,double,double,int);
  void      (HJB_SL::*METHODE_OMP)(const double*,double*,double*,double,double,int,double*);  // open MP equivalent


  double    (*funcR)          (const double* x, C u, double t);
  void      (*funcY)          (const double* x, int k, double eps, C u, double t, double h, double* res);

  //- MAY 2013: for steady equations
  double    (*discount_factor)(const double* x);

  void      euler2            (const double*,double*,double*,double,double,int,double*);
  void      heun2             (const double*,double*,double*,double,double,int,double*);
  void      midpoint2         (const double*,double*,double*,double,double,int,double*);


  double    interpolationBilineaire2d_saved (const double*, int, const double*);
  double    interpolationBilineaire3d_saved (const double*, int, const double*);
  double    interpolationBilineaire4d_saved (const double*, int, const double*);
  double    interpolationBilineaire5d_saved (const double*, int, const double*);

  double    (HJB_SL::*IB_saved) (const double*, int, const double*);
};


#endif /* HJB_SL_H */

