#ifndef _HJB_FD_H
#define _HJB_FD_H

#include "stdafx.h"
#include "HJB.h"

class HJB_FD : public HJB {

  private:
  int       type_rk;
  double    CFL;          /*! CFL condition */
  double    *Dvnum;       /*! used for function Hnum results  */
  double    *vold2;       /*! Used to stocked the previous iteration's value function [SIZE] */
  int       casedt;


  public:
  /*! constructeurs et destructeurs */

  /*!
  *  \brief default Constructor
  */
  HJB_FD(){};

  /*!
   *  \brief alternate Constructor
   *
   *  \param RegularMesh_params : includes all regular mesh parameters
   *  \param commande_params    : includes all commands/controls parameters
   *  \param commande_params2   : includes all second commands/controls parameters
   *  \param parallel_params    : includes all parallelization parameters
   *  \param general_params     : includes all general parameters for mainloop
   *  \param postprocess_params : includes all post processing parameters
   *  \param HJB_FD_params      : includes all Finite Difference method related parameters
   */
  HJB_FD(struct mrp RegularMesh_params, struct cp commande_params, struct cp commande_params2, struct pp parallel_params, struct gp general_params, struct ppp postprocess_params, struct fdp HJB_FD_params);

  /*!
   *  \brief default destructor
   */
  ~HJB_FD();

  private:
  /*!
   *  \brief copy Constructor
   */
  HJB_FD(const HJB_FD &);
  /*!
   *  \brief copy assignment
   */
  void operator=(const HJB_FD &);


  /*!
   *  \brief Lax-Friedrichs scheme
   *
   *  \param t  : current time in the mainloop algorithm
   *  \param dt : time step
   */
  void      LAFR                    (double t, double dt, double *vin, double *vout);
  void      LAFR_mpi                (double t, double dt, double *vin, double *vout);
  void      LAFR_omp                (double t, double dt, double *vin, double *vout);
  void      LAFR_mpi_omp            (double t, double dt, double *vin, double *vout);

  /*!
   *  \brief ENO2, RK1 scheme (euler forward/explicite)
   *
   *  \param t  : current time in the mainloop algorithm
   *  \param dt : time step
   */
  void      ENO2_RK1                (double t, double dt, double *vin, double *vout);
  void      ENO2_RK1_mpi            (double t, double dt, double *vin, double *vout);
  void      ENO2_RK1_omp            (double t, double dt, double *vin, double *vout);
  void      ENO2_RK1_mpi_omp        (double t, double dt, double *vin, double *vout);


  /*!
   *  \brief ENO2, RK2 scheme (euler forward/explicite)
   *
   *  \param t  : current time in the mainloop algorithm
   *  \param dt : time step
   */
  void      ENO2_RK2                (double t, double dt, double *vin, double *vout);
  void      ENO2_RK2_mpi            (double t, double dt, double *vin, double *vout);
  void      ENO2_RK2_omp            (double t, double dt, double *vin, double *vout);
  void      ENO2_RK2_mpi_omp        (double t, double dt, double *vin, double *vout);
  /*!
   *  \brief ENO2, RK3 scheme (euler forward/explicite)
   *
   *  \param t  : current time in the mainloop algorithm
   *  \param dt : time step
   *  \param vf : vector to backup the value function
   */
  void      ENO2_RK3                (double t, double dt, double *vin, double *vout);
  void      ENO2_RK3_mpi            (double t, double dt, double *vin, double *vout);
  void      ENO2_RK3_omp            (double t, double dt, double *vin, double *vout);
  void      ENO2_RK3_mpi_omp        (double t, double dt, double *vin, double *vout);

  /*!
   *  \brief function pointer where one Finite Difference iteration function is linked with
   *
   *  \param t  : current time in the mainloop algorithm
   *  \param dt : time step
   */
  void      (HJB_FD::*PROC)             (double t, double dt, double *vin, double *vout);

  /*!
   *  \brief during the mainloop, update the time step DT for the case with a non command problem and where the user has not given a specific DT
   *
   *  \param t  : current time in the mainloop algorithm
   *  \param dt : time step
   */
  double    updateDT_case0c0        (double t, double dt);

  /*!
   *  \brief during the mainloop, update the time step DT for the case with a one command problem and where the user has not given a specific DT
   *
   *  \param t  : current time in the mainloop algorithm
   *  \param dt : time step
   */
  double    updateDT_case0c1        (double t, double dt);
  double    updateDT_case0c1_mpi    (double t, double dt);
  double    updateDT_case0c1_omp0   (double t, double dt);
  double    updateDT_case0c1_omp    (double t, double dt);
  double    updateDT_case0c1_mpi_omp(double t, double dt);


  /*!
   *  \brief during the mainloop, update the time step DT for the case with a two commands problem and where the user has not given a specific DT
   *
   *  \param t  : current time in the mainloop algorithm
   *  \param dt : time step
   */
  double    updateDT_case0c2        (double t, double dt);
  double    updateDT_case0c2_omp    (double t, double dt);
  double    updateDT_case0c2_mpi    (double t, double dt);
  double    updateDT_case0c2_mpi_omp(double t, double dt);

    /*!
   *  \brief during the mainloop, update the time step DT for the case where the user has given a specific DT
   *
   *  \param t  : current time in the mainloop algorithm
   *  \param dt : time step
   */
  double    updateDT_case1          (double t, double dt);


  /*!
   *  \brief function pointer where one update time step function is linked with
   *
   *  \param t  : current time in the mainloop algorithm
   *  \param dt : time step
   */
  double    (HJB_FD::*UPDATE_DT)    (double t, double dt);

  public:

  /*!
   *  \brief execute the main algorithm of the chosen method
   */
  virtual void mainloop             ();


  private:
  double    (*HnumFunc_0C) (const double, const double*, const double, const double*);
  double    (HJB_FD::*Hnum)(const double, const double*, const double, const double*);
  double    Hnum_0C        (const double, const double*, const double, const double*);
  double    Hnum_0C_MaxMin (const double, const double*, const double, const double*);
  double    Hnum_0C_MinMax (const double, const double*, const double, const double*);
  double    Hnum_1C_Max    (const double, const double*, const double, const double*);
  double    Hnum_1C_Max_omp(const double, const double*, const double, const double*);
  double    Hnum_1C_Min    (const double, const double*, const double, const double*);
  double    Hnum_1C_Min_omp(const double, const double*, const double, const double*);
  double    Hnum_1C_MaxMin (const double, const double*, const double, const double*);
  double    Hnum_1C_MinMax (const double, const double*, const double, const double*);
  double    Hnum_2C_MaxMin (const double, const double*, const double, const double*);
  double    Hnum_2C_MinMax (const double, const double*, const double, const double*);
  void      (*computeamax) (double*, double);

  //- MAY 2013: for steady equations
  double    (*discount_factor)      (const double*);

  /*!
   *  \brief minmod function
   *
   *  \param t  : current time in the mainloop algorithm
   *  \param dt : time step
   */
  double    minmod                  (double x, double y);


  //- 2016
  int * borderSizeDims; // axis numbers for which the borderSize[d] is not nul.
  int   nbBorderDims;   // number of state variables d for which borderSize[d] is non zero.

};

#endif /* HJB_FD_H */
