#ifndef _COMMANDS_H
#define _COMMANDS_H
#include "stdafx.h"

//! comment to disable range checking
//#define RANGE_CHECK


class Commands {

  public:
  int                       dimC;   /*! commands dimension */
  int                       ncall;  /*! total number of commands */
  int                       ncall2; /*! in case of second set of commands, total number of commands for that last one */

  private:
  double**                  cdata;  /*! commands table */
  int**                     idata;  /*! commands coordinate table */
  int*                      nc;     /*! number of commands in each dimension */
  double*                   umin;   /*! lowest limits in each dimension */
  double*                   umax;   /*! highest limits in each dimension */

  public:

  /*!
  *  \brief default Constructor
  */
  Commands ()  {};

  Commands (struct cp commande_params);
  Commands (int dimensionC, const int* tailleC, const double* borneminC, const double* bornemaxC);

  /*!
   *  \brief default destructor
   */
  ~Commands();

  void init (const Commands &);
  void init (int dimensionC, const int* tailleC, const double* borneminC, const double* bornemaxC);
  void init (struct cp commande_params);
  int  init_CNES(int COMPUTE_CONTROL, int NCt);
  double& operator() (int, int);
  double  operator() (int, int) const;
  C operator[] (int) const;
  friend ostream& operator<<(ostream&, const Commands &);

  private:

  /*!
   *  \brief copy Constructor
   */
  Commands(const Commands &);
  /*!
   *  \brief copy assignment
   */
  void operator=(const Commands &);
  void range_check(int lig, int col) const;

};


inline double& Commands::operator() (int lig, int col)
{
  #ifdef RANGE_CHECK
    range_check(lig,col);
  #endif
  return cdata[lig][col];
}

inline double Commands::operator() (int lig, int col) const
{
  #ifdef RANGE_CHECK
    range_check(lig,col);
  #endif
  return cdata[lig][col];
}

inline C Commands::operator[] (int lig) const
{
  #ifdef RANGE_CHECK
    range_check(lig,1);
  #endif
  return cdata[lig];
}

#ifdef RANGE_CHECK
inline void Commands::range_check(int lig, int col) const
{
  if (lig >= (int)ncall || col >= (int)dimC)
    throw range_error("Matrix subscript out of bounds");
}
#endif


inline ostream& operator <<(ostream& f, const C & P )
{
  f << '\t';
  for(unsigned i=0;i<1;i++)
    f << P[i] << ' ';
  f << '\t';
  return f;
}

#endif /* Commands_H */

