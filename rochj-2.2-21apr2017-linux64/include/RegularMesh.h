#ifndef REGULAR_MESH_H
#define REGULAR_MESH_H

#include "stdafx.h"

class RegularMesh {
  public:
  int       nvallblock;         //- boka
  int       dim;                /*! problem dimension */
  int       PRECOMPUTE_COORDS;  /*! test if mesh point coordinates have to be precomputed and saved before the mainloop call */
  int       MESH_TYPE;          /*! if 0 : mesh points are the center of mesh cells ; 1 : otherwise and mesh's global grid size are incremented in each direction */
  int       *borderSize;        /*! ghostnodes length for every dimension */
  int       rkfirst;            /*! Mesh point's rank with integer coordinates at [0,0,..,0] in augmentated mesh */

  int       inn_nbPoints;       /*! number of points in non decomposed mesh */
  int       out_nbPoints;       /*! number of points in mesh */

  int       *periodic;          /*! mesh periodic test in each direction */

  int       *GLOBAL_gridDim;    /*! number of points in each direction in non decomposed mesh */
  int       *inn_gridDim;       /*! number of points in each direction */
  int       *out_gridDim;       /*! number of points in each direction in augmentated mesh */
  double    *GLOBAL_lengths;

  double    *GLOBAL_lowBounds;  /*! lowest  limits in each direction in non decomposed mesh */
  //double    *GENERAL_lowBounds; /*! lowest  limits in each direction in non decomposed mesh discretization with acknowledgement of 'MESH_TYPE' parameter (option for centered mesh cells) */
  double    *GLOBALMESH_lowBounds; /*! lowest  limits in each direction in non decomposed mesh discretization with acknowledgement of 'MESH_TYPE' parameter (option for centered mesh cells) */
  double    *GLOBAL_highBounds; /*! highest limits in each direction in non decomposed mesh */
  double    *lowBounds;         /*! lowest  limits in each direction */
  double    *highBounds;        /*! highest limits in each direction */
  int       *coord_start;       /*! integer coordinates of the first point in the mesh */
  int       *coord_end;         /*! integer coordinates of the last  point in the mesh */

  double    *Dx;                /*! stride in each direction */
  double    Dx_min;             /*! minimal stride */

  int       *GLOBAL_neighbors;  /*! mesh steps to get to next or previous neighbor's rank in each direction in non decomposed mesh */
  int       *inn_neighbors;     /*! mesh steps to get to next or previous neighbor's rank in each direction */
  int       *out_neighbors;     /*! mesh steps to get to next or previous neighbor's rank in each direction in augmented mesh */

  double    *returncoor;        /*! output for function getcoords call */

  double    **discretization;   /*!real coordinates of discretization nodes in each dimension */

  protected:
  double    **coordinates;      /*! real coordinates for each point in the mesh */

  /* OBSOLETE
  int       combinaison_number;
  int       **combinaison_table;
  */

  double    mj;
  int       mi;

  private:

  /*!
   *  \brief copy Constructor
   */
  RegularMesh(const RegularMesh &);
    /*!
   *  \brief copy assignment
   */
  void operator=(const RegularMesh &);



  public:

  /*!
  *  \brief default Constructor
  */
  RegularMesh() {};

  /*!
   *  \brief alternate Constructor
   *
   *  \param DIM                : problem dimension
   *  \param RegularMesh_params : includes all regular mesh parameters
   */
  RegularMesh(int DIM, struct mrp  RegularMesh_params);


  virtual void setNewBoundaries(double* XMIN, double* XMAX);

  /*!
   *  \brief default destructor
   */
  virtual ~RegularMesh();

  /*!
   *  \brief evaluate the distance between two points
   *
   *  \param p1 : index of the first point
   *  \param p2 : index of the second one
   *  \return the distance
   */
  double getDistance(int p1, int p2)const;

  /*!
   *  \brief initialize the coordinates of each point
   */
  void initializeCoordinates(int* start, int* end);

  /*!
   *  \brief get all the neighbors of a particular point
   *
   *  \param i : index of teh particular point
   *  \param neighbors[2*DIM] : indices of all the neighbors
   *      return -1 if the neighbor doesn't exit (bound)
   */
  void getNeighbors(int i,int* neighbors);

  /*!
   *  \brief get all the neighbors gridvector distance from a point
   *
   *  \param i : index of teh particular point
   *  \param neighbors[3^DIM-1] : gridvector distance of all the neighbors
   */
  void getNeighbors(int** neighbors);


  /*!
   *  \brief check if the coordinate belongs to the domain
   *
   *  \param pos: real coordinates of the point
   *  \return true if in the domain, false otherwise
   */
  bool isInDomain(double* pos);
  bool isInDomain(int* pos);

  /*!
   *  \brief function pointer where one compute real coordinates function is linked with
   *  \param i: rank in the domain
   *  \return the real coordinates of the rank in the mesh
   */
  double*   (RegularMesh::*getcoords)      (int);
  /*!
   *  \brief function pointer where one compute real coordinates function is linked with
   *  \param i: rank in the domain
   *  \param coords: real coordinates of the rank in the mesh
   */
  void      (RegularMesh::*setcoords)      (int,double*);


  /*!
   *  \brief compute the rank of the coordinates (as argument) in the domain
   *
   *  \param pos: real coordinates of the point
   *  \return the rank of the point in the domain
   */
  int       getIndex                (const int* pos) const;


  /*!
   *  \brief compute the real coordinates of the rank (as argument) in the mesh using calculation
   *
   *  \param i: rank in the domain
   *  \return the real coordinates of the rank in the mesh
   */
  double*   getPositionCompute      (int i);


  /*!
   *  \brief compute the real coordinates of the rank (as argument) in the mesh using calculation and stock it the second argument
   *
   *  \param i: rank in the domain
   *  \param coords: real coordinates of the rank in the mesh
   */
  void      setPositionCompute      (int i, double * coords);


  /*!
   *  \brief compute the real coordinates of the rank (as argument) in the mesh using the precomputed coordinates table
   *
   *  \param i: rank in the domain
   *  \return the real coordinates of the rank in the mesh
   */
  double*   getPositionTable        (int i);


  /*!
   *  \brief compute the real coordinates of the rank (as argument) in the mesh using the precomputed coordinates table and stock it the second argument
   *
   *  \param i: rank in the domain
   *  \param coords: real coordinates of the rank in the mesh

   */
  void      setPositionTable        (int i, double * coords);


  /*!
   *  \brief compute the integer coordinates of the rank (as argument) in the mesh using the precomputed coordinates table and stock it the second argument
   *
   *  \param i: rank in the domain
   */
  void      setrank                 (int i, int *returnrank);

  /*!
   *  \brief compute the rank in the augmented mesh of the rank (as argument) in the mesh
   *
   *  \param i: rank in the domain
   */
  int       getOuterRank            (int i);


  /* OBSOLETE
  void      setPositionBlocTri      (int i, int * coords);
  double    distance                (const double* p1, const double* p2) const;
  void      getBarycenterVals       (const double *coords, int * ranks, double * weights);
  void      getBarycenterVals2d     (int rk, const double *coords, int * ranks, double * weights);
  void      getBarycenterVals3d     (int rk, const double *coords, int * ranks, double * weights);

  bool      isInTriangle            (const double *p, int i);
  bool      isInTetrahedron         (const double *p, int i);
  bool      isInTri4d               (const double *p, int i);
  bool      isInTri5d               (const double *p, int i);
  bool      (RegularMesh::*isInTris)       (const double *p, int i);
  */

  friend ostream& operator<<(ostream&, const RegularMesh &);

};

/* OBSOLETE
inline double area(const double *p1, const double *p2, const double *p3)
{
  return abs(p1[0]*p2[1]+p2[0]*p3[1]+p3[0]*p1[1]-p1[0]*p3[1]-p3[0]*p2[1]-p2[0]*p1[1])*.5;
}
*/

#endif //REGULAR_MESH_H
