#ifndef _GEOMSTL_H
#define _GEOMSTL_H
#include <string>
#include <vector>
#include <cassert>
#include <fstream>
#include <iostream>
#include <sstream>
#include <streambuf>
#include "definitions.h"

/*!
 * \struct point
 * \brief Structure to hold coordinates of a point
 */

   struct point {
    float x; /*!< x coordinate */
    float y; ; /*!< y coordinate */
    float z; ; /*!< z coordinate */
    point() : x(0), y(0), z(0) {} /*!< default constructor*/
    point(float xp, float yp, float zp) : x(xp), y(yp), z(zp) {} /*!<constructor */
  };
/*!
 *\struct triangle
 * \brief structure to store normals and vertices of a triangle
 *
 */ 
  struct triangle {
    point normal; /*!< normal */ 
    point v1; /*!< coordinates of the vertex 1*/
    point v2; /*!< coordinates of the vertex 2*/
    point v3; /*!< coordinates of the vertex 3*/
    triangle(point normalp, point v1p, point v2p, point v3p) :
      normal(normalp), v1(v1p), v2(v2p), v3(v3p) {}
  };

//  std::ostream& operator<<(std::ostream& out, const triangle& t);

/*!
 *\struct stl_data
 * \brief structure to store vector of triangles read in from *.stl file
 *
 */ 
  struct stl_data {
    std::string name;
    std::vector<triangle> triangles;
    //triangle *triangles;
    stl_data(std::string namep) : name(namep) {}
  };

class Vec3 //Class for 3d Vector Magnitude and Normalization
{
public:
        Vec3(void);
        Vec3(float X, float Y, float Z);
        ~Vec3(void);
        float Length();
        Vec3 Normalize();
        Vec3 Vectors();
        float X, Y, Z;
};


class GeomSTL{

  private:
  int geom_nn;
  real *geom_xyz;
  real *__restrict__ triangle_center;
  real xyz[6];
  point centroid;
  real *distanceToBoundaries=nullptr;
  real *fitParlpd=nullptr;

  int nSTL;
  
  public:

  GeomSTL(){};
  void construct(int nSTL, real *xyz1);
  void construct(real *xyz1);

  stl_data parse_stl(const std::string& stl_path);
  float parse_float( std::ifstream &s );
  point parse_point( std::ifstream &s );
  stl_data parseSTL( const std::string &stl_path );
  void readSTLGeom(char *argv[], int n,  const real *xyz, const double init_loc[3]);
  void checkMesh( std::vector<triangle> &triangles );
  void updateDistanceToBoundaries();

  template <size_t N, typename Nvalue, size_t M, typename Mvalue >
  friend class ReblAmr;

  template <size_t N, typename Nvalue, size_t M, typename Mvalue >
  friend class ReblAmrFull;


  template <size_t N, typename Nvalue, size_t M, typename Mvalue, class T >
  friend class TemplateForest;

 ~GeomSTL();

};
 


#endif
