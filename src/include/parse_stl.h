#ifndef PARSE_STL_H
#define PARSE_STL_H

#include <string>
#include <vector>
#include <cassert>
#include <fstream>
#include <iostream>
#include <sstream>
#include <streambuf>

namespace stl {

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

  std::ostream& operator<<(std::ostream& out, const triangle& t);

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

  stl_data parse_stl(const std::string& stl_path);
  //void parse_stl(const std::string& stl_path,stl_data& info);

}

class Vector3 //Class for 3d Vector Magnitude and Normalization
{
public:
        Vector3(void);
        Vector3(float X, float Y, float Z);
        ~Vector3(void);
        float Length();
        Vector3 Normalize();
        Vector3 Vectors();
        float X, Y, Z;
};

void checkMesh(std::vector<stl::triangle> &triangles);





#endif
