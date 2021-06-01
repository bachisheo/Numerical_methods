//
// Created by kesa on 01.06.2021.
//

#ifndef NUMERICAL_METHODS_GEOMETRY_H
#define NUMERICAL_METHODS_GEOMETRY_H
#include "vect.h"
 struct vect2d {
    db x, y;
    vect2d (point<db> A, point<db> B);

};
db vectProd(vect2d a, vect2d b);

struct section{
    point<db> a, b;
    //check line segments AB and CD for intersection using vector product
    static bool intersect(section first, section second);
    section(point<db> A, point<db> B);

};


#endif //NUMERICAL_METHODS_GEOMETRY_H
