//
// Created by Olya on 01.06.2021.
//

#include "geometry.h"

bool section::intersect(section AB, section CD) {
    vect2d ab(AB.a, AB.b), ad(AB.a, CD.b), ac(AB.a, CD.a);
    if(vectProd(ab, ad) * vectProd(ab, ac) > 0)
        return false;

    vect2d cd(CD.a, CD.b), ca(CD.a, AB.a), cb(CD.a, CD.b);
    if(vectProd(cd, ca) * vectProd(cd, cb) > 0)
        return false;
    return true;
}

section::section(point<db> A, point<db> B) :a(A), b(B){}

db vectProd(vect2d a, vect2d b) {
    return a.x * b.y - a.y * b.x;
}

vect2d::vect2d(point<db> A, point<db> B) {
    x = B.x - A.x;
    y = B.y - A.y;
}
