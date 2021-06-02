//
// Created by kesa on 01.06.2021.
//

#ifndef NUMERICAL_METHODS_GEOMETRY_H
#define NUMERICAL_METHODS_GEOMETRY_H

#include "vect.h"

template<typename T>
struct point {
public:
    T x = 0, y = 0;

    point(T x, T y);

    point() = default;
};

template<typename T>
std::ostream &operator<<(ostream &fout, vector<point<T>> points);

void printPoints (ostream& out, vector<point<db>> points);
struct vect2d {
    db x, y;
    vect2d(point<db> A, point<db> B);
};

db vectProd(vect2d a, vect2d b);

struct section {
    point<db> a, b;

    //check line segments AB and CD for intersection using vector product
    static bool intersect(section first, section second);

    section(point<db> A, point<db> B);

};


///DEFINITION OF TEMPLATE METHODS
template<typename T>
point<T>::point(T x, T y) : x(x), y(y) {}

template<typename T>
std::ostream &operator<<(ostream &fout, vector<point<T>> points) {
    for (int i = 0; i < points.size(); i++)
        fout << points[i].x << " ";
    fout << endl;
    for (int i = 0; i < points.size(); i++)
        fout << points[i].y << " ";
    fout << endl;
    return fout;
}
template <typename T>
void printPoints(ostream &out, vector<point<T>> points) {
    for (int i = 0; i < points.size(); ++i)
        out << points[i].x << " ";
    out << endl;
    for (int i = 0; i < points.size(); ++i)
        out << points[i].y << " ";
    out << endl;
}

#endif //NUMERICAL_METHODS_GEOMETRY_H
