//
// Created on 29.05.2021.
//
#ifndef NUMERICAL_METHODS_VECT_H
#define NUMERICAL_METHODS_VECT_H
#include <cassert>
#include <vector>
#include <iostream>
#define ttt template <typename Type>
#define EPS 1e-6

typedef double db;

using namespace std;

template<typename Type>
class Point {
public:
    Type x = 0, y = 0;
    Point(Type x, Type y) :x(x), y(y) {}
};

template<typename Type>
class vect : public std::vector<Type> {
    std::vector<Type> v;
public:
    vect<Type>(int size) {
        v = vector<Type>(size);
    }

    vect<Type>() :vector<Type>(0){}

    vect<Type>(vector<Type> oldVector) : v(oldVector) {}

    Type min() const;

    Type max() const;

    int size() const;

    vect<db> operator*(db mult) const;

    //скалярное произведение векторов
    static db scalarProd(vect<db> a, vect<db> b);

    //найти норму вектора
    static db getNorm(vect<db> a);

    //нормирование вектора
    static vect<db> normalize(vect<db> x);

    // норма разности двух векторов
    static db diffNorm(vector<db> a, vector<db> b);
};

ttt
std::ostream & operator << (ostream& fout, vector<Point<Type>> points) {
    for (int i = 0; i < points.size(); i++)
        fout << points[i].x << " ";
    fout << endl;
    for (int i = 0; i < points.size(); i++)
        fout << points[i].y << " ";
    fout << endl;
    return fout;
}

template<typename Type>
Type vect<Type>::min() const {
    assert(v.size() > 0);
    Type mnm = v[0];
    for (int i = 0; i < v.size(); ++i)
        if (mnm > v[i])
            mnm = v[i];
    return mnm;
};

template<typename Type>
Type vect<Type>::max() const {
    assert(v.size() > 0);
    Type mxm = v[0];
    for (int i = 0; i < v.size(); ++i)
        if (mxm < v[i])
            mxm = v[i];
    return mxm;
}

template<typename Type>
int vect<Type>::size() const {
    return v.size();
}
/// product of a vector by a number
/// \param mult
/// \return

template<>
vect<db> vect<db>::operator * (db mult) const {
    vect<db> prod = vect<db>(size());
    for (int i = 0; i < size(); i++)
        prod[i] = v[i] * mult;
    return prod;
}

template<>
db vect<db>::scalarProd(vect<db> a, vect<db> b) {
    db res = 0;
    for (int i = 0; i < a.size(); i++) {
        res += a[i] * b[i];
    }
    return res;
}

template<typename Type>
db vect<Type>::getNorm(vect<db> a) {
    assert(a.size() > 0);
    db mxm = abs(a[0]);
    for (int i = 1; i < a.size(); i++)
        mxm = abs(a[i]) > mxm ? abs(a[i]): mxm;
    return mxm;
}

template<> vect<db> vect<db>::normalize(vect<db> x) {
    //найти норму вектора
    db vNorm = getNorm(x);
    //разделить вектор на его норму
    for (int i = 0; i < x.size(); i++)
        x[i] = x[i] / vNorm;
    return x;
}

// норма разности двух векторов
template<>
db vect<db>::diffNorm(vector<db> a, vector<db> b) {
    db mxm = abs(a[0] - b[0]);
    for (int i = 1; i < a.size(); i++)
        if (mxm < abs(b[i] - a[i]))
            mxm = abs(b[i] - a[i]);
    return mxm;
}


template<typename Type>
std::ostream &operator<<(std::ostream &out, const vect<Type> &v) {
    for (int i = 0; i < v.size(); ++i)
        out << v[i] << " ";
    return out;
}
#endif //NUMERICAL_METHODS_VECT_H
