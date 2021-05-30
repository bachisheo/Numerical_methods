//
// Created on 29.05.2021.
//
#ifndef NUMERICAL_METHODS_VECT_H
#define NUMERICAL_METHODS_VECT_H

#include <cassert>
#include <vector>
#include <iostream>

#define EPS 1e-6

typedef double db;
using namespace std;

template<typename T>
class Point {
public:
    T x = 0, y = 0;

    Point(T x, T y) : x(x), y(y) {}
};

template<typename T>
class vect : public std::vector<T> {
    std::vector<T> v;
public:
    vect<T>(int size) {
        v = vector<T>(size);
    }

    vect<T>() : vector<T>(0) {}

    vect<T>(vector<T> oldVector) : v(oldVector) {
        oldVector;
    }

    T min();

    T max();

    int size() const;

    vect<T> operator*(T mult) {
        vect<T> prod = vect<T>(size());
        for (int i = 0; i < size(); i++)
            prod[i] = v[i] * mult;
        return prod;
    };

    T &operator[](int index) {
        assert(index >= 0 && index < v.size());
        return v[index];
    }

    const T &operator[](int index) const {
        assert(index >= 0 && index < v.size());
        return v[index];
    }


    //скалярное произведение векторов
    static T scalarProd(vect<T> a, vect<T> b);

    //найти норму вектора
    static db getNorm(vect<db> a);

    //нормирование вектора
    static vect<db> normalize(vect<db> x);

    // норма разности двух векторов
    static db diffNorm(vector<db> a, vector<db> b);
};

template<typename T>
std::ostream &operator<<(ostream &fout, vector<Point<T>> points) {
    for (int i = 0; i < points.size(); i++)
        fout << points[i].x << " ";
    fout << endl;
    for (int i = 0; i < points.size(); i++)
        fout << points[i].y << " ";
    fout << endl;
    return fout;
}

template<typename T>
T vect<T>::min() {
    assert(v.size() > 0);
    T mnm = v[0];
    for (int i = 0; i < v.size(); ++i)
        if (mnm > v[i])
            mnm = v[i];
    return mnm;
};

template<typename T>
T vect<T>::max() {
    assert(v.size() > 0);
    T mxm = v[0];
    for (int i = 0; i < v.size(); ++i)
        if (mxm < v[i])
            mxm = v[i];
    return mxm;
}

template<typename T>
int vect<T>::size() const {
    return v.size();
}

template<typename T>
T vect<T>::scalarProd(vect<T> a, vect<T> b) {
    T res = 0;
    for (int i = 0; i < a.size(); i++) {
        res += a[i] * b[i];
    }
    return res;
}

template<typename T>
db vect<T>::getNorm(vect<db> a) {
    assert(a.size() > 0);
    db mxm = abs(a[0]);
    for (int i = 1; i < a.size(); i++)
        mxm = abs(a[i]) > mxm ? abs(a[i]) : mxm;
    return mxm;
}

template<typename T>
vect<db> vect<T>::normalize(vect<db> x) {
    //найти норму вектора
    db vNorm = getNorm(x);
    //разделить вектор на его норму
    for (int i = 0; i < x.size(); i++)
        x[i] = x[i] / vNorm;
    return x;
}

// норма разности двух векторов
template<typename T>
db vect<T>::diffNorm(vector<db> a, vector<db> b) {
    db mxm = abs(a[0] - b[0]);
    for (int i = 1; i < a.size(); i++)
        if (mxm < abs(b[i] - a[i]))
            mxm = abs(b[i] - a[i]);
    return mxm;
}


template<typename T>
std::ostream &operator<<(std::ostream &out, vect<T> v) {
    for (int i = 0; i < v.size(); ++i)
        out << v[i] << " ";
    return out;
}

#endif //NUMERICAL_METHODS_VECT_H
