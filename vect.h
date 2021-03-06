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
class vect  {
private :
    vector<T> v;
public:
    vect<T>();
    vect<T>(int size);
    vect<T>(vector<T> oldVector);

    T min();
    T max();
    int size() const;

    vect<T> operator*(T mult);
    T &operator [] (int index);


    operator vector<T>() const { return v; }

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
vect<T>::vect() {
    v = vector<T>();
}

template<typename T>
vect<T>::vect(int size) {
    v = vector<T>(size);
}

template<typename T>
vect<T>::vect(vector<T> oldVector) {
    v = oldVector;
    int a = 1;
    a++;
}
template<typename T>
T vect<T>::min() {
    assert(size() > 0);
    T mnm = this->operator[](0);
    for (int i = 0; i < size(); ++i)
        if (mnm > this->operator[](i) )
            mnm = this->operator[](i) ;
    return mnm;
};

template<typename T>
T vect<T>::max() {
    assert(size() > 0);
    T mxm = this->operator[](0) ;
    for (int i = 0; i < size(); ++i)
        if (mxm < this->operator[](i) )
            mxm = this->operator[](i) ;
    return mxm;
}

template<typename T>
int vect<T>::size() const {return  v.size();}

template<typename T>
vect<T> vect<T>::operator*(T mult) {
    vect<T> prod = vect<T>(size());
    for (int i = 0; i < size(); i++)
        prod[i] = v[i] * mult;
    return prod;
}

template<typename T>
T& vect<T>::operator[](int index) {
    return v[index];
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

