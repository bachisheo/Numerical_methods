//
// Created by kesa on 29.05.2021.
//
#ifndef NUMERICAL_METHODS_VECT_H
#define NUMERICAL_METHODS_VECT_H
#include <cassert>
#include <vector>
#include <iostream>

typedef double db;

using namespace std;

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
    static vector<db> normalize(vector<db> x);

    // норма разности двух векторов
    static db diffNorm(vector<db> a, vector<db> b);
};

template<typename Type>
std::ostream &operator<<(std::ostream &out, const vect<Type> &v);


#endif //NUMERICAL_METHODS_VECT_H
