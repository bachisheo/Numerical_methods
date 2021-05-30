//
// Created by kesa on 29.05.2021.
//
#include "vect.h"
#ifndef NUMERICAL_METHODS_NUM_METHODS_H
#define NUMERICAL_METHODS_NUM_METHODS_H

vector<db> ThomasAlgorithm(vector<db> a, vector<db> b, vector<db> c, vector<db> d);

struct Spline {
private:
    //matrix of spline's coefficient
    vector<vector<db>> SplineMatrix;
    vector<db> approxDot;
    vector<db> funcInDot;
    vector<vector<db>> getSplineMatrix(vector<db> x, vector<db> y);
public:
    Spline(vector<db> x, vector<db> y):approxDot(x),funcInDot(y),
    SplineMatrix(getSplineMatrix(x, y)) {}
    //approximate func for one dot
    db getFuncApproxInDot(db curx);
    //approximate for vector of dots
    vector<db> getFuncApproxInDots(vector<db> curx){
        vector<db> res = vector<db>(curx.size());
        for (int i = 0; i < curx.size(); ++i)
            res[i] = getFuncApproxInDot(curx[i]);
        return res;
    }
};
std::vector <db> autoGen(int n, db l, db r);



#endif //NUMERICAL_METHODS_NUM_METHODS_H
