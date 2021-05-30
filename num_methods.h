//
// Created by kesa on 29.05.2021.
//

#ifndef NUMERICAL_METHODS_NUM_METHODS_H
#define NUMERICAL_METHODS_NUM_METHODS_H
#include "vect.h"

//Абстрактный класс интеграла, содержит основные
//поля и методы
class AbstractIntegral {
private:
    ///разделить каждый интервал на два
    vector<db> divideIntervalsInto2(const vector<db> &x);

    /// <summary>
    ///автоматический выбор шага интегрирования
    /// </summary>
    /// <param name="l">текущая левая граница</param>
    /// <param name="r"></param>
    /// <param name="eps">точность</param>
    /// <param name="A">левая граница интегрирования</param>
    /// <param name="B"></param>
    /// <returns></returns>
    vector<db> calcPoints(db l, db r, db eps, db A, db B);

protected:
    //интегрируемая функция
    db (*subFunc)(db x);

//порядок точности метода
    db degreeOnInterval;
public:
    AbstractIntegral(db(*func)(db x), int degree) : subFunc(func), degreeOnInterval(degree) {};

/// вычислить значение интеграла, используя значения функции
///в точках передаваемого вектора
/// \param x
/// \return
    virtual db calcIntegral(vector<db> x) = 0;

///  автоматический выбор шага интегрирования
/// \param A левая граница
/// \param B правая граница
/// \param eps точность вычислений
/// \return вектор точек интегрирования
    vector<db> calcPoints(db A, db B, db eps) {
        return calcPoints(A, B, eps, A, B);
    }

    ///составить вектор точек, расстояние между которыми задано пользователем
    /// \param A левая граница
    /// \param B правая граница
    /// \param n количество точек
    /// \return вектор точек интегрирования
    vector<db> getPointsWithFixedStep(db A, db B, int n);
};

//интеграл, вычисляемый методом левых прямоугольников
class LeftRectangleIntegral : public AbstractIntegral {
public:
    LeftRectangleIntegral(db(*FUNC)(db x)) : AbstractIntegral(FUNC, 2) {};

    db calcIntegral(vector<db> x);
};

//интеграл, вычисляемый методом правых прямоугольников
class RightRectangleIntegral : public AbstractIntegral {
public:
    RightRectangleIntegral(db(*FUNC)(db x)) : AbstractIntegral(FUNC, 2) {};

    db calcIntegral(vector<db> x);
};

//интеграл, вычисляемый методом средних прямоугольников
class MiddleRectangleIntegral : public AbstractIntegral {
public:
    MiddleRectangleIntegral(db(*FUNC)(db x)) : AbstractIntegral(FUNC, 3) {}

    db calcIntegral(vector<db> x);
};

//интеграл, вычисляемый методом трапеций
class TrapeziumIntegral : public AbstractIntegral {
public:
    TrapeziumIntegral(db(*FUNC)(db x)) : AbstractIntegral(FUNC, 3) {}

    db calcIntegral(vector<db> x);

};

//интеграл, вычисляемый методом Симпсона
class SimpsonsIntegral : public AbstractIntegral {
public:

    db calcIntegral(vector<db> x);

    SimpsonsIntegral(db(*FUNC)(db x)) : AbstractIntegral(FUNC, 5) {}
};

/// Решение СЛАУ  с трехдиагональной матрицей
///методом прогонки (методом Томаса)
/// \param a коэффициенты элементов на диагонали, ниже главной
/// \param b коэффициенты элементов на главной диагонали
/// \param c коэффициенты элементов на диагонали, выше главной
/// \param d столбец свободных членов
/// \return решение СЛАУ
vector<db> ThomasAlgorithm(vector<db> a, vector<db> b, vector<db> c, vector<db> d);

struct Spline {
private:
    //matrix of spline's coefficient
    vector<vector<db>> SplineMatrix;
    vector<db> approxDot;
    vector<db> funcInDot;

    vector<vector<db>> getSplineMatrix(vector<db> x, vector<db> y);

public:
    Spline(){}
    Spline(vector<db> x, vector<db> y) : approxDot(x), funcInDot(y),
                                         SplineMatrix(getSplineMatrix(x, y)) {}
    //approximate func for one dot
    db getFuncApproxInDot(db curx) const;

    //approximate for vector of dots
    vect<db> getFuncApproxInDots(vector<db> curx) const{
        vect<db> res = vector<db>(curx.size());
        for (int i = 0; i < curx.size(); ++i)
            res[i] = getFuncApproxInDot(curx[i]);
        return res;
    }
};

std::vector<db> autoGen(int n, db l, db r);

//область, заданная параметрически
struct ParametricallyDefinedArea {
    vect<db> x;
    vect<db> y;
    vect<db> t;
    Spline xSpline;
    Spline ySpline;
    ParametricallyDefinedArea(vector<db> x, vector<db> y);
};
std::ostream &operator<<(std::ostream &out, const ParametricallyDefinedArea &a) ;

vector<Point<db>> createExternalRectangleArea(ParametricallyDefinedArea fig) ;
#endif //NUMERICAL_METHODS_NUM_METHODS_H
