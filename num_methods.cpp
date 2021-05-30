//
// Created by kesa on 29.05.2021.
//
#include <valarray>
#include "num_methods.h"

std::vector<db> autoGen(int n, db l, db r) {
    std::vector<db> x = std::vector<db>(n);
    for (int i = 0; i < n; i++)
        x[i] = l + i * (r - l) / (n - 1);
    return x;
}

vector<db> AbstractIntegral::divideIntervalsInto2( vector<db> &x) {
    vector<db> x2 = vector < db > ();
    x2.push_back(x[0]);
    for (int i = 1; i < x.size(); i++) {
        x2.push_back((x[i] + x[i - 1]) / 2);
        x2.push_back(x[i]);
    }
    return x2;
}

vector<db> AbstractIntegral::calcPoints(db l, db r, db eps, db A, db B) {
    //исходное разбиение отрезка на точки
    vector < db > x = {l, r};
    db h = x[1] - x[0];
    //разделить каждый интервал на два,
    //уменьшив шаг вдвое
    vector < db > x2 = divideIntervalsInto2(x);
    db intH = calcIntegral(x), intHalfH = calcIntegral(x2);
    db diff = abs(intH - intHalfH) / (pow(2, degreeOnInterval - 1) - 1);
    db mxmDiff = eps * h / (B - A);
    if (diff > mxmDiff) {
        vector < db > lPoints = calcPoints(l, (l + r) / 2., eps, A, B);
        vector < db > rPoints = calcPoints((l + r) / 2., r, eps, A, B);
        lPoints.insert(lPoints.end(), rPoints.begin() + 1, rPoints.end());
        x = lPoints;
    }
    return x;
}

///составить вектор точек, расстояние между которыми задано пользователем
vector<db> AbstractIntegral::getPointsWithFixedStep(db A, db B, int n) {
    vector < db > x = vector < db > (n);
    db h = (B - A) / db(n - 1);
    for (int i = 0; i < n - 1; i++)
        x[i] = A + h * i;
    x[n - 1] = B;
    return x;
}


db LeftRectangleIntegral::calcIntegral(vector<db> x) {
    db res = 0;
    for (int i = 1; i < x.size(); i++)
        res += subFunc(x[i - 1]) * (x[i] - x[i - 1]);
    return res;
}

db RightRectangleIntegral::calcIntegral(vector<db> x) {
    db res = 0;
    for (int i = 1; i < x.size(); i++)
        res += subFunc(x[i]) * (x[i] - x[i - 1]);
    return res;
}

db MiddleRectangleIntegral::calcIntegral(vector<db> x) {
    db res = 0;
    for (int i = 1; i < x.size(); i++) {
        db h = x[i] - x[i - 1];
        res += subFunc(x[i - 1] + h / 2.) * h;
    }
    return res;
}

db TrapeziumIntegral::calcIntegral(vector<db> x) {
    db res = 0, yiPrev = subFunc(x[0]);
    for (int i = 1; i < x.size(); i++) {
        db h = x[i] - x[i - 1];
        db yi = subFunc(x[i]);
        res += (yi + yiPrev) / 2 * h;
        yiPrev = yi;
    }
    return res;
}

db SimpsonsIntegral::calcIntegral(vector<db> x) {
    //число точек должно быть четным
    assert(x.size() % 2 == 1);
    db res = 0, yPrev = subFunc(x[0]), yi = subFunc(x[1]);
    for (int i = 1; i < x.size() - 1; i++) {
        db h = x[i + 1] - x[i - 1];
        db yNext = subFunc(x[i + 1]);
        res += (yPrev + 4 * yi + yNext) * h;
        yPrev = yi;
        yi = yNext;
    }
    res /= 6.;
    return res;
}

vector<db> ThomasAlgorithm(vector<db> a, vector<db> b, vector<db> c, vector<db> d) {
    int n = a.size();
    vector < db > alpha = vector < db > (n + 1);
    vector < db > beta = vector < db > (n + 1);
    vector < db > x = vector < db > (n);
    a[0] = c[c.size() - 1] = 0;
    //прямой ход прогонки: найти вспомогательные коэффициенты альфа и бета
    alpha[1] = c[0] / b[0];
    beta[1] = -1.0 * d[0] / b[0];
    for (int i = 1; i < n; i++) {
        alpha[i + 1] = c[i] / (b[i] - a[i] * alpha[i]);
        beta[i + 1] = (a[i] * beta[i] - d[i]) / (b[i] - a[i] * alpha[i]);
    }
    //обратный ход прогонки - найти искомые переменные
    x[n - 1] = beta[n];
    for (int i = n - 2; i > 0; i--) {
        x[i] = alpha[i + 1] * x[i + 1] + beta[i + 1];
    }
    x[0] = c[0] / b[0] * x[1] - d[1] / b[1];
    return x;
}

/// <summary>
/// Матрица коэффициентов многочленов интерполирования сплайнами: a, b, c и d
/// в 0, 1, 2 и 4 строке соответственно.
/// i-тый столбец соответствует многолену на i-том отрезке
/// </summary>
/// <param name="x">точки интерполяции</param>
/// <param name="y">значение функции в них</param>
/// <returns> матрица коэффициентов многочленов</returns>
vector<vector<db>> Spline::getSplineMatrix(vector<db> x, vector<db> y) {
    assert(x.size() == y.size());
    int n = x.size();
    vector < vector < db >> splMatr = vector < vector < db >> (4);
    //в точках интерполяции коэффициент а равен
    //значению исходной функции
    vector < db > a = y;
    // у векторов b и d нулевого элемента нет
    //вектора, соответствующие коэффициентам в сплайне
    vector < db > b = vector < db > (n);
    vector < db > c = vector < db > (n);
    vector < db > d = vector < db > (n);
    //задаем трехдиагональную матрицу 4 векторами
    //вектора, соответствующие коэффициентам в трехдиагональной матрице (задающей
    //коэффициент 'с' сплайна)
    vector < db > at = vector < db > (n, 0);
    vector < db > bt = vector < db > (n, 0);
    vector < db > ct = vector < db > (n, 0);
    //вектор свободных членов
    vector < db > dt = vector < db > (n, 0);
    //заполнить трехдиагональную матрицу и свободный член коэффициентов при c
    for (int i = 1; i < n - 1; i++) {
        db h1 = x[i] - x[i - 1];
        db h2 = x[i + 1] - x[i];
        at[i] = h1;
        bt[i] = 2 * (h1 + h2) * -1.0;
        ct[i] = h2;
        dt[i] = 6 * ((y[i + 1] - y[i]) / h2 - (y[i] - y[i - 1]) / h1);
    }
    //вычислить значение коэффициентов методом прогонки
    bt[0] = bt[n - 1] = 1;
    c = ThomasAlgorithm(at, bt, ct, dt);
    //используя формулы, вычислить остальные коэффициенты
    for (int i = 1; i < n; i++) {
        db h1 = x[i] - x[i - 1];
        d[i] = (c[i] - c[i - 1]) / h1;
        b[i] = (h1 / 2) * c[i] - ((h1 * h1) / 6) * d[i] + (y[i] - y[i - 1]) / h1;
    }
    splMatr[0] = a;
    splMatr[1] = b;
    splMatr[2] = c;
    splMatr[3] = d;
    return splMatr;
}

/// <summary>
/// Результат приближения в заданной точке
/// </summary>
/// <param name="curx">точка, в которой вычисляется приближение</param>
/// <param name="x">точки интерполяции</param>
/// <returns>результат приближения</returns>

db Spline::getFuncApproxInDot(db curx)  {
    db Sx;
    int i;
    //определить, какому отрезку принадлежит точка
    assert(curx > approxDot[0] && curx < approxDot[approxDot.size() - 1]);
    for (i = 1; i < approxDot.size(); i++)
        if (curx <= approxDot[i])
            break;
    //вычислить значение многочлена
    Sx = SplineMatrix[0][i] + SplineMatrix[1][i] * (curx - approxDot[i]) +
         (SplineMatrix[2][i] / 2) * pow((curx - approxDot[i]), 2) +
         (SplineMatrix[3][i] / 6) * pow((curx - approxDot[i]), 3);
    return Sx;
}

ParametricallyDefinedArea::ParametricallyDefinedArea(vector<db> x, vector<db>y) :x(x), y(y) {
    int xs = x.size();
    assert(xs == y.size());
    t = vect<db>(xs);
    for (int i = 0; i < xs; i++)
        t[i] = i + 1.;
    xSpline = Spline(t, x);
    ySpline = Spline(t, y);
}

std::ostream &operator<<(std::ostream &out,  ParametricallyDefinedArea &a) {
    //исходные точки
    out << a.x << endl << a.y;
    //точки для построения графика
    vector<db> printDots = autoGen(EPS * a.t[0], a.t[0], a.t[a.t.size() - 1]);
    out << a.xSpline.getFuncApproxInDots(printDots) << endl;
    out << a.ySpline.getFuncApproxInDots(printDots) << endl;
    return out;
}
vector<Point<db>> createExternalRectangleArea(ParametricallyDefinedArea fig) {
    vector<Point<db>> coord = vector<Point<db>>();
    db delta = 2;
    db minx, miny, maxx, maxy;
    minx = fig.x.min() - delta;
    miny = fig.y.min() - delta;
    maxx = fig.x.max() + delta;
    maxy = fig.y.max() + delta;
    coord.push_back(Point<db>(minx, maxy));
    coord.push_back(Point<db>(minx, miny));
    coord.push_back(Point<db>(maxx, miny));
    coord.push_back(Point<db>(maxx, maxy));
    coord.push_back(coord[0]);
    return coord;
}

