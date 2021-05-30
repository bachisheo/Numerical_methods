//
// Created by kesa on 29.05.2021.
//
#include <valarray>
#include "num_methods.h"
template <typename T>
std::vector <db> autoGen(int n, db l, db r) {
    std::vector<db> x = std::vector<db>(n);
    for (int i = 0; i < n; i++)
        x[i] = l + i * (r - l) / (n - 1);
    return x;
}

/// <summary>
/// Решение трехдиагональных СЛАУ методом прогонки
/// </summary>
/// <param name="a">коэффициенты элементов на диагонали, ниже главной</param>
/// <param name="b">коэффициенты элементов на главной диагонали</param>
/// <param name="c">коэффициенты элементов на диагонали, выше главной</param>
/// <param name="d">столбец свободных членов</param>
/// <returns>решение СЛАУ</returns>

vector<db> ThomasAlgorithm(vector<db> a, vector<db> b, vector<db> c, vector<db> d) {
    int n = a.size();
    vector<db> alpha = vector<db>(n + 1);
    vector<db> beta = vector<db>(n + 1);
    vector<db> x = vector<db>(n);
    a[0] = c[c.size() - 1] = 0;
    //прямой ход прогонки: найти вспомогательные коэффициенты альфа и бета
    alpha[1] = c[0] / b[0];
    beta[1] = -1.0 * d[0] / b[0];
    for (int i = 1; i < n; i++)
    {
        alpha[i + 1] = c[i] / (b[i] - a[i] * alpha[i]);
        beta[i + 1] = (a[i] * beta[i] - d[i]) / (b[i] - a[i] * alpha[i]);
    }
    //обратный ход прогонки - найти искомые переменные
    x[n - 1] = beta[n];
    for (int i = n - 2; i > 0; i--)
    {
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
    vector<vector<db>> splMatr = vector<vector<db>>(4);
    //в точках интерполяции коэффициент а равен
    //значению исходной функции
    vector<db> a = y;
    // у векторов b и d нулевого элемента нет
    //вектора, соответствующие коэффициентам в сплайне
    vector<db> b = vector<db>(n);
    vector<db> c = vector<db>(n);
    vector<db> d = vector<db>(n);
    //задаем трехдиагональную матрицу 4 векторами
    //вектора, соответствующие коэффициентам в трехдиагональной матрице (задающей
    //коэффициент 'с' сплайна)
    vector<db> at = vector<db>(n, 0);
    vector<db> bt = vector<db>(n, 0);
    vector<db> ct = vector<db>(n, 0);
    //вектор свободных членов
    vector<db> dt = vector<db>(n, 0);
    //заполнить трехдиагональную матрицу и свободный член коэффициентов при c
    for (int i = 1; i < n - 1; i++)
    {
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
    for (int i = 1; i < n; i++)
    {
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

db Spline::getFuncApproxInDot(db curx) {
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

