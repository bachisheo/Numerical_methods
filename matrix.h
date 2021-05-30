//
// Created on 29.05.2021.
//

#ifndef NUMERICAL_METHODS_MATRIX_H
#define NUMERICAL_METHODS_MATRIX_H
#include "vect.h"

 long long C = 10000;
using namespace std;
 double eps = 1e-7;
 int maxIterationNumber = 50000;

template<typename T>
class matrix {
protected:
public:
    vect<vect<T>> dat;

    int colNumb()  { return dat.size(); }

    int rowNumb()  { return dat.size() > 0 ? dat[0].size() : 0; }

    matrix() {}

    /// create an empty matrix of given dimensions
    /// \param cstr
    /// \param ccol
    matrix(int cstr, int ccol) {
        dat = vect<vect<T>>(cstr);
        for (int i = 0; i < ccol; i++)
            dat[i] = vect<T>(ccol);
    }

    /// конструктор матрицы из двумерного вектора
    /// \param M
    matrix(vect<vect<T>> M) : dat(M) {};

    /// Конструктор одномерной матрицы из вектора
    /// \param M
    matrix(vect<T> v) {
        dat = matrix<db>(v.size(), 1).dat;
        for (int i = 0; i < v.size(); i++)
            dat[i][0] = v[i];
    }

    /// добавить одномерный вектор
    /// в качестве еще одного столбца в матрице
    /// \param vec
    void push_back(vect<T> vec);

    /// основные алгебраические операции над матрицами
    /// \param otherMatrix
    /// \return
    // после параметров - предупреждает о том, что метод не изменяет
    //поля класса
    matrix<T> operator+( matrix<T> &otherMatrix) ;

    matrix<T> operator+(T d) ;

    matrix<T> operator*( matrix<T> &b) ;

    matrix<T> operator*(T d) ;

    vector<T> operator*(vector<T> d) ;

    vect<T> &operator[]( int rowIndex) ;

    /// получить матрицу, транспонированную к исходной
    /// \return
    matrix getTransMatrix() ;
};

template<typename T>
istream &operator>>(istream &in, matrix<T> &m);

class TriangleMatrix : public matrix<db> {
    //счетчик перестановок строк
    int numbOfTurns = 0;
    //определитель
    db det = 0;
    //измененный порядок строк
    vector<int> orderOfString;

    //поменять строки треугольной матрицы местами
    void swapString(int a, int b);

public:
    TriangleMatrix(vect<vect<db>> vec);

    //вычислить определитель матрицы
    double calcDet();

    //преобразовать вектор b к треугольному виду матрицы
    vect<db> transformVector(vect<db> b);

};

template<typename T>
std::ostream &operator<<(std::ostream &out, const matrix<T> &v) {
    for (int i = 0; i < v.rowNumb(); ++i)
        out << v[i] << endl;
    return out;
}

template<typename T>
std::istream &operator>>(std::istream &in,  vect<T> &v) {
    for (int i = 0; i < v.size(); ++i)
        in >> v[i];
    return in;
}


//обратный ход
vector<double> revStep(vector<double> b, vector<vector<double>> a);

//поиск решения СЛАУ методом Гаусса
vector<double> GaussMethod(matrix<db> B, vect<db> b);

/// вычислить невязку
/// \param A матрица коэффициентов
/// \param e1  точность
/// \param l1
/// \return
template<typename T>
vector<db> residual(matrix<db> A, vector<db> e1, db l1);
/// сумма двух матриц
/// \param otherMatrix
/// \return
template<typename T>
matrix<T> matrix<T>::operator+( matrix<T> &otherMatrix)  {
    int n = rowNumb(), m = rowNumb();
    assert(m == otherMatrix.colNumb() && n == otherMatrix.rowNumb());
    matrix c = matrix(otherMatrix);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            c.dat[i][j] = otherMatrix.dat[i][j] + dat[i][j];
    return c;
}


template<typename T>
matrix<T> matrix<T>::operator+(T d) {
    int n = rowNumb(), m = colNumb();
    matrix c = matrix(this->dat);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            c.dat[i][j] += d;
    return c;
}

template<typename T>
matrix<T> matrix<T>::operator*( matrix<T> &b)  {
    int n = rowNumb(), m = colNumb();
    matrix<T> c = matrix<T>(0, 0);
    if (m == b.rowNumb()) {
        c = matrix(n, b.colNumb());
        for (int i = 0; i < n; i++)
            for (int j = 0; j < b.colNumb(); j++) {
                c.dat[i][j] = 0;
                for (int k = 0; k < m; k++)
                    c.dat[i][j] += dat[i][k] * b.dat[k][j];
            }
    }
    return c;
}

template<typename T>
matrix<T> matrix<T>::operator*(T d)  {
    matrix c = matrix(rowNumb(), colNumb());
    for (int i = 0; i < rowNumb(); i++)
        for (int j = 0; j < colNumb(); j++)
            c.dat[i][j] = dat[i][j] * d;
    return c;
}

///умножение матрицы на вектор
/// \param v
/// \return
template<typename T>
vector<T> matrix<T>::operator*(vector<T> v)  {
    matrix<T> b = matrix<T>(v);
    return ((*this) * b).getTransMatrix().dat[0];
}

template<typename T>
vect<T> &matrix<T>::operator[]( int rowIndex)  {
    return dat[rowIndex];
}

template<typename T>
matrix<T> matrix<T>::getTransMatrix()  {
    matrix B = matrix(colNumb(), rowNumb());
    for (int i = 0; i < rowNumb(); i++)
        for (int j = 0; j < colNumb(); j++)
            B.dat[j][i] = dat[i][j];
    return B;
}

template<typename T>
void matrix<T>::push_back(vect<T> vec) {
    dat.push_back(vec);
}
#endif //NUMERICAL_METHODS_MATRIX_H
