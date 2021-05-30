//
// Created on 29.05.2021.
//

#ifndef NUMERICAL_METHODS_MATRIX_H
#define NUMERICAL_METHODS_MATRIX_H
#include "vect.h"

const long long C = 10000;
using namespace std;
const double eps = 1e-7;
const int maxIterationNumber = 50000;

template<typename Type>
class matrix {
protected:
public:
    vect<vect<Type>> dat;

    int colNumb() const { return dat.size(); }

    int rowNumb() const { return dat.size() > 0 ? dat[0].size() : 0; }

    matrix() {}

    /// create an empty matrix of given dimensions
    /// \param cstr
    /// \param ccol
    matrix(int cstr, int ccol) {
        dat = vect<vect<Type>>(cstr);
        for (int i = 0; i < ccol; i++)
            dat[i] = vect<Type>(ccol);
    }

    /// конструктор матрицы из двумерного вектора
    /// \param M
    matrix(vect<vect<Type>> M) : dat(M) {};

    /// Конструктор одномерной матрицы из вектора
    /// \param M
    matrix(vect<Type> v) {
        dat = matrix<db>(v.size(), 1).dat;
        for (int i = 0; i < v.size(); i++)
            dat[i][0] = v[i];
    }

    /// добавить одномерный вектор
    /// в качестве еще одного столбца в матрице
    /// \param vec
    void push_back(vect<Type> vec);

    /// основные алгебраические операции над матрицами
    /// \param otherMatrix
    /// \return
    //const после параметров - предупреждает о том, что метод не изменяет
    //поля класса
    matrix operator+(const matrix &otherMatrix) const;

    matrix operator+(Type d) const;

    matrix operator*(const matrix &b) const;

    matrix operator*(Type d) const;

    vector<Type> operator*(vector<Type> d) const;

    vect<Type> &operator[](const int rowIndex);

    /// получить матрицу, транспонированную к исходной
    /// \return
    matrix getTransMatrix() const;
};

template<typename Type>
istream &operator>>(istream &in, matrix<Type> &m);

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

template<typename Type>
std::ostream &operator<<(std::ostream &out, const matrix<Type> &v) {
    for (int i = 0; i < v.rowNumb(); ++i)
        out << v[i] << endl;
    return out;
}

template<typename Type>
std::istream &operator>>(std::istream &in, const vect<Type> &v) {
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
template<typename Type>
vector<db> residual(matrix<db> A, vector<db> e1, db l1);
/// сумма двух матриц
/// \param otherMatrix
/// \return
ttt
matrix<Type> matrix<Type>::operator+(const matrix<Type> &otherMatrix) const {
    int n = rowNumb(), m = rowNumb();
    assert(m == otherMatrix.colNumb() && n == otherMatrix.rowNumb());
    matrix c = matrix(otherMatrix);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            c.dat[i][j] = otherMatrix.dat[i][j] + dat[i][j];
    return c;
}


template<>
matrix<db> matrix<db>::operator+(db d) const {
    int n = rowNumb(), m = colNumb();
    matrix c = matrix(this->dat);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            c.dat[i][j] += d;
    return c;
}

template<>
matrix<db> matrix<db>::operator*(const matrix<db> &b) const {
    int n = rowNumb(), m = colNumb();
    matrix<db> c = matrix(0, 0);
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

template<>
matrix<db> matrix<db>::operator*(db d) const {
    matrix c = matrix(rowNumb(), colNumb());
    for (int i = 0; i < rowNumb(); i++)
        for (int j = 0; j < colNumb(); j++)
            c.dat[i][j] = dat[i][j] * d;
    return c;
}

///умножение матрицы на вектор
/// \param v
/// \return
template<>
vector<db> matrix<db>::operator*(vector<db> v) const {
    matrix<db> b = matrix<db>(v);
    return ((*this) * b).getTransMatrix().dat[0];
}

template<>
vect<db> &matrix<db>::operator[](const int rowIndex) {
    return dat[rowIndex];
}

ttt
matrix<Type> matrix<Type>::getTransMatrix() const {
    matrix B = matrix(colNumb(), rowNumb());
    for (int i = 0; i < rowNumb(); i++)
        for (int j = 0; j < colNumb(); j++)
            B.dat[j][i] = dat[i][j];
    return B;
}

ttt
void matrix<Type>::push_back(vect<Type> vec) {
    dat.push_back(vec);
}



#endif //NUMERICAL_METHODS_MATRIX_H
