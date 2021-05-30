//
// Created by kesa on 29.05.2021.
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

    matrix(){}

    /// create an empty matrix of given dimensions
    /// \param rowNumber
    /// \param volNumber
    matrix(int rowNumber, int volNumber);

    /// конструктор матрицы из двумерного вектора
    /// \param M
    matrix(vect<vect<Type>> M):dat(M){};

    /// Конструктор одномерной матрицы из вектора
    /// \param M
    matrix(vect<Type> v);

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

    ///
    ///умножение матрицы на вектор
    ///
    const vect<Type> operator*(vect<Type> v) const;

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
    for (int i = 0; i < v.rowNumb(); ++i)
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


#endif //NUMERICAL_METHODS_MATRIX_H
