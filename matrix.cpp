//
// Created by kesa on 29.05.2021.
//

#include "matrix.h"
#ifndef SIMPLEXMETHOD_MYMATRIX_H
#define SIMPLEXMETHOD_MYMATRIX_H
#define ttt template <typename Type>




///
/// \param cstr - количество строк
/// \param ccol количество столбцов
ttt
matrix<Type>::matrix(int cstr, int ccol) {
    dat = vect<vect<Type>>(cstr);
    for (int i = 0; i < ccol; i++)
        dat[i] = vect<Type>(ccol);
};



/// Конструктор одномерной матрицы (столбца) из вектора
/// \param M
ttt
matrix<Type>::matrix(vect<Type> v){
    dat = matrix(v.size(), 1);
    for (int i = 0; i < v.size(); i++)
        dat[i][0] = v[i];
}

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


template<> matrix<db> matrix<db>::operator+(db d) const {
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
const vector<db> matrix<db>::operator * (vector<db> v) const {
    matrix<db> b = matrix(v);
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

void TriangleMatrix::swapString(int a, int b) {
    //зафиксировать изменение в порядке строк
    swap(orderOfString[a], orderOfString[b]);
    //поменять местами содержимое строк
    vect<db> tmp = dat[a];
    dat[a] = dat[b];
    dat[b] = tmp;
    numbOfTurns++;
}

TriangleMatrix::TriangleMatrix(vect<vect<db>> M)  {
    vect<vect<db>> dat = M;
    int n = rowNumb(), m = colNumb();
    //привести матрицу к треугольному виду
    orderOfString = vector<int>(n);
    for (int i = 0; i < n; i++)
        orderOfString[i] = i;

    for (int k = 0; k < n; k++) {
        //найти главный (ведущий элемент) в столбце
        int indOfMxm = k;
        for (int p = k; p < n; p++)
            if (abs(dat[p][k]) > abs(dat[indOfMxm][k]))
                indOfMxm = p;

        //если все эелемнты 0 - пропустить итерацию
        if (abs(this->dat[indOfMxm][k]) < eps) {
            det = 0;
            continue;
        }

        //поместить ведущий элемент на главную диагональ
        if (indOfMxm != k) swapString(indOfMxm, k);
        //в каждой строке, начиная с k+1
        for (int p = k + 1; p < n; p++) {
            //найти коэффициент  С[p][k]
            double C = dat[p][k] / dat[k][k];
            //вычесть из текущей строки, k-тую строку
            //умноженную на С
            for (int l = k; l < m; l++)
                dat[p][l] = dat[p][l] - dat[k][l] * C;
            //сохранить коэффициент
            dat[p][k] = C;
        }
    }
}

//вычислить определитель матрицы
double TriangleMatrix<db>::getDet() {
    det = 1;
    //вычислить определитель треугольной матрицы
    //как произведение элементов, стоящих на главной диагонали
    for (int k = 0; k < size; k++) {
        if (abs(dat[k][k]) < eps) return 0;
        det *= dat[k][k];
    }
    //изменить знак определителя, если строки
    //были переставлены четное число раз
    if (numbOfTurns % 2) det *= (-1);
    return det;
}

//преобразовать вектор b к трекгольному виду матрицы
vect<db> TriangleMatrix::transformVector(vect<db> b) {
    vector<db> bk;
    assert(b.size() == size);
    bk = vector<double>(size);
    //переставить строки в том же порядке, что и у треугольной матрицы А
    for (int k = 0; k < size; k++)
        bk[k] = b[orderOfString[k]];

    //провести над вектором преобразовния с коэффициентами треугольной матрицы
    for (int k = 0; k < size - 1; k++)
        for (int l = k + 1; l < size; l++)
            bk[l] = bk[l] - bk[k] * dat[l][k];
    return bk;
}



//обратный ход
vector<double> revStep(vector<double> b, vector<vector<double>> a) {
    int m = b.size();
    vector<double> x = vector<double>(m);
    //найти каждый х, начиная с i-того
    for (int i = m - 1; i >= 0; i--) {
        //посчитать сумму произведений известных х
        //на соответствующие элементы матрицы А
        double sum = 0;
        for (int l = i + 1; l < m; l++) {
            sum += x[l] * a[i][l];
        }
        //вычислить текущий х
        x[i] = b[i] - sum;
        x[i] = x[i] / a[i][i];
    }
    return x;
}

//поиск решения СЛАУ методом Гаусса
vector<double> GaussMethod(Matrix B, vector<double> b) {
    TriangleMatrix A = TriangleMatrix(B.dat);
    int m = b.size();
    vector<double> x;
    //проверить определитель матрицы
    if (abs(A.getDet()) < eps) {
        return x;
    }
    //преобразовать правую часть уравнения
    vector<double> bk = A.transformVector(b);
    //обратный ход метода Гаусса
    return revStep(bk, A.dat);
}


//вычисление невязки r = Ae1 - L1e1
vector<db> residual(matrix<db> A, vect<db> e1, db l1) {
    int m = e1.size();
    vect<db> res = vect<db>(m);
    vect<db> Ae1 = A * e1;
    vect<db> L1e1 = e1 * l1;
    //для каждой строки
    for (int i = 0; i < m; i++) {
        res[i] = Ae1[i] - L1e1[i];
    }
    return res;
}









#endif