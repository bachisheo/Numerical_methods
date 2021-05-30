//
// Created by kesa on 29.05.2021.
//

#include "matrix.h"

#ifndef SIMPLEXMETHOD_MYMATRIX_H
#define SIMPLEXMETHOD_MYMATRIX_H


TriangleMatrix::TriangleMatrix(vect<vect<db>> M) {
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
};

//вычислить определитель матрицы
double TriangleMatrix::calcDet() {
    assert(rowNumb() == colNumb());
    det = 1;
    //вычислить определитель треугольной матрицы
    //как произведение элементов, стоящих на главной диагонали
    for (int k = 0; k < rowNumb(); k++) {
        if (abs(dat[k][k]) < eps)
            return 0;
        det *= dat[k][k];
    }
    //изменить знак определителя, если строки
    //были переставлены четное число раз
    if (numbOfTurns % 2)
        det *= (-1);
    return det;
}

//преобразовать вектор b к трекгольному виду матрицы
vect<db> TriangleMatrix::transformVector(vect<db> b) {
    vector<db> bk;
    int size = rowNumb();
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

void TriangleMatrix::swapString(int a, int b) {
    //зафиксировать изменение в порядке строк
    swap(orderOfString[a], orderOfString[b]);
    //поменять местами содержимое строк
    vect<db> tmp = dat[a];
    dat[a] = dat[b];
    dat[b] = tmp;
    numbOfTurns++;
}

//обратный ход
vector<double> revStep(vector<double> b, vect<vect<double>> a) {
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
vector<double> GaussMethod(matrix<db> B, vector<double> b) {
    TriangleMatrix A = TriangleMatrix(B.dat);
    int m = b.size();
    vector<double> x;
    //проверить определитель матрицы
    if (abs(A.calcDet()) < eps) {
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