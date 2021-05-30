// interpol.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <fstream>
#include <iostream>
#include <iomanip>
#include <map>
#include <cassert>
#include "num_methods.h"
using namespace  std;
#define inf 1e12
#define linp << setw(LINP) <<
int LINP = 12;

vector<db> getFunc(db(*func)(db x, int id), vector<db> x, int id);
vector<db> getFunc(db(*func)(db x), vector<db> x);

//результат приближения функции в заданной точке curx
db getFuncApproxInDotFromSpline(vector<vector<db>> SplineMatrix, vector<db> x, db curx);

//найти коэффициенты сплайна на каждом отрезке
void getDiffsFromStateH(ostream& fout, int iterCount, vector<AbstractIntegral*>& ai,
                        db A, db B, vector<vector<db>>& res) {
    //исследование зависимости точности вычисления
    //интеграла от выбора метода и шага
    db h = 1e-6;
    fout.precision(16);
    for (int i = 0; i < iterCount; i++, h += 8e-6)
    {
        fout << h << " ";
        int n = (int)(B - A) / h;
        vector<db> x = ai[0]->getPointsWithFixedStep(A, B, n);
        for (int m = 0; m < ai.size(); m++)
        {
            res[m][i] = ai[m]->calcIntegral(x);
        }
    }

}
void getPointsFromEps(ostream& fout, int iterCount, const vector<AbstractIntegral*>& ai,
                      const db A, const db B, vector<vector<db>>& res) {
    //исследование зависимости точности вычисления
//интеграла от выбора метода и шага
    db eps = 1e-7;
    fout.precision(16);
    for (int i = 0; i < iterCount; i++, eps += 1e-7)
    {
        fout << eps << " ";
        for (int m = 0; m < ai.size(); m++)
        {
            res[m][i] = ai[m]->calcPoints(A, B, eps).size();
        }
    }
}
//D-area
vector<db> dx = { 0, 0, 10, 10 };
vector<db> dy = { 10, 0, 0, 10 };
//figures dots
vector <db>fx = { 2, 2, 8, 8 };
vector <db>fy = { 7, 3, 3, 7 };
vector<db> t = { 1, 2, 3, 4 };
struct ParametricallyDefinedArea {
public:
    vector<db> x;
    vector<db> y;
    vector<db> t;
    vector<vector<db>> xSplineCoeffs;
    vector<vector<db>> ySplineCoeffs;
    ParametricallyDefinedArea(vector<db> x, vector<db>y) :x(x), y(y) {
        assert(x.size() == y.size());
        t = vector<db>(x.size());
        for (int i = 0; i < x.size(); i++)
            t[i] = i + 1.;
        xSplineCoeffs = getSplineMatrix(t, x);
        ySplineCoeffs = getSplineMatrix(t, y);
    }
    void Print(ostream& fout, db eps) {
        //исходные точки
        PrintFunc(fout, x, y);
        //точки для построения графика
        vector<db> printDots = autoGen(eps * t[0], t[0], t[t.size() - 1]);
        printVec(getSplineRes(t, x, printDots, xSplineCoeffs), fout);
        fout << endl;
        printVec(getSplineRes(t, y, printDots, ySplineCoeffs), fout);
        fout << endl;

        return;
    }
};
class Point {
public:
    db x = 0, y = 0;
    Point(db x, db y) :x(x), y(y) {}
};
vector<Point> createExternalRectangleArea(ParametricallyDefinedArea fig) {
    vector<Point> coord = vector<Point>();
    db delta = 2;
    db minx, miny, maxx, maxy;
    minx = mMin(fig.x) - delta;
    miny = mMin(fig.y) - delta;
    maxx = mMax(fig.x) + delta;
    maxy = mMax(fig.y) + delta;
    coord.push_back(Point(minx, maxy));
    coord.push_back(Point(minx, miny));
    coord.push_back(Point(maxx, miny));
    coord.push_back(Point(maxx, maxy));
    coord.push_back(coord[0]);
    return coord;
}
std::ostream& operator<< (std::ostream& out, Point& one) {
    out << one.x << " " << one.y;
    return out;
}
void PaintPoints(ostream& fout, vector<Point> points) {
    for (int i = 0; i < points.size(); i++)
        fout << points[i].x << " ";
    fout << endl;
    for (int i = 0; i < points.size(); i++)
        fout << points[i].y << " ";
    fout << endl;
}
int printGr()
{
    db eps = 1e4;
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    ParametricallyDefinedArea fig1 = ParametricallyDefinedArea({ 2, 2, 8, 8, 2 }, { 7, 3, 3, 7, 7 });
    ParametricallyDefinedArea fig2 = ParametricallyDefinedArea({ 2, 2, 2,5, 8,8, 8,5, 2 },
                                                               { 7, 5, 3,3, 3, 5, 7, 7, 7 });
    vector<Point> d = createExternalRectangleArea(fig1);
    PaintPoints(fout, d);
    fig1.Print(fout, eps);
    fout.close();
    fin.close();

    system("python integral.py");
    return 0;
}


int main() {
    setlocale(LC_ALL, "rus");
    printGr();
}
