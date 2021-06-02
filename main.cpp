// interpol.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <fstream>
#include "num_methods.h"

using namespace std;
#define inf 1e12


void printPoints (ostream& out, vector<point<db>> points) ;

//figures dots

vect<db> fx = vector<db>{0.278, 0.712, 0.898, 0.688, 0.244, 0.092, 0.302, 0.614, 0.624, 0.494, 0.336, 0.462, 0.332,
0.22, 0.308, 0.52, 0.708, 0.722, 0.586, 0.336, 0.116, 0.056, 0.278};
        //{0.336, 0.144, 0.15, 0.302, 0.522, 0.612, 0.644, 0.414, 0.336};
        //{2, 2, 2, 5, 8, 8, 8, 5, 2};
vect<db> fy = vector<db>{0.896, 0.868, 0.478, 0.116, 0.102, 0.478, 0.754, 0.662, 0.358, 0.272, 0.406, 0.522, 0.606,
        0.388, 0.21, 0.176, 0.328, 0.646, 0.802, 0.8, 0.754, 0.878, 0.896};
        //{0.82, 0.654, 0.456, 0.382, 0.352, 0.452, 0.626, 0.66, 0.82};
        //{7, 5, 3, 3, 3, 5, 7, 7, 7};
vect<db> ft = vector<db>{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23};
//{1, 2, 3, 4, 5, 6, 7, 8, 9};

int main() {
    db eps = 1e4;
    int numbDots = 1000;
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    //print origin area dots
    fout << fx << endl << fy << endl;

    ParametricallyDefinedArea fig1 = ParametricallyDefinedArea(fx, fy);

    //create dots for figures approximate
    vect<db> xGr, yGr, tGr = autoGen(1e4, fig1.t.min(), fig1.t.max());
    xGr = fig1.xSpline.getFuncApproxInDots(tGr);
    yGr = fig1.ySpline.getFuncApproxInDots(tGr);

    MonteCarlo mc = MonteCarlo(xGr, yGr);
    mc.printExternalArea(fout);



   db mcArea = mc.calcDoubleIntegral(numbDots, [](db x, db y) { return 1; });
    //print area built by Spline and all random dots
    printPoints(fout, mc.areaDots);
    printPoints(fout, mc.dotsIn);
    printPoints(fout, mc.dotsOut);




    db xDInt = mc.calcDoubleIntegral(numbDots, [](db x, db y) { return x; });
    db yDInt = mc.calcDoubleIntegral(numbDots, [](db x, db y) { return y; });
    fout  << xDInt/mcArea << endl << yDInt/mcArea << endl;

    //calculate area with SimpsonsMethod
    fout << calcParamAreabySimpson(tGr, xGr, yGr) << endl << mcArea;
    fout.close();
    fin.close();
    system("python ..\\integral.py");
    return 0;
}
void printPoints (ostream& out, vector<point<db>> points) {
    for (int i = 0; i < points.size(); ++i)
        out << points[i].x << " ";
    out << endl;
    for (int i = 0; i < points.size(); ++i)
        out << points[i].y << " ";
    out << endl;
}

