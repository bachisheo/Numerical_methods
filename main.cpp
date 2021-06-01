// interpol.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <fstream>
#include "num_methods.h"

using namespace std;
#define inf 1e12

void printPoints (ostream& out, vector<point<db>> points) ;

//figures dots

vect<db> fx = vector<db>{2, 2, 2, 5, 8, 8, 8, 5, 2};
vect<db> fy = vector<db>{7, 5, 3, 3, 3, 5, 7, 7, 7};
vect<db> ft = vector<db>{1, 2, 3, 4, 5, 6, 7, 8, 9};

int main() {
    db eps = 1e4;
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



    cout <<  mc.calcArea(10000) << endl;
    //print area built by Spline and all random dots
    printPoints(fout, mc.areaDots);
    printPoints(fout, mc.dotsIn);
    printPoints(fout, mc.dotsOut);


    //calculate area with SimpsonsMethod
    cout << calcParamAreabySimpson(tGr, xGr, yGr);
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

