// interpol.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <fstream>
#include "num_methods.h"

using namespace std;
#define inf 1e12

vector<db> dx = {2, 2, 8, 8, 2};
vector<db> dy = {7, 3, 3, 7, 7};
//figures dots

vector<db> fx{2, 2, 2, 5, 8, 8, 8, 5, 2};
vector<db> fy{7, 5, 3, 3, 3, 5, 7, 7, 7};

int main() {
    db eps = 1e4;
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    ParametricallyDefinedArea fig1 = ParametricallyDefinedArea(fx, fy);

    //paint bounding area
    vector<Point<db>> d = createExternalRectangleArea(fig1);
    fout << d;

    vect<db> dotsForGraph = autoGen(1e4, fig1.t.min(), fig1.t.max());
    fout << fig1.xSpline.getFuncApproxInDots(dotsForGraph) << endl << fig1.ySpline.getFuncApproxInDots(dotsForGraph);
    fout.close();
    fin.close();

    system("python ..\\integral.py");
    return 0;
}


