// interpol.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <fstream>
#include "num_methods.h"

using namespace std;
#define inf 1e12

//figures dots

int main() {
    ifstream fin("input.txt");
    ofstream fout("output.txt");
    /*CombainMethod cm = CombainMethod([](db x) {return cos(x) + 1./(x*x*x);},
                                     [](db x) {return -1.*sin(x) - 3./(x*x);},
                                     [](db x) {return -1.*cos(x) +6./x;});
    */CombainMethod cm = CombainMethod([](db x) {return sin(x);},
                                     [](db x) {return cos(x);},
                                     [](db x) {return -1.*sin(x);});
    cm.calcAllRoots(0, 7, 1e-8, 1, cout);

    fout.close();
    fin.close();
   // system("python ..\\integral.py");
    return 0;
}


