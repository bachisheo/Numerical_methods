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

    fout.close();
    fin.close();
    system("python ..\\integral.py");
    return 0;
}


