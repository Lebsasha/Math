#include <iostream>
#include <vector>
#include <time.h>
#include <iomanip>
#include <string>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <limits>
#include <iterator>
#include "../Lib/Math.h"
#define nullptr NULL
#define endl "\n"
#include <locale.h>
using namespace std;
void View (const vector<double>&);
vector<double> Minus (vector<double>& a, vector<double>& b)
{
    if (a.size() != b.size())
        return vector<double> (0);
    vector<double> c (a.size());
    for (vector<double>::iterator iva = a.begin(), ivb = b.begin(), ivc = c.begin(); iva < a.end(); iva++, ivb++, ivc++)
    {
        *ivc = *iva - *ivb;
    }
    return c;
}
int main ()
{
    const int N = 3;
    Matrix_SLE A;
    A.Read_from_file("A.matr_l");
    //A.View();
    vector<double> vb (N);
    vb[0] = 3;
    vb[1] = 3.8;
    vb[2] = 77;
    Matrix_SLE B(A);
    B.View();
    B = A.T();
    B.View();
    vector<double> vb1 (vb);
    View (vb1);
    vector<double> Solucion1 = A.Solve_by_Gauss_Method_v(vb);
    if (Solucion1.empty())
        cout<<"Haven't sol.";
    cout<<"Solution 1:"<<endl;
    View (Solucion1);
    vb = vb1;
    A = B;
    vb1 = A.Multi_outp_vector(Solucion1);
    vector<double> F = Minus(vb1, vb);
    setlocale (LC_ALL, "grc");
    cout<<"Δ Norm = "<<Norm(F)<<endl;
    vector<double> Solution2 = B.Solve_by_Gauss_Method_v(vb1);
    cout<<"Solution 2:"<<endl;
    View (Solution2);
    double delta = Norm (Minus(Solution2, Solucion1))/Norm(Solucion1);
    cout<<"δ delta = "<<delta;


//    A.Read_from_file("A.matr_l");//B.matr_l
//    Matrix<double> b (vb);
//    B = A.Solve(b);
    return 0;
}
void View (const vector<double>& El1)
{
    for (vector<double>::const_iterator i = El1.begin(); i < El1.end(); ++i)
        cout<<*i<<endl;
    return;
}
/*void* operator new (size_t size)
{
    Return the name of the object;
}*/
void randoe(int*, const int, const int, int=10);
/*int main()
{
tm a11;
long int i1 = 0;
time(&i1);

    string path("/home/alexander/Documents/Projects/Temp/1");
    //getline (cin, path);
    int i = path.find_last_of('/');
    string Name;
    Name.insert(0, path, i + 1, path.length() - i);
    path = "xterm -T " + Name + " -e /usr/bin/cb_console_runner LD_LIBRARY_PATH=$LD_LIBRARY_PATH:. " + path;
    //system ("gnome-terminal");
    //xterm -T 1 -e /usr/bin/cb_console_runner LD_LIBRARY_PATH=$LD_LIBRARY_PATH:. /home/alexander/Documents/Projects/Temp/1
    //system("/home/alexander/Desktop/apt-upd > qwqwqwqwwqwqwqwqww");
    //system ("rm qwqwqwqwwqwqwqwqww");
    system (path.c_str());
    //system ("read Dummy");
    return 0;
}*/
