#include <iostream>
#include "../Lib/Math.h"
using namespace std;

int main()
{
    double pa[] = {4, 2, 3, 2, 6, 5, 3, 5, 10}; //Work!
    Matrix<double> A (pa, 3, 3);
    cout<<A.If_Symmetric()<<endl;
    Matrix_SLE B (A);
    Matrix <double> lB(3, 1);
    lB[0] = 1;
    lB[1] = 2;
    lB[2] = 3;
    B.View();
    Matrix<double> C = B.Solve(lB);
    return 0;
}
