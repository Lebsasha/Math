// *** ADDED BY HEADER FIXUP ***
#include <istream>
// *** END ***
#include <iostream>
//#include "/media/alexander/loadtool/Num_Methods/Lib/Matrix.h"
#include "../Lib/Different.h"
#include "../Lib/Math.h"
using namespace std;
double Func1 (const Matrix<double>& A)
{
    try
    {
        if (A.Get_Size() != 2)
        {
            throw 1;
        }
        double* pa = A.Get_pointer();
        if (!pa)
            throw 2;
        return *pa**pa - *(++pa)**pa - 1;
    }
    catch (int i)
    {
        if (i == 1)
            cout<<"In Func1 Matrix have size "<<A.Get_N()<<':'<<A.Get_M()<<" but must have 2:2."<<endl;
        if (i == 2)
            cout<<"In Func1 Matrix have null pointer."<<endl;
        return -1;
    }
}
double Func2 (const Matrix<double>& A)
{
    try
    {
        if (A.Get_Size() != 2)
        {
            throw 1;
        }
        double* pa = A.Get_pointer();
        if (!pa)
            throw 2;
        return *pa**(++pa)**pa**pa - *pa - 3;
    }
    catch (int i)
    {
        if (i == 1)
            cout<<"In Func2 Matrix have size "<<A.Get_N()<<':'<<A.Get_M()<<" but must have 2:2."<<endl;
        if (i == 2)
            cout<<"In Func2 Matrix have null pointer."<<endl;
        return -1;
    }
}
double Func3 (const Matrix<double>& A)
{
    try
    {
        if (A.Get_Size() != 2)
        {
            return -1;
        }
        double* pa = A.Get_pointer();
        if (!pa)
            return -1;
        return 2**pa;
    }
    catch (int i)
    {
        if (i == 1)
            cout<<"In Func3 Matrix have size "<<A.Get_N()<<':'<<A.Get_M()<<" but must have 2:2."<<endl;
        if (i == 2)
            cout<<"In Func3 Matrix have null pointer."<<endl;
        return -1;
    }
}
double Func4 (const Matrix<double>& A)
{
    try
    {
        if (A.Get_Size() != 2)
        {
            throw 1;
        }
        double* pa = A.Get_pointer();
        if (!pa)
            throw 2;
        return -2**(++pa);
    }
    catch (int i)
    {
        if (i == 1)
            cout<<"In Func4 Matrix have size "<<A.Get_N()<<':'<<A.Get_M()<<" but must have 2:2."<<endl;
        if (i == 2)
            cout<<"In Func4 Matrix have null pointer."<<endl;
        return -1;
    }
}
double Func5 (const Matrix<double>& A)
{
    try
    {
        if (A.Get_Size() != 2)
        {
            throw 1;
        }
        double* pa = A.Get_pointer();
        if (!pa)
            throw 2;
        return *(++pa)**pa**pa;
    }
    catch (int i)
    {
        if (i == 1)
            cout<<"In Func5 Matrix have size "<<A.Get_N()<<':'<<A.Get_M()<<" but must have 2:2."<<endl;
        if (i == 2)
            cout<<"In Func5 Matrix have null pointer."<<endl;
        return -1;
    }
}
double Func6 (const Matrix<double>& A)
{
    try
    {
        if (A.Get_Size() != 2)
        {
            throw 1;
        }
        double* pa = A.Get_pointer();
        if (!pa)
            throw 2;
        return 3**pa**(++pa)**pa - *pa;
    }
    catch (int i)
    {
        if (i == 1)
            cout<<"In Func6 Matrix have size "<<A.Get_N()<<':'<<A.Get_M()<<" but must have 2:2."<<endl;
        if (i == 2)
            cout<<"In Func6 Matrix have null pointer."<<endl;
        return -1;
    }
}
#define next_linee cout<<endl;
int main(void)
{
    const double Eps1 = 1E-9;
    const double Eps2 = 1E-9;
    Array_of_Functions2 Ar(2);
    Ar[0] = Func1;
    Ar[1] = Func2;
    Array_of_Functions2 Der_of_Ar (4);
    Der_of_Ar[0] = Func3;
    Der_of_Ar[1] = Func4;
    Der_of_Ar[2] = Func5;
    Der_of_Ar[3] = Func6;
    Matrix<double> X1 (2, 1);
    X1[0] = 1;
    X1[1] = 1;
    Matrix<double> Y1 = Solve_SNE(Ar, Der_of_Ar, X1, Eps1, Eps2);
    X1.View();
    Y1.View();
    next_linee
    Matrix<double> X2 (2, 1);
    X2[0] = -1;
    X2[1] = -1;
    Matrix<double> Y2 = Solve_SNE(Ar, Der_of_Ar, X2, Eps1, Eps2);
    X2.View();
    Y2.View();
    cout << "Hello world!" << endl;
    Get_Pause();
    return 0;
}
