// This is a personal academic project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: http://www.viva64.com

#include <iostream>
#include <boost/test/unit_test.hpp>
#include "../Lib/Math.h"
#include "../init/labs_init.h"

using namespace std;
double Func1 (const Matrix<double>& A)
{
    try
    {
        if (A.get_size() != 2)
        {
            throw 1;
        }
        double* pa = A.data();
        if (!pa)
            throw 2;

        const double result = *pa * *pa;
        ++pa;
        return result - *pa * *pa - 1;
    }
    catch (int i)
    {
        if (i == 1)
            cout << "In Func1 Matrix have size " << A.get_n() << ':' << A.get_m() << " but must have 2:2." << endl;
        if (i == 2)
            cout<<"In Func1 Matrix have null pointer."<<endl;
        return -1;
    }
}
double Func2 (const Matrix<double>& A)
{
    try
    {
        if (A.get_size() != 2)
        {
            throw 1;
        }
        double* pa = A.data();
        if (!pa)
            throw 2;
        const double result = *pa;
        ++pa;
        return result * *pa * *pa * *pa - *pa - 3;
    }
    catch (int i)
    {
        if (i == 1)
            cout << "In Func2 Matrix have size " << A.get_n() << ':' << A.get_m() << " but must have 2:2." << endl;
        if (i == 2)
            cout<<"In Func2 Matrix have null pointer."<<endl;
        return -1;
    }
}
double Func3 (const Matrix<double>& A)
{
    try
    {
        if (A.get_size() != 2)
        {
            return -1;
        }
        double* pa = A.data();
        if (!pa)
            return -1;
        return 2 * *pa;
    }
    catch (int i)
    {
        if (i == 1)
            cout << "In Func3 Matrix have size " << A.get_n() << ':' << A.get_m() << " but must have 2:2." << endl;
        if (i == 2)
            cout<<"In Func3 Matrix have null pointer."<<endl;
        return -1;
    }
}
double Func4 (const Matrix<double>& A)
{
    try
    {
        if (A.get_size() != 2)
        {
            throw 1;
        }
        double* pa = A.data();
        if (!pa)
            throw 2;
        return -2**(++pa);
    }
    catch (int i)
    {
        if (i == 1)
            cout << "In Func4 Matrix have size " << A.get_n() << ':' << A.get_m() << " but must have 2:2." << endl;
        if (i == 2)
            cout<<"In Func4 Matrix have null pointer."<<endl;
        return -1;
    }
}
double Func5 (const Matrix<double>& A)
{
    try
    {
        if (A.get_size() != 2)
        {
            throw 1;
        }
        double* pa = A.data();
        if (!pa)
            throw 2;
        ++pa;
        return *pa * *pa * *pa;
    }
    catch (int i)
    {
        if (i == 1)
            cout << "In Func5 Matrix have size " << A.get_n() << ':' << A.get_m() << " but must have 2:2." << endl;
        if (i == 2)
            cout<<"In Func5 Matrix have null pointer."<<endl;
        return -1;
    }
}
double Func6 (const Matrix<double>& A)
{
    try
    {
        if (A.get_size() != 2)
        {
            throw 1;
        }
        double* pa = A.data();
        if (!pa)
            throw 2;
        const double result = 3 * *pa;
        ++pa;
        return result * *pa * *pa - *pa;
    }
    catch (int i)
    {
        if (i == 1)
            cout << "In Func6 Matrix have size " << A.get_n() << ':' << A.get_m() << " but must have 2:2." << endl;
        if (i == 2)
            cout<<"In Func6 Matrix have null pointer."<<endl;
        return -1;
    }
}

BOOST_AUTO_TEST_SUITE(Lab_2)
        BOOST_FIXTURE_TEST_CASE(Case_for_lab_2, labs_init)
        {
            const double Eps1 = 1E-9;
            const double Eps2 = 1E-9;
            Array_of_functions_2 Ar(2);
            Ar[0] = Function_2(Func1);
            Ar[1] = Function_2(Func2);
            Array_of_functions_2 Der_of_Ar (4);
            Der_of_Ar[0] = Function_2(Func3);
            Der_of_Ar[1] = Function_2(Func4);
            Der_of_Ar[2] = Function_2(Func5);
            Der_of_Ar[3] = Function_2(Func6);
            Matrix<double> X1 (2, 1);
            X1[0] = 1;
            X1[1] = 1;
            std::ofstream log_stream("2_Lab/output/2_Lab.log");
            Matrix<double> Y1 = solve_nonlinear_equations::solve_SNE(Ar, Der_of_Ar, X1, Eps1, Eps2, log_stream);
            X1.view();
            Y1.view();
            BOOST_CHECK_CLOSE_FRACTION(Y1[0],1.6968,1E-5);
            BOOST_CHECK_CLOSE_FRACTION(Y1[1],1.37081,1E-5);
            cout << endl;
            Matrix<double> X2 (2, 1);
            X2[0] = -1;
            X2[1] = -1;
            Matrix<double> Y2 = solve_nonlinear_equations::solve_SNE(Ar, Der_of_Ar, X2, Eps1, Eps2, log_stream);
            X2.view();
            Y2.view();
            BOOST_CHECK_CLOSE_FRACTION(Y2[0],-1.47865 ,1E-5);
            BOOST_CHECK_CLOSE_FRACTION(Y2[1],-1.08922,1E-5);
        }
BOOST_AUTO_TEST_SUITE_END()
