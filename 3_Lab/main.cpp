// This is a personal academic project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: http://www.viva64.com

#include <iostream>
#include <boost/test/unit_test.hpp>
#include "../Lib/Math.h"
using namespace std;
extern string path;

inline double Diff_u1 (const Matrix<double>& X)
{
    assert (X.get_size() == 3);
    double x = *(X.data() + 2);
    double fff = (x > 1e-10) ? sin(x)/x : 1;
    return -*X.data() * *(X.data() + 1) + fff;
}
inline double Diff_u2 (const Matrix<double>& X)
{
    assert (X.get_size() == 3);
    double x = *(X.data() + 2);
    return (-Pow(*(X.data() + 1), 2) + (2.5 + 35 / 40) * x / (1 + x * x));
}

inline double Diff_u1_sp (const Matrix<double>& X)
{
    assert (X.get_size() == 6);
    double x = *(X.data() + 4);
    double fff = (x > 1e-10) ? sin(x)/x : 1;
    return X.unsafe_index_c(0) - X.unsafe_index_c(2) - (-X.unsafe_index_c(0) * X.unsafe_index_c(1) + fff) * X.unsafe_index_c(5);
}
inline double Diff_u2_sp (const Matrix<double>& X)
{
    assert (X.get_size() == 6);
    double x = *(X.data() + 4);
    return X.unsafe_index_c(1) - X.unsafe_index_c(3) - (-Pow(X.unsafe_index_c(1), 2) + (2.5 + 35 / 40) * x / (1 + x * x)) * X.unsafe_index_c(5);
}

BOOST_AUTO_TEST_SUITE(SuItE_tests_for_Lab_3___)
BOOST_AUTO_TEST_CASE(Case_for_lab_3)
        {
            Array_of_functions_2 A (2);
            A[0] = Function_2(Diff_u1);
            A[1] = Function_2(Diff_u2);
            Matrix<double> u0 (2, 1);
            u0[0] = 0;
            u0[1] = -0.412;
            Matrix<double> Eps (2, 1);
            Eps[1] = Eps[0] = 1e-3;
            const double max_step = 0.01;
            vector<Matrix<double> > yk = solve_differential_equations::explicit_Euler_method(A, 0, 1, u0, Eps, max_step);
            ofstream oFile_Exp(path + "Explicit.txt");
            auto iEnd_e = yk[0].begin();
            bool first_not_null = true;
            for (auto iyk = yk[1].end(), itk = yk[0].end(); itk >= iEnd_e; --iyk, --itk)
            {
                if (fabs(*itk) > DBL_EPSILON || fabs(*iyk) > DBL_EPSILON || fabs(*(iyk - 1)) > DBL_EPSILON)
                {
                    if(first_not_null)
                    {
                        BOOST_CHECK_CLOSE_FRACTION(*itk, 1, max_step);
                        BOOST_CHECK_CLOSE_FRACTION(*iyk, 0.37216, 1E-5);
                        BOOST_CHECK_CLOSE_FRACTION(*(iyk-1), 0.94716, 1E-5);
                        first_not_null=!first_not_null;
                    }
                    oFile_Exp<<*itk<< ' '<<*iyk<<' ';
                    oFile_Exp<<*(--iyk)<<'\n';
                }
                else
                    --iyk;
            }
            A[0] = Function_2(Diff_u1_sp);
            A[1] = Function_2(Diff_u2_sp);
            const double t_max = (1.0 - 0) / 10;
            vector<Matrix<double> > I_all = solve_differential_equations::implicit_Euler_method(A, 0, 1, u0, Eps, 1e-2, t_max);
            ofstream oFile_Imp(path + "Implicit.txt");
            auto iEnd_i = (I_all[1]).begin();
            first_not_null=true;
            for (auto iyk = (I_all[1]).end(), itk = (I_all[0]).end(); iyk >= iEnd_i; --iyk, --itk)
            {
                if (fabs(*itk) > DBL_EPSILON || fabs(*iyk) > DBL_EPSILON || fabs(*(iyk - 1)) > DBL_EPSILON)
                {
                    if(first_not_null)
                    {
                        BOOST_CHECK_CLOSE_FRACTION(*itk, 1, t_max);
                        BOOST_CHECK_CLOSE_FRACTION(*iyk, 0.69461, 1E-5);
                        BOOST_CHECK_CLOSE_FRACTION(*(iyk-1), 1.04864, 1E-5);
                        first_not_null=!first_not_null;
                    }
                    oFile_Imp<<*itk<<' '<<*iyk<<' ';
                    oFile_Imp<<*(--iyk)<<'\n';
                }
                else
                    --iyk;
            }
            oFile_Exp.close();
            oFile_Imp.close();
        }
BOOST_AUTO_TEST_SUITE_END()

//drand48();

