#include <iostream>
#include <boost/test/unit_test.hpp>
#include "../Lib/Math.h"
using namespace std;
extern string Path;

inline double Diff_u1 (const Matrix<double>& X)
{
    assert (X.Get_Size() == 3);
    double x = *(X.Get_pointer()+2);
    double fff = (x > 1e-10) ? sin(x)/x : 1;
    return -*X.Get_pointer()**(X.Get_pointer()+1) + fff;
}
inline double Diff_u2 (const Matrix<double>& X)
{
    assert (X.Get_Size() == 3);
    double x = *(X.Get_pointer()+2);
    return (-Pow(*(X.Get_pointer()+1), 2) + (2.5+35/40)*x/(1+x*x));
}

inline double Diff_u1_sp (const Matrix<double>& X)
{
    assert (X.Get_Size() == 6);
    double x = *(X.Get_pointer()+4);
    double fff = (x < 1e-5) ? sin(x)/x : 1.0;
    return X.Unsafe_index_c(0) - X.Unsafe_index_c(2) -(-X.Unsafe_index_c(0)*X.Unsafe_index_c(1) + fff)*X.Unsafe_index_c(5);
}
inline double Diff_u2_sp (const Matrix<double>& X)
{
    assert (X.Get_Size() == 6);
    double x = *(X.Get_pointer()+4);
    return X.Unsafe_index_c(1) - X.Unsafe_index_c(3) - (-Pow(X.Unsafe_index_c(1), 2) + (2.5+35/40)*x/(1+x*x))*X.Unsafe_index_c(5);
}
int main_for_Lab_3();

BOOST_AUTO_TEST_SUITE(SuItE_tests_for_Lab_3___)
BOOST_AUTO_TEST_CASE(Case_for_lab_3)
        {
                BOOST_CHECK(main_for_Lab_3()==0);
        }
BOOST_AUTO_TEST_SUITE_END()

int main_for_Lab_3 ()
{
    Array_of_Functions2 A (2);
    A[0] = Diff_u1;
    A[1] = Diff_u2;
    Matrix<double> u0 (2, 1);
    u0[0] = 0;
    u0[1] = -0.412;
    Matrix<double> Eps (2, 1);
    Eps[1] = Eps[0] = 1e-3;
    vector<Matrix<double> > yk = Solve_Differential_Equations::Explicit_Euler_method (A, 0, 1, u0, Eps, 0.01);
    ofstream oFile_Exp(Path + "Explicit.txt");
    auto iEnd_e = yk[0].First_i();
    for (auto iyk = yk[1].Last_i(), itk = yk[0].Last_i(); itk >= iEnd_e; --iyk, --itk)
    {
        if (fabs(*itk) > DBL_EPSILON || fabs(*iyk) > DBL_EPSILON || fabs(*(iyk - 1)) > DBL_EPSILON)
            oFile_Exp<<*itk<< ' '<<*iyk<<' '<<*(iyk--)<<'\n'; //Why undefined?
        else
            --iyk;
    }
    A[0] = Diff_u1_sp;
    A[1] = Diff_u2_sp;
    vector<Matrix<double> > I_all = Solve_Differential_Equations::Implicit_Euler_method (A, 0, 1, u0, Eps, 1e-2, (1.0-0)/10);
    ofstream oFile_Imp(Path+"Implicit.txt");
    auto iEnd_i = (I_all[1]).First_i();
    for (auto iyk = (I_all[1]).Last_i(), itk = (I_all[0]).Last_i(); iyk >= iEnd_i; --iyk, --itk)
    {
        if (fabs(*itk) > DBL_EPSILON || fabs(*iyk) > DBL_EPSILON || fabs(*(iyk - 1)) > DBL_EPSILON)
        oFile_Imp<<*itk<<' '<<*iyk<<' '<<*(iyk--)<<'\n';
        else
            --iyk;
    }
    oFile_Exp.close();
    oFile_Imp.close();
//    double p[] = {1, 2, 2, 1};
//    Matrix_SLE S(p, 2, 2);
//    S.Solve(u0).View();
    return 0;
}

//drand48();

