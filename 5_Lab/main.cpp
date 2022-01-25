// This is a personal academic project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: http://www.viva64.com

#include <iostream>
#include <boost/test/unit_test.hpp>
#include "../Lib/Math.h"
using namespace std;

inline double f1(const double x)
{
    return 3*x/sqrt(1+x*x*x);
}
inline double f2 (const Matrix<double>& X)// 32
{
    return 4 - *X.data() * *X.data() - *(X.data() + 1) * *(X.data() + 1);
}

BOOST_AUTO_TEST_SUITE(Lab_5)
    BOOST_AUTO_TEST_CASE(Case_for_lab_5)
    {
        Function F (f1);
        double alg_1 = F.integral_by_trapeze_method(0, 1.075, 10e-5);
        double alg_2 = F.integral_by_Sympthon_method(0, 1.075, 10e-5);
        cout << "Trapetia " << alg_1 << "   " << "Sympthon " << alg_2 << endl;
        BOOST_CHECK_CLOSE_FRACTION(alg_1, 1.45, 1E-2);
        BOOST_CHECK_CLOSE_FRACTION(alg_2, 1.45, 1E-2);
        Function_2 F2 (f2);
        Matrix<double> from (2, 1);
        Matrix<double> to (2, 1);
        from[0] = from[1] = -1;
        to[0] = to[1] = 1;
        double twice_integral = F2.integral_by_definition(from, to);
        cout << twice_integral << endl;
        BOOST_CHECK_CLOSE_FRACTION(twice_integral, 13.33333, 1E-5);
        //1.44982 37.33
    }
BOOST_AUTO_TEST_SUITE_END()
