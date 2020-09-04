#include <math.h>
#include <boost/test/unit_test.hpp>
//#define NDEBUG 1
#include "../Lib/Math.h"


BOOST_AUTO_TEST_SUITE(SuItE_tests_for_Lab_4)
    BOOST_AUTO_TEST_CASE(Case_for_lab_4)
    {
        const int N = 6;
        double pa[] = {0, 20, 40, 60, 80, 100};
        Matrix<double> X (pa, N, 1);
        double py[] = {log10(29.5), log10(18.4), log10(11.9), log10(8.6), log10(5.0), log10(3.3)};
        Matrix<double> Y (py, N, 1);
        Matrix<double> F = Approximation::Find_Polinom_m_power(X, Y, 1);
        BOOST_CHECK_CLOSE_FRACTION(F[0], -0.90471, 1E-5);
        BOOST_CHECK(fabs(F[1] - 0.02289) < 1E-5);
        F = Approximation::Find_Polinom_m_power(X, Y, 5);
        BOOST_CHECK(fabs(F[0] - -0.36746) < 1E-5);
        BOOST_CHECK_CLOSE_FRACTION(F[1], 0.20652, 1E-5);
        BOOST_CHECK(fabs(F[2] - -0.00932) < 1E-5);
        BOOST_CHECK(fabs(F[3] - 0.00018) < 1E-5);
        BOOST_CHECK(fabs(F[4] - -1.72446e-06) < 1E-6);
        BOOST_CHECK(fabs(F[5] - 6.00188e-09) < 1E-9);
    }
BOOST_AUTO_TEST_SUITE_END()

/*int _tmain(int argc, _TCHAR* argv[])
{
	return 0;
}*/
