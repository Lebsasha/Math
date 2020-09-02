#include <math.h>
#include <boost/test/unit_test.hpp>
//#define NDEBUG 1
#include "../Lib/Math.h"

int main_for_Lab_4();

BOOST_AUTO_TEST_SUITE(SuItE_tests_for_Lab_4)
    BOOST_AUTO_TEST_CASE(Case_for_lab_4)
    {
        BOOST_CHECK(main_for_Lab_4()==0);
    }
BOOST_AUTO_TEST_SUITE_END()

int main_for_Lab_4 ()
{
    const int N = 6;
    double pa[] = {0, 20, 40, 60, 80, 100};
    Matrix<double> X (pa, N, 1);
    double py[] = {log10(29.5), log10(18.4), log10(11.9), log10(8.6), log10(5.0), log10(3.3)};
    Matrix<double> Y (py, N, 1);
    Matrix<double> F = Approximation::Find_Polinom_m_power(X, Y, 1);
    F.view();
    return 0;
}
/*int _tmain(int argc, _TCHAR* argv[])
{
	return 0;
}*/
