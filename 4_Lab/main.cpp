#include <iostream>
#include <algorithm>
#define CB 0
//#define VISUAL_STUDIO
//#define NDEBUG 1
#include <assert.h>
#include <regex>
#include <math.h>
#include "../Lib/Different.h"
#include "../Lib/Math.h"
int main()
{
    const int N = 6;
    double pa[] = {0, 20, 40, 60, 80, 100};
    Matrix<double> X (pa, N, 1);
    double py[] = {log10(29.5), log10(18.4), log10(11.9), log10(8.6), log10(5.0), log10(3.3)};
    Matrix<double> Y (py, N, 1);
    Matrix<double> F = Approximation::Find_Polinom_m_power(X, Y, 1);
    F.View();
    return 0;
}
/*int _tmain(int argc, _TCHAR* argv[])
{
	return 0;
}*/
