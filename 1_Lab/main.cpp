// This is a personal academic project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: http://www.viva64.com

#include <iostream>
#include <vector>
#include <limits>
#include <clocale>
#include <boost/test/unit_test.hpp>
#include "../Lib/Math.h"
#include "../init/labs_init.h"

#if __cplusplus  < 201103L
#define nullptr NULL
#define endl "\n"
#endif
using namespace std;
void View (const vector<double>& El1);
vector<double> Minus (const vector<double>& a, const vector<double>& b)
{
    if (a.size() != b.size())
        return vector<double> (0);
    vector<double> c (a.size());
    auto ivc = c.begin();
    for (vector<double>::const_iterator iva = a.cbegin(), ivb = b.cbegin(); iva < a.end(); ++iva, ++ivb, ++ivc)
    {
        *ivc = *iva - *ivb;
    }
    return c;
}

BOOST_AUTO_TEST_SUITE(Lab_1)
BOOST_FIXTURE_TEST_CASE(Case_for_lab_1, labs_init)
        {
            const int N = 3;
            Matrix_SLE matrix_1;
            matrix_1.read_from_file("./1_Lab/input/data.matr_l");
            matrix_1.view();
            Matrix<double> b (N,1);
            b[0] = 3;
            b[1] = 3.8;
            b[2] = 77;
            Matrix_SLE matrix_1_copy(matrix_1);
            Matrix<double> vb2 (b);
            Matrix<double> solution1 = matrix_1.solve_by_Gauss_method(b);
            BOOST_CHECK(solution1.get_size() != 0);
            for (auto p_curr_1 = matrix_1.end(), p_curr = matrix_1_copy.end(); p_curr_1 >= matrix_1.begin(); --p_curr_1, --p_curr)
                BOOST_CHECK(*p_curr_1==*p_curr);
            BOOST_CHECK(vb2[0]==b[0]);
            BOOST_CHECK(vb2[1]==b[1]);
            BOOST_CHECK(vb2[2]==b[2]);
            cout<<"Solution 1:"<<endl;
            solution1.view();
            BOOST_CHECK_CLOSE_FRACTION(solution1[0], -148.73333, 1E-5);
            BOOST_CHECK_CLOSE_FRACTION(solution1[1], 74.73333, 1E-5);
            BOOST_CHECK_CLOSE_FRACTION(solution1[2], 2.26666, 1E-5);
            b = vb2;
            vector<double> vb1 = static_cast<vector<double> >(matrix_1 * Matrix<double>(solution1));
            vector<double> F = Minus(vb1, static_cast<vector<double> >(b));
            setlocale (LC_ALL, "grc");
            const double norm_F = norm(F);
            cout << "Δ Norm = " << norm_F << endl;
            BOOST_CHECK_CLOSE_FRACTION(norm_F, -1.33226E-14, 1E-5);
            Matrix<double> solution2 = matrix_1.solve_by_Gauss_method(vb2);
            cout<<"Solution 2:"<<endl;
            solution2.view();
            BOOST_CHECK_CLOSE_FRACTION(solution2[0], -148.73333, 1E-5);
            BOOST_CHECK_CLOSE_FRACTION(solution2[1], 74.73333, 1E-5);
            BOOST_CHECK_CLOSE_FRACTION(solution2[2], 2.26666, 1E-5);

            double delta = norm(Minus(static_cast<vector<double> >(solution2), static_cast<vector<double> >(solution1))) / norm
                    (static_cast<vector<double> >(solution1));
            cout<<"δ delta = "<<delta<<endl;
            BOOST_CHECK_CLOSE_FRACTION(delta, 0, DBL_EPSILON);

            // LDLT algorithm checking
            //matrix_1.read_from_file("matrix_1.matr_l");//solution3.matr_l
            double p[] = {1, 2, 2, 1};
            Matrix_SLE C (p, 2, 2);
            p[0] = 1;
            p[1] = 2;
            Matrix<double> c (p, 2, 1);
            Matrix<double> solution3 = C.solve(c);
            solution3.view();
            BOOST_CHECK_CLOSE_FRACTION(solution3[0], 1, DBL_EPSILON);
            BOOST_CHECK_CLOSE_FRACTION(solution3[1], 0, DBL_EPSILON);
        }
BOOST_AUTO_TEST_SUITE_END()

void View (const vector<double>& El1)
{
    for (vector<double>::const_iterator i = El1.begin(); i < El1.end(); ++i)
        cout<<*i<<endl;
}
