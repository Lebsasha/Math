#include <iostream>
#include <vector>
#include <limits>
#include <clocale>
#include <boost/test/unit_test.hpp>
#include "../Lib/Math.h"
#if __cplusplus  < 201103L
#define nullptr NULL
#define endl "\n"
#endif
using namespace std;
void View (const vector<double>&);
vector<double> Minus (vector<double>& a, vector<double>& b)
{
    if (a.size() != b.size())
        return vector<double> (0);
    vector<double> c (a.size());
    for (vector<double>::iterator iva = a.begin(), ivb = b.begin(), ivc = c.begin(); iva < a.end(); iva++, ivb++, ivc++)
    {
        *ivc = *iva - *ivb;
    }
    return c;
}

BOOST_AUTO_TEST_SUITE(SuItE_tests_for_Lab_1___)
BOOST_AUTO_TEST_CASE(Case_for_lab_1)
        {
            const int N = 3;
            Matrix_SLE matrix_1;
            matrix_1.Read_from_file("../1_Lab/data.matr_l");
            vector<double> b (N);
            b[0] = 3;
            b[1] = 3.8;
            b[2] = 77;
            Matrix_SLE temp_matrix(matrix_1);
            temp_matrix.View();
            vector<double> vb1 (b);
            vector<double> Solucion1 = matrix_1.Solve_by_Gauss_Method_v(b);
            if (Solucion1.empty())
                cout<<"Haven't sol.";
            cout<<"Solution 1:"<<endl;
            View (Solucion1);
            BOOST_CHECK_CLOSE_FRACTION(Solucion1[0], -148.73333, 1E-5);
            BOOST_CHECK_CLOSE_FRACTION(Solucion1[1], 74.73333, 1E-5);
            BOOST_CHECK_CLOSE_FRACTION(Solucion1[2], 2.26666, 1E-5);
            b = vb1;
            matrix_1 = temp_matrix;
            vb1 = static_cast<vector<double> >(matrix_1 * Matrix<double>(Solucion1));
            vector<double> F = Minus(vb1, b);
            setlocale (LC_ALL, "grc");
            const double norm = Norm(F);
            cout << "Δ Norm = " << norm << endl;
            BOOST_CHECK_CLOSE_FRACTION(norm, -1.33226E-14, 1E-5);
            vector<double> Solution2 = matrix_1.Solve_by_Gauss_Method_v(vb1);
            cout<<"Solution 2:"<<endl;
            View (Solution2);
            BOOST_CHECK_CLOSE_FRACTION(Solution2[0], -148.73333, 1E-5);
            BOOST_CHECK_CLOSE_FRACTION(Solution2[1], 74.73333, 1E-5);
            BOOST_CHECK_CLOSE_FRACTION(Solution2[2], 2.26666, 1E-5);

            double delta = Norm (Minus(Solution2, Solucion1))/Norm(Solucion1);
            cout<<"δ delta = "<<delta;
            BOOST_CHECK_CLOSE_FRACTION(delta, 0, DBL_EPSILON);
            //matrix_1.read_from_file("matrix_1.matr_l");//temp_matrix.matr_l
            double p[] = {1, 2, 2, 1};
            Matrix_SLE C (p, 2, 2);
            p[0] = 1;
            p[1] = 2;
            Matrix<double> c (p, 2, 1);
            temp_matrix = C.Solve(c);
            temp_matrix.view();
            BOOST_CHECK_CLOSE_FRACTION(temp_matrix[0], 1, DBL_EPSILON);
            BOOST_CHECK_CLOSE_FRACTION(temp_matrix[1], 0, DBL_EPSILON);
        }
BOOST_AUTO_TEST_SUITE_END()

void View (const vector<double>& El1)
{
    for (vector<double>::const_iterator i = El1.begin(); i < El1.end(); ++i)
        cout<<*i<<endl;
    return;
}
/*void* operator new (size_t size)
{
    Return the name of the object;
}*/
void randoe(int*, const int, const int, int=10);
/*int main()
{
tm a11;
long int i1 = 0;
time(&i1);

    string path("/home/alexander/Documents/Projects/Temp/1");
    //getline (cin, path);
    int i = path.find_last_of('/');
    string Name;
    Name.insert(0, path, i + 1, path.length() - i);
    path = "xterm -T " + Name + " -e /usr/bin/cb_console_runner LD_LIBRARY_PATH=$LD_LIBRARY_PATH:. " + path;
    //system ("gnome-terminal");
    //xterm -T 1 -e /usr/bin/cb_console_runner LD_LIBRARY_PATH=$LD_LIBRARY_PATH:. /home/alexander/Documents/Projects/Temp/1
    //system("/home/alexander/Desktop/apt-upd > qwqwqwqwwqwqwqwqww");
    //system ("rm qwqwqwqwwqwqwqwqww");
    system (path.c_str());
    //system ("read Dummy");
    return 0;
}*/