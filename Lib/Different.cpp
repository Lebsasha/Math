// *** ADDED BY HEADER FIXUP ***
#include <ctime>
#include <istream>
// *** END ***
#include <time.h>
#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <assert.h>
using namespace std;
void randoe (int* pArray, const int Nstr, const int Nstb, int biggest_1) /*randoe (&A[0][0], n, n, 50);*/
{
    int i = 0;
    srand(static_cast<int>(time(nullptr)));
    while (i++ < Nstb*Nstr)
    {
        *(pArray++) = rand()%biggest_1;
    }
    return;
}
void randoe (double* pArray, const int Nstr, const int Nstb, int biggest_1) /*randoe (&A[0][0], n, n, 50);*/
{
    int i = 0;
    srand(static_cast<int>(time(nullptr)));
    while (i++ < Nstb*Nstr)
    {
        *(pArray) = rand()%biggest_1+rand()%100/100.0;
    }
    return;
}
void randoe (float* pArray, const int Nstr, const int Nstb, int biggest_1) /*randoe (&A[0][0], n, n, 50);*/
{
    int i = 0;
    srand(static_cast<int>(time(nullptr)));
    while (i++ < Nstb*Nstr)
    {
        *(pArray++) = rand()%biggest_1+rand()%100/100.0;
    }
    return;
}
void View (int* pArray, int Nstr, int Nstb)
{
    while (Nstr--)
    {
        for (int j = 0; j < Nstb-1; j++)
        {
            cout<<setw(3)<<(*pArray)<<" ";
            pArray++;
        }
        cout<<setw(3)<<*pArray++;
        cout<<endl;
    }
    cout<<endl;
    return;
}
void View (double* pArray, int Nstr, const int Nstb)
{
    while (Nstr--)
    {
        for (int j = 0; j < Nstb-1; j++)
        {
            cout<<setw(3)<<(*pArray)<<" ";
            pArray++;
        }
        cout<<setw(3)<<*pArray++;
        cout<<endl;
    }
    cout<<endl;
    return;
}
void View (float* pArray, int Nstr, int Nstb)
{
    while (Nstr--)
    {
        for (int j = 0; j < Nstb-1; j++)
        {
            cout<<setw(3)<<(*pArray)<<" ";
            pArray++;
        }
        cout<<setw(3)<<*pArray++;
        cout<<endl;
    }
    cout<<endl;
    return;
}
void View (double** A, const int n, int m)
{
    if (!m)
    {
        m = n;
    }
    for (int i = 0; i < n; ++i)
    {
        View (*(A++), 1, m);
    }
    return;
}
void View (float** A, const int n, int m)
{
    if (!m)
    {
        m = n;
    }
    for (int i = 0; i < n; ++i)
    {
        View (*(A++), 1, m);
    }
    return;
}
void View (int** A, const int n, int m)
{
    if (!m)
    {
        m = n;
    }
    for (int i = 0; i < n; ++i)
    {
        View (*(A++), 1, m);
    }
    return;
}
#include <vector>
#include <iterator>
void View (vector<double>& El1)
{
    for (vector<double>::const_iterator i = El1.begin(); i < El1.end() - 1; ++i)
        cout<<*i<<' ';
    cout<<*(El1.end()-1);
    return;
}
void Nul(double* p, const int n)
{
    double* pb = p;
    double* pe = p + n;
    while (pb != pe)
    {
        *(pb++) = 0;
    }
    return;
}
void Nul(int* p, const int n)
{
    int* pb = p;
    int* pe = p + n;
    while (pb != pe)
    {
        *(pb++) = 0;
    }
    return;
}
void Nul(float* p, const int n)
{
    float* pb = p;
    float* pe = p + n;
    while (pb != pe)
    {
        *(pb++) = 0;
    }
    return;
}
void Nul(char* p, const int n)
{
    char* pb = p;
    char* pe = p + n;
    while (pb != pe)
    {
        *(pb++) = 0;
    }
    return;
}
inline double Celsius_to_Farenheigts (const double C)
{
    return (1.8*C + 32);
}
inline double Farenheigts_to_Celsius (const double F)
{
    return (F-32)*5/9;
}
bool if_Simple (const int n)
{
    for (int i = 2; i <= n - 2; ++i)
    {
        if (!(n%i))
        {
            return 0;
        }
    }
    return 1;
}
double Pow (const double x, const int i, const double Prdct= 1)
{
    assert(i >= 0);
    if (i == 1)
        return Prdct*x;
    if (i == 0)
        return Prdct;
    return Pow(x, i - 1, Prdct*x);
}
void Get_Pause(void)
{
#ifdef linux
    system ("read dummy");
#endif // linux
#ifdef windows
    system("pause");
#endif // windows
    return;
}
long long Myrandnumber(void)
{
    int* pointer = new int (5);
    delete pointer;
    return (long long) pointer;
}
