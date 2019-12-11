#ifndef DIFFERENT_H
#define DIFFERENT_H
#ifdef VISUAL_STUDIO
#include "stdafx.h"
#endif // VISUAL_ST
void randoe (int*, const int, const int, int);
void randoe (double*, const int, const int, int);
void randoe (float*, const int, const int, int);
void View (int*, const int, const int);
void View (double*, const int, const int);
void View (float*, const int, const int);
void View (double**, const int, int m = 0);
void View (float**, const int, int m = 0);
void View (int**, const int, int m = 0);
#include <vector>
//using namespace std;
void View (std::vector<double>& El1);
void Nul (double*, const int);
void Nul (int*, const int);
void Nul (float*, const int);
void Nul (char*, const int);
inline double Celsius_to_Farenheigts (double a);
inline double Farenheigts_to_Celsius (double a);
bool if_Simple (int);
double Pow (const double, const int, const double = 1);
void Get_Pause (void);
long long Myrandnumber(void);
template<class T>
T Find_min_by_abs (T* pa, const int s)
{
    T a = *pa;
    for (T* p = pa + s - 1; p > pa; --p)
    {
        if (fabs(a) > fabs(*p))
            a = *p;
    }
    return a;
}
template<class T>
int Find_max_by_abs (T* pa, const int s)
{
    int i = 0;
    int indx = s - 1;
    T a = *pa;
    for (T* p = pa + indx; p > pa; --p, --indx)
    {
        if (fabs(a) < fabs(*p))
        {
            a = *p;
            i = indx;
        }
    }
    return i;
}
#endif
