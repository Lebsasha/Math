#ifndef DIFFERENT
#define DIFFERENT
#ifdef VISUAL_STUDIO
#include "stdafx.h"
#endif // VISUAL_ST
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdio>
#include <vector>
#include <random>
std::mt19937 gen(2);//time(nullptr));
template<typename T>
typename std::enable_if<std::is_integral<T>::value, void>::type fill_random_values (T* pArray, const size_t Nstr, const size_t Nstb, const T& lowest= 0, const T& biggest= std::numeric_limits<T>::max())
{
    std::uniform_int_distribution<T> dist(lowest, biggest);
    size_t i = 0;
    while (i++ < Nstb*Nstr)
    {
        *(pArray++) = dist(gen);
    }
}

template<typename T>
typename std::enable_if<std::is_floating_point<T>::value, void>::type fill_random_values (T* pArray, const size_t Nstr, const size_t Nstb, const T& lowest=0, const T& biggest=std::numeric_limits<T>::max())
{
    std::uniform_real_distribution<T> dist(lowest, biggest);
    size_t i = 0;
    while (i++ < Nstb*Nstr)
    {
        *(pArray++) = dist(gen);
    }
}
/**
 * @result Prints in matrix Nstr × Nstb form numbers.
 */
template<typename T>
void View (T* pArray, size_t Nstr, const size_t Nstb) noexcept
{
    while (Nstr--)
    {
        for (int j = 0; j < Nstb-1; ++j)
        {
            std::cout<<std::setw(3)<<(*pArray)<<" ";
            pArray++;
        }
        std::cout<<std::setw(3)<<*pArray++;
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
}
/**
 * @result Print vector elements in row like "1 2 3 4".
 */
template<typename T>
void View (const std::vector<T>& El1) noexcept
{
    for (auto i = El1.cbegin(); i < El1.cend() - 1; ++i)
        std::cout<<*i<<' ';
    std::cout<<*(El1.cend()-1)<<std::endl;
}

template<class T, size_t N>
T& Find_min_by_abs (T(&arr)[N]) noexcept
{
    T& a = *arr;
    for (T* p = arr + N - 1; p > arr; --p)
    {
        if (fabs(a) > fabs(*p))
            a = *p;
    }
    return a;
}
template<class T, size_t N>
T& Find_max_by_abs (T(&arr)[N]) noexcept
{
    T& a = *arr;
    for (T* p = arr + N - 1; p > arr; --p)
    {
        if (fabs(a) < fabs(*p))
        {
            a = *p;
        }
    }
    return a;
}
constexpr double Celsius_to_Farenheigts (const double C)
{
    return (1.8*C + 32);
}
constexpr double Farenheigts_to_Celsius (const double F)
{
    return (F-32)*5/9;
}
constexpr bool if_simple (const int n)
{
    for (int i = 2; i <= n - 2; ++i)
    {
        if (!(n%i))
        {
            return false;
        }
    }
    return true;
}
///@returns a^b
double Pow (double a, int b) noexcept;
long long my_rand_number();
inline void Get_Pause() noexcept
{
#ifdef linux
    std::cout<<"Enter something and hit \"enter\" to continue"<<std::endl;
    system ("read dummy");
#else
    #ifdef windows
    std::cout<<"Hit enter to continue"<<std::endl;
    system("pause");
#else
    cerr<<"Unsupported OS for getting pause"<<endl;
#endif // windows
#endif // linux
}
#endif
