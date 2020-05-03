#include <iostream>
#include <math.h>
//#define nullptr NULL
#include "../Lib/Math.h"
using namespace std;

inline double f1(const double x)
{
    return 3*x/sqrt(1+x*x*x);
}
inline double f2 (const Matrix<double>& X)// 32
{
    return 4 - *X.Get_pointer()**X.Get_pointer() - *(X.Get_pointer()+1)**(X.Get_pointer()+1);
}

int main()
{
    cout<<sizeof(unsigned long)<<endl;
//    int VERYIMPORTANT = 1;
    Function F (f1);
    cout<<"Trapetia "<<F.Integral_by_trapetia_method_with_epsilon(0, 1.075, 10e-5)<< "   "<<"Sympthon "<<F.Integral_by_Sympthon_method_with_epsilon(0, 1.075, 10e-5)<<endl;
    Function2 F2 (f2);
    Matrix<double> From (2, 1);
    Matrix<double> To (2, 1);
    From[0] = From[1] = -1;
    To[0] = To[1] = 1;
    cout<<F2.Integral_by_definition(From, To);
//    register int i = 0;
//    cout<<&i<<endl;//1.44982 37.33
    return 0;
}


//The compiler generates a warning if you take the address of a variable with register storage class.
//Example
//void foo(void)
//{
//    register int i;
//    int *j = &i;
//}


//For example, you might have assembly code and C code that uses the same symbol name, such as counter. Therefore, you can export a different name to be used by the assembler:
//int counter __asm__("counter_v1") = 0;
//This exports the symbol counter_v1 and not the symbol counter.


//An empty declaration, that is a semicolon with nothing before it, is permitted.
//
//Example
//; // do nothing


//#include <limits.h>
//numeric_limits<double>::infinity() = inf
//98 Standart 199711
//11 Standart 201103
//// Default placement versions of operator new.
//inline void* operator new(std::size_t, void* __p) _GLIBCXX_USE_NOEXCEPT
//{ return __p; }
//inline void* operator new[](std::size_t, void* __p) _GLIBCXX_USE_NOEXCEPT
//{ return __p; }
//
//// Default placement versions of operator delete.
//inline void operator delete  (void*, void*) _GLIBCXX_USE_NOEXCEPT { }
//inline void operator delete[](void*, void*) _GLIBCXX_USE_NOEXCEPT { }
////@}
//#if __cplusplus >= 201103L
//    /** The number of base 10 digits required to ensure that values which
//	differ are always differentiated.  */
//    static constexpr int max_digits10 = 0;
//#endif
