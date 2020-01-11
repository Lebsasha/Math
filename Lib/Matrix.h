// *** ADDED BY HEADER FIXUP ***
#include <istream>
// *** END ***
#ifdef VISUAL_STUDIO
#include "stdafx.h"
#endif // VISUAL_ST
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include <vector>
#include <iterator>
#include <math.h>
#include <assert.h>
#include <algorithm>
#if __cplusplus < 201103L
#define nullptr NULL
#endif

using namespace std;
template <class T>
class Vector
{
    T* a;
    int N;
public:
    Vector ()
    {
        N = 1;
        a = new T [N];
        *a = 0;
    }
    Vector (const T A)
    {
        N = 1;
        a = new T [N];
        *a = A;
    }
    /*Vector (const T* pA, const int N)
    {
    T* pa = pA + N -1;
    a = new T [N];
    for (--a; pA <= pa; ++pA)
    {
    *(++a) = *pA;
    }
    }*/
    Vector (Vector& A)
    {
        N = A.N;
        a = new T [N];
        --a;
        T* End = A.a + A.N-1;
        for (T* pA = A.a; pA <= End; ++pA)
        {
            *(++a) = *pA;
        }
    }
    void Rand (int Biggest_Num = 5)
    {
        ++Biggest_Num;
        for (T* pA = a + N-1; pA >= a; --pA)
            *pA = rand()%Biggest_Num;
    }
    friend ostream& operator<< (ostream& a, const Vector& v)
    {
        T* End = v.a + v.N - 1;
        for (T* pa = v.a; pa <= End; ++pa)
            a<<*pa;
        return a;
    }
    inline void Set_size (int Nt)
    {
        this->N = Nt;
        return;
    }
    inline int Get_size (void) const
    {
        return N;
    }

    T& operator[] (const int i)
    {
        if (i < N)
            return *(a+i);
        else
            exit (1);
    }
    ~Vector (void)
    {
        delete[] a;
    }
};

template <class T>
class Matrix
{
    Matrix Plus_or_Minus (const Matrix& B, bool a)
    {
            //this->View();
            //B.View();
            assert (N == B.N && M == B.M );
            T* pc = new T[N*M];
            T* pEnd = B.pa + B.N*B.M - 1;
            if (!a)
                for (T* pA = pa-1, *pB = B.pa-1, *pC = pc - 1; pB < pEnd;)
                {
                    *(++pC) = *(++pA) + *(++pB);
                }
            else
                for (T* pA = pa-1, *pB = B.pa-1, *pC = pc - 1; pB < pEnd;)
                {
                    *(++pC) = *(++pA) - *(++pB);
                }
            Matrix A (pc, N, M);
            delete[] pc;
            return A;
    }
    void Increase_Name_Count (void)
    {
        string a;
        for (int i = ++Count; i > 1; i /= 10)
        {
            a += '0' + i%10;
        }
        reverse (a.begin(), a.end());
        Name += a;
        //Name = to_string(Count);
    }
protected:
    string Name;
    int N;
    int M;
    T* pa;
    static int Count;
public:
    Matrix (bool If_null = 1): Name("Matrix "), N(1), M(1), pa (new T [N*M])
    {
        if (If_null)
        Null();
        Increase_Name_Count();
    }
    Matrix (T* pA, const int a, const int b): Name("Matrix "), N(a), M(b), pa(new T [N*M])
    {
        if (pA)
        {
            T* pt = pa+N*M;
            for (T* pB = pA + N*M - 1; pB >= pA; --pB)
            {
                *(--pt) = *pB;
            }
        }
        Increase_Name_Count();
    }
<<<<<<< HEAD
    Matrix (const int a, const int b, bool If_null = 1): Name("Matrix "), N(a >= 0 ? a : 0), M(b >= 0 ? b : 0), pa(new T [N*M])
=======
    Matrix (const int a, const int b): Name("Matrix "), N(a), M(b), pa(new T [N*M])
>>>>>>> ec95f91ed173f170c4c98998802d56dbf22b6d2e
    {
        if (If_null)
        Null();
        Increase_Name_Count();
    }
    Matrix (const Matrix& A): Name("Matrix "), N(A.N), M(A.M), pa(nullptr)
    {
        *this = A;
        Increase_Name_Count();
    }
    Matrix (vector<T> va): Name("Matrix "), N(va.size()), M(1), pa(new T [N])
    {
        typename vector<T>::reverse_iterator iva = va.rbegin();
        for (T* pb = pa + N - 1; pb >= pa; --pb, ++iva)
        {
            *pb = *iva;
        }
        Increase_Name_Count();
    }
    bool Edit_col (const int b)
    {
        if (!pa)
        {
            pa = new T [(N=1)*(M=b)];
            return pa ? 1 : 0;
        }
        T* pB = pa;
        T* pT = new T[N*b];
        if (!pT)
            return 0;
        T* pEnd = pT + N*(b-1)-1;
        --pB;
        for (T* pt = pT-1; pt <= pEnd; pB += (b >= M) ? 0 : M-b)
        {
            for (int j = 0; j < b; ++j)
            {
                if (j < M)
                    *(++pt) = *(++pB);
                else
                {
                    for (; j < b; ++j)
                        *(++pt) = 0;
                #ifdef NDEBUG
            goto AfterNullN;
#endif // NDEBUG
                }
            }
        }
        AfterNullN:
        M = b;
        delete[] pa;
        pa = pT;
        return 1;
    }
<<<<<<< HEAD
    bool If_Symmetric (void) const
    {
        if (N != M)
            return 0;
        int i = 0;
        T* pcEnd = this->First_i();
        pcEnd += M - 1;
        for (T* pb = this->Last_i(), *pc = pb; pb >= pa; pb -= i + 1, pc += M*(i-1) - 1 - 1)
        {
            for (; pc >= pcEnd; --pb, pc -= M)
            {
                if (*pb != *pc)
                {
                    return 0;
                }
            }
            ++i;
        }
        return 1;
    }
=======
>>>>>>> ec95f91ed173f170c4c98998802d56dbf22b6d2e
    bool Edit_row (const int a)
    {
        if (!pa)
        {
            pa = new T [(N=a)*(M=1)];
            return pa ? 1 : 0;
        }
        T* pT = new T [a*M];
        if (!pT)
            return 0;
        T* pEnd_T = pT + a*M - 1;
        T* pEnd_B = pa + N*M -2;
        T* pB = pa-1;
        for (T* pt = pT; pt <= pEnd_T; ++pt)
        {
            if (pB <= pEnd_B)
                *pt = *(++pB);
            else
            {
                for (; pt <= pEnd_T; ++pt)
                    *pt = 0;
#ifdef NDEBUG
            goto AfterNullR;
#endif // NDEBUG
            }
        }
AfterNullR:
        N = a;
        delete[] pa;
        pa = pT;
        return 1;
    }
    inline T* Get_pointer (void) const
    {
        return pa;
    }
    inline int Get_N (void) const
    {
        return N;
    }
    inline int Get_M (void) const
    {
        return M;
    }
    inline int Get_Size (void) const
    {
        return N*M;
    }
    inline T* Get_pa(void) const
    {
        return pa;
    }
    Matrix operator+ (const Matrix& B)
    {
        return Plus_or_Minus(B, 0);
    }
    Matrix operator- (const Matrix& B)
    {
        return Plus_or_Minus(B, 1);
    }
    Matrix operator- (const vector<T>& B)
    {
        Matrix A (B);
        return Plus_or_Minus(A, 1);
    }
    Matrix operator+ (const vector<T>& B)
    {
        Matrix A (B);
        return Plus_or_Minus(A, 1);
    }
    const Matrix& operator= (const Matrix& A)
    {
        delete[] pa;
        Name = A.Name;
        N = A.N;
        M = A.M;
        pa = new T [N*M];
        T* pT = A.pa + N*M;
        for (T* pt = pa + N*M-1; pt >= pa; --pt)
            *pt = *(--pT);
        return *this;
    }
<<<<<<< HEAD
    Matrix Multi (Matrix& B) const
=======
    Matrix& Multi (Matrix& B)
>>>>>>> ec95f91ed173f170c4c98998802d56dbf22b6d2e
    {
            assert(M == B.N);
            T* pc = new T [N*B.M];
            if (!pc)
            {
                cout<<"Haven't enough dynamic memory."<<endl;
                assert(pc);
            }
            for (T* t = pc+N*B.M-1; t >= pc; --t)
                *t = 0;
            T* pEnd = pa + N*M;
            T* pC = pc;
            for (T* pA = pa; pA < pEnd; pA += M)
                for (T* pB = B.pa; pB < B.pa + B.M; ++pB, ++pC)
                    for (T* pA1 = pA, *pB1 = pB; pA1 < pA + M; pB1 += B.M, ++pA1)
                        *pC += *pA1**pB1;
<<<<<<< HEAD
            Matrix<T> A (pc, N, B.M);
            delete[] pc;
            return A;
    }
    Matrix Multi (const vector<T>& B) const
=======
            delete[] pa;
            pa = pc;
            M = B.M;
            return *this;
        }
        catch (int i)
        {
            if (i == 1)
                cout<<"Dimentionses "<<M<<"and "<<B.N<<" must equal";
            if (i == 2)
                cout<<"Haven't enough dynamic memory."<<endl;
            return *this;
        }
    }
    Matrix& Multi (const vector<T>& B)
>>>>>>> ec95f91ed173f170c4c98998802d56dbf22b6d2e
    {
        Matrix<T> a(B);
        return Multi (a);
    }
    Matrix operator* (Matrix& B) const
    {
        return this->Multi(B);
    }
    Matrix operator* (const vector<T>& B) const
    {
        return this->Multi(B);
    }
    vector<T> Multi_outp_vector (const vector<T>& B) const
    {
//        Matrix<T> M ();
//        vector<T> a = ;
        return vector<T>(this->Multi(B));
    }
    Matrix<T> Multiple_Matrix_by_Number (const T P) const
    {
        Matrix<T> A (N, M);
        for (T* pc = A.pa + A.N*A.M - 1, *pb = pa + N*M - 1; pc >= A.pa; --pc, --pb)
        {
            *pc = *pb*P;
        }
        return A;
    }
    operator vector<T> ()
    {
        vector<T> va(0);
        if (M != 1)
        {
            cout<<"Error while trynig to create vector from "<<Name<<" with "<<N<<'*'<<M<<"dimensionses";
            return va;
        }
        T* pEnd = pa + N;
        for (T* pb = pa; pb < pEnd; ++pb)
        {
            va.push_back(*pb);
        }
        return va;
    }
    T Max_element(void) const
    {
        assert (pa != nullptr);
        T El = *(pa + N*M - 1);
        for (T* pb = pa + N*M-2; pb >= pa; --pb)
        {
            if (El < *pb)
                El = *pb;
        }
        return El;
    }
    T Max_element_fabs(void) const
    {
        assert (pa != nullptr);
        T El = fabs(*(pa + N*M - 1));
        for (T* pb = pa + N*M-2; pb >= pa; --pb)
        {
            if (El < fabs(*pb))
                El = fabs(*pb);
        }
        return El;
    }
    T Min_element(void) const
    {
        assert (pa != nullptr);
        T El = *(pa + N*M - 1);
        for (T* pb = pa + N*M-2; pb >= pa; --pb)
        {
            if (El > *pb)
                El = *pb;
        }
        return El;
    }
    T Min_element_fabs(void) const
    {
        assert (pa != nullptr);
        T El = fabs(*(pa + N*M - 1));
        for (T* pb = pa + N*M-2; pb >= pa; --pb)
        {
            if (El > fabs(*pb))
                El = fabs(*pb);
        }
        return El;
    }
    inline virtual T& Unsafe_index (const int i)
    {
        return *(pa+i);
    }
    inline virtual T Unsafe_index_c (const int i) const
    {
        return *(pa+i);
    }
    virtual T& operator[] (const int i)
    {
        assert(i < N*M && i >= 0);
//            if (i == 1)
//                cout<<"HERE"<<endl;
            return *(pa+i);
    }
    virtual T operator[] (const int i) const
    {
        assert(i < N*M && i >= 0);
            //if (i == 1)
            //cout<<"HERE"<<endl;
            return *(pa+i);
    }
#ifndef NDEBUG
    class iterator
    {
        T* p;
        T* const pa;
        int Size;
    public:
        Matrix<T>::iterator& operator= (const Matrix<T>::iterator& a)// !!! PRIVATE !!!
        {
            assert (pa == a.pa && Size == a.Size);
            p = a.p;
            return *this;
        }
    public:
        iterator (T* a, T* const pb, const int s): p(a), pa(pb), Size(s)
        {
            assert (a < pb + s);
        }
        iterator (void): p(nullptr), pa(nullptr), Size(0)
        {}
        inline Matrix<T>::iterator& operator++ (void)
        {
            assert (p <= pa + Size - 1);
#ifndef NDEBUG
            if (p == pa + Size - 1)
                cout<<"Warning! p == pa  + Size - 1 in prefix ++"<<endl;
#endif // NDEBUG
            ++p;
            return *this;
        }
        inline Matrix<T>::iterator operator++ (int)
        {
            assert (p <= pa + Size - 1);
#ifndef NDEBUG
            if (p == pa + Size - 1)
                cout<<"Warning! p == pa  + Size - 1 in postfix ++"<<endl;
#endif // NDEBUG
            return Matrix<T>::iterator(p++, pa, Size);
        }
        inline Matrix<T>::iterator& operator-- (void)
        {
            assert (p >= pa);
#ifndef NDEBUG
            if (p == pa)
                cout<<"Warning! p == pa in prefix --"<<endl;
#endif // NDEBUG
            --p;
            return *this;
        }
        inline Matrix<T>::iterator operator-- (int)
        {
            assert (p >= pa);
#ifndef NDEBUG
            if (p == pa)
                cout<<"Warning! p == pa in postfix --"<<endl;
#endif // NDEBUG
            return Matrix<T>::iterator(p--, pa, Size);
        }
        inline Matrix<T>::iterator operator-= (int i)
        {
            assert (p >= pa);
#ifndef NDEBUG
            if (p - i < pa)
                cout<<"Warning! p < pa in -="<<i<<endl;
#endif // NDEBUG
            p -= i;
            return *this;
        }
        inline Matrix<T>::iterator operator+= (int i)
        {
            assert (p < pa + Size);
#ifndef NDEBUG
            if (p + i >= pa + Size)
                cout<<"Warning! p + i >= pa + Size in +="<<i<<endl;
#endif // NDEBUG
            p += i;
            return *this;
        }
        inline Matrix<T>::iterator operator+ (const int i) const
        {
            assert (p < pa + Size);
#ifndef NDEBUG
            if (p + i >= pa + Size)
                cout<<"Warning! p + i >= pa + Size in + "<<i<<endl;
#endif // NDEBUG
            Matrix<T>::iterator a = *this;
            a.p += i;
            return a;
        }
        inline Matrix<T>::iterator operator- (const int i) const
        {
            assert (p >= pa);
#ifndef NDEBUG
            if (p - i < pa)
                cout<<"Warning! p < pa in - "<<i<<endl;
#endif // NDEBUG
            Matrix<T>::iterator a = *this;
            a.p -= i;
            return a;
        }
        inline T& operator* (void)
        {
            assert (p != nullptr);
            return *p;
        }
        inline operator T* ()
        {
            return p;
        }
        inline bool operator== (const Matrix<T>::iterator& a) const
        {
            assert (pa == a.pa);
            return this->p == a;
        }
        inline bool operator<= (const Matrix<T>::iterator& a) const
        {
            assert (pa == a.pa);
            return p <= a.p;
        }
        inline bool operator>= (const Matrix<T>::iterator& a) const
        {
            assert (pa == a.pa);
            return p >= a.p;
        }
        inline bool operator> (const Matrix<T>::iterator& a) const
        {
            assert (pa == a.pa);
            return p > a.p;
        }
        inline bool operator< (const Matrix<T>::iterator& a) const
        {
            assert (pa == a.pa);
            return p < a.p;
        }
    };
//#else
//    using Matrix<T>::iterator T*
#endif // NDEBUG
#ifndef NDEBUG
    Matrix::iterator First_i (void) const
    {
        return Matrix::iterator(pa, pa, N*M);
    }
    Matrix::iterator Last_i (void) const
    {
        return Matrix::iterator(pa + Get_Size() - 1, pa, N*M);
    }
#else // NDEBUG
    inline T* First_i (void) const
    {
        return pa;
    }
    inline T* Last_i (void) const
    {
        return pa + Get_Size() - 1;
    }
#endif // NDEBUG
    virtual T Const_Get_El (const int i) const
    {
            assert(i < N*M && i >= 0);
            return *(pa+i);
    }
    void Rand(void)
    {
        for (T* pt = pa + N*M -1; pt >= pa; --pt)
            *pt = rand()%10;
        return;
    }
    virtual const Matrix<T>& View (const int bytes_per_element = 8, const bool Show_Name = 0) const
    {
        int Nstr = N;
        const int Nstb = M;
        T* pA = pa;
        if (Show_Name)
            cout<<Name<<endl;
        while (Nstr--)
        {
            for (int j = 0; j < Nstb-1; j++)
            {
                cout<<setw(bytes_per_element)<<(*pA)<<" ";
                pA++;
            }
            cout<<setw(bytes_per_element)<<*pA++;
            cout<<endl;
        }
        return *this;
    }
    Matrix<T>& View_BAD (const int bytes_per_element = 8, const bool Show_Name = 0)
    {
        int Nstr = N;
        const int Nstb = M;
        T* pA = pa;
        if (Show_Name)
            cout<<Name<<endl;
        while (Nstr--)
        {
            for (int j = 0; j < Nstb-1; j++)
            {
                cout<<setw(bytes_per_element)<<(*pA)<<" ";
                pA++;
            }
            cout<<setw(bytes_per_element)<<*pA++;
            cout<<endl;
        }
        return *this;
    }
    virtual bool Read_from_file (string& sa)
    {
        ifstream ifsa (sa.c_str(), ios::binary);
        bool Temp = Read_from_file(ifsa);
        ifsa.close();
        return Temp;
    }
    virtual bool Read_from_file (const char* A)
    {
        ifstream ifsa (A, ios::binary);
        bool Temp = Read_from_file(ifsa);
        ifsa.close();
        return Temp;
    }
    virtual bool Read_from_file (ifstream& sa)
    {
        if (sa.eof())
            return 0;
        sa.read(reinterpret_cast<char*>(&N), sizeof(int));
        if (sa.eof())
            return 0;
        sa.read(reinterpret_cast<char*>(&M), sizeof(int));
        pa = new T [N*M];
        if (!pa)
            return 0;
        const T* pEnd = pa + N*M;
        const int Size = sizeof(T);
        for (T* pb = pa; pb < pEnd; ++pb)
        {
            if (sa.eof())
                return 0;
            sa.read(reinterpret_cast<char*>(pb), Size);
        }
        return 1;
    }
    virtual void Write_to_file (string& a) const
    {
        ofstream ofsa (a.c_str(), ios::binary);
         Write_to_file (ofsa);
        ofsa.close();
        return;
    }
    virtual void Write_to_file (const char* A) const
    {
        ofstream ofsa (A, ios::binary);
        Write_to_file (ofsa);
        ofsa.close();
        return;
    }
    virtual void Write_to_file (ofstream& sa) const
    {
        //string a = Name;
        //a.erase(0, 7);
        //int Num = a;
        //sa.write(reinterpret_cast<char*>(), sizeof(int));
        sa.write(reinterpret_cast<const char*>(&N), sizeof(int));
        sa.write(reinterpret_cast<const char*>(&M), sizeof(int));
        const int Size = sizeof(T);
        const T* pEnd = pa + N*M;
        for (T* a = pa; a < pEnd; ++a)
        {
            sa.write(reinterpret_cast<char*>(a), Size);
        }
        return;
    }
    void Null (void)
    {
        for (T* t = pa+N*M-1; t >= pa; --t)
            *t = 0;
        return;
    }
    double Norm_by_infinity (void) const
    {
        double* pEnd = pa + N*M;
        double a = *pa;
        for (double* pb = pa + 1; pb < pEnd; ++pb)
        {
            if (fabs(*pb) > fabs(a))
                a = *pb;
        }
        return a;
    }
    double Norm_by_Euclid (void) const
    {
        double a = 0;
        for (double* pb = pa + N*M - 1; pb >= pa; --pb)
        {
            a += *pb**pb;
        }
        return sqrt(a);
    }
    double Norm_by_Max (void) const
    {
        double a = 0;
        for (double* pb = pa + N*M - 1; pb >= pa; --pb)
        {
            a += fabs(*pb);
        }
        return a;
    }
    virtual ~Matrix(void)
    {
        delete[] pa;
    }
};
template <class T>
int Matrix<T>::Count = 0;





template <>
bool Matrix<double>::If_Symmetric (void) const
{
    if (N != M)
        return 0;
    int i = 1;
    double* pcEndCol = this->First_i();
    for (double* pb = this->Last_i() - 1, *pc = pb - (M - 1); pb >= pa; pb -= i, pc += (N-i)*M - 1)
    {
        for (; pc >= pcEndCol; --pb, pc -= M)
        {
            if (*pb - *pc > DBL_EPSILON)
            {
                return 0;
            }
        }
        ++i;
    }
    return 1;
}

template <>
bool Matrix<float>::If_Symmetric (void) const
{
    if (N != M)
        return 0;
    int i = 1;
    float* pcEndCol = this->First_i();
    for (float* pb = this->Last_i() - 1, *pc = pb - (M - 1); pb >= pa; pb -= i, pc += (N-i)*M - 1)
    {
        for (; pc >= pcEndCol; --pb, pc -= M)
        {
            if (*pb - *pc > DBL_EPSILON)
            {
                return 0;
            }
        }
        ++i;
    }
    return 1;
}


//class Array_of_Functions2
//{
//    int N;
//    Function2* Fn;
//public:
//    Array_of_Functions2 (void): N(1), Fn (new Function2[N])
//    {}
//    Array_of_Functions2 (const int S): N(S), Fn (new Function2[N])
//    {}
//    Array_of_Functions2 (Array_of_Functions2& Ar): N(Ar.N), Fn (new Function2[N])
//    {
//        for (Function2* pFn = Fn + N - 1, *pFa = Ar.Fn +  Ar.N - 1; pFn >= Fn; --pFn, --pFa)
//        {
//            *pFn = *pFa;
//        }
//    }
//    Function2* Get_Fn (void)
//    {
//        return Fn;
//    }
//    Function2 operator[] (const int S) const
//    {
//        if (S < N)
//            return NullF2;
//        return *(Fn + S);
//    }
//    Matrix<double> f (const Matrix<double> X) const
//    {
//        Matrix<double> A (N, 1);
//        double* pa = A.Get_pa() + N - 1;
//        for (Function2* pFn = Fn + N - 1; pFn >= Fn; --pFn, --pa)
//        {
//            *pa = pFn->f(X);
//        }
//        return A;
//    }
//    Matrix_SLE Derivative (const Matrix<double> X) const
//    {
//        int Num_of_variables = X.Get_Size();
//        Matrix<double> A (N, Num_of_variables);
//        double* pa = A.Get_pa() + A.Get_Size() - 1;
//        Matrix<double> Temp (Num_of_variables, 1);
//        double* pTemp = 0;
//        for (Function2* pFn = Fn + N*Num_of_variables - 1; pFn >= Fn; --pFn, --pa)
//        {
//            Temp = pFn->Derivative_by_definition_M_in_column(X);
//            pTemp = Temp.Get_pa() + Num_of_variables - 1;
//            for (int i = 0; i < Num_of_variables; ++i, --pa, --pTemp)
//            {
//                *pa = *pTemp;
//            }
//        }
//        return A;
//    }
//    Matrix<double> Integral (const double from_x, const double from_y, const double to_x, const double to_y) const;
//    double Max (const Matrix<double> X)
//    {
//        double m = 0;
//        double v = 0;
//        for (Function2* pA = Fn + N - 1; pA >= Fn; --pA)
//        {
//            if (m < (v = pA->f(X)))
//            m = v;
//        }
//            return m;
//    }
////    ~Array_of_Functions2 (void)
////    {
////        delete[] Fn;
////    }
//};











//using Fune = double (*)(double);
//using pFune = double (**) (double);
//class Matrix_of_Functions
//{
//    int N;
//    double (**Fn)(double);
//public:
//    Matrix_of_Functions (void): N(1), Fn (new Fune [N])
//    {}
//    Matrix_of_Functions (const int S): N(S), Fn (new Fune[N])
//    {}
//    pFune Get_Fn (void)
//    {
//        return Fn;
//    }
//    Fune operator[] (const int S)
//    {
//        if (S < N)
//        return *(Fn + S);
//        return NullF;
//    }
//    Matrix_SLE()
//    };
