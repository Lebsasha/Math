// *** ADDED BY HEADER FIXUP ***
#include <istream>
// *** END ***
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

using namespace std;
/*double** Set (const int);
void View (double** A, const int n, int m = 0);
double** set_rand1(const int n, int max_e);
double** Multi (double** M1, double** M2, const int n);
double** Summ(double** M1, double** M2, const int n);
double* set_rand2(const int n, int max_e);
double* Multi (double* M1, double* M2, const int n);
double* Summ(double* M1, double* M2, const int n);*/
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
    Matrix& Plus_or_Minus (const Matrix& B, bool a)
    {
        try
        {
            if ( N != B.N || M != B.M )
            {
                throw 1;
            }
            T* pEnd = pa + N*M;
            if (!a)
                for (T* pA = pa-1, *pB = B.pa-1; pA < pEnd - 1;)
                {
                    *(++pA) += *(++pB);
                }
            else
                for (T* pA = pa-1, *pB = B.pa-1; pA < pEnd - 1;)
                {
                    *(++pA) -= *(++pB);
                }
            return *this;
        }
        catch (...)
        {
            cout<<"Dimentionses "<<M<<'*'<<N<<"and "<<B.M<<'*'<<B.N<<" must equal";
            return *this;
        }
    }
    void Increase_Name_Count (void)//////////////////////////////////////////Integrate to constructor
    {
        for (int i = Count; i > 1; i /= 10)
        {
            Name += i%10;
        }
    }
protected:
    string Name;
    int N;
    int M;
    T* pa;
    static int Count;
public:
    Matrix (): Name("Matrix "), N(1), M(1), pa (new T [N*M])
    {
        Null();
        Name += ++Count + 48;
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
        Name += ++Count + 48;
    }
    Matrix (const int a, const int b): Name("Matrix "), N(a), M(b), pa(new T [N*M])
    {
        Null();
        Name += ++Count + 48;
    }
    Matrix (const Matrix& A): Name("Matrix "), N(A.N), M(A.M), pa(0)
    {
        *this = A;
        Name += ++Count + 48;
    }
    Matrix (vector<T> va): Name("Matrix "), N(va.size()), M(1), pa(new T [N])
    {
        typename vector<T>::reverse_iterator iva = va.rbegin();
        for (double* pb = pa + N - 1; pb >= pa; --pb, ++iva)
        {
            *pb = *iva;
        }
        Name += ++Count + 48;
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
                    for (; j < b; ++j)
                        *(++pt) = 0;
            }
        }
        //pa = const_cast<int*>(pB);
        M = b;
        delete[] pa;
        pa = pT;
        return 1;
    }
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
                for (; pt <= pEnd_T; ++pt)
                    *pt = 0;
        }
        N = a;
        delete[] pa;
        pa = pT;
        return 1;
    }
    T* Get_pointer (void) const
    {
        return pa;
    }
    int Get_N (void) const
    {
        return N;
    }
    int Get_M (void) const
    {
        return M;
    }
    int Get_Size (void) const
    {
        return N*M;
    }
    T* Get_pa(void) const
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
        if (pa)
            delete[] pa;
        N = A.N;
        M = A.M;
        pa = new T [N*M];
        T* pT = A.pa + N*M;
        for (T* pt = pa + N*M-1; pt >= pa; --pt)
            *pt = *(--pT);
        return *this;
    }
    Matrix& Multi (Matrix& B)
    {
        try
        {
            if (M != B.N)
            {
                throw 1;
            }
            T* pc = new T [N*B.M];
            if (!pc)
                throw 2;
            for (T* t = pc+N*B.M-1; t >= pc; --t)
                *t = 0;
            T* pEnd = pa + N*M;
            T* pC = pc;
            for (T* pA = pa; pA < pEnd; pA += M)
                for (T* pB = B.pa; pB < B.pa + B.M; ++pB, ++pC)
                    for (T* pA1 = pA, *pB1 = pB; pA1 < pA + M; pB1 += B.M, ++pA1)
                        *pC += *pA1**pB1;
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
    {
        Matrix<T> a(B);
        return Multi (a);
    }
    Matrix operator* (Matrix& B)
    {
        return this->Multi(B);
    }
    Matrix operator* (const vector<T>& B)
    {
        return this->Multi(B);
    }
    vector<T> Multi_outp_vector (const vector<T>& B)
    {
//        Matrix<T> M ();
//        vector<T> a = ;
        return vector<T>(this->Multi(B));
    }
    Matrix<T> Multiple_Matrix_by_Number (const T P)const
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
        if (M != 1)
        {
            cout<<"Error while trynig to create vector from "<<Name<<" with "<<N<<'*'<<M<<"dimensionses";
            return vector<T>(0);
        }
        vector<T> va(0);
        T* pEnd = pa + N;
        for (T* pb = pa; pb < pEnd; ++pb)
        {
            va.push_back(*pb);
        }
        return va;
    }
    virtual T& operator[] (const int i)
    {
        try
        {
            if (i > N*M || i < 0)
            {
                throw i;

            }
            return *(pa+i);
        }
        catch (int ind)
        {
            cout<<"You have attemted to use "<<ind<<" element in "<<Name<<". Returning first element of this matrix."<<endl;
            return *pa;
        }
    }
    virtual T Const_Get_El (const int i) const
    {
        try
        {
            if (i > N*M || i < 0)
            {
                throw i;

            }
            return *(pa+i);
        }
        catch (int ind)
        {
            cout<<"You have attemted to use "<<ind<<" element in "<<Name<<". Returning first element of this matrix."<<endl;
            return *pa;
        }
    }
    void Rand(void)
    {
        for (T* pt = pa + N*M -1; pt >= pa; --pt)
            *pt = rand()%10;
        return;
    }
    virtual const Matrix<T>& View (void) const
    {
        int Nstr = N;
        const int Nstb = M;
        T* pA = pa;
        //cout<<Name<<endl;///////////////////////////////////////////////////
        while (Nstr--)
        {
            for (int j = 0; j < Nstb-1; j++)
            {
                cout<<setw(3)<<(*pA)<<" ";
                pA++;
            }
            cout<<setw(3)<<*pA++;
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
    double Nord_by_Euclid (void) const
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







double Norm (Matrix<double> A)
{
    double* pEnd = A.Get_pointer() + A.Get_N()*A.Get_M();
    double a = *A.Get_pointer();
    for (double* pb = A.Get_pointer() + 1; pb < pEnd; ++pb)
    {
        if (fabs(*pb) > fabs(a))
            a = *pb;
    }
    return a;
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





//class Function2
//{
//    double (*Fn)(double, double);
//    const static int Num_of_variables;
//public:
//    Function2 (double (*Fn1) (double, double)): Fn (Fn1)
//    {}
//    Function2 (void): Fn (NullF2)
//    {}
////    double Derivative_by_user_def (double (*Der) (double, double), double x,)
////    {
////        return Der(x);
////    }
//    int Get_Num_of_variables (void) const
//    {
//        return Num_of_variables;
//    }
//    Matrix_SLE Derivative_by_definition_M_in_column (const double x, const double y) const
//    {
//        Matrix<double> A (Num_of_variables, 1);
//        const double Epsx = 0.1;
//        const double Epsy = 0.1;
//        double* pa = A.Get_pa() + Num_of_variables - 1;
//        double M = 0.01;
//        *pa-- = (Fn(x + M*x, y) - Fn(x, y))/Epsx;
//        *pa = (Fn(x, y + M*y) - Fn(x, y))/Epsy;
//        return A;
//    }
//    double f(const double x, const double y) const
//    {
//        return Fn(x, y);
//    }
//    double Integral_by_user_def (double (*AntiDer) (double), const double from_x, const double from_y, const double to_x, const double to_y) const;
//    double Integral_by_definition (const double from_x, const double from_y, const double to_x, const double to_y) const;
//    ~Function2 (void)
//    {
//        Fn = NullF2;
//    }
//};
//const int Function2::Num_of_variables = 2;






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
