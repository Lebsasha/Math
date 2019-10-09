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
    T* Get_pointer (void)
    {
        return pa;
    }
    int Get_N (void)
    {
        return N;
    }
    int Get_M (void)
    {
        return M;
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
    Matrix& operator* (Matrix& B)
    {
        return this->Multi(B);
    }
    Matrix& operator* (const vector<T>& B)
    {
        return this->Multi(B);
    }
    vector<T> Multi_outp_vector (const vector<T>& B)
    {
//        Matrix<T> M ();
//        vector<T> a = ;
        return vector<T>(this->Multi(B));
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
    void Rand(void)
    {
        for (T* pt = pa + N*M -1; pt >= pa; --pt)
            *pt = rand()%10;
        return;
    }
    virtual void View (void) const
    {
        int Nstr = N;
        const int Nstb = M;
        T* pA = pa;
        cout<<Name<<endl;
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
        return;
    }
    void Null (void)
    {
        for (T* t = pa+N*M-1; t >= pa; --t)
            *t = 0;
        return;
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
class Matrix_SE: public Matrix<double>
{
    bool High_accuracy;
protected:
    char Type;
public:
    Matrix_SE  (): Matrix<double>(), High_accuracy(1), Type('N')
    {}
    Matrix_SE (Matrix<double>& A): Matrix<double>(A), High_accuracy(1), Type('N')
    {}
    Matrix_SE (double* pA, const int a, const int b): Matrix<double>(pA, a, b),  High_accuracy(1), Type('N')
    {}
    Matrix_SE (Matrix_SE& A): Matrix<double>(*(reinterpret_cast<Matrix<double>*>(&A))), High_accuracy (A.High_accuracy), Type('N')
    {}
    void Set_High_accuracy (const bool a)
    {
        High_accuracy = a;
        return;
    }
    bool Get_High_accuracy (void) const
    {
        return High_accuracy;
    }
    void Multi_Row_by_Num (const int N_R, const double C)
    {
        double* pb = pa + N_R*M;
        double* pendb = pb + M;
        for (double* pt = pb; pt < pendb; ++pt)
        {
            *pt = C*(*pt);
        }
        return;
    }
    void One_Minus_Second_Row (const int one, const int second, const int C)
    {
        double* p1 = pa + one*M;
        const double* pe1 = p1 + M;
        double* p2 = pa + second*M;
        for (; p1 < pe1; ++p1, ++p2)
        {
            *p1 = *p1 - C**p2;
        }
        return;
    }
    void Replace_Rows (const int One, const int Second)
    {
        double* pe1 = pa + (One+1)*M;
        double* pe2 = pa + (Second+1)*M;
        if (!High_accuracy)
            for (double *p1 = pe1 - M, *p2 = pe2 - M; p1 < pe1; ++p1, ++p2)
            {
                *p1 = *p1 + *p2;
                *p2 = *p1 - *p2;
                *p1 = *p1 - *p2;
            }
        else
        {
            double Temp = 0;
            for (double *p1 = pe1 - M, *p2 = pe2 - M; p1 < pe1; ++p1, ++p2)
            {
                Temp = *p1;
                *p1 = *p2;
                *p2 = Temp;
            }
        }
        return;
    }
    int Find_Max_El_in_Col (const int Col) const
    {
        double* pb = pa + Col*(M+1);
        const double* pEnd = pa + N*M;
        //double* p2 = pa + Col + M;
        double El = *pb;
        int Ind_of_Max = Col;
        int Cur_Ind = Col;
        while ((pb = pb + M) < pEnd)
        {
            Cur_Ind++;
            if (*pb > El && *pb != DBL_EPSILON)
            {
                El = *pb;
                Ind_of_Max = Cur_Ind;
            }
        }
        return Ind_of_Max;
    }
    virtual void View (void) const
    {
        cout<<Type<<' ';
        Matrix::View();
    }
    double& Get_El (const int i)
    {
        if (i < N*M)
            return *(pa + i);
        double* d = new double (UINT_MAX);
        return *d;
    }
    double& operator[] (const int i)
    {
        return Get_El(i);
    }
    virtual bool Read_from_file (string& sa)
    {
        ifstream ifsa (sa.c_str(), ios::binary);
        return Read_from_file(ifsa);
    }
    virtual bool Read_from_file (const char* A)
    {
        ifstream ifsa (A, ios::binary);
        return Read_from_file(ifsa);
    }
    virtual bool Read_from_file (ifstream& sa)
    {
        delete pa;
        if (sa.eof())
            return 0;
        sa.read(&Type, sizeof(char));
        if (sa.eof())
            return 0;
        sa.read(reinterpret_cast<char*>(&N), sizeof(int));
        if (sa.eof())
            return 0;
        sa.read(reinterpret_cast<char*>(&M), sizeof(int));
        pa = new double [N*M];
        if (!pa)
            return 0;
        const double* pEnd = pa + N*M;
        const int Size = sizeof(double);
        for (double* pb = pa; pb < pEnd; ++pb)
        {
            if (sa.eof())
                return 0;
            sa.read(reinterpret_cast<char*>(pb), Size);
        }
        sa.close();
        return 1;
    }
    virtual void Write_to_file (string& a) const
    {
        ofstream ofsa (a.c_str(), ios::binary);
        return Write_to_file (ofsa);
    }
    virtual void Write_to_file (const char* A) const
    {
        ofstream ofsa (A, ios::binary);
        return Write_to_file (ofsa);
    }
    virtual void Write_to_file (ofstream& sa) const
    {
        const double* pEnd = pa + N*M;
        //string a = Name;
        //a.erase(0, 7);
        //int Num = a;
        //sa.write(reinterpret_cast<char*>(), sizeof(int));
        sa.write (&Type, sizeof(char));
        sa.write(reinterpret_cast<const char*>(&N), sizeof(int));
        sa.write(reinterpret_cast<const char*>(&M), sizeof(int));
        const int Size = sizeof(double);
        for (double* a = pa; a < pEnd; ++a)
        {
            sa.write(reinterpret_cast<char*>(a), Size);
        }
        sa.close();
    }
    virtual ~Matrix_SE (void)
    {}
};

class Matrix_SLE: public Matrix_SE
{
public:
    Matrix_SLE (void): Matrix_SE()
    {
        Type = 'L';
    }
    Matrix_SLE (Matrix<double>& A): Matrix_SE(A)
    {
        Type = 'L';
    }
    Matrix_SLE (double* pA, const int a, const int b): Matrix_SE(pA, a, b)
    {
        Type = 'L';
    }
    Matrix_SLE (Matrix_SE& A): Matrix_SE(A)
    {
        Type = 'L';
    }
    Matrix_SLE (Matrix_SLE& A): Matrix_SE(*(reinterpret_cast<Matrix_SE*>(&A)))
    {
        Type = 'L';
    }
    Matrix_SLE T (void) const
    {
        Matrix_SLE B();
        if (B.Edit_row(N) && B.Edit_col(M))
        {
        double* pb = B.pa;
        double* pbEnd = pb + N*M;
        double* paEnd = pa + N*M;
        for (double* pa1 = pa, *pb1 = pb; pa1 < paEnd; pa1+=M, ++pb1)
        {
            for (double* pa2 = pa1, *pb2 = pb1; pb2 < pbEnd; pb2+=N, ++pa2)
            {
                *pb2 = *pa2;
            }
        }
        }
        return B;
    }
    vector<double> Solve_by_Gauss_Metod (vector<double>& C)
    {
        int k = 0;
        double T = 0;
        for (int i = 0; i < N - 1; i++)
        {
            k = this->Find_Max_El_in_Col(i);
            if (Get_El(k*M+k) < 10*DBL_EPSILON)
                return vector<double> (1);
            if (i != k)
            {
                this->Replace_Rows(i, k);
                swap (C[k], C[i]);
            }
            T = Get_El(i*M+i);
            this->Multi_Row_by_Num(i, 1/T);
            C[i] = C[i]/T;
            for (k = i + 1; k < N; ++k)
            {
                T = (*this)[k*M+i];
                this->One_Minus_Second_Row(k, i, T);
                C[k] = C[k] - T*C[i];
            }
        }
        vector<double> ps (N);
        //double* pEnd = ps + N;
        //for (double* pi = 0; pi < pEnd; )
        double* pMatrix = pa + N*M - 1;
        vector<double>::reverse_iterator ivA = C.rbegin();
        vector<double>::reverse_iterator ivBeg = ps.rbegin(); //const_?
        vector<double>::const_reverse_iterator ivEnd = ps.rend();
        k = N - 1; //////////////////////////////////////////////////
        k = M;
        for (vector<double>::reverse_iterator pt = ivBeg; pt < ivEnd; ++pt, pMatrix -= (k--))
        {
            *pt = *(ivA++);
            for (vector<double>::reverse_iterator pc = ivBeg; pc < pt; --pMatrix, ++pc)
                *pt = *pt - *pMatrix*(*pc);
            *pt = *pt/(*pMatrix);
        }
        return ps;
    }
    virtual ~Matrix_SLE(void)
    {}
};
