#ifndef MATH_H
#define MATH_H
#include "../Lib/Matrix.h"
#include "Different.h"
using namespace std;
double Norm (vector<double> vect)
{
    Matrix<double> A(vect);
    double* pEnd = A.get_pointer() + A.get_n() * A.get_m();
    double a = *A.get_pointer();
    for (double* pb = A.get_pointer() + 1; pb < pEnd; ++pb)
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
    Matrix_SE (const Matrix<double>& A): Matrix<double>(A), High_accuracy(1), Type('N')
    {}
    Matrix_SE (double* pA, const int a, const int b): Matrix<double>(pA, a, b),  High_accuracy(1), Type('N')
    {}
    Matrix_SE (const Matrix_SE& A): Matrix<double>(*(reinterpret_cast<const Matrix<double>*>(&A))), High_accuracy (A.High_accuracy), Type
            ('N')
    {}
    //    Matrix_SE(Matrix_SE&& mmty) noexcept: High_accuracy(mmty.High_accuracy), Type(mmty.Type)
//    {}
    Matrix_SE (const int N1, const int M1): Matrix<double> (N1, M1), High_accuracy(1), Type('N')
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
        double* pb = p_data + N_R * M;
        double* pendb = pb + M;
        for (double* pt = pb; pt < pendb; ++pt)
        {
            *pt = C*(*pt);
        }
        return;
    }
    void One_Minus_Second_Row (const int one, const int second, const double C)
    {
        double* p1 = p_data + one * M;
        const double* pe1 = p1 + M;
        double* p2 = p_data + second * M;
        for (; p1 < pe1; ++p1, ++p2)
        {
            *p1 = *p1 - C**p2;
        }
        return;
    }
    void Replace_Rows (const int One, const int Second)
    {
        double* pe1 = p_data + (One + 1) * M;
        double* pe2 = p_data + (Second + 1) * M;
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
        double* pb = p_data + Col * (M + 1);
        const double* pEnd = p_data + N * M;
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
    virtual const Matrix_SE& View (const int i = 8, const bool Show_Name = 0) const
    {
        //cout<<Type<<' ';/////////////////////////////////////////////////////
        Matrix::view(i, Show_Name);
        return *this;
    }
    virtual bool Read_from_file (string& sa)
    {
        ifstream ifsa (sa.c_str(), ios::binary);
        bool Temp = Read_from_file(ifsa);
        ifsa.close();
        return Temp;
        return Read_from_file(ifsa);
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
        sa.read(&Type, sizeof(char));
        if (sa.eof())
            return 0;
        sa.read(reinterpret_cast<char*>(&N), sizeof(int));
        if (sa.eof())
            return 0;
        sa.read(reinterpret_cast<char*>(&M), sizeof(int));
        delete p_data;
        p_data = new double [N * M];
        if (!p_data)
            return 0;
        const double* pEnd = p_data + N * M;
        const int Size = sizeof(double);
        for (double* pb = p_data; pb < pEnd; ++pb)
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
        const double* pEnd = p_data + N * M;
        //string a = name;
        //a.erase(0, strlen("Matrix"));
        //int Num = atoi(a); ///?
        //sa.write(reinterpret_cast<char*>(&Num), sizeof(int));
        sa.write (&Type, sizeof(char));
        sa.write(reinterpret_cast<const char*>(&N), sizeof(int));
        sa.write(reinterpret_cast<const char*>(&M), sizeof(int));
        const int Size = sizeof(double);
        for (double* a = p_data; a < pEnd; ++a)
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
    Matrix_SLE (const Matrix<double>& A): Matrix_SE(A)
    {
        Type = 'L';
    }
    Matrix_SLE (double* pA, const int a, const int b): Matrix_SE(pA, a, b)
    {
        Type = 'L';
    }
    Matrix_SLE (const Matrix_SE& A): Matrix_SE(A)
    {
        Type = 'L';
    }
    Matrix_SLE (const Matrix_SLE& A): Matrix_SE(*(reinterpret_cast<const Matrix_SE*>(&A)))
    {
        Type = 'L';
    }
//    Matrix_SLE(Matrix_SLE&& A) noexcept{}
    Matrix_SLE (const int N1, const int M1): Matrix_SE(N1, M1)
    {
        Type = 'L';
    }
    Matrix_SLE T (void) const
    {
        Matrix_SLE B (M, N);
        double* pb = B.p_data;
        double* pbEnd = pb + N*M;
        double* paEnd = p_data + N * M;
        for (double* pa1 = p_data, *pb1 = pb; pa1 < paEnd; pa1+=M, ++pb1)
        {
            for (double* pa2 = pa1, *pb2 = pb1; pb2 < pbEnd; pb2+=N, ++pa2)
            {
                *pb2 = *pa2;
            }
        }
        return B;
    }
    Matrix<double> Solve (const Matrix<double>& b) const
    {
        return (if_symmetric()) ? LDLT(b) : Solve_by_Gauss_Method(b);
    }
    Matrix<double> Solve_by_Gauss_Method (Matrix_SE C) const
    {
        Matrix_SLE S (*(const_cast<Matrix_SLE*>(this)));
        int k = 0;
        double T = 0;
        for (int i = 0; i < N - 1; i++)
        {
            k = S.Find_Max_El_in_Col(i);
            if (fabs(S[k*M+k]) < 10*DBL_EPSILON)
            {
                cout<<"Ne sovmestnaya systema"<<endl;
                assert (fabs(S[k*M+k]) < 10*DBL_EPSILON);
            }
            if (i != k)
            {
                S.Replace_Rows(i, k);
                C.Replace_Rows(i, k);
            }
            T = S[i*M+i];
            S.Multi_Row_by_Num(i, 1/T);
            C.Multi_Row_by_Num(i, 1/T);
            for (k = i + 1; k < N; ++k)
            {
                T = S[k*M+i];
                S.One_Minus_Second_Row(k, i, T);
                C.One_Minus_Second_Row(k, i, T);
            }
        }
        Matrix<double> ps (N, 1);
        const double* pMatrix = S.get_pa() + S.get_n() * S.get_m() - 1;
        const double* pC = C.get_pa() + C.get_n() * 1 - 1;
        double* const pBeg = ps.get_pa() + ps.get_n() * 1 - 1; //const_?
        double* const pEnd = ps.get_pa();
        //assert(k == N - 1);
        k = M;
        for (double* pt = pBeg; pt >= pEnd; --pt, pMatrix -= (k--))
        {
            *pt = *(pC--);
            for (double* pc = pBeg; pc > pt; --pMatrix, --pc)
                *pt = *pt - *pMatrix*(*pc);
            *pt = *pt/(*pMatrix);
        }
        return ps;
    }
    Matrix<double> LDLT (const Matrix<double>& b) const
    {
        Matrix_SLE L (N, M);
        Matrix_SLE D (N, M);
        Matrix_SLE LT (N, M);
        double* plEnd = L.get_pa() + N * M;
        double* pl = L.get_pa();
        double* pd = D.get_pa();
        double* plt = LT.get_pa();
        for (int i = 0; i < N - 1; )
        {
            pd[i*M + i] = p_data[i * M + i];
            for (int k1 = 0; k1 < i; ++k1)
            {
                pd[i*M + i] -= L[i*M + k1]*L[i*M + k1]*pd[k1*M+k1];
            }
            ++i;
            for (int k = 0; k < i; ++k)
            {
                L[i*M+k] = p_data[i * M + k];
                for (int t = 0; t < k; ++t)
                    L[i*M+k] -= L[i*M+t]*L[k*M+t]*pd[t*M+t];
                L[i*M+k] /= pd[k*M+k];
            }
        }
        int i = N - 1;
        pd[i*M + i] = p_data[i * M + i];
        for (int k1 = 0; k1 < i; ++k1)
        {
            pd[i*M + i] -= L[i*M + k1]*L[i*M + k1]*pd[k1*M+k1];
        }
        for (; pl < plEnd; pl+=M+1, plt+=M+1)
        {
            *pl = 1;
            *plt = 1;
        }
        LT = L.T();
//        L.view();
//        D.view();
//        LT.view();
//        b.view();
//        (L*D*LT).view();
        Matrix_SLE Y (L.get_n(), 1);
        Y[0] = b[0];
        auto pLrow = L.first_i() + M;
        auto pLEnd = L.first_i() + M;
        for (auto pb = b.first_i() + 1, pY = Y.first_i() + 1; pb <= b.last_i(); ++pb, ++pY)
        {
            *pY = *pb;
            for (auto pL = pLrow, pY1 = pY - 1; pL >= pLEnd; --pL, --pY1)
            {
//                if (pL == pLEnd)
//                    assert (pY1 == Y.first_i());
//                assert (*pL != 0);
//                assert (*pL != 1);
                *pY -= *pL**pY1;
            }
            pLEnd += get_m();
            pLrow += get_m() + 1;
        }
//        Y.view();
        for (auto pD = D.first_i(), pY = Y.first_i(); pD <= D.last_i(); pD += M + 1, ++pY)
        {
//            assert (*pD != 0);
            *pY = *pY/(*pD);
        }
//        Y.view();
        auto pLTrow = LT.last_i() - M;
        auto pLTEnd = LT.last_i() - M;
        for (auto pY = Y.last_i() - 1; pY >= Y.first_i(); --pY)
        {
            for (auto pLT = pLTrow, pY1 = pY + 1; pLT <= pLTEnd; ++pLT, ++pY1)
            {
//                if (pLT == pLTEnd)
//                    assert (pY1 == Y.last_i());
//                assert (*pLT != 0);
//                assert (*pLT != 1);
                *pY -= *pLT**pY1;
            }
            pLTEnd -= get_m();
            pLTrow -= get_m() + 1;
        }
//        Y.view();
        return Y;
    }
    vector<double> Solve_by_Gauss_Method_v (vector<double>& C)
    {
        int k = 0;
        double T = 0;
        for (int i = 0; i < N - 1; i++)
        {
            k = this->Find_Max_El_in_Col(i);
            if (operator[](k*M+k) < 10*DBL_EPSILON)
                return vector<double> (1);
            if (i != k)
            {
                this->Replace_Rows(i, k);
                swap (C[k], C[i]);
            }
            T = operator[](i*M+i);
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
        double* pMatrix = p_data + N * M - 1;
        vector<double>::reverse_iterator ivA = C.rbegin();
        vector<double>::reverse_iterator ivBeg = ps.rbegin(); //const_?
        vector<double>::const_reverse_iterator ivEnd = ps.rend();
//        assert (k == N - 1); //////////////////////////////////////////////////
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
    virtual const Matrix_SLE& View (const int i = 8, const bool Show_Name = 0) const
    {
        Matrix_SE::View(i, Show_Name);
        return *this;
    }
    virtual ~Matrix_SLE(void)
    {}
};

class Function
{
    double (*Fn)(double);
    static inline double NullF1 (double x)
    {
        return x-x;
    }
public:
    Function (double (*Fn1) (double)): Fn (Fn1)
    {}
    Function (void): Fn (NullF1)
    {}
    Function Return_Null_Function(void)
    {
        return NullF1;
    }
    double Derivative_by_user_def (double (*Der) (double), const double x) const
    {
        return Der(x);
    }
    double Derivative_by_definition (const double x) const
    {
        double Eps = 0.1;
        return (Fn(x + Eps) - Fn(x))/Eps;
    }
    double f(const double x) const
    {
        return Fn(x);
    }
    double Integral_by_user_def (double (*AntiDer) (double, double), const double from, const double to) const;
    double Integral_by_Sympthon_method (const double from, const double to, const unsigned long long Max_N = 100000000ULL) const
    {
        unsigned long long N = 1ULL;
        auto H = N;
        double Integral = 0;
        double Length = to - from;
        while(Length/N > DBL_EPSILON && N < Max_N)// 10000000000 - very big deflection (error)!!!
        {
            N *= 2ULL;
        }
        N -=  2ULL;
        double h = Length/N;
        for (double x = from + h; x < to; x+=h, H-=2ULL)
        {
            Integral += 4*f(x);
            x+=h;
            Integral += 2*f(x);
        }
        Integral += f(from);
        Integral += f(to);
        Integral *= h/3;
        return Integral;
    }
    double Integral_by_Sympthon_method_with_epsilon(const double from, const double to, const double Epsilon) const
    {
        unsigned long long h = 1000ULL;
        double R = 0;
        double R1 = 0;
        double Epsln = 15*Epsilon;
        while (Epsln < fabs((R = Integral_by_Sympthon_method(from, to, h)) - R1))
        {
            h = 2*h;
            R1 = R;
        }
        return R;
    }
    double Integral_by_trapetia_method (const double from, const double to, const unsigned long long Max_N = 100000000ULL) const
    {
#ifndef LLONG_MAX
#define LLONG_MAX numeric_limits<long long>::max()
#endif
#ifndef ULLONG_MAX
#define ULLONG_MAX numeric_limits<unsigned long long>::max()
#endif
        unsigned long long N = 1ULL;
        double Integral = 0;
        double Length = to - from;
        //__UINT64_TYPE__
        while(Length/N > DBL_EPSILON && N < Max_N)// 10000000000 - very big deflection (error)!!! 10 - 5
        {
            N *= 2ULL;
        }
        N -=  2ULL;
        double Min = numeric_limits<double>::epsilon();
        double h = Length/N;
        unsigned long long H = N;
        for (double x = from + h; x < to; x += h, --H)
        {
            Integral += 2*f(x);
        }
        if (H > LLONG_MAX)
            H = LLONG_MAX - H;
        if (H > 0 && H < 100*Length)
        {
            while (--H != ULLONG_MAX)
            {
                Integral += 2*f(to - Min);
            }
            ++H;
        }
        assert (H == 0);
        Integral += f(from);
        Integral += f(to);
        Integral *= h/2;
        return Integral;
    }
    double Integral_by_trapetia_method_with_epsilon(const double from, const double to, const double Epsilon) const
    {
        unsigned long long h = 1000ULL;
        double R = 0;
        double R1 = 0;
        double Epsln = 3*Epsilon;
        while (Epsln < fabs((R = Integral_by_trapetia_method(from, to, h)) - R1))
        {
            h = 2*h;
            R1 = R;
        }
        return R;
    }
    friend ostream& operator<< (ostream& a, Function Fn);
    ~Function (void)
    {
        Fn = NullF1;
    }
};



class Matrix_SNE: public Matrix<Function>
{
    char Type;
public:
    Matrix_SNE (void): Matrix<Function>(), Type('N')
    {}
    Matrix_SNE (const int N1, const int M1): Matrix<Function>(N1, M1), Type('N')
    {}
    Matrix_SNE (Matrix_SNE& Ar): Matrix<Function> (*(reinterpret_cast< Matrix<Function>* > (&Ar))), Type ('N')
    {}
//    Function operator[] (const int S) const
//    {
//        if (S < N)
//            return NullF2;
//        return *(Fn + S);
//    }
    Matrix<double> f (const double x) const
    {
        Matrix<double> A (N, 1);
        double* pA = A.get_pa() + N - 1;
        for (Function* ppa = this->p_data + this->N - 1; ppa >= this->p_data; --(ppa), --pA)
        {
            *pA = ppa->f(x);
        }
        return A;
    }
//    Matrix_SLE Derivative (const Matrix<double> X) const
//    {
//        int Num_of_variables = X.get_size();
//        Matrix<double> A (N, Num_of_variables);
//        double* pa = A.get_pa() + A.get_size() - 1;
//        Matrix<double> Temp (Num_of_variables, 1);
//        double* pTemp = 0;
//        for (Function2* ppa = a + N*Num_of_variables - 1; ppa >= pa; --ppa, --pa)
//        {
//            Temp = ppa->Derivative_by_definition_M_in_column(X);
//            pTemp = Temp.get_pa() + Num_of_variables - 1;
//            for (int i = 0; i < Num_of_variables; ++i, --pa, --pTemp)
//            {
//                *pa = *pTemp;
//            }
//        }
//        return A;
//    }
//    Matrix<double> Integral (const double from_x, const double from_y, const double to_x, const double to_y) const;
    double Max (const double x)
    {
        double m = 0;
        double v = 0;
        for (Function* pA = p_data + N - 1; pA >= p_data; --pA)
        {
            if (m < (v = pA->f(x)))
                m = v;
        }
        return m;
    }
};
class Function2
{
    double (*Fn)(const Matrix<double>&);
    static inline double NullF2 (const Matrix<double>& x)
    {
        return 0;
    }
public:
    Function2 (double (*Fn1) (const Matrix<double>& x)): Fn (Fn1)
    {}
    Function2 (void): Fn (NullF2)
    {}
    Function2 Return_Null_Function(void)
    {
        return NullF2;
    }
//    double Derivative_by_user_def (double (*Der) (const Matrix<double>&), const Matrix<double> x)
//    {
//        return Der(x);
//    }
    Matrix<double> Derivative_by_definition_M_in_column (const Matrix<double>& X) const
    {
        int Num_of_variables = X.get_size();
        Matrix<double> A (Num_of_variables, 1);
        double* pa = A.get_pa();
        const double M = 0.01; //0.05   0.1////////////////////////////////////////////////////////////////////////
        double F = Fn(X);
        Matrix<double> DX(X);
        double* pEnd = DX.get_pa() + DX.get_size();
        for (double* pDX = DX.get_pa(), *pX = X.get_pa(); pDX < pEnd; ++pDX, ++pX, ++pa)
        {
            if( *pDX < 10*DBL_EPSILON)
            {
                *pDX += 1e-6;
                *pa = (Fn(DX) - F)/1e-6;
                *pDX = *pX;
            }
            else
            {
                *pDX += M*(*pDX);
                *pa = (Fn(DX) - F)/(M*(*pX));
                *pDX = *pX;
            }
        }
//        else
//            for (double* pDX = DX.Get_pa(), *pX = X.get_pa(); pDX < pEnd; ++pDX, ++pX, ++pa)
//            {
//
//            }
        return A;
    }
    double f(const Matrix<double>& X) const
    {
        return Fn(X);
    }
    double Integral_by_user_def (double (*AntiDer) (const Matrix<double>&), const Matrix<double>& From, const Matrix<double>& To) const;
    double Integral_by_definition (const Matrix<double>& From, const Matrix<double>& To, const unsigned long long x_max_steps=1000ULL, const unsigned long long y_max_steps=1000ULL) const //By Sympthon method, only for 2 variables
    {
        assert (To.get_size() == 2 && From.get_size() == 2);
        double Integral = 0;
        unsigned long long xN = 1;
        double xLength = To[0] - From[0];
        while(xLength/xN > DBL_EPSILON && xN < x_max_steps)// 10'000'000'000 - very big deflection (error)!!!
        {
            xN *= 2;
        }
        double xh = xLength/xN;
        double xBegin = From[0];
        double xEnd = To[0];
        unsigned long long yN = 1;
        double yLength = To[1] - From[1];
        while(yLength/yN > DBL_EPSILON && yN < y_max_steps)// 10000000000 - very big deflection (error)!!!
        {
            yN *= 2;
        }
        double yh = yLength/yN;
        double yEnd = To[1];
        Matrix<double> Another_Temp;
        for (Matrix<double> Temp = From; Temp.unsafe_index(1) < yEnd; Temp.unsafe_index(1)+= yh + yh)
        {
            for (Temp.unsafe_index(0) = xBegin; Temp.unsafe_index(0) < xEnd; Temp.unsafe_index(0)+= xh + xh)
            {
                Another_Temp = Temp;
                Integral += f(Another_Temp);
                Another_Temp.unsafe_index(0) += xh;
                Integral += 4*f(Another_Temp);
                Another_Temp.unsafe_index(0) += xh;
                Integral += f(Another_Temp);

                Another_Temp.unsafe_index(0) = Temp.unsafe_index(0);
                Another_Temp.unsafe_index(1) += yh;
                Integral += 4*f(Another_Temp);
                Another_Temp.unsafe_index(0) += xh;
                Integral += 16*f(Another_Temp);
                Another_Temp.unsafe_index(0) += xh;
                Integral += 4*f(Another_Temp);

                Another_Temp.unsafe_index(0) = Temp.unsafe_index(0);
                Another_Temp.unsafe_index(1) += yh;
                Integral += f(Another_Temp);
                Another_Temp.unsafe_index(0) += xh;
                Integral += 4*f(Another_Temp);
                Another_Temp.unsafe_index(0) += xh;
                Integral += f(Another_Temp);
            }
        }
        return Integral*=xh*yh/9;
    }
    ~Function2 (void)
    {
        Fn = NullF2;
    }
};
class Array_of_Functions2
{
    int N;
    Function2* Fn;
public:
    Array_of_Functions2 (void): N(1), Fn (new Function2[N])
    {}
    Array_of_Functions2 (const int S): N(S), Fn (new Function2[N])
    {}
    Array_of_Functions2 (const Array_of_Functions2& Ar): N(Ar.N), Fn (new Function2[N])
    {
        for (Function2* pFn = Fn + N - 1, *pFa = Ar.Fn +  Ar.N - 1; pFn >= Fn; --pFn, --pFa)
        {
            *pFn = *pFa;
        }
    }
    Array_of_Functions2& operator= (const Array_of_Functions2& Arr)
    {
        N = Arr.N;
        if (Fn == Arr.Fn)
        {
            return *this;
        }
        delete[] Fn;
        Function2* pFm = Arr.Fn + Arr.N - 1;
        for (Function2* pFn = Fn + N - 1; pFn >= Fn; --pFn, --pFm)
        {
            *pFn = *pFm;
        }
        return *this;
    }
    Array_of_Functions2& operator= (Array_of_Functions2& Arr)
    {
        N = Arr.N;
        if (Fn == Arr.Fn)
        {
            return *this;
        }
        delete[] Fn;
        Function2* pFm = Arr.Fn + Arr.N - 1;
        for (Function2* pFn = Fn + N - 1; pFn >= Fn; --pFn, --pFm)
        {
            *pFn = *pFm;
        }
        return *this;
    }
    inline int Get_N (void) const
    {
        return N;
    }
    inline Function2* Get_Fn (void) const
    {
        return Fn;
    }
    inline Function2* Get_Pointer (void) const
    {
        return Fn;
    }
    Function2& operator[] (const int S) const
    {
        assert (S > -1 && S < N);
        return *(Fn + S);
    }
    Matrix<double> f (const Matrix<double>& X) const //return matrix in ONE column
    {
        Matrix<double> A (N, 1);
        double* pa = A.get_pa() + N - 1;
        for (Function2* pFn = Fn + N - 1; pFn >= Fn; --pFn, --pa)
        {
            *pa = pFn->f(X);
        }
        return A;
    }
    Matrix<double> Derivative (const Matrix<double>& X) const //return matrix in Num_of_variables column and N rows
    {
        const int Num_of_variables = X.get_size();
        Matrix<double> A (N, Num_of_variables);
        double* pa = A.get_pa() + A.get_size() - 1;
        Matrix<double> Temp (Num_of_variables, 1);
        double* pTemp = nullptr;
        for (Function2* pFn = Fn + N - 1; pFn >= Fn; --pFn)
        {
            Temp = pFn->Derivative_by_definition_M_in_column(X);
            pTemp = Temp.get_pa() + Num_of_variables - 1;
            for (int i = 0; i < Num_of_variables; ++i, --pa, --pTemp)
            {
                *pa = *pTemp;
            }
        }
        return A;
    }
    Matrix<double> Integral (const Matrix<double> From, const Matrix<double> To) const //return matrix in ONE column
    {
        assert (From.get_size() == To.get_size() && (From.get_m() == To.get_m() || From.get_n() == To.get_n()));
        Matrix<double> A (N, 1);
        double* pa = A.last_i();
        for (Function2* pf = Fn +  N - 1; pf == Fn; --pf, --pa)
        {
            *pa = pf->Integral_by_definition(From, To);
        }
        return A;
    }
    double Max (const Matrix<double>& X) const
    {
        double m = 0;
        double v = 0;
        for (Function2* pA = Fn + N - 1; pA >= Fn; --pA)
        {
            if (m < (v = fabs(pA->f(X))))
                m = v;
        }
        return m;
    }
    ~Array_of_Functions2 (void)
    {
        delete[] Fn;
    }
};
namespace Solve_Nonlinear_Equations
{
double fDelta2 (const Matrix<double>& X_i, const Matrix<double>& Delta)
{
    double Delta2 = 0;
    double Temp = 0;
    for (double* x = X_i.get_pa() + X_i.get_size() - 1, *dx = Delta.get_pa() + Delta.get_size() - 1; x >= X_i.get_pa(); --x, --dx)
    {
        if (fabs(*x + *dx) < 1 && fabs(*dx) > Delta2)
        {
            Delta2 = *dx;
        }
        else if ((Temp = fabs(*dx/(*x+*dx))) > Delta2)
        {
            Delta2 = Temp;
        }
    }
    return Delta2;
}
Matrix<double> Solve_SNE_Without_Derivative (const Array_of_Functions2& Func, const Matrix<double>& xy, const double& Eps1, const double& Eps2)
{
    assert (Func.Get_Pointer() != nullptr);
    assert (xy.get_pointer() != nullptr);
    assert (xy.get_size() != 0);
    assert (Eps1 > 1e-10);
    assert (Eps2 > 1e-10);
    const int Num_it = 10000;
    int i = 0;
    Matrix<double> ixy (xy);
    Matrix<double> Delta(xy.get_n(), xy.get_m() - 1);
    Matrix_SLE Derivative;
    Matrix<double> Sol;
    double Delta1 = 0;
    double Delta2 = 0;
    do
    {
        cout<<setw(10)<<Delta1<<setw(10)<<Delta2<<setw(6)<<++i<<endl;
//        cout<<"ixy"<<endl;
//        ixy.view();
//        cout<<"Delta"<<endl;
//        Delta.view();
        ixy = ixy + Delta;
//        cout<<"Sol"<<endl;
        Sol = Func.f(ixy).multiple_matrix_by_number(-1).view();
//        cout<<"Derivative"<<endl;
        Derivative = Func.Derivative(ixy);
        Delta = Derivative.Solve(Sol);       //.multiple_matrix_by_number(-1).view();
        cout<<endl;
    }
    while ((Delta1 = Func.Max(ixy)) > Eps1 && (Delta2 = fDelta2(ixy, Delta)) > Eps2 &&  i <= Num_it);
    return ixy;
}
Matrix<double> Solve_SNE (const Array_of_Functions2& Func, const Array_of_Functions2& Der_of_Arr, const Matrix<double>& xy, const double Eps1, const double Eps2)
{
    // by Newton method
    const int Num_it = 10000;
    int i = 0;
    Matrix<double> ixy (xy);
    Matrix<double> Delta(xy.get_n(), xy.get_m());
    Matrix<double> High_Precision_ixy (xy);
    Matrix<double> High_Precision_Delta(xy.get_n(), xy.get_m());
    Matrix_SLE High_Precision_Derivative (xy.get_n(), xy.get_n()); //It is not a mistake!
    Matrix_SLE Derivative;
    Matrix<double> Sol;
    double Delta1 = 0;
    double Delta2 = 0;
    do
    {
        cout<<setw(10)<<Delta1<<setw(10)<<Delta2<<setw(6)<<++i<<endl;
//        cout<<"ixy"<<endl;
//        ixy.view();
//        cout<<"Delta"<<endl;
//        Delta.view();
        ixy = ixy + Delta;
//        cout<<"Sol"<<endl;
        Sol = Func.f(ixy).multiple_matrix_by_number(-1).view();
        Delta1 = Sol[0];
//        cout<<"Derivative"<<endl;
        Derivative = Func.Derivative(ixy);
        Delta = Derivative.Solve(Sol);
//        cout<<"High_Precision_ixy"<<endl;
//        High_Precision_ixy.view();
//        cout<<"High_Precision_Delta"<<endl;
//        High_Precision_Delta.view();
        High_Precision_ixy = High_Precision_ixy + High_Precision_Delta;
        Function2* pF = Der_of_Arr.Get_Fn() + Der_of_Arr.Get_N() - 1;
        for (double* pD = High_Precision_Derivative.get_pa() + High_Precision_Derivative.get_size() - 1; pD >=
                                                                                                         High_Precision_Derivative.get_pa(); --pD, --pF)
        {
            *pD = -pF->f(High_Precision_ixy);// Minus against write Sol = Func.f(High_Precision_ixy)
        }                         //.multiple_matrix_by_number(-1).view();
        cout<<"With user defined derivative F."<<endl;
        Sol = Func.f(High_Precision_ixy).view();
//        cout<<"High_Precision_Derivative"<<endl;
        High_Precision_Delta = High_Precision_Derivative.Solve(Sol);
        if (fabs(Delta1) > fabs(Sol[0]))
        {
            cout<<"Better"<<endl;
        }
        cout<<endl;
    }
    while ((Delta1 = Func.Max(ixy)) > Eps1 && (Delta2 = fDelta2(ixy, Delta)) > Eps2 &&  i <= Num_it);
    return ixy;
}
}


namespace Solve_Differential_Equations
{
#include<vector>


vector<Matrix<double> > Explicit_Euler_method (const Array_of_Functions2& F, const double x_From, const double x_To, const Matrix<double>& u0,const Matrix<double>& Epsilon, const double Max_Step)
{
//    static_assert (F.Get_Pointer() != nullptr, );
    assert (F.Get_Pointer() != nullptr);
    assert (u0.get_pointer() != nullptr);
    assert (Epsilon.get_pointer() != nullptr);
    assert (x_From < x_To);
    assert (F.Get_N() == Epsilon.get_size());
    assert (F.Get_N() == u0.get_size());
    assert (Max_Step > DBL_EPSILON);
    const int N = F.Get_N();
    double x = x_From;
    const int Iterations = 150000;
    Matrix<double> yk (N, Iterations);
    Matrix<double> tk (1, Iterations);
    Matrix<double> F_curr = u0;
    Matrix<double> tk_curr (1, N);
    Matrix<double> y_curr (1, N+1);
    auto tk_i = tk.first_i();
    *tk_i = x_From;
    ++tk_i;
    const auto Eps_End = Epsilon.last_i();
    auto yk_First = yk.first_i();
    auto yk_End = yk.first_i() + N;
    for (auto iu = u0.first_i() + (N - 1), iy = yk_End - 1; iu >= u0.first_i(); --iu, --iy)
    {
        *iy = *iu;
    }
    int i = 0;
    while (x < x_To)
    {
        for (auto iy = y_curr.first_i(), iyk = yk_First; iyk < yk_End; ++iy, ++iyk)
        {
            *iy = *iyk;
        }
        *y_curr.last_i() = x;
        yk_First += N;
        yk_End += N;
        F_curr = F.f(y_curr);
        for (auto it = tk_curr.first_i(), iF = F_curr.first_i(), iEps = Epsilon.first_i(); iEps <= Eps_End; ++it, ++iF, ++iEps)
        {
            *it = *iEps/(fabs(*iF)+*iEps/Max_Step);
        }
        *tk_i = tk_curr.min_element();
        for (auto it = yk_First, iF = F_curr.first_i(); it < yk_End; ++it, ++iF)
        {
            *it = *(it-N) + *tk_i**iF;
        }
        x += *tk_i;
        *tk_i = x;
        ++tk_i;
//        tk.edit_row(2*I);
        ++i;
    }
    cout<<++i<<endl;
    vector<Matrix<double> > A (2);
    A.at(0) = tk;
    A.at(1) = yk;
    return A;
}
double absol (double x)
{
    return x > 0 ? x : -x;
}
vector<Matrix<double> > Implicit_Euler_method (const Array_of_Functions2& F, const double x_From, const double x_To, const Matrix<double>& u0, const Matrix<double>& Epsilon, const double t_min, const double t_max, const bool Strategy = 0)
{
    //u0 in column, because of 25 string down!!!
    assert (F.Get_Pointer() != nullptr);
    assert (u0.get_pointer() != nullptr);
    assert (Epsilon.get_pointer() != nullptr);
    assert (F.Get_N() == Epsilon.get_size());
    assert (F.Get_N() == u0.get_size());
    assert (x_From < x_To);
    assert (t_min > DBL_EPSILON);
    assert (t_max > DBL_EPSILON);
    assert (t_min < t_max);
    for (auto ie = Epsilon.last_i(); ie >= Epsilon.first_i(); --ie)
        assert (*ie > 0);
    const int N = F.Get_N();
    int Itrtns = 20000;
    int Itrtns_curr = 0;
    double tk = t_min; //
    double tk_prv = t_min; //
    double x = x_From;
    Matrix<double> tk_curr_next (1, N); //noth
    Matrix<double> y_prv_curr (u0); //
    Matrix<double> y_curr (u0); //
    Matrix<double> y_next_curr (u0); //noth
    Matrix<double> y_next_curr_with_x (u0); //noth
    Matrix<double> Epsln_curr(Epsilon.get_n(), Epsilon.get_m()); //noth
    Matrix<double> yk (N, Itrtns);
    Matrix<double> tk_mtrx (1, Itrtns);
    y_next_curr_with_x.edit_row(N + 1); //25th string!!!!!!!!!!!
    y_next_curr_with_x.unsafe_index(N) = x_From;
    auto yk_First = yk.first_i() + N;
    auto yk_Last = yk_First;
    for (auto iyk = yk.first_i(), iu0 = u0.first_i(); iyk < yk_First; ++iyk, ++iu0)
    {
        *iyk = *iu0;
    }
    tk_mtrx.unsafe_index(0) = x_From;
    auto itk_mtrx = tk_mtrx.first_i() + 1;
    const int Num_it = 100;
    int i = 0;
    const double Eps1 = 1e-4;
    const double Eps2 = 1e-4;
    double Delta1 = 0;
    double Delta2 = 0;
    Matrix<double> Delta(u0.get_n(), u0.get_m());
    Matrix<double> F_curr;
    Matrix_SLE Drvtve;
    Matrix<double> Strange_matrix(2*N+2, 1);//x, t, k
    auto iStr = Strange_matrix.first_i();
    while (y_next_curr_with_x.unsafe_index(N) < x_To)
    {

        y_next_curr_with_x.unsafe_index(N) += tk;
        *itk_mtrx = y_next_curr_with_x.unsafe_index(N);
        ++itk_mtrx;
        iStr = Strange_matrix.last_i();
        *iStr = tk;
        --iStr;
        *iStr = *y_next_curr_with_x.last_i();
        --iStr;
        for (auto iy_curr = y_curr.last_i(); iy_curr >= y_curr.first_i(); --iy_curr, --iStr)
        {
            *iStr = *iy_curr;
        }
eps_curr_bigger_then_epsilon:
        do
        {
//            for (auto i = y_next_curr.Last_i(), ix = y_next_curr_with_x.last_i()-1; i >= y_next_curr.first_i(); --i, --ix)
//            {
//                *ix = *i;
//            }
            y_next_curr = (y_next_curr + Delta);
            iStr = Strange_matrix.first_i() + (N - 1);
            for (auto iy_next_curr = y_next_curr.last_i(); iy_next_curr >= y_next_curr.first_i(); --iy_next_curr, --iStr)
            {
                *iStr = *iy_next_curr;
            }
//            Strange_matrix.view();
//        cout<<"F_curr"<<endl;
            F_curr = F.f(Strange_matrix).multiple_matrix_by_number(-1.0)/*.view()*/;
//        cout<<"Derivative"<<endl;
            Drvtve = F.Derivative(Strange_matrix);
            Drvtve.edit_col(N);
//            Drvtve.view();
            Delta = Drvtve.Solve(F_curr);       //.multiple_matrix_by_number(-1.0).view();
//            Delta.view();
//            cout<<Delta1<<' '<<Delta2<<endl<<endl;
        }
        while (++i <= Num_it && ((Delta1 = F_curr.max_element_abs()) > Eps1) && ((Delta2 = Solve_Nonlinear_Equations::fDelta2(y_next_curr, Delta)) > Eps2));
        Delta.fill_nulls();
//        cout<<i<<endl;
        i = 0;
//        y_next_curr.view();
        for (auto ie = Epsln_curr.last_i(), iy_prv = y_prv_curr.last_i(), iy_curr = y_curr.last_i(), iy_next = y_next_curr.last_i(); ie >=
                                                                                                                                     Epsln_curr.first_i(); --ie, --iy_prv, --iy_curr, --iy_next)
        {
            *ie = tk/(tk+tk_prv)*((tk/tk_prv)*(*iy_curr - *iy_prv) + *iy_curr - *iy_next);
        }
//        Epsln_curr.view();
//        Epsilon.view();
        if (!Strategy)
        {
            for (auto ie_curr = Epsln_curr.last_i(), ie = Epsilon.last_i(); ie >= Epsilon.first_i(); --ie_curr, --ie)
            {
                if (*ie_curr > *ie)
                {
                    tk /= 2;
                    y_next_curr = y_curr;
                    for (auto inx = y_next_curr.last_i(), ix = y_next_curr_with_x.last_i() - 1; inx >= y_next_curr.first_i(); --inx, --ix)
                    {
                        *ix = *inx;
                    }
                    y_next_curr_with_x.unsafe_index(N) -= tk;
                    goto eps_curr_bigger_then_epsilon;
                }
            }
//        if (!Strategy)
            for (auto itk_curr = tk_curr_next.last_i(), iEps_curr = Epsln_curr.last_i(), iEps = Epsilon.last_i(); iEps >= Epsilon.first_i(); --itk_curr, --iEps, --iEps_curr)
            {
                *itk_curr = sqrt(*iEps/absol(*iEps_curr))*tk;
            }
        }
        else
            for (auto itk_curr = tk_curr_next.last_i(), iEps_curr = Epsln_curr.last_i(), iEps = Epsilon.last_i(); iEps >= Epsilon.first_i(); --itk_curr, --iEps, --iEps_curr)
            {
                *itk_curr = absol(*iEps_curr) > *iEps ? tk/2 : absol(*iEps_curr) > *iEps/4 ? tk : tk+tk;
            }
        tk_prv = tk;
//        tk_curr_next.view();
        tk = tk_curr_next.min_element();
        if (tk > t_max)
            tk = t_max;
//        y_curr.view();
//        y_next_curr.view();
//        cout<<tk_prv<<' '<<tk<<' '<<endl;
//        tk_curr_next.view();
        y_prv_curr = y_curr;
        y_curr = y_next_curr;
        yk_Last += N;
        for (auto iyk = yk_First, iy_curr = y_curr.first_i(); iyk < yk_Last; ++iyk, ++iy_curr)
        {
            *iyk = *iy_curr;
        }
        yk_First += N;
        ++Itrtns_curr;
    }
    cout<<Itrtns_curr<<endl;
    vector<Matrix<double> > A (2);
    A.at(0) = tk_mtrx;
    A.at(1) = yk;
    return A;
}

}

namespace Approximation
{
double Deflection (const Matrix<double>& X, const Matrix<double>& Y, const Matrix<double>& A, const int i)
{

    double S = 0;
    double S2 = 0;
    for (auto pX = X.last_i(), pY = Y.last_i(); pY >= Y.get_pointer(); --pX, --pY)
    {
        S = *pY;
        int k = i;
        for (auto pA = A.last_i(); pA >= A.get_pointer(); --pA)
        {
            if (pA == A.get_pointer())
                assert (k == 0);
            S -= *pA*Pow(*pX, k--);
        }
        S2 += S*S;
    }
    return S2;
}
Matrix<double> Find_Polinom_m_power (const Matrix<double>& X, const Matrix<double>& Y, double m)
{
    assert (X.get_size() == Y.get_size() && ((X.get_n() == 1 && Y.get_n() == 1) || (X.get_m() == 1 && Y.get_m() == 1)));
    if (m < DBL_EPSILON || m > X.get_size())
        m = X.get_size();
    const unsigned M = m;
    unsigned int i = 2*M;
    double* Sums = new double [2*M + 1];
    Sums[0] = 1;
    for (double* p = Sums + 2*M; p > Sums; --p)
    {
        m = 0;
        for(const double* pArr = X.get_pa() + X.get_size() - 1; pArr >= X.get_pa(); --pArr)
        {
            m += Pow(*pArr, i);
        }
        --i;
        *p = m;
    }
    Matrix_SLE Temp (M+1, M+1);
    for (double* p = Temp.last_i(); p >= Temp.get_pa(); p -= Temp.get_m())
    {
        double* ps = Sums + 2*M - i++;
        for (double* pln = p; /*ps > Sums + M - 2 - i*/ pln > p - M - 1; --pln, --ps)
            *pln = *ps;
    }
    i = M;
    Matrix<double> ySums (M+1, 1);
    for (double* pS = ySums.get_pa() + M; pS > ySums.get_pa(); --pS)
    {
        m = 0;
        for(const double* pArr = X.get_pa() + X.get_size() - 1, *py = Y.get_pa() + Y.get_size() - 1; pArr >= X.get_pa(); --pArr, --py)
        {
            m += Pow(*pArr, i)**py;
        }
        --i;
        *pS = m;
    }
    delete[] Sums;
    for (const double *py = Y.get_pa() + Y.get_size() - 1; py >= Y.get_pa(); --py)
    {
        *ySums.get_pa() += *py;
    }
    Matrix<double> Params = Temp.Solve(ySums);
    return Params;
}
Matrix<double> Find_Polinom (const Matrix<double> X, const Matrix<double> Y)
{
    assert (X.get_size() == Y.get_size() && ((X.get_n() == 1 && Y.get_n() == 1) || (X.get_m() == 1 && Y.get_m() == 1)));
    Matrix<double> A;
    double S2 = 0;
    double* pS = new double [X.get_size()];
    double* pSt = pS;
    for (int i = 1; i <= X.get_size(); ++i)
    {
        A = Find_Polinom_m_power(X, Y, i);
        S2 = Deflection(X, Y, A, i);
        *(pSt++) = S2/(X.get_size() - i - 1);
    }
    View(pS, X.get_size(), 1);
    return Find_Polinom_m_power(X, Y, Find_min_by_abs(pS, X.get_size()) + 1);
}
}
#endif // MATH_H
