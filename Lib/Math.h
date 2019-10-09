#include "../Lib/Matrix.h"
class Matrix_SE: public Matrix<double>
{
    bool High_accuracy;
protected:
    char Type;
public:
    Matrix_SE  (): Matrix<double>(), High_accuracy(1), Type('N')
    {}
//    Matrix_SE (const Matrix<double>& A): Matrix<double>(A), High_accuracy(1), Type('N')
//    {}
    Matrix_SE (const Matrix<double>& A): Matrix<double>(A), High_accuracy(1), Type('N')
    {}
    Matrix_SE (double* pA, const int a, const int b): Matrix<double>(pA, a, b),  High_accuracy(1), Type('N')
    {}
    Matrix_SE (Matrix_SE& A): Matrix<double>(*(reinterpret_cast<Matrix<double>*>(&A))), High_accuracy (A.High_accuracy), Type('N')
    {}
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
        double* pb = pa + N_R*M;
        double* pendb = pb + M;
        for (double* pt = pb; pt < pendb; ++pt)
        {
            *pt = C*(*pt);
        }
        return;
    }
    void One_Minus_Second_Row (const int one, const int second, const double C)
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
    virtual const Matrix_SE& View (void) const
    {
        //cout<<Type<<' ';/////////////////////////////////////////////////////
        Matrix::View();
        return *this;
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
        delete pa;
        if (sa.eof())
            return 0;
        sa.read(&Type, sizeof(char));
        return Matrix::Read_from_file(sa);
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
        sa.write (&Type, sizeof(char));
        return Matrix::Write_to_file(sa);;
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
//    Matrix_SLE (const Matrix<double>& A): Matrix_SE(A)
//    {
//        Type = 'L';
//    }
    Matrix_SLE (Matrix<double> A): Matrix_SE(A)
    {
        Type = 'L';
    }
    Matrix_SLE (double* pA, const int a, const int b): Matrix_SE(pA, a, b)
    {
        Type = 'L';
    }
//    Matrix_SLE (const Matrix_SE& A): Matrix_SE(A)
//    {
//        Type = 'L';
//    }
    Matrix_SLE (Matrix_SE A): Matrix_SE(A)
    {
        Type = 'L';
    }
    Matrix_SLE (Matrix_SLE& A): Matrix_SE(*(reinterpret_cast<Matrix_SE*>(&A)))
    {
        Type = 'L';
    }
    Matrix_SLE (const int N1, const int M1): Matrix_SE(N1, M1)
    {
        Type = 'L';
    }
    bool If_Symmetric (void) const
    {
        bool a = 1;
        return !a;
    }
    Matrix_SLE T (void) const
    {
        Matrix_SLE B;
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
    Matrix<double> Solve (const Matrix<double>& b) const
    {
        return (If_Symmetric()) ? LDLT(b) : Solve_by_Gauss_Method(b);
    }
    Matrix<double> Solve_by_Gauss_Method (Matrix_SE C) const
    {
        Matrix_SLE S (*(const_cast<Matrix_SLE*>(this)));
        int k = 0;
        double T = 0;
        for (int i = 0; i < N - 1; i++)
        {
            k = S.Find_Max_El_in_Col(i);
            if (fabs(S.Get_El(k*M+k)) < 10*DBL_EPSILON)
            {
                Matrix<double> A;
                return Matrix<double> ();
            }
            if (i != k)
            {
                S.Replace_Rows(i, k);
                C.Replace_Rows(i, k);
            }
            T = S.Get_El(i*M+i);
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
        //double* pEnd = ps + N;
        //for (double* pi = 0; pi < pEnd; )
        const double* pMatrix = S.Get_pa() + S.Get_N()*S.Get_M() - 1;
        const double* pC = C.Get_pa() + C.Get_N()*1 - 1;
        double* pBeg = ps.Get_pa() + ps.Get_N()*1 - 1; //const_?
        double* pEnd = ps.Get_pa();
        k = N - 1; //////////////////////////////////////////////////
        k = M;
        for (double* pt = pBeg; pt >= pEnd; --pt, pMatrix -= (k--))
        {
            *pt = *(pC--);
            for (double* pc = pBeg; pc > pt; --pMatrix, --pc)
                *pt = *pt - *pMatrix*(*pc);
            *pt = *pt/(*pMatrix);
        }
        //Matrix<double> D (ps.Get_pa(), N, 1);
        return ps;
    }
    Matrix<double> LDLT (Matrix<double> b) const
    {
        return Matrix<double> ();
    }
    vector<double> Solve_by_Gauss_Method_v (vector<double>& C)
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
    virtual const Matrix_SLE& View (void)
    {
        Matrix_SE::View();
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
    double Integral_by_definition (const double from, const double to) const;
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
        double* pA = A.Get_pa() + N - 1;
        for (Function* ppa = this->pa + this->N - 1; ppa >= this->pa; --(ppa), --pA)
        {
            *pA = ppa->f(x);
        }
        return A;
    }
//    Matrix_SLE Derivative (const Matrix<double> X) const
//    {
//        int Num_of_variables = X.Get_Size();
//        Matrix<double> A (N, Num_of_variables);
//        double* pa = A.Get_pa() + A.Get_Size() - 1;
//        Matrix<double> Temp (Num_of_variables, 1);
//        double* pTemp = 0;
//        for (Function2* ppa = pa + N*Num_of_variables - 1; ppa >= pa; --ppa, --pa)
//        {
//            Temp = ppa->Derivative_by_definition_M_in_column(X);
//            pTemp = Temp.Get_pa() + Num_of_variables - 1;
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
        for (Function* pA = pa + N - 1; pA >= pa; --pA)
        {
            if (m < (v = pA->f(x)))
                m = v;
        }
        return m;
    }
//    ~Matrix_SNE (void)
//    {
//        delete[] pa;
//    }
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
//    double Derivative_by_user_def (double (*Der) (double, double), double x,)
//    {
//        return Der(x);
//    }
    Matrix_SLE Derivative_by_definition_M_in_column (const Matrix<double>& X) const
    {
        int Num_of_variables = X.Get_Size();
        Matrix<double> A (Num_of_variables, 1);
        double* pa = A.Get_pa();
        const double M = 0.01; //0.05   0.1/////////////////////////////////////////////////////////////////////////////
        double F = Fn(X);
        Matrix<double> DX(X);
        double* pEnd = DX.Get_pa() + DX.Get_Size();
        for (double* pDX = DX.Get_pa(), *pX = X.Get_pa(); pDX < pEnd; ++pDX, ++pX, ++pa)
        {
            *pDX += M*(*pDX);
            *pa = (Fn(DX) - F)/(M*(*pX));
            *pDX = *pX;
        }
        //*pa-- = (Fn(x + M*x, y) - Fn(x, y))/M*;
        //*pa = (Fn(x, y + M*y) - Fn(x, y))/Epsy;
        return A;
    }
    double f(const Matrix<double>& X) const
    {
        return Fn(X);
    }
    double Integral_by_user_def (double (*AntiDer) (const Matrix<double>&), const double from_x, const double from_y, const double to_x, const double to_y) const;
    double Integral_by_definition (const double from_x, const double from_y, const double to_x, const double to_y) const;
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
//    Array_of_Functions2 (Array_of_Functions2& Ar): N(Ar.N), Fn (new Function2[N])
//    {
//        for (Function2* pFn = Fn + N - 1, *pFa = Ar.Fn +  Ar.N - 1; pFn >= Fn; --pFn, --pFa)
//        {
//            *pFn = *pFa;
//        }
//    }
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
    int Get_N (void) const
    {
        return N;
    }
    Function2* Get_Fn (void) const
    {
        return Fn;
    }
    Function2& operator[] (const int S) const
    {
        try
        {
            if (S > -1 && S < N)
            {
                return *(Fn + S);
            }
            throw 1;
        }
        catch (int i)
        {
            cout<<"Warning, Function2 haven't"<<S<<"element";
            return *Fn;
        }
    }
    Matrix<double> f (const Matrix<double>& X) const
    {
        Matrix<double> A (N, 1);
        double* pa = A.Get_pa() + N - 1;
        for (Function2* pFn = Fn + N - 1; pFn >= Fn; --pFn, --pa)
        {
            *pa = pFn->f(X);
        }
        return A;
    }
    Matrix<double> Derivative (const Matrix<double>& X) const
    {
        const int Num_of_variables = X.Get_Size();
        Matrix<double> A (N, Num_of_variables);
        double* pa = A.Get_pa() + A.Get_Size() - 1;
        Matrix<double> Temp (Num_of_variables, 1);
        double* pTemp = 0;
        for (Function2* pFn = Fn + N - 1; pFn >= Fn; --pFn)
        {
            Temp = pFn->Derivative_by_definition_M_in_column(X);
            pTemp = Temp.Get_pa() + Num_of_variables - 1;
            for (int i = 0; i < Num_of_variables; ++i, --pa, --pTemp)
            {
                *pa = *pTemp;
            }
        }
        return A;
    }
    Matrix<double> Integral (const double from_x, const double from_y, const double to_x, const double to_y) const;
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

double fDelta2 (const Matrix<double>& X_i, const Matrix<double>& Delta)
{
    double Delta2 = 0;
    double Temp = 0;
    for (double* x = X_i.Get_pa() + X_i.Get_Size() - 1, *dx = Delta.Get_pa() + Delta.Get_Size() - 1; x >= X_i.Get_pa(); --x, --dx)
    {
        if (*x + *dx < 1 && *dx > Delta2)
        {
            Delta2 = *dx;
        }
        else if ((Temp = *dx/(*x+*dx)) > Delta2)
        {
            Delta2 = Temp;
        }
    }
    return Delta2;
}
Matrix<double> Solve_SNE (const Array_of_Functions2& Func, const Array_of_Functions2& Der_of_Arr, const Matrix<double>& xy, const double Eps1, const double Eps2)
{
    const int Num_it = 10000;
    int i = 0;
    Matrix<double> ixy (xy);
    Matrix<double> Delta(xy.Get_N(), xy.Get_M());
    Matrix<double> High_Precision_ixy (xy);
    Matrix<double> High_Precision_Delta(xy.Get_N(), xy.Get_M());
    Matrix_SLE Derivative;
    Matrix_SLE High_Precision_Derivative (xy.Get_N(), xy.Get_N()); //It is not a mistake!
    Matrix<double> Sol;
    double Delta1 = 0;
    double Delta2 = 0;
    do
    {
        cout<<setw(10)<<Delta1<<setw(10)<<Delta2<<setw(6)<<++i<<endl;
//        cout<<"ixy"<<endl;
//        ixy.View();
//        cout<<"Delta"<<endl;
//        Delta.View();
        ixy = ixy + Delta;
//        cout<<"Sol"<<endl;
        Sol = Func.f(ixy).Multiple_Matrix_by_Number(-1).View();
        Delta1 = Sol[0];
//        cout<<"Derivative"<<endl;
        Derivative = Func.Derivative(ixy);
        Delta = Derivative.Solve(Sol);
//        cout<<"HP_ixy"<<endl;
//        High_Precision_ixy.View();
//        cout<<"HP_Delta"<<endl;
//        High_Precision_Delta.View();
        High_Precision_ixy = High_Precision_ixy + High_Precision_Delta;
        Function2* pF = Der_of_Arr.Get_Fn() + Der_of_Arr.Get_N() - 1;
        for (double* pD = High_Precision_Derivative.Get_pa() + High_Precision_Derivative.Get_Size() - 1; pD >= High_Precision_Derivative.Get_pa(); --pD, --pF)
        {
            *pD = -pF->f(High_Precision_ixy);// Minus against write Sol = Func.f(High_Precision_ixy)
        }                         //.Multiple_Matrix_by_Number(-1).View();
        cout<<"With user defined derivative F."<<endl;
        Sol = Func.f(High_Precision_ixy).View();
//        cout<<"HP_Derivative"<<endl;
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
