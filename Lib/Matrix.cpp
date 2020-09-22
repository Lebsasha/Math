// This is a personal academic project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: http://www.viva64.com

/*double** Set1 (const int n)
	{
        double **M = NULL;
        M = new double* [n];
        for (int i = 0; i < n; ++i)
        {
			*(M++) = new double [n];
		}
        M -= n;
        int i = 0;
        for (int i = 0; i < n; ++i)
			for (int j = 0; j < n; ++j)
				*(*(M+i)+j) = 0;
		return M;
	}
double** Multi (double** M1, double** M2, const int n)
	{
		double** Mt = Set1 (n);
		for (int i = 0; i < n; ++i)
		{
			for (int k = 0; k < n; ++k)
			{
				for (int j = 0; j < n; ++j)
				{
					*(*(Mt+i)+k) += *(*(M1+i)+j)* *(*(M2+j)+k);
				}
			}
		}
		return Mt;
	}
double** Summ(double** M1, double** M2, const int n)
	{
		double** Mt = Set1 (n);
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < n; ++j)
		{
			*(*(Mt+i)+j) = *(*(M1+i)+j) + *(*(M2+i)+j);
		}
		return Mt;
	}
double** set_rand1(const int n, int max_e)
	{
        double** Temp = Set1 (n);
		max_e++;
		double* pTemp = *Temp;
		int co = 0;
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < n; ++j)
		*(*(Temp+i)+j) = rand()%max_e-rand()%max_e;
		return Temp;
	}
double* Set2 (const int n)
	{
        double *M = new double [n*n];
		for (int j = 0; j < n*n; ++j)
			*M++ = 0;
		M -= n*n;
		return M;
	}
double* set_rand2(const int n, int max_e)
	{
		double* Temp = Set2 (n);
		for (int i = 0; i < n*n; ++i)
			*Temp++ = rand()%max_e-rand()%max_e;
		Temp -= n*n;
		return Temp;
	}
double* Multi (double* M1, double* M2, const int n)
	{
		double* Mt = Set2 (n);
		for (int i = 0; i < n; ++i)
		{
			for (int k = 0; k < n; ++k)
			{
				for (int j = 0; j < n; ++j)
				{
					*(Mt+i*n+k) += *(M1+i*n+j)**(M2+j*n+k);
				}
			}
		}
		return Mt;
	}
double* Summ(double* M1, double* M2, const int n)
	{
		double* Mt = Set2 (n);
		for (int i = 0; i < n*n; ++i)
		{
			*Mt++ = *(M1++) + *(M2++);
		}
		Mt -= n*n;
		return Mt;
	}*/

