#include <stdlib.h>
#include <string.h>
#include <string>
#include <math.h>
#include <fstream>
#include <iostream>
#include <algorithm>
#define PRES double
#define NXB 15
#define NYB 12
#define NX NXB*3+1
#define NY NYB*3+1
#define NYK2 NYB*2
#define REP 3000
#define DEL 100
#define AMAT 1.1f
#define TEM1 5.0f
#define TEM2 15.0f
#define HX 0.2f
#define HY 0.3f
using namespace std;
int main(int argc, char **argv)
{
    int i, j, k;
    int idt=0;
    int ndt=0;
    int i1 = NXB;
    int i2 = NXB*2;
    int i3 = NXB*3;
    int j1 = NYB;
    int j2 = NYB*2;
    int j3 = NYB*3;
    double T1=TEM1, T2=TEM2, h=HX, r=HY, a=AMAT, t0;
    double T[NY][NX] = {0};
    double TT[NY][NX] = {0};
    double rr = min(h, r);
    double tau = 0.25f*rr*rr/a;
    double alf_1 = -h/r;
    double alf_2 = -r/h;
    double alf_3 = 0.5f * alf_2;
    double alf_4 = 0.5f * alf_1;
    double bet_1 = a*tau / (h*r);
    double bet_2 = 2.0f*bet_1;
    double bet_3 = 4.0f*bet_1;
    double gam_1 = -2.f*(alf_1 + alf_2);
    double gam_3 = - (alf_1 + alf_2);
    double gam_4 = - (alf_3 + alf_4);
    //char filename[128];
    for (j=0; j<NY; j++)
    {
        for (i=0; i<NX; i++)
        {
            T[j][i] = 0.0f;
        }
    }
    for (j=0; j<=NYB; j++)
    {
        T[j][0] = T1;
        TT[j][0] = T1;
    }
    for (j=NYK2; j<NY; j++)
    {
        T[j][NX-1] = T2;
        TT[j][NX-1] = T2;
    }
    ofstream fout("T1.dat",ios_base::out | ios_base::trunc | ios_base::binary);
    for (j = 0; j < NY; j++)
    {
        for (i = 0; i < NX; i++)
        {
            double w = T[j][i];
            fout.write((char*)&w, sizeof w);
        }
    }
    fout.close();
    for (k=0; k<REP; k++)
    {
        for (j=0; j<NY; j++)
        {
            for (i=0; i<NX; i++)
            {
                t0 = T[j][i];
                if (i>0 && i<NX-1 && j==0)
                {
                    TT[j][i] = t0 - bet_2*(alf_3*(T[j][i-1]+T[j][i+1]) + alf_1*T[j+1][i] +gam_3*t0);
                }
                if (i>0 && i<NX-1 && j==NY-1)
                {
                    TT[j][i] = t0 - bet_2*(alf_3*(T[j][i-1]+T[j][i+1]) + alf_1*T[j-1][i] +gam_3*t0);
                }
                if (i==0 && j>NYB && j<NY-1)
                {
                    TT[j][i] = t0 - bet_2*(alf_4*(T[j-1][i]+T[j+1][i]) + alf_2*T[j][i+1] +gam_3*t0);
                }
                if (i==NX-1 && j>0 && j<NYK2)
                {
                    TT[j][i] = t0 - bet_2*(alf_4*(T[j-1][i]+T[j+1][i]) + alf_2*T[j][i-1] +gam_3*t0);
                }
                if (i==0 && j==NY-1)
                {
                    TT[j][i] = t0 - bet_3*(alf_3*T[j][i+1] + alf_4*T[j-1][i] + gam_4*t0);
                }
                if (i==NX-1 && j==0)
                {
                    TT[j][i] = t0 - bet_3*(alf_3*T[j][i-1] + alf_4*T[j+1][i] + gam_4*t0);
                }
                if (i>0 && i<NX-1 && j>0 && j<NY-1)
                {
                    TT[j][i] = t0 - bet_1*(alf_1*(T[j][i-1]+T[j][i+1]) + alf_2*(T[j-1][i]+T[j+1][i]) + gam_1*t0);
                }
            }
        }
        for (j=0; j<NY; j++)
        {
            for (i=0; i<NX; i++)
            {
                T[j][i] = TT[j][i];
            }
        }
        idt++;
        if (idt == DEL)
        {
            idt=0;
            ndt++;
            string Name;
            for (int i = ndt+1; i >= 1; i /= 10)
            {
                Name += '0' + i%10;
            }
            Name += 'T';
            reverse (Name.begin(), Name.end());
            Name += ".dat";
//            sprintf_s(filename, sizeof(filename), "T%d.dat", ndt+1);
            ofstream fout(Name.c_str(), ios_base::out | ios_base::trunc | ios_base::binary);
            for (j = 0; j < NY; j++)
            {
                for (i = 0; i < NX; i++)
                {
                    double w = T[j][i];
                    fout.write((char*)&w, sizeof w);
                }
            }
            fout.close();
        }
    }
    for (int i = 0; i <= i3; ++i)
    {
        for (int j = 0; j < j3; ++j)
            cout<<TT[j][i]<<' ';
        cout<<TT[j3][i];
        cout<<endl;
    }
    int n_x = NX;
    int n_y = NY;
    int n_k = ndt;
    ofstream fou("Param.dat",ios_base::out | ios_base::trunc | ios_base::binary);
    fou.write((char*)&n_x, sizeof n_x);
    fou.write((char*)&n_y, sizeof n_y);
    fou.write((char*)&n_k, sizeof n_y);
    fou.close();
}
