// This is a personal academic project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: http://www.viva64.com

#include <string>
#include <fstream>
#include <iostream>
#include <boost/test/unit_test.hpp>
#include <filesystem>
extern std::filesystem::path path;

//#define PRES double
#define NXB 15
#define NYB 12
#define NX NXB*3+1
#define NY NYB*3+1
//#define NYK2 NYB*2
#define REP_Lab_7 4000
#define DEL 100
#define AMAT 1.1f
#define TEM1 5.0f
#define TEM2 15.0f
#define HX 0.2f
#define HY 0.3f
using namespace std;
int main_for_Lab_7();

BOOST_AUTO_TEST_SUITE(SuItE_tests_for_Lab_7)
    BOOST_AUTO_TEST_CASE(Case_for_lab_7)
    {
        BOOST_CHECK(main_for_Lab_7()==0);
    }
BOOST_AUTO_TEST_SUITE_END()

int main_for_Lab_7 ()
{
    int i, j, k;
    int idt=0;
    int ndt=1;
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
    double bet_2 = 4.0*bet_1/3;
    double bet_3 = 2.0f*bet_1;
    double bet_4 = 4.0f*bet_1;
    double gam_1 = -2.f*(alf_1 + alf_2);
    double gam_2 = -1.5*(alf_1 + alf_2);
    double gam_3 = - (alf_1 + alf_2);
    double gam_4 = - (alf_3 + alf_4);
    double w = 0;
    string Name;
    for (j=0; j<=j3; j++)
    {
        T[j][0] = T1;
        TT[j][0] = T1;
        T[j][i3] = T2;
        TT[j][i3] = T2;
    }
    for (i = 0; i <= i1; ++i)
    {
        T[0][i] = T1;
        TT[0][i] = T1;
        T[j3][i] = T1;
        TT[j3][i] = T1;
    }
    for (i = i2; i <= i3; ++i)
    {
        T[0][i] = T1;
        TT[0][i] = T1;
        T[j3][i] = T1;
        TT[j3][i] = T1;
    }
    ofstream fout(path + "T1.dat", ios_base::out | ios_base::trunc | ios_base::binary);
    for (j = 0; j < NY; j++)
    {
        for (i = 0; i < NX; i++)
        {
            w = T[j][i];
            fout.write(reinterpret_cast<const char*>(&w), sizeof w);
        }
    }
    fout.close();
    for (k=0; k < REP_Lab_7; k++)
    {
        for (j=0; j<NY; j++)
        {
            for (i=0; i<NX; i++)
            {
                t0 = T[j][i];
                if (i>i1 && i<i2 && j==0)
                    TT[j][i] = t0 - bet_3*(alf_3*(T[j][i-1]+T[j][i+1]) + alf_1*T[j+1][i] +gam_3*t0);
                if (i>i1 && i<i2 && j==NY-1)
                    TT[j][i] = t0 - bet_3*(alf_3*(T[j][i-1]+T[j][i+1]) + alf_1*T[j-1][i] +gam_3*t0);
                if (i == i1 && j < j1)
                    TT[j][i] = t0 - bet_3*(alf_4*(T[j-1][i]+T[j+1][i]) + alf_2*T[j][i-1] + gam_3*t0);
                if (i == i1 && j > j2 && j < j3)
                    TT[j][i] = t0 - bet_3*(alf_4*(T[j-1][i]+T[j+1][i]) + alf_2*T[j][i-1] + gam_3*t0);
                if (i == i2 && j < j1)
                    TT[j][i] = t0 - bet_3*(alf_4*(T[j-1][i]+T[j+1][i]) + alf_2*T[j][i+1] + gam_3*t0);
                if (i == i2 && j > j2 && j < j3)
                    TT[j][i] = t0 - bet_3*(alf_4*(T[j-1][i]+T[j+1][i]) + alf_2*T[j][i+1] + gam_3*t0);
                if (i > i1 && i < i2 && j == j1)
                    TT[j][i] = t0 - bet_3*(alf_3*(T[j][i-1] + T[j][i+1]) + alf_1*T[j+1][i] + gam_3*t0);
                if (i > i1 && i < i2 && j == j2)
                    TT[j][i] = t0 - bet_3*(alf_3*(T[j][i-1] + T[j][i+1]) + alf_1*T[j-1][i] + gam_3*t0);

                if (i == i1 && j == j1)
                    TT[j][i] = t0 - bet_2*(alf_4*T[j-1][i] + alf_2*T[j][i-1] + alf_3*T[j][i+1] + alf_1*T[j+1][i] + gam_2*t0);
                if (i == i1 && j == j2)
                    TT[j][i] = t0 - bet_2*(alf_1*T[j-1][i] + alf_2*T[j][i-1] + alf_3*T[j][i+1] + alf_4*T[j+1][i] + gam_2*t0);
                if (i == i2 && j == j1)
                    TT[j][i] = t0 - bet_2*(alf_4*T[j-1][i] + alf_3*T[j][i-1] + alf_2*T[j][i+1] + alf_1*T[j+1][i] + gam_2*t0);
                if (i == i2 && j == j2)
                    TT[j][i] = t0 - bet_2*(alf_1*T[j-1][i] + alf_3*T[j][i-1] + alf_2*T[j][i+1] + alf_4*T[j+1][i] + gam_2*t0);

                if (i < i1 && j < j3)
                    TT[j][i] = t0 - bet_1*(alf_1*(T[j][i-1]+T[j][i+1]) + alf_2*(T[j-1][i]+T[j+1][i]) + gam_1*t0);
                if (i > i2 && i < i3 && j < j3)
                    TT[j][i] = t0 - bet_1*(alf_1*(T[j][i-1]+T[j][i+1]) + alf_2*(T[j-1][i]+T[j+1][i]) + gam_1*t0);
                if (i>=i1 && i <=i2 && j > j1 && j < j2)
                    TT[j][i] = t0 - bet_1*(alf_1*(T[j][i-1]+T[j][i+1]) + alf_2*(T[j-1][i]+T[j+1][i]) + gam_1*t0);
            }
        }
        for (j=0; j < NY; j++)
        {
            for (i=0; i < NX; i++)
            {
                T[j][i] = TT[j][i];
            }
        }
        idt++;
        if (idt == DEL)
        {
            idt=0;
            Name.clear();
            for (int is = ++ndt; is >= 1; is /= 10)
            {
                Name += '0' + is%10;
            }
            Name += 'T';
            reverse (Name.begin(), Name.end());
            Name += ".dat";
            ofstream foutn(path + Name, ios_base::out | ios_base::trunc | ios_base::binary);
            for (j = 0; j < NY; j++)
            {
                for (i = 0; i < NX; i++)
                {
                    w = T[j][i];
                    foutn.write(reinterpret_cast<const char*>(&w), sizeof w);
                }
            }
            foutn.close();
        }
//        for (int j = 0; j < j3; ++j)
//        {
//            for (int i = 0; i <= i3; ++i)
//                cout<<TT[j][i]<<' ';
//            cout<<TT[j][i3];
//            cout<<endl;
//        }
//        cout<<endl<<endl<<endl;
    }
    int n_x = NX;
    int n_y = NY;
    int n_k = --ndt;
    ofstream fou(path + "Param.dat", ios_base::out | ios_base::trunc | ios_base::binary);
    fou.write(reinterpret_cast<const char*>(&n_x), sizeof n_x);
    fou.write(reinterpret_cast<const char*>(&n_y), sizeof n_y);
    fou.write(reinterpret_cast<const char*>(&n_k), sizeof n_k);
    fou.close();
    return 0;
}
