// This is a personal academic project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: http://www.viva64.com

#include <fstream>
#include <boost/test/unit_test.hpp>
extern std::string path;

#define  PRES   double
#define  NXB     15
#define  NX        NXB*3+1
#define  NYB     12
#define  NY        NYB*3+1
#define  REP_Lab_6     3000
#define  EPSL   1.e-5
#define  LL        1.7
#define  TEM1   5.0f
#define  TEM2  15.0f
#define  HX    0.2f
#define  HY    0.3f

    using namespace std;

    void maxpvr(double *t1, double *del, double *maxdel) {
        double d = fabs(*del) / fabs(*t1);
        if (d > *maxdel) *maxdel = d;
    }

    int main_for_Lab_6();

    BOOST_AUTO_TEST_SUITE(SuItE_tests_for_Lab_6)

        BOOST_AUTO_TEST_CASE(Case_for_lab_6)
        {
            BOOST_CHECK(main_for_Lab_6() == 0);
        }

    BOOST_AUTO_TEST_SUITE_END()

    int main_for_Lab_6() {
        ofstream foutT(path + "dT.dat", ios_base::out | ios_base::trunc | ios_base::binary);
        int i1, i2, i3, j1, j2, j3, rp, i, j, k = 0;
        double T1 = TEM1, T2 = TEM2, h = HX, r = HY, tx, t0, t1, del, maxdel = 0.0f;
        double T[NY][NX];
        double lam = LL;
        double eps = EPSL;
        int prz = 1;
        int nT = 0;
        double alf_1 = -h / r;
        double alf_2 = -r / h;
        double alf_3 = alf_2 * 0.5f;
        double alf_4 = alf_1 * 0.5f;
        double gam_1 = -2.f * (alf_1 + alf_2);
        double gam_2 = -1.5f * (alf_1 + alf_2);
        double gam_3 = -(alf_1 + alf_2);
        double gam_4 = -(alf_3 + alf_4);
        i1 = NXB;
        i2 = i1 + NXB;
        i3 = NXB * 3;
        j1 = NYB;
        j2 = NYB * 2;
        j3 = NYB * 3;
        rp = REP_Lab_6;
        for (j = 0; j <= j3; j++) {
            for (i = 0; i <= i3; i++) {
                T[j][i] = 0.0f;
            }
        }
        for (j = 0; j <= j3; ++j)
            T[j][0] = T1;
        for (i = 0; i <= i1; ++i) {
            T[0][i] = T1;
            T[j3][i] = T1;
        }
//    for (i=i2; i<=i3; ++i)
//        T[j2][i] = T2;
        for (j = 0; j <= j3; ++j) {
            T[j][i3] = T2;
        }
        for (i = i2; i <= i3; ++i) {
            T[0][i] = T2;
            T[j3][i] = T2;
        }
        while (k < rp && prz == 1) {
            k++;
            for (j = 0; j < j3; j++)//???<
            {
                for (i = 1; i <= i3; i++) {
                    t0 = T[j][i];
                    if (i == i1 && ((j > 0 && j < j1) || (j > j2 && j < j3))) {
                        tx = -(alf_4 * (T[j - 1][i] + T[j + 1][i]) + alf_2 * T[j][i - 1]) / gam_3;
                        del = lam * (tx - t0);
                        t1 = t0 + del;
                        T[j][i] = t1;
                        maxpvr(&t1, &del, &maxdel);
                    }
                    if (i == i2 && ((j > 0 && j < j1) || (j > j2 && j < j3))) {
                        tx = -(alf_4 * (T[j - 1][i] + T[j + 1][i]) + alf_2 * T[j][i + 1]) / gam_3;
                        del = lam * (tx - t0);
                        t1 = t0 + del;
                        T[j][i] = t1;
                        maxpvr(&t1, &del, &maxdel);
                    }
                    if (i > i1 && i < i2 && j == j1) {
                        tx = -(alf_3 * (T[j][i - 1] + T[j][i + 1]) + alf_1 * T[j + 1][i]) / gam_3;
                        del = lam * (tx - t0);
                        t1 = t0 + del;
                        T[j][i] = t1;
                        maxpvr(&t1, &del, &maxdel);
                    }
                    if (i > i1 && i < i2 && j == j2) {
                        tx = -(alf_3 * (T[j][i - 1] + T[j][i + 1]) + alf_1 * T[j - 1][i]) / gam_3;
                        del = lam * (tx - t0);
                        t1 = t0 + del;
                        T[j][i] = t1;
                        maxpvr(&t1, &del, &maxdel);
                    }
                    if (i == i1 && j == j1) {
                        tx = -(alf_4 * T[j - 1][i] + alf_2 * T[j][i - 1] + alf_3 * T[j][i + 1] + alf_1 * T[j + 1][i]) /
                             gam_2;
                        del = lam * (tx - t0);
                        t1 = t0 + del;
                        T[j][i] = t1;
                        maxpvr(&t1, &del, &maxdel);
                    }
                    if (i == i1 && j == j2) {
                        tx = -(alf_1 * T[j - 1][i] + alf_2 * T[j][i - 1] + alf_3 * T[j][i + 1] + alf_4 * T[j + 1][i]) /
                             gam_2;
                        del = lam * (tx - t0);
                        t1 = t0 + del;
                        T[j][i] = t1;
                        maxpvr(&t1, &del, &maxdel);
                    }
                    if (i == i2 && j == j1) {
                        tx = -(alf_4 * T[j - 1][i] + alf_3 * T[j][i - 1] + alf_2 * T[j][i + 1] + alf_1 * T[j + 1][i]) /
                             gam_2;
                        del = lam * (tx - t0);
                        t1 = t0 + del;
                        T[j][i] = t1;
                        maxpvr(&t1, &del, &maxdel);
                    }
                    if (i == i2 && j == j2) {
                        tx = -(alf_1 * T[j - 1][i] + alf_3 * T[j][i - 1] + alf_2 * T[j][i + 1] + alf_4 * T[j + 1][i]) /
                             gam_2;
                        del = lam * (tx - t0);
                        t1 = t0 + del;
                        T[j][i] = t1;
                        maxpvr(&t1, &del, &maxdel);
                    }
                    if (i > 0 && i < i1 && j > 0 && j < j3) {
                        tx = -(alf_1 * (T[j - 1][i] + T[j + 1][i]) + alf_2 * (T[j][i - 1] + T[j][i + 1])) / gam_1;
                        del = lam * (tx - t0);
                        t1 = t0 + del;
                        T[j][i] = t1;
                        maxpvr(&t1, &del, &maxdel);
                    }
                    if (i >= i1 && i <= i2 && j > j1 && j < j2) {
                        tx = -(alf_1 * (T[j - 1][i] + T[j + 1][i]) + alf_2 * (T[j][i - 1] + T[j][i + 1])) / gam_1;
                        del = lam * (tx - t0);
                        t1 = t0 + del;
                        T[j][i] = t1;
                        maxpvr(&t1, &del, &maxdel);
                    }
                    if (i > i2 && i < i3 && j > 0 && j < j3) {
                        tx = -(alf_1 * (T[j - 1][i] + T[j + 1][i]) + alf_2 * (T[j][i - 1] + T[j][i + 1])) / gam_1;
                        del = lam * (tx - t0);
                        t1 = t0 + del;
                        T[j][i] = t1;
                        maxpvr(&t1, &del, &maxdel);
                    }
                }
            }
            nT++;
            double w = maxdel;
            foutT.write((char *) &w, sizeof w);
            if (maxdel < eps) prz = 0;
            maxdel = 0.0f;
        }
        foutT.close();
        ofstream fouT(path + "nT.dat", ios_base::out | ios_base::trunc | ios_base::binary);
        fouT.write((char *) &nT, sizeof nT);
        fouT.close(); // ��������� ����
        ofstream fout(path + "Pole.dat", ios_base::out | ios_base::trunc | ios_base::binary);
        for (j = 0; j < NY; j++) {
            for (i = 0; i < NX; i++) {
                double w = T[j][i];
                fout.write((char *) &w, sizeof w);
            }
        }
        fout.close();
        int n_x = NX;
        int n_y = NY;
        ofstream fou(path + "Param.dat", ios_base::out | ios_base::trunc | ios_base::binary);
        fou.write((char *) &n_x, sizeof n_x);
        fou.write((char *) &n_y, sizeof n_y);
        fou.close();
        return 0;
    }
