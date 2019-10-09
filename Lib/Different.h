#ifndef DIFFERENT_H
#define DIFFERENT_H
void randoe (int*, const int, const int, int);
void randoe (double*, const int, const int, int);
void randoe (float*, const int, const int, int);
void View (int*, const int, const int);
void View (double*, const int, const int);
void View (float*, const int, const int);
void View (double**, const int, int m = 0);
void View (float**, const int, int m = 0);
void View (int**, const int, int m = 0);
#include <vector>
//using namespace std;
void View (std::vector<double>& El1);
void Nul (double*, const int);
void Nul (int*, const int);
void Nul (float*, const int);
void Nul (char*, const int);
bool if_Simple (int);
void Get_Pause (void);
#endif
