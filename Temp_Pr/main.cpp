#include <iostream>
#define NDEBUG 1 //deactivate assert()
#include "../Lib/Math.h"
#include <cassert>
#include <time.h>
#include <stdlib.h>
#include <fstream>
using namespace std;
int main()
{
    const int N = 10;
    int pt[] = {3, 8, 1, 1, 0, 5, 3, 6, 6, 9};
    Matrix<int> A (500, 1000000);
    srand(time(nullptr));
    A.Rand(N);
    int pa[N+1] = {0};
    int S = 0;
    for (auto i = A.Last_i(); i >= A.First_i(); --i)
    {
        ++pa[*i];
    }
    View(pa, 1, N);
    ofstream a ("Data.txt");
    for (int p: pa)
    {
        S += p;
        a.write((char*)(&p), sizeof(int));
    }
    cout<<S<<endl;
    return 0;
}
