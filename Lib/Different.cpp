#include <cassert>
using namespace std;

constexpr double Pow_impl (const double x, const int i, const double Product)
{
    assert(i >= 0);
    if (i == 1)
        return Product * x;
    if (i == 0)
        return Product;
    return Pow_impl(x, i - 1, Product * x);
}
double Pow (const double x, const int i) noexcept
{
    return Pow_impl(x, i, 1);
}
