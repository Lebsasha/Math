//
// Created by alexander on 29/08/2020.
//

#include "Test_some_item.h"

#include <utility>
#include <tuple>

Test_some_item::Test_some_item(int i, std::vector<unsigned int> vv) : i(i), vv(std::move(vv))
{}

bool Test_some_item::operator==(const Test_some_item& rhs) const
{
    return i == rhs.i &&
           vv == rhs.vv;
}

Test_some_item::Test_some_item(std::vector<unsigned int> vv) : vv(std::move(vv))
{}

bool Test_some_item::operator!=(const Test_some_item& rhs) const
{
    return !(rhs == *this);
}

bool Test_some_item::operator<(const Test_some_item& rhs) const
{
    return std::tie(i, vv) < std::tie(rhs.i, rhs.vv);
}

bool Test_some_item::operator>(const Test_some_item& rhs) const
{
    return rhs < *this;
}

bool Test_some_item::operator<=(const Test_some_item& rhs) const
{
    return !(rhs < *this);
}

bool Test_some_item::operator>=(const Test_some_item& rhs) const
{
    return !(*this < rhs);
}

Test_some_item::~Test_some_item()
= default;
