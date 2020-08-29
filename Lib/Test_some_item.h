//
// Created by alexander on 29/08/2020.
//

#ifndef TEST_SOME_ITEM
#define TEST_SOME_ITEM

#include <vector>
#include <ostream>

class [[maybe_unused]] Test_some_item
{
    int i{};/// because default constructor doesn't init this value
    std::vector<unsigned int> vv;
public:
    [[nodiscard]] int get_i() const
    {
        return i;
    }

    friend std::ostream& operator<<(std::ostream& os, const Test_some_item& item)
    {
        os << "i: " << item.i << " vv: ";// << item.vv;
        return os;
    }

    bool operator==(const Test_some_item& rhs) const;

    bool operator!=(const Test_some_item& rhs) const;

    void set_i(int i_par)
    {
        Test_some_item::i = i_par;
    }

    [[nodiscard]] const std::vector<unsigned int>& get_vv() const
    {
        return vv;
    }

    void set_vv(const std::vector<unsigned int>& vv_par)
    {
        Test_some_item::vv = vv_par;
    }

    Test_some_item()
    = default;

public:
    Test_some_item(int i, std::vector<unsigned int> vv);

    explicit Test_some_item(std::vector<unsigned int> vv);

    virtual ~Test_some_item();

    bool operator<(const Test_some_item& rhs) const;

    bool operator>(const Test_some_item& rhs) const;

    bool operator<=(const Test_some_item& rhs) const;

    bool operator>=(const Test_some_item& rhs) const;
};


#endif //TEST_SOME_ITEM
