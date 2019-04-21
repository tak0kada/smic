#include <iostream>
#include <vector>
#include <utility>

template<typename T1, typename T2>
std::ostream& operator<<(std::ostream& os, std::pair<T1, T2> p)
{
    os << "pair(" << p.first << ", " << p.second << ")";
    return os;
}

template<typename T>
std::ostream& operator<<(std::ostream& os, std::vector<T> vec)
{
    os << "[";
    for (const auto& v: vec)
        os << v << ", ";
    os << "]";

    return os;
}
