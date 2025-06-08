#ifndef NFE_DIMENSIONS_H
#define NFE_DIMENSIONS_H

#include <vector>
#include <string>

namespace nfe {

struct dimension {

    int npoints;
    double origin;
    double spacing;

    double min() const throw()
    {
        return origin;
    }

    double max() const throw()
    {
        return this->operator()(npoints - 1);
    }

    double operator()(int n) const throw()
    {
        return origin + n*spacing;
    }

    void set_from_string(const std::string&);
};

} // namespace nfe

#endif // NFE_DIMENSIONS_H
