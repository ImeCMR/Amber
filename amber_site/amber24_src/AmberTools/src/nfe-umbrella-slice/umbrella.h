#ifndef NFE_UMBRELLA_H
#define NFE_UMBRELLA_H

#include <cassert>

#include "utils.h"

namespace nfe {

struct umbrella : private utils::noncopyable {
    ~umbrella() throw();

    static umbrella* load(const char*);

    inline int nextents() const throw();
    inline int extent(int) const throw();

    inline bool periodicity(int) const throw();

    inline double origin(int) const throw();
    inline double spacing(int) const throw();

    double eval(const double*) const throw();
    double eval(const double*, double*) const throw();

    static const int MAX_NEXTENTS = 6;

    double coeff(const int*) const throw();

private:

    int m_nextents;
    int m_extents[MAX_NEXTENTS];

    bool m_periodicity[MAX_NEXTENTS];

    double m_origin[MAX_NEXTENTS];
    double m_spacing[MAX_NEXTENTS];

    double* m_coeffs;

private:
    umbrella() {}
};

inline int umbrella::nextents() const throw()
{
    return m_nextents;
}

inline int umbrella::extent(int n) const throw()
{
    assert(n >= 0 && n < MAX_NEXTENTS);
    return m_extents[n];
}

inline bool umbrella::periodicity(int n) const throw()
{
    assert(n >= 0 && n < MAX_NEXTENTS);
    return m_periodicity[n];
}

inline double umbrella::origin(int n) const throw()
{
    assert(n >= 0 && n < MAX_NEXTENTS);
    return m_origin[n];
}

inline double umbrella::spacing(int n) const throw()
{
    assert(n >= 0 && n < MAX_NEXTENTS);
    return m_spacing[n];
}

} // namespace nfe

#endif // NFE_UMBRELLA_H
