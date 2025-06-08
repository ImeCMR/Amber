#ifdef HAVE_CONFIG_H
#   include "nfe.h"
#endif // HAVE_CONFIG_H

#include <new>
#include <limits>

#include <cmath>
#include <cassert>

#include <netcdf.h>

#include "utils.h"
#include "umbrella.h"

namespace {

void check_nc_status(const char* filename, const int& status)
{
    using namespace nfe;

    if (status != NC_NOERR)
        utils::fatal("'%s' : %s\n", filename, nc_strerror(status));
}

void check_att_len(const char* filename, const char* attname,
                   const int& nextents, const size_t& att_len)
{
    using namespace nfe;

    if (att_len != static_cast<size_t>(nextents))
        utils::fatal("'%s' : length of the '%s' attribute conflicts"
                     " with the 'nextents' value (%d vs %d)\n",
                     filename, attname, static_cast<int>(att_len), nextents);
}

//
// evaluates the value (f) and derivative (df)
// of the cubic B-spline centered at 0.0
// (works correctly only for -2 < x < 2)
//

void m4_vdv(const double& x, double& f, double& df) throw()
{
    if (x < -1.0) {
        const double x2 = x + 2.0;
        df = 0.5*x2*x2;
         f = (1.0/3.0)*x2*df;
    } else if (x < 0.0) {
        const double x2 = -0.5*x;
         f = x2*x*(x + 2.0) + (2.0/3.0);
        df = x2*(3.0*x + 4.0);
    } else if (x < 1.0) {
        const double x2 = 0.5*x;
         f = x2*x*(x - 2.0) + (2.0/3.0);
        df = x2*(3.0*x - 4.0);
    } else {
        const double x2 = 2.0 - x;
        df = (-0.5)*x2*x2;
         f = (-1.0/3.0)*x2*df;
    }
}

//
// evaluates the value of the cubic B-spline
// centered at 0.0 (works correctly only for -2 < x < 2)
//

double m4_v(const double& x)
{
    if (x < -1.0) {
        return (1.0/6.0)*std::pow(x + 2.0, 3);
    } else if (x < 0.0) {
        return (-0.5)*x*x*(x + 2.0) + (2.0/3.0);
    } else if (x < 1.0) {
        return (-0.5)*x*x*(2.0 - x) + (2.0/3.0);
    } else {
        return (1.0/6.0)*std::pow(2.0 - x, 3);
    }
}

//
// evaluates \int_{0}^{a}dx m4_v(x) [assumes 0 <= a <= 2]
//

double m4_i(const double& a)
{
    assert(a >= 0.0);
    assert(a <= 2.0);

    if (a < 1.0) {
        const double a2 = a*a;
        return (((1.0/8.0)*a - (1.0/3.0))*a2 + (2.0/3.0))*a;
    } else {
        const double tmp1 = 2.0 - a;
        const double tmp2 = tmp1*tmp1;
        return (1.0/2.0) - (1.0/24.0)*tmp2*tmp2;
    }
}

} // namespace

namespace nfe {

umbrella* umbrella::load(const char* filename)
{
    using utils::fatal;

    assert(filename != 0 && *filename);

    int ncid, dimid, ndims;

    size_t a_len;
    nc_type a_type;

#   define INVOKE_NC(what) check_nc_status(filename, (what));

    INVOKE_NC(nc_open(filename, NC_NOWRITE, &ncid))

    int nextents;

    INVOKE_NC(nc_inq_attlen(ncid, NC_GLOBAL, "nextents", &a_len))
    if (a_len != 1)
        fatal("'%s' : attribute 'nextents' is not scalar\n", filename);

    INVOKE_NC(nc_inq_atttype(ncid, NC_GLOBAL, "nextents", &a_type))
    if (a_type != NC_INT)
        fatal("'%s' : attribute 'nextents' is not integer\n", filename);

    INVOKE_NC(nc_get_att_int(ncid, NC_GLOBAL, "nextents", &nextents))
    if (nextents <= 0)
        fatal("'%s' : attribute 'nextents' is not positive (%d)\n",
            filename, nextents);
    if (nextents > MAX_NEXTENTS)
        fatal("'%s' : attribute 'nextents' is too big (%d > %d)\n",
            filename, nextents, MAX_NEXTENTS);

    // extents

    int extents[MAX_NEXTENTS];

    INVOKE_NC(nc_inq_attlen(ncid, NC_GLOBAL, "extents", &a_len))
    check_att_len(filename, "extents", nextents, a_len);

    INVOKE_NC(nc_inq_atttype(ncid, NC_GLOBAL, "extents", &a_type))
    if (a_type != NC_INT)
        fatal("'%s' : type of the 'extents' attribute is not integer\n",
            filename);

    INVOKE_NC(nc_get_att_int(ncid, NC_GLOBAL, "extents", extents))
    for (int i(0); i != nextents; ++i) {
        if (extents[i] <= 0)
            fatal("'%s' : value of the 'extents[%d]' is too small\n",
                filename, i);
    }

    // periodicity

    int periodicity[MAX_NEXTENTS];

    INVOKE_NC(nc_inq_attlen(ncid, NC_GLOBAL, "periodicity", &a_len))
    check_att_len(filename, "periodicity", nextents, a_len);

    INVOKE_NC(nc_inq_atttype(ncid, NC_GLOBAL, "periodicity", &a_type))
    if (a_type != NC_INT)
        fatal("'%s' : type of the 'periodicity' attribute is not integer\n",
            filename);

    INVOKE_NC(nc_get_att_int(ncid, NC_GLOBAL, "periodicity", periodicity))
    for (int i(0); i != nextents; ++i) {
        if (periodicity[i] != 0 && periodicity[i] != 1)
            fatal("'%s' : value of the 'periodicity[%d]' is not 0 or 1\n",
                filename, i);
    }

    // origin

    double origin[MAX_NEXTENTS];

    INVOKE_NC(nc_inq_attlen(ncid, NC_GLOBAL, "origin", &a_len))
    check_att_len(filename, "origin", nextents, a_len);

    INVOKE_NC(nc_inq_atttype(ncid, NC_GLOBAL, "origin", &a_type))
    if (a_type != NC_DOUBLE)
        fatal("'%s' : type of the 'origin' attribute is not double\n",
            filename);

    INVOKE_NC(nc_get_att_double(ncid, NC_GLOBAL, "origin", origin))

    // spacing

    double spacing[MAX_NEXTENTS];

    INVOKE_NC(nc_inq_attlen(ncid, NC_GLOBAL, "spacing", &a_len))
    check_att_len(filename, "spacing", nextents, a_len);

    INVOKE_NC(nc_inq_atttype(ncid, NC_GLOBAL, "spacing", &a_type))
    if (a_type != NC_DOUBLE)
        fatal("'%s' : type of the 'spacing' attribute is not double\n",
            filename);

    INVOKE_NC(nc_get_att_double(ncid, NC_GLOBAL, "spacing", spacing))
    for (int i(0); i != nextents; ++i) {
        if (spacing[i] <= 0.0)
            fatal("'%s' : value of the 'spacing[%d]' is not positive\n",
                filename, i);
    }

    // coeffs

    int varid;
    double* coeffs(0);

    INVOKE_NC(nc_inq_varid(ncid, "coeffs", &varid))

    INVOKE_NC(nc_inq_vartype(ncid, varid, &a_type))
    if (a_type != NC_DOUBLE)
        fatal("'%s' : type of the 'coeffs' variable is not double\n",
            filename);

    INVOKE_NC(nc_inq_varndims(ncid, varid, &ndims))
    if (ndims != 1)
        fatal("'%s' : 'coeffs' : number of dimensions != 1\n",
            filename);

    INVOKE_NC(nc_inq_vardimid(ncid, varid, &dimid))
    INVOKE_NC(nc_inq_dimlen(ncid, dimid, &a_len))

    size_t ncoeffs(1);
    for (int i(0); i != nextents; ++i)
        ncoeffs *= extents[i];

    if (ncoeffs != a_len)
        fatal("'%s' : length of the 'coeffs' does not"
              " correspond to the 'extents'\n", filename);

    assert(ncoeffs > 0);
    coeffs = new(std::nothrow) double[ncoeffs];
    if (!coeffs)
        fatal("'%s' : could not allocate memory for the 'coeffs'\n", filename);

    INVOKE_NC(nc_get_var_double(ncid, varid, coeffs))
    INVOKE_NC(nc_close(ncid))

    umbrella* ptr = new(std::nothrow) umbrella();
    if (!ptr)
        fatal("'%s' : could not allocate nfe::umbrella instance\n", filename);


    ptr->m_nextents = nextents;
    for (int i(0); i != nextents; ++i) {
        ptr->m_extents[i] = extents[i];
        ptr->m_periodicity[i] = (periodicity[i] == 1);
        ptr->m_origin[i] = origin[i];
        ptr->m_spacing[i] = spacing[i];
    }
    ptr->m_coeffs = coeffs;

#   undef INVOKE_NC

    return ptr;
}

umbrella::~umbrella() throw()
{
    assert(m_coeffs);
    delete[] m_coeffs;
}

double umbrella::eval(const double* x) const throw()
{
    assert(x);

    assert(m_nextents > 0);

    double m4v[MAX_NEXTENTS][4]; // values of the basis

    int gc[MAX_NEXTENTS];

    double f(0);

    for (int i(0); i != m_nextents; ++i) {
        assert(m_spacing[i] > 0.0);
        assert(m_extents[i] > 0);

        const double xs = (x[i] - m_origin[i])/m_spacing[i];
        gc[i] = static_cast<int>(std::floor(xs));

        for (int j(0); j != 4; ++j) {
            const double t = xs + static_cast<double>(1 - gc[i] - j);
            m4v[i][j] = m4_v(t);
        }
    }

    // loop over 4 x 4 x ... x 4 hypercube

    for (int n = 1 << (m_nextents << 1); n > 0; --n) {

        int p(n);
        size_t o(0);

        double m4v_prod(1);

        for (int i(0); i != m_nextents; ++i) {

            int j = gc[i] + (p & 3) - 1;

            if (j < 0 || j >= m_extents[i]) {
                if (m_periodicity[i]) {
                    j = (j < 0 ? m_extents[i] - 1 + (j + 1)%m_extents[i]
                               : j%m_extents[i]);
                } else {
                    goto next;
                }
            }

            o += (i + 1 == m_nextents ? j : m_extents[i + 1]*(o + j));

            m4v_prod *= m4v[i][p & 3];

            p >>= 2;
        }

        f += m_coeffs[o]*m4v_prod;

        next:;
    }

    return f;
}

double umbrella::eval(const double* x, double* df) const throw()
{
    assert(x);
    assert(df);

    assert(m_nextents > 0);

    double  m4v[MAX_NEXTENTS][4]; // values of the basis
    double m4dv[MAX_NEXTENTS][4]; // values of the basis' derivative

    int gc[MAX_NEXTENTS];

    double f(0);

    for (int i(0); i != m_nextents; ++i) {
        assert(m_spacing[i] > 0.0);
        assert(m_extents[i] > 0);

        df[i] = 0.0;

        const double xs = (x[i] - m_origin[i])/m_spacing[i];
        gc[i] = static_cast<int>(std::floor(xs));

        for (int j(0); j != 4; ++j) {
            const double t = xs + static_cast<double>(1 - gc[i] - j);
            m4_vdv(t, m4v[i][j], m4dv[i][j]);
        }
    }

    // loop over 4 x 4 x ... x 4 hypercube

    for (int n = 1 << (m_nextents << 1); n > 0; --n) {

        int p(n);
        size_t o(0);

        double m4v_prod(1);
        double m4dv_prod[MAX_NEXTENTS];

        for (int i(0); i != m_nextents; ++i)
            m4dv_prod[i] = 1.0;

        for (int i(0); i != m_nextents; ++i) {

            int g = gc[i] + (p & 3) - 1;

            if (g < 0 || g >= m_extents[i]) {
                if (m_periodicity[i]) {
                    g = (g < 0 ? m_extents[i] - 1 + (g + 1)%m_extents[i]
                               : g%m_extents[i]);
                } else {
                    goto next;
                }
            }

            m4v_prod *= m4v[i][p & 3];

            for (int j(0); j != i; ++j)
                m4dv_prod[j] *= m4v[i][p & 3];

            m4dv_prod[i] *= m4dv[i][p & 3];

            for (int j(i + 1); j != m_nextents; ++j)
                m4dv_prod[j] *= m4v[i][p & 3];

            o += (i + 1 == m_nextents ? g : m_extents[i + 1]*(o + g));
            p >>= 2;
        }

        f += m_coeffs[o]*m4v_prod;

        for (int i(0); i != m_nextents; ++i)
           df[i] += m_coeffs[o]*m4dv_prod[i];

        next:;
    }

    for (int i(0); i != m_nextents; ++i)
        df[i] /= m_spacing[i];

    return f;
}

double umbrella::coeff(const int* idx) const throw()
{
    assert(idx);
    assert(m_coeffs);
    assert(m_nextents > 0);

    size_t o(0);

    for (int i(0); i != m_nextents; ++i) {
        int g = idx[i];

        if (g < 0 || g >= m_extents[i]) {
            if (m_periodicity[i]) {
                g = (g < 0 ? m_extents[i] - 1 + (g + 1)%m_extents[i]
                           : g%m_extents[i]);
            } else {
                return 0.0;
            }
        }

        o += (i + 1 == m_nextents ? g : m_extents[i + 1]*(o + g));
    }

    return *(m_coeffs + o);
}

} // namespace nfe
