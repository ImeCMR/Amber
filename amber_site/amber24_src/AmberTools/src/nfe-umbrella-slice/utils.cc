#ifdef HAVE_CONFIG_H
#  include "nfe.h"
#endif // HAVE_CONFIG_H

#include <cstdio>
#include <cstdlib>
#include <cstdarg>

#include <sstream>

#include "utils.h"

namespace nfe { namespace utils {

void fatal(const char* fmt, ...)
{
    va_list ap;

    va_start(ap, fmt);
    fputs(" ** Error ** : ", stderr);
    vfprintf(stderr, fmt, ap);
    fflush(stderr);
    va_end(ap);

    exit(EXIT_FAILURE);
}

//
//  exceptionally naive :-) ... but enough for our purposes
//

std::list<std::string> tokenize
    (const std::string& source, const std::string& sep, const std::string& ign)
{
    std::list<std::string> result;

    size_t i, i0 = 0;

    for (i = 0; i < source.length(); ++i) {
        if (sep.find(source[i]) != std::string::npos) {

            if (i0 != i) {
                std::string tmp;

                for (size_t j = i0; j < i; ++j)
                    if (ign.find(source[j]) == std::string::npos)
                        tmp += source[j];

                if (!tmp.empty())
                    result.push_back(tmp);
            }
            i0 = i + 1;
        }
    }

    //
    // the last one
    //

    if (i0 != i) {
        std::string tmp;

        for (size_t j = i0; j < i; ++j)
            if (ign.find(source[j]) == std::string::npos)
                tmp += source[j];

        if (!tmp.empty())
            result.push_back(tmp);
    }

    return result;
}

int to_int(const std::string& s) throw()
{
    int res;

    std::istringstream iss(s);
    iss >> res;

    if (iss.fail() || !iss.eof())
        fatal("'%s' is not an integer number\n", s.c_str());

    return res;
}

double to_double(const std::string& s) throw()
{
    double res;

    std::istringstream iss(s);
    iss >> res;

    if (iss.fail() || !iss.eof())
        fatal("'%s' is not a real number\n", s.c_str());

    return res;
}

}} // namespace nfe::utils
