#ifdef HAVE_CONFIG_H
#   include "nfe.h"
#endif // HAVE_CONFIG_H

#include "utils.h"
#include "dimension.h"

namespace nfe {

void dimension::set_from_string(const std::string& s)
{
    std::list<std::string> toks = utils::tokenize(s, ",", " ");

    if (toks.size() == 1) {

        npoints = 1;
        origin = utils::to_double(toks.back());
        spacing = 0.0;

    } else if (toks.size() == 3) {

        npoints = utils::to_int(toks.back());
        toks.pop_back();

        if (npoints < 2)
            utils::fatal("too small number of points in '%s'\n", s.c_str());

        const double high = utils::to_double(toks.back());
        toks.pop_back();

        origin = utils::to_double(toks.back());
        toks.pop_back();

        if (origin > high)
            utils::fatal("invalid range in '%s' : %f > %f\n",
                         s.c_str(), origin, high);

        spacing = (high - origin)/(npoints - 1);

    } else {
        utils::fatal("invalid dimension specification '%s'\n", s.c_str());
    }
}

} // namespace nfe
