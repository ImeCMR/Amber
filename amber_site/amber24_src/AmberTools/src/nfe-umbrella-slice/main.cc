#ifdef HAVE_CONFIG_H
#   include "nfe.h"
#endif /* HAVE_CONFIG_H */

#include <string>
#include <iomanip>
#include <iostream>
#include <cassert>
#include <stdlib.h>

#ifdef HAVE_GETOPT_H
#  ifndef _GNU_SOURCE
#    define _GNU_SOURCE 1
#  endif // _GNU_SOURCE
#  include <unistd.h>
#  include <getopt.h>
#else
#  include <unistd.h>
#endif // HAVE_GETOPT_H

extern "C" {
    extern char *optarg;
    extern int optind, opterr, optopt;
} // extern "C"

#include "utils.h"
#include "umbrella.h"
#include "dimension.h"

using namespace std;
using namespace nfe;

namespace {

void show_usage()
{
    cout <<
    "\n"
    "Usage: nfe-umbrella-slice [options] filename\n"
    "\n"
    "Options:\n"
    "  -h --help\n"
    "  -p --pretend\n"
    "  -g --gradient (off by default)\n"
    "  -r --reset (off by default)\n"
    "  -t --translate (0.0 by default)\n"
    "  -d --dimensions\n"
    "\n"
    "Example of dimensions' argument: \"1.0:0.2,2.0,23:3.0\" \n\n"
    "(the slice along 2nd variable will be taken from 0.2 to 2.0 using\n"
    "23 points; values of the first and the last variables will be set\n"
    "to 1.0 and 3.0 correspondingly)\n"
    "\n";
}

namespace options {

static const char* umbrella_filename = NULL;
static const char* dimensions_string = NULL;

static bool pretend     = false;
static bool do_gradient = false;
static bool reset       = false;

static double translation = 0.0;

void parse(int argc, char** argv)
{
    int c;

    static const char short_options[] = "hd:gprt:";

#   ifdef HAVE_GETOPT_H
    static const struct option long_options[] = {
        {"help",       0, 0, 'h'},
        {"pretend",    0, 0, 'p'},
        {"gradient",   0, 0, 'g'},
        {"reset",      0, 0, 'r'},
        {"translate",  1, 0, 't'},
        {"dimensions", 1, 0, 'd'},
        {0,            0, 0, 0}
    };
#   endif // HAVE_GETOPT_H

    while (true) {
#       ifdef HAVE_GETOPT_H
        c = getopt_long(argc, argv, short_options, long_options, 0);
#       else
        c = getopt(argc, argv, short_options);
#       endif // HAVE_GETOPT_H

        if (c == -1)
	         break;

        switch (c) {
            case 'h':
                show_usage();
                exit(EXIT_SUCCESS);
            case 'd':
                assert(optarg);
                options::dimensions_string = (*optarg ? optarg : 0);
                break;
            case 'g':
                options::do_gradient = true;
                break;
            case 'p':
                options::pretend = true;
                break;
            case 'r':
                options::reset = true;
                break;
            case 't':
                options::translation = utils::to_double(optarg);
                break;
            case '?':
                exit(EXIT_FAILURE);
                break;
            default:
                assert(false); // should not be reached
        }
    } // while (true)

    if (optind == argc)
        utils::fatal("need more arguments (try '-h' for hints)\n");

    options::umbrella_filename = argv[optind++];
	
    if (optind < argc) {
        cerr << " ** Warning ** : ignoring unexpected command-line arguments:";
        for (; optind < argc; ++optind)
            cerr << " '" << argv[optind] << "'";
        cerr << endl;
    }
}

} // namespace options

} // namespace

int main(int argc, char** argv)
{
    options::parse(argc, argv);
    assert(options::umbrella_filename);

    // load the umbrella
    umbrella* u = umbrella::load(options::umbrella_filename);
    assert(u);
    assert(u->nextents() > 0);

    // setup dimensions
    vector<dimension> dims(u->nextents());

    if (options::dimensions_string) {

        list<string> toks =
            utils::tokenize(options::dimensions_string, ":", " ");

        if (toks.size() != dims.size())
            utils::fatal("'%s' : unexpected number of dimensions\n",
                options::dimensions_string);

        for (int n(0); n < u->nextents(); ++n) {
            dims[n].set_from_string(toks.front());
            toks.pop_front();
        }

    } else {

        for (int n(0); n < u->nextents(); ++n) {
            dims[n].npoints = u->extent(n);
            dims[n].origin = u->origin(n);
            dims[n].spacing = u->spacing(n);
        }
    }

    cout << "#\n# source = '" << options::umbrella_filename << "'\n#\n";

    cout << "# source::extents     = " << u->extent(0);
    for (int n(1); n < u->nextents(); ++n)
        cout << " x " << u->extent(n);
    cout << '\n';

    cout << "# source::origin      = " << u->origin(0);
    for (int n(1); n < u->nextents(); ++n)
        cout << ", " << u->origin(n);
    cout << '\n';

    cout << "# source::spacing     = " << u->spacing(0);
    for (int n(1); n < u->nextents(); ++n)
        cout << ", " << u->spacing(n);
    cout << '\n';

    cout << "# source::periodicity = " << (u->periodicity(0) ? "yes" : "no");
    for (int n(1); n < u->nextents(); ++n)
        cout << (u->periodicity(n) ? ", yes" : ", no");
    cout << "\n#\n";

    if (options::dimensions_string) {
        cout << "# dimensions::npoints = " << dims[0].npoints;
        for (int n(1); n < u->nextents(); ++n)
            cout << " x " << dims[n].npoints;
        cout << '\n';

        cout << "# dimensions::extents = ["
             << dims[0].min() << ", " << dims[0].max() << ']';
        for (int n(1); n < u->nextents(); ++n)
            cout << " x [" << dims[n].min() << ", " << dims[n].max() << ']';
        cout << "\n#\n";
    }

    cout << "# translation = " << options::translation
         << "\n# gradients = " << (options::do_gradient ? "on" : "off")
         << "\n#\n";

    if (options::pretend) {
        delete u;
        return EXIT_SUCCESS;
    }

    int curr[umbrella::MAX_NEXTENTS];
    for (int n(0); n < u->nextents(); ++n)
        curr[n] = 0;

    int last_running = u->nextents() - 1;
    for (; last_running >= 0; --last_running)
        if (dims[last_running].npoints != 1)
            break;

    cout << scientific << setprecision(6) << right;

    if (options::reset) {
      double max_u = 0.0;
      while (true) {

        double pos[umbrella::MAX_NEXTENTS];
        for (int n(0); n < u->nextents(); ++n) {
            pos[n] = dims[n](curr[n]);
        }
        
        if (u->eval(pos) > max_u)
           max_u = u->eval(pos);

        int k = u->nextents() - 1;
        for (; k >= 0; --k) {
            if (curr[k] != dims[k].npoints - 1)
               break;
            else
               curr[k] = 0;
        }

        if (k < 0)
            break;

        curr[k]++;
      }
      options::translation = max_u;
    }

    while (true) {

        double pos[umbrella::MAX_NEXTENTS];
        for (int n(0); n < u->nextents(); ++n) {
            pos[n] = dims[n](curr[n]);
            cout << setw(15) << pos[n];
        }

        if (options::do_gradient) {

            double grad[umbrella::MAX_NEXTENTS];

            cout << setw(15) << (options::translation - u->eval(pos, grad));
            for (int n(0); n < u->nextents(); ++n) {
                cout << setw(15) << -grad[n];
            }
            cout << '\n';

        } else {
            cout << setw(15) << (options::translation - u->eval(pos)) << '\n';
        }

        int k = u->nextents() - 1;
        for (; k >= 0; --k) {
            if (curr[k] != dims[k].npoints - 1)
               break;
            else
               curr[k] = 0;
        }

        if (k < 0)
            break;

        if (k + 1 == last_running)
            cout << '\n';

        curr[k]++;
    }

    return EXIT_SUCCESS;
}
