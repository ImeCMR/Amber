#ifndef NFE_UTILS_H
#define NFE_UTILS_H

#include <list>
#include <string>

#include "macros.h"

namespace nfe { namespace utils {

void fatal(const char*, ...) NFE_GNUC_PRINTF(1,2) NFE_GNUC_NORETURN;

struct noncopyable {
protected:
  noncopyable() {}
  ~noncopyable() {}

private:
  noncopyable(const noncopyable&);
  const noncopyable& operator=(const noncopyable&);
};

int to_int(const std::string&) throw();
double to_double(const std::string&) throw();

//
//  splits given string into tokens; i.e. for
//  s = "sda : asda as :  asda ; asda", sep=":;", ignore = " "
//  it returns 4 tokens "sda", "asdaas", "asda", "asda"
//

std::list<std::string> tokenize(const std::string& s,
    const std::string& sep, const std::string& ignore = "");

}} // namespace nfe::utils

#endif // NFE_UTILS_H
