#ifndef _Utils_hpp_
#define _Utils_hpp_

#include <string>

bool DirExists( std::string dname );
bool FileExists( std::string fname );
bool FileIsLink( std::string fname );
void CopyFile( std::string src, std::string dest );
std::string GetCwd();
bool RelativeSymLink( std::string oldfile, std::string newlink );

double SwitchOff( double const r , double const rlo , double const rhi );
double SwitchOn( double const r , double const rlo , double const rhi );

#endif
