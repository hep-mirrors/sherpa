// profiling aid
// written by Thomas Kunert

#ifndef __PROF_H
#define __PROF_H

#include <iostream>
#include <iosfwd>
#include <vector>
#include <map>



struct call_node;

struct ProfStatic;

class ProfStaticP
{
public:
  ProfStaticP( const char* s );
  ~ProfStaticP();
  ProfStatic* get() const { return p; }
private:
  ProfStatic* p;
};
  

class Prof {
public: 
    Prof( ProfStatic* ps );
    ~Prof();
private:
    ProfStatic* ps;
    int start;
    unsigned long start_mem;
    call_node* last_call_node_;
};

void print_profile( std::ostream& ost );
void set_prof();

#define PROFILE_HERE \
{}
/*
  static ProfStaticP _prof_static_(__PRETTY_FUNCTION__); \
  Prof _prof_( _prof_static_.get() );
*/

#define PROFILE_LOCAL(LOCALNAME) \
{}
/*
  static ProfStaticP _prof_static_(LOCALNAME); \
  Prof _prof_( _prof_static_.get() )
*/

#endif
