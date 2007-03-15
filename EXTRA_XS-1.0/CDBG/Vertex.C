#include "Vertex.H"

#include "Message.H"
#include "STL_Tools.H"
#ifdef PROFILE__all
#include "prof.hh"
#else
#define PROFILE_HERE
#define PROFILE_LOCAL(NAME)
#endif

using namespace EXTRAXS;
using namespace ATOOLS;

Vertex::Vertex(Current *const c): 
  p_out(c), p_a(NULL), p_b(NULL)
{
  p_out->AttachIn(this);
}

Vertex::~Vertex() 
{
}

std::ostream &EXTRAXS::operator<<(std::ostream &str,const Vertex &v)
{
  if (v.JA()!=NULL) 
    str<<'{'<<v.JA()->Type()<<','<<v.JA()->Flav()<<'}'<<v.JA()->Id();
  if (v.JB()!=NULL) 
    str<<"(+){"<<v.JB()->Type()<<','<<v.JB()->Flav()<<'}'<<v.JB()->Id();
  return str;
}
