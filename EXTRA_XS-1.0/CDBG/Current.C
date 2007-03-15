#include "Current.H"

#include "Vertex.H"
#include "Message.H"
#include "Exception.H"
#include "MyStrStream.H"
#ifdef PROFILE__all
#include "prof.hh"
#else
#define PROFILE_HERE
#define PROFILE_LOCAL(NAME)
#endif

using namespace EXTRAXS;
using namespace ATOOLS;

Current::Current(CDBG_Amplitude *const ampl,const Flavour &fl):
  p_ampl(ampl), 
  m_fl(fl) 
{
  if (m_fl.IsFermion()) m_type=ct::spinor;
  else if (m_fl.IsVector()) m_type=ct::vector;
  else if (m_fl.IsTensor()) m_type=ct::tensor;
  else THROW(fatal_error,"Invalid flavour");
}

Current::~Current()
{
  for (Vertex_Vector::const_iterator vit(m_in.begin());
       vit!=m_in.end();++vit) {
    delete *vit;
  }
}

std::ostream &EXTRAXS::operator<<(std::ostream &str,const ct::type &type)
{
  switch (type) {
  case ct::spinor: return str<<"S";
  case ct::vector: return str<<"V";
  case ct::tensor: return str<<"T";
  }
  return str;
}

std::ostream &EXTRAXS::operator<<(std::ostream &str,const Current &c)
{
  return str<<'('<<c.Type()<<','<<c.Flav()<<','
	    <<c.Id().size()<<','<<c.Key()<<')';
}

