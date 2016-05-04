#include "COMIX/Phasespace/PS_Info.H"

using namespace COMIX;

void PS_Info::Add(const CObject *c)
{
}

void PS_Info::Divide(const double &d)
{
}

void PS_Info::Multiply(const Complex &c)
{
}

void PS_Info::Conjugate()
{
}

void PS_Info::Invert()
{
}

bool PS_Info::IsZero() const
{
  return false;
}

ATOOLS::AutoDelete_Vector<PS_Info> PS_Info::s_objects;

PS_Info *PS_Info::New()
{
#ifndef USING__Threading
  if (s_objects.empty())
#endif
    return new PS_Info();
  PS_Info *v(s_objects.back());
  s_objects.pop_back();
  return v;
}

PS_Info *PS_Info::New(const PS_Info &s)
{
#ifndef USING__Threading
  if (s_objects.empty())
#endif
    return new PS_Info(s);
  PS_Info *v(s_objects.back());
  s_objects.pop_back();
  *v=s;
  return v;
}

METOOLS::CObject *PS_Info::Copy() const
{
  return PS_Info::New(*this);
}

void PS_Info::Delete()
{
#ifndef USING__Threading
  s_objects.push_back(this);
#else
  delete this;
#endif
}

std::ostream &COMIX::operator<<(std::ostream &str,const PS_Info &s)
{
  return str<<'{'<<s.H()<<';'<<s(0)<<","<<s(1)<<'|'<<s[0]<<'}';
}

