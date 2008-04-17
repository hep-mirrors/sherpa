#include "Splitting_Kernel.H"

using namespace PDF;
using namespace ATOOLS;

Splitting_Kernel::Splitting_Kernel
(const Flavour &a,const Flavour &b,const Flavour &c):
  m_fla(a), m_flb(b), m_flc(c) {}

Splitting_Kernel::~Splitting_Kernel() {}

double Splitting_Kernel::Value(const double &z) const
{
  return (*this)(z);
}
