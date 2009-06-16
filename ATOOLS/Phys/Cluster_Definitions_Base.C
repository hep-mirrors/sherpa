#include "ATOOLS/Phys/Cluster_Definitions_Base.H"

using namespace ATOOLS;

std::ostream &ATOOLS::operator<<(std::ostream &str,const CParam &cp)
{
  return str<<"CP{kt="<<sqrt(cp.m_kt2)<<",op="<<
    (cp.m_op2<0.0?"-":"")<<sqrt(dabs(cp.m_op2))<<",x="<<cp.m_x<<"}";
}

Cluster_Definitions_Base::~Cluster_Definitions_Base() 
{
}
