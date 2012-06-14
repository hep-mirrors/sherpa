#include "PDF/Main/Cluster_Definitions_Base.H"

#include "ATOOLS/Org/Exception.H"

using namespace PDF;
using namespace ATOOLS;

std::ostream &PDF::operator<<(std::ostream &str,const CParam &cp)
{
  return str<<"CP{kt="<<sqrt(cp.m_kt2)<<",op="<<
    (cp.m_op2<0.0?"-":"")<<sqrt(dabs(cp.m_op2))
	    <<",x="<<cp.m_x<<",mu="<<sqrt(cp.m_mu2)
	    <<",k="<<cp.m_kin<<"}";
}

Cluster_Definitions_Base::Cluster_Definitions_Base() : m_amode(0)
{
}

Cluster_Definitions_Base::~Cluster_Definitions_Base() 
{
}

CParam Cluster_Definitions_Base::CoreScale
(ATOOLS::Cluster_Amplitude *const ampl)
{
  THROW(fatal_error,"Invalid function call");
}
