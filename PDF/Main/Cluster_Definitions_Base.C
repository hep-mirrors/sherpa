#include "PDF/Main/Cluster_Definitions_Base.H"

#include "ATOOLS/Org/Exception.H"

using namespace PDF;
using namespace ATOOLS;

std::ostream &PDF::operator<<(std::ostream &str,const CParam &cp)
{
  return str<<"CP{kt="<<sqrt(cp.m_kt2)<<",op="<<
    (cp.m_op2<0.0?"-":"")<<sqrt(dabs(cp.m_op2))
	    <<",x="<<cp.m_x<<",mu="<<sqrt(cp.m_mu2)
	    <<",k="<<cp.m_kin<<",m="<<cp.m_mode<<"}";
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

int Cluster_Definitions_Base::ReCluster
(Cluster_Amplitude *const ampl)
{
  DEBUG_FUNC(this);
  msg_Debugging()<<*ampl<<"\n";
  for (Cluster_Amplitude *campl(ampl->Next());
       campl;campl=campl->Next()) {
    Cluster_Leg *lij(NULL);
    for (size_t ij(0);ij<campl->Legs().size();++ij)
      if (campl->Leg(ij)->K()) {
	lij=campl->Leg(ij);
	break;
      }
    if (lij==NULL) THROW(fatal_error,"Invalid amplitude");
    int i(-1), j(-1), k(-1);
    for (size_t l(0);l<campl->Prev()->Legs().size();++l) {
      Cluster_Leg *cl(campl->Prev()->Leg(l));
      if (cl->Id()&lij->Id()) {
	if (i>=0) j=l;
	else i=l;
      }
      if (cl->Id()==lij->K()) k=l;
    }
    Vec4D_Vector p=Combine
      (*campl->Prev(),i,j,k,lij->Flav(),campl->MS(),
       campl->Kin(),(lij->Stat()&2)?1:0);
    if (p.empty()) return -1;
    for (size_t m(0), n(0);n<campl->Legs().size();++m) {
      if (m==j) continue;
      if (m>n+1) THROW(fatal_error,"Invalid PS history");
      for (size_t l(0);l<campl->Legs().size();++l)
	if (campl->Leg(l)->Id()&
	    campl->Prev()->Leg(m)->Id()) {
	  campl->Leg(l)->SetMom(p[n++]);
	  break;
	}
    }
    msg_Debugging()<<*campl<<"\n";
  }
  return 1;
}
