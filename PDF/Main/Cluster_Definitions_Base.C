#include "PDF/Main/Cluster_Definitions_Base.H"

#include "ATOOLS/Org/Exception.H"

using namespace PDF;
using namespace ATOOLS;

bool Cluster_Config::operator<(const Cluster_Config &cc) const
{
  if (m_i<cc.m_i) return true;
  if (m_i>cc.m_i) return false;
  if (m_j<cc.m_j) return true;
  if (m_j>cc.m_j) return false;
  if (m_k<cc.m_k) return true;
  if (m_k>cc.m_k) return false;
  return m_mo<cc.m_mo;
}

std::ostream &PDF::operator<<(std::ostream &str,const Cluster_Config &cc)
{
  return str<<"CC{ampl="<<cc.p_ampl<<",ms="<<cc.p_ms
	    <<",i="<<cc.m_i<<",j="<<cc.m_j<<",k="<<cc.m_k<<",mo="<<cc.m_mo
	    <<",kin="<<cc.m_kin<<",mode="<<ID(cc.m_mode)<<"}";
}

std::ostream &PDF::operator<<(std::ostream &str,const Cluster_Param &cp)
{
  return str<<"CP{op="<<cp.m_op
	    <<",kt="<<(cp.m_kt2<0.0?"-":"")<<sqrt(dabs(cp.m_kt2))
	    <<",mu="<<(cp.m_mu2<0.0?"-":"")<<sqrt(dabs(cp.m_mu2))
	    <<",cpl="<<cp.m_cpl<<",kin="<<cp.m_kin
	    <<",mode="<<cp.m_mode<<",stat="<<cp.m_stat<<"}";
}

bool Cluster_Config::PureQCD() const
{
  return m_mo.Strong() &&
    p_ampl->Leg(m_i)->Flav().Strong() &&
    p_ampl->Leg(m_j)->Flav().Strong() &&
    p_ampl->Leg(m_k)->Flav().Strong();
}

Cluster_Definitions_Base::Cluster_Definitions_Base()
{
}

Cluster_Definitions_Base::~Cluster_Definitions_Base() 
{
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
	if (cl->Id()==campl->Prev()->IdNew()) j=l;
	else i=l;
      }
      if (cl->Id()==lij->K()) k=l;
    }
    if ((lij->Stat()&4) &&
	campl->Legs().size()<campl->NIn()+2) {
      for (size_t l(0);l<campl->Legs().size();++l) {
	Cluster_Leg *cl(campl->Leg(l));
	if (cl!=lij) cl->SetMom(campl->Prev()->IdLeg(cl->Id())->Mom());
      }
      lij->SetMom(campl->Prev()->Leg(i)->Mom()+
		  campl->Prev()->Leg(j)->Mom());
    }
    else {
      Cluster_Param cp=campl->CA<Cluster_Definitions_Base>()->Cluster
	(Cluster_Config(campl->Prev(),i,j,k,lij->Flav(),campl->MS(),
			NULL,campl->Kin(),((lij->Stat()&4)?1:0)|
			(campl->NLO()?16:0)));
      if (cp.m_pijt==Vec4D()) {
	cp=campl->CA<Cluster_Definitions_Base>()->Cluster
	  (Cluster_Config(campl->Prev(),j,i,k,lij->Flav(),campl->MS(),
			  NULL,campl->Kin(),((lij->Stat()&4)?1:0)|
			  (campl->NLO()?16:0)));
	if (cp.m_pijt==Vec4D()) return -1;
      }
      for (size_t n(0);n<campl->Legs().size();++n) {
	Cluster_Leg *c(campl->Leg(n));
	if (c->Id()&campl->Prev()->Leg(i)->Id()) c->SetMom(cp.m_pijt);
	else if (c->Id()&campl->Prev()->Leg(k)->Id()) c->SetMom(cp.m_pkt);
	else c->SetMom(cp.m_lam*campl->Prev()->IdLeg(c->Id())->Mom());
      }
    }
    msg_Debugging()<<*campl<<"\n";
  }
  return 1;
}
