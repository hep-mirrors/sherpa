#include "ATOOLS/Phys/Cluster_Leg.H"

#include "ATOOLS/Org/STL_Tools.H"
#include "ATOOLS/Org/MyStrStream.H"
#include <iomanip>

using namespace ATOOLS;

ClusterLeg_PVector::ClusterLeg_PVector()
{
#ifdef USING__Threading
  pthread_mutex_init(&m_mtx,NULL);
#endif
}

ClusterLeg_PVector::~ClusterLeg_PVector()
{
#ifdef USING__Threading
  pthread_mutex_destroy(&m_mtx);
#endif
  while (!empty()) {
    delete back();
    pop_back();
  }
}

ClusterLeg_PVector Cluster_Leg::s_legs;

Cluster_Leg *Cluster_Leg::New
(Cluster_Amplitude *const ampl,const Cluster_Leg &ref)
{
  s_legs.MtxLock();
  if (s_legs.empty()) {
    s_legs.MtxUnLock();
    return new Cluster_Leg(ampl,ref);
  }
  Cluster_Leg *cl(s_legs.back());
  s_legs.pop_back();
  s_legs.MtxUnLock();
  *cl=ref;
  cl->p_ampl=ampl;
  return cl;
}

Cluster_Leg *Cluster_Leg::New
(Cluster_Amplitude *const ampl,const Vec4D &p,
 const Flavour &fl,const ColorID &c)
{
  s_legs.MtxLock();
  if (s_legs.empty()) {
    s_legs.MtxUnLock();
    return new Cluster_Leg(ampl,p,fl,c);
  }
  Cluster_Leg *cl(s_legs.back());
  s_legs.pop_back();
  s_legs.MtxUnLock();
  cl->p_ampl=ampl;
  cl->m_id=0;
  cl->m_st=0;
  cl->m_n=0;
  cl->m_d=0;
  cl->m_k=0;
  cl->m_p=p;
  cl->m_fl=fl;
  cl->m_c=c;
  return cl;
}

void Cluster_Leg::Delete()
{
  s_legs.MtxLock();
  s_legs.push_back(this);
  s_legs.MtxUnLock();
}

namespace ATOOLS {

  std::ostream &operator<<(std::ostream &ostr,const ColorID &col)
  {
    return ostr<<'('<<col.m_i<<','<<col.m_j<<')';
  }

  std::ostream &operator<<(std::ostream &ostr,const Cluster_Leg &leg)
  {
    ostr<<std::right<<std::setw(12)<<ToString(ID(leg.Id()))
	<<std::setw(12)<<leg.Flav()
	<<" "<<std::left<<leg.Mom()<<" "<<leg.Col();
    ostr<<" ["<<leg.Stat()<<"|"<<leg.NMax()<<","<<leg.DMax()<<"]";
    if (leg.K()>0) ostr<<" <-> "<<ID(leg.K());
    return ostr;
  }

}
