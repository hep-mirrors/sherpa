#include "ATOOLS/Phys/NLO_Subevt.H"

#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "ATOOLS/Org/MyStrStream.H"

using namespace ATOOLS;

ME_wgtinfo::ME_wgtinfo():
  m_nx(0), m_w0(0.0), p_wx(NULL),
  m_y1(1.0), m_y2(1.0), m_mur2(0.0) {}

ME_wgtinfo::~ME_wgtinfo()
{
  if (m_nx>0) delete [] p_wx;
}

ME_wgtinfo &ME_wgtinfo::operator*=(const double &scal)
{
  m_w0*=scal;
  for (int i=0;i<m_nx;i++) p_wx[i]*=scal;
  return *this;
}

void ME_wgtinfo::Flip()
{
  std::swap<double>(m_x1,m_x2);
  std::swap<double>(m_y1,m_y2);
  if (m_nx>=10) for (int i=0;i<4;i++) std::swap<double>(p_wx[i+2],p_wx[i+6]);
  if (m_nx>=18) for (int i=0;i<4;i++) std::swap<double>(p_wx[i+10],p_wx[i+14]);
}

void ME_wgtinfo::AddMEweights(int n)
{
  m_nx=n;
  p_wx=new double[m_nx];
  for (int i=0;i<m_nx;i++) p_wx[i]=0.;
}

namespace ATOOLS {
  template <> Blob_Data<ME_wgtinfo*>::~Blob_Data() {}
  template class Blob_Data<ME_wgtinfo*>;
}

bool IDip_ID::operator<(const IDip_ID &di) const
{
  if (m_ijt<di.m_ijt) return true;
  if (m_ijt>di.m_ijt) return false;
  return m_kt<di.m_kt;
}

bool DDip_ID::operator<(const DDip_ID &di) const
{
  if (m_i<di.m_i) return true;
  if (m_i>di.m_i) return false;
  if (m_j<di.m_j) return true;
  if (m_j>di.m_j) return false;
  return m_k<di.m_k;
}

bool Dip_ID::operator<(const Dip_ID &di) const
{
  if (m_ijt<di.m_ijt) return true;
  if (m_ijt>di.m_ijt) return false;
  if (m_kt<di.m_kt) return true;
  if (m_kt>di.m_kt) return false;
  return DDip_ID::operator<(di);
}

NLO_subevt::~NLO_subevt()
{
  if (m_delete) {
    delete[] p_fl;
    delete[] p_mom;
    delete[] p_id;
  }
  if (p_ampl) p_ampl->Delete();
}

void NLO_subevt::CopyXSData(const NLO_subevt *sub)
{
  m_me=sub->m_me;
  m_mewgt=sub->m_mewgt;
  for (size_t i(0);i<m_mu2.size();++i) m_mu2[i]=sub->m_mu2[i];
  m_result=0.0;
  if (p_ampl) {
    p_ampl->Delete();
    p_ampl=NULL;
  }
  if (sub->p_ampl) {
    p_ampl=sub->p_ampl->CopyAll();
    for (Cluster_Amplitude *ampl(p_ampl);
	 ampl;ampl=ampl->Next()) ampl->SetProc(p_proc);
  }
}

void NLO_subevtlist::Mult(const double &scal)
{
  for (const_iterator it=begin();it!=end();it++) (*it)->Mult(scal);
}

void NLO_subevtlist::MultME(const double &scal)
{
  for (const_iterator it=begin();it!=end();it++) (*it)->MultME(scal);
}

void NLO_subevtlist::MultMEwgt(const double &scal)
{
  for (const_iterator it=begin();it!=end();it++) (*it)->MultMEwgt(scal);
}

ATOOLS::Particle_List *NLO_subevt::CreateParticleList() const
{
  ATOOLS::Particle_List * pl = new ATOOLS::Particle_List;
  for (size_t i=2;i<m_n;i++) {
    pl->push_back(new ATOOLS::Particle(i,p_fl[i],p_mom[i]));
  }
  return pl;
}

std::string NLO_subevt::IDString(const int mode) const
{
  std::string tag;
  bool si(mode&&(p_id[0]&1)==0);
  for (size_t i(0);i<m_n;++i) tag+=ToString(p_id[i<2&&si?1-i:i])+"_";
  tag+="_"+ToString(1<<m_i)+"_"+ToString(1<<m_j)+"_"+ToString(1<<m_k);
  return tag;
}

NLO_subevtlist &NLO_subevtlist::operator*=(const double scal)
{
  for (const_iterator it=begin();it!=end();it++) {
    (*it)->m_result*=scal;
  }
  return *this;
}

namespace ATOOLS
{
  std::ostream &operator<<(std::ostream &ostr,const IDip_ID &idi)
  {
    return ostr<<"["<<idi.m_ijt<<"]<->["<<idi.m_kt<<"]";
  }
  std::ostream &operator<<(std::ostream &ostr,const DDip_ID &ddi)
  {
    return ostr<<"("<<ddi.m_i<<","<<ddi.m_j<<")<->("<<ddi.m_k<<")";
  }
  std::ostream &operator<<(std::ostream &ostr,const Dip_ID &di)
  {
    return ostr<<"["<<di.m_ijt<<"]("<<di.m_i<<","<<di.m_j
	       <<")<->["<<di.m_kt<<"]("<<di.m_k<<")";
  }
  std::ostream &operator<<(std::ostream &ostr,const NLO_subevt &sevt)
  {
    std::vector<int> ids;
    ATOOLS::Flavour_Vector flavs;
    for (size_t i(0);i<sevt.m_n;++i) {
      flavs.push_back(sevt.p_fl[i]);
      ids.push_back(sevt.p_id[i]);
    }
    return ostr<<sevt.m_pname<<" "<<(Dip_ID)(sevt)
	       <<", idx "<<sevt.m_idx
               <<" {\n  fl: "<<flavs<<", id: "<<ids
               <<"\n  result = "<<sevt.m_result
               <<",  ME = "<<sevt.m_me<<" ("<<sevt.m_trig
               <<")\n  \\mu_Q = "<<sqrt(sevt.m_mu2[stp::res])
	       <<",  \\mu_F = "<<sqrt(sevt.m_mu2[stp::fac])
	       <<", \\mu_R = "<<sqrt(sevt.m_mu2[stp::ren])
	       <<", k_T = "<<sqrt(sevt.m_kt2)<<"\n}";
  }
  std::ostream &operator<<(std::ostream &ostr,const stp::id &scl)
  {
    switch (scl) {
    case stp::ren: return ostr<<"ren";
    case stp::fac: return ostr<<"fac";
    case stp::res: return ostr<<"res";
    case stp::size: return ostr<<"<error>";
    }
    return ostr<<"<unknown>";
  }
}

namespace ATOOLS {
  template <> Blob_Data<NLO_subevtlist*>::~Blob_Data() {}
  template class Blob_Data<NLO_subevtlist*>;
}

