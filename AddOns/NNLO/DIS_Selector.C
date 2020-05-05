#include "PHASIC++/Selectors/Selector.H"
#include "PHASIC++/Process/Process_Base.H"

namespace SHNNLO {

  class DIS_Selector : public PHASIC::Selector_Base {
    PHASIC::Process_Base *p_proc;
    double m_qtmin, m_q2;
    int m_type;
    double R2(const ATOOLS::Vec4D &p1,const ATOOLS::Vec4D &p2) const;
  public:
    DIS_Selector(const PHASIC::Selector_Key &key);
    double KT2(ATOOLS::Vec4D_Vector &moms,const int beam);
    bool Trigger(ATOOLS::Selector_List &sl);
    void BuildCuts(PHASIC::Cut_Data *) {}
  };

}// end of namespace SHNNLO

#include "PHASIC++/Main/Process_Integrator.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"

#define s_ymax std::numeric_limits<double>::max()

using namespace SHNNLO;
using namespace PHASIC;
using namespace ATOOLS;

DIS_Selector::DIS_Selector(const Selector_Key &key):
  Selector_Base("DISNNLO_Selector"), p_proc(key.p_proc)
{
  m_nin=key.p_proc->NIn();
  m_nout=key.p_proc->NOut();
  m_n=m_nin+m_nout;
  m_smin=0.0;
  m_smax=sqr(rpa->gen.Ecms());
  m_sel_log = new Selector_Log(m_name);
  Scoped_Settings s{ key.m_settings };
  const auto parameters = s.SetDefault<std::string>({}).GetVector<std::string>();
  assert(parameters[0] == "DISNNLO");
  m_qtmin=s.Interprete<double>(parameters[1]);
  m_type=m_nout-(p_proc->Info().Has(nlo_type::real)?2+1:2);
}

bool DIS_Selector::Trigger(ATOOLS::Selector_List &sl) 
{
  DEBUG_FUNC(p_proc->Name());
  Vec4D_Vector p(1,sl[1].Momentum());
  for (size_t i(3);i<m_nin+m_nout;++i)
    p.push_back(sl[i].Momentum());
  double kt2(KT2(p,1));
  msg_Debugging()<<"kT = "<<sqrt(kt2)<<"\n";
  if (kt2==-1.0) return 1-m_sel_log->Hit(1-false);
  int trig=(m_type==0 && kt2<sqr(m_qtmin)) ||
    (m_type==1 && kt2>sqr(m_qtmin));
  return 1-m_sel_log->Hit(1-trig);
}

double DIS_Selector::KT2(Vec4D_Vector &moms,const int beam)
{
  DEBUG_FUNC("beam = "<<beam);
  Vec4D qq, pp(rpa->gen.PBeam(beam));
  pp[0]=pp.PSpat();
  for (size_t i(0);i<moms.size();++i) qq+=i<1?-moms[i]:moms[i];
  Poincare cms(pp+qq);
  double Q2(-qq.Abs2()), x(Min(Q2/(2.0*pp*qq),1.0));
  msg_Debugging()<<"Q^2 = "<<Q2<<"\n";
  m_q2=Q2;
  double E(sqrt(Q2)/((beam?2.0:-2.0)*x));
  Vec4D p(sqrt(E*E+pp.Abs2()),0.0,0.0,-E);
  Vec4D q(0.0,0.0,0.0,2.0*x*E);
  cms.Boost(pp);
  cms.Boost(qq);
  Poincare zrot(pp,beam?-Vec4D::ZVEC:Vec4D::ZVEC);
  zrot.Rotate(pp);
  zrot.Rotate(qq);
  Poincare breit(p+q);
  breit.BoostBack(pp);
  breit.BoostBack(qq);
  if (Q2*(1.0-x)>x) // only if hfs energy large enough
    if (!IsEqual(pp,p,1.0e-3) || !IsEqual(qq,q,1.0e-3))
      msg_Error()<<METHOD<<"(): Boost error."<<std::endl;
  Vec4D sum;
  Vec4D_Vector jets;
  for (int i(0);i<moms.size();++i) {
    msg_Debugging()<<"p["<<i<<"] = "<<moms[i];
    cms.Boost(moms[i]);
    zrot.Rotate(moms[i]);
    breit.BoostBack(moms[i]);
    msg_Debugging()<<" -> "<<moms[i]<<"\n";
    sum+=i==0?-moms[i]:moms[i];
  }
  msg_Debugging()<<"mom sum = "<<sum<<"\n";
  int ii=0, jj=0, n=moms.size();
  std::vector<int> p_imap(n);
  std::vector<std::vector<double> > p_ktij(n,std::vector<double>(n));
  double dmin=std::numeric_limits<double>::max();
  for (int i=0;i<n;++i) {
    p_imap[i]=i;
    double di=p_ktij[i][i]=sqr(moms[i][0])*R2(moms[i],pp);
    if (di<dmin) { dmin=di; ii=i; jj=i; }
    for (int j=0;j<i;++j) {
      double dij=p_ktij[i][j]=sqr(Min(moms[i][0],moms[j][0]))*R2(moms[i],moms[j]);
      if (dij<dmin) { dmin=dij; ii=i; jj=j; }
    }
  }
  while (n>2) {
    msg_Debugging()<<"Q_{"<<n<<"->"<<n-1<<"} = "<<sqrt(dmin)
		   <<" <- "<<p_imap[jj]<<" & "<<p_imap[ii]<<"\n";
    if (ii!=jj) moms[p_imap[jj]]+=moms[p_imap[ii]];
    --n;
    for (int i=ii;i<n;++i) p_imap[i]=p_imap[i+1];
    int jjx=p_imap[jj];
    p_ktij[jjx][jjx]=sqr(moms[jjx][0])*R2(moms[jjx],pp);
    for (int j=0;j<jj;++j) p_ktij[jjx][p_imap[j]] = 
      sqr(Min(moms[jjx][0],moms[p_imap[j]][0]))
      *R2(moms[jjx],moms[p_imap[j]]);
    for (int i=jj+1;i<n;++i) p_ktij[p_imap[i]][jjx] = 
      sqr(Min(moms[jjx][0],moms[p_imap[i]][0]))
      *R2(moms[jjx],moms[p_imap[i]]);
    ii=jj=0;
    dmin=p_ktij[p_imap[0]][p_imap[0]];
    for (int i=0;i<n;++i) {
      int ix=p_imap[i]; double di=p_ktij[ix][ix];
      if (di<dmin) { dmin=di; ii=jj=i; }
      for (int j=0;j<i;++j) {
	int jx=p_imap[j]; double dij=p_ktij[ix][jx];
	if (dij<dmin) { dmin=dij; ii=i; jj=j; }
      }
    }
  }
  return dmin;
}

double DIS_Selector::R2(const Vec4D &p1, const Vec4D &p2) const
{
  return 2.0*(1.0-p1.CosTheta(p2));
}

DECLARE_GETTER(DIS_Selector,"DISNNLO",Selector_Base,Selector_Key);

Selector_Base *ATOOLS::Getter<Selector_Base,Selector_Key,DIS_Selector>::
operator()(const Selector_Key &key) const
{
  return new DIS_Selector(key);
}

void ATOOLS::Getter<Selector_Base,Selector_Key,DIS_Selector>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"NNLO selector";
}
