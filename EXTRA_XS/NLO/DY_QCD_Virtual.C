#include "PHASIC++/Process/Process_Info.H"
#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "EXTRA_XS/NLO/Logarithms.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"

using namespace PHASIC;
using namespace ATOOLS;

namespace EXTRAXS {

  class DY_QCD_Virtual : public PHASIC::Virtual_ME2_Base {
  private:
    int m_l1, m_l2;
  public:
    DY_QCD_Virtual(const Process_Info &pi,const Flavour_Vector &flavs):
      Virtual_ME2_Base(pi,flavs), m_l1(-1), m_l2(-1)
    {
      for (size_t i(0);i<flavs.size();++i)
	if (flavs[i].IsLepton() || flavs[i].IsNeutrino()) {
	  if (m_l1<0) m_l1=i;
	  else if (m_l2<0) m_l2=i;
	  else THROW(fatal_error,"Invalid flavor configuration");
	}
    }
    double Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom)
    {
      return 4.*M_PI;
    }
    void Calc(const ATOOLS::Vec4D_Vector&p)
    {
      Vec4D p1(m_l1<2?-p[m_l1]:p[m_l1]);
      Vec4D p2(m_l2<2?-p[m_l2]:p[m_l2]);
      Complex l12=LnRat(m_mur2,-(p1+p2).Abs2());
      m_res.IR2()=-2.;
      m_res.IR()=(-3.-2.*l12).real();
      m_res.Finite()=(-8.-3.*l12-sqr(l12)).real();
      m_res*=4./3.;
    }
  };// end of class DY_QCD_Virtual

  class Singlet_QCD_Virtual : public PHASIC::Virtual_ME2_Base {
  public:
    Singlet_QCD_Virtual(const Process_Info &pi,const Flavour_Vector &flavs):
      Virtual_ME2_Base(pi,flavs)
    {
    }
    double Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom)
    {
      return 4.*M_PI;
    }
    void Calc(const ATOOLS::Vec4D_Vector&p)
    {
      m_res=METOOLS::DivArrD(0.);
      for (size_t i(0);i<p_plist->size();++i) {
	for (size_t j(i+1);j<p_plist->size();++j) {
	  size_t i1((*p_plist)[i]), i2((*p_plist)[j]);
	  Vec4D p1(i1<2?-p[i1]:p[i1]);
	  Vec4D p2(i2<2?-p[i2]:p[i2]);
	  Complex l12=LnRat(m_mur2,-(p1+p2).Abs2());
	  m_res.IR2()-=-2.*(*p_dsij)[i][j];
	  m_res.IR()-=(-3.-2.*l12).real()*(*p_dsij)[i][j];
	  m_res.Finite()-=(-8.-3.*l12-sqr(l12)).real()*(*p_dsij)[i][j];
	}
      }
      m_res*=1./(*p_dsij)[0][0];
    }
  };// end of class Singlet_QCD_Virtual

}

using namespace EXTRAXS;

DECLARE_VIRTUALME2_GETTER(EXTRAXS::DY_QCD_Virtual,"DY_QCD_Virtual")
Virtual_ME2_Base *ATOOLS::Getter
<PHASIC::Virtual_ME2_Base,PHASIC::Process_Info,EXTRAXS::DY_QCD_Virtual>::
operator()(const Process_Info &pi) const
{
  DEBUG_FUNC(pi);
  if (pi.m_loopgenerator!="Internal") return NULL;
  if (pi.m_fi.m_nlotype==nlo_type::loop) {
    if (pi.m_mincpl[0]!=1. || pi.m_maxcpl[0]!=1.) return NULL;
    if (pi.m_fi.m_nlocpl[0]!=1. || pi.m_fi.m_nlocpl[1]!=0.) return NULL;
    if (pi.m_mincpl[1]!=2. || pi.m_maxcpl[1]!=2.) {
      Flavour_Vector fl=pi.ExtractFlavours();
      for (size_t i(0);i<fl.size();++i)
	if (fl[i].Strong() && !fl[i].IsQuark()) return NULL;
      return new Singlet_QCD_Virtual(pi,fl);
    }
    if (pi.m_fi.m_ps.size()!=2) return NULL;
    if (pi.m_fi.m_asscontribs!=asscontrib::none) {
      msg_Error()<<"DY_QCD_Virtual(): Error: cannot provide requested "
                 <<"associated contributions "<<pi.m_fi.m_asscontribs<<std::endl;
      return NULL;
    }
    Flavour_Vector fl=pi.ExtractFlavours();
    for (size_t i(0);i<fl.size();++i)
      if (fl[i].IsMassive()) return NULL;
    if (((fl[0].IsLepton() && fl[2].IsNeutrino()) ||
	 (fl[2].IsLepton() && fl[0].IsNeutrino())) &&
	fl[2].LeptonFamily()==fl[0].LeptonFamily()) {
      if (fl[1].IsQuark() && fl[3].QuarkFamily()==fl[1].QuarkFamily())
	return new DY_QCD_Virtual(pi,fl);
    }
    if ((fl[2].IsLepton() && fl[3].IsNeutrino() &&
	 fl[3].LeptonFamily()==fl[2].LeptonFamily())) {
      if (fl[0].IsQuark() && fl[1].QuarkFamily()==fl[0].QuarkFamily())
	return new DY_QCD_Virtual(pi,fl);
    }
    if ((fl[0].IsLepton() && fl[1]==fl[0].Bar()) ||
	(fl[0].IsNeutrino() && fl[1]==fl[0].Bar())) {
      if (fl[2].IsQuark() && fl[3]==fl[2].Bar())
	return new DY_QCD_Virtual(pi,fl);
    }
    if ((fl[0].IsLepton() && fl[2]==fl[0]) ||
	(fl[0].IsNeutrino() && fl[2]==fl[0])) {
      if (fl[1].IsQuark() && fl[3]==fl[1])
	return new DY_QCD_Virtual(pi,fl);
    }
    if ((fl[2].IsLepton() && fl[3]==fl[2].Bar()) ||
	(fl[2].IsNeutrino() && fl[3]==fl[2].Bar())) {
      if (fl[0].IsQuark() && fl[1]==fl[0].Bar())
	return new DY_QCD_Virtual(pi,fl);
    }
  }
  return NULL;
}
