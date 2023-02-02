#include "AMISIC++/Perturbative/QED_Processes.H"
#include "ATOOLS/Math/Random.H"

using namespace AMISIC;
using namespace ATOOLS;
using namespace std;

//////////////////////////////////////////////////////////////////////////
//
// only on-shell photons - quark charges handled in collection of xsecs
//
//////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////
// qg -> qgam
//////////////////////////////////////////////////////////////////////////

qg_qgamma::qg_qgamma() : XS_Base() { m_name = string("qg->qgam"); }

void qg_qgamma::Calc(const double & s,const double & t,const double & u) {
  m_lastxs = -1./3. * (s/u + u/s);
}

bool qg_qgamma::SetColours(const ATOOLS::Flavour_Vector & flavs) {
    /*
      
      q----+~~~~3
           |
           |  t
           |
    1-q====+----2
      
    */
  size_t quark = size_t(!flavs[0].IsQuark()), gluon = 1-quark;
  size_t anti = size_t(flavs[0].IsQuark()?flavs[0].IsAnti():flavs[1].IsAnti());
  m_colours[quark][anti]   = m_colours[gluon][1-anti]   = Flow::Counter();
  m_colours[2+quark][anti] = m_colours[gluon][anti]     = Flow::Counter();
  m_colours[quark][1-anti] = m_colours[2+quark][1-anti] = 0; 
  m_colours[3-quark][0]    = m_colours[3-quark][1]      = 0; 
  return true; 
}

DECLARE_XSBASE_GETTER(qg_qgamma,"XS_qg_qgamma")
XS_Base * Getter<XS_Base,Flavour_Vector,qg_qgamma>::
operator()(const ATOOLS::Flavour_Vector& flavs) const
{
  if (flavs.size()!=4) return NULL;
  if ((flavs[0].IsQuark() && flavs[1].IsGluon() && flavs[2]==flavs[0] &&
       flavs[3].IsPhoton()) ||
      (flavs[1].IsQuark() && flavs[0].IsGluon() && flavs[3]==flavs[1] &&
       flavs[2].IsPhoton())) {
    if ((flavs[0].IsQuark() && flavs[0].Mass()<1.e-6) ||
	(flavs[1].IsQuark() && flavs[1].Mass()<1.e-6)) return new qg_qgamma();
    else {
      msg_Error()<<"Error in "<<METHOD<<":\n"
		 <<"   no massive matrix element yet for Qg->Qgamma "
		 <<"with Q = "<<(flavs[0].IsQuark()?flavs[0]:flavs[1])<<".\n";
      exit(1);
    }

  }
  return NULL;
}


//////////////////////////////////////////////////////////////////////////
// qqbar_ggamma
//////////////////////////////////////////////////////////////////////////

qqbar_ggamma::qqbar_ggamma() : XS_Base() { m_name = string("qqbar->ggamma"); }

void qqbar_ggamma::Calc(const double & s,const double & t,const double & u) {
  m_lastxs = 8./9. * (t/u + u/t);
}

bool qqbar_ggamma::SetColours(const ATOOLS::Flavour_Vector & flavs) {
  size_t anti = size_t(flavs[0].IsAnti());
  m_colours[0][anti]   = m_colours[2][anti]   = Flow::Counter();
  m_colours[1][1-anti] = m_colours[2][1-anti] = Flow::Counter();
  m_colours[0][1-anti] = m_colours[1][anti] = 0;
  return true;
}

DECLARE_XSBASE_GETTER(qqbar_ggamma,"XS_qqbar_ggamma")
XS_Base * Getter<XS_Base,Flavour_Vector,qqbar_ggamma>::
operator()(const ATOOLS::Flavour_Vector& flavs) const
{
  if (flavs.size()!=4) return NULL;
  if (flavs[0].IsQuark() && flavs[1]==flavs[0].Bar() &&
      flavs[2].IsGluon() && flavs[3].IsPhoton()) {
    if (flavs[0].Mass()<=1.e-6) return new qqbar_ggamma();
    else {
      msg_Error()<<"Error in "<<METHOD<<":\n"
		 <<"   no massive matrix element yet for QQbar -> ggamma"
		 <<"with Q = "<<flavs[0]<<".\n";
      exit(1);
    }
  }
  return NULL;
}
