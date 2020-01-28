#include "CFPSHOWER++/Calculators/IF/SF_IF12.H"
#include "CFPSHOWER++/Shower/Kernel.H"
#include "PHASIC++/Channels/CSS_Kinematics.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"

using namespace CFPSHOWER;
using namespace PHASIC;
using namespace ATOOLS;

SF_IF12::SF_IF12(const Kernel_Info & info) : SF_Base(info) {}

double SF_IF12::Jacobean(const Splitting & split) const {
  double eta_old = split.eta();
  double eta_new = split.eta()/split.z();
  size_t beam    = split.GetSplitter()->Beam()-1;
  const Flavour & flav = m_flavs[m_tags[0]];
  double pdf_old = split.GetKernel()->GetXPDF(eta_old,split.t(),m_split,beam);
  double pdf_new = split.GetKernel()->GetXPDF(eta_new,split.t(),   flav,beam);
  if (dabs(pdf_old) < (split.GetKernel()->PDFMinValue() *
		       log(1.-eta_old)/log(1.-split.GetKernel()->PDFXMin())) ) {
    return 0.;
  }
  return pdf_new/pdf_old;
}

double SF_IF12::PDFEstimate(const Splitting & split) const
{
  // estimate of pdf ratios, based on current values of splitting kinematics
  double eta     = split.eta();
  double scale2  = Min(split.tstart(),split.Q2());
  size_t beam    = split.GetSplitter()->Beam()-1;
  const Flavour & flav = m_flavs[m_tags[0]];
  double pdf_old = split.GetKernel()->GetXPDF(eta,scale2,m_split,beam);
  double pdf_new = split.GetKernel()->GetXPDF(eta,scale2,flav,beam);
  // some extra treatment for heavy to light transitions in initial state
  // TODO: have to add threshold effect of heavy quarks as forced transition
  if (m_flavs[0].Mass(true)<1.0 && m_split.Mass(true)>=1.0) {
    double tcut(Max(split.tcut(),sqr(2.0*flav.Mass(true))));
    double pdf_old1 = split.GetKernel()->GetXPDF(eta,tcut,m_split,beam);
    double pdf_new1 = split.GetKernel()->GetXPDF(0.2,tcut,   flav,beam);
    if (pdf_old1!=0. && dabs(pdf_old1)<dabs(pdf_old)) pdf_old = pdf_old1;
    if (                dabs(pdf_new1)>dabs(pdf_new)) pdf_new = pdf_new1;
  }
  // estimated minimal value of pdf at splitting point.
  double pdf_min = (split.GetKernel()->PDFMinValue() *
		    log(1.0-split.eta()) /
		    log(1.0-split.GetKernel()->PDFXMin()));
  if (dabs(pdf_old)<pdf_min) return 0.;
  if (dabs(pdf_new)<pdf_min) pdf_new = pdf_old;
  msg_Out()<<METHOD<<"("<<m_split<<" <- "<<flav<<" + "<<m_flavs[m_tags[1]]<<") = "
	   <<dabs(pdf_new/pdf_old)<<".\n";
  return dabs(pdf_new/pdf_old);
}

bool SF_IF12::Construct(Splitting & split,const int & mode) {
  double y = split.t() / (split.Q2red() * (1.-split.z()));
  Kin_Args kin_args(y,split.z(),split.phi(),split.KinScheme());
  Parton * helper(NULL);
  /* 
     TODO: Have to implement different kin scheme.
     if (split.KinScheme()==0)
     for (size_t i(0);i<split.GetSplitter->Ampl()->size();++i)
     if ((*s.p_c->Ampl())[i]->Beam()==3-s.p_c->Beam()) {
     b=(*s.p_c->Ampl())[i];
     break;
     }
  */
  if (ConstructIFDipole(split.m2(0),split.m2(1),
			split.msplit2(),split.mspect2(),
			helper?split.GetKernel()->GetMSel()->Mass2(helper->Flav()):0.0,
			-split.GetSplitter()->Mom(),split.GetSpectator()->Mom(),
			helper?-helper->Mom():Vec4D(),kin_args) < 0) {
    //msg_Out()<<"    didn't construct dipole.\n";
    return false;
  }
  // TODO: Maybe have to add remnant test?
  // remember: pi is incoming particle, with flavour fl[1], resulting from splitting
  //           pj is emitted particle, with flavour fl[2]
  //           pk is spectator
  split.Set_y(y);
  m_moms[0] = -kin_args.m_pi;
  m_moms[1] =  kin_args.m_pj;
  m_specmom =  kin_args.m_pk;
  split.SetSequence(kin_args.m_lam);
  return 1;
}
