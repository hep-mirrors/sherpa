#include "CFPSHOWER++/Calculators/IF/SF_IF.H"
#include "CFPSHOWER++/Shower/Kernel.H"
#include "PHASIC++/Channels/CSS_Kinematics.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"

using namespace CFPSHOWER;
using namespace PHASIC;
using namespace ATOOLS;

SF_IF::SF_IF(const Kernel_Info & info) : SF_Base(info) {}

double SF_IF::Jacobean(const Splitting & split) const {
  double eta_old = split.Eta();
  double eta_new = split.Eta()/split.X();
  //if (eta_new>1.) return 0.;
  // ratio of pdfs s part of the Jacobean
  size_t beam    = split.GetSplitter()->Beam()-1;
  double pdf_old = split.GetKernel()->GetXPDF(eta_old,split.T(),m_flavs[0],beam);
  double pdf_new = split.GetKernel()->GetXPDF(eta_new,split.T(),m_flavs[1],beam);
  if (dabs(pdf_old) < (split.GetKernel()->PDFMinValue() *
		       log(1.-eta_old)/log(1.-split.GetKernel()->PDFXMin())) ) {
    //msg_Out()<<"*** "<<METHOD<<" yields 0 from estimate: "
    //	     <<pdf_old<<" < "<<split.GetKernel()->PDFMinValue()<<" * log(1-"<<eta_old<<")/"
    //	     <<"log(1-"<<split.GetKernel()->PDFXMin()<<")\n"
    //	     <<"    t = "<<split.T()<<", x = "<<split.X()<<", z = "<<split.Z()<<": pdfratio = "
    //	     <<pdf_new/pdf_old<<" for "<<m_flavs[0]<<" <-- "<<m_flavs[1]<<"\n";
    return 0.;
  }
  // msg_Out()<<"*** "<<METHOD<<" passed first test: "
  // 	   <<pdf_old<<" > "<<split.GetKernel()->PDFMinValue()<<" * log(1-"<<eta_old<<")/"
  // 	   <<"log(1-"<<split.GetKernel()->PDFXMin()<<")\n"
  // 	   <<"    t = "<<split.T()<<", x = "<<split.X()<<", z = "<<split.Z()<<": pdfratio =  "
  // 	   <<pdf_new/pdf_old<<" for "<<m_flavs[0]<<" <-- "<<m_flavs[1]<<"\n";
  return pdf_new/pdf_old;
}

double SF_IF::PDFEstimate(const Splitting & split) const
{
  // estimate of pdf ratios, based on current values of splitting kinematics
  double eta     = split.Eta();
  double scale2  = Min(split.T1(),split.Q2());
  size_t beam    = split.GetSpectator()->Beam()-1;
  double pdf_old = split.GetKernel()->GetXPDF(eta,scale2,m_flavs[0],beam);
  double pdf_new = split.GetKernel()->GetXPDF(eta,scale2,m_flavs[1],beam);
  // some extra treatment for heavy to lgiht transitions in initial state
  // TODO: have to add threshold effect of heavy quarks as forced transition
  if (m_flavs[1].Mass(true)<1.0 && m_flavs[0].Mass(true)>=1.0) {
    double tcut(Max(split.T0(),sqr(2.0*m_flavs[0].Mass(true))));
    double pdf_old1 = split.GetKernel()->GetXPDF(eta,tcut,m_flavs[0],beam);
    double pdf_new1 = split.GetKernel()->GetXPDF(0.2,tcut,m_flavs[1],beam);
    if (pdf_old1!=0. && dabs(pdf_old1)<dabs(pdf_old)) pdf_old = pdf_old1;
    if (                dabs(pdf_new1)>dabs(pdf_new)) pdf_new = pdf_new1;
  }
  // estimated minimal value of pdf at splitting point.
  double pdf_min = (split.GetKernel()->PDFMinValue() *
		    log(1.0-split.Eta()) /
		    log(1.0-split.GetKernel()->PDFXMin()));
  if (dabs(pdf_old)<pdf_min) return 0.;
  if (dabs(pdf_new)<pdf_min) return 1.;
  return dabs(pdf_new/pdf_old);
}


bool SF_IF::InitKinematics(Splitting & split) const {
  split.SetY(split.T() / (split.Q2() * (1.-split.Z())));
  split.SetX(split.Z());
  return true;
}

int SF_IF::Construct(Splitting & split) const {
  Kin_Args kin_if(split.Y(),split.X(),split.phi(),split.KinScheme());
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
  //msg_Out()<<"*** "<<METHOD<<"(kin = "<<split.KinScheme()<<"): "
  //	   <<"x, y, z = "<<split.X()<<", "<<split.Y()<<", "<<split.Z()<<"\n";
  if (ConstructIFDipole(split.mi2(),split.mj2(),split.mij2(),split.mk2(),
			helper?split.GetKernel()->GetMSel()->Mass2(helper->Flav()):0.0,
			-split.GetSplitter()->Mom(),split.GetSpectator()->Mom(),
			helper?-helper->Mom():Vec4D(),kin_if) < 0) {
    //msg_Out()<<"    didn't construct dipole.\n";
    return -1;
  }
  // TODO: Maybe have to add remnant test?
  // remember: pi is incoming particle, with flavour fl[1], resulting from splitting
  //           pj is emitted particle, with flavour fl[2]
  //           pk is spectator
  split.SetMom(0, -kin_if.m_pi);
  split.SetMom(1,  kin_if.m_pj);
  split.SetSpecMom(kin_if.m_pk);
  split.SetSequence(kin_if.m_lam);
  //msg_Out()<<"    worked out ok.\n";
  return 1;
}
