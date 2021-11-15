#include "CFPSHOWER++/Kinematics/Kinematics_Base.H"
#include "CFPSHOWER++/Shower/Kernel.H"
#include "PHASIC++/Channels/CSS_Kinematics.H"
#include "ATOOLS/Org/Message.H"

namespace CFPSHOWER {
  class FI_CS_Kinematics : public Kinematics_Base {
  private:
    bool   KinCheck(Splitting & split);
    double CalculateY(Splitting & split);
  public:
    FI_CS_Kinematics(const Kernel_Info & info);
    ~FI_CS_Kinematics() {}
    bool Init(Splitting & split, Configuration & config,
	      const ATOOLS::Mass_Selector * msel);
    bool operator()(Splitting & split, Configuration & config);
    bool InitLimits(Splitting & split, Configuration & config);
    void CalculateJacobean(Splitting & split, Configuration & config);
    bool UpdateSystem(Splitting & split, Configuration & config);
  };
}

using namespace CFPSHOWER;
using namespace PDF;
using namespace ATOOLS;

FI_CS_Kinematics::FI_CS_Kinematics(const Kernel_Info & info)  :
  Kinematics_Base(info)
{
  m_name   = std::string("FI: CS Kinematics");
  m_scheme = kin_type::CataniSeymour;
}

bool FI_CS_Kinematics::Init(Splitting & split, Configuration & config,
			    const ATOOLS::Mass_Selector * msel) {
  size_t beam    = split.GetSpectator()->Beam();
  PDF_Base * pdf = p_kernel->GetPDF(beam);
  if (beam<1 || beam>2 || pdf==NULL) {
    msg_Error()<<"Error in "<<METHOD<<"(beam = "<<beam<<"):\n"
	       <<"   We should not arrive here.  Will exit for debugging reasons.\n";
    exit(1);
    return false;
  }
  double xmin = Max(1.e-6,pdf->XMin());
  double xmax = Min(0.999999,pdf->XMax());
  double x    = split.GetSpectator()->XB();
  if (x>=xmax || x<=xmin || split.Q2()>=pdf->Q2Max() || split.Q2()<=pdf->Q2Min() || 
      split.Tcut()*x > split.Q2()*(1.-x)) return false;
  double disc = 1.-4.*Min(1.,x/(1-x))*split.Tcut()/split.Q2();
  if (disc<0.) return false;
  double delta = sqrt(disc), zmin = (1.-delta)/2., zmax = (1.+delta)/2.;
  if (zmin>zmax) return false;
  split.SetZmin(zmin);
  split.SetZmax(zmax);
  InitSystem(split,msel);
  return true;
}

bool FI_CS_Kinematics::KinCheck(Splitting & split) {
  split.SetY(m_y = CalculateY(split));
  return (split.Y()>=0.0 && split.Y()<=(1.0-split.GetSpectator()->XB()) &&
	  split.Q2()>split.Mass2(0)+split.Mass2(1)+m_mspect2);
}

double FI_CS_Kinematics::CalculateY(Splitting & split) {
  double arg = ( ( split.T() / (split.Z()*(1.-split.Z())) +
		   split.Mass2(0) * (1.-split.Z())/split.Z() +
		   split.Mass2(1) * split.Z()/(1.-split.Z()) ) /
		 ( split.Q2()-split.Mass2(0)-split.Mass2(1)-m_mspect2) );
  return 1./(1.-arg); 
}

bool FI_CS_Kinematics::operator()(Splitting & split, Configuration & config) {
  if (!KinCheck(split)) return false;
  // have to add mode=8 to kinargs to ensure the correct y is taken for the kinematics
  PHASIC::Kin_Args kinargs(split.Y(),split.Z(),split.Phi(),8); 
  if (PHASIC::ConstructFIDipole(split.Mass2(0),split.Mass2(1),
				m_msplit2,m_mspect2,
				m_psplit,m_pspect,kinargs) < 0) return false;
  split.SetMom(0,kinargs.m_pi);
  split.SetMom(1,kinargs.m_pj);
  split.SetMom(2,kinargs.m_pk);
  split.SetKinSpect(kinargs.m_pk);
  //msg_Out()<<METHOD
  //	   <<"(z = "<<split.Z()<<" vs "<<m_z[0]<<" from "<<split.Ztest()<<", "
  //	   <<"y = "<<m_y<<", phi = "<<split.Phi()<<"):\n"
  //	   <<m_psplit<<" + "<<m_pspect<<" --> \n"
  //	   <<kinargs.m_pi<<" + "<<kinargs.m_pj<<" + "<<kinargs.m_pk<<"\n";
  return true;
}

void FI_CS_Kinematics::CalculateJacobean(Splitting & split,
					 Configuration & config) {
  if (split.T()<0.) m_weight = 1.;
  else {
    m_weight = 0.;
    const double newscale = split.T(), oldscale = split.T();
    const double newpdf   = p_kernel->GetXPDF(newscale,split.Eta()/(1.-split.Y()),
					      split.GetSpectator()->Flav(),
					      split.Beam());
    const double oldpdf   = p_kernel->GetXPDF(oldscale,split.Eta(),
					      split.GetSpectator()->Flav(),
					      split.Beam());
    if (newpdf > 0. && oldpdf > 0. &&
	std::abs(oldpdf) > p_kernel->DynamicPDFThreshold(split.Eta()))
      m_weight = (1.-split.Y()) * newpdf/oldpdf;
  }
}

bool FI_CS_Kinematics::UpdateSystem(Splitting & split,Configuration & config) {
  Vec4D check = ( (split.GetSpectator()->Mom()-split.GetSplitter()->Mom()) -
		  (split.Mom(2)-split.Mom(0)+split.Mom(1)) );
  if (dabs(check.Abs2())>1.e-6 || dabs(check[0])>1.e-3) {
    msg_Error()<<"\n"<<METHOD<<" throws momentum conservation error: "
	       <<"check = "<<check<<"\n";
    return false;
  }
  split.GetSpectator()->SetMom(split.Mom(2));
  return true;
}



DECLARE_GETTER(FI_CS_Kinematics,"Kinematics_FI_Catani-Seymour",
	       Kinematics_Base,Kernel_Info);

Kinematics_Base * ATOOLS::Getter<Kinematics_Base,Kernel_Info,FI_CS_Kinematics>::
operator()(const Parameter_Type & info) const
{
  if (info.Type()==kernel_type::FI &&
      info.KinType()==kin_type::CataniSeymour) {
    return new FI_CS_Kinematics(info);
  }
  return NULL;
}

void ATOOLS::Getter<Kinematics_Base,Kernel_Info,FI_CS_Kinematics>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"Catani-Seymour Kinematics (FI)";
}


