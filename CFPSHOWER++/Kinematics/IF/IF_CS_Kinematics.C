#include "CFPSHOWER++/Kinematics/Kinematics_Base.H"
#include "CFPSHOWER++/Shower/Kernel.H"
#include "PHASIC++/Channels/CSS_Kinematics.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"

namespace CFPSHOWER {
  class IF_CS_Kinematics : public Kinematics_Base {
  private:
    double        m_otherISmass2;
    ATOOLS::Vec4D m_otherISmom;
    bool   KinCheck(Splitting & split);
    double CalculateY(Splitting & split);
  public:
    IF_CS_Kinematics(const Kernel_Info & info);
    ~IF_CS_Kinematics() {}
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

IF_CS_Kinematics::IF_CS_Kinematics(const Kernel_Info & info)  :
  Kinematics_Base(info)
{
  m_name   = std::string("IF: CS Kinematics");
  m_scheme = kin_type::CataniSeymour;
}

bool IF_CS_Kinematics::Init(Splitting & split, Configuration & config,
			    const ATOOLS::Mass_Selector * msel) {
  size_t beam    = split.GetSplitter()->Beam();
  PDF_Base * pdf = p_kernel->GetPDF(beam);
  if (beam<1 || beam>2 || pdf==NULL) {
    msg_Error()<<"Error in "<<METHOD<<"(beam = "<<beam<<"):\n"
	       <<"   We should not arrive here.  Will exit for debugging reasons.\n";
    exit(1);
    return false;
  }
  double xmin = Max(1.e-6,pdf->XMin());
  double xmax = Min(0.999999,pdf->XMax());
  double x    = split.GetSplitter()->XB();
  if (x>=xmax || x<=xmin || split.Q2()>=pdf->Q2Max() || split.Q2()<=pdf->Q2Min() ||
      split.Tcut()>split.Q2()) return false;
  double zmin = x/xmax, zmax = split.Q2()/(split.Q2()+split.Tcut());
  if (zmin>zmax) return false;
  split.SetZmin(zmin);
  split.SetZmax(zmax);
  InitSystem(split,msel);
  for (Parton_List::iterator pit=config.begin();pit!=config.end();pit++) {
    if ((*pit)->On() && (*pit)->Beam()>0 && (*pit)!=split.GetSplitter()) {
      m_otherISmass2 = sqr(msel->Mass((*pit)->Flav().Mass()));
      m_otherISmom   = (*pit)->Mom();
      return true;
    }
  }
  THROW(fatal_error,"No 2nd IS parton found.");
  return false;
}

bool IF_CS_Kinematics::KinCheck(Splitting & split) {
  split.SetY(m_y = CalculateY(split));
  return (split.Z()>=0.0 && split.Z()<=1.0 &&
	  split.Q2()>(split.Mass2(0)+split.Mass2(1)+m_mspect2));
}

double IF_CS_Kinematics::CalculateY(Splitting & split) {
  return -(split.Z()/( split.Q2()-split.Mass2(0)-split.Mass2(1)-m_mspect2) *
	   ((split.T()+split.Mass2(0))/(1.-split.Z()) +
	    (1.-split.Z())*m_msplit2)); 
}

bool IF_CS_Kinematics::operator()(Splitting & split, Configuration & config) {
  if (!KinCheck(split)) return false;
  PHASIC::Kin_Args kinargs(split.Y(),split.Z(),split.Phi());
  if (PHASIC::ConstructIFDipole(split.Mass2(0),split.Mass2(1),
				m_msplit2,m_mspect2,m_otherISmass2,
				m_psplit,m_pspect,m_otherISmom,kinargs) < 0) return false;
  split.SetPoincare(kinargs.m_lam);
  kinargs.m_lam.Invert();
  split.SetMom(0,kinargs.m_lam * kinargs.m_pi);
  split.SetMom(1,kinargs.m_lam * kinargs.m_pj);
  split.SetMom(2,kinargs.m_lam * kinargs.m_pk);
  split.SetKinSpect(kinargs.m_lam * kinargs.m_pk);
  //msg_Out()<<METHOD
  //	   <<"(z = "<<split.Z()<<" vs "<<m_z[0]<<" from "<<split.Ztest()<<", "
  //	   <<"y = "<<m_y<<", phi = "<<split.Phi()<<"):\n"
  //	   <<m_psplit<<" + "<<m_pspect<<" --> \n"
  //	   <<kinargs.m_pi<<" + "<<kinargs.m_pj<<" + "<<kinargs.m_pk<<"\n";
  return true;
}

void IF_CS_Kinematics::CalculateJacobean(Splitting & split,
					 Configuration & config) {
  if (split.T()<0.) m_weight = 1./split.Z();
  else {
    m_weight = 0.;
    const double newscale = split.T(), oldscale = split.T();
    const double newpdf   = p_kernel->GetXPDF(newscale,split.Eta()/split.Z(),
					      split.GetKernel()->GetFlavs()[0],
					      split.Beam());
    const double oldpdf   = p_kernel->GetXPDF(oldscale,split.Eta(),
					      split.GetKernel()->GetSplit(),
					      split.Beam());
    if (newpdf > 0. && oldpdf > 0. &&
	std::abs(oldpdf) > p_kernel->DynamicPDFThreshold(split.Eta()))
      m_weight = newpdf/oldpdf;
  }
}

bool IF_CS_Kinematics::UpdateSystem(Splitting & split,Configuration & config) {
  // need to check in how far there are more boosts here.
  exit(1);
  Vec4D check = ( (split.GetSplitter()->Mom()-split.GetSpectator()->Mom()) -
		  (split.Mom(0)-split.Mom(1)+split.Mom(2)) );
  if (dabs(check.Abs2())>1.e-6 || dabs(check[0])>1.e-3) {
    msg_Error()<<"\n"<<METHOD<<" throws momentum conservation error: "
	       <<"check = "<<check<<"\n";
    return false;
  }
  split.GetSpectator()->SetMom(split.Mom(2));
  return true;
}



DECLARE_GETTER(IF_CS_Kinematics,"Kinematics_IF_Catani-Seymour",
	       Kinematics_Base,Kernel_Info);

Kinematics_Base * ATOOLS::Getter<Kinematics_Base,Kernel_Info,IF_CS_Kinematics>::
operator()(const Parameter_Type & info) const
{
  if (info.Type()==kernel_type::IF &&
      info.KinType()==kin_type::CataniSeymour) {
    return new IF_CS_Kinematics(info);
  }
  return NULL;
}

void ATOOLS::Getter<Kinematics_Base,Kernel_Info,IF_CS_Kinematics>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"Catani-Seymour Kinematics (IF)";
}


