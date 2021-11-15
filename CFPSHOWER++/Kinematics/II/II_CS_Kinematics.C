#include "CFPSHOWER++/Kinematics/Kinematics_Base.H"
#include "CFPSHOWER++/Shower/Kernel.H"
#include "PHASIC++/Channels/CSS_Kinematics.H"
#include "ATOOLS/Org/Message.H"

namespace CFPSHOWER {
  class II_CS_Kinematics : public Kinematics_Base {
  private:
    bool   KinCheck(Splitting & split);
    double CalculateY(Splitting & split);
  public:
    II_CS_Kinematics(const Kernel_Info & info);
    ~II_CS_Kinematics() {}
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

II_CS_Kinematics::II_CS_Kinematics(const Kernel_Info & info)  :
  Kinematics_Base(info)
{
  m_name   = std::string("II: CS Kinematics");
  m_scheme = kin_type::CataniSeymour;
}

bool II_CS_Kinematics::Init(Splitting & split, Configuration & config,
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
  return true;
}

bool II_CS_Kinematics::KinCheck(Splitting & split) {
  split.SetY(m_y = CalculateY(split));
  return (m_y>=0.0 && m_y<=1.0);
}

double II_CS_Kinematics::CalculateY(Splitting & split) {
  return ( split.Z()/( split.Q2()-split.Mass2(0)-split.Mass2(1)-m_mspect2) *
	   ((split.T()+split.Mass2(0))/(1.-split.Z()) +
	    (1.-split.Z())*m_msplit2)); 
}

bool II_CS_Kinematics::operator()(Splitting & split, Configuration & config) {
  if (!KinCheck(split)) return false;
  PHASIC::Kin_Args kinargs(split.Y(),split.Z(),split.Phi());
  if (PHASIC::ConstructIIDipole(split.Mass2(0),split.Mass2(1),
				m_msplit2,m_mspect2,
				m_psplit,m_pspect,kinargs) < 0) return false;
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

void II_CS_Kinematics::CalculateJacobean(Splitting & split,
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

bool II_CS_Kinematics::UpdateSystem(Splitting & split,Configuration & config) {
  // need to check in how far there are more boosts here ad/or signs.
  exit(1);
  Vec4D check = ( (split.GetSplitter()->Mom()+split.GetSpectator()->Mom()) -
		  (split.Mom(0)+split.Mom(1)+split.Mom(2)) );
  if (dabs(check.Abs2())>1.e-6 || dabs(check[0])>1.e-3) {
    msg_Error()<<"\n"<<METHOD<<" throws momentum conservation error: "
	       <<"check = "<<check<<"\n";
    return false;
  }
  split.GetSpectator()->SetMom(split.Mom(2));
  return true;
}



DECLARE_GETTER(II_CS_Kinematics,"Kinematics_II_Catani-Seymour",
	       Kinematics_Base,Kernel_Info);

Kinematics_Base * ATOOLS::Getter<Kinematics_Base,Kernel_Info,II_CS_Kinematics>::
operator()(const Parameter_Type & info) const
{
  if (info.Type()==kernel_type::II &&
      info.KinType()==kin_type::CataniSeymour) {
    return new II_CS_Kinematics(info);
  }
  return NULL;
}

void ATOOLS::Getter<Kinematics_Base,Kernel_Info,II_CS_Kinematics>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"Catani-Seymour Kinematics (II)";
}


