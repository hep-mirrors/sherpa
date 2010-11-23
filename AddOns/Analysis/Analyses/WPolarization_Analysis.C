#include "AddOns/Analysis/Analyses/Analysis_Base.H"

namespace ANALYSIS {

  class WPolarization_Analysis: public Analysis_Base {  
  public:

    WPolarization_Analysis(const std::string &listname);

    void Evaluate(double weight, double ncount,int mode);
    Primitive_Observable_Base * Copy() const;

  };// end of class WPolarization_Analysis

}// end of namespace ANALYSIS

#include "AddOns/Analysis/Main/Primitive_Analysis.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include <algorithm>

using namespace ANALYSIS;
using namespace ATOOLS;

WPolarization_Analysis::WPolarization_Analysis(const std::string &listname):
  Analysis_Base(listname)
{
  m_name+="_WPolarization";
  // m_dists.resize(2,NULL);
  // m_dists[0] = new Normalized_Observable(1,0.0,10.0,10);
  m_histos.resize(2,NULL);
  // d\sigma/d\cos(theta^*)
  m_histos[0] = new Histogram(1,-1.0,1.0,100,"CosThetaStar");
  // d\sigma/d\phi^*
  m_histos[1] = new Histogram(1,0.0,360.0,90,"PhiStar");
}

void WPolarization_Analysis::Evaluate(double weight,double ncount,int mode)
{
  DEBUG_FUNC("");
  Particle_List all(*p_ana->GetParticleList(m_listname));
  Particle *l(NULL), *nu(NULL);
  for (size_t i(0);i<all.size();++i) {
    if (Flavour(kf_lepton).Includes(all[i]->Flav())) {
      if (l) THROW(fatal_error,"More than one lepton found");
      l=all[i];
    }
    if (Flavour(kf_neutrino).Includes(all[i]->Flav())) {
      if (nu) THROW(fatal_error,"More than one neutrino found");
      nu=all[i];
    }
  }
  if (l==NULL || nu==NULL) AddZero(weight,ncount);
  Vec4D pb1(rpa.gen.PBeam(0)), pb2(rpa.gen.PBeam(1));
  Vec4D pl(l->Momentum()), plnu(pl+nu->Momentum());
  Poincare cms(plnu), zrot(plnu,Vec4D::ZVEC);
  cms.Boost(pl);
  cms.Boost(pb1);
  cms.Boost(pb2);
  zrot.Rotate(pl);
  zrot.Rotate(pb1);
  zrot.Rotate(pb2);
  if (pb1.PPerp2()==0.0 && pb2.PPerp2()==0.0) {
    msg_Debugging()<<"using arbitrary x axis\n";
  }
  else {
    Vec4D xref(pb1.PPerp2()/pb1.PSpat2()>
	       pb2.PPerp2()/pb2.PSpat2()?pb1:pb2);
    Poincare xrot(xref.Perp(),Vec4D::XVEC);
    xrot.Rotate(pl);
    xrot.Rotate(pb1);
    xrot.Rotate(pb2);
    if (pb2.Theta(-Vec4D::XVEC)<pb1.Theta(Vec4D::XVEC))
      msg_Error()<<METHOD<<"(): \\theta_1 = "<<pb1.Theta(Vec4D::XVEC)
		 <<" > \\theta_2 = "<<pb2.Theta(-Vec4D::XVEC)<<std::endl;
  }
  double costhetas(pl.CosTheta()), phis(pl.Phi());
  if (phis<0.0) phis+=2.0*M_PI;
  msg_Debugging()<<"cos\\theta^* = "<<costhetas
		 <<", \\phi^* = "<<phis*180.0/M_PI<<"\n";
  FillHisto(0,costhetas,weight,ncount,mode);
  FillHisto(1,phis*180.0/M_PI,weight,ncount,mode);
}

Primitive_Observable_Base *WPolarization_Analysis::Copy() const 
{
  return new WPolarization_Analysis(m_listname);
}

DECLARE_GETTER(WPolarization_Getter,"WPolarization",
	       Primitive_Observable_Base,Argument_Matrix);

Primitive_Observable_Base *
WPolarization_Getter::operator()(const Argument_Matrix &parameters) const
{
  if (parameters.size()==0 || parameters[0].size()<1) return NULL;
  return new WPolarization_Analysis(parameters[0][0]);
}

void WPolarization_Getter::PrintInfo
(std::ostream &str,const size_t width) const
{ 
  str<<"list"; 
}

