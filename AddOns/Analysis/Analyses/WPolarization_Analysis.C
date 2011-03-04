#include "AddOns/Analysis/Analyses/Analysis_Base.H"

namespace ANALYSIS {

  class WPolarization_Analysis: public Analysis_Base {
  private:

    Argument_Matrix m_params;

    double m_ptwc;

  public:

    WPolarization_Analysis(const Argument_Matrix &params);

    void Evaluate(double weight, double ncount,int mode);
    Primitive_Observable_Base * Copy() const;

  };// end of class WPolarization_Analysis

}// end of namespace ANALYSIS

#include "AddOns/Analysis/Main/Primitive_Analysis.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include <algorithm>

using namespace ANALYSIS;
using namespace ATOOLS;

WPolarization_Analysis::WPolarization_Analysis(const Argument_Matrix &params):
  Analysis_Base(params[0][0]), m_params(params), m_ptwc(20.0)
{
  Data_Reader reader(",",";","!","=");
  Algebra_Interpreter *ip=reader.Interpreter();
  for (size_t i(1);i<params.size();++i) {
    if (params.size()<2) continue;
    if (params[i][0]=="PT_W_Min") ToType<double>(ip->Interprete(params[i][1]));
  }
  m_name+="_WPolarization";
  m_dists.resize(11,NULL);
  m_dists[0] = new Normalized_Observable(1,0.0,1000.0,100,"A0",1);
  m_dists[1] = new Normalized_Observable(1,0.0,1000.0,100,"A1",1);
  m_dists[2] = new Normalized_Observable(1,0.0,1000.0,100,"A2",1);
  m_dists[3] = new Normalized_Observable(1,0.0,1000.0,100,"A3",1);
  m_dists[4] = new Normalized_Observable(1,0.0,1000.0,100,"A4",1);
  m_dists[5] = new Normalized_Observable(1,0.0,1000.0,100,"A5",1);
  m_dists[6] = new Normalized_Observable(1,0.0,1000.0,100,"A6",1);
  m_dists[7] = new Normalized_Observable(1,0.0,1000.0,100,"A7",1);
  m_dists[8] = new Normalized_Observable(1,0.0,1000.0,100,"fL",1);
  m_dists[9] = new Normalized_Observable(1,0.0,1000.0,100,"fR",1);
  m_dists[10] = new Normalized_Observable(1,0.0,1000.0,100,"f0",1);
  m_histos.resize(6,NULL);
  m_histos[0] = new Histogram(1,-1.0,1.0,100,"CosThetaStar");
  m_histos[1] = new Histogram(1,0.0,360.0,90,"PhiStar");
  m_histos[2] = new Histogram(1,-1.0,1.0,100,"CosThetaStar_PTW"+ToString(m_ptwc));
  m_histos[3] = new Histogram(1,0.0,360.0,90,"PhiStar_PTW"+ToString(m_ptwc));
  m_histos[4] = new Histogram(1,0.0,1000.0,100,"PTW");
  m_histos[5] = new Histogram(10,0.1,1000.0,100,"logPTW");
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
  if (l==NULL || nu==NULL) AddZero(ncount,mode);
  Vec4D pb1(rpa.gen.PBeam(0)), pb2(rpa.gen.PBeam(1));
  Vec4D pl(l->Momentum()), plnu(pl+nu->Momentum());
  Poincare cms(plnu), zrot(plnu,Vec4D::ZVEC);
  cms.Boost(pl);
  cms.Boost(pb1);
  cms.Boost(pb2);
  zrot.Rotate(pl);
  zrot.Rotate(pb1);
  zrot.Rotate(pb2);
  Vec4D xref(pb1.CosTheta()>pb2.CosTheta()?pb1:pb2);
  Poincare xrot(xref.Perp(),Vec4D::XVEC);
  xrot.Rotate(pl);
  double ptw((l->Momentum()+nu->Momentum()).PPerp());
  FillHisto(4,ptw,weight,ncount,mode);
  FillHisto(5,ptw,weight,ncount,mode);
  double thetas(pl.Theta()), phis(pl.Phi());
  double costhetas(cos(thetas)), sinthetas(sin(thetas));
  double cosphis(cos(phis)), sinphis(sin(phis));
  if (phis<0.0) phis+=2.0*M_PI;
  msg_Debugging()<<"cos\\theta^* = "<<costhetas
		 <<", \\phi^* = "<<phis*180.0/M_PI<<"\n";
  FillHisto(0,costhetas,weight,ncount,mode);
  FillHisto(1,phis*180.0/M_PI,weight,ncount,mode);
  double ptcweight(ptw>m_ptwc?weight:0.0);
  FillHisto(2,costhetas,ptcweight,ncount,mode);
  FillHisto(3,phis*180.0/M_PI,ptcweight,ncount,mode);
  FillDist(0,ptw,10.0/3.0*(1.0-3.0*sqr(costhetas))+2.0/3.0,weight,ncount,mode);
  FillDist(1,ptw,10.0*sinthetas*costhetas*cosphis,weight,ncount,mode);
  FillDist(2,ptw,10.0*sqr(sinthetas)*(sqr(cosphis)-sqr(sinphis)),weight,ncount,mode);
  FillDist(3,ptw,4.0*sinthetas*cosphis,weight,ncount,mode);
  FillDist(4,ptw,4.0*costhetas,weight,ncount,mode);
  FillDist(5,ptw,4.0*sinthetas*sinphis,weight,ncount,mode);
  FillDist(6,ptw,10.0*costhetas*sinthetas*sinphis,weight,ncount,mode);
  FillDist(7,ptw,10.0*sqr(sinthetas)*cosphis*sinphis,weight,ncount,mode);
  FillDist(8,ptw,0.5*sqr(1.0-costhetas)-(1.0-2.0*sqr(costhetas)),weight,ncount,mode);
  FillDist(9,ptw,0.5*sqr(1.0+costhetas)-(1.0-2.0*sqr(costhetas)),weight,ncount,mode);
  FillDist(10,ptw,5.0*sqr(sinthetas)-3.0,weight,ncount,mode);
}

Primitive_Observable_Base *WPolarization_Analysis::Copy() const 
{
  return new WPolarization_Analysis(m_params);
}

namespace ANALYSIS {

DECLARE_GETTER(WPolarization_Getter,"WPolarization",
	       Primitive_Observable_Base,Argument_Matrix);

Primitive_Observable_Base *
WPolarization_Getter::operator()(const Argument_Matrix &parameters) const
{
  if (parameters.size()==0 || parameters[0].size()<1) return NULL;
  return new WPolarization_Analysis(parameters);
}

void WPolarization_Getter::PrintInfo
(std::ostream &str,const size_t width) const
{ 
  str<<"list"; 
}

}
