#include "AddOns/Analysis/Analyses/Analysis_Base.H"

namespace ANALYSIS {

  class BFKL_Analysis: public Analysis_Base {  
  private:

    std::string m_reflist, m_jetslist;

  public:

    BFKL_Analysis(const std::string &listname,
		  const std::string &jetsname,
		  const std::string &refname);

    void Evaluate(double weight,double ncount,int mode);
    Primitive_Observable_Base * Copy() const;

  };// end of class BFKL_Analysis

}// end of namespace ANALYSIS

#include "AddOns/Analysis/Main/Primitive_Analysis.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Message.H"
#include <algorithm>

using namespace ANALYSIS;
using namespace ATOOLS;

BFKL_Analysis::BFKL_Analysis(const std::string &listname,
			     const std::string &jetsname,
			     const std::string &refname):
  Analysis_Base(m_listname)
{
  m_jetslist=jetsname;
  m_reflist=refname;
  m_name+="_BFKL";
  m_dists.resize(2,NULL);
  // N_{jet}(\Delta y_{1n})
  m_dists[0] = new Normalized_Observable(4,0.0,10.0,40);
  // \cos(\pi-\Delta\phi_{1n})(\Delta y_{1n})
  m_dists[1] = new Normalized_Observable(4,0.0,10.0,40);
  m_histos.resize(6,NULL);
  // d\sigma/d\Delta y_{1n}
  m_histos[0] = new Histogram(1,0.0,10.0,40);
  // d\sigma/d\Delta y_{12}
  m_histos[1] = new Histogram(1,0.0,10.0,40);
  // d\sigma/dp_{T,H}
  m_histos[2] = new Histogram(1,0.0,500.0,50);
  // d\sigma/dy_c
  m_histos[3] = new Histogram(1,0.0,4.0,40);
  // d\sigma/dy_c 30
  m_histos[4] = new Histogram(1,0.0,4.0,40);
  // d\sigma/dy_c 20
  m_histos[5] = new Histogram(1,0.0,4.0,40);
}

class JS_Order_PT {
public:
  bool operator()(const Particle *a,const Particle *b)
  { return a->Momentum().PPerp2()>b->Momentum().PPerp2(); }
};// end of class JS_Order_PT

class JS_Order_Y {
public:
  bool operator()(const Particle *a,const Particle *b)
  { return a->Momentum().Y()>b->Momentum().Y(); }
};// end of class JS_Order_Y

void BFKL_Analysis::Evaluate(double weight, double ncount,int mode)
{
  Particle_List jets(*p_ana->GetParticleList(m_listname));
  Particle_List alljets(*p_ana->GetParticleList(m_jetslist));
  Particle_List ref(*p_ana->GetParticleList(m_reflist));
  if (jets.size()<2 || ref.size()!=1) AddZero(ncount,mode);
  if (dabs(ref.front()->Momentum().Y())>4.5) AddZero(ncount,mode);
  // sort pt
  std::sort(jets.begin(),jets.end(),JS_Order_PT());
  FillHisto(1,jets[0]->Momentum().Y()-
	    jets[1]->Momentum().Y(),weight,ncount,mode);
  // sort y
  std::sort(jets.begin(),jets.end(),JS_Order_Y());
  double dy(jets.front()->Momentum().Y()-
	    jets.back()->Momentum().Y());
  double cdphi(cos(M_PI-jets.front()->Momentum().DPhi
		   (jets.back()->Momentum())));
  FillHisto(0,dy,weight,ncount,mode);
  FillDist(0,dy,jets.size(),weight,ncount,mode);
  FillDist(1,dy,cdphi,weight,ncount,mode);
  if (dy>4.0) {
    m_histos[2]->Insert
      (ref.front()->Momentum().PPerp(),weight,ncount);
    double ymean=0.5*(jets.front()->Momentum().Y()+
		      jets.back()->Momentum().Y());
    double yc[3]={100.0,100.0,100.0};
    for (size_t i(0);i<alljets.size();++i) {
      if (alljets[i]->Momentum()==jets.front()->Momentum() ||
          alljets[i]->Momentum()==jets.back()->Momentum()) continue;
      double pt(alljets[i]->Momentum().PPerp());
      double dy(alljets[i]->Momentum().Y()-ymean);
      if (pt>40.0) yc[0]=Min(yc[0],dabs(dy));
      if (pt>30.0) yc[1]=Min(yc[1],dabs(dy));
      if (pt>20.0) yc[2]=Min(yc[2],dabs(dy));
    }
    for (int j(0);j<3;++j) {
      for (double y(0.05);y<=yc[j];y+=0.1)
	FillHisto(3+j,y,weight,0,mode);
      FillHisto(3+j,0.0,0.0,ncount,mode);
    }
  }
  else {
    FillHisto(2,0.0,0.0,ncount,mode);
    FillHisto(3,0.0,0.0,ncount,mode);
    FillHisto(4,0.0,0.0,ncount,mode);
    FillHisto(5,0.0,0.0,ncount,mode);
  }
}

Primitive_Observable_Base *BFKL_Analysis::Copy() const 
{
  return new BFKL_Analysis(m_listname,m_jetslist,m_reflist);
}

DECLARE_ND_GETTER(BFKL_Analysis_Getter,"BFKL",
		  Primitive_Observable_Base,Argument_Matrix,false);

Primitive_Observable_Base *
BFKL_Analysis_Getter::operator()(const Argument_Matrix &parameters) const
{
  if (parameters.size()==0 || parameters[0].size()<3) return NULL;
  return new BFKL_Analysis
    (parameters[0][0],parameters[0][1],parameters[0][2]);
}

void BFKL_Analysis_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"list jetlist reflist"; 
}

