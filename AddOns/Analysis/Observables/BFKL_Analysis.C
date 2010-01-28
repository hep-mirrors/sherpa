#include "AddOns/Analysis/Observables/Normalized_Observable.H"

namespace ANALYSIS {

  class BFKL_Analysis: public Primitive_Observable_Base {  
  private:

    std::string m_reflist, m_jetslist;

    std::vector<Normalized_Observable*> m_dists;
    std::vector<ATOOLS::Histogram*>     m_histos;

  public:

    BFKL_Analysis(const std::string &listname,
		  const std::string &jetsname,
		  const std::string &refname);

    ~BFKL_Analysis();
    
    Analysis_Object &operator+=(const Analysis_Object &obj);
    void EndEvaluation(double scale=1.0);
    void Restore(double scale=1.0);
    void Output(const std::string & pname);
    void Evaluate(const ATOOLS::Particle_List & pl, double weight, double ncount);
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
  Primitive_Observable_Base(1,0.,1.,100)
{
  m_listname=listname;
  m_jetslist=jetsname;
  m_reflist=refname;
  m_name=m_listname+"_JS";
  m_dists.resize(2,NULL);
  // N_{jet}(\Delta y_{1n})
  m_dists[0] = new Normalized_Observable(1,0.0,10.0,40);
  // \cos(\pi-\Delta\phi_{1n})(\Delta y_{1n})
  m_dists[1] = new Normalized_Observable(1,0.0,10.0,40);
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

#define AddZero(weight,ncount)			\
  {						\
    for (size_t i(0);i<m_dists.size();++i)	\
      m_dists[i]->Fill(0.0,0.0,0.0,ncount);	\
    for (size_t i(0);i<m_histos.size();++i)	\
      m_histos[i]->Insert(0.0,0.0,ncount);	\
    return;					\
  }

void BFKL_Analysis::Evaluate(const ATOOLS::Particle_List& pl,
			     double weight, double ncount)
{
  Particle_List jets(*p_ana->GetParticleList(m_listname));
  Particle_List alljets(*p_ana->GetParticleList(m_jetslist));
  Particle_List ref(*p_ana->GetParticleList(m_reflist));
  if (jets.size()<2 || ref.size()!=1) AddZero(weight,ncount);
  if (dabs(ref.front()->Momentum().Y())>4.5) AddZero(weight,ncount);
  // sort pt
  std::sort(jets.begin(),jets.end(),JS_Order_PT());
  m_histos[1]->Insert(jets[0]->Momentum().Y()-
		      jets[1]->Momentum().Y(),weight,ncount);
  // sort y
  std::sort(jets.begin(),jets.end(),JS_Order_Y());
  double dy(jets.front()->Momentum().Y()-
	    jets.back()->Momentum().Y());
  double cdphi(cos(M_PI-jets.front()->Momentum().DPhi
		   (jets.back()->Momentum())));
  m_histos[0]->Insert(dy,weight,ncount);
  m_dists[0]->Fill(dy,jets.size(),weight,ncount);
  m_dists[1]->Fill(dy,cdphi,weight,ncount);
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
	m_histos[3+j]->Insert(y,weight,0);
      m_histos[3+j]->Insert(0.0,0.0,ncount);
    }
  }
  else {
    m_histos[2]->Insert(0.0,0.0,ncount);
    m_histos[3]->Insert(0.0,0.0,ncount);
    m_histos[4]->Insert(0.0,0.0,ncount);
    m_histos[5]->Insert(0.0,0.0,ncount);
  }
}

Analysis_Object &BFKL_Analysis::operator+=
(const Analysis_Object &obj)
{
  const BFKL_Analysis *vob((const BFKL_Analysis*)&obj);
  for (size_t i(0);i<m_dists.size();++i) 
    *m_dists[i]+=*vob->m_dists[i];
  for (size_t i(0);i<m_histos.size();++i) 
    *m_histos[i]+=*vob->m_histos[i];
  return *this;
}

void BFKL_Analysis::EndEvaluation(double scale) 
{
  for (size_t i(0);i<m_dists.size();++i) 
    m_dists[i]->EndEvaluation(scale);
  for (size_t i(0);i<m_histos.size();++i) {
    m_histos[i]->Finalize();
    m_histos[i]->Scale(scale);
  }
}

void BFKL_Analysis::Restore(double scale) 
{
  for (size_t i(0);i<m_dists.size();++i) 
    m_dists[i]->Restore(scale);
  for (size_t i(0);i<m_histos.size();++i) {
    m_histos[i]->Scale(scale);
    m_histos[i]->Restore();
  }
}

void BFKL_Analysis::Output(const std::string & pname) 
{
  msg_Debugging()<<METHOD<<"(): {\n";
  for (size_t i(0);i<m_dists.size();++i) {
    m_dists[i]->SetName(m_name+"_f"+ToString(i)+".dat");
    m_dists[i]->Output(pname);
  }
  for (size_t i(0);i<m_histos.size();++i) {
    m_histos[i]->Output
      (pname+"/"+m_name+"_h"+ToString(i)+".dat");
  }
  msg_Debugging()<<"}\n";
}

BFKL_Analysis::~BFKL_Analysis()
{
  while (m_dists.size()) {
    delete m_dists.back();
    m_dists.pop_back();
  }
  while (m_histos.size()) {
    delete m_histos.back();
    m_histos.pop_back();
  }
}

Primitive_Observable_Base *BFKL_Analysis::Copy() const 
{
  return new BFKL_Analysis(m_listname,m_jetslist,m_reflist);
}

DECLARE_ND_GETTER(JS_Getter,"BFKL",
		  Primitive_Observable_Base,Argument_Matrix,false);

Primitive_Observable_Base *
JS_Getter::operator()(const Argument_Matrix &parameters) const
{
  if (parameters.size()==0 || parameters[0].size()<3) return NULL;
  return new BFKL_Analysis
    (parameters[0][0],parameters[0][1],parameters[0][2]);
}

void JS_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"list jetlist reflist"; 
}

