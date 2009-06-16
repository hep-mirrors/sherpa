#include "AddOns/Analysis/Observables/Normalized_Observable.H"

namespace ANALYSIS {

  class BFKL_Analysis: public Primitive_Observable_Base {  
  private:

    std::string m_reflist;

    std::vector<Normalized_Observable*> m_dists;
    std::vector<ATOOLS::Histogram*>     m_histos;

    double m_dy, m_ptmin;

  public:

    BFKL_Analysis(const std::string & listname,
		  const double &dy,const double &ptmin);
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

BFKL_Analysis::BFKL_Analysis
(const std::string & listname,
 const double &dy,const double &ptmin) :
  Primitive_Observable_Base(1,0.,1.,100),
  m_dy(dy), m_ptmin(ptmin)
{
  m_listname=listname;
  m_name=m_listname+"_JS_Dy"+ToString(m_dy)
    +"_pT"+ToString(m_ptmin);
  m_dists.resize(9,NULL);
  // N_{jet}(\Delta y_{1n})
  m_dists[0] = new Normalized_Observable(1,0.0,10.0,100);
  // \Delta\phi_{1n}(\Delta y_{1n})
  m_dists[1] = new Normalized_Observable(1,0.0,10.0,100);
  // 2-jets(\Delta y_{1n})
  m_dists[2] = new Normalized_Observable(1,0.0,10.0,100);
  // 3-jets(\Delta y_{1n})
  m_dists[3] = new Normalized_Observable(1,0.0,10.0,100);
  // 4-jets(\Delta y_{1n})
  m_dists[4] = new Normalized_Observable(1,0.0,10.0,100);
  // 5-jets(\Delta y_{1n})
  m_dists[5] = new Normalized_Observable(1,0.0,10.0,100);
  // p_{T,mid lead}(\Delta y_{1n})
  m_dists[6] = new Normalized_Observable(1,0.0,10.0,100);
  // H_T(\Delta y_{1n})
  m_dists[7] = new Normalized_Observable(1,0.0,10.0,100);
  // H_T/N_{jet}(\Delta y_{1n})
  m_dists[8] = new Normalized_Observable(1,0.0,10.0,100);
}

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
  for (Particle_List::iterator 
	 pit(jets.begin());pit!=jets.end();) {
    // jets only
    if ((*pit)->Flav().Kfcode()!=kf_jet) pit=jets.erase(pit);
    // y cut
    else if (dabs((*pit)->Momentum().Y())>m_dy) pit=jets.erase(pit);
    // pt cut
    else if ((*pit)->Momentum().PPerp()<m_ptmin) pit=jets.erase(pit);
    else ++pit;
  }
  if (jets.size()<2) AddZero(weight,ncount);
  // sort
  std::sort(jets.begin(),jets.end(),JS_Order_Y());
  double dy(jets.front()->Momentum().Y()-
	    jets.back()->Momentum().Y());
  double cdphi(cos(M_PI-jets.front()->Momentum().DPhi
		   (jets.back()->Momentum())));
  m_dists[0]->Fill(dy,jets.size(),weight,ncount);
  m_dists[1]->Fill(dy,cdphi,weight,ncount);
  for (size_t i(2);i<6;++i) {
    if (i==jets.size()) m_dists[i]->Fill(dy,1.0,weight,ncount);
    else m_dists[i]->Fill(dy,0.0,weight,ncount);
  }
  double ptml(0.0), ht(0.0);
  for (size_t i(1);i<jets.size()-1;++i) {
    double pt(jets[i]->Momentum().PPerp());
    if (pt>ptml) ptml=pt;
    ht+=pt;
  }
  ht+=jets.front()->Momentum().PPerp();
  ht+=jets.back()->Momentum().PPerp();
  if (ptml>0.0) m_dists[6]->Fill(dy,ptml,weight,ncount);
  else m_dists[6]->Fill(dy,0.0,0.0,ncount);
  m_dists[7]->Fill(dy,ht,weight,ncount);
  m_dists[8]->Fill(dy,ht/jets.size(),weight,ncount);
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
  return new BFKL_Analysis(m_listname,m_dy,m_ptmin);
}

DECLARE_ND_GETTER(JS_Getter,"BFKL",
		  Primitive_Observable_Base,Argument_Matrix,false);

Primitive_Observable_Base *
JS_Getter::operator()(const Argument_Matrix &parameters) const
{
  if (parameters.size()==0 || parameters[0].size()<3) return NULL;
  return new BFKL_Analysis(parameters[0][2],
			   ToType<double>(parameters[0][0]),
			   ToType<double>(parameters[0][1]));
}

void JS_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"deltay ptmin list"; 
}

