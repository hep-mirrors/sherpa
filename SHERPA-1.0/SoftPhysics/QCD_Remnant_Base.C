#include "QCD_Remnant_Base.H"

#include "Exception.H"
#include "Random.H"

#ifdef DEBUG__QCD_Remnant_Base
#include "Run_Parameter.H"
#define EVENT 0
#endif

using namespace SHERPA;

QCD_Remnant_Base::QCD_Remnant_Base(PDF::ISR_Handler *isrhandler,const unsigned int beam,
				   const double scale,const Remnant_Base::TypeID type):
  Remnant_Base(type,beam),
  m_deltax(0.0125), 
  m_scale(scale),
  m_ecms(sqrt(isrhandler->Pole())),
  m_xscheme(1), 
  m_maxtrials(100)
{
  if (isrhandler==NULL) {
    throw(ATOOLS::Exception(ATOOLS::ex::fatal_error,"QCD remnant needs ISR Handler.",
			    "QCD_Remnant_Base","QCD_Remnant_Base"));
  }
  p_pdfbase=isrhandler->PDF(m_beam)->GetBasicPDF();
}

void QCD_Remnant_Base::Clear()
{
  m_undo.clear();
  Remnant_Base::Clear();
}

bool QCD_Remnant_Base::TestColours(ATOOLS::Particle *particle,int oldc,int newc,
				 bool singlet,bool force,int anti)
{
  if (particle->GetFlow(anti)==oldc && 
      (m_adjusted.find(particle)==m_adjusted.end() || 
       (force && m_singlet.find(particle)==m_singlet.end()))) {
    return true;
  }
  return false;
}

bool QCD_Remnant_Base::AdjustColours(ATOOLS::Particle *particle,int oldc,int newc,
				   bool &singlet,bool force,int anti,bool forward)
{
  if (m_adjusted.find(particle)!=m_adjusted.end()) m_singlet.insert(particle);
  m_adjusted.insert(particle);
#ifdef DEBUG__QCD_Remnant_Base
  if (ATOOLS::rpa.gen.NumberOfDicedEvents()==EVENT) {
    std::cout<<"["<<m_adjusted.size()<<"] ("<<oldc<<" -> "<<newc<<")"<<particle<<std::endl;
  }
#endif
  if (!force && (particle->GetFlow(1)==newc || particle->GetFlow(2)==newc)) {
    ATOOLS::msg.Tracking()<<"QCD_Remnant_Base::AdjustColours(..): "
			  <<"Created colour singlet. Retry."<<std::endl;
#ifdef DEBUG__QCD_Remnant_Base
    if (ATOOLS::rpa.gen.NumberOfDicedEvents()==EVENT) {
      std::cout<<"QCD_Remnant_Base::AdjustColours(..): "
	       <<"Created colour singlet. Retry."<<std::endl;
    }
#endif
    return singlet=true;
  }
  if ((forward && particle->DecayBlob()==NULL) ||
      (!forward && particle->ProductionBlob()==NULL)) {
    return true;
  }
  if (m_adjusted.size()>100) {
    ATOOLS::msg.Error()<<"QCD_Remnant_Base::AdjustColours(..): "
		       <<"Colour nesting is too deep (more than "<<m_adjusted.size()-1
		       <<" levels)."<<std::endl
		       <<"   Cannot adjust colours completely. "
		       <<"Result might be unreliable."<<std::endl;
    return false;
  }
  ATOOLS::Blob *cur=particle->DecayBlob();
  int newanti=anti;
  bool newforward=forward;
  if (!forward) {
    newanti=3-anti;
    newforward=!forward;
    cur=particle->ProductionBlob();
  }
  for (int i=0;i<cur->NOutP();++i) {
    ATOOLS::Particle *help=cur->OutParticle(i);
    if (TestColours(help,oldc,newc,singlet,force,newanti)) { 
      if (!AdjustColours(help,oldc,newc,singlet,force,newanti,newforward)) return false;
      if (!singlet) help->SetFlow(newanti,newc);
      return true;
    }
  }
  newanti=3-newanti;
  newforward=!newforward;
  for (int i=0;i<cur->NInP();++i) {
    ATOOLS::Particle *help=cur->InParticle(i);
    if (TestColours(help,oldc,newc,singlet,force,newanti)) { 
      if (!AdjustColours(help,oldc,newc,singlet,force,newanti,newforward)) return false;
      if (!singlet) help->SetFlow(newanti,newc);
      return true;
    }
  }
  return true;
}

bool QCD_Remnant_Base::AdjustColours(ATOOLS::Particle *particle,int oldc,int newc,
				   bool &singlet,bool force)
{
  m_singlet.clear();
  m_adjusted.clear();
  if (oldc==newc) return true; 
  size_t i=1;
  for (;i<3;++i) if (particle->GetFlow(i)==oldc) break;
  bool result=AdjustColours(particle,oldc,newc,singlet,force,i,true);
  if (result && !singlet) { 
    particle->SetFlow(i,newc);
    m_undo.push_back(std::pair<ATOOLS::Particle*,
		     std::pair<int,int> >(particle,std::pair<int,int>(i,oldc)));
  }
  return result;
}

ATOOLS::Particle *QCD_Remnant_Base::FindConnected(ATOOLS::Particle *particle,bool same,int orig) 
{
  if (!particle->Flav().IsGluon()) {
    if (particle->Flav().IsAnti()^particle->Flav().IsDiQuark()) orig=2; else orig=1;
  }
  int comp=3-orig;
  if (same) comp=orig;
  for (unsigned int set=0;set<2;++set) {
    for (ATOOLS::Particle_Iterator pit=m_parton[set].begin();pit!=m_parton[set].end();++pit) {
      if ((*pit)->GetFlow(comp)==particle->GetFlow(orig)) return *pit;
    }
  }
  return NULL;
}

ATOOLS::Particle *QCD_Remnant_Base::SelectCompanion(ATOOLS::Particle *cur) 
{
  size_t trials=0;
  bool tested=true;
  ATOOLS::Particle *companion=cur;
  while (tested && trials<m_maxtrials/10) {
    ++trials;
    tested=false;
    size_t set=2*(size_t)ATOOLS::ran.Get();
    companion=m_parton[set][(size_t)(ATOOLS::ran.Get()*((int)m_parton[set].size()-1))];
    if (m_companions.find(companion)!=m_companions.end()) tested=true;
  }
  if (trials==m_maxtrials/10) {
    for (size_t i=0;i<3;++i) {
      for (size_t j=0;j<m_parton[i].size();++j) {
	if (m_companions.find(m_parton[i][j])==m_companions.end()) companion=m_parton[i][j];
      }
    }
  }
  m_companions.insert(companion);
  return companion;
}

void QCD_Remnant_Base::UnDo() 
{
  ATOOLS::msg.Tracking()<<"QCD_Remnant_Base::UnDo(): Undoing changes on blob list."<<std::endl;
#ifdef DEBUG__QCD_Remnant_Base
  if (ATOOLS::rpa.gen.NumberOfDicedEvents()==EVENT) {
    std::cout<<"QCD_Remnant_Base::UnDo(): ["<<m_errors<<"] Undoing changes on blob list."<<std::endl;
  }
#endif
  for (int i=(int)m_undo.size()-1;i>=0;--i) {
    bool singlet=false;
    AdjustColours(m_undo[i].first,m_undo[i].first->GetFlow(m_undo[i].second.first),
		  m_undo[i].second.second,singlet,true);
    p_beamblob->RemoveOutParticle(m_undo[i].first);
  }
  m_undo.clear();
  Remnant_Base::UnDo();
}
