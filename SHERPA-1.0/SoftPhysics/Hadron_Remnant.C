#include "Hadron_Remnant.H"

#include "Exception.H"
#include "Random.H"

#ifdef DEBUG__Hadron_Remnant
#include "Run_Parameter.H"
#define EVENT 753
#endif

using namespace SHERPA;

Hadron_Remnant::Hadron_Remnant(PDF::ISR_Handler *isrhandler,const unsigned int _m_beam,double _m_scale):
  Remnant_Base(Hadron_Remnant::Hadron,_m_beam),
  m_deltax(0.0125), 
  m_scale(-_m_scale),
  m_ecms(sqrt(isrhandler->Pole())),
  m_xscheme(1), 
  m_maxtrials(100)
{
  if (isrhandler==NULL) {
    throw(ATOOLS::Exception(ATOOLS::ex::fatal_error,"Hadron remnant needs ISR Handler.",
			    "Hadron_Remnant","Hadron_Remnant"));
  }
  p_pdfbase=isrhandler->PDF(m_beam)->GetBasicPDF();
  GetConstituents(isrhandler->Flav(m_beam));
}

void Hadron_Remnant::GetConstituents(const ATOOLS::Flavour flav) 
{
  int hadint=(flav.Kfcode()-flav.Kfcode()/10000)/10;
  if ((hadint>100)&&(hadint<1000)) {
    m_constit.resize(3);
    m_constit[0]=ATOOLS::Flavour(ATOOLS::kf::code(hadint)/100);
    m_constit[1]=ATOOLS::Flavour(ATOOLS::kf::code((hadint-(hadint/100)*100)/10));
    m_constit[2]=ATOOLS::Flavour(ATOOLS::kf::code(hadint-(hadint/10)*10));
    if (flav.IsAnti()) {
      for(int i=0;i<3;i++) m_constit[i]=m_constit[i].Bar();
    }
    ATOOLS::msg.Tracking()<<"Hadron_Remnant::FindConstituents("<<flav<<"): "
			  <<"Hadron is baryon."<<std::endl<<"   Constituents are ["
			  <<m_constit[0]<<","<<m_constit[1]<<","<<m_constit[2]<<"]."<<std::endl;
    return;
  }
  if ((hadint>10)&&(hadint<100)) {
    m_constit.resize(2);
    m_constit[0]=ATOOLS::Flavour(ATOOLS::kf::code(hadint)/10);
    m_constit[1]=ATOOLS::Flavour(ATOOLS::kf::code(hadint-(hadint/10)*10));
    if (flav.IsAnti()) {
      for(int i=0;i<2;i++) m_constit[i]=m_constit[i].Bar();
    }
    ATOOLS::msg.Tracking()<<"Hadron_Remnant::FindConstituents("<<flav<<"): "
			  <<"Hadron is meson."<<std::endl<<"   Constituents are ["
			  <<m_constit[0]<<","<<m_constit[1]<<"]."<<std::endl;
    return;
  }
  throw(ATOOLS::Exception(ATOOLS::ex::critical_error,"Cannot determine constituents.",
			  "Hadron_Remnant","GetConstituents"));
}

bool Hadron_Remnant::FillBlob(ATOOLS::Blob *beamblob,ATOOLS::Particle_List *particlelist)
{
  if (p_partner==NULL) {
    throw(ATOOLS::Exception(ATOOLS::ex::critical_error,"Partner Remnant not set.",
			    "Hadron_Remnant","GetConstituents"));
  }
  p_beamblob=beamblob;
  m_pbeam=beamblob->InParticle(0)->Momentum();
  m_undo.clear();
  // decompose hadron
  m_hardpt=ATOOLS::Vec4D();
  for (size_t i=0;i<m_parton[1].size();++i) {
    ATOOLS::Particle *cur=m_parton[1][i];
    m_hardpt+=cur->Momentum();
    if (i>0) {
      if (cur->Flav().Kfcode()==ATOOLS::kf::gluon) {
	if (!TreatGluon(cur)) {
	  UnDo();
	  p_partner->UnDo();
	  return false;
	}
      }
      else {
	if (!TreatQuark(cur)) {
	  UnDo();
	  p_partner->UnDo();
	  return false;
	}
      }
    }
    else {
      if (cur->Flav().Kfcode()==ATOOLS::kf::gluon) TreatFirstGluon(cur);
      else TreatFirstQuark(cur);
    }
#ifdef DEBUG__Hadron_Remnant
    if (ATOOLS::rpa.gen.NumberOfDicedEvents()==EVENT) {
      std::cout<<"particle ["<<i<<"]("<<m_parton[1].size()<<") {\n";
      for (size_t j=0;j<m_parton[0].size();++j) std::cout<<m_parton[0][j]<<std::endl;
      std::cout<<"\n";
      for (size_t j=0;j<m_parton[2].size();++j) std::cout<<m_parton[2][j]<<std::endl;
      std::cout<<"}\n";
    }
#endif
  }
  // select x's according to pdf
  DiceKinematics();
  // fill blob
  for (int i=1;i>=0;--i) {
    for (size_t j=0;j<m_parton[i].size();++j) {
      if (i==0) {
	m_parton[i][j]->SetNumber((long int)m_parton[i][j]);
	m_parton[i][j]->SetInfo('F');
      }
      // if (i==1) beamblob->AddToOutParticles(m_parton[i][j]);
      if (particlelist!=NULL) {
	m_parton[i][j]->SetNumber(particlelist->size());
	particlelist->push_back(m_parton[i][j]);
      }
    }
  }
  return true;
}

void Hadron_Remnant::DiceKinematics()
{
  unsigned int trials;
  ATOOLS::Vec4D ptot=m_pbeam;
  double m_xtot=1.0;
  for (unsigned int i=0;i<m_parton[1].size();++i) m_xtot-=m_parton[1][i]->Momentum()[0]/ptot[0];
  trials=0;
  m_xscheme=1;
  double xtot;
  std::map<ATOOLS::Particle*,double> xmap;
  do {
    ++trials;
    xtot=m_xtot;
    p_pdfbase->Reset();
    for(unsigned int i=0;i<m_parton[1].size();++i) {
      p_pdfbase->Extract(m_parton[1][i]->Flav(),
			 2.*m_parton[1][i]->Momentum()[0]/m_ecms);
    }
    for (unsigned int i=0;i<m_parton[0].size();++i) {
      if (!m_parton[0][i]->Flav().IsDiQuark()) {
	double value=1.0;
	for (unsigned int j=0;(xtot-value<m_deltax)&&(j<m_maxtrials/10);++j) { 
	  value=GetXPDF(m_parton[0][i]->Flav(),m_scale); 
      	}
	p_pdfbase->Extract(m_parton[0][i]->Flav(),value);
	xmap[m_parton[0][i]]=value;
 	xtot-=value;
      }
      else {
	p_last[0]=m_parton[0][i];
      }
    }
    xmap[p_last[0]]=xtot;
    if (trials>m_maxtrials) {
      ATOOLS::msg.Tracking()<<"Hadron_Remnant::DiceKinematics(): "
			    <<"Too many trials to find appropriate x values for partons."<<std::endl
			    <<"   Using naive distribution instead."<<std::endl;
      m_xscheme=0;
    }
  } while ((xtot<m_deltax)&&(m_xscheme!=0)&&(xmap[p_last[0]]<p_last[0]->Flav().PSMass()));
  p_pdfbase->Reset();
  if (m_xscheme==0) {
    xtot=0.;
    for (std::map<ATOOLS::Particle*,double>::iterator it=xmap.begin();it!=xmap.end(); ++it) {
      xtot+=it->second=it->first->Flav().PSMass()/m_pbeam[0];
    }
    for (std::map<ATOOLS::Particle*,double>::iterator it=xmap.begin();it!=xmap.end(); ++it) {
      it->second*=m_xtot/xtot;
    }
  }
  for (unsigned int i=0;i<m_parton[1].size();++i) ptot-=m_parton[1][i]->Momentum();
  for (unsigned int j=0;j<m_parton[0].size();++j) {
    double E=xmap[m_parton[0][j]]*m_pbeam[0];
    double m=m_parton[0][j]->Flav().PSMass();
    // crude, to be changed soon ...
    if (m>E) m=0.;
    double pz=ATOOLS::Sign(m_pbeam[3])*sqrt(E*E-ATOOLS::sqr(m)-m_hardpt.PPerp2()/ATOOLS::sqr(m_parton[0].size()));
    m_parton[0][j]->SetMomentum(ATOOLS::Vec4D(E,-m_hardpt[1]/m_parton[0].size(),
					      -m_hardpt[2]/m_parton[0].size(),pz));
    // the brackets are necessary for 'nan'-values
    if (!(E>0.) || (!(pz>0.) && !(pz<=0.))) {
      ATOOLS::msg.Error()<<"Hadron_Remnant::DiceKinematics(): "                 
                         <<"Parton ("<<(long int)m_parton[0][j]<<") has non-positive momentum: p = "
                         <<m_parton[0][j]->Momentum()<<" m_{"<<m_parton[0][j]->Flav()<<"} = "
                         <<m_parton[0][j]->Flav().PSMass()<<" <- "<<m_xscheme<<std::endl;
    }
  }
}

bool Hadron_Remnant::TreatFirstGluon(ATOOLS::Particle *cur) 
{
  size_t single=(size_t)(ATOOLS::ran.Get()*3.); 
  ATOOLS::Flavour difl, fl=m_constit[single];
  int di[2];
  for (unsigned int j=0, i=0;i<3;i++) {
    if (i!=single) di[j++]=m_constit[i].Kfcode();
  }
  if (di[0]!=di[1]) {
    if (ATOOLS::ran.Get()<0.25) difl=ATOOLS::Flavour(ATOOLS::kf::code(abs(di[0])*1000+abs(di[1])*100+1));
    else difl=ATOOLS::Flavour(ATOOLS::kf::code(abs(di[0])*1000+abs(di[1])*100+3));
  }
  else {
    difl=ATOOLS::Flavour(ATOOLS::kf::code(abs(di[0])*1100+3));
  }
  if (m_constit[0].IsAnti()) difl=difl.Bar();
  ATOOLS::Particle *newpart[2];
  newpart[0] = new ATOOLS::Particle(-1,difl); 
  newpart[1] = new ATOOLS::Particle(-1,fl); 
  int anti=difl.IsAnti();
  newpart[0]->SetFlow(2-anti,cur->GetFlow(1+anti)); 
  newpart[0]->SetFlow(1+anti,0);
  newpart[1]->SetFlow(1+anti,cur->GetFlow(2-anti));
  newpart[1]->SetFlow(2-anti,0); 
  for (unsigned int i=0;i<2;++i) {
    m_parton[0].push_back(newpart[i]);
    p_beamblob->AddToOutParticles(newpart[i]);
  }
  m_parton[2].push_back(cur);
  p_beamblob->AddToOutParticles(cur);
  p_last[0]=newpart[0];
  return true;
}

bool Hadron_Remnant::TreatFirstQuark(ATOOLS::Particle *cur) 
{
  bool found=false;
  int di[3];
  for (size_t j=0,i=0;i<m_constit.size();++i) {
    if (found||cur->Flav()!=m_constit[i]) di[j++]=m_constit[i];
    else found=true;
  }
  if (found) {
    ATOOLS::Flavour difl;
    if (di[0]!=di[1]) {
      if (ATOOLS::ran.Get()<0.25) difl=ATOOLS::Flavour(ATOOLS::kf::code(abs(di[0])*1000+abs(di[1])*100+1));
      else difl=ATOOLS::Flavour(ATOOLS::kf::code(abs(di[0])*1000+abs(di[1])*100+3));
    }
    else {
      difl=ATOOLS::Flavour(ATOOLS::kf::code(abs(di[0])*1100+3));
    }
    if (m_constit[0].IsAnti()) difl=difl.Bar();
    ATOOLS::Particle *newpart = new ATOOLS::Particle(-1,difl); 
    newpart->SetFlow(1,cur->GetFlow(2));
    newpart->SetFlow(2,cur->GetFlow(1));
    m_parton[0].push_back(newpart);
    p_beamblob->AddToOutParticles(newpart);
    p_last[0]=newpart;
  }
  else {
    unsigned int single=(unsigned int)(ATOOLS::ran.Get()*3.0); 
    ATOOLS::Flavour difl, fl=m_constit[single];
    for (unsigned int j=0,i=0;i<3;i++) {
      if (i!=single) di[j++]=m_constit[i].Kfcode();
    }
    if (di[0]!=di[1]) {
      if (ATOOLS::ran.Get()<0.25) difl=ATOOLS::Flavour(ATOOLS::kf::code(di[0]*1000+di[1]*100+1));
      else difl=ATOOLS::Flavour(ATOOLS::kf::code(di[0]*1000+di[1]*100+3));
    }
    else difl=ATOOLS::Flavour(ATOOLS::kf::code(di[0]*1100+3));
    if (m_constit[0].IsAnti()) difl=difl.Bar();
    ATOOLS::Particle *newpart[3];
    newpart[0] = new ATOOLS::Particle(-1,difl); 
    newpart[1] = new ATOOLS::Particle(-1,(cur->Flav()).Bar()); 
    newpart[2] = new ATOOLS::Particle(-1,fl); 
    if (cur->Flav().IsAnti()^difl.IsAnti()) {
      ATOOLS::Particle *rem=newpart[0]; newpart[0]=newpart[2]; newpart[2]=rem;
    }
    int anti=cur->Flav().IsAnti();
    newpart[0]->SetFlow(2-anti,cur->GetFlow(1+anti));
    newpart[1]->SetFlow(2-anti,-1);
    newpart[1]->SetFlow(1+anti,0);
    newpart[2]->SetFlow(1+anti,newpart[1]->GetFlow(2-anti));
    newpart[2]->SetFlow(2-anti,0);
    for (unsigned int i=0;i<3;++i) {
      m_parton[0].push_back(newpart[i]);
      p_beamblob->AddToOutParticles(newpart[i]);
    }
    p_last[0]=newpart[0];
  }
  m_parton[2].push_back(cur);
  p_beamblob->AddToOutParticles(cur);
  return true;
}

bool Hadron_Remnant::TestColours(ATOOLS::Particle *particle,int oldc,int newc,
				 bool singlet,bool force,int anti)
{
  if (particle->GetFlow(anti)==oldc && 
      (m_adjusted.find(particle)==m_adjusted.end() || 
       (force && m_singlet.find(particle)==m_singlet.end()))) {
    return true;
  }
  return false;
}

bool Hadron_Remnant::AdjustColours(ATOOLS::Particle *particle,int oldc,int newc,
				   bool &singlet,bool force,int anti,bool forward)
{
  if (m_adjusted.find(particle)!=m_adjusted.end()) m_singlet.insert(particle);
  m_adjusted.insert(particle);
#ifdef DEBUG__Hadron_Remnant
  if (ATOOLS::rpa.gen.NumberOfDicedEvents()==EVENT) {
    std::cout<<"["<<m_adjusted.size()<<"] ("<<oldc<<" -> "<<newc<<")"<<particle<<std::endl;
  }
#endif
  if (!force && (particle->GetFlow(1)==newc || particle->GetFlow(2)==newc)) {
    ATOOLS::msg.Tracking()<<"Hadron_Remnant::AdjustColours(..): "
			  <<"Created colour singlet. Retry."<<std::endl;
#ifdef DEBUG__Hadron_Remnant
    if (ATOOLS::rpa.gen.NumberOfDicedEvents()==EVENT) {
      std::cout<<"Hadron_Remnant::AdjustColours(..): "
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
    ATOOLS::msg.Error()<<"Hadron_Remnant::AdjustColours(..): "
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

bool Hadron_Remnant::AdjustColours(ATOOLS::Particle *particle,int oldc,int newc,
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

ATOOLS::Particle *Hadron_Remnant::FindConnected(ATOOLS::Particle *particle,bool same,int orig) 
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

ATOOLS::Particle *Hadron_Remnant::SelectCompanion(ATOOLS::Particle *cur) 
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

bool Hadron_Remnant::TreatQuark(ATOOLS::Particle *cur) 
{
  unsigned int trials=0;
  bool success=true, singlet=true;
  ATOOLS::Particle *comp[2];
  ATOOLS::Particle *newpart = new ATOOLS::Particle(-1,cur->Flav().Bar());
  newpart->SetFlow(1+(int)newpart->Flav().IsAnti(),-1);
  m_companions.clear();
  while (singlet && success && trials<m_maxtrials/10) {
    ++trials;
    singlet=false;
    do { comp[0]=SelectCompanion(cur); } while ((comp[1]=FindConnected(comp[0]))==NULL);
    bool anticol=comp[0]->Flav().IsAnti()^comp[0]->Flav().IsDiQuark();
    if (!(cur->Flav().IsAnti()^anticol)) {
      ATOOLS::Particle *rem=comp[0]; comp[0]=comp[1]; comp[1]=rem;
    }
    if (cur->Flav().IsAnti()) {
      int old=comp[1]->GetFlow(2);
      success=success&&AdjustColours(comp[1],old,newpart->GetFlow(1),
				     singlet,m_errors>=m_maxtrials/10);
      if (singlet) continue;
      success=success&&AdjustColours(cur,cur->GetFlow(2),comp[0]->GetFlow(1),
				     singlet,m_errors>=m_maxtrials/10);
      if (singlet) {
	AdjustColours(comp[1],comp[1]->GetFlow(2),old,singlet=false,false);
	singlet=true;
	continue;
      }
    }
    else {
      int old=comp[1]->GetFlow(1);
      success=success&&AdjustColours(comp[1],old,newpart->GetFlow(2),
				     singlet,m_errors>=m_maxtrials/10);
      if (singlet) continue;
      success=success&&AdjustColours(cur,cur->GetFlow(1),comp[0]->GetFlow(2),
				     singlet,m_errors>=m_maxtrials/10);
      if (singlet) {
	AdjustColours(comp[1],comp[1]->GetFlow(1),old,singlet=false,false);
	singlet=true;
	continue;
      }
    }
  }
  if (trials==m_maxtrials/10 && m_errors<m_maxtrials/10) success=false;
  m_parton[0].push_back(newpart);
  p_beamblob->AddToOutParticles(newpart);
  m_parton[2].push_back(cur);
  p_beamblob->AddToOutParticles(cur);
  return success;
}

bool Hadron_Remnant::TreatGluon(ATOOLS::Particle *cur) 
{
  unsigned int trials=0;
  bool success=true, singlet=true;
  ATOOLS::Particle *comp[2];
  m_companions.clear();
  while (singlet && success && trials<m_maxtrials/10) {
    ++trials;
    singlet=false;
    do { comp[0]=SelectCompanion(cur); } while ((comp[1]=FindConnected(comp[0]))==NULL);
    if (comp[0]->Flav().IsAnti()^comp[0]->Flav().IsDiQuark()) {
      int old=comp[0]->GetFlow(2);
      success=success&&AdjustColours(comp[0],old,cur->GetFlow(1),
				     singlet,m_errors>=m_maxtrials/10);
      if (singlet) continue;
      success=success&&AdjustColours(cur,cur->GetFlow(2),comp[1]->GetFlow(1),
				     singlet,m_errors>=m_maxtrials/10);
      if (singlet) {
	AdjustColours(comp[0],comp[0]->GetFlow(2),old,singlet=false,false);
	singlet=true;
	continue;
      }
    }
    else {
      int old=comp[0]->GetFlow(1);
      success=success&&AdjustColours(comp[0],old,cur->GetFlow(2),
				     singlet,m_errors>=m_maxtrials/10);
      if (singlet) continue;
      success=success&&AdjustColours(cur,cur->GetFlow(1),comp[1]->GetFlow(2),
				     singlet,m_errors>=m_maxtrials/10);
      if (singlet) {
	AdjustColours(comp[0],comp[0]->GetFlow(1),old,singlet=false,false);
	singlet=true;
	continue;
      }
    }
  }
  if (trials==m_maxtrials/10 && m_errors<m_maxtrials/10) success=false;
#ifdef DEBUG__Hadron_Remnant
  if (ATOOLS::rpa.gen.NumberOfDicedEvents()==EVENT) {
    std::cout<<"exited with trials = "<<trials<<", singlet = "<<singlet
	     <<", m_errors = "<<m_errors<<", m_maxtrials = "<<m_maxtrials<<std::endl;
  }
#endif
  m_parton[2].push_back(cur);
  p_beamblob->AddToOutParticles(cur);
  return success;
}

double Hadron_Remnant::GetXPDF(ATOOLS::Flavour flavour,double scale) 
{
  double cut, x;
  cut=2.0*(flavour.PSMass()+m_hardpt.PPerp2()/ATOOLS::sqr(m_parton[0].size()))/sqrt(scale);
  if (scale<p_pdfbase->Q2Min()) {
    ATOOLS::msg.Error()<<"Hadron_Remnant::GetXPDF("<<flavour<<","<<scale<<"): "
		       <<"Scale under-runs minimum as given by PDF: "
		       <<scale<<" < "<<p_pdfbase->Q2Min()<<std::endl;
    return cut;
  } 
  if (cut>0.49) return 0.5;
  unsigned int xtrials, pdftrials=0;
  while (true) {
    ++pdftrials;
    xtrials=0;
    do { 
      ++xtrials;
      x=ATOOLS::ran.Get(); 
      if (xtrials>=m_maxtrials) x=cut;
    } while (x<cut);
    p_pdfbase->Calculate(x,0.,0.,scale);
    if (pdftrials>=m_maxtrials) { m_xscheme=0; return 0.01; }
    if (p_pdfbase->GetXPDF(flavour)/x>ATOOLS::ran.Get()) return x;
  } 
  return 0.0;
}

void Hadron_Remnant::UnDo() 
{
  ATOOLS::msg.Tracking()<<"Hadron_Remnant::UnDo(): Undoing changes on blob list."<<std::endl;
#ifdef DEBUG__Hadron_Remnant
  if (ATOOLS::rpa.gen.NumberOfDicedEvents()==EVENT) {
    std::cout<<"Hadron_Remnant::UnDo(): ["<<m_errors<<"] Undoing changes on blob list."<<std::endl;
  }
#endif
  for (int i=(int)m_undo.size()-1;i>=0;--i) {
    bool singlet=false;
    AdjustColours(m_undo[i].first,m_undo[i].first->GetFlow(m_undo[i].second.first),
		  m_undo[i].second.second,singlet,true);
    p_beamblob->RemoveOutParticle(m_undo[i].first);
  }
  Remnant_Base::UnDo();
}
