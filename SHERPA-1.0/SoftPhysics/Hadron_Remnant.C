#include "Hadron_Remnant.H"

#include "Exception.H"
#include "Random.H"

using namespace SHERPA;

Hadron_Remnant::Hadron_Remnant(PDF::ISR_Handler *isrhandler,const unsigned int _m_beam,double _m_scale):
  Remnant_Base(Hadron_Remnant::Hadron,_m_beam),
  m_deltax(0.0125), 
  m_scale(-_m_scale),
  m_xscheme(1), 
  m_maxtrials(100)
{
  if (isrhandler==NULL) {
    throw(ATOOLS::Exception(ATOOLS::ex::fatal_error,"Hadron remnant needs ISR Handler.",
			    "Hadron_Remnant","Hadron_Remnant"));
  }
  p_pdfbase=isrhandler->PDF(m_beam);
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
  // decompose hadron
  m_hardpt=ATOOLS::Vec4D();
  for (size_t i=0;i<m_parton[1].size();++i) {
    ATOOLS::Particle *cur=m_parton[1][i];
    m_hardpt+=cur->Momentum();
    if (i>0) {
      if (cur->Flav().Kfcode()==ATOOLS::kf::gluon) TreatGluon(cur);
      else TreatQuark(cur);
    }
    else {
      if (cur->Flav().Kfcode()==ATOOLS::kf::gluon) TreatFirstGluon(cur);
      else TreatFirstQuark(cur);
    }
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
      beamblob->AddToOutParticles(m_parton[i][j]);
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
    for (unsigned int i=0;i<m_parton[0].size();++i) {
      if (!m_parton[0][i]->Flav().IsDiQuark()) {
	double value=1.0;
	for (unsigned int j=0;(xtot-value<m_deltax)&&(j<m_maxtrials/100);++j) { 
	  value=GetXPDF(m_parton[0][i]->Flav(),m_scale); 
	}
	xmap[m_parton[0][i]]=value;
 	xtot-=value;
      }
      else {
	p_last[0]=m_parton[0][i];
      }
    }
    xmap[p_last[0]]=xtot;
    if (trials>m_maxtrials) {
      ATOOLS::msg.Tracking()<<"Hadron_Remnant::TreatHadron: "
			    <<"Too many trials to find appropriate x values for partons."<<std::endl
			    <<"   Using naive distribution instead."<<std::endl;
      m_xscheme=0;
    }
  } while ((xtot<m_deltax)&&(m_xscheme!=0));
  if (m_xscheme==0) {
    for (std::map<ATOOLS::Particle*,double>::iterator it=xmap.begin();it!=xmap.end(); ++it) {
      it->second=m_xtot/xmap.size();
    }
  }
  for (unsigned int i=0;i<m_parton[1].size();++i) ptot-=m_parton[1][i]->Momentum();
  for (unsigned int j=0;j<m_parton[0].size();++j) {
    double E=xmap[m_parton[0][j]]*m_pbeam[0];
    double pz=ATOOLS::Sign(m_pbeam[3])*sqrt(E*E-ATOOLS::sqr(m_parton[0][j]->Flav().PSMass())
					    -m_hardpt.PPerp2()/ATOOLS::sqr(m_parton[0].size()));
    m_parton[0][j]->SetMomentum(ATOOLS::Vec4D(E,-m_hardpt[1]/m_parton[0].size(),
					      -m_hardpt[2]/m_parton[0].size(),pz));
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
  for (unsigned int i=0;i<2;++i) m_parton[0].push_back(newpart[i]);
  m_parton[2].push_back(cur);
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
    for (unsigned int i=0;i<3;++i) m_parton[0].push_back(newpart[i]);
    p_last[0]=newpart[0];
  }
  m_parton[2].push_back(cur);
  return true;
}

bool Hadron_Remnant::AdjustColours(ATOOLS::Particle *particle,unsigned int oldc,unsigned int newc,
				   bool anti,bool forward,unsigned int catcher)
{
  if (oldc==newc) return true;
  if (++catcher>100) {
    ATOOLS::msg.Tracking()<<"Hadron_Remnant::AdjustColours(..): "
			  <<"Colour nesting is too deep!"<<std::endl
			  <<"   Cannot adjust colours completely. Result might be unreliable."<<std::endl;
    return false;
  }
  ATOOLS::Blob *cur;
  if (forward) cur=particle->DecayBlob();
  else cur=particle->ProductionBlob();
  if (cur!=NULL) {
    for (int i=0;i<cur->NOutP();++i) {
      for (int j=1;j<3;++j) {
	ATOOLS::Particle *help=cur->OutParticle(i);
	if (help->GetFlow(j)==(int)oldc && 
	    (help->Flav().IsAnti()!=anti || help->Flav().IsGluon())) {
	  help->SetFlow(j,newc);
	  return AdjustColours(help,oldc,newc,anti,true,catcher);
	}
      }
    }
    for (int i=0;i<cur->NInP();++i) {
      for (int j=1;j<3;++j) {
	ATOOLS::Particle *help=cur->InParticle(i);
	if (help->GetFlow(j)==(int)oldc && 
	    (help->Flav().IsAnti()!=anti || help->Flav().IsGluon())) {
	  help->SetFlow(j,newc);
	  return AdjustColours(help,oldc,newc,anti,true,catcher);
	}
      }
    }
  }
  else {
    return true;
  }
  return false;
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
  size_t set=2*(size_t)ATOOLS::ran.Get();
  return m_parton[set][(size_t)(ATOOLS::ran.Get()*((int)m_parton[set].size()-1))];
}

bool Hadron_Remnant::TreatQuark(ATOOLS::Particle *cur) 
{
  bool success=true;
  ATOOLS::Particle *comp[2];
  do { comp[0]=SelectCompanion(cur); } while ((comp[1]=FindConnected(comp[0]))==NULL);
  ATOOLS::Particle *newpart = new ATOOLS::Particle(-1,cur->Flav().Bar());
  newpart->SetFlow(1+(int)newpart->Flav().IsAnti(),-1);
  bool anticol=comp[0]->Flav().IsAnti()^comp[0]->Flav().IsDiQuark();
  if (!(cur->Flav().IsAnti()^anticol)) {
      ATOOLS::Particle *rem=comp[0]; comp[0]=comp[1]; comp[1]=rem;
  }
  if (cur->Flav().IsAnti()) {
    int catcher=0;
    success=success&&AdjustColours(comp[1],comp[1]->GetFlow(2),newpart->GetFlow(1),
				   false,true,catcher);
    success=success&&AdjustColours(cur,cur->GetFlow(2),comp[0]->GetFlow(1),
				   false,true,catcher);
    comp[1]->SetFlow(2,newpart->GetFlow(1));
    cur->SetFlow(2,comp[0]->GetFlow(1));
  }
  else {
    int catcher=0;
    success=success&&AdjustColours(comp[1],comp[1]->GetFlow(1),newpart->GetFlow(2),
				   true,true,catcher);
    success=success&&AdjustColours(cur,cur->GetFlow(1),comp[0]->GetFlow(2),
				   true,true,catcher);
    comp[1]->SetFlow(1,newpart->GetFlow(2));
    cur->SetFlow(1,comp[0]->GetFlow(2));
  }
  m_parton[0].push_back(newpart);
  m_parton[2].push_back(cur);
  return success;
}

bool Hadron_Remnant::TreatGluon(ATOOLS::Particle *cur) 
{
  bool success=true;
  ATOOLS::Particle *comp[2];
  do { comp[0]=SelectCompanion(cur); } while ((comp[1]=FindConnected(comp[0]))==NULL);
  if (comp[0]->Flav().IsAnti()^comp[0]->Flav().IsDiQuark()) {
    int catcher=0;
    success=success&&AdjustColours(comp[0],comp[0]->GetFlow(2),cur->GetFlow(1),
				   false,true,catcher);
    success=success&&AdjustColours(comp[1],comp[1]->GetFlow(1),cur->GetFlow(2),
				   true,true,catcher);
    comp[0]->SetFlow(2,cur->GetFlow(1));
    comp[1]->SetFlow(1,cur->GetFlow(2));
  }
  else {
    int catcher=0;
    success=success&&AdjustColours(comp[0],comp[0]->GetFlow(1),cur->GetFlow(2),
				   true,true,catcher);
    success=success&&AdjustColours(comp[1],comp[1]->GetFlow(2),cur->GetFlow(1),
				   false,true,catcher);
    comp[0]->SetFlow(1,cur->GetFlow(2));
    comp[1]->SetFlow(2,cur->GetFlow(1));
  }
  m_parton[2].push_back(cur);
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
    p_pdfbase->GetBasicPDF()->Calculate(x,0.,0.,scale);
    if (pdftrials>=m_maxtrials) { m_xscheme=0; return 0.01; }
    if (p_pdfbase->GetBasicPDF()->GetXPDF(flavour)/x>ATOOLS::ran.Get()) return x;
  } 
  return 0.0;
}


