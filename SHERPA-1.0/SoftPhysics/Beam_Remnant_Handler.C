#include "Beam_Remnant_Handler.H"

#include "Primordial_KPerp.H"
#include "Running_AlphaS.H"
#include "Run_Parameter.H"
#include "Data_Read.H"
#include "Random.H"

using namespace SHERPA;
using namespace ATOOLS;

std::set<Blob*> checked_blobs;

void SumMomenta(Blob * bl, Vec4D & inisum, Vec4D & finsum,bool iterate=true) 
{
  // check if caught in loop
  std::set<Blob*>::const_iterator bit=checked_blobs.find(bl);
  if (bit!=checked_blobs.end()) return;
  checked_blobs.insert(bl);
  for (int i=0;i<bl->NInP();++i) {
    Particle * p =bl->InParticle(i);
    if (p->ProductionBlob()==NULL) inisum+=p->Momentum();
    else if (iterate) SumMomenta(p->ProductionBlob(),inisum,finsum);
    else inisum+=p->Momentum();
  }
  for (int i=0;i<bl->NOutP();++i) {
    Particle * p =bl->OutParticle(i);
    if (p->DecayBlob()==NULL) finsum+=p->Momentum();
    else if (iterate) SumMomenta(p->DecayBlob(),inisum,finsum);
    else finsum+=p->Momentum();
  }
}

bool SumMomenta(Blob * bl)
{
  Vec4D inisum,finsum;
  checked_blobs.clear();
  SumMomenta(bl,inisum,finsum);
  std::cout<<" inisum="<<inisum<<std::endl;
  std::cout<<" finsum="<<finsum<<std::endl;
  return (inisum==finsum);
}

Beam_Remnant_Handler::Beam_Remnant_Handler(std::string _m_path,std::string _m_file,
					   PDF::ISR_Handler * _p_isr,BEAM::Beam_Spectra_Handler * _p_beam,
					   double _m_scale) :
  p_isr(_p_isr), 
  p_beam(_p_beam),
  p_constituents(NULL), 
  p_numberofconstituents(NULL), 
  m_path(_m_path), 
  m_file(_m_file), 
  m_deltax(0.0125), 
  m_scale(-_m_scale),
  m_xscheme(1), 
  m_maxtrials(100)
{
  p_numberofconstituents = new int[2];
  for (int i=0;i<2;i++) p_numberofconstituents[i] = 0;
  if (p_isr->Flav(0).IsHadron() || p_isr->Flav(1).IsHadron()) {
    p_constituents = new Flavour *[2];
    for (int i=0;i<2;i++) {
      if (p_isr->Flav(i).IsHadron()) 
	p_numberofconstituents[i] = Constituents(p_isr->Flav(i),p_constituents[i]);
      else p_constituents[i] = NULL;
    }
    m_fill = 1;
  }
  else if ((p_isr->Flav(0).Kfcode()==kf::e && (p_isr->On()==1 || p_isr->On()==3)) ||
	   (p_isr->Flav(1).Kfcode()==kf::e && (p_isr->On()==2 || p_isr->On()==3))) {
    m_fill = 1;
  }
  else m_fill = 0;
  p_particle[0] = new ATOOLS::Particle_List[3];
  p_particle[1] = new ATOOLS::Particle_List[3];
  p_x = new std::map<ATOOLS::Particle*,double>();

  m_q2min = 1.;

  p_kperp = new Primordial_KPerp(_m_path,_m_file);
}

int Beam_Remnant_Handler::Constituents(Flavour _had,Flavour *& _flavs) {
  int hadint = (_had.Kfcode() - _had.Kfcode()/10000)/10;
  
  if ((hadint > 100) && (hadint < 1000)) {
    _flavs    = new Flavour[3];
    _flavs[0] = Flavour(kf::code(hadint)/100);
    _flavs[1] = Flavour(kf::code((hadint-(hadint/100)*100)/10));
    _flavs[2] = Flavour(kf::code(hadint-(hadint/10)*10));
    if (_had.IsAnti()) {
      for(int i=0;i<3;i++) _flavs[i]=_flavs[i].Bar();
    }
    msg.Tracking()<<"Beam_Remnant_Handler::Constituents for "<<_had<<" is baryon :"
		  <<hadint<<std::endl<<"   "<<_flavs[0]<<", "<<_flavs[1]<<", "<<_flavs[2]<<std::endl;
    return 3;
  }
  if ((hadint > 10) && (hadint < 100)) {
    _flavs    = new Flavour[2];
    _flavs[0] = Flavour(kf::code(hadint)/10);
    _flavs[1] = Flavour(kf::code(hadint-(hadint/10)*10));
    if (_had.IsAnti()) {
      for(int i=0;i<2;i++) _flavs[i]=_flavs[i].Bar();
    }
    msg.Tracking()<<"Beam_Remnant_Handler::Constituents for "<<_had<<" is meson :"
		  <<hadint<<std::endl<<"   "<<_flavs[0]<<", "<<_flavs[1]<<std::endl;
    return 2;
  }
  
  msg.Error()<<"Beam_Remnant_Handler::Constituents :"
	     <<"No idea how to handle this case: "<<_had<<". Abort."<<std::endl;
  abort();
}      

Beam_Remnant_Handler::~Beam_Remnant_Handler() {  
  if (p_constituents) {
    for (int i=0;i<2;i++) {
      if (p_constituents[i]) delete [] p_constituents[i];
    }
    delete [] p_constituents; p_constituents = NULL;
  }
  
  if (p_numberofconstituents) delete [] p_numberofconstituents;
  delete [] p_particle[0];
  delete [] p_particle[1];
  delete p_kperp;
  delete p_x;
}


bool Beam_Remnant_Handler::FillBunchBlobs(Blob_List * _bloblist,Particle_List * _particlelist)
{
  Blob_Iterator endblob = _bloblist->end(); 
  Blob * blob;
  bool flag=false;
  int  pos,pos1,pos2,number;
  Particle * p;
  for (int i=0;i<2;i++) {
    for (Blob_Iterator biter = _bloblist->begin();biter != endblob;++biter) {
      pos1 = (*biter)->Type().find(std::string("Beam Remnant"));
      pos2 = (*biter)->Type().find(std::string("IS"));
      pos  = Max(pos1,pos2);
      if ((*biter)->Status()==1 && (*biter)->Beam()==i && pos>-1) {
	(*biter)->SetStatus(2);
	blob = new ATOOLS::Blob();
	blob->SetId(_bloblist->size());
	blob->SetType(std::string("Bunch"));
	blob->SetBeam(i);
	blob->AddToOutParticles((*biter)->InParticle(0));
	(*biter)->InParticle(0)->SetProductionBlob(blob);
	if ((*biter)->InParticle(0)->Flav()==p_beam->GetBeam(i)->Beam() &&
	    IsEqual((*biter)->InParticle(0)->E(),p_beam->GetBeam(i)->InMomentum()[0])) {
	  p = new Particle((*biter)->InParticle(0));
	  if (_particlelist) number = _particlelist->size();
	  else number = (long int)p;
	  p->SetNumber(number);
	  p->SetDecayBlob(blob);
	  p->SetProductionBlob(NULL);
	  blob->AddToInParticles(p);
	}
	else {
	  p = new Particle(-1,p_beam->GetBeam(i)->Beam(),p_beam->GetBeam(i)->InMomentum());
	  if (_particlelist) number = _particlelist->size();
	  else number = (long int)p;
	  p->SetNumber(number);	  
	  p->SetStatus(2);
	  p->SetDecayBlob(blob);
	  blob->AddToInParticles(p);
	  Particle * p = new Particle(-1,p_beam->GetBeam(i)->Remnant(),
				      p_beam->GetBeam(i)->InMomentum()+(-1.)*(*biter)->InParticle(0)->Momentum());
	  if (_particlelist) number = _particlelist->size();
	  else number = (long int)p;
	  p->SetNumber(number);
	  blob->AddToOutParticles(p);
	  p->SetProductionBlob(blob);
	}
	_bloblist->insert(_bloblist->begin(),blob);
	flag=true;
      }
    }
  }
//   for (Blob_Iterator biter = _bloblist->begin();biter != endblob;++biter) {
//     cout<<(*biter);
//   }
//  SumMomenta(_bloblist->front());
  return flag;
}

bool Beam_Remnant_Handler::FillBeamBlobs(Blob_List * _bloblist,Particle_List * _particlelist)
{ 
  if (!m_fill) return false;
  Blob_Iterator endblob = _bloblist->end(); 
  Blob * blob;
  bool okay=false, hadron[2];
  hadron[1]=hadron[0]=false;
  int pos,number;
  for (int i=0;i<2;i++) {
    bool flag=true,treat=false;
    for (Blob_Iterator biter = _bloblist->begin();biter != endblob;++biter) {
      pos = (*biter)->Type().find(std::string("IS"));
      if (m_fill) {
	if ((*biter)->Beam()==i && (*biter)->Status()==1 && pos>-1) { 
	  if (p_numberofconstituents[i]>0) {
	    if (flag) {
	      blob = new ATOOLS::Blob();
	      blob->SetId(_bloblist->size());
	      blob->SetType(std::string("Beam Remnant"));
	      blob->SetBeam(i);
	      blob->SetStatus(1);
	      p_last[i] = new Particle(-1,p_isr->Flav(i),p_beam->GetBeam(i)->OutMomentum());
	      blob->AddToInParticles(p_last[i]);
	      p_last[i]->SetStatus(2);
	      p_last[i]->SetDecayBlob(blob);
	      _bloblist->insert(_bloblist->begin(),blob);
	      p_blob[i]=blob;
	      flag=false;
	      okay=true;
	    }
	    blob->AddToOutParticles((*biter)->InParticle(0));
	    (*biter)->InParticle(0)->SetProductionBlob(blob);
	    (*biter)->SetStatus(2);
	    treat=true;
	  }
	  else if (!IsEqual((*biter)->InParticle(0)->E(),p_beam->GetBeam(i)->OutMomentum()[0])) {
	    (*biter)->SetStatus(2);
	    blob = new ATOOLS::Blob();
	    blob->SetId(_bloblist->size());
	    blob->SetType(std::string("Beam Remnant"));
	    blob->SetBeam(i);
	    blob->SetStatus(1);
	    blob->AddToOutParticles((*biter)->InParticle(0));
	    (*biter)->InParticle(0)->SetProductionBlob(blob);
	    p_last[i] = new Particle(-1,p_isr->Flav(i),p_beam->GetBeam(i)->OutMomentum());
	    if (_particlelist) number = _particlelist->size();
	    else number = (long int)p_last[i];
	    p_last[i]->SetNumber(number);
	    p_last[i]->SetStatus(2);
	    p_last[i]->SetDecayBlob(blob);
	    blob->AddToInParticles(p_last[i]);
	    Particle *p = new Particle(-1,Flavour(kf::photon),
				       p_beam->GetBeam(i)->OutMomentum()+(-1.)*(*biter)->InParticle(0)->Momentum());
	    if (_particlelist) number = _particlelist->size();
	    else number = (long int)p;
	    p->SetNumber(number);
	    blob->AddToOutParticles(p);
	    p->SetProductionBlob(blob);
	    _bloblist->insert(_bloblist->begin(),blob);
	    flag=false;
	    okay=true;
	  }
	}
      }
    }
    if ((p_isr->Flav(i).IsHadron())&&(treat)) {
      okay==okay&&FillHadron(blob);
      hadron[i]=true;
    }
  }
  for (size_t beam=0;beam<2;++beam) if (hadron[beam]) {
    FillHadronBlob(_particlelist,beam);
  }
  bool success;
  size_t trials=0;
  if (hadron[0] && hadron[1]) 
  do {
    ++trials;
    success=true;
    for (m_beam=0;m_beam<2;++m_beam) if (hadron[m_beam]) {
      if (!FillHadron()) success=false;
    }
    m_Erem=rpa.gen.Ecms(); m_pzrem=0.0;
    for (size_t j=0;j<2;++j) p_kperp->CreateKPerp(p_blob[j],p_blob[j]->NOutP());
    for (size_t j=0;j<2;++j) {
      for (size_t i=0;i<p_particle[j][2].size();++i) {
	if (!p_kperp->FillKPerp(p_particle[j][2][i],j)) success=false;
	m_Erem-=p_particle[j][2][i]->Momentum()[0]; 
	m_pzrem-=p_particle[j][2][i]->Momentum()[3];
      }
      for (size_t i=0;i<p_particle[j][1].size();++i) {
	if (!p_kperp->FillKPerp(p_particle[j][1][i],j)) success=false;
	if ((p_particle[j][1][i]!=p_last[0])&&
	    (p_particle[j][1][i]!=p_last[1])) {
	  m_Erem-=p_particle[j][1][i]->Momentum()[0]; 
	  m_pzrem-=p_particle[j][1][i]->Momentum()[3];
	}
      }
    }
    if (hadron[0]||hadron[1]) FillLastPartons();
    if (trials==m_maxtrials) {
      for (size_t i=0;i<2;++i) {
	p_kperp->SetKPerpMean(p_kperp->KPerpMean(i)/10.0,i);
	if (IsZero(p_kperp->KPerpMean(i))) p_kperp->SetKPerpMean(0.0,i);
      }
      trials=0;
    }
  } while (!success);
  //  SumMomenta(_bloblist->front());
  return okay;
}

void Beam_Remnant_Handler::FillHadronBlob(ATOOLS::Particle_List *pl,unsigned int beam)
{
  ATOOLS::Blob *blob=p_blob[beam];
  for (size_t j=0;j<p_particle[beam][1].size();++j) {
    p_particle[beam][1][j]->SetNumber((long int)p_particle[beam][1][j]);
    p_particle[beam][1][j]->SetInfo('F');
    blob->AddToOutParticles(p_particle[beam][1][j]);
    p_particle[beam][1][j]->SetProductionBlob(blob);
    if (pl!=NULL) {
      p_particle[beam][1][j]->SetNumber(pl->size());
      pl->push_back(p_particle[beam][1][j]);
    }
  }
}

bool Beam_Remnant_Handler::FillLastPartons()
{
  // DIS scenario still missing!
  Vec4D pr1=p_last[0]->Momentum(), pr2=p_last[1]->Momentum();
  double Erem=m_Erem;
  double pzrem=m_pzrem;
  double sp, sp1, sp2, lam2, c1, c2, c3, E1, E2, pz1, pz2;
  sp1=sqr(p_last[0]->Flav().PSMass())+sqr(pr1[1])+sqr(pr1[2]);
  sp2=sqr(p_last[1]->Flav().PSMass())+sqr(pr2[1])+sqr(pr2[2]);
  sp=Erem*Erem-pzrem*pzrem;
  lam2=(sp-sp1-sp2)*(sp-sp1-sp2)/4.0-sp1*sp2;
  c1=0.5*(sp-sp1+sp2); c2=0.5*(sp+sp1-sp2); 
  c3=0.5*(sp-sp1-sp2)*Erem*Erem-lam2;
  double spn, ytn, yto=(Erem+pzrem)/(Erem-pzrem);
  for (double sign=1.0;sign>=-1.0;sign-=2.0) {
    E1=Erem*c2/(c1+c2)*(1.0+sign*sqrt(1.0+(c3/(Erem*Erem)-c2)/(c2*c2)*(c1+c2)));
    E2=Erem-E1;
    pz1=Sign(pr1[3])*sqrt(E1*E1-sp1); pz2=Sign(pr2[3])*sqrt(E2*E2-sp2);
    spn=sqr(Erem)-sqr(pz1+pz2);
    if (!IsEqual(spn,sp)) pz1*=-1.0; { spn=sqr(E1+E2)-sqr(pz1+pz2); }
    if (!IsEqual(spn,sp)) continue;
    ytn=(E1+E2+pz1+pz2)/(E1+E2-pz1-pz2);
    if (!IsEqual(ytn,yto)) { pz1*=-1.0; pz2*=-1.0; ytn=(E1+E2+pz1+pz2)/(E1+E2-pz1-pz2); }
    if (!IsEqual(ytn,yto)) continue;
  }
  pr1=Vec4D(E1,pr1[1],pr1[2],pz1); p_last[0]->SetMomentum(pr1);
  pr2=Vec4D(E2,pr2[1],pr2[2],pz2); p_last[1]->SetMomentum(pr2);
  return true;
}

bool Beam_Remnant_Handler::FillHadron(ATOOLS::Blob *blob)
{
  m_beam=blob->Beam();
  for (unsigned int i=0;i<3;++i) p_particle[m_beam][i].clear();
  for (int i=0;i<blob->NOutP();++i) {
    Particle *cur=blob->OutParticle(i);
    if (i>0) {
      if (cur->Flav().Kfcode()==kf::gluon) TreatGluon(cur,blob->Beam());
      else TreatQuark(cur,blob->Beam());
    }
    else {
      if (cur->Flav().Kfcode()==kf::gluon) TreatFirstGluon(cur,blob->Beam());
      else TreatFirstQuark(cur,blob->Beam());
    }
  }
//   p_kperp->CreateKPerp(blob,p_particle[m_beam][1].size()+blob->NOutP());
  return true;
}

bool Beam_Remnant_Handler::FillHadron()
{
  unsigned int trials;
  Vec4D ptot=p_beam->GetBeam(m_beam)->OutMomentum();
  m_xtot=1.0;
  for (unsigned int i=0;i<p_particle[m_beam][2].size();++i) {
    m_xtot-=p_particle[m_beam][2][i]->Momentum()[0]/ptot[0];
  }
  trials=0;
  m_xscheme=1;
  double xtot;
  p_x->clear();
  do {
    ++trials;
    xtot=m_xtot;
    for (unsigned int i=0;i<p_particle[m_beam][1].size();++i) {
      if (!p_particle[m_beam][1][i]->Flav().IsDiQuark()) {
	double value=1.0;
	for (unsigned int j=0;(xtot-value<m_deltax)&&(j<m_maxtrials/100);++j) { 
	  value=GetXPDF(p_particle[m_beam][1][i]->Flav(),m_scale,m_beam); 
	}
	(*p_x)[p_particle[m_beam][1][i]]=value;
 	xtot-=value;
      }
      else {
	p_last[m_beam]=p_particle[m_beam][1][i];
      }
    }
    (*p_x)[p_last[m_beam]]=xtot;
    if (trials>m_maxtrials) {
      ATOOLS::msg.Tracking()<<"Beam_Remnant_Handler::TreatHadron: "
			    <<"Too many trials to find appropriate x values for partons."<<std::endl
			    <<"   Using naive distribution instead."<<std::endl;
      m_xscheme=0;
    }
  } while ((xtot<m_deltax)&&(m_xscheme!=0));
  if (m_xscheme==0) {
    for (std::map<ATOOLS::Particle*,double>::iterator it=p_x->begin();it!=p_x->end(); ++it) {
      it->second=m_xtot/p_x->size();
    }
  }
  for (unsigned int i=0;i<p_particle[m_beam][2].size();++i) {
    ptot-=p_particle[m_beam][2][i]->Momentum();
  }
  double pztot=ptot[3];
  for (unsigned int j=0;j<p_particle[m_beam][1].size();++j) {
    double E=(*p_x)[p_particle[m_beam][1][j]]/m_xtot*ptot[0];
    double pz=Sign(pztot)*sqrt(E*E-sqr(p_particle[m_beam][1][j]->Flav().PSMass()));
    if (p_particle[m_beam][1][j]!=p_last[m_beam]) pztot-=pz;
    p_particle[m_beam][1][j]->SetMomentum(ATOOLS::Vec4D(E,0.0,0.0,pz));
    if (p_last[m_beam]==p_particle[m_beam][1][j]) continue;
  }
  double E=(*p_x)[p_last[m_beam]]/m_xtot*ptot[0];
  p_last[m_beam]->SetMomentum(ATOOLS::Vec4D(E,0.0,0.0,pztot));
  return true;
}

bool Beam_Remnant_Handler::TreatFirstGluon(ATOOLS::Particle *cur,unsigned int beam) 
{
  unsigned int single=(unsigned int)(ran.Get()*3.); 
  ATOOLS::Flavour difl, fl=p_constituents[beam][single];
  int di[2];
  for (unsigned int j=0, i=0;i<3;i++) {
    if (i!=single) di[j++] = p_constituents[beam][i].Kfcode();
  }
  if (di[0]!=di[1]) {
    if (ran.Get()<0.25) difl=Flavour(kf::code(abs(di[0])*1000+abs(di[1])*100+1));
    else difl=Flavour(kf::code(abs(di[0])*1000+abs(di[1])*100+3));
    if (di[0]<0) difl=difl.Bar();
  }
  else {
    difl=Flavour(kf::code(abs(di[0])*1100+3));
    if (di[0]<0) difl=difl.Bar();
  }
  Particle *newpart[2];
  newpart[1] = new Particle(-1,fl); 
  newpart[0] = new Particle(-1,difl); 
  newpart[1]->SetInfo('F');
  newpart[0]->SetInfo('F');
  if (newpart[0]->Flav().IsAnti()) {
    newpart[1]->SetFlow(1,0);
    newpart[1]->SetFlow(2,cur->GetFlow(1));
    newpart[0]->SetFlow(1,cur->GetFlow(2));
    newpart[0]->SetFlow(2,0);
  }
  else {
    newpart[1]->SetFlow(1,cur->GetFlow(2));
    newpart[1]->SetFlow(2,0);
    newpart[0]->SetFlow(1,0);
    newpart[0]->SetFlow(2,cur->GetFlow(1));
  }
  for (unsigned int i=0;i<2;++i) p_particle[m_beam][1].push_back(newpart[i]);
  p_particle[m_beam][2].push_back(cur);
  return true;
}

bool Beam_Remnant_Handler::TreatFirstQuark(ATOOLS::Particle *cur,unsigned int beam) 
{
  bool found=false;
  int di[3];
  for (int j=0,i=0;i<p_numberofconstituents[beam];++i) {
    if (found||cur->Flav()!=p_constituents[beam][i]) di[j++]=p_constituents[beam][i];
    else found=true;
  }
  if (found) {
    ATOOLS::Flavour difl;
    if (di[0]!=di[1]) {
      if (ran.Get()<0.25) difl=Flavour(kf::code(abs(di[0])*1000+abs(di[1])*100+1));
      else difl=ATOOLS::Flavour(kf::code(abs(di[0])*1000+abs(di[1])*100+3));
      if (di[0]<0) difl=difl.Bar();
    }
    else {
      difl=Flavour(kf::code(abs(di[0])*1100+3));
      if (di[0]<0) difl=difl.Bar();
    }
    Particle *newpart = new Particle(-1,difl); 
    newpart->SetFlow(1,cur->GetFlow(2));
    newpart->SetFlow(2,cur->GetFlow(1));
    newpart->SetInfo('F');
    p_particle[m_beam][1].push_back(newpart);
    p_particle[m_beam][2].push_back(cur);
  }
  else {
    unsigned int single=(unsigned int)(ran.Get()*3.0); 
    ATOOLS::Flavour difl, fl=p_constituents[beam][single];
    for (unsigned int j=0, i=0;i<3;i++) {
      if (i!=single) di[j++] = p_constituents[beam][i].Kfcode();
    }
    if (di[0]!=di[1]) {
      if (ran.Get()<0.25) difl = Flavour(kf::code(di[0]*1000+di[1]*100+1));
      else difl = Flavour(kf::code(di[0]*1000+di[1]*100+3));
    }
    else difl = Flavour(kf::code(di[0]*1100+3));
    if (p_constituents[beam][0].IsAnti()) difl = difl.Bar();
    ATOOLS::Particle *newpart[3];
    newpart[2] = new Particle(-1,fl); 
    newpart[1] = new Particle(-1,(cur->Flav()).Bar()); 
    newpart[0] = new Particle(-1,difl); 
    newpart[2]->SetInfo('F');
    newpart[1]->SetInfo('F');
    newpart[0]->SetInfo('F');
    if ( (fl.IsAnti() && !((cur->Flav()).IsAnti()) ) ||
	 (!(fl.IsAnti()) && (cur->Flav()).IsAnti() ) ) {
      if (newpart[0]->Flav().IsAnti()) {
	newpart[2]->SetFlow(2,cur->GetFlow(1));
	newpart[1]->SetFlow(1,0);
	newpart[1]->SetFlow(2,-1);
	newpart[0]->SetFlow(1,newpart[1]->GetFlow(2));
	newpart[0]->SetFlow(2,0);
      }
      else {
	newpart[2]->SetFlow(1,cur->GetFlow(2));
	newpart[1]->SetFlow(1,-1);
	newpart[1]->SetFlow(2,0);
	newpart[0]->SetFlow(1,0);
	newpart[0]->SetFlow(2,newpart[1]->GetFlow(1));
      }
    }
    else {
      newpart[0]->SetFlow(1,cur->GetFlow(2));
      newpart[0]->SetFlow(2,cur->GetFlow(1));
      if (fl.IsAnti()) {
	newpart[1]->SetFlow(1,-1);
	newpart[1]->SetFlow(2,0);
	newpart[2]->SetFlow(1,0);
	newpart[2]->SetFlow(2,newpart[1]->GetFlow(1));
      }
      else {
	newpart[1]->SetFlow(1,0);
	newpart[1]->SetFlow(2,-1);
	newpart[2]->SetFlow(1,newpart[1]->GetFlow(2));
	newpart[2]->SetFlow(2,0);
      }
    }
    for (unsigned int i=0;i<3;++i) p_particle[m_beam][1].push_back(newpart[i]);
    p_particle[m_beam][2].push_back(cur);
  }
  return true;
}

bool Beam_Remnant_Handler::AdjustColours(ATOOLS::Particle *particle,unsigned int oldc,unsigned int newc,
					 unsigned int catcher,bool forward)
{
  if (oldc==newc) return true;
  if (++catcher>100) {
    ATOOLS::msg.Tracking()<<"Beam_Remnant_Handler::AdjustColours(..): "
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
	if (cur->OutParticle(i)->GetFlow(j)==(int)oldc) {
	  cur->OutParticle(i)->SetFlow(j,newc);
	  return AdjustColours(cur->OutParticle(i),oldc,newc,catcher,true);
	}
      }
    }
    for (int i=0;i<cur->NInP();++i) {
      for (int j=1;j<3;++j) {
	if (cur->InParticle(i)->GetFlow(j)==(int)oldc) {
	  cur->InParticle(i)->SetFlow(j,newc);
	  return AdjustColours(cur->InParticle(i),oldc,newc,catcher,false);
	}
      }
    }
  }
  else {
    return true;
  }
  return false;
}

ATOOLS::Particle *Beam_Remnant_Handler::FindConnected(ATOOLS::Particle *particle,bool same,int orig) 
{
  if (particle->Flav().IsQuark()) {
    if (particle->Flav().IsAnti()) orig=2; else orig=1; 
  }
  else if (particle->Flav().IsDiQuark()) orig=2;
  int comp=3-orig;
  if (same) comp=orig;
  for (unsigned int set=1;set<3;++set) {
    for (ATOOLS::Particle_Iterator pit=p_particle[m_beam][set].begin();pit!=p_particle[m_beam][set].end();++pit) {
      if ((*pit)->GetFlow(comp)==particle->GetFlow(orig)) return *pit;
    }
  }
  return NULL;
}

ATOOLS::Particle *Beam_Remnant_Handler::SelectCompanion(ATOOLS::Particle *cur) 
{
  ATOOLS::Particle *selected;
  if (true) {
    unsigned int set=(unsigned int)ATOOLS::ran.Get();
    selected=p_particle[m_beam][1+set][(unsigned int)(ATOOLS::ran.Get()*((int)p_particle[m_beam][1+set].size()-1))];
  }
  return selected;
}

bool Beam_Remnant_Handler::TreatQuark(ATOOLS::Particle *cur,unsigned int beam) 
{
  bool success=true;
  ATOOLS::Particle *comp[2];
  do { comp[0]=SelectCompanion(cur); } while ((comp[1]=FindConnected(comp[0]))==NULL);
  ATOOLS::Particle *newpart = new ATOOLS::Particle(-1,cur->Flav().Bar());
  newpart->SetFlow(1+(int)newpart->Flav().IsAnti(),-1);
  bool quark[2], antiquark[2], diquark[2];
  quark[0]=comp[0]->Flav().IsQuark()&&!comp[0]->Flav().IsAnti();
  quark[1]=comp[1]->Flav().IsQuark()&&!comp[1]->Flav().IsAnti();
  antiquark[0]=comp[0]->Flav().IsQuark()&&comp[0]->Flav().IsAnti();
  antiquark[1]=comp[1]->Flav().IsQuark()&&comp[1]->Flav().IsAnti();
  diquark[0]=comp[0]->Flav().IsDiQuark(); diquark[1]=comp[1]->Flav().IsDiQuark();
  if (!cur->Flav().IsAnti()) {
    if (quark[0]||antiquark[1]||diquark[1]) {
      ATOOLS::Particle *rem=comp[0]; comp[0]=comp[1]; comp[1]=rem;
    }
  }
  else {
    if (antiquark[0]||diquark[0]||quark[1]) {
      ATOOLS::Particle *rem=comp[0]; comp[0]=comp[1]; comp[1]=rem;
    }
  }
  if (cur->Flav().IsAnti()) {
    int catcher=0;
    success=success&&AdjustColours(comp[1],comp[1]->GetFlow(2),newpart->GetFlow(1),catcher);
    success=success&&AdjustColours(cur,cur->GetFlow(2),comp[0]->GetFlow(1),catcher);
    comp[1]->SetFlow(2,newpart->GetFlow(1));
    cur->SetFlow(2,comp[0]->GetFlow(1));
  }
  else {
    int catcher=0;
    success=success&&AdjustColours(comp[1],comp[1]->GetFlow(1),newpart->GetFlow(2),catcher);
    success=success&&AdjustColours(cur,cur->GetFlow(1),comp[0]->GetFlow(2),catcher);
    comp[1]->SetFlow(1,newpart->GetFlow(2));
    cur->SetFlow(1,comp[0]->GetFlow(2));
  }
  p_particle[m_beam][1].push_back(newpart);
  p_particle[m_beam][2].push_back(cur);
  return success;
}

bool Beam_Remnant_Handler::TreatGluon(ATOOLS::Particle *cur,unsigned int beam) 
{
  bool success=true;
  ATOOLS::Particle *comp[2];
  do { comp[0]=SelectCompanion(cur); } while ((comp[1]=FindConnected(comp[0]))==NULL);
  if (comp[0]->Flav().IsAnti()||comp[0]->Flav().IsDiQuark()) {
    int catcher=0;
    success=success&&AdjustColours(comp[0],comp[0]->GetFlow(2),cur->GetFlow(1),catcher);
    success=success&&AdjustColours(comp[1],comp[1]->GetFlow(1),cur->GetFlow(2),catcher);
    comp[0]->SetFlow(2,cur->GetFlow(1));
    comp[1]->SetFlow(1,cur->GetFlow(2));
  }
  else {
    int catcher=0;
    success=success&&AdjustColours(comp[0],comp[0]->GetFlow(1),cur->GetFlow(2),catcher);
    success=success&&AdjustColours(comp[1],comp[1]->GetFlow(2),cur->GetFlow(1),catcher);
    comp[0]->SetFlow(1,cur->GetFlow(2));
    comp[1]->SetFlow(2,cur->GetFlow(1));
  }
  p_particle[m_beam][2].push_back(cur);
  return success;
}

double Beam_Remnant_Handler::GetXPDF(Flavour flav,double scale,int beam) 
{
  double cut, x;
  PDF::PDF_Base *pdf=p_isr->PDF(beam);
  cut=2.0*flav.PSMass()/sqrt(scale);
  if (scale<pdf->GetQ2Min()) {
    ATOOLS::msg.Error()<<"Beam_Remnant_Handler::GetXPDF("<<flav<<","<<scale<<","<<beam<<"): "
		       <<"Scale under-runs minimum as given by PDF: "
		       <<scale<<" < "<<pdf->GetQ2Min()<<std::endl;
    return cut;
  } 
  if (cut>0.49) return 0.5;
  unsigned int xtrials, pdftrials=0;
  while (true) {
    ++pdftrials;
    xtrials=0;
    do { 
      ++xtrials;
      x=ran.Get(); 
      if (xtrials>=m_maxtrials) x=cut;
    } while (x<cut);
    pdf->Calculate(x,scale);
    if (pdftrials>=m_maxtrials) { m_xscheme=0; return 0.01; }
    if (pdf->GetXPDF(flav)/x>ran.Get()) return x;
  } 
  return 0.0;
}

