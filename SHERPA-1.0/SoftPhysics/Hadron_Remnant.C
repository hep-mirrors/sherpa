#include "Hadron_Remnant.H"

#include "Run_Parameter.H"
#include "Exception.H"
#include "Random.H"

#ifdef PROFILE__all
#define PROFILE__Hadron_Remnant
#endif
#ifdef PROFILE__Hadron_Remnant
#include "prof.hh" 
#else
#define PROFILE_HERE
#endif

using namespace SHERPA;

Hadron_Remnant::Hadron_Remnant(PDF::ISR_Handler *isrhandler,
			       const unsigned int beam,const double scale):
  QCD_Remnant_Base(isrhandler,beam,-scale,rtp::hadron)
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
    msg_Tracking()<<"Hadron_Remnant::FindConstituents("<<flav<<"): "
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
    msg_Tracking()<<"Hadron_Remnant::FindConstituents("<<flav<<"): "
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
  m_hardpt=ATOOLS::Vec4D();
  for (size_t i=0;i<m_parton[1].size();++i) {
    p_beamblob->AddToOutParticles(m_parton[1][i]);
    m_hardpt+=m_parton[1][i]->Momentum();
  }
  // decompose hadron
  SortRemnants();
  msg_Debugging()<<*p_beamblob<<std::endl;
  if (!SelectCompanions()) return false;
  if (m_initial>1) if (!ConnectRemnants()) return false;
  if (!AttachLastRemnants()) return false;
  // select x's according to pdf
  DiceKinematics();
  // fill blob
  for (int i=1;i>=0;--i) {
    for (size_t j=0;j<m_parton[i].size();++j) {
      if (i==0) {
	m_parton[i][j]->SetNumber(1);
	m_parton[i][j]->SetInfo('F');
      }
      if (particlelist!=NULL) {
	m_parton[i][j]->SetNumber(-particlelist->size());
	particlelist->push_back(m_parton[i][j]);
      }
    }
  }
  msg_Debugging()<<*p_beamblob<<p_beamblob->CheckMomentumConservation()<<std::endl;
  return true;
}

void Hadron_Remnant::DiceKinematics()
{
  PROFILE_HERE;
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
	for (unsigned int j=0;xtot-value<m_deltax && j<m_maxtrials/10;++j) { 
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
      msg_Tracking()<<"Hadron_Remnant::DiceKinematics(): "
			    <<"Too many trials to find appropriate x values for partons."<<std::endl
			    <<"   Using naive distribution instead."<<std::endl;
      m_xscheme=0;
    }
  } while (xtot<m_deltax && m_xscheme!=0 &&
	   xmap[p_last[0]]*m_pbeam[0]<=p_last[0]->Flav().PSMass());
  p_pdfbase->Reset();
  if (m_xscheme==0) {
    xtot=0.;
    for (std::map<ATOOLS::Particle*,double>::iterator it=xmap.begin();it!=xmap.end(); ++it) {
      double x=it->first->Flav().PSMass()/m_pbeam[0];
      if (x==0.) x=10.*ATOOLS::rpa.gen.Accu();
      xtot+=it->second=x;
    }
    for (std::map<ATOOLS::Particle*,double>::iterator it=xmap.begin();it!=xmap.end(); ++it) {
      it->second*=m_xtot/xtot;
    }
  }
  double xperptot=0.9999999999;
  for (unsigned int i=0;i<m_parton[1].size();++i) ptot-=m_parton[1][i]->Momentum();
  for (unsigned int j=0;j<m_parton[0].size();++j) {
    double E=xmap[m_parton[0][j]]*ptot[0];
    double m=m_parton[0][j]->Flav().PSMass();
    double pmax=0.9999999999*sqrt(E*E-m*m);
    double xperp=ATOOLS::Min(xperptot,pmax/ptot.PPerp());
    if (j==m_parton[0].size()-1) xperp=xperptot;
    xperptot-=xperp;
    ATOOLS::Vec4D p=xperp*ptot;
    p[0]=E;
    p[3]=ATOOLS::Sign(m_pbeam[3])*sqrt(E*E-p.PPerp2()-ATOOLS::sqr(m));
    m_parton[0][j]->SetMomentum(p);
    if (!(E>0.) || (!(p[3]>0.) && !(p[3]<=0.))) {
      ATOOLS::msg.Error()<<"Hadron_Remnant::DiceKinematics(): "                 
                         <<"Parton ("<<m_parton[0][j]<<") has non-positive momentum: p = "
                         <<m_parton[0][j]->Momentum()<<" m_{"<<m_parton[0][j]->Flav()<<"} = "
                         <<m_parton[0][j]->Flav().PSMass()<<" <- "<<m_xscheme<<std::endl;
    }
  }
}

bool Hadron_Remnant::AttachLastRemnants() 
{
  PROFILE_HERE;
  m_last=0;
  short int pos=-1;
  QCD_Remnant_Info *constit=NULL;
  for (size_t i=0;i<m_sorted.size();++i) {
    for (size_t j=0;j<m_constit.size();++j) {
      if ((*m_sorted[i])->Flav()==m_constit[j]) {
	pos=j;
	constit=m_sorted[i];
      }
    }
  }
  if (pos<0) pos=(short int)(ATOOLS::ran.Get()*3.); 
  short int rem[2];
  for (short unsigned int i=0,j=0;i<3;++i) {
    if (i!=pos) rem[j++]=m_constit[i].Kfcode();
  }
  ATOOLS::Flavour part=m_constit[pos];
  ATOOLS::Flavour anti=ATOOLS::Flavour(ATOOLS::kf::code(abs(rem[0])*1000+abs(rem[1])*100+3));
  if (rem[0]!=rem[1]) {
    if (ATOOLS::ran.Get()<0.25) anti=ATOOLS::Flavour(ATOOLS::kf::code(abs(rem[0])*1000+abs(rem[1])*100+1));
  }
  else {
    anti=ATOOLS::Flavour(ATOOLS::kf::code(abs(rem[0])*1100+3));
  }
  if (m_constit[0].IsAnti()) anti=anti.Bar();
  if (constit==NULL) {
    constit=m_sorted[m_sorted.size()-1];
    if (constit!=NULL) msg_Debugging()<<"c<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n"<<*constit<<std::endl;
    ATOOLS::Particle *newpart[3];
    newpart[0] = new ATOOLS::Particle(-1,part); 
    newpart[1] = new ATOOLS::Particle(-1,anti); 
    if (newpart[0]->Flav().IsAnti()) std::swap<ATOOLS::Particle*>(newpart[0],newpart[1]);
    newpart[0]->SetFlow(1,ATOOLS::Flow::Counter());
//     newpart[0]->SetFlow(1,(*constit)(1)->GetFlow(2));
    p_beamblob->RemoveOutParticle((*constit)(1));
    msg_Debugging()<<"-> "<<*newpart[0]<<"\n-> "<<*newpart[1]<<"\n-> "<<*(*(*constit)[1])(0)<<"\n-> "<<*(*(*constit)[0])(1)<<std::endl;
    bool singlet=false;
    unsigned int old=(*constit)(1)->GetFlow(2);
    if (!AdjustColours((*constit)(1),old,newpart[0]->GetFlow(1),singlet,false)) return false;
    msg_Debugging()<<"singlet: "<<old<<" "<<newpart[0]->GetFlow(1)<<" "<<singlet<<" "<<*(*constit)(1)<<std::endl;
    if (singlet) {
      newpart[2] = new ATOOLS::Particle(-1,ATOOLS::kf::gluon);
      m_parton[0].push_back(newpart[2]);
      newpart[2]->SetFlow(1,old);
      newpart[2]->SetFlow(2,newpart[0]->GetFlow(1));
      p_beamblob->AddToOutParticles(newpart[2]);
      ++m_last;
      msg_Debugging()<<"inserted "<<newpart[2]<<std::endl;
    }
    newpart[1]->SetFlow(2,(*(*constit)[1])(0)->GetFlow(1)); 
    p_beamblob->AddToOutParticles((*constit)(1));
    for (unsigned int i=0;i<2;++i) {
      m_parton[0].push_back(newpart[i]);
      p_beamblob->AddToOutParticles(newpart[i]);
    }
    p_last[0]=newpart[1];
    m_last+=2;
  }
  else {
    //  if (constit!=NULL) msg_Debugging()<<"f<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n"<<*constit<<std::endl;
    unsigned int pos=anti.IsAnti();
    (*constit)(1-pos)->SetFlav(anti);
    if (m_initial==1) (*constit)(1-pos)->SetFlow(2-pos,(*constit)(pos)->GetFlow(1+pos));
    p_last[0]=(*constit)(1-pos);
  }
  return true;
}

double Hadron_Remnant::GetXPDF(ATOOLS::Flavour flavour,double scale) 
{
  PROFILE_HERE;
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

double Hadron_Remnant::MinimalEnergy(const ATOOLS::Flavour &flavour) 
{
  if (!m_initialized) {
    m_initialized=true;
    if (flavour.IsGluon()) {
      size_t single=(size_t)(ATOOLS::ran.Get()*3.); 
      ATOOLS::Flavour difl, fl=m_constit[single];
      int di[2];
      for (unsigned int j=0, i=0;i<3;i++) if (i!=single) di[j++]=m_constit[i].Kfcode();
      if (di[0]!=di[1]) {
	ATOOLS::Flavour singlet=ATOOLS::Flavour(ATOOLS::kf::code(abs(di[0])*1000+abs(di[1])*100+1));
	ATOOLS::Flavour triplet=ATOOLS::Flavour(ATOOLS::kf::code(abs(di[0])*1000+abs(di[1])*100+3));
	if (singlet.PSMass()>triplet.PSMass()) difl=singlet;
	else difl=triplet;
      }
      else {
	difl=ATOOLS::Flavour(ATOOLS::kf::code(abs(di[0])*1100+3));
      }
      if (m_constit[0].IsAnti()) difl=difl.Bar();
      return difl.PSMass()+fl.PSMass();
    }
    else if (flavour.IsQuark()) {
      bool found=false;
      int di[3];
      for (size_t j=0,i=0;i<m_constit.size();++i) {
	if (found||flavour!=m_constit[i]) di[j++]=m_constit[i];
	else found=true;
      }
      if (found) {
	ATOOLS::Flavour difl;
	if (di[0]!=di[1]) {
	  ATOOLS::Flavour singlet=ATOOLS::Flavour(ATOOLS::kf::code(abs(di[0])*1000+abs(di[1])*100+1));
	  ATOOLS::Flavour triplet=ATOOLS::Flavour(ATOOLS::kf::code(abs(di[0])*1000+abs(di[1])*100+3));
	  if (singlet.PSMass()>triplet.PSMass()) difl=singlet;
	  else difl=triplet;
	}
	else {
	  difl=ATOOLS::Flavour(ATOOLS::kf::code(abs(di[0])*1100+3));
	}
	if (m_constit[0].IsAnti()) difl=difl.Bar();
	return difl.PSMass();
      }
      else {
	unsigned int single=(unsigned int)(ATOOLS::ran.Get()*3.0); 
	ATOOLS::Flavour difl, fl=m_constit[single];
	for (unsigned int j=0,i=0;i<3;i++) {
	  if (i!=single) di[j++]=m_constit[i].Kfcode();
	}
	if (di[0]!=di[1]) {
	  ATOOLS::Flavour singlet=ATOOLS::Flavour(ATOOLS::kf::code(abs(di[0])*1000+abs(di[1])*100+1));
	  ATOOLS::Flavour triplet=ATOOLS::Flavour(ATOOLS::kf::code(abs(di[0])*1000+abs(di[1])*100+3));
	  if (singlet.PSMass()>triplet.PSMass()) difl=singlet;
	  else difl=triplet;
	}
	else {
	  difl=ATOOLS::Flavour(ATOOLS::kf::code(di[0]*1100+3));
	}
	if (m_constit[0].IsAnti()) difl=difl.Bar();
	return difl.PSMass()+fl.PSMass()+flavour.Bar().PSMass();
      }
    }
  }
  else {
    if (flavour.IsQuark()) return flavour.Bar().PSMass();
  }
  return 0.;
}

void Hadron_Remnant::UnDo() 
{
//   for (size_t i=0;i<m_last;++i) {
//     p_beamblob->DeleteOutParticle(m_parton[0][m_parton[0].size()-i-1]);
//   }
//   m_parton[0].resize(m_parton[0].size()-m_last);
//   QCD_Remnant_Base::UnDo();
}
