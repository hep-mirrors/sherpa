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
    throw(ATOOLS::Exception(ATOOLS::ex::fatal_error,
			    "Hadron remnant needs ISR Handler.",
			    "Hadron_Remnant","Hadron_Remnant"));
  }
  GetConstituents(isrhandler->Flav(m_beam));
}

void Hadron_Remnant::GetConstituents(const ATOOLS::Flavour flav) 
{
  int hadint=(flav.Kfcode()-flav.Kfcode()/10000)/10;
  if ((hadint>100)&&(hadint<1000)) {
    m_constit.resize(3);
    m_constit[0]=ATOOLS::Flavour(ATOOLS::kf::code(hadint)/100);
    m_constit[1]=ATOOLS::Flavour(ATOOLS::kf::code((hadint-
						   (hadint/100)*100)/10));
    m_constit[2]=ATOOLS::Flavour(ATOOLS::kf::code(hadint-(hadint/10)*10));
    if (flav.IsAnti()) {
      for(int i=0;i<3;i++) m_constit[i]=m_constit[i].Bar();
    }
    msg_Tracking()<<"Hadron_Remnant::FindConstituents("<<flav<<"): "
		  <<"Hadron is baryon."<<std::endl<<"   Constituents are ["
		  <<m_constit[0]<<","<<m_constit[1]<<","
		  <<m_constit[2]<<"]."<<std::endl;
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
  throw(ATOOLS::Exception(ATOOLS::ex::critical_error,
			  "Cannot determine constituents.",
			  "Hadron_Remnant","GetConstituents"));
}

bool Hadron_Remnant::FillBlob(ATOOLS::Blob *beamblob,
			      ATOOLS::Particle_List *particlelist)
{
  p_beamblob=beamblob;
  m_pbeam=beamblob->InParticle(0)->Momentum();
  m_hardpt=ATOOLS::Vec4D();
  for (size_t i=0;i<m_extracted.size();++i) {
    m_hardpt+=m_extracted[i]->Momentum();
  }
  bool success=true;
  if (!DecomposeHadron()) success=false;
  AssignRemnants();
  if (!ConnectRemnants()) success=false;
  FillRemnants();
  if (!DiceKinematics()) success=false;
  for (size_t j=0;j<m_extracted.size();++j) {
    if (particlelist!=NULL) {
      m_extracted[j]->SetNumber(-particlelist->size());
      particlelist->push_back(m_extracted[j]);
    }
  }
  for (size_t j=0;j<m_companions.size();++j) {
    m_companions[j]->SetNumber(1);
    m_companions[j]->SetInfo('F');
    if (particlelist!=NULL) {
      m_companions[j]->SetNumber(-particlelist->size());
      particlelist->push_back(m_companions[j]);
    }
  }
  return success;
}

bool Hadron_Remnant::DiceKinematics()
{
  PROFILE_HERE;
  unsigned int trials;
  ATOOLS::Vec4D ptot=m_pbeam;
  double m_xtot=1.0;
  for (unsigned int i=0;i<m_extracted.size();++i) 
    m_xtot-=m_extracted[i]->Momentum()[0]/ptot[0];
  trials=0;
  m_xscheme=1;
  double xtot;
  std::map<ATOOLS::Particle*,double> xmap;
  do {
    ++trials;
    xtot=m_xtot;
    p_pdfbase->Reset();
    for(unsigned int i=0;i<m_extracted.size();++i) {
      p_pdfbase->Extract(m_extracted[i]->Flav(),
			 2.*m_extracted[i]->Momentum()[0]/m_ecms);
    }
    for (unsigned int i=0;i<m_companions.size();++i) {
      if (!m_companions[i]->Flav().IsDiQuark()) {
	double value=1.0;
	for (unsigned int j=0;xtot-value<m_deltax && j<m_maxtrials/10;++j) { 
	  value=GetXPDF(m_companions[i]->Flav(),m_scale);
      	}
	p_pdfbase->Extract(m_companions[i]->Flav(),value);
	xmap[m_companions[i]]=value;
 	xtot-=value;
      }
      else {
	p_last[0]=m_companions[i];
      }
    }
    xmap[p_last[0]]=xtot;
    if (trials>m_maxtrials) {
      msg_Tracking()<<"Hadron_Remnant::DiceKinematics(): "
		    <<"Too many trials to find appropriate x values.\n"
		    <<"   Using naive distribution instead."<<std::endl;
      m_xscheme=0;
    }
  } while (xtot<m_deltax && m_xscheme!=0 &&
	   xmap[p_last[0]]*m_pbeam[0]<=p_last[0]->Flav().PSMass());
  p_pdfbase->Reset();
  if (m_xscheme==0) {
    xtot=0.;
    for (std::map<ATOOLS::Particle*,double>::iterator it=xmap.begin();
	 it!=xmap.end(); ++it) {
      double x=it->first->Flav().PSMass()/m_pbeam[0];
      if (x==0.) x=10.*ATOOLS::rpa.gen.Accu();
      xtot+=it->second=x;
    }
    for (std::map<ATOOLS::Particle*,double>::iterator it=xmap.begin();
	 it!=xmap.end(); ++it) {
      it->second*=m_xtot/xtot;
    }
  }
  double xperptot=0.9999999999;
  for (unsigned int i=0;i<m_extracted.size();++i) 
    ptot-=m_extracted[i]->Momentum();
  for (unsigned int j=0;j<m_companions.size();++j) {
    double E=xmap[m_companions[j]]*ptot[0];
    double m=m_companions[j]->Flav().PSMass();
    double pmax=0.9999999999*sqrt(E*E-m*m);
    double xperp=ATOOLS::Min(xperptot,pmax/ptot.PPerp());
    if (m_companions[j]->Flav().IsDiQuark() && j<m_companions.size()-1) {
      xperp=0.;
    }
    if (j==m_companions.size()-1) xperp=xperptot;
    xperptot-=xperp;
    ATOOLS::Vec4D p=xperp*ptot;
    p[0]=E;
    p[3]=ATOOLS::Sign(m_pbeam[3])*sqrt(E*E-p.PPerp2()-ATOOLS::sqr(m));
    m_companions[j]->SetMomentum(p);
    if (!(E>0.) || (!(p[3]>0.) && !(p[3]<=0.))) {
      if (!m_dupdf) {
	ATOOLS::msg.Error()<<"Hadron_Remnant::DiceKinematics(): "                 			   <<"Parton ("<<m_companions[j]<<") "
			   <<" has non-positive momentum: p = "
			   <<m_companions[j]->Momentum()<<" m_{"
			   <<m_companions[j]->Flav()<<"} = "
			   <<m_companions[j]->Flav().PSMass()<<" <- "
			   <<m_xscheme<<std::endl;
      }
      return false;
    }
  }
  return true;
}

bool Hadron_Remnant::ValenceQuark(ATOOLS::Particle *const quark) 
{
  double x=2.0*quark->Momentum()[0]/m_ecms;
  p_pdfbase->Calculate(x,0.,0.,m_scale);
  double val=p_pdfbase->GetXPDF(quark->Flav());
  return val>(p_pdfbase->GetXPDF(quark->Flav().Bar())+val)*ATOOLS::ran.Get();
}

ATOOLS::Flavour Hadron_Remnant::Opposite(ATOOLS::Flavour flav) const
{
  bool found=false;
  ATOOLS::kf::code rem[2];
  for (short unsigned int i=0,j=0;i<3;++i) {
    if (m_constit[i]==flav && !found) found=true;
    else rem[j++]=m_constit[i].Kfcode();
  }
  ATOOLS::Flavour anti=ATOOLS::Flavour(ATOOLS::kf::code(abs(rem[0])*1000+
							abs(rem[1])*100+3));
  if (rem[0]!=rem[1]) {
    if (ATOOLS::ran.Get()<0.25) 
      anti=ATOOLS::Flavour(ATOOLS::kf::code(abs(rem[0])*1000+
					    abs(rem[1])*100+1));
  }
  else {
    anti=ATOOLS::Flavour(ATOOLS::kf::code(abs(rem[0])*1100+3));
  }
  if (flav.IsAnti()) anti=anti.Bar();
  return anti;
}

bool Hadron_Remnant::DecomposeHadron() 
{
  PROFILE_HERE;
  for (ATOOLS::Particle_List::iterator pit=m_extracted.begin();
       pit!=m_extracted.end();++pit) {
    for (size_t j=0;j<m_constit.size();++j) {
      if ((*pit)->Flav()==m_constit[j]) {
	if (ValenceQuark(*pit)) {
	  p_start = new Color_Dipole(*pit,&m_companions);  
	  p_start->Begin(ANTI((*pit)->Flav().IsAnti()))->
	    SetFlav(Opposite((*pit)->Flav()));
// 	  std::cout<<"val "<<ATOOLS::om::red<<*p_start
// 		   <<ATOOLS::om::reset<<std::endl;
	  return true;
	}
      }
    }
  }
  ATOOLS::Flavour flav=m_constit[(size_t)(ATOOLS::ran.Get()*3.)];
  ATOOLS::Particle *part = new ATOOLS::Particle(-1,flav); 
  part->SetFlow(COLOR((qri::type)(flav.IsAnti())),ATOOLS::Flow::Counter());
  p_start = new Color_Dipole(part,&m_companions);  
  p_start->Begin(ANTI(flav.IsAnti()))->SetFlav(Opposite(flav));
  m_companions.push_back(part);
//   std::cout<<"sea "<<ATOOLS::om::red<<*p_start
// 	   <<ATOOLS::om::reset<<std::endl;
  return true;
}

double Hadron_Remnant::GetXPDF(ATOOLS::Flavour flavour,double scale) 
{
  PROFILE_HERE;
  double cut, x;
  cut=2.0*(flavour.PSMass()+m_hardpt.PPerp2()/
	   ATOOLS::sqr(m_companions.size()))/sqrt(scale);
  if (scale<p_pdfbase->Q2Min()) {
    ATOOLS::msg.Error()<<"Hadron_Remnant::GetXPDF("<<flavour<<","<<scale<<"): "
		       <<"Scale under-runs minimum given by PDF: "
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
      for (unsigned int j=0, i=0;i<3;i++) 
	if (i!=single) di[j++]=m_constit[i].Kfcode();
      if (di[0]!=di[1]) {
	ATOOLS::Flavour singlet=
	  ATOOLS::Flavour(ATOOLS::kf::code(abs(di[0])*1000+abs(di[1])*100+1));
	ATOOLS::Flavour triplet=
	  ATOOLS::Flavour(ATOOLS::kf::code(abs(di[0])*1000+abs(di[1])*100+3));
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
	  ATOOLS::Flavour singlet=
	    ATOOLS::Flavour(ATOOLS::kf::code(abs(di[0])*1000+
					     abs(di[1])*100+1));
	  ATOOLS::Flavour triplet=
	    ATOOLS::Flavour(ATOOLS::kf::code(abs(di[0])*1000+
					     abs(di[1])*100+3));
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
	  ATOOLS::Flavour singlet=
	    ATOOLS::Flavour(ATOOLS::kf::code(abs(di[0])*1000+
					     abs(di[1])*100+1));
	  ATOOLS::Flavour triplet=
	    ATOOLS::Flavour(ATOOLS::kf::code(abs(di[0])*1000+
					     abs(di[1])*100+3));
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

