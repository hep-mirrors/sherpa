#include "PDF/Remnant/Hadron_Remnant.H"

#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "BEAM/Main/Beam_Base.H"
#include "ATOOLS/Math/Random.H"

using namespace PDF;
using namespace ATOOLS;

Hadron_Remnant::Hadron_Remnant(PDF::ISR_Handler *isrhandler,
			       const unsigned int beam):
  QCD_Remnant_Base(isrhandler,beam,rtp::hadron)
{
  if (isrhandler==NULL) {
    THROW(fatal_error,"Hadron remnant needs ISR Handler.");
  }
  GetConstituents(isrhandler->Flav(m_beam));
  m_emin=0.0;
}

const std::vector<ATOOLS::Flavour> &Hadron_Remnant::
GetConstituents(const ATOOLS::Flavour flav) 
{
  if (m_constit.size()>0) return m_constit;
  int hadint=(flav.Kfcode()-flav.Kfcode()/10000)/10;
  if ((hadint>100)&&(hadint<1000)) {
    m_constit.resize(3);
    m_constit[0]=Flavour((kf_code)(hadint)/100);
    m_constit[1]=Flavour((kf_code)((hadint-(hadint/100)*100)/10));
    m_constit[2]=Flavour((kf_code)(hadint-(hadint/10)*10));
    if (flav.IsAnti()) {
      for(int i=0;i<3;i++) m_constit[i]=m_constit[i].Bar();
    }
    msg_Tracking()<<"Hadron_Remnant::FindConstituents("<<flav<<"): "
		  <<"Hadron is baryon."<<std::endl<<"   Constituents are ["
		  <<m_constit[0]<<","<<m_constit[1]<<","
		  <<m_constit[2]<<"]."<<std::endl;
    return m_constit;
  }
  if ((hadint>10)&&(hadint<100)) {
    m_constit.resize(2);
    m_constit[0]=Flavour((kf_code)(hadint)/10);
    m_constit[1]=Flavour((kf_code)(hadint-(hadint/10)*10));
    if (flav.IsAnti()) {
      for(int i=0;i<2;i++) m_constit[i]=m_constit[i].Bar();
    }
    msg_Tracking()<<"Hadron_Remnant::FindConstituents("<<flav<<"): "
		  <<"Hadron is meson."<<std::endl<<"   Constituents are ["
		  <<m_constit[0]<<","<<m_constit[1]<<"]."<<std::endl;
    return m_constit;
  }
  THROW(critical_error,"Cannot determine constituents.");
  return m_constit;
}

bool Hadron_Remnant::FillBlob(ATOOLS::Blob *beamblob,
			      ATOOLS::Particle_List *particlelist)
{
  p_beamblob=beamblob;
  m_pbeam=beamblob->InParticle(0)->Momentum();
  m_hardpt=Vec4D();
  for (size_t i=0;i<m_extracted.size();++i) {
    m_hardpt+=m_extracted[i]->Momentum();
  }
  bool success=true;
  if (!DecomposeHadron()) success=false;
  AssignRemnants();
  FillRemnants();
  if (!GenerateKinematics()) success=false;
  for (size_t j=0;j<m_extracted.size();++j) {
    if (particlelist!=NULL) {
      m_extracted[j]->SetNumber(-particlelist->size());
      particlelist->push_back(m_extracted[j]);
    }
    else m_extracted[j]->SetNumber(0);
  }
  for (size_t j=0;j<m_companions.size();++j) {
    m_companions[j]->SetNumber(1);
    m_companions[j]->SetInfo('B');
    if (particlelist!=NULL) {
      m_companions[j]->SetNumber(-particlelist->size());
      particlelist->push_back(m_companions[j]);
    }
    else m_companions[j]->SetNumber(0);
  }
  return success;
}

bool Hadron_Remnant::GenerateKinematics()
{
  unsigned int trials;
  Vec4D ptot=m_pbeam;
  double m_xtot=1.0;
  for (unsigned int i=0;i<m_extracted.size();++i) 
    m_xtot-=m_extracted[i]->Momentum()[0]/ptot[0];
  trials=0;
  m_xscheme=1;
  std::map<Particle*,double> xmap;
  do {
    ++trials;
    m_xrem=m_xtot;
    p_pdfbase->Reset();
    for(unsigned int i=0;i<m_extracted.size();++i) {
      p_pdfbase->Extract(m_extracted[i]->Flav(),
			 m_extracted[i]->Momentum()[0]/p_beam->Energy());
    }
    for (unsigned int i=0;i<m_companions.size();++i) {
      if (!m_companions[i]->Flav().IsDiQuark()) {
	double value=1.0;
	for (unsigned int j=0;m_xrem-value<m_deltax && j<m_maxtrials/10;++j) { 
	  value=GetXPDF(m_companions[i]->Flav(),m_scale);
      	}
	p_pdfbase->Extract(m_companions[i]->Flav(),value);
	xmap[m_companions[i]]=value;
 	m_xrem-=value;
      }
      else {
	p_last[0]=m_companions[i];
      }
    }
    xmap[p_last[0]]=m_xrem;
    if (trials>m_maxtrials) m_xscheme=0;
  } while (m_xrem<m_deltax && m_xscheme!=0 &&
	   xmap[p_last[0]]*m_pbeam[0]<=p_last[0]->Flav().Mass());
  p_pdfbase->Reset();
  if (m_xscheme==0) {
    double xtot=0.;
    for (std::map<Particle*,double>::iterator it=xmap.begin();
	 it!=xmap.end(); ++it) {
      double x=1.1*it->first->Flav().Mass()/m_pbeam[0];
      if (x==0.0) x=10.*rpa->gen.Accu();
      xtot+=it->second=x;
    }
    for (std::map<Particle*,double>::iterator it=xmap.begin();
	 it!=xmap.end(); ++it) {
      it->second*=m_xtot/xtot;
    }
  }
  double xperptot=1.0-1.0e-12;
  for (unsigned int i=0;i<m_extracted.size();++i) 
    ptot-=m_extracted[i]->Momentum();
  for (unsigned int j=0;j<m_companions.size();++j) {
    double E=xmap[m_companions[j]]*m_pbeam[0];
    double m=m_companions[j]->Flav().Mass();
    if (p_partner->Type()==rtp::electron) m=0.0;
    double pmax=(1.0-1.0e-12)*sqrt(E*E-m*m);
    double xperp=Min(xperptot,pmax/ptot.PPerp());
    xperptot-=xperp;
    Vec4D p=xperp*ptot;
    p[0]=E;
    p[3]=Sign(m_pbeam[3])*sqrt(E*E-p.PPerp2()-sqr(m));
    m_companions[j]->SetMomentum(p);
    if (!(E>0.) || (!(p[3]>0.) && !(p[3]<=0.))) {
      msg_Tracking()<<"Hadron_Remnant::GenerateKinematics(): "
			 <<"Parton ("<<m_companions[j]<<") "
			 <<" has non-positive momentum: p = "
			 <<m_companions[j]->Momentum()<<" m_{"
			 <<m_companions[j]->Flav()<<"} = "
			 <<m_companions[j]->Flav().Mass()<<" <- "
			 <<m_xscheme<<std::endl;
      return false;
    }
  }
  return true;
}

bool Hadron_Remnant::ValenceQuark(Particle *const quark) 
{
  double x=quark->Momentum()[0]/p_beam->Energy();
  if (x>1.) {
    msg_Out()<<" WARNING in Hadron_Remnant::ValenceQuark \n"
	     <<" (x-1)="<<x-1<<std::endl;
    x = 1.;
  }
  if (x<p_pdfbase->XMin() || x>p_pdfbase->XMax()) return false;
  if (m_scale<p_pdfbase->Q2Min()) m_scale=1.001*p_pdfbase->Q2Min();
  p_pdfbase->Calculate(x,m_scale);
  double val=p_pdfbase->GetXPDF(quark->Flav());
  return val>(p_pdfbase->GetXPDF(quark->Flav().Bar())+val)*ran->Get();
}

ATOOLS::Flavour Hadron_Remnant::Opposite(ATOOLS::Flavour flav) const
{
  bool found=false;
  kf_code rem[2];
  for (short unsigned int i=0,j=0;i<3;++i) {
    if (m_constit[i]==flav && !found) found=true;
    else rem[j++]=m_constit[i].Kfcode();
  }
  Flavour anti=Flavour((kf_code)(abs(rem[0])*1000+abs(rem[1])*100+3));
  if (rem[0]!=rem[1]) {
    if (ran->Get()<0.25) 
      anti=Flavour((kf_code)(abs(rem[0])*1000+abs(rem[1])*100+1));
  }
  else {
    anti=Flavour((kf_code)(abs(rem[0])*1100+3));
  }
  if (flav.IsAnti()) anti=anti.Bar();
  return anti;
}

bool Hadron_Remnant::DecomposeHadron() 
{
  bool success=true;
  double Eb(p_beam->Energy());
  for (Particle_List::iterator pit=m_extracted.begin();
       pit!=m_extracted.end();++pit) {
    if ((*pit)->Momentum()[0]>Eb || (*pit)->Momentum()[0]<0.0) {
      msg_Error()<<"Hadron_Remnant::DecomposeHadron(): "
			 <<"Constituent energy out of range. \n   E_"
			 <<(*pit)->Flav()<<" = "<<(*pit)->Momentum()[0]
			 <<"."<<std::endl;
      success=false;
    }
    for (size_t j=0;j<m_constit.size();++j) {
      if ((*pit)->Flav()==m_constit[j]) {
	//std::cout<<METHOD<<" "<<success<<":"<<(*pit)->Flav()
	//	 <<" ("<<ValenceQuark(*pit)<<")"<<std::endl;
	if (success && ValenceQuark(*pit)) {
	  p_start = new Color_Dipole(*pit,&m_companions);  
	  p_start->Begin(ANTI((*pit)->Flav().IsAnti()))->
	    SetFlav(Opposite((*pit)->Flav()));
	  return success;
	}
      }
    }
  }
  Flavour    flav = m_constit[(size_t)(ran->Get()*3.)];
  Particle * part = new Particle(-1,flav); 
  part->SetStatus(part_status::active);
  part->SetFinalMass(flav.Mass());
  part->SetFlow(COLOR((qri::type)(flav.IsAnti())),Flow::Counter());
  //std::cout<<METHOD<<":"<<flav<<std::endl
  //	   <<"  "<<(*part)<<std::endl;
  p_start = new Color_Dipole(part,&m_companions);  
  p_start->Begin(ANTI(flav.IsAnti()))->SetFlav(Opposite(flav));
  m_companions.push_back(part);
  return success;
}

double Hadron_Remnant::GetXPDF(ATOOLS::Flavour flavour,double scale) 
{
  double cut, x;
  cut=2.0*(flavour.HadMass()+m_hardpt.PPerp()/
	   sqr(m_companions.size()))/p_beam->OutMomentum()[0];
  // assume heavy flavours have been pair-produced
  // => scale should be approximately (2m)^2
  scale=Max(scale,4.0*sqr(flavour.Mass()));
  if (scale<p_pdfbase->Q2Min()) {
    msg_Error()<<"Hadron_Remnant::GetXPDF("<<flavour<<","<<scale<<"): "
		       <<"Scale under-runs minimum given by PDF: "
		       <<scale<<" < "<<p_pdfbase->Q2Min()<<std::endl;
    scale=1.001*p_pdfbase->Q2Min();
  } 
  unsigned int xtrials, pdftrials=0;
  while (true) {
    ++pdftrials;
    xtrials=0;
    do { 
      ++xtrials;
      x=m_xrem*ran->Get();
      if (xtrials>=m_maxtrials) {
	x=Min(cut,0.999999*p_pdfbase->RescaleFactor());
	break;
      }
      if (x>p_pdfbase->RescaleFactor()) continue;
    } while (x<cut);	
    if (x>p_pdfbase->XMin() && x<p_pdfbase->XMax()) {
      p_pdfbase->Calculate(x,scale);
    }
    else {
      m_xscheme=0; return 0.01;
    }
    if (pdftrials>=m_maxtrials) { m_xscheme=0; return 0.01; }
    if (p_pdfbase->GetXPDF(flavour)/x>ran->Get()) return x;
  } 
  return 0.0;
}

double Hadron_Remnant::MinimalEnergy(const ATOOLS::Flavour &flavour) 
{
  if (!m_initialized) {
    if (!flavour.Strong()) {
      return p_beam->Beam().HadMass();
    }
    bool found(false);
    kf_code di[3] = {0, 0, 0};
    if (flavour.IsQuark()) {
      short unsigned int j=0;
      for (Flavour_Vector::const_iterator flit(m_constit.begin());
           flit!=m_constit.end();
           ++flit) {
        if (found || flavour!=*flit) di[j++]=flit->Kfcode();
        else found=true;
      }
    }
    Flavour difl;
    if (!found || flavour.IsGluon()) {
      int single=-1;
      for (size_t i=0;i<m_constit.size();++i) {
        for (size_t j=i+1;j<m_constit.size();++j) {
          if (m_constit[i]==m_constit[j]) { 
            single=j; 
            break; 
          }
        }
        if (single>0) break;
      }
      Flavour fl(m_constit[single]);
      for (short unsigned int j=0, i=0;i<3;i++) 
        if (i!=single) di[j++]=m_constit[i].Kfcode();
      if (di[0]>di[1]) difl=Flavour((kf_code)(di[0]*1000+di[1]*100+1));
      else if (di[1]>di[0]) difl=Flavour((kf_code)(di[1]*1000+di[0]*100+1));
      else difl=Flavour((kf_code)(di[0]*1100+3));
      if (m_constit[single].IsAnti()) difl=difl.Bar();
      return difl.Mass()+fl.Mass()+flavour.Bar().Mass();
    }
    if (di[0]>di[1]) difl=Flavour((kf_code)(di[0]*1000+di[1]*100+1));
    else if (di[1]>di[0]) difl=Flavour((kf_code)(di[1]*1000+di[0]*100+1));
    else difl=Flavour((kf_code)(di[0]*1100+3));
    if (m_constit[0].IsAnti()) difl=difl.Bar();
    return difl.Mass();
  } else {
    if (flavour.IsQuark()) return flavour.Bar().Mass();
  }
  return 0.;
}

