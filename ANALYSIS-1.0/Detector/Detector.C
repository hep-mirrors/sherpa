#include "Detector.H"
#include "Rambo.H"
#include "Random.H"
#include "MyStrStream.H"
#include "Message.H"
#include "CXXFLAGS.H"

#ifdef USING__ROOT
#include "Scaling.H"
#include "TH2D.h"
#include "My_Root.H"
#endif 

using namespace ANALYSIS;
using namespace ATOOLS;

Detector::Detector(Primitive_Analysis * ana) : 
  p_ana(ana), m_inlist("FinalState"), m_outlist("Detected_FS")
{
  m_name = "Detector";
}

Detector::~Detector() {}

void Detector::AddDetectorElement(Detector_Element * element) {
  std::string name = element->Name();
  if (m_elements.find(name)!=m_elements.end()) {
    msg.Error()<<"Error in "<<METHOD<<" : "<<std::endl
	       <<"   Element of type '"<<name<<"' already in detector."<<std::endl
	       <<"   Abort the run."<<std::endl;
    abort();
  }
  m_elements[name] = element;
}

Detector_Element * Detector::GetElement(std::string name) {
  if (m_elements.find(name)==m_elements.end()) {
    msg.Error()<<"Error in "<<METHOD<<" : "<<std::endl
	       <<"   Element of type '"<<name<<"' does not exist in detector."<<std::endl
	       <<"   Return 'NULL' and hope for the best."<<std::endl;
    return NULL;
  }
  return m_elements[name];
}

void Detector::AddParticleSmearer(Particle_Smearer_Base * smearer) {
  std::string name = smearer->Name();
  if (m_smearers.find(name)!=m_smearers.end()) {
    msg.Error()<<"Error in "<<METHOD<<" : "<<std::endl
	       <<"   Smearer of type '"<<name<<"' already in detector."<<std::endl
	       <<"   Abort the run."<<std::endl;
    abort();
  }
  m_smearers[name] = smearer;
} 

Particle_Smearer_Base * Detector::GetParticleSmearer(std::string name) {
  if (m_smearers.find(name)==m_smearers.end()) {
    msg.Error()<<"Error in "<<METHOD<<" : "<<std::endl
	       <<"   Element of type '"<<name<<"' does not exist in detector."<<std::endl
	       <<"   Return 'NULL' and hope for the best."<<std::endl;
    return NULL;
  }
  return m_smearers[name];
}

void Detector::AddObjectDefinition(Object_Definition_Base * definition) {
  std::string name = definition->Name();
  if (m_definitions.find(name)!=m_definitions.end()) {
    msg.Error()<<"Error in "<<METHOD<<" : "<<std::endl
	       <<"   Definition of type '"<<name<<"' already in detector."<<std::endl
	       <<"   Abort the run."<<std::endl;
    abort();
  }
  m_definitions[name] = definition;
} 

Object_Definition_Base * Detector::GetObjectDefinition(std::string name) {
  if (m_definitions.find(name)==m_definitions.end()) {
    msg.Error()<<"Error in "<<METHOD<<" : "<<std::endl
	       <<"   Element of type '"<<name<<"' does not exist in detector."<<std::endl
	       <<"   Return 'NULL' and hope for the best."<<std::endl;
    return NULL;
  }
  return m_definitions[name];
}

void Detector::Evaluate(const ATOOLS::Blob_List &bl, double weight, int ncount)
{
  std::cout<<METHOD<<" : "<<p_ana<<std::endl;
  //Print();
  Particle_List *outlist(new Particle_List);
  Particle_List *inlist(p_ana->GetParticleList(m_inlist));
  //std::cout<<METHOD<<": After list definitions."<<std::endl;
  if (inlist==NULL) {
    msg.Error()<<METHOD<<"(): List '"<<m_inlist<<"' not found."<<std::endl;
    p_ana->AddParticleList(m_outlist,outlist);
    return;
  }
  Fill(inlist);
  ReconstructObjects(outlist);
  p_ana->AddParticleList(m_outlist,outlist);
  //std::cout<<METHOD<<" : "<<outlist->size()<<std::endl
  //	   <<"==========================================================="<<std::endl;
  Reset();
}


void Detector::Fill(Particle_List * plist) {
  //std::cout<<METHOD<<std::endl<<(*plist)<<std::endl;
  int ehits=0,hhits=0;
  Particle * part(NULL);
  double eta,phi,E(0);
  for (size_t i=0;i<plist->size();i++) {
    part = (*plist)[i];
    for (std::map<std::string,Particle_Smearer_Base *>::iterator smit=m_smearers.begin();
	 smit!=m_smearers.end();smit++) {
      if (smit->second->TreatParticle(part)) {
	if (smit->second->GivesTrack(eta,phi)) 
	  m_elements["Tracker"]->Fill(E,eta,phi,part);
	E = smit->second->EnergyInECal(eta,phi);
	if (E>0) { ehits+=int(m_elements["ECal"]->Fill(E,eta,phi,part)); }
	E = smit->second->EnergyInHCal(eta,phi);
	if (E>0) { hhits+=int(m_elements["HCal"]->Fill(E,eta,phi,part)); }
	if (smit->second->GivesMuon(eta,phi)) 
	  m_elements["Muon_Chambers"]->Fill(part->Momentum()[0],eta,phi,part);
	continue;
      }
    }
  }
  std::cout<<METHOD<<" ehits = "<<ehits<<", hhits = "<<hhits<<", total = "<<(ehits+hhits)<<std::endl;
}

void Detector::ReconstructObjects(Particle_List *& plist) {
  for (std::map<std::string,Object_Definition_Base *>::iterator defit=m_definitions.begin();
       defit!=m_definitions.end();defit++) {
    //std::cout<<METHOD<<" for "<<defit->first<<std::endl;
    defit->second->ReconstructObjects(plist);
  }
}

void Detector::Reset() {
  for (std::map<std::string,Detector_Element *>::iterator elit=m_elements.begin();
       elit!=m_elements.end();elit++) elit->second->Reset();
  for (std::map<std::string,Object_Definition_Base *>::iterator defit=m_definitions.begin();
       defit!=m_definitions.end();defit++) defit->second->Reset();
}

void Detector::Print() {
  msg_Out()<<METHOD
	   <<" for "<<m_elements.size()<<" detector elements and "
	   <<m_definitions.size()<<" definitions."<<std::endl;
  for (std::map<std::string,Detector_Element *>::iterator elit=m_elements.begin();
       elit!=m_elements.end();elit++) {
    msg_Out()<<"  Detector Element  : "<<elit->first<<" "<<elit->second<<std::endl;
  }
  for (std::map<std::string,Object_Definition_Base *>::iterator defit=m_definitions.begin();
       defit!=m_definitions.end();defit++) {
    msg_Out()<<"  Object Definition : "<<defit->first<<" "<<defit->second<<std::endl;
  }
}


/*-----------------------------------------------------------------------------------------
 *
 *
 *  TESTING ROUTINES
 *
 *
 *
 --------------------------------------------------------------------------------------------*/

void Detector::Test(const int mode) {
  int testevents(1),number(20);
  Flavour flav(Flavour(kf::none));
  Print();
  InitHistograms(mode);
  switch (mode) {
  case 1:
    for (long int i=0;i<testevents;i++) TestRandomIsotropicEvent(number);
    break;
  case 2:
    testevents = 100; number = 100;
    for (int j=0;j<10;j++) {
      flav = Flavour(kf::e);
      msg_Out()<<"================================================================"<<std::endl
  	       <<"================================================================"<<std::endl
  	       <<"================================================================"<<std::endl
  	       <<"Testing "<<m_elements.size()<<" detector elements "
  	       <<"with "<<flav<<" at "<<Energy(j)<<" GeV."<<std::endl;
      for (long int i=0;i<testevents;i++) TestIsotropicEvent(number,j,flav);
      flav = Flavour(kf::photon);
      msg_Out()<<"================================================================"<<std::endl
 	       <<"================================================================"<<std::endl
	       <<"================================================================"<<std::endl
	       <<"Testing "<<m_elements.size()<<" detector elements "
	       <<"with "<<flav<<" at "<<Energy(j)<<" GeV."<<std::endl;
      for (long int i=0;i<testevents;i++) TestIsotropicEvent(number,j,Flavour(kf::photon));
      flav = Flavour(kf::pi);
      msg_Out()<<"================================================================"<<std::endl
 	       <<"================================================================"<<std::endl
 	       <<"================================================================"<<std::endl
 	       <<"Testing "<<m_elements.size()<<" detector elements "
  	       <<"with "<<flav<<" at "<<Energy(j)<<" GeV."<<std::endl;
      for (long int i=0;i<testevents;i++) TestIsotropicEvent(number,j,Flavour(kf::pi));
    } 
    break;
  case 3:
    for (int j=0;j<10;j++) {
      flav = Flavour(kf::e);
      msg_Out()<<"================================================================"<<std::endl
  	       <<"================================================================"<<std::endl
  	       <<"================================================================"<<std::endl
  	       <<"Testing "<<m_definitions.size()<<" reconstruction codes "
  	       <<"with "<<flav<<" at "<<Energy(j)<<" GeV."<<std::endl;
      for (long int i=0;i<testevents;i++) TestReconstructionCodes(number,j,flav,false);
    }
  case 4:
    for (int j=0;j<10;j++) {
      flav = Flavour(kf::e);
      msg_Out()<<"================================================================"<<std::endl
  	       <<"================================================================"<<std::endl
  	       <<"================================================================"<<std::endl
  	       <<"Testing "<<m_definitions.size()<<" reconstruction codes "
  	       <<"with "<<flav<<" at "<<Energy(j)<<" GeV."<<std::endl;
      for (long int i=0;i<testevents;i++) TestReconstructionCodes(number,j,flav,true);
    }
    break;
  }
  msg_Out()<<"================================================================"<<std::endl
	   <<"================================================================"<<std::endl
	   <<"================================================================"<<std::endl;
}

void Detector::TestRandomIsotropicEvent(const int number) {
  msg_Out()<<"================================================================"<<std::endl
	   <<"================================================================"<<std::endl
	   <<"================================================================"<<std::endl;
  for (std::map<std::string,Detector_Element *>::iterator elit=m_elements.begin();
       elit!=m_elements.end();elit++) elit->second->Reset();
  Particle_List * plist = ProduceParticleList("Rambo",number,1000,Flavour(kf::none));
  Fill(plist);
  for (std::map<std::string,Detector_Element *>::iterator elit=m_elements.begin();
       elit!=m_elements.end();elit++) elit->second->PrintHits();
  plist->Clear();
  delete plist;
}

void Detector::TestIsotropicEvent(const int number,const int j,const ATOOLS::Flavour flav) {
  for (std::map<std::string,Detector_Element *>::iterator elit=m_elements.begin();
       elit!=m_elements.end();elit++) elit->second->Reset();
  double E(Energy(j));
  Particle_List * plist = ProduceParticleList("FixedEnergy",number,E,flav);
  Fill(plist);
  FillHistograms(j,flav);
  //for (std::map<std::string,Detector_Element *>::iterator elit=m_elements.begin();
  //    elit!=m_elements.end();elit++) elit->second->PrintHits();
  plist->Clear();
  delete plist;
}

void Detector::TestReconstructionCodes(const int number,const int j,Flavour flav,
				       const bool background) {
  msg_Out()<<"================================================================"<<std::endl
	   <<"================================================================"<<std::endl
	   <<"================================================================"<<std::endl;
  for (std::map<std::string,Detector_Element *>::iterator elit=m_elements.begin();
       elit!=m_elements.end();elit++) elit->second->Reset();
  double E(Energy(j));
  Particle_List * plist(NULL),* qlist(NULL);
  if (!background) plist = ProduceParticleList("Radiation",number,E,flav);
              else plist = ProduceParticleList("Radiation+",number,E,flav);
  Fill(plist);
  //for (std::map<std::string,Detector_Element *>::iterator elit=m_elements.begin();
  //    elit!=m_elements.end();elit++) elit->second->PrintHits();
  ReconstructObjects(qlist);
  plist->Clear();
  qlist->Clear();
  delete plist;
  delete qlist;
}

double Detector::Energy(const int j) {
  switch (j) {
  case 0:  return 0.4;
  case 1:  return 1.;
  case 2:  return 2.;
  case 3:  return 5.;
  case 4:  return 10.;
  case 5:  return 15.;
  case 6:  return 20.;
  case 7:  return 50.;
  case 8:  return 100.;
  case 9:  return 200.;
  }
  return 500.;
}

void Detector::InitHistograms(const int mode) {
#ifdef USING__ROOT
  if (mode==2) {
    DS_List * segments;
    double etamin,etamax;
    std::string tag(std::string("")),name;
    segments = m_elements["ECal"]->GetDetectorSegments();
    
    for (int j=0;j<10;j++) {
      for (std::set<Detector_Segment *,DS_Order>::iterator ds=segments->begin();
	   ds!=segments->end();ds++) {
	(*ds)->Dimensions(etamin,etamax);
	tag  = std::string("pion");
	name = 
	  std::string("rhad_")+tag+std::string("_")+
	  ATOOLS::ToString(Energy(j))+
	  std::string("_")+ToString(etamin)+
	  std::string("_")+ToString(etamax);
	(*MYROOT::myroot)(new TH2D(name.c_str(),name.c_str(),70,0.,1.4,50,0.,1.0),name);
	name = 
	  std::string("rdep_")+tag+std::string("_")+
	  ATOOLS::ToString(Energy(j))+
	  std::string("_")+ToString(etamin)+
	  std::string("_")+ToString(etamax);
	(*MYROOT::myroot)(new TH1D(name.c_str(),name.c_str(),70,0.,1.4),name);
	tag  = std::string("electron");
	name = 
	  std::string("rdep_")+tag+std::string("_")+
	  ATOOLS::ToString(Energy(j))+
	  std::string("_")+ToString(etamin)+
	  std::string("_")+ToString(etamax);
	(*MYROOT::myroot)(new TH1D(name.c_str(),name.c_str(),70,0.,1.4),name);
	tag  = std::string("photon");
	name = 
	  std::string("rdep_")+tag+std::string("_")+
	  ATOOLS::ToString(Energy(j))+
	  std::string("_")+ToString(etamin)+
	  std::string("_")+ToString(etamax);
	(*MYROOT::myroot)(new TH1D(name.c_str(),name.c_str(),70,0.,1.4),name);
      }
    }
  }
#endif
}

void Detector::FillHistograms(const int j,const Flavour flav) {
#ifdef USING__ROOT
  std::list<Cell *> * cells(m_elements["ECal"]->GetHitCells());

  double E(0.),etamin(0.),etamax(0.);
  std::string tag(std::string("")),name1,name2;

  for (std::list<Cell *>::iterator cit=cells->begin();cit!=cells->end();cit++) {
    Detector_Segment * segment = (*cit)->GetEtastrip()->GetSegment();
    segment->Dimensions(etamin,etamax);
    if (flav.Kfcode()==kf::e)           tag = std::string("electron");
    else if (flav.Kfcode()==kf::photon) tag = std::string("photon");
    else if (flav.Kfcode()==kf::pi)     tag = std::string("pion");
    else return;

    if ((*cit)->ParticleEntries()->size()>1) continue;
    E = (*cit)->ParticleEntries()->begin()->first->Momentum()[0];

    double eta,phi;
    (*cit)->Centroid(eta,phi);
    Cell * hcell(m_elements["HCal"]->GetCell(eta,phi));
    if (hcell->ParticleEntries()->size()>1) continue;
    name1 = 
      std::string("rdep_")+tag+std::string("_")+
      ATOOLS::ToString(Energy(j))+
      std::string("_")+ToString(etamin)+
      std::string("_")+ToString(etamax);
    name2 = 
      std::string("rhad_")+tag+std::string("_")+
      ATOOLS::ToString(Energy(j))+
      std::string("_")+ToString(etamin)+
      std::string("_")+ToString(etamax);
    if (tag==std::string("pion")) {
      if (hcell!=0) {
	((TH2D*)(*MYROOT::myroot)[name2])
	  ->Fill(((*cit)->TotalDeposit()+hcell->TotalDeposit())/E,
		 hcell->TotalDeposit()/((*cit)->TotalDeposit()+hcell->TotalDeposit()),1.);
	((TH1D*)(*MYROOT::myroot)[name1])->Fill(((*cit)->TotalDeposit()+hcell->TotalDeposit())/E,1.);
      }
      else {
	((TH2D*)(*MYROOT::myroot)[name2])->Fill((*cit)->TotalDeposit()/E,0.,1.);
	((TH1D*)(*MYROOT::myroot)[name1])->Fill(((*cit)->TotalDeposit())/E,1.);
      }
    }
    else {
      ((TH1D*)(*MYROOT::myroot)[name1])->Fill(((*cit)->TotalDeposit())/E,1.);
    }  
  }
#endif
}

Particle_List * Detector::ProduceParticleList(const std::string mode,
					      int n,const double E,Flavour flav) {
  if (mode=="Radiation+") n=3*n;
  Vec4D    * moms  = new Vec4D[n+1];
  Flavour  * flavs = new Flavour[n+1];
  int code;
  for (int i=1;i<n+1;i++) {
    if (flav==Flavour(kf::none)) {
      code = int((ran.Get()>0.5?-1:1)*(0.5+3.49999*ran.Get()));
      switch (code) {
      case -3: flavs[i] = Flavour(kf::pi_plus).Bar(); break;
      case -2: flavs[i] = Flavour(kf::mu).Bar();      break;
      case -1: flavs[i] = Flavour(kf::e).Bar();       break;
      case  1: flavs[i] = Flavour(kf::e);             break;
      case  2: flavs[i] = Flavour(kf::mu);            break;
      case  3: flavs[i] = Flavour(kf::pi_plus);       break;
      case  0: 
      default: flavs[i] = Flavour(kf::photon);        break;
      }
    }
    else if (mode=="Radiation+") {
      code = int((ran.Get()>0.5?-1:1)*(0.5+1.49999*ran.Get()));
      switch (code) {
      case -1: flavs[i] = Flavour(kf::pi_plus).Bar(); break;
      case  1: flavs[i] = Flavour(kf::pi_plus);       break;
      case  0: 
      default: flavs[i] = Flavour(kf::pi);            break;
      }
    }
    else flavs[i] = flav;
  }

  if (mode=="Rambo") {
    moms[0] = Vec4D(E,0.,0.,0.);
    PHASIC::Rambo rambo(1,n,flavs,true);
    rambo.GeneratePoint(moms);
  }
  else if (mode=="FixedEnergy") {
    double p,phi,costheta,sintheta;
    for (int i=1;i<n+1;i++) {
      p        = sqrt(sqr(E)-sqr(flavs[i].PSMass()));
      phi      = 2.*M_PI*ran.Get();
      costheta = -1.+2.*ran.Get(),sintheta=sqrt(1.-sqr(costheta));
      moms[i]  = Vec4D(E,p*sintheta*cos(phi),p*sintheta*sin(phi),p*costheta);
    }
  }
  else if ((mode=="Radiation" || mode=="Radiation+") &&
	   (flav==Flavour(kf::e) || flav==Flavour(kf::e).Bar())) {
    int count(1), number;
    if (mode=="Radiation+") {
      moms[0] = Vec4D(n*E,0.,0.,0.);
      PHASIC::Rambo rambo(1,n,flavs,true);
      rambo.GeneratePoint(moms);
      count += 2*n/3;
    }
    double e(E),p,phi,costheta,sintheta;
    double lambda(0.01),epsilon(0.005),omega,z,coszeta,sinzeta,xi;
    double rana,dummy;
    Vec4D  dir;
    do {
      number = int(pow(E/lambda,0.5)*ran.Get());
      if (count+number>n) number = n-count;
      phi      = 2.*M_PI*ran.Get();
      costheta = -1.+2.*ran.Get(),sintheta=sqrt(1.-sqr(costheta));
      dir      = Vec4D(1.,sintheta*cos(phi),sintheta*sin(phi),costheta);
      Poincare rotate(Vec4D(1.,0.,0.,1.),dir);
      for (int i=1;i<number;i++) {
	do {
	  z       = 1.-pow(1.-lambda/E,ran.Get());
	  omega   = z*e;
	} while (omega>e/2. && (e-omega)>2.*lambda);
	e      -= omega;
	do {
	  ran.Gaussian(rana,dummy);
	  coszeta = 1.-epsilon*dabs(rana);
	} while (coszeta>1. || coszeta<-1.);
	sinzeta = sqrt(1.-sqr(coszeta));
	xi      = 2.*M_PI*ran.Get();
	moms[count+i] = omega*Vec4D(1.,sinzeta*cos(xi),sinzeta*sin(xi),coszeta);
	rotate.Rotate(moms[count+i]);
	flavs[count+i] = Flavour(kf::photon);
	//double delta_eta = moms[count+i].Eta()-dir.Eta(),
	//       delta_phi = moms[count+i].Phi()-dir.Phi(),
	//       R = sqrt(sqr(delta_eta)+sqr(delta_phi));
	//if (omega>0.005 && (dabs(delta_eta)<0.075 && dabs(delta_phi)<0.075)) 
	//  std::cout<<"Photon with E = "<<omega<<" in R = "<<R
	//	   <<" ("<<delta_eta<<","<<delta_phi<<")"<<std::endl;
	//std::cout<<"   "<<(count+i)<<" : "<<flavs[count+i]<<" "<<moms[count+i]<<std::endl;
      }
      p = sqrt(sqr(e)-sqr(flav.PSMass()));
      moms[count]  = Vec4D(E,p*sintheta*cos(phi),p*sintheta*sin(phi),p*costheta);
      flavs[count] = flav;
      //std::cout<<"Electron : "<<flavs[count]<<" "<<moms[count]
      //	       <<" --> "<<moms[count].Eta()<<"/"<<moms[count].Phi()<<std::endl;
      count+=number;
    } while (count<n);
  }
  Particle_List * plist = new Particle_List;
  Particle * part(NULL);
  for (int i=1;i<n+1;i++) {
    part = new Particle(i,flavs[i],moms[i],'t');
    plist->push_back(part);
  }
  return plist;
}
