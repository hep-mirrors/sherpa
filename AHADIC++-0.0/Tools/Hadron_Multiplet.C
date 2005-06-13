#include "Hadron_Multiplet.H"
#include "Hadronisation_Parameters.H"
#include "Message.H"

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;



All_Hadron_Multiplets::All_Hadron_Multiplets() :
  p_wavefunctions(NULL), p_multiplets(NULL)
{
  ConstructWaveFunctions();
  ConstructAntiWaveFunctions();
  CreateMultiplets();
  AddMultipletWeights();
}

All_Hadron_Multiplets::~All_Hadron_Multiplets() 
{
  if (p_wavefunctions!=NULL && !p_wavefunctions->empty()) {
    for (Hadron_WF_Miter wfm=p_wavefunctions->begin();wfm!=p_wavefunctions->end();wfm++) {
      if (wfm->second!=NULL) { delete wfm->second; wfm->second=NULL; }
    }
    p_wavefunctions->clear();
    delete p_wavefunctions;
  }
  if (p_multiplets!=NULL && !p_multiplets->empty()) {
    for (Hadron_Multiplet_Miter mp=p_multiplets->begin();mp!=p_multiplets->end();mp++) {
      if (mp->second!=NULL) { delete mp->second; mp->second=NULL; }
    }
    p_multiplets->clear();
    delete p_multiplets;
  }
}

void All_Hadron_Multiplets::ConstructWaveFunctions() 
{
  p_wavefunctions = new Hadron_WF_Map; 
  Flavour hadron;
  int     kfcode;
  Hadron_Wave_Function * wavefunction;

  for (int i=0;i<kf_table.Size();i++) {
    kfcode = kf_table.FromInt(i);
    hadron = Flavour(kf::code(kfcode));
    if (!hadron.IsHadron() || !hadron.IsOn()) continue;
    
    if (int(kfcode/1000)-int(kfcode/10000)*10==0) {
      wavefunction = ConstructMesonWaveFunction(kfcode-10*int(kfcode/10),
						int(kfcode/10000),
						int(kfcode/10)-int(kfcode/100)*10,
						int(kfcode/100)-int(kfcode/1000)*10); 
    }
    else {
      wavefunction = ConstructBaryonWaveFunction(kfcode-10*int(kfcode/10),
						 int(kfcode/10000),
						 int(kfcode/10)-int(kfcode/100)*10,
						 int(kfcode/100)-int(kfcode/1000)*10,
						 int(kfcode/1000)-int(kfcode/10000)*10); 
    }
    if (wavefunction!=NULL) {
      wavefunction->SetFlavour(hadron);
      wavefunction->SetKfCode(kfcode);
      (*p_wavefunctions)[hadron] = wavefunction;
    }
  }
}


Hadron_Wave_Function * All_Hadron_Multiplets::ConstructMesonWaveFunction(int spin,int lp,int fl1,int fl2) 
{
  if (spin==0) return NULL;
  Hadron_Wave_Function * wavefunction=NULL;
  double                 costh, sinth, weight;
  Flavour                flavs[2];
  FlavPair             * pair;

  flavs[0]     = Flavour(kf::code(fl1)).Bar();
  flavs[1]     = Flavour(kf::code(fl2));
  pair         = new FlavPair;
  pair->first  = flavs[1];
  pair->second = flavs[0];
  if (fl1!=fl2 && (flavs[0].Charge()+flavs[1].Charge()!=0.)) {
    wavefunction = new Hadron_Wave_Function;
    wavefunction->AddToWaves(pair,1.);
  }
  else {
    if (fl1==fl2) {
      LookUpAngles(lp,spin,costh,sinth);
      if (fl1==1) {
	wavefunction = new Hadron_Wave_Function;
	wavefunction->AddToWaves(pair,-1./sqrt(2.));
	flavs[0]     = Flavour(kf::code(fl1+1));
	pair         = new FlavPair;
	pair->first  = flavs[0];
	pair->second = flavs[0].Bar();
	wavefunction->AddToWaves(pair,+1./sqrt(2.));
      } 
      if (fl1==2) {
	weight         = costh/sqrt(3.)+sinth/sqrt(6.);
	if (weight>1.e-3) {
	  wavefunction = new Hadron_Wave_Function;
	  wavefunction->AddToWaves(pair,weight);
	  flavs[0]     = Flavour(kf::code(1));
	  pair         = new FlavPair;
	  pair->first  = flavs[0];
	  pair->second = flavs[0].Bar();
	  wavefunction->AddToWaves(pair,weight);
	}
	weight         = costh/sqrt(3.)-2.*sinth/sqrt(6.);
	if (weight>1.e-3) {
	  flavs[0]     = Flavour(kf::code(3));
	  pair         = new FlavPair;
	  pair->first  = flavs[0];
	  pair->second = flavs[0].Bar();
	  if (wavefunction==NULL) wavefunction = new Hadron_Wave_Function;
	  wavefunction->AddToWaves(pair,weight);
	}
      } 
      if (fl1==3) {
	weight         = 2.*costh/sqrt(6.)+sinth/sqrt(3.);
	if (weight>1.e-3) {
	  wavefunction = new Hadron_Wave_Function;
	  wavefunction->AddToWaves(pair,weight);
	}
	weight         = costh/sqrt(6.)-sinth/sqrt(3.);
	if (weight>1.e-3) {
	  flavs[0]     = Flavour(kf::code(1));
	  pair         = new FlavPair;
	  pair->first  = flavs[0];
	  pair->second = flavs[0].Bar();
	  if (wavefunction==NULL) wavefunction = new Hadron_Wave_Function;
	  wavefunction->AddToWaves(pair,weight);
	  flavs[0]     = Flavour(kf::code(2));
	  pair         = new FlavPair;
	  pair->first  = flavs[0];
	  pair->second = flavs[0].Bar();
	  wavefunction->AddToWaves(pair,weight);
	}
      } 
    }
    else if (flavs[0].Charge()+flavs[1].Charge()==0.) {
      wavefunction = new Hadron_Wave_Function;
      wavefunction->AddToWaves(pair,1.);
    }
  }
  wavefunction->SetSpin(spin);
  return wavefunction;
}

Hadron_Wave_Function * All_Hadron_Multiplets::ConstructBaryonWaveFunction(int spin,int lp,
									  int fl1,int fl2,int fl3) 
{
  // Baryon wave functions according to Lichtenberg, Namgung, Wills & Predazzi
  if (spin==0) return NULL;

  int wf=-1;
  int pos1,pos2,pos3,di;

  if (spin==2 && (fl1<fl2 || fl2<fl3)) {
    // Octet
    if (fl1==fl2 && fl1<fl3) {
      if (fl3>3) wf = 810;
            else wf = 81;
      pos1 = fl1; pos2 = fl2; pos3 = fl3; 
    }
    else if (fl1<fl2 && fl2==fl3) {
      wf   = 81;
      pos1 = fl2; pos2 = fl3; pos3 = fl1; 
    }
    else if (fl1<fl2 && fl2<fl3) {
      if (fl3>3) wf = 820;
            else wf = 82;
      pos1 = fl1; pos2 = fl2; pos3 = fl3; 
    }
    else if (fl1>fl2 && fl2<fl3) {
      if (fl3>3) wf = 830;
            else wf = 83;
      pos1 = fl1; pos2 = fl2; pos3 = fl3; 
    }
  }
  else if (spin==4) {
    // Decuplet
    if (fl1==fl2 && fl2==fl3) {
      wf   = 101;
      pos1 = pos2 = pos3 = fl1;
    }
    else if (fl1<fl2 && fl2==fl3) {
      wf   = 102;
      pos1 = fl1; pos2 = fl2; pos3 = fl3;
    }
    else if (fl1==fl2 && fl2<fl3) {
      if (fl3>3) wf = 1020;
            else wf = 102;
      pos1 = fl3; pos2 = fl1; pos3 = fl2;
    }
    else if (fl1<fl2 && fl2<fl3) {
      if (fl3>3) wf = 1030;
            else wf = 103;
      pos1 = fl3; pos2 = fl2; pos3 = fl1;
    }
  }
  //cout<<"Check "<<spin<<" "<<fl1<<" "<<fl2<<" "<<fl3<<" -> "<<wf<<endl;

  Hadron_Wave_Function * wavefunction = new Hadron_Wave_Function;
  FlavPair             * pair;
  switch (wf) {
  case 1020:
    pair         = new FlavPair;
    pair->first  = Flavour(kf::code(pos1));
    pair->second = Flavour(kf::code(pos2*1000+pos3*100+3));
    wavefunction->AddToWaves(pair,1.);
    break;
  case 1030:
    pair         = new FlavPair;
    pair->first  = Flavour(kf::code(pos1));
    pair->second = Flavour(kf::code(pos2*1000+pos3*100+3));
    wavefunction->AddToWaves(pair,1.);
    break;
  case 810:
    pair         = new FlavPair;
    pair->first  = Flavour(kf::code(pos3));
    pair->second = Flavour(kf::code(pos1*1000+pos2*100+3));
    wavefunction->AddToWaves(pair,+1.);
    break;
  case 820:
  case 830:
    if (wf==82) di=+1;
           else di=-1;
    pair         = new FlavPair;
    pair->first  = Flavour(kf::code(pos3));
    pair->second = (pos1>pos2)? Flavour(kf::code(pos1*1000+pos2*100+2+di)) :
                                Flavour(kf::code(pos2*1000+pos1*100+2+di));
    wavefunction->AddToWaves(pair,+1.);
    break;
  case 100:
    pair         = new FlavPair;
    pair->first  = Flavour(kf::code(pos1));
    pair->second = Flavour(kf::code(pos3*1000+pos2*100+3));
    wavefunction->AddToWaves(pair,+1./sqrt(3.));
    pair         = new FlavPair;
    pair->first  = Flavour(kf::code(pos2));
    pair->second = (pos1>pos3)? Flavour(kf::code(pos1*1000+pos3*100+3)) : 
                                Flavour(kf::code(pos3*1000+pos1*100+3)); 
    wavefunction->AddToWaves(pair,+1./sqrt(3.));
    pair         = new FlavPair;
    pair->first  = Flavour(kf::code(pos3));
    pair->second = Flavour(kf::code(pos1*1000+pos2*100+3)); 
    wavefunction->AddToWaves(pair,+1./sqrt(3.));
    break;    
  case 101:
    pair         = new FlavPair;
    pair->first  = Flavour(kf::code(pos1));
    pair->second = Flavour(kf::code(pos3*1000+pos2*100+3));
    wavefunction->AddToWaves(pair,1.);
    break;
  case 102:
    pair         = new FlavPair;
    pair->first  = Flavour(kf::code(pos1));
    pair->second = Flavour(kf::code(pos2*1000+pos3*100+3));
    wavefunction->AddToWaves(pair,1./sqrt(3.));
    pair         = new FlavPair;
    pair->first  = Flavour(kf::code(pos2));
    pair->second = (pos1>pos3) ? Flavour(kf::code(pos1*1000+pos3*100+3)) :
                                 Flavour(kf::code(pos3*1000+pos1*100+3));
    wavefunction->AddToWaves(pair,1./sqrt(3.));
    pair         = new FlavPair;
    pair->first  = Flavour(kf::code(pos3));
    pair->second = (pos1>pos2) ? Flavour(kf::code(pos1*1000+pos2*100+3)) :
                                 Flavour(kf::code(pos2*1000+pos1*100+3));
    wavefunction->AddToWaves(pair,1./sqrt(3.));
    break;
  case 103:
    pair         = new FlavPair;
    pair->first  = Flavour(kf::code(pos1));
    pair->second = Flavour(kf::code(pos2*1000+pos3*100+3));
    wavefunction->AddToWaves(pair,1./sqrt(3.));
    pair         = new FlavPair;
    pair->first  = Flavour(kf::code(pos2));
    pair->second = Flavour(kf::code(pos1*1000+pos3*100+3));
    wavefunction->AddToWaves(pair,1./sqrt(3.));
    pair         = new FlavPair;
    pair->first  = Flavour(kf::code(pos3));
    pair->second = Flavour(kf::code(pos1*1000+pos2*100+3));
    wavefunction->AddToWaves(pair,1./sqrt(3.));
    break;
  case 81:
    pair         = new FlavPair;
    pair->first  = Flavour(kf::code(pos3));
    pair->second = Flavour(kf::code(pos1*1000+pos2*100+3));
    wavefunction->AddToWaves(pair,+1./sqrt(3.));
    pair         = new FlavPair;
    pair->first  = Flavour(kf::code(pos1));
    pair->second = (pos2>pos3)? Flavour(kf::code(pos2*1000+pos3*100+3)) : 
                                Flavour(kf::code(pos3*1000+pos2*100+3)); 
    wavefunction->AddToWaves(pair,+1./sqrt(6.));
    pair         = new FlavPair;
    pair->first  = Flavour(kf::code(pos1));
    pair->second = (pos2>pos3)? Flavour(kf::code(pos2*1000+pos3*100+1)) : 
                                Flavour(kf::code(pos3*1000+pos2*100+1)); 
    wavefunction->AddToWaves(pair,+1./sqrt(2.));
    break;
  case 82:
  case 83:
    if (wf==82) di=+1;
           else di=-1;
    pair         = new FlavPair;
    pair->first  = Flavour(kf::code(pos3));
    pair->second = (pos1>pos2)? Flavour(kf::code(pos1*1000+pos2*100+2+di)) :
                                Flavour(kf::code(pos2*1000+pos1*100+2+di));
    wavefunction->AddToWaves(pair,+1./sqrt(3.));
    pair         = new FlavPair;
    pair->first  = Flavour(kf::code(pos1));
    pair->second = (pos2>pos3)? Flavour(kf::code(pos2*1000+pos3*100+2+di)) : 
                                Flavour(kf::code(pos3*1000+pos2*100+2+di)); 
    wavefunction->AddToWaves(pair,+1./sqrt(12.));
    pair         = new FlavPair;
    pair->first  = Flavour(kf::code(pos2));
    pair->second = (pos1>pos3)? Flavour(kf::code(pos1*1000+pos3*100+2+di)) : 
                                Flavour(kf::code(pos3*1000+pos1*100+2+di)); 
    wavefunction->AddToWaves(pair,+1./sqrt(12.));
    pair         = new FlavPair;
    pair->first  = Flavour(kf::code(pos1));
    pair->second = (pos2>pos3)? Flavour(kf::code(pos2*1000+pos3*100+2-di)) : 
                                Flavour(kf::code(pos3*1000+pos2*100+2-di)); 
    wavefunction->AddToWaves(pair,+1./sqrt(4.));
    pair         = new FlavPair;
    pair->first  = Flavour(kf::code(pos2));
    pair->second = (pos1>pos3)? Flavour(kf::code(pos1*1000+pos3*100+2-di)) : 
                                Flavour(kf::code(pos3*1000+pos1*100+2-di)); 
    wavefunction->AddToWaves(pair,+1./sqrt(4.));
    break;
  }
  wavefunction->SetSpin(spin);
  return wavefunction;
}


void All_Hadron_Multiplets::ConstructAntiWaveFunctions() 
{
  Hadron_Wave_Function * anti; 
  for (Hadron_WF_Miter wfm=p_wavefunctions->begin();wfm!=p_wavefunctions->end();wfm++) {
    anti = wfm->second->Anti();
    //cout<<"Anti wave function for "<<wfm->first<<" : "<<anti<<endl;
    if (anti!=NULL) (*p_wavefunctions)[wfm->first.Bar()] = anti;
  } 
}

void All_Hadron_Multiplets::LookUpAngles(const int angular,const int spin,double & costh,double & sinth)
{
  double angle=0.;
  switch (spin) {
  case 7 : angle = hadpars.Get("Mixing_Angle_3+"); break;
  case 5 : angle = hadpars.Get("Mixing_Angle_2-"); break;
  case 3 : angle = hadpars.Get("Mixing_Angle_1+"); break;
  case 1 : angle = hadpars.Get("Mixing_Angle_0-"); break; 
  default: break;  
  }
  costh = cos(angle); sinth = sin(angle);
}

void All_Hadron_Multiplets::CreateMultiplets()
{
  p_multiplets           = new Hadron_Multiplet_Map;
  int                      kfcode,spin;
  Hadron_Multiplet_Miter   mpletiter;
  Hadron_Multiplet       * multiplet;
  for (Hadron_WF_Miter wfm=p_wavefunctions->begin();wfm!=p_wavefunctions->end();wfm++) {
    kfcode = wfm->second->KfCode();
    spin   = wfm->second->Spin();
    if (spin%2==0) {
      if (abs(kfcode)>4000 && abs(kfcode)<10000) spin += 10;
      if (abs(kfcode)>5000 && abs(kfcode)<10000) spin += 10;
      if (wfm->first.IsAnti()) spin = -spin;
    }
    else {
      if (abs(kfcode)>400 && abs(kfcode)<1000) spin += 10;
      if (abs(kfcode)>500 && abs(kfcode)<1000) spin += 10;
    }
    mpletiter = p_multiplets->find(spin);
    //cout<<"Look for : "<<wfm->first<<" to "<<spin<<" / "<<kfcode<<endl;
    if (mpletiter!=p_multiplets->end()) {
      //msg.Debugging()<<"Add "<<wfm->first<<" to "<<mpletiter->second->m_name<<endl;
      mpletiter->second->AddToElements(wfm->first);
    }
    else {
      multiplet = new Hadron_Multiplet;
      if (dabs(kfcode)<1000) {
	switch (spin) {
	case 1:  multiplet->SetName(string("Pseudoscalars (0-)")); break;
	case 3:  multiplet->SetName(string("Vectors       (1+)")); break;
	case 5:  multiplet->SetName(string("Tensors       (2-)")); break;
	case 7:  multiplet->SetName(string("Tensors       (3+)")); break;
	default: multiplet->SetName(string("Don't know        ")); break;
	}
	multiplet->AddToElements(wfm->first);
	(*p_multiplets)[spin] = multiplet;
      }
      else if (dabs(kfcode)<10000) {
	//cout<<"New multiplet: "<<wfm->first<<" to "<<spin<<endl;
	switch (spin) {
	case   2:  multiplet->SetName(string("Nucleons        (1/2)")); break;
	case   4:  multiplet->SetName(string("Decuplet        (3/2)")); break;
	case  12:  multiplet->SetName(string("Nucleons_c      (1/2)")); break;
	case  14:  multiplet->SetName(string("Decuplet_c      (3/2)")); break;
	case  22:  multiplet->SetName(string("Nucleons_b      (1/2)")); break;
	case  24:  multiplet->SetName(string("Decuplet_b      (3/2)")); break;
	case  -2:  multiplet->SetName(string("Anti-Nucleons   (1/2)")); break;
	case  -4:  multiplet->SetName(string("Anti-Decuplet   (3/2)")); break;
	case -12:  multiplet->SetName(string("Anti-Nucleons_c (1/2)")); break;
	case -14:  multiplet->SetName(string("Anti-Decuplet_c (3/2)")); break;
	case -22:  multiplet->SetName(string("Anti-Nucleons_c (1/2)")); break;
	case -24:  multiplet->SetName(string("Anti-Decuplet_c (3/2)")); break;
	default:   multiplet->SetName(string("Don't know     (anti)")); break;
	}
	multiplet->AddToElements(wfm->first);
	(*p_multiplets)[spin] = multiplet;
      }
    }
  }
} 

void All_Hadron_Multiplets::AddMultipletWeights()
{
  Hadron_Multiplet     * mplet;
  Hadron_Wave_Function * wave;
  int spin;
  for (Hadron_Multiplet_Miter miter=p_multiplets->begin();
       miter!=p_multiplets->end();miter++) {
    spin  = abs(miter->first);
    if (spin>20) spin-=20;
    if (spin>10) spin-=10;
    mplet = miter->second;
    miter->second->SetSpinWeight(double(spin));
    miter->second->SetExtraWeight(1.);
    miter->second->SetWeight();
    for (FlSetIter flit=miter->second->GetElements()->begin();
	 flit!=miter->second->GetElements()->end();flit++) {
      wave = GetWaveFunction((*flit));
      if (wave!=NULL) wave->SetMultipletWeight(miter->second->Weight());
      else {
	msg.Error()<<"ERROR in All_Hadron_Multiplets::AddMultipletWeights():"<<endl
		   <<"   No wave function found for "<<(*flit)
		   <<",continue and hope for the best."<<endl;
      }
    }
  }
}

Hadron_Wave_Function * All_Hadron_Multiplets::GetWaveFunction(ATOOLS::Flavour flav) 
{
  Hadron_WF_Miter wfm=p_wavefunctions->find(flav);
  if (wfm!=p_wavefunctions->end()) return wfm->second;
  return NULL;
}

void All_Hadron_Multiplets::PrintWaveFunctions() 
{
  Hadron_WF_Miter        wfm;
  Hadron_Wave_Function * wf;
  map<Flavour,double> checkit;
  for (Hadron_Multiplet_Miter mplet=p_multiplets->begin();mplet!=p_multiplets->end();mplet++) {
    msg.Out()<<"-----------------------------------------------"<<endl
	     <<" "<<mplet->second->Name()<<" with "<<mplet->second->Size()<<" elements: "<<endl;
    for (FlSetIter fl=mplet->second->GetElements()->begin();
	 fl!=mplet->second->GetElements()->end();fl++) {
      wfm = p_wavefunctions->find((*fl));
      if (wfm!=p_wavefunctions->end()) { wf = wfm->second; msg.Out()<<(*wf); }
      
      if (mplet->second->Name()==string("Nucleons        (1/2)")) {
       	WFcomponent * wfc = wf->GetWaves();
       	for (WFcompiter cit = wfc->begin();cit!=wfc->end();cit++) {
       	  if (checkit.find(cit->first->first)==checkit.end()) 
       	    checkit[cit->first->first] = sqr(cit->second);
       	  else checkit[cit->first->first] += sqr(cit->second);
	  if (checkit.find(cit->first->second)==checkit.end()) 
       	    checkit[cit->first->second] = sqr(cit->second);
	  else checkit[cit->first->second] += sqr(cit->second);
       	}
      }
    } 
    if (mplet->second->Name()==string("Nucleons        (1/2)")) {
      msg.Out()<<"Weight of individual components in "<<mplet->second->Name()<<" : "<<endl;
      for (map<Flavour,double>::iterator it=checkit.begin();it!=checkit.end();it++) 
       	msg.Out()<<"     -> "<<it->first<<" : "<<it->second<<endl;
    }
    msg.Out()<<"-----------------------------------------------"<<endl;
  }
  msg.Out()<<endl
	   <<"-------- END OF ALL_HADRON_MULTIPLETS ------"<<endl;  
}

