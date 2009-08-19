#include "AHADIC++/Tools/Hadron_Multiplet.H"
#include "AHADIC++/Tools/Hadronisation_Parameters.H"
#include "ATOOLS/Org/Message.H"

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;



All_Hadron_Multiplets::All_Hadron_Multiplets() :
  m_singletsuppression(hadpars.Get("Singlet_Suppression")),
  p_wavefunctions(NULL), p_multiplets(NULL)
{
  ConstructWaveFunctions();
  ConstructAntiWaveFunctions();
  CreateMultiplets();
  AddMultipletWeights();
  //PrintMultiplets();
  //PrintWaveFunctions();
  //abort();
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
  int     fl1,fl2,fl3, help;
  Hadron_Wave_Function * wavefunction;

  for(KFCode_ParticleInfo_Map::const_iterator kfit(s_kftable.begin());
      kfit!=s_kftable.end();++kfit) {
    kf_code kfcode = kfit->first;
    hadron = Flavour(kfcode);
    if (!hadron.IsHadron() || !hadron.IsOn()) continue;
    fl1 = int(kfcode/10)-int(kfcode/100)*10;
    fl2 = int(kfcode/100)-int(kfcode/1000)*10;
    fl3 = int(kfcode/1000)-int(kfcode/10000)*10;

    //std::cout<<"Construct wave function for "<<hadron<<" --> "
    //	     <<fl1<<" "<<fl2<<" "<<fl3<<" spin = "<<kfcode-10*int(kfcode/10)<<"."<<std::endl;
    if (fl3==0) {
      if (int(fl2/2.)!=fl2/2.) { help = fl1; fl1 = fl2; fl2 = help; }
      wavefunction = ConstructMesonWaveFunction(int(kfcode/9000000),
						int(kfcode/100000),
						int(kfcode/10000),
						kfcode-10*int(kfcode/10),
						fl1,fl2); 
    }
    else {
      wavefunction = ConstructBaryonWaveFunction(int(kfcode/10000),
						 kfcode-10*int(kfcode/10),
						 fl1,fl2,fl3); 
    }
    if (wavefunction!=NULL) {
      wavefunction->SetFlavour(hadron);
      wavefunction->SetKfCode(kfcode);
      (*p_wavefunctions)[hadron] = wavefunction;
    }
  }
}


Hadron_Wave_Function * 
All_Hadron_Multiplets::ConstructMesonWaveFunction(const int iso0,const int rp,const int lp,
						  const int spin,const int fl1,const int fl2) 
{
  if (spin==0) return NULL;
  Hadron_Wave_Function * wavefunction=NULL;
  double                 costh, sinth, weight;
  Flavour                flavs[2];
  Flavour_Pair         * pair;

  flavs[0]     = Flavour((kf_code)(fl1)).Bar();
  flavs[1]     = Flavour((kf_code)(fl2));
  pair         = new Flavour_Pair;
  pair->first  = flavs[1];
  pair->second = flavs[0];
  if (fl1!=fl2 && (flavs[0].Charge()+flavs[1].Charge()!=0.)) {
    wavefunction = new Hadron_Wave_Function;
    wavefunction->AddToWaves(pair,1.);
  }
  else if (fl1==fl2) {
    if (fl1==1) {
      LookUpAngles(lp,spin,costh,sinth);
      wavefunction = new Hadron_Wave_Function;
      wavefunction->AddToWaves(pair,-1./sqrt(2.));
      flavs[0]     = Flavour((kf_code)(2));
      pair         = new Flavour_Pair;
      pair->first  = flavs[0];
      pair->second = flavs[0].Bar();
      wavefunction->AddToWaves(pair,+1./sqrt(2.));
    } 
    else if (fl1==2 && iso0==1) {
      weight       = 1/sqrt(3.);
      wavefunction = new Hadron_Wave_Function;
      wavefunction->AddToWaves(pair,weight);
      flavs[0]     = Flavour((kf_code)(1));
      pair         = new Flavour_Pair;
      pair->first  = flavs[0];
      pair->second = flavs[0].Bar();
      wavefunction->AddToWaves(pair,weight);
      flavs[0]     = Flavour((kf_code)(3));
      pair         = new Flavour_Pair;
      pair->first  = flavs[0];
      pair->second = flavs[0].Bar();
      wavefunction->AddToWaves(pair,weight);
    }
    else if ((fl1==2 && spin!=1) || (fl1==3 && spin==1)) {
      weight         = (sinth/sqrt(6.)+costh/sqrt(3.))*sqrt(1.-sqr(m_singletsuppression));
      if (dabs(weight)>1.e-3) {
	wavefunction = new Hadron_Wave_Function;
	wavefunction->AddToWaves(pair,weight);
	flavs[0]     = Flavour((kf_code)(1));
	pair         = new Flavour_Pair;
	pair->first  = flavs[0];
	pair->second = flavs[0].Bar();
	wavefunction->AddToWaves(pair,weight);
      }
      weight         = (-2.*sinth/sqrt(6.)+costh/sqrt(3.))*sqrt(1.-sqr(m_singletsuppression));
      if (dabs(weight)>1.e-3) {
	flavs[0]     = Flavour((kf_code)(3));
	pair         = new Flavour_Pair;
	pair->first  = flavs[0];
	pair->second = flavs[0].Bar();
	if (wavefunction==NULL) wavefunction = new Hadron_Wave_Function;
	wavefunction->AddToWaves(pair,weight);
      }
    }
    else if ((fl1==2 && spin==1) || (fl1==3 && spin!=1)) {
      weight         = (-2.*costh/sqrt(6.)-sinth/sqrt(3.))*m_singletsuppression;
      if (dabs(weight)>1.e-3) {
	wavefunction = new Hadron_Wave_Function;
	wavefunction->AddToWaves(pair,weight);
      }
      weight         = (costh/sqrt(6.)-sinth/sqrt(3.))*m_singletsuppression;
      if (dabs(weight)>1.e-3) {
	flavs[0]     = Flavour((kf_code)(1));
	pair         = new Flavour_Pair;
	pair->first  = flavs[0];
	pair->second = flavs[0].Bar();
	if (wavefunction==NULL) wavefunction = new Hadron_Wave_Function;
	wavefunction->AddToWaves(pair,weight);
	flavs[0]     = Flavour((kf_code)(2));
	pair         = new Flavour_Pair;
	pair->first  = flavs[0];
	pair->second = flavs[0].Bar();
	wavefunction->AddToWaves(pair,weight);
      }
    }
    else if (fl1==fl2 && (fl1==4 || fl1==5)) {
      wavefunction = new Hadron_Wave_Function;
      wavefunction->AddToWaves(pair,1.);
    }
  }
  else if (flavs[0].Charge()+flavs[1].Charge()==0.) {
    wavefunction = new Hadron_Wave_Function;
    wavefunction->AddToWaves(pair,1.);
  }
  if (wavefunction) wavefunction->SetSpin(spin);
  return wavefunction;
}

Hadron_Wave_Function * 
All_Hadron_Multiplets::ConstructBaryonWaveFunction(int lp,int spin,
						   int fl1,int fl2,int fl3) 
{
  // Baryon wave functions according to Lichtenberg, Namgung, Wills & Predazzi
  if (spin==0) return NULL;

  int wf(-1),pos1(-1),pos2(-1),pos3(-1),di(-1);

  if ((spin==2 || lp!=0) && (fl1<fl2 || fl2<fl3)) {
    // Octet
    // proton- and neutron-like  
    if (fl1==fl2 && fl1<fl3) {
      if (fl3>3) wf = 810;
            else wf = 81;
      pos1 = fl1; pos2 = fl2; pos3 = fl3; 
    }
    // proton- and neutron-like  
    else if (fl1<fl2 && fl2==fl3) {
      wf   = 81;
      pos1 = fl2; pos2 = fl3; pos3 = fl1; 
    }
    // sigma0-like  
    else if (fl1<fl2 && fl2<fl3) {
      if (fl3>3) wf = 820;
            else wf = 82;
      pos1 = fl1; pos2 = fl2; pos3 = fl3; 
    }
    // lambda-like  
    else if (fl1>fl2 && fl2<fl3) {
      if (fl3>3) wf = 830;
            else wf = 83;
      pos1 = fl1; pos2 = fl2; pos3 = fl3; 
    }
  }
  else if (spin==4) {
    // Decuplet
    // delta++-like
    if (fl1==fl2 && fl2==fl3) {
      wf   = 101;
      pos1 = pos2 = pos3 = fl1;
    }
    // delta+/--like
    else if (fl1<fl2 && fl2==fl3) {
      wf   = 102;
      pos1 = fl1; pos2 = fl2; pos3 = fl3;
    }
    // delta+/--like
    else if (fl1==fl2 && fl2<fl3) {
      if (fl3>3) wf = 1020;
            else wf = 102;
      pos1 = fl3; pos2 = fl1; pos3 = fl2;
    }
    // lambda*0-like  
    else if (fl1<fl2 && fl2<fl3) {
      if (fl3>3) wf = 1030;
            else wf = 103;
      pos1 = fl3; pos2 = fl2; pos3 = fl1;
    }
  }

  if (pos1<0 || pos2<0 || pos3<0 || wf<0) {
    msg_Tracking()<<"Warning in "<<METHOD<<"("<<lp<<","<<spin<<","<<fl1<<","<<fl2<<","<<fl3<<"):"
		  <<std::endl
		  <<"   leads to unknown combination of wf type and positions."<<std::endl
		  <<"   wf = "<<wf<<", pos123 = "<<pos1<<" "<<pos2<<" "<<pos3<<";"
		  <<" will ignore this hadron."<<std::endl;
    return NULL;
  }

  Hadron_Wave_Function * wavefunction = new Hadron_Wave_Function;
  Flavour_Pair         * pair;

  switch (wf) {
  case 1020:
    pair         = new Flavour_Pair;
    pair->first  = Flavour((kf_code)(pos1));
    pair->second = Flavour((kf_code)(pos2*1000+pos3*100+3));
    wavefunction->AddToWaves(pair,1.);
    break;
  case 1030:
    pair         = new Flavour_Pair;
    pair->first  = Flavour((kf_code)(pos1));
    pair->second = Flavour((kf_code)(pos2*1000+pos3*100+3));
    wavefunction->AddToWaves(pair,1.);
    break;
  case 810:
    pair         = new Flavour_Pair;
    pair->first  = Flavour((kf_code)(pos3));
    pair->second = Flavour((kf_code)(pos1*1000+pos2*100+3));
    wavefunction->AddToWaves(pair,+1.);
    break;
  case 820:
  case 830:
    if (wf==82 || wf==820) di=+1;
                      else di=-1;
    pair         = new Flavour_Pair;
    pair->first  = Flavour((kf_code)(pos3));
    pair->second = (pos1>pos2)? Flavour((kf_code)(pos1*1000+pos2*100+2+di)) :
                                Flavour((kf_code)(pos2*1000+pos1*100+2+di));
    wavefunction->AddToWaves(pair,+1.);
    break;
  case 100:
    pair         = new Flavour_Pair;
    pair->first  = Flavour((kf_code)(pos1));
    pair->second = Flavour((kf_code)(pos3*1000+pos2*100+3));
    wavefunction->AddToWaves(pair,+1./sqrt(3.));
    pair         = new Flavour_Pair;
    pair->first  = Flavour((kf_code)(pos2));
    pair->second = (pos1>pos3)? Flavour((kf_code)(pos1*1000+pos3*100+3)) : 
                                Flavour((kf_code)(pos3*1000+pos1*100+3)); 
    wavefunction->AddToWaves(pair,+1./sqrt(3.));
    pair         = new Flavour_Pair;
    pair->first  = Flavour((kf_code)(pos3));
    pair->second = Flavour((kf_code)(pos1*1000+pos2*100+3)); 
    wavefunction->AddToWaves(pair,+1./sqrt(3.));
    break;    
  case 101:
    pair         = new Flavour_Pair;
    pair->first  = Flavour((kf_code)(pos1));
    pair->second = Flavour((kf_code)(pos3*1000+pos2*100+3));
    wavefunction->AddToWaves(pair,1.);
    break;
  case 102:
    pair         = new Flavour_Pair;
    pair->first  = Flavour((kf_code)(pos1));
    pair->second = Flavour((kf_code)(pos2*1000+pos3*100+3));
    wavefunction->AddToWaves(pair,1./sqrt(3.));
    pair         = new Flavour_Pair;
    pair->first  = Flavour((kf_code)(pos2));
    pair->second = (pos1>pos3) ? Flavour((kf_code)(pos1*1000+pos3*100+3)) :
                                 Flavour((kf_code)(pos3*1000+pos1*100+3));
    wavefunction->AddToWaves(pair,sqrt(2./3.));
    break;
  case 103:
    pair         = new Flavour_Pair;
    pair->first  = Flavour((kf_code)(pos1));
    pair->second = Flavour((kf_code)(pos2*1000+pos3*100+3));
    wavefunction->AddToWaves(pair,1./sqrt(3.));
    pair         = new Flavour_Pair;
    pair->first  = Flavour((kf_code)(pos2));
    pair->second = Flavour((kf_code)(pos1*1000+pos3*100+3));
    wavefunction->AddToWaves(pair,1./sqrt(3.));
    pair         = new Flavour_Pair;
    pair->first  = Flavour((kf_code)(pos3));
    pair->second = Flavour((kf_code)(pos1*1000+pos2*100+3));
    wavefunction->AddToWaves(pair,1./sqrt(3.));
    break;
  case 81:
    pair         = new Flavour_Pair;
    pair->first  = Flavour((kf_code)(pos3));
    pair->second = Flavour((kf_code)(pos1*1000+pos2*100+3));
    wavefunction->AddToWaves(pair,+1./sqrt(3.));
    pair         = new Flavour_Pair;
    pair->first  = Flavour((kf_code)(pos1));
    pair->second = (pos2>pos3)? Flavour((kf_code)(pos2*1000+pos3*100+3)) : 
                                Flavour((kf_code)(pos3*1000+pos2*100+3)); 
    wavefunction->AddToWaves(pair,+1./sqrt(6.));
    pair         = new Flavour_Pair;
    pair->first  = Flavour((kf_code)(pos1));
    pair->second = (pos2>pos3)? Flavour((kf_code)(pos2*1000+pos3*100+1)) : 
                                Flavour((kf_code)(pos3*1000+pos2*100+1)); 
    wavefunction->AddToWaves(pair,+1./sqrt(2.));
    break;
  case 82:
  case 83:
    if (wf==82) di=+1;
           else di=-1;
    pair         = new Flavour_Pair;
    pair->first  = Flavour((kf_code)(pos3));
    pair->second = (pos1>pos2)? Flavour((kf_code)(pos1*1000+pos2*100+2+di)) :
                                Flavour((kf_code)(pos2*1000+pos1*100+2+di));
    wavefunction->AddToWaves(pair,+1./sqrt(3.));
    pair         = new Flavour_Pair;
    pair->first  = Flavour((kf_code)(pos1));
    pair->second = (pos2>pos3)? Flavour((kf_code)(pos2*1000+pos3*100+2+di)) : 
                                Flavour((kf_code)(pos3*1000+pos2*100+2+di)); 
    wavefunction->AddToWaves(pair,+1./sqrt(12.));
    pair         = new Flavour_Pair;
    pair->first  = Flavour((kf_code)(pos2));
    pair->second = (pos1>pos3)? Flavour((kf_code)(pos1*1000+pos3*100+2+di)) : 
                                Flavour((kf_code)(pos3*1000+pos1*100+2+di)); 
    wavefunction->AddToWaves(pair,+1./sqrt(12.));
    //pair         = new Flavour_Pair;
    //pair->first  = Flavour((kf_code)(pos3));
    //pair->second = (pos1>pos2)? Flavour((kf_code)(pos1*1000+pos2*100+2-di)) :
    //                            Flavour((kf_code)(pos2*1000+pos1*100+2-di));
    //wavefunction->AddToWaves(pair,+1./sqrt(4.));
    pair         = new Flavour_Pair;
    pair->first  = Flavour((kf_code)(pos1));
    pair->second = (pos2>pos3)? Flavour((kf_code)(pos2*1000+pos3*100+2-di)) : 
                                Flavour((kf_code)(pos3*1000+pos2*100+2-di)); 
    wavefunction->AddToWaves(pair,+1./sqrt(4.));
    pair         = new Flavour_Pair;
    pair->first  = Flavour((kf_code)(pos2));
    pair->second = (pos1>pos3)? Flavour((kf_code)(pos1*1000+pos3*100+2-di)) : 
                                Flavour((kf_code)(pos3*1000+pos1*100+2-di)); 
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
    if (anti!=NULL) (*p_wavefunctions)[wfm->first.Bar()] = anti;
  } 
}

void All_Hadron_Multiplets::LookUpAngles(const int angular,const int spin,double & costh,double & sinth)
{
  double angle=0.;
  switch (spin) {
  case 9 : angle = hadpars.Get("Mixing_Angle_4+"); break;
  case 7 : angle = hadpars.Get("Mixing_Angle_3-"); break;
  case 5 : angle = hadpars.Get("Mixing_Angle_2+"); break;
  case 3 : angle = hadpars.Get("Mixing_Angle_1-"); break;
  case 1 : angle = hadpars.Get("Mixing_Angle_0+"); break; 
  default: break;  
  }
  costh = cos(angle); sinth = sin(angle);
}

void All_Hadron_Multiplets::CreateMultiplets()
{
  p_multiplets           = new Hadron_Multiplet_Map;
  int                      kfcode,spin,orbital,radial,iso0;
  Hadron_Multiplet_Miter   mpletiter;
  Hadron_Multiplet       * multiplet;
  std::string              prefix;
  for (Hadron_WF_Miter wfm=p_wavefunctions->begin();wfm!=p_wavefunctions->end();wfm++) {
    kfcode  = abs(wfm->second->KfCode());
    spin    = wfm->second->Spin();
    iso0    = int(kfcode/9000000);
    kfcode -= iso0*9000000;
    radial  = int(kfcode/100000);
    kfcode -= radial*100000;
    orbital = int(kfcode/10000);
    kfcode -= orbital*10000;
    spin   += iso0*10000+radial*1000+orbital*100;
    if (spin%2==0) { if (wfm->first.IsAnti()) spin = -spin; }
    Flavour flav = wfm->first;
    mpletiter = p_multiplets->find(spin);
    if (mpletiter!=p_multiplets->end()) mpletiter->second->AddToElements(flav);
    else {
      multiplet = new Hadron_Multiplet;
      if (kfcode<1000) {
	switch (spin) {
	case     1:  
	  multiplet->SetName(string("Pseudoscalars         (0-+)")); 
	  multiplet->SetExtraWeight(hadpars.Get("Multiplet_Meson_L0R0S0"));
	  break;
	case     3:  
	  multiplet->SetName(string("Vectors               (1--)"));
	  multiplet->SetExtraWeight(hadpars.Get("Multiplet_Meson_L0R0S1"));
	  break;
	case     5:  
	  multiplet->SetName(string("Tensors               (2++)"));
	  multiplet->SetExtraWeight(hadpars.Get("Multiplet_Meson_L0R0S2"));
	  break;
	case     7:  
	  multiplet->SetName(string("Tensors               (3--)"));
	  multiplet->SetExtraWeight(hadpars.Get("Multiplet_Meson_L0R0S3"));
	  break;
	case     9:  
	  multiplet->SetName(string("Tensors               (4++)"));
	  multiplet->SetExtraWeight(hadpars.Get("Multiplet_Meson_L0R0S4"));
	  break;
	case   101:  
	  multiplet->SetName(string("L=1 Pseudoscalars     (0++)"));
	  multiplet->SetExtraWeight(hadpars.Get("Multiplet_Meson_L1R0S0"));
	  break;
	case   103:  
	  multiplet->SetName(string("L=1 Vectors           (1+-)"));
	  multiplet->SetExtraWeight(hadpars.Get("Multiplet_Meson_L1R0S1"));
	  break;
	case   105:  
	  multiplet->SetName(string("L=1 Tensors           (2-+)"));
	  multiplet->SetExtraWeight(hadpars.Get("Multiplet_Meson_L1R0S2"));
	  break;
	case   107:  
	  multiplet->SetName(string("L=1 Tensors           (3+-)"));
	  multiplet->SetExtraWeight(hadpars.Get("Multiplet_Meson_L1R0S3"));
	  break;
	case   109:  
	  multiplet->SetName(string("L=1 Tensors           (4-+)"));
	  multiplet->SetExtraWeight(hadpars.Get("Multiplet_Meson_L1R0S4"));
	  break;
	case   203:  
	  multiplet->SetName(string("L=2 Vectors           (1++)"));
	  multiplet->SetExtraWeight(hadpars.Get("Multiplet_Meson_L2R0S1"));
	  break;
	case   205:  
	  multiplet->SetName(string("L=2 Tensors           (2--)"));
	  multiplet->SetExtraWeight(hadpars.Get("Multiplet_Meson_L2R0S2"));
	  break;
	case   207:  
	  multiplet->SetName(string("L=2 Tensors           (3++)"));
	  multiplet->SetExtraWeight(hadpars.Get("Multiplet_Meson_L2R0S3"));
	  break;
	case   209:  
	  multiplet->SetName(string("L=2 Tensors           (4--)"));
	  multiplet->SetExtraWeight(hadpars.Get("Multiplet_Meson_L2R0S4"));
	  break;
	case   303:  
	  multiplet->SetName(string("L=3 Vectors           (1--)"));
	  multiplet->SetExtraWeight(hadpars.Get("Multiplet_Meson_L3R0S1"));
	  break;
	case   305:  
	  multiplet->SetName(string("L=3 Tensors           (2++)"));
	  multiplet->SetExtraWeight(hadpars.Get("Multiplet_Meson_L3R0S2"));
	  break;
	case   307:  
	  multiplet->SetName(string("L=3 Tensors           (3--)"));
	  multiplet->SetExtraWeight(hadpars.Get("Multiplet_Meson_L3R0S3"));
	  break;
	case   309:  
	  multiplet->SetName(string("L=3 Tensors           (4++)"));
	  multiplet->SetExtraWeight(hadpars.Get("Multiplet_Meson_L3R0S4"));
	  break;
	case  1001:  
	  multiplet->SetName(string("R=1 Pseudoscalars     (0-+)"));
	  multiplet->SetExtraWeight(hadpars.Get("Multiplet_Meson_L0R1S0"));
	  break;
	case  1003:  
	  multiplet->SetName(string("R=1 Vectors           (1--)"));
	  multiplet->SetExtraWeight(hadpars.Get("Multiplet_Meson_L0R1S1"));
	  break;
	case 10001:  
	  multiplet->SetName(string("Isoscalar Pseudoscalars     (0-+)"));
	  multiplet->SetExtraWeight(hadpars.Get("Multiplet_Meson_Iso_L0R0S0"));
	  break;
	case 10003:  
	  multiplet->SetName(string("Isoscalar Vectors           (1--)"));
	  multiplet->SetExtraWeight(hadpars.Get("Multiplet_Meson_Iso_L0R0S1"));
	  break;
	case 10101:  
	  multiplet->SetName(string("Isoscalar L=1 Pseudoscalars (0-+)"));
	  multiplet->SetExtraWeight(hadpars.Get("Multiplet_Meson_Iso_L1R0S0"));
	  break;
	case 10103:  
	  multiplet->SetName(string("Isoscalar L=1 Vectors       (1--)"));
	  multiplet->SetExtraWeight(hadpars.Get("Multiplet_Meson_Iso_L1R0S1"));
	  break;
	default: multiplet->SetName(string("Don't know        ")); break;
	}
      }
      else if (kfcode<10000) {
	switch (spin) {
	case    2:  
	  multiplet->SetName(string("Nucleons                  (1/2)"));
	  multiplet->SetExtraWeight(hadpars.Get("Multiplet_Nucleon_L0R0S1/2"));
	  break;
	case   -2:  
	  multiplet->SetName(string("Anti-Nucleons             (1/2)"));
	  multiplet->SetExtraWeight(hadpars.Get("Multiplet_Nucleon_L0R0S1/2"));
	  break;
	case  102:                   
	  multiplet->SetName(string("L=0 excited Nucleons      (1/2)"));
	  multiplet->SetExtraWeight(hadpars.Get("Multiplet_exc_Nucleon_L0R0S1/2"));
	  break;
	case -102: 
	  multiplet->SetName(string("L=0 excited Anti-Nucleons  (1/2)"));
	  multiplet->SetExtraWeight(hadpars.Get("Multiplet_exc_Nucleon_L0R0S1/2"));
	  break;
	case  202:                   
	  multiplet->SetName(string("L=1 excited Nucleons      (1/2)"));
	  multiplet->SetExtraWeight(hadpars.Get("Multiplet_exc_Nucleon_L1R0S1/2"));
	  break;
	case -202: 
	  multiplet->SetName(string("L=1 excited Anti-Nucleons  (1/2)"));
	  multiplet->SetExtraWeight(hadpars.Get("Multiplet_exc_Nucleon_L1R0S1/2"));
	  break;
	case  104:                   
	  multiplet->SetName(string("L=1 excited Nucleons      (3/2)"));
	  multiplet->SetExtraWeight(hadpars.Get("Multiplet_exc_Nucleon_L1R0S3/2"));
	  break;
	case -104: 
	  multiplet->SetName(string("L=1 excited Anti-Nucleons  (3/2)"));
	  multiplet->SetExtraWeight(hadpars.Get("Multiplet_exc_Nucleon_L1R0S3/2"));
	  break;
	case    4:  
	  multiplet->SetName(string("Decuplet                  (3/2)"));
	  multiplet->SetExtraWeight(hadpars.Get("Multiplet_Delta_L0R0S3/2"));
	  break;
	case   -4:  
	  multiplet->SetName(string("Anti-Decuplet             (3/2)"));
	  multiplet->SetExtraWeight(hadpars.Get("Multiplet_Delta_L0R0S3/2"));
	  break;
	case  10104:
	  multiplet->SetName(string("L=1 excited Decuplet      (3/2)"));
	  multiplet->SetExtraWeight(hadpars.Get("Multiplet_exc_Delta_L1R0S3/2"));
	  break;
	case -10104:
	  multiplet->SetName(string("L=1 excited anti-Decuplet (3/2)"));
	  multiplet->SetExtraWeight(hadpars.Get("Multiplet_exc_Delta_L1R0S3/2"));
	  break;
	default:   multiplet->SetName(string("Don't know                (anti)")); break;
	}
      }
      multiplet->GetElements()->insert(flav);
      (*p_multiplets)[spin] = multiplet;
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
    spin = abs(miter->first);
    spin -= int(spin/1000)*1000;
    spin -= int(spin/100)*100;
    spin -= int(spin/10)*10;
    mplet = miter->second;
    miter->second->SetSpinWeight(double(spin));
    miter->second->SetWeight();
    for (FlSetIter flit=miter->second->GetElements()->begin();
 	 flit!=miter->second->GetElements()->end();flit++) {
      wave = GetWaveFunction((*flit));
      if (wave!=NULL) {
	wave->SetMultipletWeight(miter->second->Weight());
      }
      else {
	if (wave==NULL) {
	  msg_Error()<<"ERROR in All_Hadron_Multiplets::AddMultipletWeights():"<<endl
		     <<"   No wave function found for "<<(*flit)
		     <<",continue and hope for the best."<<endl;
	}
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
  map<Flavour,double> checkit;
  for (Hadron_Multiplet_Miter mplet=p_multiplets->begin();mplet!=p_multiplets->end();mplet++) {
    if (mplet->second->Weight()<=0.) continue;
    msg_Out()<<"-----------------------------------------------"<<endl
	     <<" "<<mplet->second->Name()<<" with "<<mplet->second->Size()<<" elements: "<<endl;
    for (FlSetIter fl=mplet->second->GetElements()->begin();
	 fl!=mplet->second->GetElements()->end();fl++) {
      wfm = p_wavefunctions->find((*fl));
      if (wfm!=p_wavefunctions->end()) msg_Out()<<(*wfm->second); 
      if (mplet->second->Name()==string("Nucleons        (1/2)")) {
       	WFcomponent * wfc = wfm->second->GetWaves();
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
      msg_Out()<<"Weight of individual components in "<<mplet->second->Name()<<" : "<<endl;
      for (map<Flavour,double>::iterator it=checkit.begin();it!=checkit.end();it++) 
       	msg_Out()<<"     -> "<<it->first<<" : "<<it->second<<endl;
    }
    msg_Out()<<"-----------------------------------------------"<<endl;
  }
  msg_Out()<<endl
	   <<"-------- END OF ALL_HADRON_MULTIPLETS ------"<<endl;  
}


void All_Hadron_Multiplets::PrintMultiplets() 
{
  for (Hadron_Multiplet_Miter miter=p_multiplets->begin();
       miter!=p_multiplets->end();miter++) {
    msg_Out()<<"* "<<miter->first<<" "<<miter->second->Name()<<" : "
	     <<"spin weight = "<<miter->second->Weight()<<", "
	     <<"extra weight = "<<miter->second->ExtraWeight()<<endl;
    for (FlSetIter flit=miter->second->GetElements()->begin();
 	 flit!=miter->second->GetElements()->end();flit++) {
      msg_Out()<<"  "<<(*flit);
    }
    msg_Out()<<endl<<endl;
  }
}
