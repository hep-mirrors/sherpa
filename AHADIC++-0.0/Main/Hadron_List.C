#include "Hadron_List.H"
#include "Flavour.H"
#include "MathTools.H"
#include "Data_Read.H"
#include "Message.H"

using namespace AHADIC;
using namespace ATOOLS;

Hadron_Parameters AHADIC::hadpars;

//##############################################################################
//##############################################################################
//##############################################################################
//##############################################################################
//##############################################################################


ConstituentMasses::ConstituentMasses(bool no_diquarks) {
  ConstituentCharacteristic * cc;
  double flwt=0., spwt=0., sm = hadpars.Get("AngularSmearing");

  // Gluon
  cc = new ConstituentCharacteristic(hadpars.Get("Mass_glue"),2,flwt,spwt);
  CCMap[Flavour(kf::gluon)] = cc;

  // Light quarks
  spwt = 2.;
  flwt = (1.-hadpars.Get("Baryon_supression"));
  cc = new ConstituentCharacteristic(hadpars.Get("Mass_down"),1,flwt,spwt,sm);
  CCMap[Flavour(kf::d)] = cc;
  cc = new ConstituentCharacteristic(hadpars.Get("Mass_up"),1,flwt,spwt,sm);
  CCMap[Flavour(kf::u)] = cc;  
  flwt = (1.-hadpars.Get("Baryon_supression"))*hadpars.Get("Strange_supression");
  cc = new ConstituentCharacteristic(hadpars.Get("Mass_strange"),1,flwt,spwt,sm);
  CCMap[Flavour(kf::s)] = cc;

  if (no_diquarks) return;

  // Light Di-quarks, spin 0
  spwt = 1.;
  flwt = hadpars.Get("Baryon_supression");
  cc = new ConstituentCharacteristic(hadpars.Get("Mass_ud0"),0,flwt,spwt,sm);
  CCMap[Flavour(kf::ud_0)] = cc;
  flwt = hadpars.Get("Baryon_supression")*hadpars.Get("P_qs_by_P_qq");;
  cc = new ConstituentCharacteristic(hadpars.Get("Mass_sd0"),0,flwt,spwt,sm);
  CCMap[Flavour(kf::sd_0)] = cc;
  cc = new ConstituentCharacteristic(hadpars.Get("Mass_su0"),0,flwt,spwt,sm);
  CCMap[Flavour(kf::su_0)] = cc;

  // Light Di-quarks, spin 1
  spwt = 3.*hadpars.Get("P_di_1_by_P_di_0");
  flwt = hadpars.Get("Baryon_supression");
  cc = new ConstituentCharacteristic(hadpars.Get("Mass_dd1"),2,flwt,spwt,sm);
  CCMap[Flavour(kf::dd_1)] = cc;
  cc = new ConstituentCharacteristic(hadpars.Get("Mass_ud1"),2,flwt,spwt,sm);
  CCMap[Flavour(kf::ud_1)] = cc;
  cc = new ConstituentCharacteristic(hadpars.Get("Mass_uu1"),2,flwt,spwt,sm);
  CCMap[Flavour(kf::uu_1)] = cc;
  flwt = hadpars.Get("Baryon_supression")*hadpars.Get("P_qs_by_P_qq");;
  cc = new ConstituentCharacteristic(hadpars.Get("Mass_sd1"),2,flwt,spwt,sm);
  CCMap[Flavour(kf::sd_1)] = cc;
  cc = new ConstituentCharacteristic(hadpars.Get("Mass_su1"),2,flwt,spwt,sm);
  CCMap[Flavour(kf::su_1)] = cc;
  flwt = hadpars.Get("Baryon_supression")*hadpars.Get("P_ss_by_P_qq");;
  cc = new ConstituentCharacteristic(hadpars.Get("Mass_ss1"),2,flwt,spwt,sm);
  CCMap[Flavour(kf::ss_1)] = cc;


  // Heavy Quark flavours : Won't show up in cluster break up => flwt = spwt = 0.
  flwt = spwt = 0.;
  cc = new ConstituentCharacteristic(hadpars.Get("Mass_charm"),1,flwt,spwt,sm);
  CCMap[Flavour(kf::c)] = cc;
  cc = new ConstituentCharacteristic(hadpars.Get("Mass_bottom"),1,flwt,spwt,sm);
  CCMap[Flavour(kf::b)] = cc;
  cc = new ConstituentCharacteristic(hadpars.Get("Mass_cd0"),0,flwt,spwt,sm);
  CCMap[Flavour(kf::cd_0)] = cc;
  cc = new ConstituentCharacteristic(hadpars.Get("Mass_cu0"),0,flwt,spwt,sm);
  CCMap[Flavour(kf::cu_0)] = cc;
  cc = new ConstituentCharacteristic(hadpars.Get("Mass_cs0"),0,flwt,spwt,sm);
  CCMap[Flavour(kf::cs_0)] = cc;
  cc = new ConstituentCharacteristic(hadpars.Get("Mass_cd1"),0,flwt,spwt,sm);
  CCMap[Flavour(kf::cd_1)] = cc;
  cc = new ConstituentCharacteristic(hadpars.Get("Mass_cu1"),0,flwt,spwt,sm);
  CCMap[Flavour(kf::cu_1)] = cc;
  cc = new ConstituentCharacteristic(hadpars.Get("Mass_cs1"),0,flwt,spwt,sm);
  CCMap[Flavour(kf::cs_1)] = cc;
  cc = new ConstituentCharacteristic(hadpars.Get("Mass_cc1"),0,flwt,spwt,sm);
  CCMap[Flavour(kf::cc_1)] = cc;
  cc = new ConstituentCharacteristic(hadpars.Get("Mass_bd0"),0,flwt,spwt,sm);
  CCMap[Flavour(kf::bd_0)] = cc;
  cc = new ConstituentCharacteristic(hadpars.Get("Mass_bu0"),0,flwt,spwt,sm);
  CCMap[Flavour(kf::bu_0)] = cc;
  cc = new ConstituentCharacteristic(hadpars.Get("Mass_bs0"),0,flwt,spwt,sm);
  CCMap[Flavour(kf::bs_0)] = cc;
  cc = new ConstituentCharacteristic(hadpars.Get("Mass_bc0"),0,flwt,spwt,sm);
  CCMap[Flavour(kf::bc_0)] = cc;
  cc = new ConstituentCharacteristic(hadpars.Get("Mass_bd1"),0,flwt,spwt,sm);
  CCMap[Flavour(kf::bd_1)] = cc;
  cc = new ConstituentCharacteristic(hadpars.Get("Mass_bu1"),0,flwt,spwt,sm);
  CCMap[Flavour(kf::bu_1)] = cc;
  cc = new ConstituentCharacteristic(hadpars.Get("Mass_bs1"),0,flwt,spwt,sm);
  CCMap[Flavour(kf::bs_1)] = cc;
  cc = new ConstituentCharacteristic(hadpars.Get("Mass_bc1"),0,flwt,spwt,sm);
  CCMap[Flavour(kf::bc_1)] = cc;
  cc = new ConstituentCharacteristic(hadpars.Get("Mass_bb1"),0,flwt,spwt,sm);
  CCMap[Flavour(kf::bb_1)] = cc;
}

ConstituentMasses::~ConstituentMasses() {
  for (FlavCCMap_Iterator cmit=CCMap.begin(); cmit!=CCMap.end();cmit++) {
    if (cmit->second) { delete cmit->second; cmit->second=NULL; }
  }
  CCMap.clear();
}

double ConstituentMasses::Mass(Flavour & flav) {
  FlavCCMap_Iterator     cmit = CCMap.find(flav);
  if (cmit==CCMap.end()) cmit = CCMap.find(flav.Bar());
  return (cmit!=CCMap.end())? cmit->second->Mass() : 1.e36;
}

double ConstituentMasses::FlWeight(ATOOLS::Flavour & flav) {
  FlavCCMap_Iterator     cmit = CCMap.find(flav);
  if (cmit==CCMap.end()) cmit = CCMap.find(flav.Bar());
  return (cmit!=CCMap.end())? cmit->second->FlWeight() : 0.;
}

double ConstituentMasses::SpWeight(ATOOLS::Flavour & flav) {
  FlavCCMap_Iterator     cmit = CCMap.find(flav);
  if (cmit==CCMap.end()) cmit = CCMap.find(flav.Bar());
  return (cmit!=CCMap.end())? cmit->second->SpWeight() : 0.;
}

double ConstituentMasses::TotWeight(ATOOLS::Flavour & flav) {
  FlavCCMap_Iterator     cmit = CCMap.find(flav);
  if (cmit==CCMap.end()) cmit = CCMap.find(flav.Bar());
  return (cmit!=CCMap.end())? cmit->second->TotWeight() : 0.;
}

double ConstituentMasses::Smearing(Flavour & flav) {
  FlavCCMap_Iterator     cmit = CCMap.find(flav);
  if (cmit==CCMap.end()) cmit = CCMap.find(flav.Bar());
  return (cmit!=CCMap.end())? cmit->second->Smearing() : 0.;
}

int ConstituentMasses::ISpin(Flavour & flav) {
  FlavCCMap_Iterator     cmit = CCMap.find(flav);
  if (cmit==CCMap.end()) cmit = CCMap.find(flav.Bar());
  return (cmit!=CCMap.end())? cmit->second->ISpin() : 0;
}

void ConstituentMasses::PrintConstituents() {
  for (FlavCCMap_Iterator cmit=CCMap.begin();cmit!=CCMap.end();cmit++) {
    msg.Out()<<cmit->first<<" : "<<cmit->second->m_mass<<" GeV, "
	     <<"Spin = "<<double(cmit->second->m_ispin/2.)<<", "
	     <<"Weights = "<<cmit->second->m_flweight<<", "
	     <<cmit->second->m_spweight<<std::endl;
  }
}

//##############################################################################
//##############################################################################
//##############################################################################
//##############################################################################
//##############################################################################


Wave_Function::Wave_Function() :
  m_hadron(Flavour(kf::none)), m_spin(-1), m_barrable(false)
{ }

Wave_Function::Wave_Function(const ATOOLS::Flavour & _hadron) :
  m_hadron(_hadron), m_spin(-1), m_barrable(false)
{ }

void Wave_Function::AddToWaves(FlavourPair pair,double weight)
{
  if (m_waves.find(pair)==m_waves.end()) m_waves[pair] = weight;
  else {
    msg.Error()<<"Potential error in Wave_Function::AddToWaves"<<std::endl
	       <<"   Pair "<<pair.first<<"/"<<pair.second<<" already in map."<<std::endl
	       <<"   Will ignore this and continue."<<std::endl;
    return;
  }
  if (pair.first!=pair.second.Bar()) m_barrable = true;
}

std::ostream & AHADIC::operator<<(std::ostream & s, Wave_Function & wf) 
{
  s<<" "<<wf.m_hadron<<" ("<<wf.m_kfcode<<"), spin = "<<(wf.m_spin-1)/2.<<std::endl;
  WFcomponent * waves = wf.GetWaves();
  for (WFcompiter wfc=waves->begin();wfc!=waves->end();wfc++) {
    s<<"     "<<wfc->first.first<<" "<<wfc->first.second<<" : "<<wfc->second<<std::endl;
  }
  return s;
}

Wave_Function * Wave_Function::Anti() {
  Wave_Function * wf = NULL;
  if (m_barrable) {
    wf = new Wave_Function(m_hadron.Bar());
    wf->SetSpin(m_spin);
    wf->SetKfCode(-m_kfcode);
    FlavourPair * pair;
    for (WFcompiter wfc=m_waves.begin();wfc!=m_waves.end();wfc++) {
      pair         = new FlavourPair;
      pair->first  = wfc->first.second.Bar();
      pair->second = wfc->first.first.Bar();
      wf->AddToWaves((*pair),wfc->second);
    }
  }
  return wf;
}

std::ostream & AHADIC::operator<<(std::ostream & s, const CChannel_Element & ce) {
  s<<"  ("<<ce.m_pairs[0].first<<","<<ce.m_pairs[0].second
   <<" ; "<<ce.m_pairs[1].first<<","<<ce.m_pairs[1].second
   <<") -> ("<<ce.m_hadrons[0]<<","<<ce.m_hadrons[1]<<")"
   <<"   "<<ce.m_mass<<"  "<<ce.m_weight<<std::endl;
  return s;
}

std::ostream & AHADIC::operator<<(std::ostream & s, const ChannelSet & cs) {
  for (Channeliter ch=cs.begin();ch!=cs.end();ch++) 
    s<<"  ("<<(*ch)->m_pairs[0].first<<","<<(*ch)->m_pairs[0].second
     <<" ; "<<(*ch)->m_pairs[1].first<<","<<(*ch)->m_pairs[1].second
     <<") -> ("<<(*ch)->m_hadrons[0]<<","<<(*ch)->m_hadrons[1]<<")"
     <<"   "<<(*ch)->m_mass<<"  "<<(*ch)->m_weight<<std::endl;
  return s;
}



//##############################################################################
//##############################################################################
//##############################################################################
//##############################################################################
//##############################################################################






Hadron_List::Hadron_List() :
  p_wavefunctions(NULL), p_multiplets(NULL)
{
  ConstructWaveFunctions();
  ConstructAntiWaveFunctions();
  CreateMultiplets();
  PrintWaveFunctions();
}

Hadron_List::~Hadron_List() 
{
  if (p_wavefunctions!=NULL && !p_wavefunctions->empty()) {
    for (WFmapiter wfm=p_wavefunctions->begin();wfm!=p_wavefunctions->end();wfm++) {
      if (wfm->second!=NULL) { delete wfm->second; wfm->second=NULL; }
    }
    p_wavefunctions->clear();
    delete p_wavefunctions;
  }
  if (p_multiplets!=NULL && !p_multiplets->empty()) {
    for (Mpletmapiter mp=p_multiplets->begin();mp!=p_multiplets->end();mp++) {
      if (mp->second!=NULL) { delete mp->second; mp->second=NULL; }
    }
    p_multiplets->clear();
    delete p_multiplets;
  }
}

void Hadron_List::ConstructWaveFunctions() 
{
  Flavour         hadron;
  Wave_Function * wavefunction;
  int             spin, lp, kfcode, fl1, fl2, fl3;
  double          costh, sinth, weight;
  Flavour         trip,antitrip,flavs[3];
  FlavourPair   * pair;

  p_wavefunctions = new WFMap; 
  for (int i=0;i<kf_table.Size();i++) {
    kfcode = kf_table.FromInt(i);
    hadron = Flavour(kf::code(kfcode));
    if (!hadron.IsHadron() || !hadron.IsOn()) continue;
    wavefunction = new Wave_Function(hadron);
    spin = (kfcode-10*int(kfcode/10));
    lp   = int(kfcode/10000);
    fl1  = int(kfcode/10)-int(kfcode/100)*10;       
    fl2  = int(kfcode/100)-int(kfcode/1000)*10;       
    fl3  = int(kfcode/1000)-int(kfcode/10000)*10;       
    

    if (fl3==0) {
      flavs[0]     = Flavour(kf::code(fl1)).Bar();
      flavs[1]     = Flavour(kf::code(fl2));
      pair         = new FlavourPair;
      pair->first  = flavs[1];
      pair->second = flavs[0];
      if (fl1!=fl2 && (flavs[0].Charge()+flavs[1].Charge()!=0.)) 
	wavefunction->AddToWaves((*pair),1.);
      else {
	if (fl1==fl2) {
	  LookUpAngles(lp,spin,costh,sinth);
	  if (fl1==1) {
	    wavefunction->AddToWaves((*pair),-1./sqrt(2.));
	    flavs[0]     = Flavour(kf::code(fl1+1));
	    pair         = new FlavourPair;
	    pair->first  = flavs[0];
	    pair->second = flavs[0].Bar();
	    wavefunction->AddToWaves((*pair),+1./sqrt(2.));
	  } 
	  if (fl1==2) {
	    weight         = costh/sqrt(3.)+sinth/sqrt(6.);
	    if (weight>1.e-3) {
	      wavefunction->AddToWaves((*pair),weight);
	    }
	    weight         = costh/sqrt(3.)+sinth/sqrt(6.);
	    if (weight>1.e-3) {
	      flavs[0]     = Flavour(kf::code(1));
	      pair         = new FlavourPair;
	      pair->first  = flavs[0];
	      pair->second = flavs[0].Bar();
	      wavefunction->AddToWaves((*pair),weight);
	    }
	    weight         = costh/sqrt(3.)-2.*sinth/sqrt(6.);
	    if (weight>1.e-3) {
	      flavs[0]     = Flavour(kf::code(3));
	      pair         = new FlavourPair;
	      pair->first  = flavs[0];
	      pair->second = flavs[0].Bar();
	      wavefunction->AddToWaves((*pair),weight);
	    }
	  } 
	  if (fl1==3) {
	    weight         = 2.*costh/sqrt(6.)+sinth/sqrt(3.);
	    if (weight>1.e-3) {
	      wavefunction->AddToWaves((*pair),weight);
	    }
	    weight         = costh/sqrt(6.)-sinth/sqrt(3.);
	    if (weight>1.e-3) {
	      flavs[0]     = Flavour(kf::code(1));
	      pair         = new FlavourPair;
	      pair->first  = flavs[0];
	      pair->second = flavs[0].Bar();
	      wavefunction->AddToWaves((*pair),weight);
	    }
	    weight         = costh/sqrt(6.)-sinth/sqrt(3.);
	    if (weight>1.e-3) {
	      flavs[0]     = Flavour(kf::code(2));
	      pair         = new FlavourPair;
	      pair->first  = flavs[0];
	      pair->second = flavs[0].Bar();
	      wavefunction->AddToWaves((*pair),weight);
	    }
	  } 
	}
	else if (flavs[0].Charge()+flavs[1].Charge()==0.) {
	  wavefunction->AddToWaves((*pair),1.);
	}
      }
    }

    /*
      std::cout<<i<<":"<<hadron<<":"<<kfcode<<"->"<<spin
      <<"("<<fl1<<","<<fl2<<","<<fl3<<")"
      <<"-> ("<<flavs[0]<<","<<flavs[1]<<")"<<std::endl;
    */
    if (spin>0) {
      wavefunction->SetSpin(spin);
      wavefunction->SetKfCode(kfcode);
      (*p_wavefunctions)[hadron] = wavefunction;
    }
    else delete wavefunction;
  }
}

void Hadron_List::ConstructAntiWaveFunctions() 
{
  Wave_Function * anti; 
  for (WFmapiter wfm=p_wavefunctions->begin();wfm!=p_wavefunctions->end();wfm++) {
    anti = wfm->second->Anti();
    if (anti!=NULL) (*p_wavefunctions)[wfm->first.Bar()] = anti;
  } 
}

void Hadron_List::LookUpAngles(const int angular,const int spin,double & costh,double & sinth)
{
  double angle=0.;
  switch (spin) {
  case 7 : hadpars.Get("Mixing_Angle_3+"); break;
  case 5 : hadpars.Get("Mixing_Angle_2-"); break;
  case 3 : hadpars.Get("Mixing_Angle_1+"); break;
  case 1 : hadpars.Get("Mixing_Angle_0-"); break; 
  default: break;  
  }
  costh = cos(angle); sinth = sin(angle);
}

void Hadron_List::CreateMultiplets()
{
  p_multiplets = new MultipletMap;
  int                kfcode,spin;
  Mpletmapiter       mpletiter;
  Hadron_Multiplet * multiplet;
  for (WFmapiter wfm=p_wavefunctions->begin();wfm!=p_wavefunctions->end();wfm++) {
    kfcode = wfm->second->KfCode();
    spin   = wfm->second->Spin();
    if (kfcode<1000) {
      mpletiter = p_multiplets->find(spin);
      if (mpletiter!=p_multiplets->end()) {
	//msg.Debugging()<<"Add "<<wfm->first<<" to "<<mpletiter->second->m_name<<std::endl;
	mpletiter->second->m_elements.insert(wfm->first);
      }
      else {
	multiplet = new Hadron_Multiplet;
	multiplet->m_weight = 1.;
	switch (spin) {
	case 1:  multiplet->m_name=std::string("Pseudoscalars (0-)"); break;
	case 3:  multiplet->m_name=std::string("Vectors       (1+)"); break;
	case 5:  multiplet->m_name=std::string("Tensors       (2-)"); break;
	case 7:  multiplet->m_name=std::string("Tensors       (3+)"); break;
	default: multiplet->m_name=std::string("Don't know   ");      break;
	}
	multiplet->m_elements.insert(wfm->first);
	//msg.Debugging()<<"Add "<<wfm->first<<" to new "<<multiplet->m_name<<std::endl;
	(*p_multiplets)[spin] = multiplet;
      }
    }
  }
} 

void Hadron_List::PrintWaveFunctions() 
{
  WFmapiter       wfm;
  Wave_Function * wf;
  for (Mpletmapiter mplet=p_multiplets->begin();mplet!=p_multiplets->end();mplet++) {
    msg.Out()<<"-----------------------------------------------"<<std::endl
	     <<" "<<mplet->second->m_name<<" with "
	     <<mplet->second->m_elements.size()<<" elements: "<<std::endl;
    for (FlSetIter fl=mplet->second->m_elements.begin();
	 fl!=mplet->second->m_elements.end();fl++) {
      wfm = p_wavefunctions->find((*fl));
      if (wfm!=p_wavefunctions->end()) { wf = wfm->second; msg.Out()<<(*wf); }
    }
    msg.Out()<<"-----------------------------------------------"<<std::endl;
  }
}


//##############################################################################
//##############################################################################
//##############################################################################
//##############################################################################
//##############################################################################




Hadron_Parameters::Hadron_Parameters() :
  p_constituents(NULL), p_hadrons(NULL), 
  p_allweights(NULL), p_allchannels(NULL)
{ }

Hadron_Parameters::~Hadron_Parameters() {
  if (p_constituents!=NULL) { 
    delete p_constituents; p_constituents=NULL; 
  }
  if (p_allweights!=NULL && !p_allweights->empty()) {
    for (Allweightiter all=p_allweights->begin();all!=p_allweights->end();all++) {
      if (all->second!=NULL) { delete all->second; all->second=NULL; }
    }
    p_allweights->clear();
    delete p_allweights;
  }
  if (p_allchannels!=NULL && !p_allchannels->empty()) {
    for (Allchanneliter all=p_allchannels->begin();all!=p_allchannels->end();all++) {
      if (all->second!=NULL) { delete all->second; all->second=NULL; }
    }
    p_allchannels->clear();
    delete p_allchannels;
  }
}

void Hadron_Parameters::Init(std::string dir,std::string file)
{
  std::cout<<"In Hadron_Parameters::Init("<<dir<<file<<")"<<std::endl;
  ReadParameters(dir,file);
  p_constituents = new ConstituentMasses(true);
  p_constituents->PrintConstituents();
  p_hadrons      = new Hadron_List;
  CreateWeightLists();  
  PrintWeightLists();  
  CreateChannelLists();  
  PrintChannelLists();  
}
  
void Hadron_Parameters::ReadParameters(std::string dir,std::string file)
{
  Data_Read dataread(dir+file);
  m_parametermap[std::string("Strange_supression")] =
    dataread.GetValue<double>("STRANGE_SUPRESSION",0.2);      
  m_parametermap[std::string("Baryon_supression")]  = 
    dataread.GetValue<double>("BARYON_SUPRESSION",0.2);
  m_parametermap[std::string("Q_breakup")]          = 
    dataread.GetValue<double>("Q_BREAKUP",1.);
  m_parametermap[std::string("AngularSmearing")]    = 
    dataread.GetValue<double>("AngularSmearing",Get("Q_breakup"));
  m_parametermap[std::string("Max_Prod_Mass")]      = 
    dataread.GetValue<double>("Max_Prod_Mass",Get("Q_breakup"));
  m_parametermap[std::string("P_qs_by_P_qq")]       = 
    dataread.GetValue<double>("P_{QS}/P_{QQ}",Get("Strange_supression"));
  m_parametermap[std::string("P_ss_by_P_qq")]       = 
    dataread.GetValue<double>("P_{SS}/P_{QQ}",sqr(Get("Strange_supression")));    
  m_parametermap[std::string("P_di_1_by_P_di_0")]   = 
    dataread.GetValue<double>("P_{QQ_1}/P_{QQ_0}",1.);
  m_parametermap[std::string("Mass_glue")]          = 
    dataread.GetValue<double>("M_GLUE",0.75);
  m_parametermap[std::string("Mass_down")]          = 
    dataread.GetValue<double>("M_DOWN",0.32);
  m_parametermap[std::string("Mass_up")]            = 
    dataread.GetValue<double>("M_UP",0.32);
  m_parametermap[std::string("Mass_strange")]       = 
    dataread.GetValue<double>("M_STRANGE",0.45);
  m_parametermap[std::string("Mass_charm")]         = 
    dataread.GetValue<double>("M_CHARM",1.55);
  m_parametermap[std::string("Mass_bottom")]        = 
    dataread.GetValue<double>("M_BOTTOM",4.55);
  m_parametermap[std::string("Mass_dd1")]           = 
    dataread.GetValue<double>("M_DD_1",2.*Get("Mass_down"));
  m_parametermap[std::string("Mass_ud0")]           = 
    dataread.GetValue<double>("M_UD_0",Get("Mass_up")+Get("Mass_down"));
  m_parametermap[std::string("Mass_ud1")]           = 
    dataread.GetValue<double>("M_UD_1",Get("Mass_up")+Get("Mass_down"));
  m_parametermap[std::string("Mass_uu1")]           = 
    dataread.GetValue<double>("M_UU_1",2.*Get("Mass_up"));
  m_parametermap[std::string("Mass_sd0")]           = 
    dataread.GetValue<double>("M_SD_0",Get("Mass_strange")+Get("Mass_down"));
  m_parametermap[std::string("Mass_sd1")]           = 
    dataread.GetValue<double>("M_SD_1",Get("Mass_strange")+Get("Mass_down"));
  m_parametermap[std::string("Mass_su0")]           = 
    dataread.GetValue<double>("M_SU_0",Get("Mass_strange")+Get("Mass_up"));
  m_parametermap[std::string("Mass_su1")]           = 
    dataread.GetValue<double>("M_SU_1",Get("Mass_strange")+Get("Mass_up"));
  m_parametermap[std::string("Mass_ss1")]           = 
    dataread.GetValue<double>("M_SS_1",2.*Get("Mass_strange"));
  m_parametermap[std::string("Mass_cd0")]           = 
    dataread.GetValue<double>("M_CD_0",Get("Mass_charm")+Get("Mass_down"));
  m_parametermap[std::string("Mass_cd1")]           = 
    dataread.GetValue<double>("M_CD_1",Get("Mass_charm")+Get("Mass_down"));
  m_parametermap[std::string("Mass_cu0")]           = 
    dataread.GetValue<double>("M_CU_0",Get("Mass_charm")+Get("Mass_up"));
  m_parametermap[std::string("Mass_cu1")]           = 
    dataread.GetValue<double>("M_CU_1",Get("Mass_charm")+Get("Mass_up"));
  m_parametermap[std::string("Mass_cs0")]           = 
    dataread.GetValue<double>("M_CS_0",Get("Mass_charm")+Get("Mass_strange"));
  m_parametermap[std::string("Mass_cs1")]           = 
    dataread.GetValue<double>("M_CS_1",Get("Mass_charm")+Get("Mass_strange"));
  m_parametermap[std::string("Mass_cc1")]           = 
    dataread.GetValue<double>("M_CC_1",2.*Get("Mass_charm"));
  m_parametermap[std::string("Mass_bd0")]           = 
    dataread.GetValue<double>("M_BD_0",Get("Mass_bottom")+Get("Mass_down"));
  m_parametermap[std::string("Mass_bd1")]           = 
    dataread.GetValue<double>("M_BD_1",Get("Mass_bottom")+Get("Mass_down"));
  m_parametermap[std::string("Mass_bu0")]           = 
    dataread.GetValue<double>("M_BU_0",Get("Mass_bottom")+Get("Mass_up"));
  m_parametermap[std::string("Mass_bu1")]           = 
    dataread.GetValue<double>("M_BU_1",Get("Mass_bottom")+Get("Mass_up"));
  m_parametermap[std::string("Mass_bs0")]           = 
    dataread.GetValue<double>("M_BS_0",Get("Mass_bottom")+Get("Mass_strange"));
  m_parametermap[std::string("Mass_bs1")]           = 
    dataread.GetValue<double>("M_BS_1",Get("Mass_bottom")+Get("Mass_strange"));
  m_parametermap[std::string("Mass_bc0")]           = 
    dataread.GetValue<double>("M_BC_0",Get("Mass_bottom")+Get("Mass_charm"));
  m_parametermap[std::string("Mass_bc1")]           = 
    dataread.GetValue<double>("M_BC_1",Get("Mass_bottom")+Get("Mass_charm"));
  m_parametermap[std::string("Mass_bb1")]           = 
    dataread.GetValue<double>("M_BB_1",2.*Get("Mass_bottom"));
  m_parametermap[std::string("Mixing_Angle_0-")]    = 
    dataread.GetValue<double>("Mixing_0-",-0.3491);
  m_parametermap[std::string("Mixing_Angle_1+")]    = 
    dataread.GetValue<double>("Mixing_0-",0.6155);
}

double Hadron_Parameters::Get(std::string keyword) 
{
  m_piter = m_parametermap.find(keyword);
  if (m_piter!=m_parametermap.end()) return m_piter->second;
  msg.Error()<<"Error in Hadron_Parameters::Get("<<keyword<<") : "<<std::endl
	     <<"   Keyword not found. Return 0 and hope for the best."<<std::endl;
  return 0.;
}

double Hadron_Parameters::MinimalMass(FlavourPair * pair)
{
  Allchanneliter all = p_allchannels->find((*pair));
  if (all!=p_allchannels->end()) return (*all->second->begin())->m_mass;

  msg.Error()<<"Error in Hadron_List::MinimalMass("<<pair->first<<","<<pair->second<<") :"
	     <<"   Pair not found. Abort the run."<<std::endl;
  abort();
}

void Hadron_Parameters::PrintWeightLists()
{
  Weightiter one;
  for (Allweightiter all=p_allweights->begin();all!=p_allweights->end();all++) {
    msg.Out()<<"-----------------------------------------------"<<std::endl
	     <<" "<<all->first.first<<"+"<<all->first.second<<" with "
	     <<all->second->size()<<" elements: "<<std::endl;
    for (one=all->second->begin();one!=all->second->end();one++) {
      msg.Out()<<"   "<<one->first<<"   : "<<one->second<<std::endl;
    }
    msg.Out()<<"-----------------------------------------------"<<std::endl;
  }
}

void Hadron_Parameters::PrintChannelLists()
{
  for (Allchanneliter all=p_allchannels->begin();all!=p_allchannels->end();all++) {
    msg.Out()<<"-----------------------------------------------"<<std::endl
	     <<" "<<all->first.first<<"+"<<all->first.second<<" with "
	     <<all->second->size()<<" elements: "<<std::endl<<(*all->second);
    msg.Out()<<"-----------------------------------------------"<<std::endl;
  }
}


void Hadron_Parameters::CreateWeightLists() 
{
  p_allweights   = new AllWeights;
  WFcomponent   * wf;
  for (WFmapiter wfm=p_hadrons->GetWaveFunctions()->begin();
       wfm!=p_hadrons->GetWaveFunctions()->end();wfm++) {
    wf     = wfm->second->GetWaves();
    for (WFcompiter wfc=wf->begin();wfc!=wf->end();wfc++) {
      AddToWeights(wfc->first,wfm->first,wfc->second);
    }
  }
}

void Hadron_Parameters::AddToWeights(FlavourPair flpair,Flavour had,double wt) 
{
  Allweightiter all = p_allweights->find(flpair);
  if (all!=p_allweights->end()) all->second->insert(std::make_pair(had,wt));
  else {
    WeightMap * wmap = new WeightMap;
    wmap->insert(std::make_pair(had,wt));
    (*p_allweights)[flpair] = wmap;
  }
}

void Hadron_Parameters::CreateChannelLists() 
{
  p_allchannels    = new AllChannels;
  FlavCCMap          flavmassmap = p_constituents->CCMap;
  Allchanneliter     chan;
  ChannelSet       * channels;
  CChannel_Element * channel;
  FlavourPair        pair,pair1,pair2;
  Flavour            flav;
  Allweightiter      all, all1, all2;
  Weightiter         wgt1, wgt2;
  double             mass;

  for (all=p_allweights->begin();all!=p_allweights->end();all++) {
    pair         = all->first;
    pair1.first  = pair.first;
    pair2.second = pair.second;
    channels     = new ChannelSet;
    chan         = p_allchannels->find(pair);
    if (chan!=p_allchannels->end()) {
      msg.Error()<<"Error in Hadron_List::CreateChannelLists :"<<std::endl
		 <<"   Pair already found in channel lists : "
		 <<"{"<<pair.first<<","<<pair.second<<"}"<<std::endl
		 <<"   Will abort the run."<<std::endl;
      abort();
    }
    p_allchannels->insert(std::make_pair(pair, channels));
    for (FlavCCMap_Iterator fliter=flavmassmap.begin();
	 fliter!=flavmassmap.end();fliter++) {
      flav = fliter->first;
      if ((flav.IsQuark() || flav.IsDiQuark()) && 
	  fliter->second->Mass()<Get("Max_Prod_Mass")) {
	pair1.second = fliter->first.Bar();
	pair2.first  = fliter->first;
      }
      else continue;
      all1 = p_allweights->find(pair1);
      all2 = p_allweights->find(pair2);
      if (all1==p_allweights->end() || all2==p_allweights->end()) {
	msg.Error()<<"Error in Hadron_List::CreateChannelLists :"<<std::endl
		   <<"   Pairs not found in weight lists : "<<std::endl
		   <<"{"<<pair.first<<","<<pair.second<<"} -> "
		   <<"{"<<pair1.first<<","<<pair1.second<<"}, "
		   <<"{"<<pair2.first<<","<<pair2.second<<"}"<<std::endl
		   <<"   Will continue and neglect this pair."<<std::endl;
	continue;
      }
      for (wgt1=all1->second->begin();wgt1!=all1->second->end();wgt1++) {
	for (wgt2=all2->second->begin();wgt2!=all2->second->end();wgt2++) {
	  channel               = new CChannel_Element;
	  channel->m_pairs[0]   = pair1;
	  channel->m_pairs[1]   = pair2;
	  channel->m_hadrons[0] = wgt1->first;
	  channel->m_hadrons[1] = wgt2->first;
	  mass                  = Max(wgt1->first.Mass()+wgt2->first.Mass(),
				      p_constituents->Mass(pair1.first)+
				      p_constituents->Mass(pair1.second)+
				      p_constituents->Mass(pair2.first)+
				      p_constituents->Mass(pair2.second));
	  channel->m_mass       = 1.1*mass;
	  channel->m_weight     = wgt1->second*wgt2->second;
	  channels->insert(channel);
	}
      }
    }
  }
}




