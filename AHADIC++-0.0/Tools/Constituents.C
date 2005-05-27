#include "Hadronisation_Parameters.H"
#include "Constituents.H"
#include "Data_Read.H"
#include "Message.H"

using namespace AHADIC;
using namespace ATOOLS;

Constituents::Constituents(bool no_diquarks) {
  bool no_heavies = true;

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

  if (no_heavies) return;

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

Constituents::~Constituents() {
  for (FlavCCMap_Iterator cmit=CCMap.begin(); cmit!=CCMap.end();cmit++) {
    if (cmit->second) { delete cmit->second; cmit->second=NULL; }
  }
  CCMap.clear();
}

double Constituents::Mass(Flavour & flav) {
  FlavCCMap_Iterator     cmit = CCMap.find(flav);
  if (cmit==CCMap.end()) cmit = CCMap.find(flav.Bar());
  return (cmit!=CCMap.end())? cmit->second->Mass() : 1.e36;
}

double Constituents::FlWeight(ATOOLS::Flavour & flav) {
  FlavCCMap_Iterator     cmit = CCMap.find(flav);
  if (cmit==CCMap.end()) cmit = CCMap.find(flav.Bar());
  return (cmit!=CCMap.end())? cmit->second->FlWeight() : 0.;
}

double Constituents::SpWeight(ATOOLS::Flavour & flav) {
  FlavCCMap_Iterator     cmit = CCMap.find(flav);
  if (cmit==CCMap.end()) cmit = CCMap.find(flav.Bar());
  return (cmit!=CCMap.end())? cmit->second->SpWeight() : 0.;
}

double Constituents::TotWeight(ATOOLS::Flavour & flav) {
  FlavCCMap_Iterator     cmit = CCMap.find(flav);
  if (cmit==CCMap.end()) cmit = CCMap.find(flav.Bar());
  return (cmit!=CCMap.end())? cmit->second->TotWeight() : 0.;
}

double Constituents::Smearing(Flavour & flav) {
  FlavCCMap_Iterator     cmit = CCMap.find(flav);
  if (cmit==CCMap.end()) cmit = CCMap.find(flav.Bar());
  return (cmit!=CCMap.end())? cmit->second->Smearing() : 0.;
}

int Constituents::ISpin(Flavour & flav) {
  FlavCCMap_Iterator     cmit = CCMap.find(flav);
  if (cmit==CCMap.end()) cmit = CCMap.find(flav.Bar());
  return (cmit!=CCMap.end())? cmit->second->ISpin() : 0;
}

void Constituents::PrintConstituents() {
  for (FlavCCMap_Iterator cmit=CCMap.begin();cmit!=CCMap.end();cmit++) {
    msg.Out()<<cmit->first<<" : "<<cmit->second->m_mass<<" GeV, "
	     <<"Spin = "<<double(cmit->second->m_ispin/2.)<<", "
	     <<"Weights = "<<cmit->second->m_flweight<<", "
	     <<cmit->second->m_spweight<<std::endl;
  }
  msg.Out()<<"------------- END OF CONSTITUENTS ---------------"<<std::endl;
}

