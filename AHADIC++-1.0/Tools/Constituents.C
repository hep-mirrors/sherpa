#include "Hadronisation_Parameters.H"
#include "Constituents.H"
#include "Data_Read.H"
#include "Message.H"

using namespace AHADIC;
using namespace ATOOLS;

Constituents::Constituents(bool no_diquarks) :
  m_minmass(100.)
{
  bool no_heavies = true;

  ConstituentCharacteristic * cc;
  double flwt=0., spwt=0., sm = hadpars.Get("AngularSmearing");

  // Gluon
  cc = new ConstituentCharacteristic(hadpars.Get("Mass_glue"),2,flwt,spwt);
  CCMap[Flavour(kf::gluon)] = cc;

  // Light quarks
  double sfrac(hadpars.Get("Strange_fraction")), bfrac(hadpars.Get("Baryon_fraction"));
  double wt(1.-bfrac);
  cc = new ConstituentCharacteristic(hadpars.Get("Mass_down"),1,wt*(1.-sfrac)/2.,sm);
  CCMap[Flavour(kf::d)] = cc;
  cc = new ConstituentCharacteristic(hadpars.Get("Mass_up"),1,wt*(1.-sfrac)/2.,sm);
  CCMap[Flavour(kf::u)] = cc;  
  cc = new ConstituentCharacteristic(hadpars.Get("Mass_strange"),1,wt*sfrac,sm);
  CCMap[Flavour(kf::s)] = cc;

  cc = new ConstituentCharacteristic(hadpars.Get("Mass_charm"),1,0.,sm);
  CCMap[Flavour(kf::c)] = cc;
  cc = new ConstituentCharacteristic(hadpars.Get("Mass_bottom"),1,0.,sm);
  CCMap[Flavour(kf::b)] = cc;

  if (!no_diquarks) {
    // Light Di-quarks, spin 0
    double qssup(hadpars.Get("P_qs_by_P_qq"));
    double sssup(hadpars.Get("P_ss_by_P_qq"));
    double sp1sup(hadpars.Get("P_di_1_by_P_di_0"));
    wt = bfrac/(1.*(1.+2.*qssup) + 3.*sp1sup*(3.+2.*qssup+sssup));
    
    cc = new ConstituentCharacteristic(hadpars.Get("Mass_ud0"),0,wt,sm);
    CCMap[Flavour(kf::ud_0)] = cc;
    cc = new ConstituentCharacteristic(hadpars.Get("Mass_sd0"),0,wt*qssup,sm);
    CCMap[Flavour(kf::sd_0)] = cc;
    cc = new ConstituentCharacteristic(hadpars.Get("Mass_su0"),0,wt*qssup,sm);
    CCMap[Flavour(kf::su_0)] = cc;
    
    // Light Di-quarks, spin 1
    cc = new ConstituentCharacteristic(hadpars.Get("Mass_dd1"),2,3.*wt*sp1sup,sm);
    CCMap[Flavour(kf::dd_1)] = cc;
    cc = new ConstituentCharacteristic(hadpars.Get("Mass_ud1"),2,3.*wt*sp1sup,sm);
    CCMap[Flavour(kf::ud_1)] = cc;
    cc = new ConstituentCharacteristic(hadpars.Get("Mass_uu1"),2,3.*wt*sp1sup,sm);
    CCMap[Flavour(kf::uu_1)] = cc;
    cc = new ConstituentCharacteristic(hadpars.Get("Mass_sd1"),2,3.*wt*sp1sup*qssup,sm);
    CCMap[Flavour(kf::sd_1)] = cc;
    cc = new ConstituentCharacteristic(hadpars.Get("Mass_su1"),2,3.*wt*sp1sup*qssup,sm);
    CCMap[Flavour(kf::su_1)] = cc;
    cc = new ConstituentCharacteristic(hadpars.Get("Mass_ss1"),2,3.*wt*sp1sup*sssup,sm);
    CCMap[Flavour(kf::ss_1)] = cc;

    if (!no_heavies) {
      cc = new ConstituentCharacteristic(hadpars.Get("Mass_cd0"),0,0.,sm);
      CCMap[Flavour(kf::cd_0)] = cc;
      cc = new ConstituentCharacteristic(hadpars.Get("Mass_cu0"),0,0.,sm);
      CCMap[Flavour(kf::cu_0)] = cc;
      cc = new ConstituentCharacteristic(hadpars.Get("Mass_cs0"),0,0.,sm);
      CCMap[Flavour(kf::cs_0)] = cc;
      cc = new ConstituentCharacteristic(hadpars.Get("Mass_cd1"),0,0.,sm);
      CCMap[Flavour(kf::cd_1)] = cc;
      cc = new ConstituentCharacteristic(hadpars.Get("Mass_cu1"),0,0.,sm);
      CCMap[Flavour(kf::cu_1)] = cc;
      cc = new ConstituentCharacteristic(hadpars.Get("Mass_cs1"),0,0.,sm);
      CCMap[Flavour(kf::cs_1)] = cc;
      cc = new ConstituentCharacteristic(hadpars.Get("Mass_cc1"),0,0.,sm);
      CCMap[Flavour(kf::cc_1)] = cc;
      cc = new ConstituentCharacteristic(hadpars.Get("Mass_bd0"),0,0.,sm);
      CCMap[Flavour(kf::bd_0)] = cc;
      cc = new ConstituentCharacteristic(hadpars.Get("Mass_bu0"),0,0.,sm);
      CCMap[Flavour(kf::bu_0)] = cc;
      cc = new ConstituentCharacteristic(hadpars.Get("Mass_bs0"),0,0.,sm);
      CCMap[Flavour(kf::bs_0)] = cc;
      cc = new ConstituentCharacteristic(hadpars.Get("Mass_bc0"),0,0.,sm);
      CCMap[Flavour(kf::bc_0)] = cc;
      cc = new ConstituentCharacteristic(hadpars.Get("Mass_bd1"),0,0.,sm);
      CCMap[Flavour(kf::bd_1)] = cc;
      cc = new ConstituentCharacteristic(hadpars.Get("Mass_bu1"),0,0.,sm);
      CCMap[Flavour(kf::bu_1)] = cc;
      cc = new ConstituentCharacteristic(hadpars.Get("Mass_bs1"),0,0.,sm);
      CCMap[Flavour(kf::bs_1)] = cc;
      cc = new ConstituentCharacteristic(hadpars.Get("Mass_bc1"),0,0.,sm);
      CCMap[Flavour(kf::bc_1)] = cc;
      cc = new ConstituentCharacteristic(hadpars.Get("Mass_bb1"),0,0.,sm);
      CCMap[Flavour(kf::bb_1)] = cc;
    }
  }
  for (FlavCCMap_Iterator cmit=CCMap.begin(); cmit!=CCMap.end();cmit++) {
    if (cmit->second->Mass()<m_minmass) m_minmass = cmit->second->Mass();
  }
}

Constituents::~Constituents() {
  for (FlavCCMap_Iterator cmit=CCMap.begin(); cmit!=CCMap.end();cmit++) {
    if (cmit->second) { delete cmit->second; cmit->second=NULL; }
  }
  CCMap.clear();
}

double Constituents::MinMass() { return m_minmass; }

double Constituents::Mass(Flavour & flav) {
  FlavCCMap_Iterator     cmit = CCMap.find(flav);
  if (cmit==CCMap.end()) cmit = CCMap.find(flav.Bar());
  return (cmit!=CCMap.end())? cmit->second->Mass() : flav.Mass();
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
  double wt(0.);
  for (FlavCCMap_Iterator cmit=CCMap.begin();cmit!=CCMap.end();cmit++) {
    wt+=cmit->second->m_weight;
    msg_Out()<<cmit->first<<" : "<<cmit->second->m_mass<<" GeV, "
	     <<"Spin = "<<double(cmit->second->m_ispin/2.)<<", "
	     <<"Weight = "<<cmit->second->m_weight<<std::endl;
  }
  msg_Out()<<"Total weight : "<<wt<<std::endl
	   <<"------------- END OF CONSTITUENTS ---------------"<<std::endl;
}

