#include "AHADIC++/Tools/Hadronisation_Parameters.H"
#include "AHADIC++/Tools/Constituents.H"
#include "ATOOLS/Org/Message.H"

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
  CCMap[Flavour(kf_gluon)] = cc;

  // Light quarks
  double sfrac(hadpars.Get("Strange_fraction")), bfrac(hadpars.Get("Baryon_fraction"));
  double wt(1.-bfrac);
  cc = new ConstituentCharacteristic(hadpars.Get("Mass_down"),1,1.,sm);
  CCMap[Flavour(kf_d)] = cc;
  cc = new ConstituentCharacteristic(hadpars.Get("Mass_up"),1,1.,sm);
  CCMap[Flavour(kf_u)] = cc;  
  cc = new ConstituentCharacteristic(hadpars.Get("Mass_strange"),1,sfrac,sm);
  CCMap[Flavour(kf_s)] = cc;

  cc = new ConstituentCharacteristic(hadpars.Get("Mass_charm"),1,0.,sm);
  CCMap[Flavour(kf_c)] = cc;
  cc = new ConstituentCharacteristic(hadpars.Get("Mass_bottom"),1,0.,sm);
  CCMap[Flavour(kf_b)] = cc;

  if (!no_diquarks) {
    // Light Di-quarks, spin 0
    double qssup(hadpars.Get("P_qs_by_P_qq"));
    double sssup(hadpars.Get("P_ss_by_P_qq"));
    double sp1sup(hadpars.Get("P_di_1_by_P_di_0"));

    cc = new ConstituentCharacteristic(hadpars.Get("Mass_ud0"),0,bfrac,sm);
    CCMap[Flavour(kf_ud_0)] = cc;
    cc = new ConstituentCharacteristic(hadpars.Get("Mass_sd0"),0,bfrac*qssup,sm);
    CCMap[Flavour(kf_sd_0)] = cc;
    cc = new ConstituentCharacteristic(hadpars.Get("Mass_su0"),0,bfrac*qssup,sm);
    CCMap[Flavour(kf_su_0)] = cc;
    
    // Light Di-quarks, spin 1
    cc = new ConstituentCharacteristic(hadpars.Get("Mass_dd1"),2,3.*bfrac*sp1sup,sm);
    CCMap[Flavour(kf_dd_1)] = cc;
    cc = new ConstituentCharacteristic(hadpars.Get("Mass_ud1"),2,3.*bfrac*sp1sup,sm);
    CCMap[Flavour(kf_ud_1)] = cc;
    cc = new ConstituentCharacteristic(hadpars.Get("Mass_uu1"),2,3.*bfrac*sp1sup,sm);
    CCMap[Flavour(kf_uu_1)] = cc;
    cc = new ConstituentCharacteristic(hadpars.Get("Mass_sd1"),2,3.*bfrac*sp1sup*qssup,sm);
    CCMap[Flavour(kf_sd_1)] = cc;
    cc = new ConstituentCharacteristic(hadpars.Get("Mass_su1"),2,3.*bfrac*sp1sup*qssup,sm);
    CCMap[Flavour(kf_su_1)] = cc;
    cc = new ConstituentCharacteristic(hadpars.Get("Mass_ss1"),2,3.*bfrac*sp1sup*sssup,sm);
    CCMap[Flavour(kf_ss_1)] = cc;

    if (!no_heavies) {
      cc = new ConstituentCharacteristic(hadpars.Get("Mass_cd0"),0,0.,sm);
      CCMap[Flavour(kf_cd_0)] = cc;
      cc = new ConstituentCharacteristic(hadpars.Get("Mass_cu0"),0,0.,sm);
      CCMap[Flavour(kf_cu_0)] = cc;
      cc = new ConstituentCharacteristic(hadpars.Get("Mass_cs0"),0,0.,sm);
      CCMap[Flavour(kf_cs_0)] = cc;
      cc = new ConstituentCharacteristic(hadpars.Get("Mass_cd1"),0,0.,sm);
      CCMap[Flavour(kf_cd_1)] = cc;
      cc = new ConstituentCharacteristic(hadpars.Get("Mass_cu1"),0,0.,sm);
      CCMap[Flavour(kf_cu_1)] = cc;
      cc = new ConstituentCharacteristic(hadpars.Get("Mass_cs1"),0,0.,sm);
      CCMap[Flavour(kf_cs_1)] = cc;
      cc = new ConstituentCharacteristic(hadpars.Get("Mass_cc1"),0,0.,sm);
      CCMap[Flavour(kf_cc_1)] = cc;
      cc = new ConstituentCharacteristic(hadpars.Get("Mass_bd0"),0,0.,sm);
      CCMap[Flavour(kf_bd_0)] = cc;
      cc = new ConstituentCharacteristic(hadpars.Get("Mass_bu0"),0,0.,sm);
      CCMap[Flavour(kf_bu_0)] = cc;
      cc = new ConstituentCharacteristic(hadpars.Get("Mass_bs0"),0,0.,sm);
      CCMap[Flavour(kf_bs_0)] = cc;
      cc = new ConstituentCharacteristic(hadpars.Get("Mass_bc0"),0,0.,sm);
      CCMap[Flavour(kf_bc_0)] = cc;
      cc = new ConstituentCharacteristic(hadpars.Get("Mass_bd1"),0,0.,sm);
      CCMap[Flavour(kf_bd_1)] = cc;
      cc = new ConstituentCharacteristic(hadpars.Get("Mass_bu1"),0,0.,sm);
      CCMap[Flavour(kf_bu_1)] = cc;
      cc = new ConstituentCharacteristic(hadpars.Get("Mass_bs1"),0,0.,sm);
      CCMap[Flavour(kf_bs_1)] = cc;
      cc = new ConstituentCharacteristic(hadpars.Get("Mass_bc1"),0,0.,sm);
      CCMap[Flavour(kf_bc_1)] = cc;
      cc = new ConstituentCharacteristic(hadpars.Get("Mass_bb1"),0,0.,sm);
      CCMap[Flavour(kf_bb_1)] = cc;
    }
  }
  for (FlavCCMap_Iterator cmit=CCMap.begin(); cmit!=CCMap.end();cmit++) {
    if (cmit->first!=Flavour(kf_gluon) &&
	cmit->second->Mass()<m_minmass) m_minmass = cmit->second->Mass();
  }
}

Constituents::~Constituents() {
  for (FlavCCMap_Iterator cmit=CCMap.begin(); cmit!=CCMap.end();cmit++) {
    if (cmit->second) { delete cmit->second; cmit->second=NULL; }
  }
  CCMap.clear();
}

double Constituents::MinMass() { return m_minmass; }

double Constituents::Mass(const Flavour & flav) {
  FlavCCMap_Iterator     cmit = CCMap.find(flav);
  if (cmit==CCMap.end()) cmit = CCMap.find(flav.Bar());
  return (cmit!=CCMap.end())? cmit->second->Mass() : flav.HadMass();
}

double Constituents::TotWeight(const ATOOLS::Flavour & flav) {
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

