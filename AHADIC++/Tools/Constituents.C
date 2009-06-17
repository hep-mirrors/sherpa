#include "AHADIC++/Tools/Hadronisation_Parameters.H"
#include "AHADIC++/Tools/Constituents.H"
#include "ATOOLS/Org/Message.H"

using namespace AHADIC;
using namespace ATOOLS;

Constituents::Constituents(bool no_diquarks) :
  m_minmass(100.)
{
  // Light quarks and diquarks
  double total(0.),udfrac(1.), ud0(1.);
  double sfrac(hadpars.Get("Strange_fraction")); 
  double bfrac(hadpars.Get("Baryon_fraction"));
  double qssup(hadpars.Get("P_qs_by_P_qq"));
  double sssup(hadpars.Get("P_ss_by_P_qq"));
  double sp1sup(hadpars.Get("P_di_1_by_P_di_0"));

  total  = 2.*(2.*udfrac+sfrac);
  total += bfrac*ud0*(1.+2.*qssup+3.*sp1sup*(3.+2.*qssup+sssup));
  double norm = 1./total;

  Flavour flav(Flavour(kf_d));
  CCMap[flav] = new ConstituentCharacteristic(flav.HadMass(),1,udfrac*norm);
  flav = Flavour(kf_u);
  CCMap[flav] = new ConstituentCharacteristic(flav.HadMass(),1,udfrac*norm);
  flav = Flavour(kf_s);
  CCMap[flav] = new ConstituentCharacteristic(flav.HadMass(),1,sfrac*norm);
  flav = Flavour(kf_c);
  CCMap[flav] = new ConstituentCharacteristic(flav.HadMass(),1,0.);
  flav = Flavour(kf_b);
  CCMap[flav] = new ConstituentCharacteristic(flav.HadMass(),1,0.);


  if (!no_diquarks && bfrac>0.) {
    // Light Di-quarks, spin 0
    flav = Flavour(kf_ud_0);
    CCMap[flav] = new ConstituentCharacteristic(flav.HadMass(),0,ud0*norm);
    flav = Flavour(kf_sd_0);
    CCMap[flav] = new ConstituentCharacteristic(flav.HadMass(),0,qssup*norm);
    flav = Flavour(kf_su_0);
    CCMap[flav] = new ConstituentCharacteristic(flav.HadMass(),0,qssup*norm);

    // Light Di-quarks, spin 1
    flav = Flavour(kf_uu_1);
    CCMap[flav] = new ConstituentCharacteristic(flav.HadMass(),2,3.*sp1sup*norm);
    flav = Flavour(kf_ud_1);
    CCMap[flav] = new ConstituentCharacteristic(flav.HadMass(),2,3.*sp1sup*norm);
    flav = Flavour(kf_dd_1);
    CCMap[flav] = new ConstituentCharacteristic(flav.HadMass(),2,3.*sp1sup*norm);
    flav = Flavour(kf_su_1);
    CCMap[flav] = new ConstituentCharacteristic(flav.HadMass(),2,3.*sp1sup*qssup*norm);
    flav = Flavour(kf_sd_1);
    CCMap[flav] = new ConstituentCharacteristic(flav.HadMass(),2,3.*sp1sup*qssup*norm);
    flav = Flavour(kf_ss_1);
    CCMap[flav] = new ConstituentCharacteristic(flav.HadMass(),2,3.*sp1sup*sssup*norm);

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

int Constituents::ISpin(Flavour & flav) {
  FlavCCMap_Iterator     cmit = CCMap.find(flav);
  if (cmit==CCMap.end()) cmit = CCMap.find(flav.Bar());
  return (cmit!=CCMap.end())? cmit->second->ISpin() : 0;
}

void Constituents::PrintConstituents() {
  double wt(0.),wtq(0.),wtd(0.);
  for (FlavCCMap_Iterator cmit=CCMap.begin();cmit!=CCMap.end();cmit++) {
    wt+=cmit->second->m_weight;
    if (cmit->first.IsQuark()) wtq+=cmit->second->m_weight;
    else wtd+=cmit->second->m_weight;
    msg_Out()<<cmit->first<<" : "<<cmit->second->m_mass<<" GeV, "
	     <<"Spin = "<<double(cmit->second->m_ispin/2.)<<", "
	     <<"Weight = "<<cmit->second->m_weight<<std::endl;
  }
  msg_Out()<<"Total weight : "<<wt<<" (quarks = "<<wtq<<", diquarks = "<<wtd<<")."<<std::endl
	   <<"------------- END OF CONSTITUENTS ---------------"<<std::endl;
}

