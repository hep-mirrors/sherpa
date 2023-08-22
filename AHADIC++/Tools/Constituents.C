#include "AHADIC++/Tools/Hadronisation_Parameters.H"
#include "AHADIC++/Tools/Constituents.H"
#include "ATOOLS/Org/Message.H"

using namespace AHADIC;
using namespace ATOOLS;

Constituents::Constituents(bool diquarks) :
  m_minmass(100.),m_maxmass(0.)
{
  // Light quarks and diquarks
  auto v_sfrac(hadpars->GetVec("Strange_fraction"));
  auto v_bfrac(hadpars->GetVec("Baryon_fraction"));
  auto v_qssup(hadpars->GetVec("P_qs_by_P_qq"));
  auto v_sssup(hadpars->GetVec("P_ss_by_P_qq"));
  auto v_sp1sup(hadpars->GetVec("P_di_1_by_P_di_0"));

  std::vector<Flavour> quarks = {
    Flavour(kf_d), Flavour(kf_u), Flavour(kf_s), Flavour(kf_c), Flavour(kf_b),
  };

  std::vector<Flavour> diquarks_spin_zero = {
    Flavour(kf_ud_0), Flavour(kf_sd_0), Flavour(kf_su_0),
  };

  std::vector<Flavour> diquarks_spin_one = {
    Flavour(kf_uu_1), Flavour(kf_ud_1), Flavour(kf_dd_1),
    Flavour(kf_su_1), Flavour(kf_sd_1), Flavour(kf_ss_1),
  };

  // init CCMap
  for (auto f : quarks) {
    CCMap[f] = new ConstituentCharacteristic(f.HadMass());
    CCMap[f]->m_ispin = 1;
    CCMap[f]->m_weight.resize(0);
  }

  for (auto f : diquarks_spin_zero) {
    CCMap[f] = new ConstituentCharacteristic(f.HadMass());
    CCMap[f]->m_ispin = 0;
    CCMap[f]->m_weight.resize(0);
  }

  for (auto f : diquarks_spin_one) {
    CCMap[f] = new ConstituentCharacteristic(f.HadMass());
    CCMap[f]->m_ispin = 2;
    CCMap[f]->m_weight.resize(0);
  }

  double total(0.),udfrac(1.), ud0(1.), norm;
  m_nvars = v_sfrac.size();
  for(int i{0}; i<v_sfrac.size(); ++i) {
    const double sfrac  = v_sfrac[i];
    const double bfrac  = v_bfrac[i];
    const double qssup  = v_qssup[i];
    const double sssup  = v_sssup[i];
    const double sp1sup = v_sp1sup[i];

    total  = 2.*(2.*udfrac+sfrac);
    total += bfrac*ud0*(1.+2.*qssup+3.*sp1sup*(3.+2.*qssup+sssup));
    norm = 1./total;

    CCMap[Flavour(kf_d)]->m_weight.push_back(2.*udfrac*norm);
    CCMap[Flavour(kf_u)]->m_weight.push_back(2.*udfrac*norm);
    CCMap[Flavour(kf_s)]->m_weight.push_back(2.*sfrac*norm);
    CCMap[Flavour(kf_c)]->m_weight.push_back(0.);
    CCMap[Flavour(kf_b)]->m_weight.push_back(0.);

    if (diquarks && bfrac>0.) {
      // Light Di-quarks, spin 0
      CCMap[Flavour(kf_ud_0)]->m_weight.push_back(bfrac*ud0*norm);
      CCMap[Flavour(kf_sd_0)]->m_weight.push_back(bfrac*qssup*norm);
      CCMap[Flavour(kf_su_0)]->m_weight.push_back(bfrac*qssup*norm);

      // Light Di-quarks, spin 1
      CCMap[Flavour(kf_uu_1)]->m_weight.push_back(3.*bfrac*sp1sup*norm);
      CCMap[Flavour(kf_ud_1)]->m_weight.push_back(3.*bfrac*sp1sup*norm);
      CCMap[Flavour(kf_dd_1)]->m_weight.push_back(3.*bfrac*sp1sup*norm);
      CCMap[Flavour(kf_su_1)]->m_weight.push_back(3.*bfrac*sp1sup*qssup*norm);
      CCMap[Flavour(kf_sd_1)]->m_weight.push_back(3.*bfrac*sp1sup*qssup*norm);
      CCMap[Flavour(kf_ss_1)]->m_weight.push_back(3.*bfrac*sp1sup*sssup*norm);
    }
  }

  //CCMap = FlavourMaps[0];
  double massoffset(hadpars->Get("minmass2"));
  for (FlavCCMap_Iterator cmit=CCMap.begin(); cmit!=CCMap.end();cmit++) {
    if (cmit->first==Flavour(kf_gluon)) continue;
    if (cmit->second->Mass()+massoffset<m_minmass) {
      m_minmass  = cmit->second->Mass()+massoffset;
      m_lightest = cmit->first;
    }
    if (cmit->second->Mass()+massoffset>m_maxmass)
      m_maxmass = cmit->second->Mass()+massoffset;
  }
}

Constituents::~Constituents() {
  for (FlavCCMap_Iterator cmit=CCMap.begin(); cmit!=CCMap.end();cmit++) {
    if (cmit->second) { delete cmit->second; cmit->second=NULL; }
  }
  CCMap.clear();
}

double Constituents::MinMass() { return m_minmass; }
double Constituents::MaxMass() { return m_maxmass; }

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

std::vector<double> Constituents::Weights(const ATOOLS::Flavour & flav) {
  FlavCCMap_Iterator     cmit = CCMap.find(flav);
  if (cmit==CCMap.end()) cmit = CCMap.find(flav.Bar());
  if (cmit!=CCMap.end())
    return cmit->second->m_weight;
  return {};
}

int Constituents::ISpin(Flavour & flav) {
  FlavCCMap_Iterator     cmit = CCMap.find(flav);
  if (cmit==CCMap.end()) cmit = CCMap.find(flav.Bar());
  return (cmit!=CCMap.end())? cmit->second->ISpin() : 0;
}

void Constituents::PrintConstituents() {
  double wt(0.),wtq(0.),wtd(0.);
  for (FlavCCMap_Iterator cmit=CCMap.begin();cmit!=CCMap.end();cmit++) {
    wt+=cmit->second->TotWeight();
    if (cmit->first.IsQuark()) wtq+=cmit->second->TotWeight();
    else wtd+=cmit->second->TotWeight();
    msg_Out()<<cmit->first<<" : "<<cmit->second->m_mass<<" GeV, "
	     <<"Spin = "<<double(cmit->second->m_ispin/2.)<<", "
	     <<"Weight = "<<cmit->second->TotWeight()<<", "
	     <<"VarWeight = "<<cmit->second->m_weight
	     <<std::endl;
  }
  msg_Out()<<"Total weight : "<<wt<<" (quarks = "<<wtq<<", diquarks = "<<wtd<<")."<<std::endl
	   <<"------------- END OF CONSTITUENTS ---------------"<<std::endl;
}

