#include "AHADIC++/Tools/Multiplet_Constructor.H"
#include "AHADIC++/Tools/Hadronisation_Parameters.H"
#include "ATOOLS/Org/Message.H"
#include <stdio.h>

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;

namespace AHADIC {
  ostream & operator<<(ostream & s,const HadInfo & info) 
  {
    s<<" "<<info.flav<<" ["<<info.multiname<<"]:"
     <<info.iso<<" ex_{r,l} = "<<info.exr<<","<<info.exl<<" [";
    if (info.fl3!=0) s<<info.fl3<<",";
    s<<info.fl2<<","<<info.fl1<<"], spin = "<<(info.spin2-1.)/2.;
    return s;
  }
}

// To test the wave function and multiplet construction, together
// with the transitions, set test=true.
// Default is test=false.

Multiplet_Constructor::Multiplet_Constructor(bool test) :
  m_test(test),
  m_singletsuppression(hadpars->Get("Singlet_Suppression"))
{
  CreateMultiplets();
}

Multiplet_Constructor::~Multiplet_Constructor() {
  while (!m_wavefunctions.empty()) {
    delete m_wavefunctions.begin()->second;
    m_wavefunctions.erase(m_wavefunctions.begin());
  }
  m_wavefunctions.clear();
}

void Multiplet_Constructor::CreateMultiplets() {
  // iterate over all hadrons in Hadron_Init:
  // - extract flavour, spin and excitation information,
  // - ignore for the moment the 9xxxxxx hadrons
  // - attach hard-coded quark-diquark wavefunctions
  // - add "anti-mesons" in same multiplet
  // - add anti-baryons in separate multiplets
  // - put the hadrons in multiplets.
  for(KFCode_ParticleInfo_Map::const_iterator kfit(s_kftable.begin());
      kfit!=s_kftable.end();++kfit) {
    if (!ExtractInfo(kfit->first) || m_info.iso>0 || m_info.multiwt<1.e-6 ||
	!ConstructWaveFunction()) continue;
    m_multiplets[m_info.multiname].insert(m_info.flav);
    if (m_info.barrable && !m_info.flav.IsBaryon()) {
      m_multiplets[m_info.multiname].insert(m_info.flav.Bar());
    }
    if (m_info.flav.IsBaryon()) {
      if (!ConstructAntiBaryonWaveFunction(m_info.flav)) {
	msg_Error()<<METHOD<<" throws error:\n"
		   <<"   Could not derieve anti-particle wave function for "
		   <<m_info.flav<<", will exit the run.\n";
	exit(1);
      }
      string antiname = string("Anti-")+m_info.multiname;
      m_multiplets[antiname].insert(m_info.flav.Bar());
    }
  }
}

bool Multiplet_Constructor::ExtractInfo(const kf_code & kfc) {
  Flavour hadron   = Flavour(kfc);
  if (!hadron.IsHadron() || !hadron.IsOn()) return false;
  m_info.flav  = hadron;
  int kf     = int(kfc);
  m_info.iso   = int(kf/9000000);
  if (m_info.iso>0) return false;
  kf        -= m_info.iso*9000000;
  m_info.exr   = int(kf/100000);
  kf        -= m_info.exr*100000;
  m_info.exl   = int(kf/10000);
  kf        -= m_info.exl*10000;
  m_info.fl3   = int(kf/1000);
  kf        -= m_info.fl3*1000;
  m_info.fl2   = int(kf/100);
  kf        -= m_info.fl2*100;
  m_info.fl1   = int(kf/10);
  kf        -= m_info.fl1*10;
  m_info.spin2 = int(kf);
  if (m_info.spin2<=0) return false;
  m_info.barrable  = m_info.flav.IsBaryon() || (m_info.fl2!=m_info.fl1);
  m_info.multiname = MultipletName();
  FillMultipletWeights();
  return (m_info.multiname!=string(""));
}

std::string Multiplet_Constructor::MultipletName() {
  string name = string("");
  if (m_info.exr!=0)   name += string("R=")+ToString(m_info.exr)+string("_");
  if (m_info.exl!=0)   name += string("L=")+ToString(m_info.exl)+string("_");
  if (m_info.fl3==0) {
    if (m_info.spin2==1) name += "Scalars";
    if (m_info.spin2==3) name += "Vectors";
    if (m_info.spin2==5) name += "Tensors";
  }
  else {
    if (m_info.spin2==2) name += "Octet";
    if (m_info.spin2==4) {
      if (m_info.exr==0) name += "Decuplet";
      if (m_info.exr==1) name += "1_Octet";
    }
  }
  return name;
}

void Multiplet_Constructor::FillMultipletWeights() {
  m_info.multiwt = 0.;
  if (m_info.multiname==string("Scalars"))
    m_info.multiwt = hadpars->Get("Multiplet_Meson_R0L0S0");
  if (m_info.multiname==string("Vectors"))
    m_info.multiwt = hadpars->Get("Multiplet_Meson_R0L0S1");
  if (m_info.multiname==string("Tensors"))
    m_info.multiwt = hadpars->Get("Multiplet_Meson_R0L0S2");
  if (m_info.multiname==string("L=1_Scalars"))
    m_info.multiwt = hadpars->Get("Multiplet_Meson_R0L1S0");
  if (m_info.multiname==string("L=1_Vectors"))
    m_info.multiwt = hadpars->Get("Multiplet_Meson_R0L1S1");
  if (m_info.multiname==string("L=2_Vectors"))
    m_info.multiwt = hadpars->Get("Multiplet_Meson_R0L2S2");
  if (m_info.multiname==string("Octet"))
    m_info.multiwt = hadpars->Get("Multiplet_Baryon_R0L0S1/2");
  if (m_info.multiname==string("Decuplet"))
    m_info.multiwt = hadpars->Get("Multiplet_Baryon_R0L0S3/2");
  if (m_info.multiname==string("R=1_Octet"))
    m_info.multiwt = hadpars->Get("Multiplet_Baryon_R1L0S1/2");
  if (m_info.multiname==string("R=1_1_Octet"))
    m_info.multiwt = hadpars->Get("Multiplet_Baryon_R1_1L0S1/2");
  if (m_info.multiname==string("R=2_Octet"))
    m_info.multiwt = hadpars->Get("Multiplet_Baryon_R0L0S1/2");
  m_info.spinwt  = double(m_info.spin2);
  m_info.extrawt = 1.;
  if (m_info.flav.IsBaryon() && m_info.fl3>3)
    m_info.extrawt *= hadpars->Get("Heavy_Baryon_Enhancement");
}

bool Multiplet_Constructor::ConstructWaveFunction() 
{
  bool constructed = (m_info.fl3==0?
		      ConstructMesonWaveFunction():
		      ConstructBaryonWaveFunction());
  if (m_test && constructed) return true;
  if (constructed && m_wavefunctions.find(m_info.flav)!=m_wavefunctions.end()) {
    m_wavefunctions[m_info.flav]->SetMultipletWeight(m_info.multiwt);
    m_wavefunctions[m_info.flav]->SetSpin(m_info.spin2);
    m_wavefunctions[m_info.flav]->SetExtraWeight(m_info.extrawt);
    if (m_info.fl3==0 && m_info.barrable) {
      m_wavefunctions[m_info.flav.Bar()]->SetMultipletWeight(m_info.multiwt);
      m_wavefunctions[m_info.flav.Bar()]->SetSpin(m_info.spin2);
      m_wavefunctions[m_info.flav.Bar()]->SetExtraWeight(m_info.extrawt);
    }
    return true;
  }
  return false;
}

bool Multiplet_Constructor::ConstructMesonWaveFunction() 
{
  // these are the "funny mesons" ... a0(980) and friends ...
  // no idea (yet) how to deal with them.
  if (m_info.iso>0) return false;
  
  if ((m_info.fl1!=m_info.fl2) ||
      (m_info.fl1==m_info.fl2 && (m_info.fl1==4 || m_info.fl1==5)) ||
      (m_info.fl1==4 && m_info.fl2==5)) {
    m_wavefunctions[m_info.flav] = TrivialMesonWaveFunction();
    if (m_info.barrable) {
      m_wavefunctions[m_info.flav.Bar()] =
	m_wavefunctions[m_info.flav]->GetAnti();
    }
  }
  else if (m_info.fl1==m_info.fl2 && m_info.fl1==1) 
    m_wavefunctions[m_info.flav] = Pi0WaveFunction();
  else if ((m_info.fl1==m_info.fl2 && m_info.fl1==2 && m_info.spin2==1) ||
	   (m_info.fl1==m_info.fl2 && m_info.fl1==3 && m_info.spin2!=1)) {
    double theta   = MixingAngle(), costh = cos(theta), sinth = sin(theta);
    m_info.extrawt = costh*costh+sinth*sinth*m_singletsuppression; 
    m_wavefunctions[m_info.flav] = OctetMesonWaveFunction();
  }
  else if ((m_info.fl1==m_info.fl2 && m_info.fl1==3 && m_info.spin2==1) ||
	   (m_info.fl1==m_info.fl2 && m_info.fl1==2 && m_info.spin2!=1)) {
    double theta   = MixingAngle(), costh = cos(theta), sinth = sin(theta);
    m_info.extrawt = costh*costh*m_singletsuppression+sinth*sinth; 
    m_wavefunctions[m_info.flav] = SingletMesonWaveFunction();
  }
  return (m_wavefunctions.find(m_info.flav)!=m_wavefunctions.end());
}

Wave_Function * Multiplet_Constructor::TrivialMesonWaveFunction() {
  bool even(m_info.fl2/2==m_info.fl2/2.);
  Flavour fl1   = Flavour((kf_code)(m_info.fl1));
  Flavour fl2   = Flavour((kf_code)(m_info.fl2));
  Flavour_Pair  * pair         = new Flavour_Pair;
  pair->first   = even?fl2:fl1;
  pair->second  = even?fl1.Bar():fl2.Bar();
  Wave_Function * wavefunction = new Wave_Function(m_info.flav.Bar());
  wavefunction->AddToWaves(pair,1.);
  return wavefunction;
}

Wave_Function * Multiplet_Constructor::Pi0WaveFunction() {
  // Wave function for neutral state with isospin = strangeness = 0
  // Essentially pi^0, rho^0 ... .  Fixed by isospin.
  Wave_Function * wavefunction = new Wave_Function(m_info.flav);
  Flavour_Pair  * pair = new Flavour_Pair;
  pair->first  = Flavour(kf_d);
  pair->second = Flavour(kf_d).Bar();
  wavefunction->AddToWaves(pair,-1./sqrt(2.));
  pair = new Flavour_Pair;
  pair->first  = Flavour(kf_u);
  pair->second = Flavour(kf_u).Bar();
  wavefunction->AddToWaves(pair,+1./sqrt(2.));
  return wavefunction;
}

Wave_Function * Multiplet_Constructor::OctetMesonWaveFunction() {
  // Trivially, for mixing angle = 0, this is the octet state, i.e.
  // 1/sqrt(6) [d dbar + u ubar - 2 s sbar]
  // Unfortunately for all but the pseudoscalars this is the heavier state
  double theta  = MixingAngle(), sinth=sin(theta), costh=cos(theta);
  double weight = costh/sqrt(6.)-sinth/sqrt(3.);
  Wave_Function * wavefunction = new Wave_Function(m_info.flav);
  Flavour_Pair  * pair;
  if (dabs(weight)>1.e-3) {
    pair = new Flavour_Pair;
    pair->first  = Flavour(kf_d);
    pair->second = Flavour(kf_d).Bar();
    wavefunction->AddToWaves(pair,weight);
    pair = new Flavour_Pair;
    pair->first  = Flavour(kf_u);
    pair->second = Flavour(kf_u).Bar();
    wavefunction->AddToWaves(pair,weight);
  }
  weight         = -2.*costh/sqrt(6.)-sinth/sqrt(3.);
  if (dabs(weight)>1.e-3) {
    pair = new Flavour_Pair;
    pair->first  = Flavour(kf_s);
    pair->second = Flavour(kf_s).Bar();
    wavefunction->AddToWaves(pair,weight);
  }
  return wavefunction;
}

Wave_Function * Multiplet_Constructor::SingletMesonWaveFunction() {
  // Trivially, for mixing angle = 0, this is the singlet state, i.e.
  // 1/sqrt(3) [d dbar + u ubar + s sbar], up to a phase.
  double theta  = MixingAngle(), sinth = sin(theta), costh = cos(theta);
  double weight = sinth/sqrt(6.)+costh/sqrt(3.);
  Wave_Function * wavefunction = new Wave_Function(m_info.flav);
  Flavour_Pair  * pair;
  if (dabs(weight)>1.e-3) {
    pair = new Flavour_Pair;
    pair->first  = Flavour(kf_d);
    pair->second = Flavour(kf_d).Bar();
    wavefunction->AddToWaves(pair,weight);
    pair = new Flavour_Pair;
    pair->first  = Flavour(kf_u);
    pair->second = Flavour(kf_u).Bar();
    wavefunction->AddToWaves(pair,weight);
  }
  weight  = -2.*sinth/sqrt(6.)+costh/sqrt(3.);
  if (dabs(weight)>1.e-3) {
    pair = new Flavour_Pair;
    pair->first  = Flavour(kf_s);
    pair->second = Flavour(kf_s).Bar();
    wavefunction->AddToWaves(pair,weight);
  }
  return wavefunction;
}

double Multiplet_Constructor::MixingAngle()
{
  switch (m_info.spin2) {
  case 5 : return hadpars->Get("Mixing_Angle_2+");
  case 3 : return hadpars->Get("Mixing_Angle_1-");
  case 1 : return hadpars->Get("Mixing_Angle_0+");
  default: break;  
  }
  return 0.;
}

bool Multiplet_Constructor::ConstructBaryonWaveFunction() 
{
  // SU(3) baryon wave functions according to
  // Lichtenberg, Namgung, Wills & Predazzi
  //
  // Octet like -- there are two kinds of multiplets here:
  // the 56-plet and the 70-plet: the former have one Lambda
  // in SU(3), the latter have two Lambda's, the extra on
  // being a singlet, with no or little mising with the
  // "normal" Sigma and Lambda - quite often the heavies
  // are unknown
  if (m_info.spin2==2 || (m_info.spin2==4 && m_info.exr==1)) {
    if (m_info.fl3<4) {
      if (m_info.fl3>m_info.fl2 && m_info.fl3>m_info.fl1) {
	if (m_info.fl2==m_info.fl1)
	  m_wavefunctions[m_info.flav] = NeutronWaveFunction();
	else if (m_info.fl2>m_info.fl1)
	  m_wavefunctions[m_info.flav] = SigmaWaveFunction();
	else if (m_info.fl2<m_info.fl1 && m_info.fl1<4)
	  m_wavefunctions[m_info.flav] = LambdaWaveFunction();
	else if (m_info.fl2<m_info.fl1 && m_info.fl1>3) {
	  // heavy singlet Lambda - have to exchnage positions 1 and 3
	  int help   = m_info.fl1;
	  m_info.fl1 = m_info.fl3;
	  m_info.fl3 = help;
	  m_wavefunctions[m_info.flav] = LambdaHWaveFunction();
	}
      }
      else if (m_info.fl3==m_info.fl2 && m_info.fl2>m_info.fl1) 
	m_wavefunctions[m_info.flav] = ProtonWaveFunction();
      else if (m_info.fl3==2 && m_info.fl2==1 && m_info.fl1==3)
	m_wavefunctions[m_info.flav] = Lambda1WaveFunction();
    }
    else if (m_info.fl3>3 && m_info.fl2<4 && m_info.fl1<4) {
      if (m_info.fl2>=m_info.fl1)
	m_wavefunctions[m_info.flav] = SigmaHWaveFunction();
      else if (m_info.fl2<m_info.fl1)
	m_wavefunctions[m_info.flav] = LambdaHWaveFunction();
    }
    return true;
  }
  // Decuplet-like - the heavies practically are unknown
  else if (m_info.spin2==4) {
    if (m_info.fl3<4) {
      if (m_info.fl3==m_info.fl2 && m_info.fl2==m_info.fl1) {
	m_wavefunctions[m_info.flav] = DeltaPPWaveFunction();
      }
      if (m_info.fl3==m_info.fl2 && m_info.fl2>m_info.fl1) {
	m_wavefunctions[m_info.flav] = DeltaPWaveFunction();
      }
      if (m_info.fl3>m_info.fl2 && m_info.fl2==m_info.fl1) {
	m_wavefunctions[m_info.flav] = Delta0WaveFunction();
      }
      if (m_info.fl3>m_info.fl2 && m_info.fl2>m_info.fl1) {
	m_wavefunctions[m_info.flav] = DecupletSigmaWaveFunction();
      }
    }
    else if (m_info.fl3>3 && m_info.fl2<4 && m_info.fl1<4) {
      m_wavefunctions[m_info.flav] = HeavyDecupletWaveFunction();
    }
    return true;
  }
  msg_Error()<<"Error in "<<METHOD<<":\n"
	     <<"   No wavefunction (yet) for "<<m_info.flav<<".\n";
  return false;
}

bool Multiplet_Constructor::
ConstructAntiBaryonWaveFunction(Flavour & flav) {
  m_wavefunctions[flav.Bar()] = m_wavefunctions[flav]->GetAnti();
  return (m_wavefunctions[flav.Bar()]!=NULL);
}

Wave_Function * Multiplet_Constructor::NeutronWaveFunction() {
  // kfcode = 2112
  // 1/sqrt{3}[u+ (dd)_1] + 1/sqrt{6}[d+(ud)_1] + 1/\sqrt{2}[d+(ud)_0]
  Wave_Function * wavefunction = new Wave_Function(m_info.flav);
  Flavour_Pair * pair = new Flavour_Pair;
  pair->first  = Flavour((kf_code)(m_info.fl3));
  pair->second = Flavour((kf_code)(m_info.fl1*1000+m_info.fl2*100+3));
  wavefunction->AddToWaves(pair,+1./sqrt(3.));
  pair         = new Flavour_Pair;
  pair->first  = Flavour((kf_code)(m_info.fl2));
  pair->second = Flavour((kf_code)(m_info.fl3*1000+m_info.fl1*100+3)); 
  wavefunction->AddToWaves(pair,+1./sqrt(6.));
  pair         = new Flavour_Pair;
  pair->first  = Flavour((kf_code)(m_info.fl2));
  pair->second = Flavour((kf_code)(m_info.fl3*1000+m_info.fl1*100+1)); 
  wavefunction->AddToWaves(pair,+1./sqrt(2.));
  return wavefunction;
}

Wave_Function * Multiplet_Constructor::ProtonWaveFunction() {
  // kfcode = 2212
  // 1/sqrt{3}[d+ (uu)_1] + 1/sqrt{6}[u+(ud)_1] + 1/\sqrt{2}[u+(ud)_0]
  Wave_Function * wavefunction = new Wave_Function(m_info.flav);
  Flavour_Pair * pair = new Flavour_Pair;
  pair->first  = Flavour((kf_code)(m_info.fl1));
  pair->second = Flavour((kf_code)(m_info.fl2*1000+m_info.fl3*100+3));
  wavefunction->AddToWaves(pair,+1./sqrt(3.));
  pair         = new Flavour_Pair;
  pair->first  = Flavour((kf_code)(m_info.fl2));
  pair->second = Flavour((kf_code)(m_info.fl3*1000+m_info.fl1*100+3)); 
  wavefunction->AddToWaves(pair,+1./sqrt(6.));
  pair         = new Flavour_Pair;
  pair->first  = Flavour((kf_code)(m_info.fl2));
  pair->second = Flavour((kf_code)(m_info.fl3*1000+m_info.fl1*100+1));  
  wavefunction->AddToWaves(pair,+1./sqrt(2.));
  return wavefunction;
}

Wave_Function * Multiplet_Constructor::SigmaWaveFunction() {
  // kfcode = 3212
  // 1/sqrt{3}[s+ (ud)_1] + 1/sqrt{12}[d+(su)_1] + 1/\sqrt{4}[d+(su)_0]
  //                      + 1/sqrt{12}[u+(sd)_1] + 1/\sqrt{4}[u+(sd)_0]
  Wave_Function * wavefunction = new Wave_Function(m_info.flav);
  Flavour_Pair * pair = new Flavour_Pair;
  pair->first  = Flavour((kf_code)(m_info.fl3));
  pair->second = Flavour((kf_code)(m_info.fl2*1000+m_info.fl1*100+3));
  wavefunction->AddToWaves(pair,+1./sqrt(3.));
  pair         = new Flavour_Pair;
  pair->first  = Flavour((kf_code)(m_info.fl1));
  pair->second = Flavour((kf_code)(m_info.fl3*1000+m_info.fl2*100+3));
  wavefunction->AddToWaves(pair,+1./sqrt(12.));
  pair         = new Flavour_Pair;
  pair->first  = Flavour((kf_code)(m_info.fl1));
  pair->second = Flavour((kf_code)(m_info.fl3*1000+m_info.fl2*100+1));  
  wavefunction->AddToWaves(pair,+1./sqrt(4.));
  pair         = new Flavour_Pair;
  pair->first  = Flavour((kf_code)(m_info.fl2));
  pair->second = Flavour((kf_code)(m_info.fl3*1000+m_info.fl1*100+3)); 
  wavefunction->AddToWaves(pair,+1./sqrt(12.));
  pair         = new Flavour_Pair;
  pair->first  = Flavour((kf_code)(m_info.fl2));
  pair->second = Flavour((kf_code)(m_info.fl3*1000+m_info.fl1*100+1)); 
  wavefunction->AddToWaves(pair,+1./sqrt(4.));
  return wavefunction;
}

Wave_Function * Multiplet_Constructor::LambdaWaveFunction() {
  // kfcode = 3122
  // 1/sqrt{3}[s+ (ud)_0] + 1/sqrt{12}[d+(su)_0] + 1/\sqrt{4}[d+(su)_1]
  //                      + 1/sqrt{12}[u+(sd)_0] + 1/\sqrt{4}[u+(sd)_1]
  Wave_Function * wavefunction = new Wave_Function(m_info.flav);
  Flavour_Pair * pair = new Flavour_Pair;
  pair->first  = Flavour((kf_code)(m_info.fl3));
  pair->second = Flavour((kf_code)(m_info.fl1*1000+m_info.fl2*100+1));
  wavefunction->AddToWaves(pair,+1./sqrt(3.));
  pair         = new Flavour_Pair;
  pair->first  = Flavour((kf_code)(m_info.fl1));
  pair->second = Flavour((kf_code)(m_info.fl3*1000+m_info.fl2*100+1));
  wavefunction->AddToWaves(pair,+1./sqrt(12.));
  pair         = new Flavour_Pair;
  pair->first  = Flavour((kf_code)(m_info.fl1));
  pair->second = Flavour((kf_code)(m_info.fl3*1000+m_info.fl2*100+3));
  wavefunction->AddToWaves(pair,+1./sqrt(4.));
  pair         = new Flavour_Pair;
  pair->first  = Flavour((kf_code)(m_info.fl2));
  pair->second = Flavour((kf_code)(m_info.fl3*1000+m_info.fl1*100+1)); 
  wavefunction->AddToWaves(pair,+1./sqrt(12.));
  pair         = new Flavour_Pair;
  pair->first  = Flavour((kf_code)(m_info.fl2));
  pair->second = Flavour((kf_code)(m_info.fl3*1000+m_info.fl1*100+3)); 
  wavefunction->AddToWaves(pair,+1./sqrt(4.));
  return wavefunction;
}

Wave_Function * Multiplet_Constructor::Lambda1WaveFunction() {
  // kfcode = 2132
  // 1/sqrt{3}[s+ (ud)_0] + 1/sqrt{3}[d+(su)_0] + 1/sqrt{3}[u+(sd)_0]
  Wave_Function * wavefunction = new Wave_Function(m_info.flav);
  Flavour_Pair * pair = new Flavour_Pair;
  pair->first  = Flavour(kf_d);
  pair->second = Flavour(kf_su_0); 
  wavefunction->AddToWaves(pair,+1./sqrt(3.));
  pair = new Flavour_Pair;
  pair->first  = Flavour(kf_u);
  pair->second = Flavour(kf_sd_0); 
  wavefunction->AddToWaves(pair,+1./sqrt(3.));
  pair = new Flavour_Pair;
  pair->first  = Flavour(kf_s);
  pair->second = Flavour(kf_ud_0); 
  wavefunction->AddToWaves(pair,+1./sqrt(3.));
  return wavefunction;
}

Wave_Function * Multiplet_Constructor::SigmaHWaveFunction() {
  // kfcode = 4212 or 4312 or 4322
  // [c+ (ud)_1] or [c + (su)_1] or [c + (sd)_1] (wild guess ...)
  Wave_Function * wavefunction = new Wave_Function(m_info.flav);
  Flavour_Pair * pair = new Flavour_Pair;
  pair->first  = Flavour((kf_code)(m_info.fl3));
  pair->second = Flavour((kf_code)(m_info.fl2*1000+m_info.fl1*100+3));
  wavefunction->AddToWaves(pair,1.);
  return wavefunction;
}

Wave_Function * Multiplet_Constructor::LambdaHWaveFunction() {
  // kfcode = 4122 or 4132 or 4232
  // [c+ (ud)_0] or [c + (su)_0] or [c + (sd)_0]
  Wave_Function * wavefunction = new Wave_Function(m_info.flav);
  Flavour_Pair * pair = new Flavour_Pair;
  pair->first  = Flavour((kf_code)(m_info.fl3));
  pair->second = Flavour((kf_code)(m_info.fl1*1000+m_info.fl2*100+1));
  wavefunction->AddToWaves(pair,1.);
  return wavefunction;
}

Wave_Function * Multiplet_Constructor::DeltaPPWaveFunction() {
  // kfcode = 1114, 2224, or 3334
  // [q (qq)_1]
  Wave_Function * wavefunction = new Wave_Function(m_info.flav);
  Flavour_Pair * pair = new Flavour_Pair;
  pair->first  = Flavour((kf_code)(m_info.fl3));
  pair->second = Flavour((kf_code)(m_info.fl3*1100+3));
  wavefunction->AddToWaves(pair,1.);
  return wavefunction;
}

Wave_Function * Multiplet_Constructor::DeltaPWaveFunction() {
  // kfcode = 2214, 3314, or 3324
  // [q (qq)_1]
  Wave_Function * wavefunction = new Wave_Function(m_info.flav);
  Flavour_Pair * pair = new Flavour_Pair;
  pair->first  = Flavour((kf_code)(m_info.fl3));
  pair->second = Flavour((kf_code)(m_info.fl2*1000+m_info.fl1*100+3));
  wavefunction->AddToWaves(pair,sqrt(2./3.));
  pair = new Flavour_Pair;
  pair->first  = Flavour((kf_code)(m_info.fl1));
  pair->second = Flavour((kf_code)(m_info.fl3*1100+3));
  wavefunction->AddToWaves(pair,1./sqrt(3.));
  return wavefunction;
}

Wave_Function * Multiplet_Constructor::Delta0WaveFunction() {
  // kfcode = 2114, 3114, or 3224
  // [q (qq)_1]
  Wave_Function * wavefunction = new Wave_Function(m_info.flav);
  Flavour_Pair * pair = new Flavour_Pair;
  pair->first  = Flavour((kf_code)(m_info.fl3));
  pair->second = Flavour((kf_code)(m_info.fl2*1100+3));
  wavefunction->AddToWaves(pair,1./sqrt(3.));
  pair = new Flavour_Pair;
  pair->first  = Flavour((kf_code)(m_info.fl1));
  pair->second = Flavour((kf_code)(m_info.fl3*1000+m_info.fl2*100+3));
  wavefunction->AddToWaves(pair,sqrt(2./3.));
  return wavefunction;
}

Wave_Function * Multiplet_Constructor::DecupletSigmaWaveFunction() {
  // kfcode = 3214
  // 1/sqrt{3}[s+ (ud)_1] + 1/sqrt{3}[d+(su)_1] + 1/sqrt{3}[u+(sd)_1]
  Wave_Function * wavefunction = new Wave_Function(m_info.flav);
  Flavour_Pair * pair = new Flavour_Pair;
  pair->first  = Flavour((kf_code)(m_info.fl3));
  pair->second = Flavour((kf_code)(m_info.fl2*1000+m_info.fl1*100+3));
  wavefunction->AddToWaves(pair,+1./sqrt(3.));
  pair = new Flavour_Pair;
  pair->first  = Flavour((kf_code)(m_info.fl2));
  pair->second = Flavour((kf_code)(m_info.fl3*1000+m_info.fl1*100+3));
  wavefunction->AddToWaves(pair,+1./sqrt(3.));
  pair = new Flavour_Pair;
  pair->first  = Flavour((kf_code)(m_info.fl1));
  pair->second = Flavour((kf_code)(m_info.fl3*1000+m_info.fl2*100+3));
  wavefunction->AddToWaves(pair,+1./sqrt(3.));
  return wavefunction;
}

Wave_Function * Multiplet_Constructor::HeavyDecupletWaveFunction() {
  // kfcode = 4114, 4214, 4224, 4314, 4324, 4334, or similar for beauties
  // [c+ (ud)_1] or [c + (su)_1] or [c + (sd)_1] 
  Wave_Function * wavefunction = new Wave_Function(m_info.flav);
  Flavour_Pair * pair = new Flavour_Pair;
  pair->first  = Flavour((kf_code)(m_info.fl3));
  pair->second = Flavour((kf_code)(m_info.fl2*1000+m_info.fl1*100+3));
  wavefunction->AddToWaves(pair,1.);
  return wavefunction;
}



void Multiplet_Constructor::PrintWaveFunctions(bool checkonly) 
{
  map<Flavour,double> checkit;
  for (map<string,set<Flavour> >::iterator mplet=m_multiplets.begin();
       mplet!=m_multiplets.end();mplet++) {
    checkit.clear();
    msg_Out()<<"-----------------------------------------------\n"
	     <<" "<<mplet->first<<" with "
	     <<mplet->second.size()<<" elements: "<<endl;
    for (set<Flavour>::iterator flit=mplet->second.begin();
	 flit!=mplet->second.end();flit++) {
      if (m_wavefunctions.find((*flit))==m_wavefunctions.end()) {
	msg_Out()<<"   no wave functions found for "<<(*flit)
		 <<" ["<<flit->Kfcode()<<"].\n";
      }
      Wave_Function * wavefunction = m_wavefunctions.find(*flit)->second;
      if (!checkonly) msg_Out()<<(*wavefunction);
      for (WaveComponents::iterator wit=wavefunction->GetWaves()->begin();
	   wit!=wavefunction->GetWaves()->end();wit++) {
	if (checkit.find(wit->first->first)!=checkit.end())
	  checkit[wit->first->first] += sqr(wit->second);
	else
	  checkit[wit->first->first] = sqr(wit->second);
	if (checkit.find(wit->first->second)!=checkit.end())
	  checkit[wit->first->second] += sqr(wit->second);
	else
	  checkit[wit->first->second] = sqr(wit->second);
      }
    }
    msg_Out()<<"-----------------------------------------------\n"
	     <<"-- Total weights for the constituents: --------\n"
	     <<"-----------------------------------------------\n";
    for (map<Flavour,double>::iterator cit=checkit.begin();
	 cit!=checkit.end();cit++)
      msg_Out()<<"   "<<cit->first<<" : "<<cit->second<<"\n";
    msg_Out()<<"-----------------------------------------------\n"
	     <<"-----------------------------------------------\n";
  }
}


void Multiplet_Constructor::PrintMultiplets() 
{
  msg_Out()<<"**********************************************************\n";
  for (map<string,set<Flavour> >::iterator multi=m_multiplets.begin();
       multi!=m_multiplets.end();multi++) {
    msg_Out()<<"*** "<<multi->first<<"  ["<<multi->second.size()<<"]:\n   ";
    for (set<Flavour>::iterator flit=multi->second.begin();
	 flit!=multi->second.end();flit++)
      msg_Out()<<(*flit)<<" ";
    msg_Out()<<"\n\n";
  }
}
