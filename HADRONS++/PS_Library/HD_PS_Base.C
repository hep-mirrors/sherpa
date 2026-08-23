// ============================================================
// PROPOSED PATCH to HD_PS_Base.C
// ============================================================
//
// Root cause: the "Dalitz" (nout==3) branch has an explicit
// string->kf_code lookup table and derives BOTH the settings key
// (Flavour(kfres).IDName()) AND the default mass/width
// (Flavour(kfres).HadMass()/.Width()) from the RESOLVED kf-code - so
// an unrecognised or newly-introduced resonance name still gets a
// physically sensible default (whatever kfres falls back to), and a
// recognised one gets the CORRECT resonance's mass/width automatically,
// with no YAML override required.
//
// The "TwoResonances" (nout==4) branch has no such lookup: it builds
// SimpleResonanceFlavour directly from the raw parsed name string, and
// defaults SILENTLY to a1(1260)'s mass/width (res_a) or rho(770)'s
// mass/width (res_v) for ANY name not explicitly overridden via
// "Mass_<name>"/"Width_<name>" in Settings - regardless of what
// resonance was actually requested. This is the logical inconsistency:
// every genuinely new resonance name used in a TwoResonances phase-
// space channel (K(1)(1270)+, K(1)(1400)+, K*(1410)+, rho(1450)+, ...)
// silently gets the WRONG mass/width unless a matching override is
// supplied by hand, with no error or warning either way.
//
// Fix: factor the name->kf_code lookup out into one shared helper,
// used by BOTH branches, and extended to cover the resonance names
// actually in use. This keeps the "recognised name -> correct default
// from Sherpa's own particle data" property Dalitz already has, and
// gives TwoResonances the same property instead of requiring manual
// YAML overrides for every resonance.
// ============================================================

#include "HADRONS++/PS_Library/HD_PS_Base.H"
#include "HADRONS++/PS_Library/Two_Body_PSs.H"
#include "HADRONS++/PS_Library/Three_Body_PSs.H"
#include "HADRONS++/PS_Library/Four_Body_PSs.H"
#include "HADRONS++/PS_Library/Five_Body_PSs.H"
#include "PHASIC++/Channels/Rambo.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "HADRONS++/PS_Library/ResonanceFlavour.H"

using namespace HADRONS;
using namespace PHASIC;
using namespace ATOOLS;
using namespace std;


// ------------------------------------------------------------------
// Shared resonance-name -> kf_code lookup, used by BOTH the Dalitz
// (single-resonance) and TwoResonances (two-resonance) phase-space
// channel constructors below. Extend this table - not either
// constructor separately - whenever a NEW resonance name is used in a
// phase-space channel string anywhere in the decay YAML; this keeps
// both code paths using exactly the same mapping and avoids the two
// silently drifting apart again.
//
// Falls back to the given default kf-code (with a warning) for any
// unrecognised name, rather than silently defaulting with no
// indication at all - this makes a genuinely NEW/misspelled resonance
// name visible in the log instead of quietly sampling the wrong peak.
// ------------------------------------------------------------------
static kf_code ResonanceNameToKfCode(const std::string & name,
				      kf_code fallback) {
  if (name=="photon")        return kf_photon;
  if (name=="rho(770)+")     return kf_rho_770_plus;
  if (name=="rho(1450)+")    return kf_rho_1450_plus;
  if (name=="rho(1700)+")    return kf_rho_1700_plus;
  if (name=="K*(892)+")      return kf_K_star_892_plus;
  if (name=="K*(1410)+")     return kf_K_star_1410_plus;
  if (name=="K*(1680)+")     return kf_K_star_1680_plus;
  if (name=="a(1)(1260)+")   return kf_a_1_1260_plus;
  if (name=="K(1)(1270)+")   return kf_K_1_1270_plus;
  if (name=="K(1)(1400)+")   return kf_K_1_1400_plus;
  // Added on request: f0(600)/sigma. Confirmed real, pre-registered
  // kf-code (9000221, matching the general "9000xxx = non-qq-bar/
  // exotic scalar nonet" PDG MC-ID pattern also used for a0(980)/
  // a0(1450) below) - NOT a guess. Sherpa's own particle database
  // mass/width for this state will differ somewhat from the CLEO 3pi
  // fit's own value (M=0.860, Gamma=0.880 GeV) - confirmed acceptable
  // for phase-space importance-sampling purposes (the actual
  // amplitude calculation is unaffected either way; this only shapes
  // how efficiently events are generated near this resonance).
  if (name=="f(0)(600)")     return kf_f_0_600;
  // K_0*(700)/kappa - NEWLY registered this session (kf 9000311/
  // 9000321, see K_0_800_Decays.H/.C), needed for the Kpipi KS_f0/
  // CLEO_f0 phase-space channels' Kpi-subsystem scalar term (not yet
  // used in any live channel below, but added proactively).
  if (name=="K(0)(800)")     return kf_K_0_800;
  if (name=="K(0)(800)+")    return kf_K_0_800_plus;
  // Added proactively (not yet used in any live phase-space channel,
  // but properly registered Total_Width_Base lineshapes exist for
  // both - see A0_Decays.H/.C) in case a future eta-eta-pi/K-eta-eta
  // current (currently unmodeled) or similar ever uses either as an
  // intermediate resonance.
  if (name=="a(0)(980)+")    return kf_a_0_980_plus;
  if (name=="a(0)(1450)+")   return kf_a_0_1450_plus;
  // NOTE: the CLEO 3pi current's dormant "a1'" term (beta=0 by
  // default, currently unused - see m_betaCLEOa1p in FF_0_PPP.C) has
  // NO entry here and none was added: no kf-code for a higher-mass
  // a1 state (e.g. a1(1640)) exists anywhere in Flavour_Tags.H, so
  // there is nothing concrete to map it to. If that term is ever
  // activated, this table needs a genuinely new kf-code registered
  // first (out of scope for this file alone) before it can be
  // properly named in a phase-space channel string.
  if (name=="J/psi(1S)")     return kf_J_psi_1S;
  if (name=="psi(2S)")       return kf_psi_2S;
  if (name=="psi(4040)")     return kf_psi_4040;
  msg_Error()<<"Error in ResonanceNameToKfCode: unrecognised resonance "
	     <<"name \""<<name<<"\" in a phase-space channel descriptor - "
	     <<"falling back to kf-code "<<fallback<<"'s mass/width "
	     <<"(silently using the WRONG resonance otherwise). Add this "
	     <<"name to the lookup table above if it is a genuinely new "
	     <<"resonance, or check for a typo in the decay YAML.\n";
  return fallback;
}


bool HD_Channel_Selector::DecomposeChannel(string name, ChannelInformation & ci)
{
  ci.name = "noname";
  ci.a=0; ci.b=0; ci.c=0; ci.d=0;
  ci.res1 = "no res";
  ci.res2 = "no res";
  ci.res3 = "no res";
  
  Data_Reader reader("_",";","#");
  vector<string> exploded;
  reader.SetString(name);
  reader.VectorFromString(exploded);
  
  if(exploded.size() < 1) return false;
  
  if(exploded[0]=="Isotropic" || exploded[0]=="Iso2") {
    ci.name = exploded[0];
    ci.nRes = 0;
  }
  else if(exploded[0]=="Dalitz" && exploded.size()==3) {
    ci.name = exploded[0];
    ci.res1 = exploded[1];
    int ab = ToType<int>(exploded[2]);
    ci.b=ab%10; ci.a=ab/10;   // int/int !
    ci.nRes = 1;
  }
  else if(exploded[0]=="TwoResonances" && exploded.size()==5) {
    ci.name = exploded[0];
    ci.res1 = exploded[1]; 
    ci.a    = ToType<int>(exploded[2]);
    ci.res2 = exploded[3]; 
    int bc = ToType<int>(exploded[4]);
    ci.c=bc%10; ci.b=bc/10;   // int/int !
    ci.nRes = 2;
  }
  // ThreeResonances_res1_k_res2_l_res3_ij - direct extension of
  // TwoResonances by one more sequential resonance (see
  // Five_Body_PSs.H's ThreeResonances class). 7 tokens: name, res1,
  // k, res2, l, res3, ij (ij combined two-digit exactly like
  // TwoResonances' own trailing index pair).
  else if(exploded[0]=="ThreeResonances" && exploded.size()==7) {
    ci.name = exploded[0];
    ci.res1 = exploded[1];
    ci.a    = ToType<int>(exploded[2]);
    ci.res2 = exploded[3];
    ci.b    = ToType<int>(exploded[4]);
    ci.res3 = exploded[5];
    int cd  = ToType<int>(exploded[6]);
    ci.d=cd%10; ci.c=cd/10;   // int/int !
    ci.nRes = 3;
  }
  else if(exploded[0]=="IsotropicSpectator" && exploded.size()==2) {
    ci.name = exploded[0];
    ci.a = ToType<int>(exploded[1]); // spectator index
  }
  
  if( ci.name==string("noname") ) return false;
  return true;
}

Single_Channel * HD_Channel_Selector::GetChannel( 
    int nin, 
    int nout, 
    const Flavour * flavs, 
    string name,
    ATOOLS::Scoped_Settings& s,
    const ATOOLS::Mass_Selector* ms)
{
  if ( nin>1 || nout<1 ) {
    msg_Error()<<METHOD<<": Error: "<<endl
           <<"   No PS for channel ("<<nin<<" -> "<<nout<<" )"<<endl
           <<"   Return nothing and hope for the best."<<endl;
    return NULL;
  }
  ChannelInformation ci;
  if( DecomposeChannel( name, ci ) ) {
    if (ci.name==string("Isotropic")) {
      if ( nout == 2 ) return new Iso2Channel(flavs);
      if ( nout == 1 ) return new Iso1Channel(flavs);
      return new Rambo(1,nout,flavs,ms);
    }
  }
  if (ci.name==string("Iso2") || nout==2 ) return new Iso2Channel(flavs);
  if (nout==3) {
    if (ci.name==string("Dalitz")) {
      // Unchanged in spirit, now routed through the shared helper -
      // "W" keeps its special-cased width (2.06 GeV, not from the
      // particle table) exactly as before.
      kf_code kfres = ResonanceNameToKfCode(ci.res1, kf_rho_770_plus);
      double width = Flavour(kfres).Width();
      if( ci.res1==string("W") ) {
        kfres = kf_Wplus;
        width = 2.06;
      }
      SimpleResonanceFlavour res(
          Flavour(kfres).IDName(),
          s["Mass_"+Flavour(kfres).IDName()].SetDefault(Flavour(kfres).HadMass()).Get<double>(),
          s["Width_"+Flavour(kfres).IDName()].SetDefault(width).Get<double>());
      return new Dalitz(flavs,res,ci.a,ci.b);
    }
    if( ci.name==string("IsotropicSpectator") ) {
      return new IsotropicSpectator( flavs, nout, ci.a, ms );
    }
  }
  if (nout==4) {
    if( ci.name==string("TwoResonances") ) {
      // FIX: previously built SimpleResonanceFlavour directly from the
      // raw parsed strings ci.res1/ci.res2, with the DEFAULT mass/
      // width hardcoded to a1(1260) (res_a) / rho(770) (res_v)
      // regardless of what resonance was actually named - meaning
      // every resonance other than those exact two silently got the
      // WRONG default unless a matching "Mass_<name>"/"Width_<name>"
      // override was supplied by hand. Now resolved through the same
      // shared lookup Dalitz already uses, so any recognised name
      // gets the CORRECT mass/width from Sherpa's own particle
      // database automatically, and an unrecognised one is flagged
      // with an explicit error instead of failing silently.
      kf_code kfres_a = ResonanceNameToKfCode(ci.res1, kf_a_1_1260_plus);
      kf_code kfres_v = ResonanceNameToKfCode(ci.res2, kf_rho_770_plus);
      SimpleResonanceFlavour res_a(
          Flavour(kfres_a).IDName(),
          s["Mass_"+Flavour(kfres_a).IDName()].SetDefault(Flavour(kfres_a).HadMass()).Get<double>(),
          s["Width_"+Flavour(kfres_a).IDName()].SetDefault(Flavour(kfres_a).Width()).Get<double>());
      SimpleResonanceFlavour res_v(
          Flavour(kfres_v).IDName(),
          s["Mass_"+Flavour(kfres_v).IDName()].SetDefault(Flavour(kfres_v).HadMass()).Get<double>(),
          s["Width_"+Flavour(kfres_v).IDName()].SetDefault(Flavour(kfres_v).Width()).Get<double>());
      return new TwoResonances( flavs, res_a, ci.a, res_v, ci.b, ci.c );
    }
    if( ci.name==string("IsotropicSpectator") ) {
      return new IsotropicSpectator( flavs, nout, ci.a, ms );
    }
  }
  if (nout==5) {
    if( ci.name==string("ThreeResonances") ) {
      // Direct extension of the TwoResonances (nout==4) branch by one
      // more resonance - same lookup mechanism, same default-fallback
      // reasoning (res1/res2 default to a1(1260)/rho(770) as before;
      // res3, being the innermost resonance and thus most often
      // vector-like in practice, defaults to rho(770) too - override
      // via Mass_<name>/Width_<name> as needed, same as the others).
      kf_code kfres_a = ResonanceNameToKfCode(ci.res1, kf_a_1_1260_plus);
      kf_code kfres_b = ResonanceNameToKfCode(ci.res2, kf_rho_770_plus);
      kf_code kfres_c = ResonanceNameToKfCode(ci.res3, kf_rho_770_plus);
      SimpleResonanceFlavour res_a(
          Flavour(kfres_a).IDName(),
          s["Mass_"+Flavour(kfres_a).IDName()].SetDefault(Flavour(kfres_a).HadMass()).Get<double>(),
          s["Width_"+Flavour(kfres_a).IDName()].SetDefault(Flavour(kfres_a).Width()).Get<double>());
      SimpleResonanceFlavour res_b(
          Flavour(kfres_b).IDName(),
          s["Mass_"+Flavour(kfres_b).IDName()].SetDefault(Flavour(kfres_b).HadMass()).Get<double>(),
          s["Width_"+Flavour(kfres_b).IDName()].SetDefault(Flavour(kfres_b).Width()).Get<double>());
      SimpleResonanceFlavour res_c(
          Flavour(kfres_c).IDName(),
          s["Mass_"+Flavour(kfres_c).IDName()].SetDefault(Flavour(kfres_c).HadMass()).Get<double>(),
          s["Width_"+Flavour(kfres_c).IDName()].SetDefault(Flavour(kfres_c).Width()).Get<double>());
      return new ThreeResonances( flavs, res_a, ci.a, res_b, ci.b, res_c, ci.c, ci.d );
    }
    if( ci.name==string("IsotropicSpectator") ) {
      return new IsotropicSpectator( flavs, nout, ci.a, ms );
    }
  }

  msg_Error()<<METHOD<<": Error: "<<endl
    <<"   No channel for ("<<nin<<" -> "<<nout<<") with name "<<name<<endl
         <<"   Return nothing and hope for the best."<<endl;
  return NULL;
}
