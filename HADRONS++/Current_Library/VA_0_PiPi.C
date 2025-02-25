#include "HADRONS++/Current_Library/VA_0_PiPi.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"

using namespace HADRONS;
using namespace ATOOLS;
using namespace std;


///////////////////////////////////////////////////////////////////////////
//
// Comment to be deleted once we have it all checked ....
// Please take a look at a summary of what is implemented in Tauola
// - https://arxiv.org/pdf/1509.09140
// - https://arxiv.org/abs/1609.04617
// and "harvest" the references within.
//
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
//
// Form factors for charged pi pi, K K, K pi final-state currents from:
// - RChT, 2: Resonance Chiral Perturbation Theory
//   * pi pi/KK fit: Eur.Phys.J.C 79 (2019) 5, 436
//     (https://doi.org/10.1140/epjc/s10052-019-6943-9)
//           - We could also have some alternative dispersive approach, 
//             maybe as form factor model 3, based on RchT
//             (https://arxiv.org/pdf/1112.0962)
//             === Zara - maybe this would be a nice way to promote
//                 yourself from "debugging" to "coding"?
//   * K pi: Phys.Lett.B 640 (2006) 176
//     (https://doi.org/10.1016/j.physletb.2006.06.058)
//           - not implemented (but maybe worth it, especially because
//             of the treatment of eta and eta'?):
//             (https://arxiv.org/pdf/1201.1794)
// - KS, 1 (Kuehn-Santamaria model):
//   * pi pi (original version): Z.Phys.C 48 (1990) 445-452
//     (https://doi.org/10.1007/BF01572024)
//   * pi pi (KS with Gounaris-Sakurai form factors):
//     from Belle measurement & analysis:
//     (https://arxiv.org/pdf/0805.3773)
//   * pi eta (KS inspired): Phys.Rev.D 82 (2010) 057301
//     (https://doi.org/10.1103/PhysRevD.82.057301)
//   * pi eta' (KS inspired): Phys.Rev.D 84 (2011) 017302
//     (https://doi.org/10.1103/PhysRevD.84.017302)
//   * K pi: Z.Phys.C 69 (1996) 243 (vector form factor)
//     (https://doi.org/10.1007/s002880050024)
//   * K pi: Z.Phys.C 72 (1996) 619 (scalar form factor)
//     (https://doi.org/10.1007/s002880050284)
// - none, 0: no form factor
//
///////////////////////////////////////////////////////////////////////////

VA_0_PiPi::VA_0_PiPi(const ATOOLS::Flavour_Vector& flavs,
		     const std::vector<int>& indices,
		     const std::string& name) :
  Current_Base(flavs, indices, name),
  m_restype(resonance_type::running),
  m_PSmode(PSmode::unknown),
  m_ffmodel(ffmodel::KS),
  m_norm(1.), m_deltaM2(0.),
  m_m2_pi(sqr(Flavour(kf_pi_plus).Mass(true))),
  m_m2_K(sqr(Flavour(kf_K_plus).Mass(true))),
  m_m2_eta(sqr(Flavour(kf_eta).Mass(true))),
  m_mu2(sqr(Flavour(kf_rho_770_plus).Mass(true)))
{
  msg_Out()<<"==========================================================\n";
}

VA_0_PiPi::~VA_0_PiPi() {
  while (!m_vectors.empty()) {
    delete m_vectors.front().second;
    m_vectors.pop_front();
  }
  while (!m_scalars.empty()) {
    delete m_scalars.front().second;
    m_scalars.pop_front();
  }
}


void VA_0_PiPi::Calc(const ATOOLS::Vec4D_Vector& moms, bool m_anti)
{
  Vec4D  q  = moms[p_i[1]]+moms[p_i[0]];
  double Q2 = q.Abs2();
  msg_Out()<<METHOD<<"(Q = "<<sqrt(Q2)<<"): "
	   <<"V = "<<VectorFF(Q2)<<", S = "<<ScalarFF(Q2)<<"\n";
  Insert(m_norm *
	 ( VectorFF(Q2) * ((moms[p_i[1]]-moms[p_i[0]]) - m_deltaM2/Q2 * q ) +
	   ScalarFF(Q2) * q
	   ),   0);
}

void VA_0_PiPi::SetModelParameters(struct GeneralModel model)
{
  FixMode();
  FixNorm(model);
  //m_fa0(0.015),
  //m_GS(0.), m_eps_etapi(0.0134), m_eps_etaprimepi(0.003),

  m_mu2     = sqr(model("mu", Flavour(kf_rho_770_plus).Mass(true)));
  m_fpi     = model("fpi", 0.13041 );
  m_GV      = model("GV",  0.0619 );
  m_ffmodel = ffmodel(model("FORM_FACTOR",1));
  m_deltaM2 = m_flavs[p_i[1]].Mass(true)-m_flavs[p_i[0]].Mass(true);
  if (m_ffmodel==ffmodel::RChT) {
    m_fpi   = 0.092316;
    if (m_PSmode==PSmode::pipi || m_PSmode==PSmode::KK) {
      m_RChTpref = -1./(96.*sqr(M_PI*m_fpi));
    }
    else if (m_PSmode==PSmode::Kpi || m_PSmode==PSmode::Keta) {
      m_RChTpref = (sqr(m_GV)*Flavour(kf_K_star_892_plus).Mass(true)/
		    (64.*M_PI*pow(m_fpi,4)));
    }
  }
  else if (m_ffmodel==ffmodel::KS) {
    if (m_PSmode==PSmode::etapi) {
      m_GS = m_fa0*m_g_a0_eta_pi/m_deltaM2;
    }
    if (m_PSmode==PSmode::etapi) {
      m_GS = m_fa0*m_g_a0_etaprime_pi/m_deltaM2;
    }
  }
  switch (int(model("RUNNING_WIDTH",1))) {
  case 10: m_restype = resonance_type::bespoke; break; 
  case 2:  m_restype = resonance_type::GS;      break; 
  case 1:  m_restype = resonance_type::running; break;
  case 0:  m_restype = resonance_type::fixed;   break;
  default: m_restype = resonance_type::running; break;
  }
  if (m_ffmodel!=ffmodel::none) { 
    list<kf_code> tags;
    if (SelectResonances(tags,false)) {
      for (list<kf_code>::iterator kfit=tags.begin();kfit!=tags.end();kfit++)
	InitResonance(model, *kfit, false);
    }
    tags.clear();
    if ( (m_PSmode==PSmode::Kpi ||
	  m_PSmode==PSmode::etapi || m_PSmode==PSmode::etaprimepi ||
	  m_PSmode==PSmode::Keta  || m_PSmode==PSmode::Ketaprime) &&
	  SelectResonances(tags,true)) {
      for (list<kf_code>::iterator kfit=tags.begin();kfit!=tags.end();kfit++)
	InitResonance(model, *kfit, true);
    }
    if (m_vectors.empty() && m_scalars.empty()) {
      THROW(fatal_error,"Could not find suitable vector and scalar resonances.");
    }
  }
}

void VA_0_PiPi::FixMode() {
  if (m_flavs[p_i[0]].Kfcode()==kf_pi &&
      m_flavs[p_i[1]].Kfcode()==kf_pi_plus)             m_PSmode = PSmode::pipi;
  else if ( (m_flavs[p_i[0]].Kfcode()==kf_K_L ||
	     m_flavs[p_i[0]].Kfcode()==kf_K_S ||
	     m_flavs[p_i[0]].Kfcode()==kf_K) &&
	    m_flavs[p_i[1]].Kfcode()==kf_K_plus)        m_PSmode = PSmode::KK;
  else if ( (m_flavs[p_i[0]].Kfcode()==kf_pi_plus && 
	     (m_flavs[p_i[1]].Kfcode()==kf_K_S ||
	      m_flavs[p_i[1]].Kfcode()==kf_K_L ||
	      m_flavs[p_i[1]].Kfcode()==kf_K) ) ||
	    (m_flavs[p_i[0]].Kfcode()==kf_pi && 
	     (m_flavs[p_i[1]].Kfcode()==kf_K_plus) ) )  m_PSmode = PSmode::Kpi;
  else if ( m_flavs[p_i[0]].Kfcode()==kf_pi_plus &&
	    m_flavs[p_i[1]].Kfcode()==kf_eta )          m_PSmode = PSmode::etapi;
  else if ( m_flavs[p_i[0]].Kfcode()==kf_pi_plus &&
	    m_flavs[p_i[1]].Kfcode()==kf_eta_prime_958) m_PSmode = PSmode::etaprimepi;
  else if ( m_flavs[p_i[0]].Kfcode()==kf_eta &&
	    m_flavs[p_i[1]].Kfcode()==kf_K_plus)        m_PSmode = PSmode::Keta;
  else if ( m_flavs[p_i[0]].Kfcode()==kf_eta_prime_958 &&
	    m_flavs[p_i[1]].Kfcode()==kf_pi_plus)       m_PSmode = PSmode::Ketaprime;
  if (m_PSmode==PSmode::unknown)
    THROW(fatal_error,"Current called for illegal flavour combination.");
}

void VA_0_PiPi::FixNorm(struct GeneralModel & model) {
  // global pre-factor: 1/sqrt(2) for pi_0 wave-function, V_ud for the
  // quark-level coupling producing a rho (or rho-resonance), 1/sqrt(2)
  // for the overall normalisation.
  double iso = 0.;
  switch (int(m_PSmode)) {
  case int(PSmode::KK):
    iso = 1./sqrt(2.);
    break;
  case int(PSmode::Kpi):
    if (m_flavs[p_i[0]].Kfcode()==kf_pi_plus) iso = 1./2.;
    else if (m_flavs[p_i[0]].Kfcode()==kf_pi) iso = 1./sqrt(2.);
    break;
  case int(PSmode::pipi):
    iso = 1.;
    break;
  default: break;
  }
  double CKM = 0.;
  if (m_PSmode==PSmode::pipi  || m_PSmode==PSmode::KK ||
      m_PSmode==PSmode::etapi || m_PSmode==PSmode::etaprimepi)
    CKM = model("Vud", Tools::Vud);
  else if (m_PSmode==PSmode::Kpi || m_PSmode==PSmode::Keta ||
	   m_PSmode==PSmode::Ketaprime) 
    CKM = model("Vus", Tools::Vus);
  m_norm = iso * CKM / sqrt(0.5);
  if (m_norm<=0.) THROW(fatal_error,"Current with zero norm.");
}
  
bool VA_0_PiPi::SelectResonances(list<kf_code> & tags,const bool & isScalar) {
  if (isScalar) {
    if (m_PSmode==PSmode::etapi   ||
	m_PSmode==PSmode::etaprimepi) {
      tags.push_back(kf_a_0_980_plus);
      tags.push_back(kf_a_0_1450_plus);
      return true;
    }
    else if (m_PSmode==PSmode::Kpi   ||
	     m_PSmode==PSmode::Keta  ||
	     m_PSmode==PSmode::Ketaprime) {
      tags.push_back(kf_K_0_star_1430_plus);
      return true;
    }
  }
  else {
    if (m_PSmode==PSmode::pipi    ||
	m_PSmode==PSmode::KK      ||
	m_PSmode==PSmode::etapi   ||
	m_PSmode==PSmode::etaprimepi) {
      tags.push_back(kf_rho_770_plus);
      tags.push_back(kf_rho_1450_plus);
      tags.push_back(kf_rho_1700_plus);
      return true;
    }
    else if (m_PSmode==PSmode::Kpi   ||
	     m_PSmode==PSmode::Keta  ||
	     m_PSmode==PSmode::Ketaprime) {
      tags.push_back(kf_K_star_892_plus);
      tags.push_back(kf_K_star_1410_plus);
      return true;
    }
  }
  return false;
}

bool VA_0_PiPi::FillDefaults(const kf_code & tag,vector<double> & defs) {
  // Default values as vector of 4 doubles:
  // weight, mass, width, phase
  switch (m_ffmodel) {
  case ffmodel::KS:
    return KSDefaults(tag,defs);
  case ffmodel::RChT:
    return RChTDefaults(tag,defs);
  }
  THROW(fatal_error,"unknown form factor model.");
}

bool VA_0_PiPi::KSDefaults(const kf_code & tag,vector<double> & defs) {
  if (m_PSmode==PSmode::pipi) {
    switch (tag) {
    case kf_rho_1700_plus:
      if (m_restype==resonance_type::running)
	defs = { -0.0370, 1.7000, 0.2350, 0.0000 };
      else if (m_restype==resonance_type::GS)
	defs = {  0.0280, 1.6940, 0.1350, -3./180.*M_PI };
      break;
    case kf_rho_1450_plus: 
      if (m_restype==resonance_type::running)
	defs = { -0.1030, 1.3200, 0.3900, 0.0000 };
      else if (m_restype==resonance_type::GS)
	defs = {  0.1300, 1.4280, 0.4130, 197./180.*M_PI };
      break;
    case kf_rho_770_plus:  
      if (m_restype==resonance_type::running)
	defs = {  1.0000, 0.7749, 0.1480, 0.0000 };
      else if (m_restype==resonance_type::GS)
	defs = {  1.0000, 0.7740, 0.1440, 0.0000 };
      break;
    }
  }
  else if (m_PSmode==PSmode::etapi || m_PSmode==PSmode::etaprimepi) {
    switch (tag) {
    case kf_rho_1700_plus:
      defs = { -0.0370, 1.7000, 0.2350, 0.0000 }; break;
    case kf_rho_1450_plus: 
      defs = { -0.1030, 1.3200, 0.3900, 0.0000 }; break;
    case kf_rho_770_plus:  
      defs = {  1.0000, 0.7749, 0.1480, 0.0000 }; break;
    case kf_a_0_1450_plus:  
      defs = {  0.0000, 0.1450, 0.2650, 0.0000 }; break;
    case kf_a_0_980_plus:  
      defs = {  1.0000, 0.9800, 0.1000, 0.0000 }; break;
    }
  }
  else if (m_PSmode==PSmode::KK && m_restype==resonance_type::running) {
    switch (tag) {
    case kf_rho_1700_plus:
      defs = {  0.0500, 1.7000, 0.2350, 0.0000 }; break;
    case kf_rho_1450_plus:
      defs = { -0.1670, 1.3630, 0.3100, 0.0000 }; break;
    case kf_rho_770_plus:  
      defs = {  1.0000, 0.7749, 0.1480, 0.0000 }; break;
    }
  }
  else if ((m_PSmode==PSmode::Kpi  ||
	    m_PSmode==PSmode::Keta || m_PSmode==PSmode::Ketaprime) &&
	   m_restype==resonance_type::running) {
    switch (tag) {
    case kf_K_star_1410_plus:
      defs = { -0.1350, 1.4120, 0.2270, 0.0000 }; break;
    case kf_K_star_892_plus:      
      defs = {  1.0000, 0.8920, 0.0500, 0.0000 }; break;
    case kf_K_0_star_1430_plus:      
      defs = {  1.7000, 1.4300, 0.2870, 0.0000 }; break;
    }
  }
  if (!defs.empty()) return true;
  msg_Error()<<METHOD<<" will throw error:\n"
	     <<"Combination of decay mode = "<<int(m_PSmode)<<", "
	     <<"resonancve type = "<<int(m_restype)<<" and "
	     <<"flavour = "<<tag<<" not yet implemented.\n";
  THROW(fatal_error,"unknown combination in KS defaults");
  return false;
}

bool VA_0_PiPi::RChTDefaults(const kf_code & tag,vector<double> & defs) {
  if (m_PSmode==PSmode::pipi || m_PSmode==PSmode::KK) {
    switch (tag) {
    case kf_rho_1700_plus:
      defs ={ -0.1200, 1.7540, 0.4120, -0.0200 }; break;
    case kf_rho_1450_plus: 
      defs = { 0.1800, 1.4380, 0.5350, -0.3600 }; break;
    case kf_rho_770_plus:  
      defs = { 1.0000, 0.7752, 0.0000,  0.0000 }; break;
    }
  }
  else if (m_PSmode==PSmode::Kpi) {
    switch (tag) {
    case kf_K_star_1410_plus:
      defs = {  0.0130, 1.4140, 0.2320, 0.0000 }; break;
    case kf_K_star_892_plus:      
      defs = {  1.0000, 0.8917, 0.0500, 0.0000 }; break;
    case kf_K_0_star_1430_plus:      
      defs = {  1.7000, 1.4300, 0.2870, 0.0000 }; break;
    }
  }
  if (!defs.empty()) return true;
  msg_Error()<<METHOD<<" will throw error:\n"
	     <<"Combination of decay mode = "<<int(m_PSmode)<<", "
	     <<"resonancve type = "<<int(m_restype)<<" and "
	     <<"flavour = "<<tag<<" not yet implemented.\n";
  THROW(fatal_error,"unknown combination in KS defaults");
  return false;
}
    
void VA_0_PiPi::InitResonance(struct GeneralModel & model,const kf_code & tag,
			      const bool & isScalar)
{
  msg_Out()<<METHOD<<"(mode = "<<int(m_PSmode)<<" for "
	   <<int(tag)<<", scalar = "<<isScalar<<")\n";
  vector<double> defaults;
  if (!FillDefaults(tag,defaults) || dabs(defaults[0])<1.e-8) return;
  string id = string("");
  switch (tag) {
  case kf_rho_1700_plus:      id = string("rho(1700)+");  break;
  case kf_rho_1450_plus:      id = string("rho(1450)+");  break;
  case kf_rho_770_plus:       id = string("rho(770)+");   break;
  case kf_K_star_1410_plus:   id = string("K*(1410)+");   break;
  case kf_K_star_892_plus:    id = string("K*(892)+");    break;
  case kf_K_0_star_1430_plus: id = string("K*_0(1430)+"); break;
  case kf_a_0_980_plus:       id = string("a_0(980)+");   break;
  case kf_a_0_1450_plus:      id = string("a_0(1450)+");  break;
  default: break;
  }
  vector<Flavour> outflavs;
  if (m_PSmode==PSmode::pipi || m_PSmode==PSmode::KK) {
    outflavs.push_back(Flavour(kf_pi_plus));
    outflavs.push_back(Flavour(kf_pi));
  }
  else if (m_PSmode==PSmode::Kpi ||
	   m_PSmode==PSmode::etapi || m_PSmode==PSmode::etaprimepi ||
	   m_PSmode==PSmode::Keta  || m_PSmode==PSmode::Ketaprime) {
    outflavs.push_back(m_flavs[p_i[0]]);
    outflavs.push_back(m_flavs[p_i[1]]);
  }
  Flavour flav = Flavour(tag);
  double weight = model(string("Weight_")+id, defaults[0]);
  double mass   = model(string("Mass_")  +id, defaults[1]);
  double width  = model(string("Width_") +id, defaults[2]);
  double phase  = model(string("Phase_") +id, defaults[3]);

  Resonance_Parameters params(flav,outflavs,m_restype,mass,width,phase);
  Resonance_Base * res = NULL;
  switch (m_restype) {
  case resonance_type::fixed:
    res = new FixedWidth_Resonance(params);
    break;
  case resonance_type::GS:
    res = new GS_Resonance(params);
    break;
  case resonance_type::running:
  default:
    res = new RunningWidth_Resonance(params);
    break;
  }
  if (!isScalar) m_vectors.push_back(make_pair(weight,res));
  else           m_scalars.push_back(make_pair(weight,res));
}

Complex VA_0_PiPi::VectorFF(const double & s) {
  if (m_vectors.empty() && m_ffmodel!=ffmodel::none) return Complex(0.,0.); 
  if (m_PSmode==PSmode::pipi && s<4.*m_m2_pi) return Complex(0.,0.);
  if (m_PSmode==PSmode::KK   && s<4.*m_m2_K)  return Complex(0.,0.);
  switch (int(m_ffmodel)) {
  case int(ffmodel::RChT): return VectorRChT(s);
  case int(ffmodel::KS):   return VectorKS(s);
  case int(ffmodel::none): 
  default: break;
  }
  return Complex(1.,0.);
}

Complex VA_0_PiPi::ScalarFF(const double & s) {
  if (m_scalars.empty() && m_ffmodel!=ffmodel::none)  return Complex(0.,0.); 
  if (m_PSmode==PSmode::pipi || m_PSmode==PSmode::KK) return Complex(0.,0.);
  switch (int(m_ffmodel)) {
  case int(ffmodel::RChT): return ScalarRChT(s);
  case int(ffmodel::KS):   return ScalarKS(s);
  case int(ffmodel::none): 
  default: break;
  }
  return Complex(1.,0.);
}

Complex VA_0_PiPi::VectorKS(const double & s)
{
  Complex ff = Complex(0.,0.), norm = Complex(0.,0.), pref;
  for (list<pair<double, Resonance_Base * > >::iterator rtit=m_vectors.begin();
       rtit!=m_vectors.end();rtit++) {
    pref  = rtit->first * Complex(cos(rtit->second->Phase()),
				  sin(rtit->second->Phase()));
    ff   += pref * rtit->second->BreitWigner(s);
    norm += pref;
  }
  return ff/norm;
}

Complex VA_0_PiPi::ScalarKS(const double & s)
{
  Complex ff = Complex(0.,0.), pref;
  for (list<pair<double, Resonance_Base * > >::iterator rtit=m_vectors.begin();
       rtit!=m_vectors.end();rtit++) {
    pref  = ( rtit->first * m_deltaM2/rtit->second->Mass2() *
	      Complex(cos(rtit->second->Phase()),sin(rtit->second->Phase())) );
    ff   += pref * rtit->second->BreitWigner(s);
  }
  return ff;
}
  
Complex VA_0_PiPi::VectorRChT(const double & s) {
  switch (m_PSmode) {
  case PSmode::pipi:
  case PSmode::KK:   return VectorRChT_pipi(s);
  case PSmode::Kpi:  return VectorRChT_Kpi(s);
  default: break;
  }
  return Complex(0.,0.);
}

Complex VA_0_PiPi::VectorRChT_pipi(const double & s)
{
  Complex ff      = Complex(0.,0.);
  Complex ExtraBW = Complex(0.,0.);
  double  RhoExpo = 0., Gamma = 0., pref;
  for (list<pair<double, Resonance_Base * > >::iterator
	 rtit=m_vectors.begin();rtit!=m_vectors.end();rtit++) {
    if (rtit->second->Flav()==Flavour(kf_rho_770_plus)) {
      Complex Ampl = ( Loop(m_m2_pi, s, m_mu2) +
		       Loop(m_m2_K,  s, m_mu2)/2. );
      // remember that  m_RChTpref = -1./(96.*sqr(M_PI*m_fpi))
      RhoExpo = exp(m_RChTpref*s*Ampl.real());
      // rho width repressed through the imaginar part of loop integral:
      // automatically captures the mass-shell conditions
      Gamma   = m_RChTpref*rtit->second->Mass()*s*Ampl.imag();
      ExtraBW = rtit->second->AltBreitWigner(s,rtit->second->Mass2(),
					     -rtit->second->Mass()*Gamma);
      ff  += (rtit->first *
	      rtit->second->BreitWigner(s,rtit->second->Mass2(),
					-rtit->second->Mass()*Gamma) *
	      RhoExpo );
    }
    else {
      Gamma = (s>4.*m_m2_pi ? rtit->second->Width(s) : 0.);
      pref = - ( s*rtit->second->Width()/M_PI/
		 pow(2.*rtit->second->Lambda(rtit->second->Mass2(),m_m2_pi,m_m2_pi), 3) );
      ff += - (rtit->first * exp(Complex(0.,rtit->second->Phase())) *
	       ( rtit->second->AltBreitWigner(s,rtit->second->Mass2(),
					      -rtit->second->Mass()*Gamma) *
		 exp(pref * Loop(m_m2_pi, s, m_mu2).real()) -
		 ExtraBW * RhoExpo) );
    }
  }
  return ff;
}

Complex VA_0_PiPi::VectorRChT_Kpi(const double & s)
{
  Complex ff      = Complex(0.,0.);
  Complex ExtraBW = Complex(0.,0.);
  double  Expo    = exp(3/2.*(Htilde(s,m_m2_K,m_m2_pi)+
			      Htilde(s,m_m2_K,m_m2_eta)).real());
  double  Gamma;
  for (list<pair<double, Resonance_Base * > >::iterator
	 rtit=m_vectors.begin();rtit!=m_vectors.end();rtit++) {
    if (rtit->second->Flav()==Flavour(kf_K_star_892_plus)) {
      Gamma   = (m_RChTpref * 8./sqrt(s) *
		 (pow(rtit->second->Lambda(s,m_m2_K,m_m2_pi),3)+
		  pow(rtit->second->Lambda(s,m_m2_K,m_m2_eta),3)));
      ExtraBW = rtit->second->AltBreitWigner(s,rtit->second->Mass2(),
					     -rtit->second->Mass()*Gamma);
      ff     += (rtit->first *
		 rtit->second->BreitWigner(s,rtit->second->Mass2(),
					   -rtit->second->Mass()*Gamma) );
    }
    else {
      Gamma = (s>sqr(sqrt(m_m2_pi)+sqrt(m_m2_K)) ? rtit->second->Width(s) : 0.);
      ff     += (rtit->first *
		 rtit->second->AltBreitWigner(s,rtit->second->Mass2(),
					      -rtit->second->Mass()*Gamma) );
      ff     -= rtit->first * ExtraBW;
    }
	
  }
  return ff * Expo;
}


Complex VA_0_PiPi::ScalarRChT(const double & s) {
  switch (m_PSmode) {
  case PSmode::pipi:
  case PSmode::KK:   return Complex(0.,0.);
  case PSmode::Kpi:  // return ScalarRChT_Kpi(s);
  default: break;
  }
  return Complex(0.,0.);
}



Complex VA_0_PiPi::Loop(const double & m2, const double & s,const double & mu2) {
  Complex sigma = csqrt(1.-4.*m2/s);
  return ( log(m2/mu2) + 8.*m2/s - 5./3. +
	   pow(sigma,3.)*log((sigma+1.)/(sigma-1.)) );
}
 
Complex VA_0_PiPi::Htilde(const double & s, const double & MP2,const double & MQ2) {
  return (s*Mr(s,MP2,MQ2)-L(s,MP2,MQ2))/sqr(m_fpi);
}

Complex VA_0_PiPi::L(const double & s, const double & MP2,const double & MQ2) {
  return sqr(MP2-MQ2)/(4.*s) * JBar(s,MP2,MQ2);
}

Complex VA_0_PiPi::Mr(const double & s, const double & MP2,const double & MQ2) {
  double sigma = MP2+MQ2, delta = MP2-MQ2, delta2 = sqr(delta);
  return ( (s-2.*sigma)/(12.*s) * JBar(s,MP2,MQ2) +
	   delta2/(3.*sqr(s)) * JBarBar(s,MP2,MQ2) -
	   1./6. * k(MP2,MQ2) + 1./(288.*sqr(M_PI)) );
}

Complex VA_0_PiPi::k(const double & MP2,const double & MQ2) {
  return (MP2*log(MP2/m_mu2)-MQ2*log(MQ2/m_mu2))/(32.*sqr(M_PI)*(MP2-MQ2));
}
    
Complex VA_0_PiPi::JBar(const double & s, const double & MP2,const double & MQ2) {
  double disc  = sqr(s-MP2-MQ2)-4.*MP2*MQ2;
  if (disc<0.) return 0.;
  double nu    = sqrt(disc), nuPs = sqr(nu+s), nuMs = sqr(nu-s);
  double sigma = MP2+MQ2, delta = MP2-MQ2, delta2 = sqr(delta);
  return ( (2. + (delta/s-sigma/delta)*log(MQ2/MP2) -
	    nu/s*log((nuPs-delta2)/(nuMs-delta2)))/
	   (32.*sqr(M_PI)) );
}

Complex VA_0_PiPi::JBarBar(const double & s, const double & MP2,const double & MQ2) {
  double sigma    = MP2+MQ2, delta = MP2-MQ2, delta2 = sqr(delta);
  double Jprime_0 = (sigma + 2.*MP2*MQ2/delta*log(MQ2/MP2))/(32.*sqr(M_PI * delta));
  return JBar(s,MP2,MQ2) - Jprime_0;
}


// need to update the references once we have the tau's under control
DEFINE_CURRENT_GETTER(HADRONS::VA_0_PiPi,"VA_0_PiPi")

void ATOOLS::Getter<HADRONS::Current_Base,
		    HADRONS::ME_Parameters,HADRONS::VA_0_PiPi>::
PrintInfo(std::ostream &st,const size_t width) const {
  st<<"Example: $ 0 \\rightarrow \\pi \\pi $ \n\n"
    <<"Order: 0 = $\\pi^0$, 1 = $\\pi^\\pm$ \n\n"
    <<"Available form factors: \n "
    <<"  \\begin{itemize} \n"
    <<"    \\item {\\tt FORM\\_FACTOR = 1 :} Kuehn-Santamaria \n"
    <<"    \\item {\\tt FORM\\_FACTOR = 2 :} Resonance Chiral Theory \n"
    <<"  \\end{itemize} \n"
    <<"Reference: https://sherpa.hepforge.org/olddokuwiki/data/media/publications/theses/diplom\\__laubrich.pdf \n"
    <<std::endl;
}
