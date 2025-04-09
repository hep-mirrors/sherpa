#include "HADRONS++/Current_Library/VA_0_PiPiPi.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"

using namespace HADRONS;
using namespace ATOOLS;
using namespace std;

///////////////////////////////////////////////////////////////////////////
//
// Form factors for pi pi pi, K K, K pi final-state currents from:
// - KS, 1 (Kuehn-Santamaria model):
//   * pi pi pi (original version): Z.Phys.C 48 (1990) 445-452
//     (https://doi.org/10.1007/BF01572024)
// - RChT, 2 (Resonance Chiral Theory):
//   * pi pi pi (original version): Phys.Lett.B 685 (2010) 158
//     (https://doi.org/10.1016/j.physletb.2010.01.059)
//     + addition of sigma(600) in: Phys.Rev.D 88 (2013) 093012
//     (https://doi.org/10.1103/PhysRevD.88.093012)
// - none, 0: no form factor
//
///////////////////////////////////////////////////////////////////////////

VA_0_PiPiPi::VA_0_PiPiPi(const ATOOLS::Flavour_Vector& flavs,
		       const std::vector<int>& indices,
		       const std::string& name) :
  Current_Base(flavs, indices, name),
  m_restype(resonance_type::running),
  m_PSmode(PSmode::unknown),
  m_ffmodel(ffmodel::KS),
  m_fpi(0.13041), m_m2_pi(sqr(Flavour(kf_pi_plus).Mass(true))),
  m_mu2(sqr(Flavour(kf_rho_770_plus).Mass(true))),
  m_c1(Complex(0.,0.)),m_c2(Complex(0.,0.)),
  m_c4(Complex(0.,0.)),m_c5(Complex(0.,0.))
{
  msg_Out()<<METHOD<<": "
  	   <<"{"<<m_flavs[p_i[0]].Kfcode()<<", "
  	   <<m_flavs[p_i[1]].Kfcode()<<", "<<m_flavs[p_i[2]].Kfcode()<<"}\n";
}

VA_0_PiPiPi::~VA_0_PiPiPi() {
}

void VA_0_PiPiPi::Calc(const ATOOLS::Vec4D_Vector& moms, bool m_anti)
{
  Vec4D  qijk = moms[p_i[m_i]]+moms[p_i[m_j]]+moms[p_i[m_k]];
  Vec4D  qij  = moms[p_i[m_i]]+moms[p_i[m_j]], dij = moms[p_i[m_i]]-moms[p_i[m_j]];
  Vec4D  qkj  = moms[p_i[m_k]]+moms[p_i[m_j]], dkj = moms[p_i[m_k]]-moms[p_i[m_j]];
  double sijk = qijk.Abs2(), sij = qij.Abs2(), skj = qkj.Abs2();

  /*
  msg_Out()<<METHOD<<"("<<int(m_PSmode)<<"):\n"
	   <<m_flavs[p_i[m_i]]<<": "<<moms[p_i[m_i]]
	   <<" ("<<sqrt(moms[p_i[m_i]].Abs2())<<")\n" 
	   <<m_flavs[p_i[m_j]]<<": "<<moms[p_i[m_j]]
	   <<" ("<<sqrt(moms[p_i[m_j]].Abs2())<<")\n" 
	   <<m_flavs[p_i[m_k]]<<": "<<moms[p_i[m_k]]
	   <<" ("<<sqrt(moms[p_i[m_k]].Abs2())<<")\n";

  msg_Out()<<" V1("<<sqrt(sij)<<") = "<<FF_V1(sijk,sij)<<", "
	   <<" V2("<<sqrt(skj)<<") = "<<FF_V2(sijk,skj)<<"\n";
  */
  //msg_Out()<<METHOD<<": "
  //	   <<m_c1<<" * "<<dij<<"\n"
  //	   <<m_c2<<" * "<<dkj<<"\n";

  Insert( m_norm * (   m_c1 * FF_V1(sijk,sij) * dij +
		       m_c2 * FF_V2(sijk,skj) * dkj +
		       m_c4 * FF_A(sijk)      * qijk      +
		     ( m_c5 * FF_P(sijk)    *
		       Vec4C(cross(moms[p_i[m_i]],moms[p_i[m_j]],moms[p_i[m_k]])) )
		     ), 0);
}

void VA_0_PiPiPi::SetModelParameters(struct GeneralModel model)
{
  FixMode();
  FixMappings();
  FixNorm(model);
  m_fpi     = model("fpi", 0.0933 );
  m_ffmodel = ffmodel(model("FORM_FACTOR",1));
  switch (int(model("RUNNING_WIDTH",1))) {
  case 10: m_restype = resonance_type::bespoke; break; 
  case 1:  m_restype = resonance_type::running; break;
  case 0:  m_restype = resonance_type::fixed;   break;
  default: m_restype = resonance_type::running; break;
  }
  if (m_ffmodel!=ffmodel::none) {
    SelectResonances(model);
  }
  if (m_PSmode==PSmode::unknown)
    THROW(fatal_error,"Current called for illegal flavour combination.");
}

void VA_0_PiPiPi::FixMode() {
  if (m_flavs[p_i[0]].Kfcode()==kf_pi_plus &&
      m_flavs[p_i[1]].Kfcode()==kf_pi_plus &&
      m_flavs[p_i[2]].Kfcode()==kf_pi_plus)  {
    m_c1 = m_c2 = Complex(0.,-2.*sqrt(2.)/3./m_fpi);
    m_PSmode = PSmode::pipipi_pmp;
  }
  else if (m_flavs[p_i[0]].Kfcode()==kf_pi &&
	   m_flavs[p_i[1]].Kfcode()==kf_pi &&
	   m_flavs[p_i[2]].Kfcode()==kf_pi_plus) m_PSmode = PSmode::pipipi_p00;
  else if (m_flavs[p_i[0]].Kfcode()==kf_pi_plus &&
	   m_flavs[p_i[1]].Kfcode()==kf_pi_plus &&
	   m_flavs[p_i[2]].Kfcode()==kf_K_plus)  m_PSmode = PSmode::Kpipi_pmp;
  else if (m_flavs[p_i[0]].Kfcode()==kf_pi &&
	   m_flavs[p_i[1]].Kfcode()==kf_pi &&
	   m_flavs[p_i[2]].Kfcode()==kf_K_plus)  m_PSmode = PSmode::Kpipi_p00;
  else if (m_flavs[p_i[0]].Kfcode()==kf_pi_plus &&
	   (m_flavs[p_i[1]].Kfcode()==kf_K_S ||
	    m_flavs[p_i[1]].Kfcode()==kf_K_L ||
	    m_flavs[p_i[1]].Kfcode()==kf_K ) &&
	   m_flavs[p_i[2]].Kfcode()==kf_pi)      m_PSmode = PSmode::piKpi_p00;
  else if (m_flavs[p_i[0]].Kfcode()==kf_pi_plus &&
	   m_flavs[p_i[1]].Kfcode()==kf_K_plus  &&
	   m_flavs[p_i[2]].Kfcode()==kf_pi_plus) m_PSmode = PSmode::piKpi_ppm;
}

void VA_0_PiPiPi::FixMappings() {
  size_t pos=0, neu=0, neg=0;
  for (size_t i=0;i<3;i++) {
    if (m_flavs[p_i[i]].Charge()>0)  pos++;
    if (m_flavs[p_i[i]].Charge()==0) neu++;
    if (m_flavs[p_i[i]].Charge()<0)  neg++;
  }
  int total = pos-neg;
  m_i = m_j = m_k = 3;
  for (size_t i=0;i<3;i++) {
    int q = m_flavs[p_i[i]].Charge();
    if (total==1 && pos==2 && neg==1) {
      if (q<0)         m_j = i;
      else if (m_i==3) m_i = i;
      else             m_k = i;
    }
    else if (total==-1 && neg==2 && pos==1) {
      if (q>0)         m_j = i;
      else if (m_i==3) m_i = i;
      else             m_k = i;
    }
    else if (total==1 && pos==1 && neu==2) {
      if (q>0)         m_j = i;
      else if (m_i==3) m_i = i;
      else             m_k = i;
    }
    else if (total==-1 && neg==1 && neu==2) {
      if (q<0)         m_j = i;
      else if (m_i==3) m_i = i;
      else             m_k = i;
    }
  }
  msg_Out()<<METHOD<<"(tot = "<<total<<" from "<<pos<<"/"<<neu<<"/"<<neg<<") "
	   <<"makes two resonances: "
	   <<"["<<m_i<<m_j<<"] and ["<<m_k<<m_j<<"]\n";
}

void VA_0_PiPiPi::FixNorm(struct GeneralModel & model) {
  double iso = 0., CKM = 0.;
  switch (int(m_PSmode)) {
  case int(PSmode::pipipi_pmp):
    iso = 1./sqrt(2.);
    CKM = model("Vud", Tools::Vud);
    break;
  case int(PSmode::pipipi_p00):
    iso = 1./sqrt(2.);
    CKM = model("Vud", Tools::Vud);
    break;
  case int(PSmode::Kpipi_pmp):
    iso = 1./sqrt(2.);
    CKM = model("Vus", Tools::Vus);
    break;
  case int(PSmode::Kpipi_p00):
    iso = 1./sqrt(2.);
    CKM = model("Vus", Tools::Vus);
    break;
  case int(PSmode::piKpi_p00):
    iso = 1./sqrt(2.);
    CKM = model("Vus", Tools::Vus);
    break;
  default: break;
  }
  m_norm = iso * CKM / sqrt(0.5);
  if (m_norm<=0.) THROW(fatal_error,"Current with zero norm.");
}

bool VA_0_PiPiPi::SelectResonances(struct GeneralModel & model) {
  msg_Out()<<METHOD<<"(mode = "<<int(m_PSmode)<<")\n";
  vector<Flavour> ijk, ij,kj;
  for (size_t i=0;i<3;i++) ijk.push_back(m_flavs[p_i[i]]);
  ij.push_back(m_flavs[p_i[m_i]]);ij.push_back(m_flavs[p_i[m_j]]);
  kj.push_back(m_flavs[p_i[m_k]]);kj.push_back(m_flavs[p_i[m_j]]);
  msg_Out()<<"ijk : "<<ijk<<"\n"<<"ij  : "<<ij<<"\n"<<"kj  : "<<kj<<"\n";
  if (m_PSmode==PSmode::pipipi_pmp || m_PSmode==PSmode::pipipi_p00) {
    // a bit of a cheat - the resonances should be neutral or charged,
    // but they are the same mass etc, hence I treat them as exchangeable.
    MakeTree(model,fftag::ffV_1,
	     Flavour(kf_a_1_1260_plus),ijk,
	     { Flavour(kf_rho_770_plus),
	       Flavour(kf_rho_1450_plus),
	       Flavour(kf_rho_1700_plus) },
	     ij);
    MakeTree(model,fftag::ffV_2,
	     Flavour(kf_a_1_1260_plus),ijk,
	     { Flavour(kf_rho_770_plus),
	       Flavour(kf_rho_1450_plus),
	       Flavour(kf_rho_1700_plus) },
	     kj);
    return true;
  }
  if (m_PSmode==PSmode::Kpipi_pmp || m_PSmode==PSmode::Kpipi_p00 ||
      m_PSmode==PSmode::piKpi_ppm || m_PSmode==PSmode::piKpi_p00) {
    return true;
  }
  return false;
}

void VA_0_PiPiPi::MakeTree(struct GeneralModel & model,fftag ff,
			   const Flavour & start,const vector<Flavour> & ijk,
			   const list<Flavour> & rlist,const vector<Flavour> & ij) {
  Resonance_Base * R = MakeResonance(model,start,ijk);
  if (R) {
    msg_Out()<<METHOD<<" inits new tree for "<<int(ff)<<" with "<<start<<".\n";
    m_treemap[ff] = new Resonance_Tree(R,R->Weight());
    for (list<Flavour>::const_iterator fit=rlist.begin();fit!=rlist.end();fit++) {
      R = MakeResonance(model,(*fit),ij);
      if (R) {
	msg_Out()<<"    adding new resonance: "<<(*fit)<<", "
		 <<R->Mass()<<" +/- "<<R->Width()
		 <<" with wt = "<<R->Weight()<<"\n";
	m_treemap[ff]->Add(R);
      }
    }
  }
}

Resonance_Base * VA_0_PiPiPi::
MakeResonance(struct GeneralModel & model,const Flavour & flav,const vector<Flavour> & ij) {
  msg_Out()<<METHOD<<"(mode = "<<int(m_PSmode)<<" for "<<flav<<")\n";
  vector<double> defaults;
  if (!FillDefaults(flav.Kfcode(),defaults)) return NULL;
  Res_Params params(flav,ij);
  string id = string("");
  switch (flav.Kfcode()) {
  case kf_a_1_1260_plus:      id = string("a_1(1260)+");  break;
  case kf_K_1_1400_plus:      id = string("K_1(1400)+");  break;
  case kf_K_1_1270_plus:      id = string("K_1(127)+");  break;
  case kf_rho_1700_plus:      id = string("rho(1700)+");  break;
  case kf_rho_1450_plus:      id = string("rho(1450)+");  break;
  case kf_rho_770_plus:       id = string("rho(770)+");   break;
  case kf_rho_1700:           id = string("rho(1700)");  break;
  case kf_rho_1450:           id = string("rho(1450)");  break;
  case kf_rho_770:            id = string("rho(770)");   break;
  case kf_K_star_1410_plus:   id = string("K*(1410)+");   break;
  case kf_K_star_892_plus:    id = string("K*(892)+");    break;
  case kf_K_star_1410:        id = string("K*(1410)");   break;
  case kf_K_star_892:         id = string("K*(892)");    break;
  case kf_K_0_star_1430_plus: id = string("K*_0(1430)+"); break;
  case kf_K_0_star_1430:      id = string("K*_0(1430)"); break;
  case kf_a_0_980_plus:       id = string("a_0(980)+");   break;
  case kf_a_0_1450_plus:      id = string("a_0(1450)+");  break;
  default: break;
  }
  params.m_weight  = model(string("Weight_")+id, defaults[0]);
  msg_Out()<<METHOD<<": "<<defaults[0]<<", "<<defaults[1]<<", "
	   <<defaults[2]<<", "<<defaults[3]<<" -- "
	   <<params.m_outflavs[0]<<"/"<<params.m_outflavs[1]<<"\n";
  Resonance_Base * R = NULL;
  if (dabs(params.m_weight)>1.e-8) {
    params.m_OSmass  = model(string("Mass_")  +id, defaults[1]);
    params.m_OSwidth = model(string("Width_") +id, defaults[2]);
    params.m_phase   = model(string("Phase_") +id, defaults[3]);
    switch (m_restype) {
    case resonance_type::fixed:
      R = new FixedWidth_Resonance(params);
      break;
    case resonance_type::GS:
      R = new GS_Resonance(params);
      break;
    case resonance_type::running:
      if (flav.Kfcode()==kf_a_1_1260_plus ||
	  flav.Kfcode()==kf_K_1_1400_plus ||
	  flav.Kfcode()==kf_K_1_1270_plus) {
	R = new RunningWidth3_Resonance(params);
	break;
      }
      R = new RunningWidth_Resonance(params);
      break;
    default:
      break;
    }
    if (R) m_initialised[flav] = R;
  }
  return R;
}

bool VA_0_PiPiPi::FillDefaults(const kf_code & tag,vector<double> & defs) {
  // Default values as vector of 4 doubles:
  // weight, mass, width, phase
  switch (m_ffmodel) {
  case ffmodel::KS:
    return KSDefaults(tag,defs);
  }
  THROW(fatal_error,"unknown form factor model.");
}

bool VA_0_PiPiPi::KSDefaults(const kf_code & tag,vector<double> & defs) {
  string error = string("Resonance type: ")+to_string(int(m_restype));
  if (s_KSdefaults.find(int(m_restype))==
      s_KSdefaults.end())
    THROW(fatal_error,string("Run will terminate.")+error+string(" not found."));
  error += string(", decay mode: ")+to_string(int(m_PSmode));
  if (s_KSdefaults[int(m_restype)].find(int(m_PSmode))==
      s_KSdefaults[int(m_restype)].end())
    THROW(fatal_error,string("Run will terminate.")+error+string(" not found."));
  error += string(", resonance: ")+to_string((long unsigned int)(tag));
  if (s_KSdefaults[int(m_restype)][int(m_PSmode)].find(tag)==
      s_KSdefaults[int(m_restype)][int(m_PSmode)].end())
    THROW(fatal_error,string("Run will terminate.")+error+string(" not found."));
  defs = s_KSdefaults[int(m_restype)][int(m_PSmode)][tag];
  return true;
}

ResDefaults VA_0_PiPiPi::s_KSdefaults = {
  { int(resonance_type::running), {
      { int(VA_0_PiPiPi::PSmode::pipipi_pmp),
	{
	  { (long unsigned int)(kf_a_1_1260_plus),      {  1.0000, 1.2510, 0.5990, 0.0000 } },
	  { (long unsigned int)(kf_rho_1700_plus),      { -0.0370, 1.7000, 0.2350, 0.0000 } },
	  { (long unsigned int)(kf_rho_1450_plus),      { -0.1030, 1.3200, 0.3900, 0.0000 } },
	  { (long unsigned int)(kf_rho_770_plus),       {  1.0000, 0.7749, 0.1480, 0.0000 } }
	}
      },
      { int(VA_0_PiPiPi::PSmode::pipipi_p00),
	{
	  { (long unsigned int)(kf_a_1_1260_plus),      {  1.0000, 1.2510, 0.5990, 0.0000 } },
	  { (long unsigned int)(kf_rho_1700_plus),      { -0.0370, 1.7000, 0.2350, 0.0000 } },
	  { (long unsigned int)(kf_rho_1450_plus),      { -0.1030, 1.3200, 0.3900, 0.0000 } },
	  { (long unsigned int)(kf_rho_770_plus),       {  1.0000, 0.7749, 0.1480, 0.0000 } }
	}
      }
    }
  }
};

Complex VA_0_PiPiPi::FF_V1(const double & sijk,const double & sij) {
  if (m_ffmodel==ffmodel::KS) {
    if (m_treemap.find(fftag::ffV_1)!=m_treemap.end()) {
      Complex result = (*m_treemap[fftag::ffV_1])(sijk,sij);
      //msg_Out()<<METHOD<<"("<<sqrt(sijk)<<", "<<sqrt(sij)<<") = "<<result<<"\n";
      return result;
    }
  }
  return Complex(0.,0.);
}

Complex VA_0_PiPiPi::FF_V2(const double & sijk,const double & skj) {
  if (m_ffmodel==ffmodel::KS) {
    if (m_treemap.find(fftag::ffV_2)!=m_treemap.end()) {
      Complex result = (*m_treemap[fftag::ffV_2])(sijk,skj); 
      return result;
    }
  }
  return Complex(0.,0.);
}

Complex VA_0_PiPiPi::FF_A(const double & s) {
  return Complex(0.,0.);
}

Complex VA_0_PiPiPi::FF_P(const double & s) {
  return Complex(0.,0.);
}


// need to update the references once we have the tau's under control
DEFINE_CURRENT_GETTER(HADRONS::VA_0_PiPiPi,"VA_0_PiPiPi")

void ATOOLS::Getter<HADRONS::Current_Base,
		    HADRONS::ME_Parameters,HADRONS::VA_0_PiPiPi>::
PrintInfo(std::ostream &st,const size_t width) const {
  st<<"Example: $ 0 \\rightarrow \\pi \\pi \\pi $ \n\n"
    <<"Order: 0 = $\\pi^0$, 1 = $\\pi^\\pm$ \n\n"
    <<"Available form factors: \n "
    <<"  \\begin{itemize} \n"
    <<"    \\item {\\tt FORM\\_FACTOR = 1 :} Kuehn-Santamaria \n"
    <<"    \\item {\\tt FORM\\_FACTOR = 2 :} Resonance Chiral Theory \n"
    <<"  \\end{itemize} \n"
    <<"Reference: https://sherpa.hepforge.org/olddokuwiki/data/media/publications/theses/diplom\\__laubrich.pdf \n"
    <<std::endl;
}
