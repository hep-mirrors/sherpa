#include "METOOLS/HadronCurrents/FormFactors/FF_0_PPP.H"
#include "METOOLS/HadronCurrents/Tools.H"

using namespace METOOLS;
using namespace ATOOLS;
using namespace std;


void FF_0_PPP_Base::FixMode() {
  if (m_flavs[m_pi[0]].Kfcode()==kf_pi &&
      m_flavs[m_pi[1]].Kfcode()==kf_pi &&
      m_flavs[m_pi[2]].Kfcode()==kf_pi_plus)      m_mode = FF_0_PPP_mode::piP_pi0_pi0;
  else if (m_flavs[m_pi[0]].Kfcode()==kf_pi_plus &&
	   m_flavs[m_pi[1]].Kfcode()==kf_pi_plus &&
	   m_flavs[m_pi[2]].Kfcode()==kf_pi_plus) m_mode = FF_0_PPP_mode::piM_piP_piP;
}

Complex FF_0_PPP_Base::operator()(const ATOOLS::Vec4D_Vector& moms) {
  Vec4D p1    = moms[m_pi[0]],  p2 = moms[m_pi[1]],   p3  = moms[m_pi[2]];
  Vec4D q     = p1+p2+p3;
  double s123 = q.Abs2(),      s12 = (p1+p2).Abs2(), s13 = (p1+p3).Abs2();
  msg_Out()<<p1<<"("<<p1.Abs2()<<") + "<<p2<<"("<<p2.Abs2()<<") + "<<p3<<"("<<p3.Abs2()<<")\n";
  switch (m_ffmodel) {
  case ff_model::none:    return Complex(1.,0.);
  case ff_model::KS:      return m_norm * FF_KS(s123,s12,s13);
    //case ff_model::RChiPT:  return m_norm * FF_RChiPT(s123,s12,s13);
  case ff_model::unknown:
  default:
    break;
  }
  return Complex(0.,0.);
}

Complex FF_0_PPP_Base::
FF_KS(const double & s123,const double & s1,const double & s2) {
  msg_Error()<<"Error in "<<METHOD<<": RChiPT not available for "
	     <<m_flavs[m_pi[0]]<<"+"<<m_flavs[m_pi[1]]<<"+"<<m_flavs[m_pi[2]]
	     <<" form factor.\n"
	     <<"   Will exit the run.\n";
  THROW(critical_error,"No form factor for 3 pseudoscalars.")
}

Complex FF_0_PPP_Base::
FF_RChiPT(const double & s123,const double & s1,const double & s2) {
  msg_Error()<<"Error in "<<METHOD<<": RChiPT not available for "
	     <<m_flavs[m_pi[0]]<<"+"<<m_flavs[m_pi[1]]<<"+"<<m_flavs[m_pi[2]]
	     <<" form factor.\n"
	     <<"   Will exit the run.\n";
  THROW(critical_error,"No form factor for 3 pseudoscalars.")
}


//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//
// {pi pi pi}^+, {K pi pi}^+, {K K pi}^+  form factors.
// Todos:
//  - add more theory - i.e. Gounaris-Sakurai forms and variations of RChiPT
//  - mirror parameters to run card/decay yaml files
//
// Form factors from:
// - KS
//     * pi pi pi from
// - none, 0: no form factor
//
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

class F1_0_PiPlusPiZeroPiZero : public FF_0_PPP_Base {
protected:
  bool    m_isF2;
  double  m_fpi;
  Complex m_alpha, m_gamma, m_delta;
  
  Summed_Propagator * p_a1s, * p_rhos;
  
  void    FixParameters(const FF_Parameters & params);
  void    Construct();
  Complex FF_KS(const double & s123,const double & s1,const double & s2);
  Complex FF_RChiPT(const double & s123,const double & s1,const double & s2);
  /*
    Complex A(const double & m2,const double & s,const double & mu2);
    double  Gamma_V(const double & s);
    double  Gamma_Vp(const double & s);
    double  Gamma_Vpp(const double & s);
  */
public :
  F1_0_PiPlusPiZeroPiZero(const FF_Parameters & params);
  ~F1_0_PiPlusPiZeroPiZero();
};

F1_0_PiPlusPiZeroPiZero::F1_0_PiPlusPiZeroPiZero(const FF_Parameters & params)  :
  FF_0_PPP_Base(params),
  p_a1s(NULL), p_rhos(NULL),
  m_isF2(false), 
  m_fpi((*params.p_model)("fpi",0.1307)/sqrt(2.))
{
  if (params.m_name=="F2_0_PPP") m_isF2 = true;
  FixParameters(params);
  Construct();
}

F1_0_PiPlusPiZeroPiZero::~F1_0_PiPlusPiZeroPiZero() {
  if (p_a1s)  { delete p_a1s;  p_a1s  = NULL; }
  if (p_rhos) { delete p_rhos; p_rhos = NULL; }
}

void F1_0_PiPlusPiZeroPiZero::FixParameters(const FF_Parameters & params)  {
  if (m_mode==FF_0_PPP_mode::piP_pi0_pi0 ||
      m_mode==FF_0_PPP_mode::piM_piP_piP) {
    m_norm    = Complex(0., -((2.*sqrt(2)*(*params.p_model)("Vud", Tools::Vud)) /
			      (3.*m_fpi) ));
    if (m_ffmodel==ff_model::KS) {
      m_gamma = Complex(-0.14500,0.0000);
      m_delta = Complex( 0.00000,0.0000);
      m_alpha = Complex( 0.00185,0.0000);
    }
  }
}

void F1_0_PiPlusPiZeroPiZero::Construct() {
  if (m_ffmodel==ff_model::KS) {
    if (m_mode==FF_0_PPP_mode::piP_pi0_pi0) {
      Propagator_Base * a11260  =
	new BreitWigner(LineShapes->Get(Flavour(kf_a_1_1260_plus)));
      Propagator_Base * rho770  =
	new BreitWigner(LineShapes->Get(Flavour(kf_rho_770_plus)));
      Propagator_Base * rho1450 =
	new BreitWigner(LineShapes->Get(Flavour(kf_rho_1450_plus)));
      p_a1s  = new Summed_Propagator();
      p_a1s->Add(a11260,  Complex(1., 0.));
      p_rhos = new Summed_Propagator();
      p_rhos->Add(rho770,  Complex(1.,0.));
      p_rhos->Add(rho1450, m_gamma);
    }
    else if (m_mode==FF_0_PPP_mode::piM_piP_piP) {
      Propagator_Base * a11260  =
	new BreitWigner(LineShapes->Get(Flavour(kf_a_1_1260_plus)));
      Propagator_Base * rho770  =
	new BreitWigner(LineShapes->Get(Flavour(kf_rho_770)));
      Propagator_Base * rho770_1  =
	new BreitWigner(LineShapes->Get(Flavour(kf_rho_770)));
      Propagator_Base * omega782  =
	new BreitWigner(LineShapes->Get(Flavour(kf_omega_782)));
      Propagator_Base * rho1450 =
	new BreitWigner(LineShapes->Get(Flavour(kf_rho_1450_plus)));
      p_a1s   = new Summed_Propagator();
      p_a1s->Add(a11260,  Complex(1., 0.));
      Summed_Propagator     * rho   = new Summed_Propagator();
      rho->Add(rho770,  Complex(1.,0.));
      Multiplied_Propagator * inter = new Multiplied_Propagator(); 
      inter->Add(rho770_1, Complex(1.,0.));
      inter->Add(omega782, Complex(1.,0.));
      rho->Add(inter,  m_alpha);
      p_rhos  = new Summed_Propagator();
      p_rhos->Add(rho,     Complex(1.,0.));
      p_rhos->Add(rho1450, m_gamma);
    }
  }
}

Complex F1_0_PiPlusPiZeroPiZero::
FF_KS(const double & s123,const double & s12,const double & s13) {
  msg_Out()<<METHOD<<"(F2 = "<<m_isF2<<", s = "<<sqrt(s123)<<", "
	   <<"s12 = "<<sqrt(s12)<<", s13 = "<<sqrt(s13)<<") ";
  if (p_a1s==NULL || p_rhos==NULL) return Complex(0.,0.);
  if (m_isF2) {
    msg_Out()<<m_norm<<" * "<<(*p_a1s)(s123)<<" * "<<(*p_rhos)(s12)<<"\n";
    return m_norm * (*p_a1s)(s123) * (*p_rhos)(s12);
  }
  msg_Out()<<m_norm<<" * "<<(*p_a1s)(s123)<<" * "<<(*p_rhos)(s13)<<"\n";
  return m_norm * (*p_a1s)(s123) * (*p_rhos)(s13);
}

Complex F1_0_PiPlusPiZeroPiZero::
FF_RChiPT(const double & s123,const double & s1,const double & s2) {
  return Complex(0.,0.);
}

//////////////////////////////////////////////////////////////////////////////////
//
// Form factors for
//   * pi pi pi vanish due to identical masses
// - Kuehn-Santamaria
//   * K pi: Z.Phys.C 72 (1996) 619 (scalar form factor)
//     (https://doi.org/10.1007/s002880050284)
//
// Todo: implement all scalar form factors here!
//
//////////////////////////////////////////////////////////////////////////////////

class F3_0_PiPlusPiZeroPiZero : public FF_0_PPP_Base {
  double  m_fpi;
  Complex m_gamma, m_delta;
  
  void    FixParameters(const FF_Parameters & params);
  void    Construct();
  Complex FF_KS(const double & s123,const double & s1,const double & s2) {
    return complex(0.,0.);
  }
  Complex FF_RChiPT(const double & s123,const double & s1,const double & s2) {
    return complex(0.,0.);
  }
public :
  F3_0_PiPlusPiZeroPiZero(const FF_Parameters & params)  :
    FF_0_PPP_Base(params),
    m_fpi((*params.p_model)("fpi",0.1307)/sqrt(2.))
  {
    msg_Out()<<METHOD<<"\n";
    FixParameters(params);
    Construct();
  }
};

void F3_0_PiPlusPiZeroPiZero::FixParameters(const FF_Parameters & params)  {}

void F3_0_PiPlusPiZeroPiZero::Construct()  {}

//////////////////////////////////////////////////////////////////////////////////
//
// Form factors for
//   * pi pi pi vanish due to identical masses
// - Kuehn-Santamaria
//   * K pi: Z.Phys.C 72 (1996) 619 (scalar form factor)
//     (https://doi.org/10.1007/s002880050284)
//
// Todo: implement all scalar form factors here!
//
//////////////////////////////////////////////////////////////////////////////////

class FS_0_PiPlusPiZeroPiZero : public FF_0_PPP_Base {
  double  m_fpi;
  Complex m_gamma, m_delta;
  
  void    FixParameters(const FF_Parameters & params);
  void    Construct();
  Complex FF_KS(const double & s123,const double & s1,const double & s2) {
    return complex(0.,0.);
  }
  Complex FF_RChiPT(const double & s123,const double & s1,const double & s2) {
    return complex(0.,0.);
  }
public :
  FS_0_PiPlusPiZeroPiZero(const FF_Parameters & params)  :
    FF_0_PPP_Base(params),
    m_fpi((*params.p_model)("fpi",0.1307)/sqrt(2.))
  {
    msg_Out()<<METHOD<<"\n";
    FixParameters(params);
    Construct();
  }
};

void FS_0_PiPlusPiZeroPiZero::FixParameters(const FF_Parameters & params)  {}

void FS_0_PiPlusPiZeroPiZero::Construct()  {}


DECLARE_FF_GETTER(FF_0_PPP_Base,"FF_0_PPP")

FormFactor_Base * ATOOLS::Getter<FormFactor_Base,FF_Parameters,
				 FF_0_PPP_Base>:: 
operator()(const METOOLS::FF_Parameters &params) const
{
  msg_Out()<<METHOD<<"("<<params.m_name<<", N_f = "<<params.m_flavs.size()<<"):\n";
  size_t Nmesons = 0;
  for (size_t i=0;i<params.m_pi.size();i++) {
    msg_Out()<<"    *  i = "<<i<<": "<<params.m_pi[i]<<"\n";
    if (params.m_flavs[params.m_pi[i]].IsMeson()) Nmesons++;
  }
  if (Nmesons!=3) return NULL;
  // Below a first round of decays/currents for which we have both
  // Kuehn-Santamaria and RChiPT parametrizations
  if (//   pi^+ pi^0 pi^0
      (params.m_flavs[params.m_pi[0]].Kfcode()==kf_pi_plus &&
       params.m_flavs[params.m_pi[1]].Kfcode()==kf_pi &&
       params.m_flavs[params.m_pi[2]].Kfcode()==kf_pi) ||
      //   pi^- pi^+ pi^+
      (params.m_flavs[params.m_pi[0]].Kfcode()==kf_pi_plus &&
       params.m_flavs[params.m_pi[1]].Kfcode()==kf_pi_plus &&
       params.m_flavs[params.m_pi[2]].Kfcode()==kf_pi_plus)      
      )
    {
    if (params.m_name=="F1_0_PPP") return new F1_0_PiPlusPiZeroPiZero(params);
    if (params.m_name=="F2_0_PPP") return new F1_0_PiPlusPiZeroPiZero(params);
    if (params.m_name=="F3_0_PPP") return new F3_0_PiPlusPiZeroPiZero(params);
    if (params.m_name=="FS_0_PPP") return new FS_0_PiPlusPiZeroPiZero(params);
  }
  return NULL;
}




