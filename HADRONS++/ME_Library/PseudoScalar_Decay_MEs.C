#include "HADRONS++/ME_Library/PseudoScalar_Decay_MEs.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "HELICITIES/Main/XYZFuncs.H"
#include "HELICITIES/Main/Polarization_Tools.H"

using namespace HADRONS;
using namespace ATOOLS;

void P_PP::SetModelParameters( GeneralModel md )
{
  double g   = md("g", 1.);
  double phi = md("phi", 0.);
  m_global   = g*exp(Complex(0.,1.)*phi);
}
 
void P_PP::operator()(const Vec4D * p,Spin_Amplitudes * amps,int k0_n)
{
  std::vector<std::pair<int,int> > spins;
  spins.push_back(std::pair<int,int>(p_i[0],0));
  spins.push_back(std::pair<int,int>(p_i[1],0));
  spins.push_back(std::pair<int,int>(p_i[2],0));
  amps->Add(m_global,spins);
}

DEFINE_ME_GETTER(P_PP,P_PP_Getter,"P_PP")

void P_PP_Getter::PrintInfo(std::ostream &st,const size_t width) const {
  st<<"Example: $K_S \\rightarrow \\pi\\pi$ \n\n"
    <<"Order: 0 = IS (Pseudo-)scalar, 1, 2 = FS (Pseudo-)scalar \n\n"
    <<"\\[ \\mathcal{M}=g\\exp(i\\phi) \\] \n\n"
    <<std::endl;
}





void P_LNu::SetModelParameters( GeneralModel md )
{  
  m_global   = SQRT_05*md("GF",rpa.gen.ScalarConstant(std::string("GF")))*
    md("f_P", 1.)*md("V_CKM", 1.)*(m_flavs[p_i[1]].HadMass());
}
 
void P_LNu::operator()(const Vec4D * p,Spin_Amplitudes * amps,int k0_n)
{
  Complex ampl(0.,0.);
  Complex i(0.0,1.0);
  XYZFunc F(m_n,p,m_flavs,k0_n,m_anti,p_i);
  for( int hfermion1=0; hfermion1<2; hfermion1++ ) {
    for( int hfermion2=0; hfermion2<2; hfermion2++ ) {
      ampl = m_global*F.X(1,hfermion1,p[p_i[0]],2,hfermion2,1.,-1.);
      
      std::vector<std::pair<int,int> > spins;
      spins.push_back(std::pair<int,int>(p_i[0],0));
      spins.push_back(std::pair<int,int>(p_i[1],hfermion1));
      spins.push_back(std::pair<int,int>(p_i[2],hfermion2));
      amps->Add(ampl,spins);
    }
  }
}

DEFINE_ME_GETTER(P_LNu,P_LNu_Getter,"P_LNu")

void P_LNu_Getter::PrintInfo(std::ostream &st,const size_t width) const {
  st<<"Example: $\\pi^+ \\rightarrow \\ell^+\\nu_\\ell$ \n\n"
    <<"Order: 0 = IS Pseudoscalar, 1, 2 = FS Fermion \n\n"
    <<"\\[ \\mathcal{M}=\\frac{G_F}{\\sqrt{2}} f_P m_P V_{\rm CKM}"
    <<"\\bar u_\\ell(1-\\gamma_5)u_{\\bar\\nu}\\] \n\n"
    <<std::endl;
}





void P_PV::SetModelParameters( GeneralModel md )
{
  m_global   = md("g", 1.);
  if (m_flavs[p_i[2]]==Flavour(kf_photon)) m_npol = 2;
}

void P_PV::operator()(const Vec4D * p,Spin_Amplitudes * amps,int k0_n)
{
  Complex ampl(0.,0.);
  Vec4D summom = p[p_i[0]]+p[p_i[1]];
  for( int hvector=0; hvector<m_npol; hvector++ ) {
    Vec4C eps = Polarization_Tools::ComplexBosonPolarizationVector(p[p_i[2]],hvector);
    ampl = m_global*eps*summom;
    std::vector<std::pair<int,int> > spins;
    spins.push_back(std::pair<int,int>(p_i[0],0));
    spins.push_back(std::pair<int,int>(p_i[1],0));
    spins.push_back(std::pair<int,int>(p_i[2],hvector));
    amps->Add(ampl,spins);
  }
}

DEFINE_ME_GETTER(P_PV,P_PV_Getter,"P_PV")

void P_PV_Getter::PrintInfo(std::ostream &st,const size_t width) const {
  st<<"Example: $ B^{+} \\rightarrow \\bar{D}D_{s}^{*}  $ \n\n"
    <<"Order: 0 = Scalar, 1 = Scalar, 2 = Vector \n\n"
    <<"\\[ \\mathcal{M} = g\\epsilon^\\mu_V(p_P+p_{P'})_\\mu \\]"
    <<std::endl;
}





void P_VV::SetModelParameters( GeneralModel md )
{
  m_global   = md("g", 1.);
  if (m_flavs[p_i[1]]==Flavour(kf_photon)) m_npol1 = 2;
  if (m_flavs[p_i[2]]==Flavour(kf_photon)) m_npol2 = 2;
}

void P_VV::operator()(const Vec4D * p,Spin_Amplitudes * amps,int k0_n)
{
  Complex ampl(0.,0.);
  for(int h1=0;h1<m_npol1;h1++) {
    Vec4C eps1 = Polarization_Tools::ComplexBosonPolarizationVector(p[p_i[1]],h1);
    for(int h2=0;h2<m_npol2;h2++) {
      Vec4C eps2 = Polarization_Tools::ComplexBosonPolarizationVector(p[p_i[2]],h2);
      ampl = m_global*p[p_i[0]]*cross(eps1,p[p_i[1]],eps2);
      
      std::vector<std::pair<int,int> > spins;
      spins.push_back(std::pair<int,int>(p_i[0],0));
      spins.push_back(std::pair<int,int>(p_i[1],h1));
      spins.push_back(std::pair<int,int>(p_i[2],h2));
      amps->Add(ampl,spins);
    }
  }
}

DEFINE_ME_GETTER(P_VV,P_VV_Getter,"P_VV")

void P_VV_Getter::PrintInfo(std::ostream &st,const size_t width) const {
  st<<"Example: $ B^{+} \\rightarrow D_{s1}(2536)^{+} \\bar{D}^{*}(2007) $ \n\n"
    <<"Order: 0 = Pseudoscalar, 1, 2 = Vector \n\n"
    <<"\\[ \\mathcal{M} = g\\varepsilon_{\\mu \\nu \\sigma \\rho}  "
    <<"p_{0}^{\\mu} \\epsilon^{\\nu}_{1} p^{\\sigma}_{1} \\epsilon^{\\rho}_{2} \\]"
    <<std::endl;
}






void S_VV::SetModelParameters( GeneralModel md )
{
  m_global   = md("g", 1.);
  if (m_flavs[p_i[1]]==Flavour(kf_photon)) m_npol1 = 2;
  if (m_flavs[p_i[2]]==Flavour(kf_photon)) m_npol2 = 2;
}

void S_VV::operator()(const Vec4D * p,Spin_Amplitudes * amps,int k0_n)
{
  Complex ampl(0.,0.);
  for(int h1=0;h1<m_npol1;h1++) {
    Vec4C eps1 = Polarization_Tools::ComplexBosonPolarizationVector(p[p_i[1]],h1);
    for(int h2=0;h2<m_npol2;h2++) {
      Vec4C eps2 = Polarization_Tools::ComplexBosonPolarizationVector(p[p_i[2]],h2);
      ampl = m_global*eps1*eps2;
      
      std::vector<std::pair<int,int> > spins;
      spins.push_back(std::pair<int,int>(p_i[0],0));
      spins.push_back(std::pair<int,int>(p_i[1],h1));
      spins.push_back(std::pair<int,int>(p_i[2],h2));
      amps->Add(ampl,spins);
    }
  }
}

DEFINE_ME_GETTER(S_VV,S_VV_Getter,"S_VV")

void S_VV_Getter::PrintInfo(std::ostream &st,const size_t width) const {
  st<<"Order: 0 = Scalar, 1, 2 = Vector \n\n"
    <<"\\[ \\mathcal{M} = g\\epsilon^{\\mu}_{1}\\epsilon_{\\mu,2} \\]"
    <<std::endl;
}




void P_PT::SetModelParameters( GeneralModel md )
{
  m_global   = md("g", 1.);
}

void P_PT::operator()(const Vec4D * p,Spin_Amplitudes * amps,int k0_n)
{
  Complex ampl(0.,0.);
  for( int h=0; h<5; h++ ) {
    CMatrix eps = Polarization_Tools::ComplexSpin2BosonPolarizationTensor(p[p_i[1]], h);
    ampl = m_global*p[p_i[0]]*(eps*p[p_i[0]]);

    std::vector<std::pair<int,int> > spins;
    spins.push_back(std::pair<int,int>(p_i[0],0));
    spins.push_back(std::pair<int,int>(p_i[1],h));
    spins.push_back(std::pair<int,int>(p_i[2],0));
    amps->Add(ampl,spins);
  }
}

DEFINE_ME_GETTER(P_PT,P_PT_Getter,"P_PT")

void P_PT_Getter::PrintInfo(std::ostream &st,const size_t width) const {
  st<<"Example: $ B^{+} \\rightarrow \\chi_{c2}(1P) K^{+} $ \n\n"
    <<"Order: 0 = Scalar, 1 = Tensor, 2 = Scalar \n\n"
    <<"\\[ \\mathcal{M} = \\frac{m_{1}^{2}}{m_{0}\\vec{p}_{1}^{2}} "
    <<"p_{0}^{\\mu} p_{0}^{\\nu} \\epsilon_{\\mu \\nu} \\]"
    <<std::endl;
}





void P_PPP::SetModelParameters( GeneralModel md )
{
  m_global = md("A",1.);
  m_ff     = int(md("Formfactor",0));
  switch (m_ff) {
  case 1:
    // formfactor = 1+g u + h u^2 + j v + k v^2 + f uv
    m_g    = md("g",0.);
    m_h    = md("h",0.);
    m_j    = md("j",0.);
    m_k    = md("k",0.);
    m_f    = md("f",0.);
  case 0:
  default:
    break;
  }
}

void P_PPP::operator()(const Vec4D * p,Spin_Amplitudes * amps,int k0_n)
{
  Complex ampl(m_global*Formfactor(p));
  std::vector<std::pair<int,int> > spins;
  spins.push_back(std::pair<int,int>(p_i[0],0));
  spins.push_back(std::pair<int,int>(p_i[1],0));
  spins.push_back(std::pair<int,int>(p_i[2],0));
  spins.push_back(std::pair<int,int>(p_i[3],0));
  amps->Add(ampl,spins);
}

double P_PPP::Formfactor(const Vec4D * p) {
  double s1((p[p_i[0]]-p[p_i[1]]).Abs2()),s2((p[p_i[0]]-p[p_i[2]]).Abs2()),s3((p[p_i[0]]-p[p_i[3]]).Abs2());
  double s0(1./3.*(s1+s2+s3)), mav2(1./3.*(p[p_i[1]].Abs2()+p[p_i[2]].Abs2()+p[p_i[3]].Abs2()));
  double u((s3-s0)/mav2),v((s2-s1)/mav2);
  switch (m_ff) {
  case 1: 
    return 1 + m_g*u + m_h*sqr(u) + m_j*v + m_k*sqr(v) + m_f*u*v;
  case 0: 
  default: 
    return 1.;
  }
}

DEFINE_ME_GETTER(P_PPP,P_PPP_Getter,"P_PPP")
  
void P_PPP_Getter::PrintInfo(std::ostream &st,const size_t width) const {
  st<<"Example: $ K \\rightarrow \\pi\\ell\\bar\\nu_\\ell $ with form factors\n\n"
    <<"Order: 0,1,2,3 = Pseudoscalars\n\n"
    <<"\\[ \\mathcal{M} = A*Formfactor(p). \\]"
    <<std::endl;
}





void P_VFF::SetModelParameters( GeneralModel md )
{
  m_global   = md("g", 1.);
  if (md("VDM",false)) {
    m_VDM       = true;
    m_VDM_mass  = md("mass",Flavour(kf_rho_770).HadMass());
    m_VDM_width = md("width",Flavour(kf_rho_770).Width());
  }
  if (m_flavs[p_i[1]]==Flavour(kf_photon)) m_npol = 2;
}

void P_VFF::operator()(const Vec4D * p,Spin_Amplitudes * amps,int k0_n)
{
  Complex ampl(0.,0.);
  Complex i(0.0,1.0);
  XYZFunc F(m_n,p,m_flavs,k0_n,m_anti,p_i);
  Vec4D   q(p[p_i[2]]+p[p_i[3]]);
  for( int h=0; h<m_npol; h++ ) {
    Vec4C eps = Polarization_Tools::ComplexBosonPolarizationVector(p[p_i[1]], h);
    for( int hfermion1=0; hfermion1<2; hfermion1++ ) {
      for( int hfermion2=0; hfermion2<2; hfermion2++ ) {
	ampl = m_global*F.L(3,hfermion2,2,hfermion1,1.0,1.0)*
          cross(p[p_i[1]],eps,q)/q.Abs2();
	if (m_VDM) {
	  ampl *= -(Complex(0.,1.)*m_VDM_mass*m_VDM_width)/
	    (q.Abs2()-sqr(m_VDM_mass)+Complex(0.,1.)*m_VDM_mass*m_VDM_width);
	}
	std::vector<std::pair<int,int> > spins;
	spins.push_back(std::pair<int,int>(p_i[0],0));
	spins.push_back(std::pair<int,int>(p_i[1],h));
	spins.push_back(std::pair<int,int>(p_i[2],hfermion1));
	spins.push_back(std::pair<int,int>(p_i[3],hfermion2));
	amps->Add(ampl,spins);
      }
    }
  }
}

DEFINE_ME_GETTER(P_VFF,P_VFF_Getter,"P_VFF")

void P_VFF_Getter::PrintInfo(std::ostream &st,const size_t width) const {
  st<<"Example: $ B \\rightarrow D^{*}\\ell\\bar\\nu_\\ell $ without form factors\n\n"
    <<"Order: 0 = Pseudoscalar, 1 = Vector, 2, 3 = Fermions \n\n"
    <<"\\[ \\mathcal{M} = \\frac{g}{q^2_{\\ell\\bar\\ell}}\\epsilon_{\\mu\\nu\\rho\\sigma}"
    <<"p^\\mu_V\\epsilon^\\nu_Vq^\\rho_{\\ell\\bar\\ell}"
    <<"(\\bar u_\\ell\\gamma^\\sigma u_{\\bar\\ell}) \\]"
    <<std::endl;
}





void P_PLNu::SetModelParameters( GeneralModel md )
{
  m_fP     = md("f_P", 1.);
  m_global = SQRT_05*md("GF",rpa.gen.ScalarConstant(std::string("GF")))*md("V_CKM",1.)*m_fP;
  m_ff     = int(md("Formfactor",0));
  switch (m_ff) {
  case 1: 
    m_Norm    = md("N_PP'",1.);
    m_aplus0  = md("a_+^0",-1.);
    m_aplus1  = md("a_+^1", 0.);
    m_azero0  = md("a_0^0", 0.);
    m_azero1  = md("a_0^1", 0.);
    m_masssqr_diff = sqr(m_flavs[p_i[0]].HadMass())-sqr(m_flavs[p_i[1]].HadMass());
    break;
  case 0:
  default:
    break;
  }
}

void P_PLNu::operator()(const Vec4D * p,Spin_Amplitudes * amps,int k0_n)
{
  Complex ampl(0.,0.);
  XYZFunc F(m_n,p,m_flavs,k0_n,m_anti,p_i);
  Vec4D   q(p[p_i[2]]+p[p_i[3]]);
  for( int hfermion1=0; hfermion1<2; hfermion1++ ) {
    for( int hfermion2=0; hfermion2<2; hfermion2++ ) {
      ampl = m_global*F.L(2,hfermion1,3,hfermion2,1.0,-1.0)*
	(fplus(p[p_i[0]].Abs2(),p[p_i[1]].Abs2(),q.Abs2())*(p[p_i[0]]+p[p_i[1]]-m_masssqr_diff/q.Abs2()*q)+
	 fzero(p[p_i[0]].Abs2(),p[p_i[1]].Abs2(),q.Abs2())*m_masssqr_diff/q.Abs2()*q);
      std::vector<std::pair<int,int> > spins;
      spins.push_back(std::pair<int,int>(p_i[0],0));
      spins.push_back(std::pair<int,int>(p_i[1],0));
      spins.push_back(std::pair<int,int>(p_i[2],hfermion1));
      spins.push_back(std::pair<int,int>(p_i[3],hfermion2));
      amps->Add(ampl,spins);
    }
  }
}

double P_PLNu::fplus(const double p2,const double pp2,const double q2) {
  switch (m_ff) {
  case 1:  return m_aplus0 + m_aplus1*q2/sqr(m_fP);
  case 0: 
  default: return 1.;
  }
}

double P_PLNu::fzero(const double p2,const double pp2,const double q2) {
  switch (m_ff) {
  case 1:  return m_azero0 + m_azero1*(p2-pp2)/sqr(m_fP);
  case 0: 
  default: return 0.;
  }
}

DEFINE_ME_GETTER(P_PLNu,P_PLNu_Getter,"P_PLNu")

void P_PLNu_Getter::PrintInfo(std::ostream &st,const size_t width) const {
  st<<"Example: $ K \\rightarrow \\pi\\ell\\bar\\nu_\\ell $ with form factors\n\n"
    <<"Order: 0,1 = Pseudoscalars, 2, 3 = Fermions \n\n"
    <<"\\[ \\mathcal{M} = \\frac{G_F}{\\sqrt{2}}\\frac{1}{N_{PP'}}"
    <<"\\left(f_+(q_{\\ell\\nu}^2)(p+p')_\\mu+f_-(q_{\\ell\\nu}^2)(p-p')_\\mu\\right)"
    <<"\\bar u_{\\bar\\ell}\\gamma^\\mu(1-\\gamma_5)u_\\nu\\,,\\]"<<std::endl
    <<"   where the form factors are given through the expansion"<<std::endl
    <<"\\[f_+(p^2,{p'}^2,q^2) = a_+^0 + a_+^1\\frac{q^2}{f_{P'}^2} \\]"<<std::endl
    <<"\\[f_-(p^2,{p'}^2,q^2) = a_-^0 + a_-^1\\frac{p^2-{p'}^2}{f_{P'}^2}. \\]"
    <<std::endl;
}





void P_PPLNu::SetModelParameters( GeneralModel md )
{
  m_fP     = md("f_P", 1.);
  m_global = SQRT_05*md("GF",rpa.gen.ScalarConstant(std::string("GF")))*md("V_CKM",1.)*m_fP;
  m_ff     = int(md("Formfactor",0));
  switch (m_ff) {
  case 1:
    m_ff      = 1;
    m_f1_0    = md("f1_0",1.);
    m_lambda1 = md("lambda_1",0.);
    m_f2_0    = md("f2_0",1.);
    m_lambda2 = md("lambda_2",0.);
    m_g_0     = md("g_0",1.);
    m_kappa   = md("kappa",0.);
    break;
  case 0:
  default :
    break;
  }
}

void P_PPLNu::operator()(const Vec4D * p,Spin_Amplitudes * amps,int k0_n)
{
  Complex ampl(0.,0.);
  XYZFunc F(m_n,p,m_flavs,k0_n,m_anti,p_i);
  double  q2((p[p_i[2]]+p[p_i[3]]).Abs2());
  for( int hfermion1=0; hfermion1<2; hfermion1++ ) {
    for( int hfermion2=0; hfermion2<2; hfermion2++ ) {
      ampl = m_global*F.L(3,hfermion1,4,hfermion2,1.0,-1.0)*
	(f1(0.,0.,q2)*(p[p_i[1]]+p[p_i[2]]) + f2(0.,0.,q2)*(p[p_i[1]]-p[p_i[2]]) - 
	 g(0.,0.,q2)*cross(p[p_i[0]],p[p_i[1]],p[p_i[2]]));
      std::vector<std::pair<int,int> > spins;
      spins.push_back(std::pair<int,int>(p_i[0],0));
      spins.push_back(std::pair<int,int>(p_i[1],0));
      spins.push_back(std::pair<int,int>(p_i[2],0));
      spins.push_back(std::pair<int,int>(p_i[3],hfermion1));
      spins.push_back(std::pair<int,int>(p_i[4],hfermion2));
      amps->Add(ampl,spins);
    }
  }
}

double P_PPLNu::f1(const double p2,const double pp2,const double q2) {
  switch (m_ff) {
  case 1:  return m_f1_0 + m_lambda1*q2/sqr(m_fP);
  case 0: 
  default: return 1.;
  }
}

double P_PPLNu::f2(const double p2,const double pp2,const double q2) {
  switch (m_ff) {
  case 1:  return m_f2_0 + m_lambda2*q2/sqr(m_fP);
  case 0: 
  default: return 1.;
  }
}

double P_PPLNu::g(const double p2,const double pp2,const double q2) {
  switch (m_ff) {
  case 1:  return m_g_0 + m_kappa*q2/sqr(m_fP);
  case 0: 
  default: return 1.;
  }
}

DEFINE_ME_GETTER(P_PPLNu,P_PPLNu_Getter,"P_PPLNu")

void P_PPLNu_Getter::PrintInfo(std::ostream &st,const size_t width) const {
  st<<"Example: $ K \\rightarrow \\pi\\ell\\bar\\nu_\\ell $ with form factors\n\n"
    <<"Order: 0,1 = Pseudoscalars, 2, 3 = Fermions \n\n"
    <<"\\[ \\mathcal{M} = \\frac{G_F}{\\sqrt{2}}\\frac{1}{N_{PP'}}"
    <<"\\left(f_+(q_{\\ell\\nu}^2)(p+p')_\\mu+f_-(q_{\\ell\\nu}^2)(p-p')_\\mu\\right)"
    <<"\\bar u_{\\bar\\ell}\\gamma^\\mu(1-\\gamma_5)u_\\nu\\,, \\]"<<std::endl
    <<"   where the form factors are given through the expansion"<<std::endl
    <<"\\[   f_+(p^2,{p'}^2,q^2) = a_+^0 + a_+^1\\frac{q^2}{f_{P'}^2} \\]"<<std::endl
    <<"\\[   f_-(p^2,{p'}^2,q^2) = a_-^0 + a_-^1\\frac{p^2-{p'}^2}{f_{P'}^2}. \\]"
    <<std::endl;
}





void P_FFFF::SetModelParameters( GeneralModel md )
{
  m_global   = md("g", 1.);
  if (m_flavs[p_i[1]]==m_flavs[p_i[3]]) m_exchange = true;
}

void P_FFFF::operator()(const Vec4D * p,Spin_Amplitudes * amps,int k0_n)
{
  Complex ampl(0.,0.);
  Complex i(0.0,1.0);
  XYZFunc F(m_n,p,m_flavs,k0_n,m_anti,p_i);
  Vec4D   q12(p[p_i[1]]+p[p_i[2]]);
  Vec4D   q34(p[p_i[1]]+p[p_i[2]]);
  Vec4D   q14(p[p_i[1]]+p[p_i[2]]);
  Vec4D   q32(p[p_i[1]]+p[p_i[2]]);
  for( int hfermion1=0; hfermion1<2; hfermion1++ ) {
    for( int hfermion2=0; hfermion2<2; hfermion2++ ) {
      for( int hfermion3=0; hfermion3<2; hfermion3++ ) {
	for( int hfermion4=0; hfermion4<2; hfermion4++ ) {
	  Vec4C help12 = F.L(1,hfermion1,2,hfermion2,1.0,1.0);
	  Vec4C help34 = F.L(3,hfermion3,4,hfermion4,1.0,1.0);
	  ampl = m_global*help34*cross(q12,help12,q34)/(q12.Abs2()*q34.Abs2());
	  if (m_exchange) {
	    Vec4C help14 = F.L(1,hfermion1,4,hfermion4,1.0,1.0);
	    Vec4C help32 = F.L(3,hfermion3,2,hfermion2,1.0,1.0);
	    ampl -= m_global*help32*cross(q14,help14,q32)/(q14.Abs2()*q32.Abs2());
	  }
	  std::vector<std::pair<int,int> > spins;
	  spins.push_back(std::pair<int,int>(p_i[0],0));
	  spins.push_back(std::pair<int,int>(p_i[1],hfermion1));
	  spins.push_back(std::pair<int,int>(p_i[2],hfermion2));
	  spins.push_back(std::pair<int,int>(p_i[3],hfermion3));
	  spins.push_back(std::pair<int,int>(p_i[4],hfermion4));
	  amps->Add(ampl,spins);
	}
      }
    }
  }
}

DEFINE_ME_GETTER(P_FFFF,P_FFFF_Getter,"P_FFFF")

void P_FFFF_Getter::PrintInfo(std::ostream &st,const size_t width) const {
  st<<"Example: $ \\pi^0 \\rightarrow \\ell\\bar\\ell\\ell'\\bar\\ell' $ without form factors\n\n"
    <<"Order: 0 = Pseudoscalar, 2, 3, 4, 5 = Fermions \n\n"
    <<"\\[ \\mathcal{M} = \\frac{g}{q^2_{\\ell\\bar\\ell}}\\epsilon_{\\mu\\nu\\rho\\sigma}"
    <<"q^\\mu_{\\ell\\bar\\ell}(\\bar u_\\ell\\gamma^\\nu u_{\\bar\\ell}) "
    <<"q^\\rho_{\\ell'\\bar\\ell'}(\\bar u_\\ell'\\gamma^\\sigma u_{\\bar\\ell'}) \\]"
    <<std::endl;
}
