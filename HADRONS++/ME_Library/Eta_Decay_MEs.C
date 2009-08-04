#include "HADRONS++/ME_Library/Eta_Decay_MEs.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "METOOLS/Main/Polarization_Tools.H"
#include "MODEL/Main/Model_Base.H"

using namespace HADRONS;
using namespace ATOOLS;
using namespace METOOLS;
using namespace MODEL;

void Eta_PPV::SetModelParameters(GeneralModel md)
{
  if (m_flavs[p_i[3]]==Flavour(kf_photon)) m_npol = 2;
  m_ff             = int(md("Formfactor",0));
  m_fP             = md("f_pi",0.130)/sqrt(2.);  
  // extra factor in pion decay constant due to convention in hep-ph/0112150
  double f8_by_fpi = md("f_8/f_pi",1.30); 
  double f0_by_fpi = md("f_8/f_pi",1.04);
  double theta     = md("Theta",-20./180.*M_PI);
  double e(sqrt(4.*M_PI*s_model->ScalarConstant(std::string("alpha_QED(0)"))));
  m_global = 3.*e/(12.*sqrt(3.)*sqr(M_PI)*pow(m_fP,3));
  if (m_flavs[p_i[0]]==Flavour(kf_eta)) 
    m_global *= (1./f8_by_fpi*cos(theta)-sqrt(2.)/f0_by_fpi*sin(theta));
  else if (m_flavs[p_i[0]]==Flavour(kf_eta_prime_958)) 
    m_global *= (1./f8_by_fpi*sin(theta)+sqrt(2.)/f0_by_fpi*cos(theta));

  m_VDM_mass  = md("M_Rho",Flavour(kf_rho_770).HadMass());
  m_VDM_width = md("Gamma_Rho",Flavour(kf_rho_770).Width());

  switch (m_ff) {
  case 2:
    // Omnes-type formfactor
    m_mpipi2 = sqr(m_flavs[p_i[1]].HadMass()+m_flavs[p_i[2]].HadMass());
    m_mrho2  = sqr(m_VDM_mass);
    m_pref   = m_mpipi2/(96.*sqr(M_PI*m_fP));
    break;
  case 1:
  case 0:
  default:
    break;
  }
}
 
void Eta_PPV::operator()(const Vec4D * p,Spin_Amplitudes * amps)
{
  Complex ampl(0.,0.),pref(m_global*Formfactor((p[p_i[1]]+p[p_i[2]]).Abs2()));
  Polarization_Vector pol(p[p_i[3]], sqr(m_flavs[p_i[3]].HadMass()));
  for (int hvector=0;hvector<m_npol;hvector++) {
    Vec4C eps = pol[hvector];
    ampl = pref*eps*cross(p[p_i[1]],p[p_i[2]],p[p_i[3]]);
    std::vector<std::pair<int,int> > spins;
    spins.push_back(std::pair<int,int>(p_i[0],0));
    spins.push_back(std::pair<int,int>(p_i[1],0));
    spins.push_back(std::pair<int,int>(p_i[2],0));
    spins.push_back(std::pair<int,int>(p_i[3],hvector));
    amps->Add(ampl,spins);
  }
}

Complex Eta_PPV::Formfactor(const double s) {
  Complex i(Complex(0.,1.));
  double runwidth = 
    (pow(lambdaNorm(sqrt(s),m_flavs[p_i[1]].HadMass(),m_flavs[p_i[2]].HadMass()),3.)*m_mrho2)/
    (pow(lambdaNorm(m_VDM_mass,m_flavs[p_i[1]].HadMass(),m_flavs[p_i[2]].HadMass()),3.)*s)*m_VDM_width;
  switch (m_ff) {
  case 2: // Omnes form
    return Omnes_Formfactor(s,runwidth);
  case 1: //VDM model
    return 1.-1.5*s/(s-m_mrho2+i*m_VDM_mass*runwidth);
  case 0:
  default:
    return 1.;
  }
}

Complex Eta_PPV::Omnes_Formfactor(const double s,const double runwidth) {
  double  c(1.),pabs;
  Complex i(Complex(0.,1.)),ff;
  if(s>m_mpipi2) {
    pabs = sqrt(1.-m_mpipi2/s);
    ff   = (1.-s/m_mpipi2)*pabs*log((1.+pabs)/(1.-pabs))-2.;
    //ff  += i*runwidth;
  }
  else {
    pabs = sqrt(m_mpipi2/s-1.);
    ff   = 2.*(1.-s/m_mpipi2)*pabs*atan(1./pabs)-2.;
  }
  Complex D = 1.-s/m_mrho2-s/(96.*sqr(M_PI*m_fP))*log(4.*m_mrho2/m_mpipi2)-m_pref*ff;

  //std::cout<<"Check this "<<D<<"("<<m_pref*ff<<"), "<<((1.+s/(2.*m_mrho2))/D)<<" vs. "
  //	   <<(1.-1.5*s/(s-m_mrho2+i*m_VDM_mass*runwidth))<<" for "<<s<<std::endl;
  return 1.-c+c*(1.+s/(2.*(m_mrho2-i*m_VDM_mass*m_VDM_width)))/D;
}

DEFINE_ME_GETTER(Eta_PPV,Eta_PPV_Getter,"Eta_PPV")

void Eta_PPV_Getter::PrintInfo(std::ostream &st,const size_t width) const {
  st<<"Example: $\\eta \\rightarrow \\pi\\pi\\gamma$ \n\n"
    <<"Order: 0 = $\\eta$, 1, 2 = $\\pi$, 3 = $\\gamma$ \n\n"
    <<"\\[ \\mathcal{M}=gB(s,t,u)\\epsilon_{\\mu\\nu\\rho\\sigma}"
    <<"\\epsilon^\\mu p_+^\\nu p_-^\\rho k_\\gamma^\\sigma\\] \n\n"
    <<std::endl;
}





void Eta_PVV::SetModelParameters(GeneralModel md)
{
  if (m_flavs[p_i[2]]==Flavour(kf_photon)) m_npol1 = 2;
  if (m_flavs[p_i[3]]==Flavour(kf_photon)) m_npol2 = 2;
  m_ff             = int(md("Formfactor",0));
  m_fP             = md("f_pi",0.130)/sqrt(2.);  
  // extra factor in pion decay constant due to convention in hep-ph/0112150
  double f8_by_fpi = md("f_8/f_pi",1.30); 
  double f0_by_fpi = md("f_8/f_pi",1.04);
  double theta     = md("Theta",-20./180.*M_PI);

  m_mrho           = md("M_Rho",Flavour(kf_rho_770).HadMass());
  m_Grho           = md("Gamma_Rho",Flavour(kf_rho_770).Width());
  m_momega         = md("M_Rho",Flavour(kf_omega_782).HadMass());
  m_Gomega         = md("Gamma_Rho",Flavour(kf_omega_782).Width());
  m_mrho2          = sqr(m_mrho);
  m_momega2        = sqr(m_momega);
  double alpha     = s_model->ScalarConstant(std::string("alpha_QED(0)"));
  double g2        = m_mrho2/(2.*sqr(m_fP)); 
  m_g_rhoG         = 2.*sqrt(4.*M_PI*alpha*g2)*sqr(m_fP);

  double help      = 3.*m_momega2/(64.*pow(M_PI,3)*m_fP*alpha);
  double g_omrhopi = md("g_omrhopi",help);

  m_global         = 2.*sqrt(3.)/9.*sqr(g_omrhopi);
  if (m_ff==0) m_global/=m_mrho2;

  std::cout<<"g = "<<g_omrhopi<<" ("<<m_momega2<<" "<<m_fP<<")"<<std::endl;

  if (m_flavs[p_i[0]]==Flavour(kf_eta)) 
    m_global *= (1./f8_by_fpi*cos(theta)-sqrt(2.)/f0_by_fpi*sin(theta));
  else if (m_flavs[p_i[0]]==Flavour(kf_eta_prime_958)) 
    m_global *= (1./f8_by_fpi*sin(theta)+sqrt(2.)/f0_by_fpi*cos(theta));
}
 
void Eta_PVV::operator()(const Vec4D * p,Spin_Amplitudes * amps)
{
  Complex ampl(0.,0.);
  Complex FormD = D(p), FormE = E(p);

  //std::cout<<" pref,D,E = "<<m_global<<" "<<FormD<<" "<<FormE<<std::endl;
  Polarization_Vector pol1(p[p_i[2]], sqr(m_flavs[p_i[2]].HadMass()));
  Polarization_Vector pol2(p[p_i[3]], sqr(m_flavs[p_i[3]].HadMass()));
  for (int hvector1=0;hvector1<m_npol1;hvector1++) {
    Vec4C eps1 = pol1[hvector1];
    for (int hvector2=0;hvector2<m_npol2;hvector2++) {
      Vec4C eps2 = pol2[hvector2];

      Complex eps12(eps1*eps2),q12(p[p_i[2]]*p[p_i[3]]),eps1q2(eps1*p[p_i[3]]),eps2q1(eps2*p[p_i[2]]);
      Complex pq1(p[p_i[0]]*p[p_i[2]]),pq2(p[p_i[0]]*p[p_i[3]]),peps1(p[p_i[0]]*eps1),peps2(p[p_i[0]]*eps2);

      ampl = m_global* (FormD * (eps12*q12-eps1q2*eps2q1) +
			FormE * (eps12*pq1*pq2 + peps1*peps2*q12 - 
				 eps1q2*peps2*pq1 - eps2q1*peps1*pq2));
      std::vector<std::pair<int,int> > spins;
      spins.push_back(std::pair<int,int>(p_i[0],0));
      spins.push_back(std::pair<int,int>(p_i[1],0));
      spins.push_back(std::pair<int,int>(p_i[2],hvector1));
      spins.push_back(std::pair<int,int>(p_i[3],hvector2));
      amps->Add(ampl,spins);
    }
  }
}

Complex Eta_PVV::D(const Vec4D * p) {
  Complex i(Complex(0.,1.)),help(1.,0.);
  double meta2(p[p_i[0]].Abs2()),p02(p[p_i[0]]*p[p_i[2]]),p03(p[p_i[0]]*p[p_i[3]]),
    t((p[p_i[1]]+p[p_i[3]]).Abs2()),u((p[p_i[1]]+p[p_i[2]]).Abs2());
  double runwidth_rho_t = 
    (pow(lambdaNorm(sqrt(t),m_flavs[p_i[1]].HadMass(),m_flavs[p_i[2]].HadMass()),3.)*m_mrho2)/
    (pow(lambdaNorm(m_mrho,m_flavs[p_i[1]].HadMass(),m_flavs[p_i[2]].HadMass()),3.)*t)*m_Grho;
  double runwidth_rho_u = 
    (pow(lambdaNorm(sqrt(u),m_flavs[p_i[1]].HadMass(),m_flavs[p_i[3]].HadMass()),3.)*m_mrho2)/
    (pow(lambdaNorm(m_mrho,m_flavs[p_i[1]].HadMass(),m_flavs[p_i[3]].HadMass()),3.)*u)*m_Grho;
  switch (m_ff) {
  case 1: //VDM model
    help = 
      sqr(m_g_rhoG/m_mrho2)*((meta2-p03)/(t-m_mrho2+i*m_mrho*runwidth_rho_t)+
			     (meta2-p02)/(u-m_mrho2+i*m_mrho*runwidth_rho_u))+
      sqr(m_g_rhoG/m_momega2)*(1./3.)*((meta2-p03)/(t-m_momega2+i*m_momega*m_Gomega)+
				       (meta2-p02)/(u-m_momega2+i*m_momega*m_Gomega));
    break;
  case 0:
  default:
    break;
  }
  return help;
}

Complex Eta_PVV::E(const Vec4D * p) {
  Complex i(Complex(0.,1.)),help(1.,0.);
  double  t((p[p_i[1]]+p[p_i[3]]).Abs2()),u((p[p_i[1]]+p[p_i[2]]).Abs2());
  double runwidth_rho_t = 
    (pow(lambdaNorm(sqrt(t),m_flavs[p_i[1]].HadMass(),m_flavs[p_i[2]].HadMass()),3.)*m_mrho2)/
    (pow(lambdaNorm(m_mrho,m_flavs[p_i[1]].HadMass(),m_flavs[p_i[2]].HadMass()),3.)*t)*m_Grho;
  double runwidth_rho_u = 
    (pow(lambdaNorm(sqrt(u),m_flavs[p_i[1]].HadMass(),m_flavs[p_i[3]].HadMass()),3.)*m_mrho2)/
    (pow(lambdaNorm(m_mrho,m_flavs[p_i[1]].HadMass(),m_flavs[p_i[3]].HadMass()),3.)*u)*m_Grho;
  switch (m_ff) {
  case 1: //VDM model
    help = 
      m_g_rhoG/m_mrho2*(1./(t-m_mrho2+i*m_mrho*runwidth_rho_t)+
			1./(u-m_mrho2+i*m_mrho*runwidth_rho_u))+
      m_g_rhoG/m_momega2*(1./3.)*(1./(t-m_momega2+i*m_momega*m_Gomega)+
				  1./(u-m_momega2+i*m_momega*m_Gomega));
    break;
  case 0:
  default:
    break;
  }
  return help;
}

DEFINE_ME_GETTER(Eta_PVV,Eta_PVV_Getter,"Eta_PVV")

void Eta_PVV_Getter::PrintInfo(std::ostream &st,const size_t width) const {
  st<<"Example: $\\eta \\rightarrow \\pi\\pi\\gamma$ \n\n"
    <<"Order: 0 = $\\eta$, 1, 2 = $\\pi$, 3 = $\\gamma$ \n\n"
    <<"\\[ \\mathcal{M}=gB(s,t,u)\\epsilon_{\\mu\\nu\\rho\\sigma}"
    <<"\\epsilon^\\mu p_+^\\nu p_-^\\rho k_\\gamma^\\sigma\\] \n\n"
    <<std::endl;
}






void Eta_PPP::SetModelParameters(GeneralModel md)
{
}
 
void Eta_PPP::operator()(const Vec4D * p,Spin_Amplitudes * amps)
{
}

DEFINE_ME_GETTER(Eta_PPP,Eta_PPP_Getter,"Eta_PPP")

void Eta_PPP_Getter::PrintInfo(std::ostream &st,const size_t width) const {
  st<<"Example: $\\eta \\rightarrow \\pi\\pi\\pi$ \n\n"
    <<"Order: 0 = $\\eta$, 1, 2, 3 = $\\pi^{+,-,0}$\n\n"
    <<std::endl;
}

