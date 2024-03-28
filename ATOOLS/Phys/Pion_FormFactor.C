#include "ATOOLS/Phys/Pion_FormFactor.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "MODEL/Main/Running_AlphaQED.H"

using namespace ATOOLS;


Pion_FormFactor::Pion_FormFactor(){
	Scoped_Settings s{ Settings::GetMainSettings()["Pion"] };
    m_form_mode = s["Form_Factor"].SetDefault(0).Get<int>();
    m_crho    = 1;
    m_crhop   = 0.14104;
    m_crhopp  = 0.0614;
    m_crhoppp = 0.0047;

    m_comega = 0.00158;
    m_cphi   = 0.00045;

    m_prho    = 0;
    m_prhop   = 3.7797;
    m_prhopp  = 1.429;
    m_prhoppp = 0.921;

    m_pomega = 0.075;
    m_pphi   = 2.888;

    m_I = Complex(0,1);
    m_flv =  Flavour(kf_pi_plus);
}

Pion_FormFactor::~Pion_FormFactor(){}


double Pion_FormFactor::ppi(const double &q2){
    return 0.5*sqrt(q2-4*sqr(m_flv.Mass()));
  }

  double Pion_FormFactor::d(const double &m){
    double L = log( (m+2.*ppi(m*m))/2./m_flv.Mass());
    return 3./M_PI*sqr(m_flv.Mass()/ppi(m*m))*L+0.5*m/ppi(m*m)/M_PI-m*sqr(m_flv.Mass())/M_PI/pow(ppi(m*m),3);
  }

  double Pion_FormFactor::h(const double &q2){
    double L = log( (sqrt(q2) + 2.*ppi(q2) )/2./m_flv.Mass());
    return 2./M_PI*ppi(q2)/sqrt(q2)*L;
  }

  double Pion_FormFactor::dh(const double &q2){
    return h(q2)/8*(1./sqr(ppi(q2))-4./q2)+0.5/M_PI/q2;
  }

  double Pion_FormFactor::f(const double &q2, const double &m, const double &g){
    double t1 = sqr(ppi(q2))*(h(q2)-h(m*m))+sqr(ppi(m*m))*(m*m-q2)*dh(m*m);
    return t1*g*m*m/pow(ppi(m*m),3);
  }

  double Pion_FormFactor::gamma(const double &q2, const double &m, const double &g){
    return g*m/sqrt(q2)*pow(ppi(q2)/ppi(m*m),3);
  }

  Complex Pion_FormFactor::b(const double &q2, const Flavour &fl){
    double m = fl.Mass();
    double w = fl.Width();
    Complex I(0., 1.);

    return m*m/(m*m-q2-m*w*I);
  }

  Complex Pion_FormFactor::GS(const double &q2, const Flavour &fl){
    //  Gounaris-Sakurai (GS) function https://inspirehep.net/literature/53152
    double mf = fl.Mass();
    double wf = fl.Width();
    Complex I(0, 1.);

    Complex num = mf*mf+d(mf)*mf*wf;
    Complex den = mf*mf-q2+f(q2,mf,wf)-I*mf*gamma(q2,mf,wf);
    if(IsBad(num) || IsBad(den)){
      PRINT_VAR(fl);
      msg_Error()<<"NaN in "<<METHOD<<std::endl
                 <<"num = "<<num<<std::endl
                 <<"den = "<<num<<std::endl
                 <<"mf = "<<mf<<std::endl
                 <<"wf = "<<wf<<std::endl
                 <<"d(mf) = "<<d(mf)<<std::endl
                 <<"f(q2,mf,wf) = "<<f(q2,mf,wf)<<std::endl
                 <<"gamma(q2,mf,wf) = "<<gamma(q2,mf,wf)<<std::endl;
    }
    return num/den;
  }

  double Pion_FormFactor::Eval(const double &q2){
    if(m_form_mode==0) return 1;
    // https://gitlab.com/strong2020/monte-carlo-results/-/tree/root/pion-formfactor?ref_type=heads
    Particle_Info rhoi(kf_rho_770,  0.77456, 0, 0.14832,3,0,0,0,"rho","rho");
    
    Particle_Info rhopi(kf_rho_770, 1.4859,  0,  0.37360, -3,0,0,1,"rho pp","rho");
    
    Particle_Info rhoppi(kf_rho_770, 1.8668,  0, 0.30334, -3,0,0,1,"rho ppp","rho");
    
    Particle_Info rhopppi(kf_rho_770, 2.2645,  0, 0.11327, -3,0,0,1,"rho ppp","rho");
    
    Particle_Info omegapi(kf_omega_782, 0.78248, 0, 0.00855, -3,0,0,1,"omega","omega");
    
    Particle_Info phipi(kf_phi_1020, 1.01947, 0, 0.00425, -3,0,0,1,"phi","phi");
    
    m_rho   = Flavour(rhoi);
    m_rhop  = Flavour(rhopi);
    m_rhopp = Flavour(rhoppi);
    m_rhoppp = Flavour(rhopppi);

    m_omega = Flavour(omegapi);
    m_phi   = Flavour(phipi);

    Complex num1 = 1.+( q2/sqr(m_omega.Mass())*m_comega*exp(m_I*m_pomega))*b(q2,m_omega);
    num1 += q2/(sqr(m_phi.Mass()))*m_cphi*exp(m_I*m_pphi)*b(q2,m_phi);
    num1 *= GS(q2,m_rho);
    Complex den = 1.+m_crhop*exp(m_I*m_prhop)+m_crhopp*exp(m_I*m_prhopp)+m_crhoppp*exp(m_I*m_prhoppp);

    Complex num2 = GS(q2,m_rhop) * m_crhop*exp(m_I*m_prhop);
    num2 += GS(q2,m_rhopp) * m_crhopp*exp(m_I*m_prhopp);
    num2 += GS(q2,m_rhoppp) * m_crhoppp*exp(m_I*m_prhoppp);
    if(IsBad(num1) || IsBad(num2) || IsBad(den)){
      msg_Error()<<"NaN in pion form-factor\n"
                 <<"Numertor 1 = "<<num1<<std::endl
                 <<"Numertor 2 = "<<num2<<std::endl
                 <<"Den 1 = "<<den<<std::endl;
    }
    Complex Form = (num1+num2)/den;
    return (Form*conj(Form)).real();
  }
