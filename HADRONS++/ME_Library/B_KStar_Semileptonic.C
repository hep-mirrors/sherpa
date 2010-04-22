#include "HADRONS++/ME_Library/B_KStar_Semileptonic.H"
#include "ATOOLS/Org/Message.H"
#include "METOOLS/Main/XYZFuncs.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Math/MathTools.H"
#include "METOOLS/Main/Polarization_Tools.H"

using namespace HADRONS;
using namespace ATOOLS;
using namespace METOOLS;
using namespace MODEL;
using namespace std;

void B_KStar_Semileptonic::SetModelParameters( GeneralModel _md )
{
  m_mB     = p_masses[0];
  m_mKhat  = p_masses[1]/m_mB;
  double Vts  = _md("Vts", Tools::Vts);
  m_Vtb  = _md("Vtb", Tools::Vtb);
  m_Vub = _md("Vub", Tools::Vub);
  double GF   = _md("GF", s_model->ScalarConstant(string("GF")) );
  double alpha   = _md("alpha", s_model->ScalarConstant(string("alpha_QED(0)")) );
  m_global = GF * alpha * Vts * m_Vtb * m_mB / 2.0 / sqrt(2.0) / M_PI;
  m_LD    = bool(int(_md("LD",1.0)+0.5));
  m_C1    = _md("C1",   -0.248);
  m_C2    = _md("C2",    1.107);
  m_C3    = _md("C3",    0.011);
  m_C4    = _md("C4",   -0.026);
  m_C5    = _md("C5",    0.007);
  m_C6    = _md("C6",   -0.031);
  m_C7eff = _md("C7eff",-0.313);
  m_C9    = _md("C9",    4.344);
  m_C10   = _md("C10",  -4.669);

  m_A1_0  = _md("A1_0",0.337);
  m_A2_0  = _md("A2_0",0.282);
  m_A0_0  = _md("A0_0",0.471);
  m_V_0   = _md("V_0",0.457);
  m_T1_0  = _md("T1_0",0.379);
  m_T2_0  = _md("T2_0",0.379);
  m_T3_0  = _md("T3_0",0.260);
  m_A1_c1 = _md("A1_c1",0.602);
  m_A2_c1 = _md("A2_c1",1.172);
  m_A0_c1 = _md("A0_c1",1.505);
  m_V_c1  = _md("V_c1",1.482);
  m_T1_c1 = _md("T1_c1",1.519);
  m_T2_c1 = _md("T2_c1",0.517);
  m_T3_c1 = _md("T3_c1",1.129);
  m_A1_c2 = _md("A1_c2",0.258);
  m_A2_c2 = _md("A2_c2",0.567);
  m_A0_c2 = _md("A0_c2",0.710);
  m_V_c2  = _md("V_c2",1.015);
  m_T1_c2 = _md("T1_c2",1.030);
  m_T2_c2 = _md("T2_c2",0.426);
  m_T3_c2 = _md("T3_c2",1.128);
  m_A1_c3 = _md("A1_c3",0.0);
  m_A2_c3 = _md("A2_c3",0.0);
  m_A0_c3 = _md("A0_c3",0.0);
  m_V_c3  = _md("V_c3",0.0);
  m_T1_c3 = _md("T1_c3",0.0);
  m_T2_c3 = _md("T2_c3",0.0);
  m_T3_c3 = _md("T3_c3",0.0);
  m_mc  = _md("mc",1.4);
  m_ms  = _md("ms",0.2); 
}

void B_KStar_Semileptonic::operator()( 
    const Vec4D             * _p,
    METOOLS::Spin_Amplitudes * amps)
{
  double s     = (_p[p_i[2]]+_p[p_i[3]]).Abs2();
  double shat  = s/(m_mB*m_mB);
  Vec4D pBhat = _p[p_i[0]]/(m_mB);
  Vec4D pKhat = _p[p_i[1]]/(m_mB);
  Vec4D phat  = pBhat + pKhat;
  Vec4D qhat  = pBhat - pKhat;
  
  // get m_b^hat at scale mu
  double mbpole = 4.8;
  double mu = mbpole; // fixme: correct scale?
  double alphas = s_model->ScalarFunction("alpha_S",mu);
  double m_mb = mbpole*(1.0-4.0/3.0*alphas/M_PI);
  double m_mbhat = m_mb/m_mB;

  // calculate formfactors
  double shat2=shat*shat;
  double shat3=shat2*shat;
  double A1 = m_A1_0 * exp(m_A1_c1*shat + m_A1_c2*shat2 + m_A1_c3*shat3);
  double A2 = m_A2_0 * exp(m_A2_c1*shat + m_A2_c2*shat2 + m_A2_c3*shat3);
  double A0 = m_A0_0 * exp(m_A0_c1*shat + m_A0_c2*shat2 + m_A0_c3*shat3);
  double V  = m_V_0  * exp(m_V_c1 *shat + m_V_c2 *shat2 + m_V_c3*shat3);
  double T1 = m_T1_0 * exp(m_T1_c1*shat + m_T1_c2*shat2 + m_T1_c3*shat3);
  double T2 = m_T2_0 * exp(m_T2_c1*shat + m_T2_c2*shat2 + m_T2_c3*shat3);
  double T3 = m_T3_0 * exp(m_T3_c1*shat + m_T3_c2*shat2 + m_T3_c3*shat3);

  // calculate coefficients
  Flavour cquark(kf_c); Flavour squark(kf_s);
  Complex C9eff;
  if(m_LD==true) C9eff = C9sehgal(shat) + sehgalld(shat);
  else C9eff = m_C9 + gSD(m_mc/p_masses[0],shat) * (3.*m_C1+m_C2+3.*m_C3+m_C4+3.*m_C5+m_C6)
		- 0.5 * gSD(m_mb/p_masses[0],shat) * (4.*m_C3+4.*m_C4+3.*m_C5+m_C6)
		- 0.5 * gSD(m_ms/p_masses[0],shat) * (m_C3+3.*m_C4)
		+ (2./9.) * (3.*m_C3+m_C4+3.*m_C5+m_C6);
  Complex A = 2.0/(1.0+m_mKhat)*C9eff*V + 4.0*m_mbhat/shat*m_C7eff*T1;
  Complex B = (1.0+m_mKhat) * ( C9eff*A1+2.0*m_mbhat/shat*(1.0-m_mKhat)*m_C7eff*T2 ); 
  Complex T3T2term = T3+(1.0-m_mKhat*m_mKhat)/shat*T2;
  Complex C = 1.0/(1.0-m_mKhat*m_mKhat) * ( (1.0-m_mKhat)*C9eff*A2+2.0*m_mbhat*m_C7eff*T3T2term );
  Complex A1A2A0term = (1.0+m_mKhat)*A1 - (1.0-m_mKhat)*A2 - 2.0*m_mKhat*A0;
  Complex D = 1.0/shat * ( C9eff*A1A2A0term - 2.0*m_mbhat*m_C7eff*T3 );
  Complex E = 2.0/(1.0+m_mKhat) * m_C10 * V;
  Complex F = (1.0+m_mKhat) * m_C10 * A1;
  Complex G = 1.0/(1.0+m_mKhat)*m_C10*A2;
  Complex H = 1.0/shat*m_C10 * ( (1.0+m_mKhat)*A1-(1.0-m_mKhat)*A2-2.0*m_mKhat*A0 );  

  XYZFunc Func(m_n,_p,m_flavs,Tools::k0,m_anti,p_i);
  Complex i(0.0,1.0);
  Polarization_Vector pol(_p[p_i[1]], sqr(m_flavs[p_i[1]].HadMass()));
  for(int hhad = 0; hhad <3; hhad++) {
    Vec4C eps = conj(pol[hhad]);
    Vec4C T1vec = A*cross(eps,pBhat,pKhat)
                       - i*B*eps
                       + i*(eps*pBhat)*(C*phat+D*qhat);
    Vec4C T2vec = E*cross(eps,pBhat,pKhat)
                       - i*F*eps
                       + i*(eps*pBhat)*(G*phat+H*qhat);
    for( int hlm=0; hlm<2; hlm++ ) {
      for( int hlp=0; hlp<2; hlp++ ) {
	Complex M(0.0,0.0);
	M += Func.X( 2,hlm, T1vec, 3,hlp, Complex(1.0,0.0), Complex(1.0,0.0) );
	M += Func.X( 2,hlm, T2vec, 3,hlp, Complex(1.0,0.0), Complex(-1.0,0.0) );
        vector<pair<int,int> > spins;
        spins.push_back(make_pair(p_i[0],0));    // B
        spins.push_back(make_pair(p_i[1],hhad)); // K*
        spins.push_back(make_pair(p_i[2],hlm));  // lepton-
        spins.push_back(make_pair(p_i[3],hlp));  // lepton+
        amps->Add(m_global * M, spins);
      }
    }
  }
}
 
Complex B_KStar_Semileptonic::C9sehgal(double sHat) {
  double LUT = m_Vub/m_Vtb;
  return m_C9 + gc(sHat)*(3.0*m_C1+m_C2+3.0*m_C3+m_C4+3.0*m_C5+m_C6)
    -1.0/2.0*g0(sHat)*(m_C3+3.0*m_C4)-1.0/2.0*g(sHat)
    *(4.0*m_C3+4.0*m_C4+3.0*m_C5+m_C6)
    +2.0/9.0*(3.0*m_C3+m_C4+3.0*m_C5+m_C6) 
    -LUT*(3.0*m_C1+m_C2)*(g0(sHat)-gc(sHat));
}

Complex B_KStar_Semileptonic::sehgalld(double sHat) {
  Complex i = Complex(0.0,1.0);
  double LUT = m_Vub/m_Vtb;
  double alphaQED=1.0/129.0;
  double rho=9.0/sqr(alphaQED);
  double md = 0.135;
  double mdh = md/p_masses[0];
  double fump=0.875;
  double DELTA = 0.0;
  double DELTAI = 0.0;
  double mcc[6] = {3.097, 3.686, 3.770, 4.040, 4.159, 4.415};
  double mcch[6];
  double BR[6] = {6.02e-2, 8.1e-3, 1.12e-5, 1.4e-5, 1.0e-5, 1.1e-5};
  double GammaTot[6] = {0.087e-3, 0.277e-3, 23.6e-3, 52.0e-3, 78.0e-3, 43.0e-3};
  double GammaToth[6];
  for(int i=0; i<6; i++){
    mcch[i] = mcc[i]/p_masses[0]; 
    GammaToth[i] = GammaTot[i]/p_masses[0]; 
    double ZZ = M_PI*(sHat-sqr(mcch[i]))+2.0*(sqr(mcch[i])-sHat)*
      atan((4.0*sqr(mdh)-sqr(mcch[i]))/(GammaToth[i]*mcch[i]))-
      GammaToth[i]*mcch[i]*log(16.0*pow(mdh,4)+sqr(GammaToth[i])*sqr(mcch[i])-
      8.0*sqr(mdh)*sqr(mcch[i])+pow(mcch[i],4))+
      2.0*GammaToth[i]*mcch[i]*log(abs(4.0*sqr(mdh)-sHat));
    DELTA += BR[i]/mcch[i]*GammaToth[i]*ZZ/
      (sqr(sqr(mcch[i])-sHat)+sqr(mcch[i])*sqr(GammaToth[i]));
    DELTAI += sHat*BR[i]*sqr(GammaToth[i])/
      (sqr(sqr(mcch[i])-sHat)+sqr(mcch[i])*sqr(GammaToth[i]));
  }
  return fump*
    (sHat/3.0*(-rho/2.0)*DELTA+i*M_PI/3.0*rho*DELTAI)*(1.0+LUT);
}

Complex B_KStar_Semileptonic::g(double shat)
{
  Complex i = Complex(0.0,1.0);
  double mu = 4.8;
  double y = 4.0/shat;  
  return -8.0/9.0*log(p_masses[0]/mu)+8.0/27.0+4.0/9.0*y-2.0/9.0*(2.0+y)
    *sqrt(abs(1.0-y))*(Theta(1.0-y)*(log(abs(
    (1.0+sqrt(abs(1.0-y) ) )/(1.0-sqrt(abs(1.0-y)))))
    -i*M_PI )
    +Theta(y-1.0)*2.0*atan(1.0/sqrt(abs(y-1.0))));
}

Complex B_KStar_Semileptonic::g0(double shat)
{
  Complex i = Complex(0.0,1.0);
  double mu = 4.8;
  return 8.0/27.0-4.0/9.0*log(shat*p_masses[0]/mu)+4.0/9.0*i*M_PI;
}

Complex B_KStar_Semileptonic::gc(double shat)
{
  Complex i = Complex(0.0,1.0);
  double a = -6.8;
  double b = 11.33;
  double c = 1.02;
  double xl = 0.69;
  double md = 1.8693;
  double mdh = md/p_masses[0];
  double mc = 1.4;
  return -8.0/9.0*log(mc/p_masses[0])-4.0/9.0+shat/3.0*(
     (c-a)/shat*log(xl)+a/shat*log(4.0*sqr(mdh))+
     (a+b*shat-c)/shat*log(abs(xl-shat))-
     (a+b*shat)/shat*log(abs(4.0*sqr(mdh)-shat)))+i*
     M_PI/3.0*(Theta(xl-shat)*Theta(shat-0.6)*(a+b*shat)+
     c*Theta(shat-xl)*Theta(1.0-shat)); 
}

double B_KStar_Semileptonic::Theta(double x){
  if(x>0) return 1.0;
  else return 0.0;
}


Complex B_KStar_Semileptonic::gSD(double mhat, double shat)
{
  Complex i = Complex(0.0,1.0);
  double y = 4.0*mhat*mhat/shat;
  Complex theta;
  if (1 > y) theta = log((1.0+sqrt(1.-y))/(1.0-sqrt(1.0-y)))-i*M_PI;
  else       theta = 2.0*atan(1.0/sqrt(y-1.0));
  return -8./9.*log(mhat) +
          8./27. +
          4./9.*y -
          2./9.*(2.+y)*sqrt(fabs(1.-y))*theta;
}

DEFINE_ME_GETTER(B_KStar_Semileptonic,B_KStar_Semileptonic_Getter,"B_KStar_Semileptonic")

void B_KStar_Semileptonic_Getter::PrintInfo(std::ostream &st,const size_t width) const {
  st<<"Example: $ B \\rightarrow K^* \\; l^- \\; l^+ $ \n \n"
    <<"Order: 0 = $B$, 1 = $K^*$, 2 = $l^-$, 3 = $l^+$ \n\n"
    <<"For matrix element and form factors: hep-ph/9910221 \n"
    <<endl;
}
