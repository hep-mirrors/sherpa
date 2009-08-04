#include "HADRONS++/ME_Library/B_K_Semileptonic.H"
#include "ATOOLS/Org/Message.H"
#include "METOOLS/Main/XYZFuncs.H"
#include "ATOOLS/Math/Histogram.H"
#include "MODEL/Main/Model_Base.H"

using namespace HADRONS;
using namespace ATOOLS;
using namespace METOOLS;
using namespace MODEL;
using namespace std;

void B_K_Semileptonic::SetModelParameters( GeneralModel _md )
{ 
  double GF = _md("GF",s_model->ScalarConstant(string("GF")));
  double alpha = _md("alpha",1./129.);
  double Vts = _md("Vts",Tools::Vts);
  m_Vtb = _md("Vtb",Tools::Vtb);
  m_Vub = _md("Vub",Tools::Vub);
  m_global = GF*alpha*Vts*m_Vtb/(2.*M_PI)*SQRT_05;
  m_LD    = bool(int(_md("LD",1.0)+0.5));
  m_cR_T1 = Complex(1.0 ,0.0); 
  m_cR_T2 = Complex(1.0 ,0.0);
  m_cL_T1 = Complex(1.0 ,0.0);
  m_cL_T2 = Complex(-1.0,0.0);

  m_C1  = _md("C1",-0.248);
  m_C2  = _md("C2",1.107);
  m_C3  = _md("C3",0.011);
  m_C4  = _md("C4",-0.026);
  m_C5  = _md("C5",0.007);
  m_C6  = _md("C6",-0.031);
  m_C7eff = _md("C7eff",-0.313);
  m_C9  = _md("C9",4.344);
  m_C10 = _md("C10",-4.669);

  m_mc  = _md("mc",1.4);
  m_ms  = _md("ms",0.2); 
}


void B_K_Semileptonic::operator()( 
    const Vec4D             * _p,
    METOOLS::Spin_Amplitudes * amps)
{
  double s     = (_p[p_i[2]]+_p[p_i[3]]).Abs2();
  double shat  = s/sqr(p_masses[0]);
  
  double m_bpole = 4.8;
  double mu = m_bpole;
  double alphas = s_model->ScalarFunction("alpha_S",mu);
  m_mb = m_bpole*(1.-4./3.*alphas/M_PI);

  // calculate formfactors
  double fplus  = 0.319 * exp(1.465*shat+0.372*sqr(shat)+0.782*pow(shat,3));
  double fzero  = 0.319 * exp(0.633*shat-0.095*sqr(shat)+0.591*pow(shat,3));
  double fminus = (fzero -fplus)/(shat*sqr(p_masses[0]))*(sqr(p_masses[0])-sqr(p_masses[1]));
  double fT     = 0.355 * exp(1.478*shat+0.373*sqr(shat)+0.7*pow(shat,3));
  
  // calculate coefficients
  Complex C9eff;
  if(m_LD==true) C9eff = C9sehgal(shat) + sehgalld(shat);
  else C9eff = m_C9 + gSD(m_mc/p_masses[0],shat) * (3.*m_C1+m_C2+3.*m_C3+m_C4+3.*m_C5+m_C6)
		- 0.5 * gSD(m_mb/p_masses[0],shat) * (4.*m_C3+4.*m_C4+3.*m_C5+m_C6)
		- 0.5 * gSD(m_ms/p_masses[0],shat) * (m_C3+3.*m_C4)
		+ (2./9.) * (3.*m_C3+m_C4+3.*m_C5+m_C6);
  Complex Aprime =  C9eff * fplus +
                    2.*m_mb/p_masses[0]/(1+p_masses[1]/p_masses[0]) * m_C7eff * fT;
  
  Complex Bprime =  C9eff * fminus -
		    2.*m_mb/p_masses[0]/shat*(1.-p_masses[1]/p_masses[0]) * m_C7eff * fT;
  Complex Cprime =  m_C10 * fplus;
  Complex Dprime =  m_C10 * fminus;
  
  XYZFunc F(m_n,_p,m_flavs,Tools::k0,m_anti,p_i);
  for( int hlm=0; hlm<2; hlm++ ) {
    for( int hlp=0; hlp<2; hlp++ ) {
      Complex amplitude = Aprime*F.X(2,hlm, _p[p_i[0]]+_p[p_i[1]], 3,hlp, m_cR_T1, m_cL_T1)
			+ Bprime*F.X(2,hlm, _p[p_i[0]]-_p[p_i[1]], 3,hlp, m_cR_T1, m_cL_T1)
			+ Cprime*F.X(2,hlm, _p[p_i[0]]+_p[p_i[1]], 3,hlp, m_cR_T2, m_cL_T2)
			+ Dprime*F.X(2,hlm, _p[p_i[0]]-_p[p_i[1]], 3,hlp, m_cR_T2, m_cL_T2);
      vector<pair<int,int> > spins;
      spins.push_back(make_pair(p_i[0],0));    // B
      spins.push_back(make_pair(p_i[1],0));    // K
      spins.push_back(make_pair(p_i[2],hlm));  // lepton-
      spins.push_back(make_pair(p_i[3],hlp));  // lepton+
      amps->Add(amplitude*m_global*p_masses[0], spins);
    }
  }
}
 
Complex B_K_Semileptonic::C9sehgal(double sHat) {
  double LUT = m_Vub/m_Vtb;
  return m_C9 + gc(sHat)*(3.0*m_C1+m_C2+3.0*m_C3+m_C4+3.0*m_C5+m_C6)
    -1.0/2.0*g0(sHat)*(m_C3+3.0*m_C4)-1.0/2.0*g(sHat)
    *(4.0*m_C3+4.0*m_C4+3.0*m_C5+m_C6)
    +2.0/9.0*(3.0*m_C3+m_C4+3.0*m_C5+m_C6) 
    -LUT*(3.0*m_C1+m_C2)*(g0(sHat)-gc(sHat));
}

Complex B_K_Semileptonic::sehgalld(double sHat) {
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

Complex B_K_Semileptonic::g(double shat)
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

Complex B_K_Semileptonic::g0(double shat)
{
  Complex i = Complex(0.0,1.0);
  double mu = 4.8;
  return 8.0/27.0-4.0/9.0*log(shat*p_masses[0]/mu)+4.0/9.0*i*M_PI;
}

Complex B_K_Semileptonic::gc(double shat)
{
  Complex i = Complex(0.0,1.0);
  double a = -6.8;
  double b = 11.33;
  double c = 1.02;
  double xl = 0.69;
  double md = 1.8693;
  double mdh = md/p_masses[0];
  return -8.0/9.0*log(m_mc/p_masses[0])-4.0/9.0+shat/3.0*(
     (c-a)/shat*log(xl)+a/shat*log(4.0*sqr(mdh))+
     (a+b*shat-c)/shat*log(abs(xl-shat))-
     (a+b*shat)/shat*log(abs(4.0*sqr(mdh)-shat)))+i*
     M_PI/3.0*(Theta(xl-shat)*Theta(shat-0.6)*(a+b*shat)+
     c*Theta(shat-xl)*Theta(1.0-shat)); 
}

double B_K_Semileptonic::Theta(double x){
  if(x>0) return 1.0;
  else return 0.0;
}

Complex B_K_Semileptonic::gSD(double mhat, double shat)
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

DEFINE_ME_GETTER(B_K_Semileptonic,B_K_Semileptonic_Getter,"B_K_Semileptonic")

void B_K_Semileptonic_Getter::PrintInfo(std::ostream &st,const size_t width) const {
  st<<"Example: $ B \\rightarrow K \\; l^- \\; l^+ $ \n \n"
    <<"Order: 0 = $B$, 1 = $K$, 2 = $l^-$, 3 = $l^+$ \n\n"
    <<"For matrix element and form factors: hep-ph/9910221 \n"
    <<endl;
}
