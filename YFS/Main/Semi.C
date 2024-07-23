#include "YFS/Main/Semi.H"

#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/My_Limits.H"
#include "ATOOLS/Phys/Blob.H"
#include "MODEL/Main/Model_Base.H"
#include "MODEL/Main/Running_AlphaQED.H"

#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
//#ifdef PHOTONS_DEBUG
//#include "ATOOLS/Math/Histogram_2D.H"
//#include "ATOOLS/Org/Shell_Tools.H"
//#endif

using namespace YFS;
using namespace ATOOLS;
using namespace MODEL;
std::ofstream outfile;

#define PARAMETER_TYPE PHASIC::Process_Info

// member functions of class Photon
Semi_YFS::Semi_YFS() :
m_name("YFS test")
{  
    // Data_Reader     reader(" ,",";","#","=");       
    // m_BeamEnergy1   = reader.GetValue<double>("BEAM_ENERGY_1",1);
    // m_BeamEnergy2   = reader.GetValue<double>("BEAM_ENERGY_2",1);
    // m_yfsmode       = reader.GetValue<int>("YFS_ISR_MODE",0);
    // m_s             = pow(rpa->gen.Ecms(),2);
    // mass_f          = 0.000511;
    // // mass_f                       = rpa->gen.Beam1().Mass();
    // m_alqed         = 1./137.;
    // m_epsilon       = 1e-4;         
    // // m_epsilon            = reader.GetValue<double>("YFS_IR_CUTOFF",0.001);               
    // m_Kmin           = reader.GetValue<double>("YFS_IR_CUTOFF",0.01);
    // m_Kmax           = (m_BeamEnergy1+m_BeamEnergy2)/2*0.1; //ask Frank about this value
    // m_var   = ran->Get();
    // m_Emass = 1.;
    // Process_Info pi;
    // Flavour_Vector fl;
    // msg_Out()<<"Process Info = "<<pi<<"\n";
    // msg_Out()<<"Flavour = "<<fl<<"\n";
}

void Semi_YFS::initial() 
{
}


double Semi_YFS::bornxs(double s,double theta){
    double Qf     = -1;
    double sin2th =  0.22293;
    double vmu    = -0.5+2*sin2th;
    double t3mu   = -0.5;
    double Amu    = t3mu;
    double Mz     = 91.18;
    double alpha  = 0.00776384;
    double pi     = M_PI;
    double gammaz = 2.49;
    double GF     = 1.116e-5;
    double mz2    = pow(Mz,2.);
    double kappa  = pow(2.,0.5)*GF*mz2/(16.*pi*alpha); // c.f eq 3.2 pg 55 
    double chi1   = kappa*(s*(s-mz2)/(pow((s-mz2),2.)+gammaz*gammaz*mz2));
    double chi2   = pow(kappa,2.)*(pow(s,2.)/(pow((s-mz2),2.)+pow(gammaz,2.)*mz2));
    
    double t1     = pow(Qf,2.)-2.*Qf*vmu*vmu*chi1;
    double t2     = pow((pow(Amu,2.)+pow(vmu,2.)),2.)*chi2;
    double t3     = cos(theta)*(-4.*Qf*Amu*Amu*chi1+8.*Amu*Amu*vmu*vmu*chi2);
    double dsigma = pi*alpha*alpha/(2.*s)*pow((1.+cos(theta)),2)*(t1+t2+t3);
    return dsigma;
}

double Semi_YFS::sigma_crude(double s,double v)
{
    double mass_f = 0.000511;
    double alpha  = 1./137;
    double L      = log((s*s)/(mass_f*mass_f));
    double pi     = M_PI;
    double gamma  = 2.*alpha/pi*(L-1.);
    double gammap = 2.*alpha/pi*(L);
    double yfsfac = exp(alpha/pi*(0.5*L-1.+(pi*pi)/3.));
    double Jac    = 0.5*(1.+pow((1.-v),-0.5));
    double E      = 1.+gammap/gamma*pow(v/0.01,(gammap-gamma));
    double t1     = yfsfac*gamma;
    double t2     = pow(v,(gamma-1.))*Jac*E;
    double t3     = bornxs(s*s*(1.-v),pi/2.);
    double tot    = t1*t2*t3;
// #     tot    = dxs(s**2*(1.-v),np.pi/2)
    return tot*3.89379e8;
}


double Semi_YFS::Weight( double v, double alp, const int &semiyfs)
{   
    // if (x<=(1.-0.999999)) return 0;
    // if(x<1e-06 || x>1.-1e-6) return 0;
    double m_rho;
    double m_alpha = alp;
    // alp = 1.0/137;
    double x = 1.-v;
    double m_s = pow(rpa->gen.Ecms(),2.0);
    double m_sp = m_s*(1.-x);
    double mass_f  = rpa->gen.Beam1().Mass();
    double mass_f2 = mass_f*mass_f;
    double L       = log(m_s/mass_f2);
    double beta_l = 2*alp/M_PI*(L-1.);
    double gamma_f = exp(Gammln(1. + beta_l));
    double C = 0.5772156649;
    double rho0 = exp(beta_l*(1.0/4. - GAMMA_E) + alp/M_PI*(pow(M_PI,2.)/3.-0.5))/gamma_f*beta_l*pow(1-x,(beta_l-1.)); //see arxiv 1607.03210 Appendix
    // rho0 = rho0*(1.0 + beta_l*0.5 -0.5*(1.0-x*x));
    double diljac0 = 0.5*(1.+1./sqrt(1.-x));
    if(semiyfs == 1){
        m_rho = rho0;//*(1.-0.25*beta_l*log(1.-x)-0/5*m_alpha*pow(log(1.-x),2.));
    }
  
    else if (semiyfs == 2){
        m_rho = rho0*(1.0  + beta_l*0.5 +x*(-1.0+x/2.)+beta_l*(-x*x/2.+0.25*(-1+4.*x-2.*x*x)*log(1.-x)));
    }
    else if (semiyfs == 3){
        m_rho = rho0*(1.0 + beta_l*0.5 + 0.125*beta_l*beta_l-0.5*(1.0-x*x)+ beta_l*(-0.5*(1.-x)-0.125*(1.0+3.*x*x)*log(x)));
    }
    else if (semiyfs == 5){
        int i = 1;
        double Li = 0.;
        double del= 1.;
        double del2=1.;
        for (i=1;del2>rpa->gen.Accu();++i) {  
        del *=x;
        del2 = del/double(i*i);
        Li  += del2;
        }
        
        double t1 = 0.5*(1.0+x*x) + m_alpha/M_PI*beta_l*(3./32. - M_PI*M_PI + 3./2.*(1.202057)); // zeta(3) ~ 1.202057
        double t2 = 0.25*beta_l*(-0.5*(1.+3.*x*x)*log(x) - pow(1.-x,2));
        
        double st1 = -(1.+3.*x*x)*log(x)*log(x);
        double st2 = 4.*(1.+x*x)*(Li*(1.-x)+log(x)*log(1-x));
        double st3 = 2.*(1.-x)*(3.-2.*x) + 2.*(3. + 2.*x + x*x)*log(x);

        double t3 = 0.125*m_alpha/M_PI*(st1 + st2 + st3);
        double t4 = 0.125*beta_l*beta_l*(0.5*(3.0*x*x - 4*x + 1)*log(x)+1./12.*(1+7*x*x)*log(x)*log(x) + (1.-x*x)*Li*(1.-x)
                    + pow(1.-x,2));
        
        m_rho = rho0*(t1 + t2 + t3 + t4);
    }
    else if (semiyfs == 6){
        m_rho = m_alpha/M_PI*pow(1-x,(beta_l-1.))*((1+L)*log(4*0.0001)+L*(1+0.5*L)+M_PI*M_PI/3.);
    }
    // if((1.-x)<1e-04) m_rho*= pow(100.,beta_l/2)/(pow(100.,beta_l/2)-1.);
    return m_rho;
}
