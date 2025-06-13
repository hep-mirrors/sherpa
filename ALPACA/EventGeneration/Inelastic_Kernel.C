#include "ALPACA/EventGeneration/Inelastic_Kernel.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Math/Random.H"


#include <string>
#include <vector>
#include <cmath>
#include <complex>
#include<stdio.h>

using namespace ALPACA;
using namespace ATOOLS;
using namespace std;

Inelastic_Kernel::Inelastic_Kernel(shared_ptr<Dynamic_Quantities_Handler> ptr_dyn_quant_handler, shared_ptr<Analysis_Handler> ptr_analysis_handler):
  p_dyn_quant_handler(ptr_dyn_quant_handler), p_analysis_handler(ptr_analysis_handler),
  m_fixed_gamma(HIPars.FixedGamma()), m_gaussian_kT2(HIPars.GaussiankT2()), 
  m_kT2_reg(HIPars.kT2Reg()), m_formation_time(HIPars.FormationTime())
{};


Inelastic_Kernel::~Inelastic_Kernel() {
};

double Inelastic_Kernel::FindGamma(double p, const ATOOLS::Flavour flavour_a, const ATOOLS::Flavour flavour_b, const ATOOLS::Flavour flavour_c, 
                                  double x, double Tstar, double mg2, double mf2){
    /*
      Find the effective splitting kernel \gamma^a_{bc}(p,xp,[1-x]p) as defined in AMY
    */

  double alpha_s =  p_dyn_quant_handler->GetAlphaS();
  double gamma, prefactor, mu2_old, mu2_new;
  double CA = 3, dA = 8.;
  double CF =  4./3., dF = 3.;
  double N_c = 3;
  double lambda = 4.*M_PI*alpha_s*N_c;
  double eta = x*(1.-x)*lambda*Tstar*p/mg2;
  double xi = 9.09916;
  double mu2_interp, mu2_BH, mu2_LPM;
  double Mhat2;
  double m2s1, m2s2, m2s3;
  const complex<double> i(0.,1.);

  if(flavour_a.IsGluon()){
    if(flavour_b.IsGluon() || flavour_c.IsGluon()){ // g <-> gg
      prefactor = sqrt(2.)*dA*CA*alpha_s/(pow(2.*M_PI,4));
      prefactor *= (1.+pow(x,4)+pow(1.-x,4))/(x*x*pow(1.-x,2));
      m2s1 = mg2;
      m2s2 = mg2;
      m2s3 = mg2;
      Mhat2 = 1. - x + x*x; //TEMPORARY
    } else{ //g <-> qqb
      prefactor = sqrt(2.)*dF*CF*alpha_s/(pow(2.*M_PI,4));
      prefactor *= (x*x+pow(1.-x,2))/(x*(1-x));
      m2s1 = mg2;
      m2s2 = mf2;
      m2s3 = mf2;
      Mhat2 = 4./9.-x*(1.-x); //TEMPORARY
    }
  } else if(flavour_a.IsQuark()){ // q <-> qg
    prefactor = sqrt(2.)*dF*CF*alpha_s/(pow(2.*M_PI,4));
    prefactor *= (1.+pow(1.-x,2))/(x*x*(1.-x));
    m2s1 = mf2;
    m2s2 = mf2;
    m2s3 = mg2;
    Mhat2 = 1.-x+(4./9.)*x*x;
  } else THROW(fatal_error,"Did not find process");

  double C_s1s2s3 = FindCs1s2s3(flavour_a, flavour_b, flavour_c);
  double C_s2s3s1 = FindCs1s2s3(flavour_b, flavour_c, flavour_a);
  double C_s3s1s2 = FindCs1s2s3(flavour_c, flavour_a, flavour_b);

  // #### Solve for mu^2_BH ####
  //Mhat2 = x*(1-x)*( m2s1/(x*(1.-x)) - m2s2/x - m2s3/(1.-x) )/mg2;
  //msg_Out() << "Mhat2 = " << Mhat2 << endl;
  //msg_Out() << "m2s1/(x*(1.-x)) = " << m2s1/(x*(1.-x)) << ", m2s2/x = " << m2s2/x << ", m2s3/(1.-x) = " << m2s3/(1.-x) << endl;
  if(Mhat2 < pow(10,-10)){
    mu2_BH = 0;
  } else{
    mu2_BH = (eta*4.*M_PI/sqrt(2))*( C_s1s2s3*QBH(2./Mhat2) + C_s2s3s1*QBH(x*x*2./Mhat2) + C_s3s1s2*QBH((1.-x)*(1.-x)*2./Mhat2) );
    mu2_BH *= mg2;
  }


  // #### Solve for mu^2_LPM ####
  // Iterative solution to solve for mu^2_LPM, start with mu^2 = sqrt(2eta)
  int N_counter = 0;
  int N_limit = 0; //0 to match Mathematica data, anything larger produces too small results
  mu2_new = sqrt(2.*eta);
  while(N_counter <= N_limit){
    mu2_old = mu2_new;
    //mu2_new = sqrt((eta/(2.*M_PI))*( C_s1s2s3*log(xi*mu2_old/mg2+1) + C_s2s3s1*x*x*log(xi*mu2_old/(x*x*mg2)+1) + C_s3s1s2*(1-x)*(1-x)*log(xi*mu2_old/((1-x)*(1-x)*mg2)+1) ));
    //Initial guess for mu2_new should actually be for mu2/mg2 according to mathematica doc
    mu2_new = sqrt((eta/(2.*M_PI))*( C_s1s2s3*log(xi*mu2_old+1.) + C_s2s3s1*x*x*log(xi*mu2_old/(x*x)+1.) + C_s3s1s2*(1.-x)*(1.-x)*log(xi*mu2_old/((1.-x)*(1.-x))+1.) ));
    N_counter += 1;
  }

  mu2_LPM = mg2*mu2_new;

  mu2_interp = (sqrt(eta+1.)-1.)*(2.*mu2_BH/(eta*(1.+eta)) + eta*mu2_LPM/(sqrt(eta)*(1.+eta)));

  if(std::isnan(mu2_interp) || std::isinf(mu2_interp)){
    mu2_interp = 0.;
  }

  //msg_Out() << "C_s1s2s3 = " << C_s1s2s3 << ", C_s2s3s1 = " << C_s2s3s1 << ", C_s3s1s2 = " << C_s3s1s2 << endl;
  //msg_Out() << "Mhat2 = " << Mhat2 << endl;
  //msg_Out() << "mu2_BH = " << mu2_BH << ", mu2_LPM = " << mu2_LPM << ", mu2_interp = " << mu2_interp << endl;


  gamma  = prefactor * mu2_interp;

  return gamma;
}


double Inelastic_Kernel::FindCs1s2s3(ATOOLS::Flavour s1, ATOOLS::Flavour s2, ATOOLS::Flavour s3){
  /*
    Find factors C(s1, s2, s3) = (C_s2+C_s3-C_s1)/C_A used in \mu^2_BH and \mu2_LPM in \gamma (collinear splitting)
  */
  double C_s1, C_s2, C_s3;
  if(s1.IsQuark()){
    C_s1 = 4./3.;
  } else{
    C_s1 = 3;
  }
  if(s2.IsQuark()){
    C_s2 = 4./3.;
  } else{
    C_s2 = 3;
  }
  if(s3.IsQuark()){
    C_s3 = 4./3.;
  } else{
    C_s3 = 3;
  }
  return (C_s2+C_s3-C_s1)/3.;
}

double Inelastic_Kernel::QBH(double r){
  /*
    Find Q_BH(r) used in \mu^2_BH in \gamma (collinear splitting)
  */

  const complex<double> i(0.,1.);
  complex<double> Q(0.,0.);

  complex<double> r_plus;
  complex<double> r_minus;

  if(r == 0){
    Q = 0.;
  } else if(r == 2. || r == 4.){
    Q = (2. - log(r))/(8.*M_PI*M_PI);
  } else{
    double tol = 0.01;
    int N_max = 100;

    if(r < 4.){
      r_plus = 1. - r/2. + i*sqrt((4.-r)*r)/2.;
      r_minus = 1. - r/2. - i*sqrt((4.-r)*r)/2.;
      Q = ( i*(r-2.)*(Li2(r_minus,tol,N_max)-Li2(r_plus,tol,N_max))/(sqrt((4.-r)*r)) + 2. - log(r) )/(8.*M_PI*M_PI);
    } else{
      r_plus = 1. -r/2. - sqrt((r-4.)*r)/2.;
      r_minus = 1. -r/2. + sqrt((r-4.)*r)/2.;
      Q = ( (r-2.)*(Li2(r_minus,tol,N_max)-Li2(r_plus,tol,N_max))/(sqrt((r-4.)*r)) + 2. - log(r) )/(8.*M_PI*M_PI);
    }

    
  }

  //msg_Out() << "Q_BH = " << Q << " for r = " << r << endl; 
  if(!Q.imag()){
    return Q.real();
  } else{
    msg_Out() << "#### ERROR: QBH.imag() != 0" << " for r = " << r << endl;
    msg_Out() << "r_plus = " << r_plus << ", r_minus = " << r_minus << endl;
    msg_Out() << "  QBH = " << Q << ", returning 0." << endl;
    msg_Out() << "Li2(r_minus,tol,N_max) = " << Li2(r_minus,0.01,100) << ", Li2(r_plus,tol,N_max) = " << Li2(r_plus,0.01,100) << endl;
    return Q.real();
    //return 0.;
  }

}


std::complex<double> Inelastic_Kernel::Li2(std::complex<double> z, double tol, int N_max){
  /*
    Function to compute the dilograithm.
    Input:
      z     - Complex number for Li_2(z)
      tol   - Tolerance for when to end numerical calculation, when r_new/r_old and theta_old/theta_new are within 1+-tol
      N_max - Cutoff, if not within 1+-tol after N_max iterations, break
  */
  complex<double> dilog_new(0.,0.);
  complex<double> dilog_old(0.,0.);
  complex<double> z_sum(0.,0.);
  double sign = 1.;

  complex<double> real_zero(0.,0.);
  complex<double> real_one(1.,0.);
  const complex<double> i(0.,1.);

  bool condition_r = true, condition_theta = true; //For sum expansion, make sure both amplitude and complex angle ratios converges within 1+-tol

  double r = abs(z);
  if(z == real_zero){         //If z = 0+0i, set dilog = 0
    dilog_new = real_zero;
  } else if(z == real_one){  //If z = 1+0i, set dilog = pi^2/6
    dilog_new = M_PI*M_PI/6;
  } else{                     //If non of the special cases above, calculate from sum expansion
    //Select case for sum
    if(r < 1){                        //For |z| < 1
      z_sum = z;
    } else if(!z.imag() && z.real() > 1){    //For z real and z > 1
      if(r <= 2){                         //For z real and 1 < z <= 2
        // Should be this according to Morris, but does not match wolfram alpha for e.g. Li_2(1.5)                
        //z_sum = real_one - 1./z;
        //dilog_new = M_PI*M_PI/6. - log(z)*log(real_one-z) - log(z)/2. - i*M_PI*log(z);
        //
        z_sum = 1./z;
        sign = -1.;
        dilog_new = M_PI*M_PI/3. - pow(log(z),2)/2. - i*M_PI*log(z);
      } else {                            //For z real and 2 < z
        z_sum = 1./z;
        sign = -1.;
        dilog_new = M_PI*M_PI/3. - pow(log(z),2)/2. - i*M_PI*log(z);
      }
    } else{                           //For z not along the branch cut [1,infinity), and |z|>1
      z_sum = 1./z;
      sign= -1.;
      dilog_new = -M_PI*M_PI/6.-pow(log(-z),2)/2.;
    }
    //Calculate sum part
    int k = 1;
    while(condition_r || condition_theta){
      if(k-1 >= N_max){ //If series has not converged within 1+-tol after N_max iterations, break
        msg_Out() << "#### Warning ####" << endl;
        msg_Out() << "  Li_2(z) did not converge after " << k-1 << " iterations, ends iteration." << endl;
        msg_Out() << "  z = " << z << ", dilog_old = " << dilog_old << ", dilog_new = " << dilog_new << endl;
        msg_Out() << "  abs(dilog_new)/abs(dilog_old) = " << abs(dilog_new)/abs(dilog_old) << ", arg(dilog_new)/arg(dilog_old) = " << arg(dilog_new)/arg(dilog_old) << endl;
        break;
      }
      dilog_old = dilog_new;
      dilog_new += sign*pow(z_sum,k)/double(k*k);
      k += 1;
      condition_r = !(1.-tol <= abs(dilog_new)/abs(dilog_old) && abs(dilog_new)/abs(dilog_old) <= 1.+tol);
      condition_theta = !(1.-tol <= arg(dilog_new)/arg(dilog_old) && arg(dilog_new)/arg(dilog_old) <= 1.+tol);
    }
  }

  //msg_Out() << "dilog = " << dilog_new << " for z = " << z << endl;
  return dilog_new;
}


double Inelastic_Kernel::FindFormationTime(double E_in, double T_star, double mg2, double kT2){
  double t_form = 0.; 
  if(m_formation_time > 0.){
    double tau_soft = 1./(pow(2.*M_PI*(p_dyn_quant_handler->GetAlphaS()),2.)*T_star);
    double N_form = sqrt(1.+E_in/(mg2*tau_soft));
    if(kT2 < 0.){
      t_form = E_in/(N_form*mg2);
    } else{
      t_form = E_in/(N_form*kT2);
    }
    t_form = m_formation_time*t_form;
  }
  return t_form;
}


double Inelastic_Kernel::SamplekT2(double kT2_min, double kT2_max){
  double kT2;
  if(m_gaussian_kT2){
    //Accept/reject with constant overestimate
    double weight;
    int N = 0;
    bool found_kT2 = false;
    while(!found_kT2){
      kT2 = kT2_min + (kT2_max - kT2_min)*ran->Get();
      weight = exp(-pow(kT2/m_kT2_reg, 2.)/2.); //Sampling with constant over estimate of the Gaussian

      if(ran->Get() < weight){
        found_kT2 = true;
      }

      N += 1;
      if(N > pow(10., 7)){
        msg_Out() << METHOD << ": ERROR. Cannot find sample of k_T^2, exit" << endl;
        exit(1.);
      }
    }
  } else{
    //Primitive function known, sample directly
    kT2 = pow((kT2_max+m_kT2_reg)/(kT2_min+m_kT2_reg),ran->Get())*(kT2_min+m_kT2_reg) - m_kT2_reg;
  }

  return kT2;
}


double Inelastic_Kernel::FindMergingXsec(std::shared_ptr<Parton> part_in_1, std::shared_ptr<Parton> part_in_2, double T_star, double mg2, double mf2, double x, double kT2, double Ep, double Ek, double p_out_boost, double kT2_max){
  /*
    Function to find the cross section for 2->1 merging for part_in_1 and part_in_2. Returns xsec_merge.
    NOTE: does not include Bose/Pauli-factors
  */
  Flavour flav_out;
  double shat = (part_in_1->Momentum()+part_in_2->Momentum()).Abs2();
  double xsec_merge = 0.;
  double gamma = 0.;
  
  double p_in_1 = part_in_1->Abs3Momentum();
  double p_in_2 = part_in_2->Abs3Momentum();
  Vec4D p_temp_1 = part_in_1->Momentum();
  Vec4D p_temp_2 = part_in_2->Momentum();
  double p_out = p_in_1 + p_in_2;
  Flavour flav_in_1 = part_in_1->Flav();
  Flavour flav_in_2 = part_in_2->Flav();

  double nu_a, nu_b;
  if(flav_in_1.IsGluon() && flav_in_2.IsGluon()){ // gg -> g
    flav_out = kf_gluon;
    nu_a = 16.;
    nu_b = 16.;
    if(m_fixed_gamma.first){
      gamma = m_fixed_gamma.second;
    } else{
      gamma = FindGamma(p_out_boost, flav_out, flav_in_1, flav_in_2, x, T_star, mg2, mf2);
    }
  } else if( (flav_in_1.IsFermion() && flav_in_2.IsGluon()) || (flav_in_1.IsGluon() && flav_in_2.IsFermion()) ){ // gq -> q
    nu_a = 16.;
    nu_b = 6.;
    Flavour temp_flav_in;
    if(flav_in_1.IsFermion()){
      flav_out = flav_in_1;
      temp_flav_in = flav_in_1; //Make sure particle b (second position in gamma, which gets energy fraction x, is a gluon)
    } else{
      flav_out = flav_in_2;
      temp_flav_in = flav_in_2; //Make sure particle b (second position in gamma, which gets energy fraction x, is a gluon)
    }
    if(m_fixed_gamma.first){
      gamma = m_fixed_gamma.second;
    } else{
      gamma = 2.*FindGamma(p_out_boost, flav_out, flav_in_1, flav_in_2, x, T_star, mg2, mf2); //Factor *2 for q->gq and q->qg
    }
  } else if( (flav_in_1.IsAnti() && !flav_in_2.IsAnti()) || (!flav_in_1.IsAnti() && flav_in_2.IsAnti()) ){ // qqb -> g
    flav_out = kf_gluon;
    nu_a = 6.;
    nu_b = 6.;
    if(m_fixed_gamma.first){
      gamma = m_fixed_gamma.second;
    } else{
      gamma = FindGamma(p_out_boost, flav_out, flav_in_1, flav_in_2, x, T_star, mg2, mf2);
    }
  }

  if(gamma > 0){ //If not rejected, calculate cross section for merging
    double dirac_factor = Ek/Ep + Ep/Ek + 2.;
    double kT2_dist_factor;

    if(m_gaussian_kT2){
      kT2_dist_factor = 1./( (2./(sqrt(2*M_PI)*m_kT2_reg))*exp(-pow(kT2/m_kT2_reg, 2.)/2.) );
    } else{
      kT2_dist_factor = (kT2+m_kT2_reg)*log((kT2_max+m_kT2_reg)/m_kT2_reg);
    }

    xsec_merge = 2.*pow(2.*M_PI,5)*gamma/(nu_a*nu_b*shat*dirac_factor*kT2_dist_factor);
  }

  return xsec_merge;
}