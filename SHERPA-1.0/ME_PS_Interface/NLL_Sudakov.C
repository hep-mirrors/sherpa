#include "NLL_Sudakov.H"
#include "Message.H"
#include "MathTools.H"
#include <stdio.h>

using namespace MOCAIC;
using namespace AMEGIC;
using namespace APHYTOOLS;
using namespace AORGTOOLS;
using namespace AMATOOLS;

NLL_Sudakov::NLL_Sudakov(double _tmax,double _tmin,bool _run) :
  tmax(_tmax), tmin(_tmin), run_as(_run), gauss(this), inte(NLL::analytic) {

  qmin   = sqrt(tmin);
  qmax   = sqrt(tmax);

  Nc     = 3;
  CA     = Nc;
  CF     = (Nc*Nc-1)/(2.*Nc);
  Nf     = 5;
  TR     =  1./2.;
 
  beta0  = (11.*CA-2.*Nf)/3.;
  beta1  = (17.*CA*CA- 3.*CF*Nf-5.*CA*Nf)/3.;
  FixLambda2();
  nlo    = 0;
  K      = 0;             // LO
  // NLO:  
  nlo    = 0;
  // with     
  if (nlo) K  =  CA*(67./18.-M_PI*M_PI/6.)-10./9.*TR*Nf;

  msg.Debugging()<<"Init the NLL_Sudakov :"<<std::endl
		 <<"  Nf = "<<Nf<<", TR = "<<TR<<", CA = "<<CA<<", CF = "<<CF<<std::endl
		 <<"  lambda = "<<sqrt(lambda2)<<"   --->  "<<std::endl
		 <<"  alphaS test: as(mue  ="<<sqrt(mu2)<<") = "<<asmu<<std::endl
		 <<"               as(qmin = "<<qmin<<") = "<<AlphaS(qmin*qmin)<<std::endl
		 <<"               as(qmax = "<<qmax<<") = "<<AlphaS(qmax*qmax)<<std::endl;
};


void NLL_Sudakov::FixLambda2() 
{
  mu2     = sqr(Flavour(kf::Z).Mass());
  asmu    = (*as)(mu2);
  lambda2 = mu2 * exp(1./(beta0 * asmu));
}                 

double NLL_Sudakov::AlphaS(double t)                
{ 
  if (t<0.) t = -t;
  // - lambda - parametrisation:
  //  return 4.*M_PI/(beta0*log(t/lambda2));
  // - first order:
  double   w = 1.-beta0*asmu/(4.*M_PI)*log(mu2/t);
  double   a = asmu/w;
  if (nlo) a *= (1. - (beta1*asmu)/(beta0*2.*M_PI*w)*log(w));
  return a;
}





// nomenclature Frank-Preprint
double NLL_Sudakov::GammaQ(double q, double Q) {
  double a = AlphaS(q*q);
  return 2.*CF/M_PI * a/q * ((1.+a/(2.*M_PI)*K)*log(Q/q)-3./4.); 
}

double NLL_Sudakov::GammaG(double q, double Q) {
  double a = AlphaS(q*q);
  return 2.*CA/M_PI * a/q * ((1.+a/(2.*M_PI)*K)*log(Q/q)-11./12.);
}

double NLL_Sudakov::GammaF(double q) {
  return Nf/(3.*M_PI) * AlphaS(q*q)/q;
}

double NLL_Sudakov::DeltaQ(double Q,double q) 
{ 
  if ((Q-q)/q<1.e-8) return 1.;
  /*
    if (inte==NLL::use_table) {
    if (q != qmin) {
    msg.Error()<<" ERROR: wrong q!!!!"<<std::endl;
    // WARNING:  Fast_Func has to be extended to deal with two variables!!!!
    //           in case qmin changes!!! This is indeed the case, see e+e- with
    //           ISR. Therefore I commented this out.
    msg.Error()<<" q=    " <<q<<std::endl;
    msg.Error()<<" qmin= " <<qmin<<std::endl;
    } 
    return exp(-log_delta_q(Q));
    }
  */
  double dq = exp(-IntGammaQ(q,Q));
  if (dq<=1.) return dq;
  return 1.;
}

double NLL_Sudakov::DeltaG(double Q,double q) 
{ 
  if ((Q-q)/q<1.e-8) return 1.;
  /*
    if (inte==NLL::use_table) {
    if (q != qmin) {
    msg.Error()<<" ERROR: wrong q!!!!"<<std::endl;
    // WARNING:  Fast_Func has to be extended to deal with two variables!!!!
    //           in case qmin changes!!! This is indeed the case, see e+e- with
    //           ISR. Therefore I commented this out.
    msg.Error()<<" q=    " <<q<<std::endl;
    msg.Error()<<" qmin= " <<qmin<<std::endl;
    } 
    return exp(-log_delta_g(Q));
    }
  */
  double dg= exp(-IntGammaG(q,Q)-IntGammaF(q,Q));
  if (dg<=1.) return  dg;
  return 1.;
}

double NLL_Sudakov::DeltaF(double Q,double q) 
{ 
  return sqr(DeltaQ(Q,q))/DeltaG(Q,q);
}



double NLL_Sudakov::IntGammaQ(double Q0, double Q) {
  switch (inte) {
  case NLL::use_table : 
    if (Q0!=qmin) {
      msg.Error()<<" ERROR: wrong q!!!!"<<std::endl;
      // WARNING:  Fast_Func has to be extended to deal with two variables!!!!
      //           in case qmin changes!!!
      msg.Error()<<" q0=   " <<Q0<<std::endl;
      msg.Error()<<" qmin= " <<qmin<<std::endl;
    } 
    return log_delta_q(Q);
    // Integrate NLL branching (using nummeric gauss integration);
  case NLL::inte_self : 
    {
      NLL::code savetype = type;
      double    saveqmin = qmin;
      double    saveqmax = qmax;

      type       = NLL::int_G_q;
      qmin       = Q0;
      qmax       = Q;
      double sum = (*this)(Q);

      qmin = saveqmin;
      qmax = saveqmax;
      type = savetype;

      return sum;
    }
    // Integrate NLL branching (using analytic NLO formula);
  case NLL::analytic :
    {
      if ((nlo==0.) && (K==0)) {
	// ensure that the integrated splitting function is positive!
	//	double Q0_eff = exp(3./4.) * Q0;
	double Q0_eff = Q0;

	double balpi = beta0*asmu/(4.*M_PI) ;
	double eta0  = balpi * log(sqr(Q0_eff)/mu2);
	double eta1  = balpi * log(sqr(Q)/mu2);
	return 8.*M_PI * CF/(sqr(beta0)* asmu) *
	  ( balpi * log(sqr(Q0_eff/Q)) + 
	    (1.+ balpi*(log(sqr(Q)/mu2) - 3./2.))*log((1+eta1)/(1+eta0)));
      } 
      else {
	msg.Error()<<"ERROR in NLL_Sudakov::IntGammaQ."<<std::endl
		   <<"    Analytic version for higher orders not implemented yet."<<std::endl;
	return 0.;
      }
    }
    // Integrate NLL branching (using analytic NLO formula with alphaS);
  case NLL::alphas_repl :
    {
      double balpi = beta0*asmu/(4.*M_PI) ;
      double eta0  = balpi * log(sqr(Q0)/mu2);
      double eta1  = balpi * log(sqr(Q)/mu2);
      return 8.*M_PI * CF/(sqr(beta0)* asmu) *
	       ( balpi * log(sqr(Q0/Q)) + 
		 (1.+ balpi*(log(sqr(Q)/mu2) - 3./2.)) * log(AlphaS(Q0*Q0)/AlphaS(Q*Q))  );
    }
  }

  msg.Error()<<" ERROR wrong type specified in NLL_Sudakov::IntGammaF  "<<std::endl;
  return 0.;
}

double NLL_Sudakov::IntGammaG(double Q0, double Q) {
  switch (inte) {
    // Integrate NLL branching (using nummeric gauss integration);
  case NLL::inte_self : 
    {
      NLL::code savetype = type;
      double    saveqmin = qmin;
      double    saveqmax = qmax;

      type       = NLL::int_G_g;
      qmin       = Q0;
      qmax       = Q;
      double sum = (*this)(Q);

      qmin = saveqmin;
      qmax = saveqmax;
      type = savetype;

      return sum;
    }
    // Integrate NLL branching (using analytic NLO formula);
  case NLL::analytic :
    {
      if ((nlo==0.) && (K==0)) {
	// ensure that the integrated splitting function is positive!
	//	double Q0_eff = exp(11./12.) * Q0;
	double Q0_eff = Q0;

	double balpi = beta0 * asmu/(4.*M_PI) ;
	double eta0  = balpi * log(sqr(Q0_eff)/mu2);
	double eta1  = balpi * log(sqr(Q)/mu2);
	return 8.* M_PI * CA/(sqr(beta0)* asmu) *
	  ( balpi * log(sqr(Q0_eff/Q)) + 
	    (1.+ balpi*(log(sqr(Q)/mu2) - 11./6.))*log((1+eta1)/(1+eta0)));
      }
      else {
	msg.Error()<<"ERROR in NLL_Sudakov::IntGammaQ."<<std::endl
		   <<"    Analytic version for higher orders not implemented yet."<<std::endl;
	return 0.;
      }
    }
    // Integrate NLL branching (using analytic NLO formula with alphaS);
  case NLL::alphas_repl :
    {
      double balpi = beta0*asmu/(4.*M_PI) ;
      double eta0  = balpi * log(sqr(Q0)/mu2);
      double eta1  = balpi * log(sqr(Q)/mu2);
      return 8.*M_PI * CA/(sqr(beta0)* asmu) *
	       ( balpi * log(sqr(Q0/Q)) + 
		 (1.+ balpi*(log(sqr(Q)/mu2) - 11./6.)) * log(AlphaS(Q0*Q0)/AlphaS(Q*Q)) );
    }
  }

  msg.Error()<<" ERROR wrong type specified in NLL_Sudakov::IntGammaF  "<<std::endl;
  return 0.;
}


double NLL_Sudakov::IntGammaF(double Q0, double Q) {
  switch (inte) {
    // Integrate NLL branching (using nummeric gauss integration);
  case NLL::inte_self : 
    {
      NLL::code savetype = type;
      double    saveqmin = qmin;
      double    saveqmax = qmax;

      type       = NLL::int_G_f;
      qmin       = Q0;
      qmax       = Q;
      double sum = (*this)(Q);

      qmin = saveqmin;
      qmax = saveqmax;
      type = savetype;

      return sum;
    }
    // Integrate NLL branching (using analytic LO formula with fixed alphaS);
  case NLL::analytic :
    {
      if ((nlo==0.) && (K==0)) {
	double eta0 = beta0*asmu/(4.*M_PI) * log(sqr(Q0)/mu2);
	double eta1 = beta0*asmu/(4.*M_PI) * log(sqr(Q)/mu2);
	return 2.*Nf/(3.*beta0) * log((1+eta1)/(1+eta0));
      }
      else {
	msg.Error()<<"ERROR in NLL_Sudakov::IntGammaQ."<<std::endl
		   <<"    Analytic version for higher orders not implemented yet."<<std::endl;
	return 0.;
      }
    }
    // Integrate NLL branching (using analytic LO formula with running alphaS);
  case NLL::alphas_repl :
    return 2.*Nf/(3.*beta0) * log(AlphaS(Q0*Q0)/AlphaS(Q*Q));
  }

  msg.Error()<<" ERROR wrong type specified in NLL_Sudakov::IntGammaF  "<<std::endl;
  return 0.;
}








double NLL_Sudakov::R2(double Q0, double Q)
{
  return sqr(DeltaQ(Q,Q0));
}

double NLL_Sudakov::R3(double Q0, double Q)
{
  qmin=Q0;
  qmax=Q;
  type=NLL::three_jet;
  double sum=gauss.Integrate(Q0,Q,1.e-5,1); 
  sum*=2.*R2(Q0,Q);
  return sum;
} 

double NLL_Sudakov::R4(double Q0, double Q) 
{
  qmin=Q0;
  qmax=Q;
  double fac=2.*R2(Q0,Q);
  type=NLL::three_jet;
  double sum1=gauss.Integrate(Q0,Q,1.e-5,1); 
  type=NLL::four_jet_1;
  double sum2=gauss.Integrate(Q0,Q,1.e-5,1); 
  return fac*(sqr(sum1)+sum2);
}

double NLL_Sudakov::R5(double Q0, double Q) 
{
  msg.Error()<<" ERROR: Not implemented !!"<<std::endl;
  return 0.;
}




double NLL_Sudakov::operator()(double x) {
  double save,sum;
  double zmin;

  switch (type) {
    // NLL branching fractions
  case NLL::G_g_G_f : return GammaG(x,qmax) + GammaF(x);
  case NLL::G_q     : return GammaQ(x,qmax);
  case NLL::G_g     : return GammaG(x,qmax);
  case NLL::G_f     : return GammaF(x);

    // integrals of branching fractions (exponents of Sudakovs)
  case NLL::int_G_q : 
    save=qmax;
    qmax=x;
    type=NLL::G_q;
    sum=gauss.Integrate(qmin,x,1.e-5,1); 
    type=NLL::int_G_q;
    qmax=save;
    return sum;
  case NLL::int_G_g_G_f :
    save=qmax;
    qmax=x;
    type=NLL::G_g_G_f;
    sum=gauss.Integrate(qmin,x,1.e-5,1);
    type=NLL::int_G_g_G_f;
    qmax=save;
    return sum;
  case NLL::int_G_g : 
    save=qmax;
    qmax=x;
    type=NLL::G_g;
    sum=gauss.Integrate(qmin,x,1.e-5,1); 
    type=NLL::int_G_g;
    qmax=save;
    return sum;
  case NLL::int_G_f : 
    save=qmax;
    qmax=x;
    type=NLL::G_f;
    sum=gauss.Integrate(qmin,x,1.e-5,1); 
    type=NLL::int_G_f;
    qmax=save;
    return sum;

    // LLA branching fractions - direct calculated from LLA Splitting functions
  case NLL::G_g_lla : 
    type  = NLL::P_g_lla;
    qstat = x;
    zmin  = 0.5-0.5*sqrt(1.-sqr(qmin/qstat));
    // factor 2 for transformation q->q2
    sum   = 2.*gauss.Integrate(zmin,1.-zmin,1.e-5,5)/qstat/(2.*M_PI);
    type  = NLL::G_g_lla;
    return sum;
  case NLL::G_q_lla : 
    type  = NLL::P_q_lla;
    qstat = x;
    zmin  = 0.5-0.5*sqrt(1.-sqr(qmin/qstat));
    // factor 2 for transformation q->q2
    sum   = 2.*gauss.Integrate(zmin,1.-zmin,1.e-5,5)/qstat/(2.*M_PI);
    type  = NLL::G_q_lla;
    return sum;

    // MLLA branching fractions - direct calculateted from MLLA Splitting functions
  case NLL::G_g_mlla: 
    type = NLL::P_g_mlla;
    qstat = x;
    zmin = qmin/(2.*qstat);
    // factor 2 for transformation q->q2
    sum = 2.*gauss.Integrate(zmin,1.-zmin,1.e-5,5)/qstat/(2.*M_PI);
    type = NLL::G_g_mlla;
    return sum;
  case NLL::G_q_mlla: 
    type = NLL::P_q_mlla;
    qstat = x;
    zmin = qmin/(2.*qstat);
    // factor 2 for transformation q->q2
    sum = 2.*gauss.Integrate(zmin,1.-zmin,1.e-5,5)/qstat/(2.*M_PI);
    type = NLL::G_q_mlla;
    return sum;

    // splitting functions
  case NLL::P_g_lla:
    return P_g_LLA(qstat,x);
  case NLL::P_q_lla: 
    return P_q_LLA(qstat,x);
  case NLL::P_g_mlla: 
    return P_g_MLLA(qstat,x);
  case NLL::P_q_mlla: 
    return P_q_MLLA(qstat,x);



    // NLL - Jetrates (in parts)
  case NLL::three_jet:
    return GammaQ(x,qmax)*DeltaQ(x,qmin);

  case NLL::four_jet_1:
    save=qmax;
    qmax=x;
    type=NLL::four_jet_2;
    sum=gauss.Integrate(qmin,x,1.e-5,1); 
    sum*=GammaQ(x,qmax)*DeltaG(x,qmin);
    type=NLL::four_jet_1;
    qmax=save;
    return sum;
  case NLL::four_jet_2:
    return GammaG(x,qmax)*DeltaG(x,qmin) +  GammaF(x)*DeltaF(x,qmin);
  default           : 
    msg.Error()<<" ERROR wrong type specified in NLL_Sudakov::operator() !!!!!  "<<std::endl;
  }
}














double NLL_Sudakov::P_q_MLLA(double q,double z)
{
  // factor 2 for q->qg + q->gq
  //MLLA
  return 2.*(*as)(z*z*(1.-z)*(1.-z)*q*q)*CF*(1.+z*z)/(1.-z);
}

double NLL_Sudakov::P_q_LLA(double q,double z)
{
  // factor 2 for q->qg + q->gq
  //LLA
  return 2.*(*as)((z-z*z)*q*q)*CF*(1.+z*z)/(1.-z);
}

double NLL_Sudakov::P_g_MLLA(double q,double z)
{  
 return 
   //MLLA
   // Nc = 3
   (*as)(z*z*(1.-z)*(1.-z)*q*q)*(6.*(1.-z+z*z)*(1.-z+z*z)/(z-z*z)
   // Nf=5
   +2.5*(2.*z*z+1.-2.*z));
}

double NLL_Sudakov::P_g_LLA(double q,double z)
{  
 return 
   //LLA
   // Nc = 3
   (*as)((z-z*z)*q*q)*(6.*(1.-z+z*z)*(1.-z+z*z)/(z-z*z)
   // Nf=5
   +2.5*(2.*z*z+1.-2.*z));
}

Kabbala NLL_Sudakov::DeltaQ(int j, double Q, double q)
{
  double value = DeltaQ(Q,q);
  msg.Tracking()<<" called NLL_Sudakov::DeltaQ("<<j<<", "<<Q<<", "<<q<<")"<<std::endl;
  msg.Tracking()<<"            result ="<<value<<std::endl<<std::endl;

  char help[10];
  sprintf(help,"Dq[%i]",j);
  return Kabbala(std::string(help),value);
}

Kabbala NLL_Sudakov::DeltaG(int j, double Q, double q)
{
  double value = DeltaG(Q,q);
  msg.Tracking()<<" called NLL_Sudakov::DeltaG("<<j<<", "<<Q<<", "<<q<<")"<<std::endl;
  msg.Tracking()<<"            result ="<<value<<std::endl<<std::endl;

  char help[10];
  sprintf(help,"Dg[%i]",j);
  return Kabbala(std::string(help),value);
}


Kabbala NLL_Sudakov::DeltaF(int j, double Q, double q)
{
  return DeltaQ(j,Q,q) * DeltaQ(j,Q,q)/DeltaG(j,Q,q);
}

double NLL_Sudakov::IntGammaI(double q, double i0, double i1, double i2, double i3) {
  /* original from mathematica
  double sum=(a2*i0*(108*Power(a0,2)*Power(a1,3)*a3 - 
      108*Power(a0,2)*Power(a1,3)*a2*i1 - 162*a0*Power(a1,2)*a2*a3*i1 + 
      19*a1*a2*Power(a3,2)*i1 - 108*Power(a0,2)*Power(a1,2)*a3*i2 + 
      108*Power(a0,2)*Power(a1,2)*a2*i1*i2 - 54*a0*a1*a2*a3*i1*i2 + 
      8*a2*Power(a3,2)*i1*i2 + 108*Power(a0,2)*Power(a1,2)*a3*i3 - 
      108*Power(a0,3)*Power(Log(q),4) + 
      108*Power(a0,2)*Power(a1,3)*(-(a2*i1) + a0*(a1 - i2 + i3))*
       Log(a1 - Log(q)) + 108*Power(a0,2)*Power(a1,3)*a3*
       Log(a0*(a1 - Log(q))) - 
      108*a0*Power(a1,2)*a2*a3*i1*Log(a0*(a1 - Log(q))) + 
      30*a1*a2*Power(a3,2)*i1*Log(a0*(a1 - Log(q))) - 
      108*Power(a0,2)*Power(a1,2)*a3*i2*Log(a0*(a1 - Log(q))) - 
      108*a0*a1*a2*a3*i1*i2*Log(a0*(a1 - Log(q))) + 
      24*a2*Power(a3,2)*i1*i2*Log(a0*(a1 - Log(q))) + 
      108*Power(a0,2)*Power(a1,2)*a3*i3*Log(a0*(a1 - Log(q))) + 
      54*Power(a0,2)*Power(a1,3)*a3*Power(Log(a0*(a1 - Log(q))),2) + 
      18*a1*a2*Power(a3,2)*i1*Power(Log(a0*(a1 - Log(q))),2) + 
      36*a2*Power(a3,2)*i1*i2*Power(Log(a0*(a1 - Log(q))),2) - 
      54*Power(a0,2)*Power(Log(q),3)*
       (-6*a0*a1 + 2*(-(a2*i1) + a0*(a1 - i2 + i3))*Log(a1 - Log(q)) + 
         a3*Power(Log(a0*(a1 - Log(q))),2)) + 
      54*a0*Power(Log(q),2)*(2*
          (-3*Power(a0,2)*Power(a1,2) - 2*a2*a3*i1 + 
            a0*(a1*a3 - a1*a2*i1 - a3*i2 + a2*i1*i2 + a3*i3)) + 
         6*a0*a1*(-(a2*i1) + a0*(a1 - i2 + i3))*Log(a1 - Log(q)) + 
         2*a3*(-2*a2*i1 + a0*(a1 - i2 + i3))*Log(a0*(a1 - Log(q))) + 
         3*a0*a1*a3*Power(Log(a0*(a1 - Log(q))),2)) - 
      27*Log(q)*(-4*Power(a0,3)*Power(a1,3) + 
         8*Power(a0,2)*Power(a1,2)*a3 - 8*Power(a0,2)*Power(a1,2)*a2*i1 - 
         14*a0*a1*a2*a3*i1 + a2*Power(a3,2)*i1 - 8*Power(a0,2)*a1*a3*i2 + 
         8*Power(a0,2)*a1*a2*i1*i2 - 2*a0*a2*a3*i1*i2 + 
         8*Power(a0,2)*a1*a3*i3 + 
         12*Power(a0,2)*Power(a1,2)*(-(a2*i1) + a0*(a1 - i2 + i3))*
          Log(a1 - Log(q)) + 2*a3*
          (a2*a3*i1 - 2*a0*a2*i1*(3*a1 + i2) + 
            4*Power(a0,2)*a1*(a1 - i2 + i3))*Log(a0*(a1 - Log(q))) + 
         2*a3*(3*Power(a0,2)*Power(a1,2) + a2*a3*i1)*
          Power(Log(a0*(a1 - Log(q))),2))))/
  (108.*Power(a0,4)*Power(a1 - Log(q),3))
  */
  /*
  double a02=a0*a0;
  double a03=a0*a02;
  double a12=a1*a1;
  double a13=a1*a12;
  double a32=a3*a3;
  double lnq=Log(q);
  double lnq2=lnq*lnq;
  double lnq3=lnq*lnq2;
  double lnq4=lnq*lnq3;
  double w = a0*(a1 - lnq);
  double lnw =Log(w);
  double lnw2=lnw*lnw; 

  double sum=(a2*i0*(108*a02*a13*a3 - 
      108*a02*a13*a2*i1 - 162*a0*a12*a2*a3*i1 + 
      19*a1*a2*a32*i1 - 108*a02*a12*a3*i2 + 
      108*a02*a12*a2*i1*i2 - 54*a0*a1*a2*a3*i1*i2 + 
      8*a2*a32*i1*i2 + 108*a02*a12*a3*i3 - 
      108*a03*lnq4 + 
      108*a02*a13*(-(a2*i1) + a0*(a1 - i2 + i3))*Log(a1 - lnq) + 
      108*a02*a13*a3*lnw - 
      108*a0*a12*a2*a3*i1*lnw + 
      30*a1*a2*a32*i1*lnw - 
      108*a02*a12*a3*i2*lnw - 
      108*a0*a1*a2*a3*i1*i2*lnw + 
      24*a2*a32*i1*i2*lnw + 
      108*a02*a12*a3*i3*lnw + 
      54*a02*a13*a3*lnw2 + 
      18*a1*a2*a32*i1*lnw2 + 
      36*a2*a32*i1*i2*lnw2 - 
      54*a02*lnq3*
       (-6*a0*a1 + 2*(-(a2*i1) + a0*(a1 - i2 + i3))*Log(a1 - lnq) + 
         a3*lnw2) + 
      54*a0*lnq2*(2*
          (-3*a02*a12 - 2*a2*a3*i1 + 
            a0*(a1*a3 - a1*a2*i1 - a3*i2 + a2*i1*i2 + a3*i3)) + 
         6*a0*a1*(-(a2*i1) + a0*(a1 - i2 + i3))*Log(a1 - lnq) + 
         2*a3*(-2*a2*i1 + a0*(a1 - i2 + i3))*lnw + 
         3*a0*a1*a3*lnw2) - 
      27*lnq*(-4*a03*a13 + 
         8*a02*a12*a3 - 8*a02*a12*a2*i1 - 
         14*a0*a1*a2*a3*i1 + a2*a32*i1 - 8*a02*a1*a3*i2 + 
         8*a02*a1*a2*i1*i2 - 2*a0*a2*a3*i1*i2 + 
         8*a02*a1*a3*i3 + 
         12*a02*a12*(-(a2*i1) + a0*(a1 - i2 + i3))*
          Log(a1 - lnq) + 2*a3*
          (a2*a3*i1 - 2*a0*a2*i1*(3*a1 + i2) + 
            4*a02*a1*(a1 - i2 + i3))*lnw + 
         2*a3*(3*a02*a12 + a2*a3*i1)*
          lnw2)))/
    (108.*a0*Power(w,3))
    return sum;
*/


  // new a_i definition!!!
/* denominator of alpha 
   beta0  = (11.*CA-2.*Nf)/3.;
   beta1  = (17.*CA*CA- 3.*CF*Nf-5.*CA*Nf)/3.;

  Nc     = 3;
  CA     = Nc;
  CF     = (Nc*Nc-1)/(2.*Nc);
  Nf     = 5;


    a0 = - beta0/2 
    a1 = Log[mu]

   w= (a1 -  Log[q])*a0*a2;

   alpha_s/pi 
    a2 = asmu/PI
    nlo* a3 = nlo * beta1/beta0 /2 
   alpi= a2/w * (1 - a2*a3 * Log[w]/w); 
*/ 
// WARNING! overlay Definitons of CA and Nf
  const double CA=3.;
  const double Nf=5.;
  const double a0=(2.*Nf-11.*CA)/6.;
  const double a3=0.5 *(17.*CA*CA- 3.*CF*Nf-5.*CA*Nf)/(11.*CA-2.*Nf);

  double a1=0.5*log(mu2);
  double a2=asmu/M_PI;

  const  double a02=a0*a0;
  const  double a03=a0*a02;
  double a12=a1*a1;
  double a13=a1*a12;
  const  double a32=a3*a3;
  double lnq=log(q);
  double lnq2=lnq*lnq;
  double lnq3=lnq*lnq2;
  double lnq4=lnq*lnq3;
  double w = a0*a2*(a1 - lnq);
  double lnw =log(w);
  double lnw2=lnw*lnw; 

  double sum= (-108*a03*i0*lnq4 - 
    54*a02*i0*lnq3*
     (-6*a0*a1 + 2*(-i1 + a0*(a1 - i2 + i3))*log(a1 - lnq) + nlo*a3*lnw2) + 
    54*a0*i0*lnq2*(
       2*(-3*a02*a12 - 2*nlo*a3*i1 + a0*(a1*nlo*a3 - a1*i1 - nlo*a3*i2 + i1*i2 + nlo*a3*i3)) + 
       6*a0*a1*(-i1 + a0*(a1 - i2 + i3))*log(a1 - lnq) + 
       2*nlo*a3*(-2*i1 + a0*(a1 - i2 + i3))*lnw + 
       3*a0*a1*nlo*a3*lnw2) - 
    27*i0*lnq*(
       -4*a03*a13 + 8*a02*a12*nlo*a3 - 8*a02*a12*i1 - 
       14*a0*a1*nlo*a3*i1 + nlo*a32*i1 - 8*a02*a1*nlo*a3*i2 + 
        8*a02*a1*i1*i2 - 2*a0*nlo*a3*i1*i2  + 8*a02*a1*nlo*a3*i3 + 
        12*a02*a12*(-i1 + a0*(a1 - i2 + i3))*log(a1 - lnq) + 
        2*nlo*a3*(nlo*a3*i1 - 2*a0*i1*(3*a1 + i2) + 4*a02*a1*(a1 - i2 + i3))*lnw + 
        2*nlo*a3*(3*a02*a12 + nlo*a3*i1)*lnw2) 
    + i0*(108*a02*a13*nlo*a3 - 108*a02*a13*i1 - 162*a0*a12*nlo*a3*i1 + 19*a1*nlo*a32*i1 - 
       108*a02*a12*nlo*a3*i2 + 108*a02*a12*i1*i2 - 54*a0*a1*nlo*a3*i1*i2 + 
       8*nlo*a32*i1*i2 + 108*a02*a12*nlo*a3*i3 + 
       108*a02*a13*(-i1 + a0*(a1 - i2 + i3))*log(a1 - lnq) + 
       6*nlo*a3*(-18*a0*a1*i1*(a1 + i2) + nlo*a3*i1*(5*a1 + 4*i2) + 
       18*a02*a12*(a1 - i2 + i3))*lnw + 
	  18*nlo*a3*(3*a02*a13 + nlo*a3*i1*(a1 + 2*i2))*lnw2))/(108.*a0*pow(w,3));
  return sum;
}
