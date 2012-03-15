#include "SHRiMPS/Eikonals/Form_Factors.H"
#include "SHRiMPS/Tools/Special_Functions.H"
#include "ATOOLS/Math/Gauss_Integrator.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/Histogram.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"


using namespace SHRIMPS;
using namespace ATOOLS;

Special_Functions SHRIMPS::SF;


inline double Form_Factor::Norm_Argument::operator()(double q) { 
  return 2.*M_PI*q*(*p_ff)(q); 
}

inline double Form_Factor::FT_Argument::operator()(double q) { 
  return 2.*M_PI*q*SF.Jn(0,q*m_b)*(*p_ff)(q); 
}




Form_Factor::Form_Factor(const int & number,const int & test) :
  m_form(MBpars.FF_Form()), m_number(number), 
  m_Lambda2(0.), m_beta(0.), m_kappa(0.), m_xi(0.), 
  m_prefactor(0.), m_ftnorm(4.*M_PI*M_PI), m_norm(0.),
  m_bmax(0.), m_deltab(1.), m_bsteps(100), 
  m_ffmin(0.), m_ffmax(0.), m_accu(0.0001), m_test(test)
{
  switch (m_test) {
  case 1: 
    m_form = ff_form::Gauss;
    break;
  case -1:
    m_form = ff_form::dipole;
    break;
  default:
    break;
  }
}

Form_Factor::Form_Factor(const ff_form::code & fff,const int & number,
			 const int & test) :
  m_form(fff), m_number(number), 
  m_Lambda2(0.), m_beta(0.), m_kappa(0.), m_xi(0.), 
  m_prefactor(0.), m_ftnorm(4.*M_PI*M_PI), m_norm(0.),
  m_bmax(0.), m_deltab(1.), m_bsteps(100), 
  m_ffmin(0.), m_ffmax(0.), m_accu(0.0001), m_test(test)
{
  switch (m_test) {
  case 1: 
    m_form = ff_form::Gauss;
    break;
  case -1:
    m_form = ff_form::dipole;
    break;
  default:
    break;
  }
}

void Form_Factor::Initialise(const std::vector<double> & params) {
  m_prefactor = params[0]; 
  if ((m_form==ff_form::dipole && params.size()<7) ||
      (m_form==ff_form::Gauss && params.size()<6)) {
    msg_Error()<<"Error in "<<METHOD<<":"<<std::endl
	       <<"    Wrong number of parameters ("<<params.size()
	       <<") for form factor form "<<m_form<<", will abort."<<std::endl;
    exit(1);
  }

  switch (m_form) {
  case ff_form::dipole: 
    m_xi      = params[4];
  case ff_form::Gauss:
    m_Lambda2 = params[1];
    m_beta    = params[2];
    m_kappa   = params[3];
    break;
  default:
    msg_Error()<<"Error in "<<METHOD<<":"<<std::endl
	       <<"    Form factor form "<<m_form
	       <<" not known, will abort."<<std::endl;
    exit(1);
  }
  m_norm   = NormAnalytical();
  m_bmax   = ATOOLS::Min(16.,params[params.size()-2]);
  m_deltab = (m_bmax/double(m_bsteps));
  m_accu   = params[params.size()-1];
  m_ftarg  = FT_Argument(this);

  FillFTGrid();

  msg_Info()<<"Initialised form factor "<<m_number<<".\n";
  msg_Tracking()
    <<"  "<<m_form<<"(Lambda2 = "<<m_Lambda2<<", "
    <<"kappa = "<<m_kappa<<", xi = "<<m_xi<<")"<<std::endl
    <<"  beta = "<<m_beta<<" --> Norm = "<<m_norm<<","<<std::endl
    <<"  evaluate in (naively) "<<m_bsteps
    <<" steps up to b = "<<m_bmax<<", "
    <<"accuracy goal = "<<m_accu<<"."<<std::endl;
  
  if (m_test!=0) TestFormFactor();
}

void Form_Factor::FillFTGrid() {
  msg_Tracking()<<"In "<<METHOD
		<<" for interval [0"<<", "<<m_bmax<<"] "
		<<"with delta b = "<<m_deltab<<"."<<std::endl;
  m_ffmin = 1.e12;
  m_ffmax = 0.;
  double b(0.), value;
  bool   run(true);
  while (run) {
    while (b<=m_bmax) {
      value = CalculateFourierTransform(b);
      if (dabs(value)<1.e-6 || 
	  (m_ffmax>0. && dabs(value/m_ffmax)<1.e-7)) value  = 0.;
      if (value<m_ffmin) m_ffmin = value;
      if (value>m_ffmax) m_ffmax = value;
      msg_Debugging()<<METHOD<<": fill in FT("<<b<<") = "<<value
		     <<" ("<<m_values.size()<<")."<<std::endl;
      m_values.push_back(value);
      b += m_deltab;
    }
    run = false;
    //if (dabs(m_ffmin/m_ffmax)>1.e-3 && m_bmax<12.) {
    //  msg_Tracking()<<"   - does not meet accuracy goal in b_max = "
    //<<m_bmax<<", double it."<<std::endl
    //		    <<"     Use now "<<m_bsteps<<" steps --> delta_b = "
    //<<m_deltab
    //		    <<"     (min, max = "<<m_ffmin<<", "<<m_ffmax<<" -> "
    //<<m_ffmin/m_ffmax<<")."<<std::endl;
    //  m_bsteps *= 2;
    //  m_bmax   *= 2.;
    //  run = true;
    //}
    if (!run) {
      for (size_t i=1;i<m_values.size()-1;i++) {
	double btest((i+0.5)*m_deltab);
	double fit(FourierTransform(btest));
	double exact(Max(0.,CalculateFourierTransform(btest)));
	if (exact/m_ffmax>1.e-6 && dabs(fit/exact-1.)>m_accu/5. && 
	    dabs(fit/m_ffmax)>m_accu) {
	  msg_Tracking()<<"   - does not meet accuracy goal yet "
			<<"("<<(dabs(fit/exact-1.)>0.01)<<") "
			<<"in "<<m_bsteps<<" steps:"<<std::endl
			<<"     i = "<<i<<", "
			<<"b = "<<btest<<": "<<dabs(fit/exact-1.)
			<<" from : exact = "<<exact<<" vs. grid = "<<fit<<"."
			<<std::endl
			<<"     Use now "<<(2*m_bsteps)
			<<" steps --> delta_b = "<<(m_deltab/2.)<<"."
			<<std::endl;
	  m_ffmin   = 1.e12;
	  m_ffmax   = -1.;
	  m_bsteps *= 2;
	  m_deltab /= 2.;
	  b         = 0.;
	  run       = true;
	  break;
	}
      }
      if (run) m_values.clear();
    }
  }
  msg_Tracking()<<"Out "<<METHOD<<": accuracy goal ("<<m_accu<<") reached."
		<<std::endl
		<<"   B-Interval: [0"<<", "<<m_bmax<<"] in "
		<<m_bsteps<<"steps, "
		<<"yields form factors in ["<<m_ffmin<<", "<<m_ffmax<<"]."
		<<std::endl
		<<"   Set the value for b_max  ==>  F(b_max) = 0."<<std::endl;
  m_values[m_values.size()-1] = 0.;
}

double Form_Factor::CalculateFourierTransform(const double & b) {
  double ft(0.), diff(1.), qmin(0.), qmax(10.);
  m_ftarg.SetB(b);
  Gauss_Integrator gauss(&m_ftarg);
  while (dabs(diff)>1.e-8) {
    diff  = gauss.Integrate(qmin,qmax,sqr(m_accu));
    ft   += diff;
    qmin  = qmax;
    qmax *= 2.;
  }
  if (dabs(ft)<1.e-6) ft = 0.;
  return ft/m_ftnorm;  
}

double Form_Factor::FourierTransform(const double & b) const 
{
  double ft(0.);
  if (b<0. || b>m_bmax) {
    msg_Error()<<"Warning in "<<METHOD<<":"<<std::endl
	       <<"   Impact parameter b = "<<b<<" outside interval"
	       <<" [0"<<", "<<m_bmax<<"];"<<std::endl
	       <<"   Will return 0 and hope for the best."<<std::endl;
    return 0.;
    abort();
  }
  else {
    size_t bbin(size_t(b/m_deltab));
    if (bbin>=0 && bbin<m_bsteps) {
      if (dabs(b-bbin*m_deltab)/m_deltab<1.e-3) ft = m_values[bbin];
      else if (bbin>=1 && bbin<m_values.size()-2) {
	double ft1(m_values[bbin-1]), b1=(bbin-1)*m_deltab;
	double ft2(m_values[bbin+0]), b2=(bbin+0)*m_deltab;
	double ft3(m_values[bbin+1]), b3=(bbin+1)*m_deltab;
	double ft4(m_values[bbin+2]), b4=(bbin+2)*m_deltab;
	ft =
	  ft1 * (b-b2)*(b-b3)*(b-b4)/((b1-b2)*(b1-b3)*(b1-b4)) +
	  ft2 * (b-b1)*(b-b3)*(b-b4)/((b2-b1)*(b2-b3)*(b2-b4)) +
	  ft3 * (b-b1)*(b-b2)*(b-b4)/((b3-b1)*(b3-b2)*(b3-b4)) +
	  ft4 * (b-b1)*(b-b2)*(b-b3)/((b4-b1)*(b4-b2)*(b4-b3));	
      }
      else if (bbin>=0 && bbin<m_values.size()-1) {
	double ft1(m_values[bbin]),   b1=bbin*m_deltab;
	double ft2(m_values[bbin+1]), b2=(bbin+1)*m_deltab;
	ft = (ft1*(b2-b) + ft2*(b-b1))/m_deltab;
      }
    }
    if (ft<0.) ft = 0.;
    return ft;
  }
  return 0.;
}

double Form_Factor::ImpactParameter(const double & val) const 
{
  // assuming a monotonously decreasing function
  if (val>m_values.front()) {
    msg_Error()<<"Warning in "<<METHOD<<":"<<std::endl
	       <<"   Fourier Transform ft = "<<val<<" outside interval"
	       <<" ["<<m_values.front()<<", "<<m_values.back()<<"]."<<std::endl
	       <<"   Will return 0 and hope for the best."<<std::endl;
    return 0.;
  }
  if (val<m_values.back()) {
    msg_Error()<<"Warning in "<<METHOD<<":"<<std::endl
	       <<"   Fourier Transform ft = "<<val<<" outside interval"
	       <<" ["<<m_values.front()<<", "<<m_values.back()<<"]."<<std::endl
	       <<"   Will return b_max = "<<m_bmax
	       <<" and hope for the best."<<std::endl;
    return m_bmax;
  }
  size_t i;
  for (i=0;i<m_bsteps;i++) { if (m_values[i]<val) break; }
  double b2(i*m_deltab), b1(b2-m_deltab); 
  double val2(m_values[i]), val1(m_values[i-1]);
  return b1 * (val-val2)/(val1-val2) + b2 * (val-val1)/(val2-val1);
}


double Form_Factor::operator()(const double q) {
  double pref(sqr(m_beta)*(1.+m_kappa));
  double q2tilde_Lambda2((1.+m_kappa)*(q*q/m_Lambda2)), ff(0.);
  switch (m_form) {
  case ff_form::Gauss:
    ff = exp(-q2tilde_Lambda2);
    break;
  case ff_form::dipole:
    ff = exp(-m_xi*q2tilde_Lambda2)/sqr(1.+q2tilde_Lambda2);
    break;
  default:
    break;
  }
  if (ff<1.e-6) ff=0.;
  return pref * ff;
}

double Form_Factor::
SelectQT2(const double & qt2max,const double & qt2min) const {
  double qt2(0.), pref(m_Lambda2/(1.+m_kappa)), effexp(1./pref),random;
  switch (m_form) {
  case ff_form::Gauss:
    do {
      random = ATOOLS::ran->Get();
      qt2 = -pref*log(1.-random*(1.-exp(-qt2max/pref)));
    } while (qt2<qt2min);
    break;
  case ff_form::dipole:
    do {
      random = ATOOLS::ran->Get();
      //Original form factor
      //exp(-xi* ...)/(1+q^2/Lambda^2)^2
      qt2 = (pref*qt2max*random)/(qt2max*(1.-random)+pref);
    } while (qt2<qt2min || exp(-m_xi*effexp*qt2)<ATOOLS::ran->Get());
    break;
  default:
    break;
  }
  return qt2;
}

double Form_Factor::NormAnalytical() {
  double norm(m_beta*m_beta*M_PI*m_Lambda2/m_ftnorm);
  switch (m_form) {
  case ff_form::Gauss:
    break;
  case ff_form::dipole:
    norm *= (1.-exp(m_xi)*m_xi*SF.IncompleteGamma(0,m_xi));
    break;
  default: 
    norm = 0.;
    break;
  }
  return norm;
}

double Form_Factor::AnalyticalFourierTransform(const double & b) {
  double kernel(0.), pref(m_beta*m_beta*m_Lambda2*M_PI/m_ftnorm), help(0.);
  switch (m_form) {
  case ff_form::Gauss:
    kernel = exp(-b*b*m_Lambda2/(4.*(1.+m_kappa)));
    break;
  case ff_form::dipole:
    if (b<=1.e-8) kernel = 1.;
    else {
      help   = sqrt((1.+m_kappa)/m_Lambda2);
      kernel = b/help * SF.Kn(1,b/help);
    }
    break;
  default:
    break;
  }
  return pref * kernel;
}

double Form_Factor::Norm() {
  Norm_Argument normarg(this);
  double norm(0.), diff(0.), rel(1.), qmin(0.), qmax(1.);
  Gauss_Integrator gauss(&normarg);
  while (rel>m_accu) {
    diff  = gauss.Integrate(qmin,qmax,m_accu,1);
    norm += diff;
    rel   = dabs(diff/norm);
    qmin  = qmax;
    qmax *= 2.;
  }
  return norm/m_ftnorm;
}


void Form_Factor::TestFormFactor() {
  if (m_test!=0) {
    msg_Out()<<"In "<<METHOD<<"("<<m_test<<") : "<<std::endl
	     <<"   Formfactor(0) = "<<operator()(0.)<<" "
	     <<"vs. "<<(m_beta*m_beta*(1.+m_kappa))<<"."<<std::endl
	     <<"   Norm = 1/(2 Pi)^2 Int_0^Infinity dq [q 2 Pi f(q)] = "
	     <<Norm()<<" vs. analytical = "<<NormAnalytical()<<std::endl
	     <<"                                        vs. estimate = "
	     <<AnalyticalFourierTransform(0.)
	     <<" from approximate FT(0)."<<std::endl
	     <<"   Fourier transform for b = 0  : exact = "
	     <<CalculateFourierTransform(0.)<<","
	     <<" approximately = "<<AnalyticalFourierTransform(0.)<<","
	     <<std::endl
	     <<"                     for b = 1  : exact = "
	     <<CalculateFourierTransform(1.)<<","
	     <<" approximately = "<<AnalyticalFourierTransform(1.)<<","
	     <<std::endl
	     <<"                     for b = 10 : exact = "
	     <<CalculateFourierTransform(10.)<<","
	     <<" approximately = "<<AnalyticalFourierTransform(10.)<<"."
	     <<std::endl
	     <<"   Grid in impact parameter space: "<<m_bsteps<<" bins "
	     <<"up to bmax = "<<m_bmax<<"."<<std::endl;
    double qt2max(20.),q,b,ff;
    std::ofstream was;
    std::string filename = std::string("InclusiveQuantities/formfactors-ana.dat");
    was.open(filename.c_str());
    was<<"# q     form factor"<<std::endl;
    for (int i=0;i<100;i++) { 
      q  = sqrt(qt2max*double(i)/100.);
      ff = (*this)(q);
      was<<" "<<q<<"  "<<ff<<std::endl;
    }
    was<<std::endl<<std::endl;
    was<<"# b     FT of form factor num      ana"<<std::endl;
    for (int i=0;i<100;i++) { 
      b  = m_bmax*double(i)/100.;
      was<<" "<<b<<"   "<<CalculateFourierTransform(b)<<"   "<<AnalyticalFourierTransform(b)<<std::endl;
    }
    was.close();
    if (m_test==-1) {
      PrintFFGrids(0);
      exit(1);
    }
    else if (m_test==-2) {
      msg_Out()<<"In "<<METHOD<<": "
	       <<"testing analytical functions, will exit afterwards."
	       <<std::endl;
      double b;
      for (int t=1;t<10;t++) {
	b = ran->Get();
	msg_Out()<<"   Exp[LnGamma("<<t<<") = "
		 <<exp(SF.LnGamma(double(t)))<<"   "
		 <<"Gamma(0, "<<b<<") = "<<SF.IncompleteGamma(0.,b)<<"."
		 <<std::endl;
      }
      SF.TestBessel();
      exit(1);
    }
    else if (m_test==-3) {
//       double qt2max(20.),q,ff;
      ATOOLS::Histogram histo1(0,0.0,qt2max,100);
      for (int i=0;i<100000;i++) 
	histo1.Insert(SelectQT2(qt2max)); 
      histo1.Finalize();
      histo1.Output(std::string("SelectQt2.dat"));
      std::ofstream was;
      std::string filename = std::string("Analytical.dat");
      was.open(filename.c_str());
      for (int i=0;i<100;i++) { 
	q  = sqrt(qt2max*double(i)/100.);
	ff = (*this)(q);
	was<<" "<<q*q<<"  "<<ff<<std::endl;
      }
      was.close();
      exit(1);
    }
  }
}

void Form_Factor::PrintFFGrids(const int & mode) {
  std::ofstream was;
  std::string filename, tag;

  tag = std::string("all");
  Form_Factor dipana(ff_form::dipole,10,0);
  Form_Factor diporig(ff_form::dipole,11,0);
  Form_Factor diphalf(ff_form::dipole,12,0);
  Form_Factor dipGauss(ff_form::Gauss,13,0);
  Form_Factor dipmin(ff_form::dipole,14,0);
  Form_Factor dipzero(ff_form::dipole,15,0);
  std::vector<double> params;
  params.push_back(m_prefactor);
  params.push_back(m_Lambda2);
  params.push_back(m_beta);
  params.push_back(m_kappa);
  params.push_back(m_xi);
  params.push_back(m_bmax);
  params.push_back(m_accu);
  diporig.Initialise(params);
  dipGauss.Initialise(params);
  params[4] = 1.e-6;
  dipana.Initialise(params);
  params[4] = 0.5;
  diphalf.Initialise(params);
  params[3] = -m_kappa;
  params[4] = m_xi;
  dipmin.Initialise(params);
  params[3] = 0.;
  dipzero.Initialise(params);
  
  double q(0.);
  msg_Out()<<"   Form factor in Q space: "<<std::endl;
  filename = std::string("Form_Factor_In_Qspace.")+tag+std::string(".dat");
  was.open(filename.c_str());
  for (int qstep=0;qstep<10000;qstep++) {
    q = qstep*1./1000.;
    was<<" "<<q<<"   "<<diporig(q)<<"   "<<dipana(q)<<"   "
       <<diphalf(q)<<"   "<<dipGauss(q)<<std::endl;
  }
  was.close();
  
  msg_Out()<<"   Form factor in B space: "<<std::endl
	   <<" b   orig   xi->0  analytic   Gauss   analytic"<<std::endl;  
  filename = std::string("Form_Factor_In_Bspace.")+tag+std::string(".dat");
  was.open(filename.c_str());
  double b(0.), val1, val2, val2a, val3, val3a;
  while (b<=8.) {
    val1  = diporig.FourierTransform(b);
    val2  = dipana.FourierTransform(b);
    val2a = dipana.AnalyticalFourierTransform(b);
    val3  = dipGauss.FourierTransform(b);
    val3a = dipGauss.AnalyticalFourierTransform(b);
    was<<" "<<b<<"   "<<val1
       <<"   "<<val2<<"   "<<val2a<<" ("<<(100.*(1.-val2/val2a))<<")    "
       <<"   "<<val3<<"   "<<val3a<<" ("<<(100.*(1.-val3/val3a))<<")"
       <<std::endl;
    if (b<=4.) {
      msg_Out()<<" "<<b<<"   "<<val1
	       <<"   "<<val2<<"   "<<val2a<<" "
	       <<"("<<(100.*(1.-val2/val2a))<<"%)    "
	       <<"   "<<val3<<"   "<<val3a<<" "
	       <<"("<<(100.*(1.-val3/val3a))<<"%)"<<std::endl;
      b += m_deltab/10.;
    }
    else b+= m_deltab/5.;
  }
  was.close();


  msg_Out()<<"   Form factor in B space, dependence on kappa: "<<std::endl
	   <<" b   orig   kappa=0   kappa->-kappa"<<std::endl;  
  filename = std::string("Form_Factor_In_Bspace_kappa.")+tag+
    std::string(".dat");
  was.open(filename.c_str());
  b = 0.;
  while (b<=8.) {
    val1  = diporig.FourierTransform(b);
    val2  = dipzero.FourierTransform(b);
    val3  = dipmin.FourierTransform(b);
    was<<" "<<b<<"   "<<val1<<"   "<<val2<<"   "<<val3<<std::endl;
    if (b<=4.) {
      msg_Out()<<" "<<b<<"   "<<val1<<"   "<<val2<<"   "<<val3<<std::endl;
      b += m_deltab/10.;
    }
    else b+= m_deltab/5.;
  }
  was.close();
}





