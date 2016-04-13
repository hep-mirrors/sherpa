#include "SHRiMPS/Eikonals/Form_Factors.H"
#include "SHRiMPS/Tools/MinBias_Parameters.H"
#include "SHRiMPS/Tools/Special_Functions.H"
#include "ATOOLS/Math/Gauss_Integrator.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/Histogram.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"


using namespace SHRIMPS;
using namespace ATOOLS;

Special_Functions SHRIMPS::SF;

double Form_Factor::Norm_Argument::operator()(double q) { 
  return 2.*M_PI*q*(*p_ff)(q); 
}

double Form_Factor::FT_Argument::operator()(double q) { 
  return 2.*M_PI*q*SF.Jn(0,q*m_b)*(*p_ff)(q); 
}


Form_Factor::Form_Factor(const FormFactor_Parameters & params) :
  m_ftarg(this), 
  m_number(params.number), m_form(params.form), 
  m_prefactor(params.norm), m_beta(sqrt(params.beta02)),
  m_Lambda2(params.Lambda2), m_kappa(params.kappa), m_xi(params.xi), 
  m_bmax(params.bmax), m_bsteps(16), m_deltab(m_bmax/double(m_bsteps)), 
  m_accu(10.*params.accu), m_ffmin(1.e-8), m_ffmax(0.), 
  m_ftnorm(4.*M_PI*M_PI)
{ }

void Form_Factor::Initialise() {  
  m_norm   = NormAnalytical();
  m_ftarg  = FT_Argument(this);
  FillFourierTransformGrid();
}

void Form_Factor::FillFourierTransformGrid() {
  msg_Tracking()<<METHOD<<"(bmax = "<<m_bmax<<") "
		<<"with "<<m_bsteps<<" initial steps,\n"
		<<"   accuracy goal = "<<m_accu<<";\n"
		<<"   start evaluating in (naively) "<<(2*m_bsteps)
		<<" steps up to b = "<<m_bmax<<".\n";
  do {
    m_bsteps *= 2;
    m_deltab /= 2.;
    FillTestGrid();
  } while (!GridGood());
  m_values[m_values.size()-1] = 0.;
  m_ffmax = CalculateFourierTransform(0.);
}

void Form_Factor::FillTestGrid() {
  msg_Tracking()<<METHOD<<" for "<<m_bsteps
		<<" steps with size = "<<m_deltab<<"\n";
  m_values.clear();
  double b(0.), value;
  while (b<=m_bmax) {
    value = CalculateFourierTransform(b);
    if (m_ffmax>0. && dabs(value/m_ffmax)<m_ffmin) value  = 0.;
    if (value>m_ffmax) m_ffmax = value;
    m_values.push_back(value);
    b += m_deltab;
  }
}

bool Form_Factor::GridGood() {
  for (size_t i=1;i<m_values.size()-1;i++) {
    double btest((i+0.5)*m_deltab);
    double fit(FourierTransform(btest));
    double exact(Max(0.,CalculateFourierTransform(btest)));
    if (exact/m_ffmax>1.e-6 && dabs(fit)/m_ffmax>1.e-6 &&
	dabs(fit/exact-1.)>m_accu/5.) {
      msg_Tracking()<<"   - does not meet accuracy goal yet "
	       <<"("<<(dabs(fit/exact-1.)>0.01)<<") "
	       <<"in "<<m_bsteps<<" steps:\n"
	       <<"     i = "<<i<<", "
	       <<"b = "<<btest<<": "<<dabs(fit/exact-1.)
	       <<" from : exact = "<<exact<<" vs. grid = "<<fit<<".\n"
	       <<"     Use now "<<(2*m_bsteps)
	       <<" steps --> delta_b = "<<(m_deltab/2.)<<".\n";
      return false;
    }
  }
  msg_Tracking()<<"   - did meet accuracy goal.\n";
  return true;
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
  if (ft<m_ffmin) ft = 0.;
  return ft/m_ftnorm;
}

double Form_Factor::FourierTransform(const double & b) const 
{
  double ft(0.);
  if (b<0. || b>m_bmax) return 0.;
  size_t bbin(size_t(b/m_deltab));
  if (bbin<m_bsteps) {
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
    else if (bbin<m_values.size()-1) {
      double ft1(m_values[bbin]),   b1=bbin*m_deltab;
      double ft2(m_values[bbin+1]), b2=(bbin+1)*m_deltab;
      ft = (ft1*(b2-b) + ft2*(b-b1))/m_deltab;
    }
  }
  if (ft<0.) ft = 0.;
  return ft;
}

double Form_Factor::ImpactParameter(const double & val) const 
{
  if (val>m_values.front()) return 0.;
  if (val<m_values.back())  m_bmax;

  size_t i;
  for (i=0;i<m_bsteps;i++) { if (m_values[i]<val) break; }
  double b2(i*m_deltab),    b1(b2-m_deltab); 
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


void Form_Factor::Test(const std::string & dirname) {
  TestSpecialFunctions(dirname);
  WriteOutFF_Q(dirname);
  WriteOutFF_B(dirname);
  TestNormAndSpecificBs(dirname);
  TestQ2Selection(dirname);
}

void Form_Factor::TestSpecialFunctions(const std::string & dirname) {  
  std::ofstream was1,was2;
  std::string filename1 = dirname+std::string("/LnGamma.dat");
  was1.open(filename1.c_str());
  std::string filename2 = dirname+std::string("/IncompleteGamma.dat");
  was2.open(filename2.c_str());
  for (int t=1;t<51;t++) {
    was1<<t<<"  "<<SF.LnGamma(double(t))<<"\n";
    was2<<(double(t)/50.)<<"   "
	<<SF.IncompleteGamma(0.,double(t)/50.)<<"\n";
  }
  was1.close();
  was2.close();
  SF.TestBessel(dirname);
}

void Form_Factor::WriteOutFF_Q(const std::string & dirname) {  
  double qt2max(20.),q,ff;
  std::ofstream was;
  std::string filename = dirname+std::string("/FormFactor_Q.dat");
  was.open(filename.c_str());
  was<<"# q     FF(q^2)\n";
  for (int i=0;i<100;i++) { 
    q  = sqrt(qt2max*double(i)/100.);
    ff = (*this)(q);
    was<<" "<<q<<"  "<<ff<<"\n";
  }
  was.close();
}  

void Form_Factor::WriteOutFF_B(const std::string & dirname) {  
  double b;
  std::ofstream was;
  std::string filename = dirname+std::string("/FormFactor_B.dat");
  was.open(filename.c_str());
  was<<"# b     FT of form factor num      ana"<<std::endl;
  for (int i=0;i<100;i++) { 
    b  = m_bmax*double(i)/100.;
    was<<" "<<b<<"   "<<CalculateFourierTransform(b)<<"   "
       <<AnalyticalFourierTransform(b)<<"\n";
  }
  was.close();
}

void Form_Factor::TestNormAndSpecificBs(const std::string & dirname) {
  std::ofstream was;
  std::string filename = dirname+std::string("/FormFactor_Summary.dat");
  was.open(filename.c_str());
  was<<"Formfactor(0, "<<m_form<<") = "<<operator()(0.)<<" "
     <<"vs. "<<(m_beta*m_beta*(1.+m_kappa))<<".\n"
     <<"Norm = 1/(2 Pi)^2 Int_0^Infinity dq [q 2 Pi f(q)] = "
     <<Norm()<<"\n"
     <<"                                   vs. analytical = "
     <<NormAnalytical()<<"\n"
     <<"                                   vs. estimate   = "
     <<AnalyticalFourierTransform(0.)<<" from approximate FT(0).\n";
  was<<"Fourier transform for b : exact : analytical : interpolated\n";
  for (size_t i=0;i<11;i++) {
    const double b = double(i);
    was<<"  "<<b<<"  "<<CalculateFourierTransform(b)<<"  "
       <<AnalyticalFourierTransform(b)<<"  "
       <<FourierTransform(b)<<"\n";
  }
   was<<"Grid in impact parameter space: "<<m_bsteps<<" bins "
      <<"up to bmax = "<<m_bmax<<", will be in separate file for plotting.\n";
  was.close();
}
  

void Form_Factor::TestQ2Selection(const std::string & dirname) {
  double qt2max(16.);
  ATOOLS::Histogram histo1(0,0.0,qt2max,100);
  for (int i=0;i<100000;i++) 
    histo1.Insert(SelectQT2(qt2max)); 
  histo1.Finalize();
  histo1.Output(dirname+std::string("/SelectQt2.dat"));
  std::ofstream was;
  std::string filename = dirname+std::string("/FF_Q2_Analytical.dat");
  was.open(filename.c_str());
  double q,ff;
  for (int i=0;i<100;i++) { 
    q  = sqrt(qt2max*double(i)/100.);
    ff = (*this)(q);
    was<<" "<<q*q<<"  "<<ff<<std::endl;
  }
  was.close();
}

/*
void Form_Factor::PrintFFGrids(const int & mode) {
  std::string tag("all");
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
  std::string filename = dirname+
    std::string("Form_Factor_In_Qspace.")+tag+std::string(".dat");
  std::ofstream was;
  was.open(filename.c_str());
  was<<"#   Form factor in Q space: \n";
  for (int qstep=0;qstep<10000;qstep++) {
    q = qstep*1./1000.;
    was<<" "<<q<<"   "<<diporig(q)<<"   "<<dipana(q)<<"   "
       <<diphalf(q)<<"   "<<dipGauss(q)<<"\n";
  }
  was.close();
  
  filename = dirname+
    std::string("/Form_Factor_In_Bspace.")+tag+std::string(".dat");
  was.open(filename.c_str());
  was<<"#   Form factor in B space: \n"
     <<"# b   orig   xi->0  analytic     xi=0.5    Gauss   analytic  \n";  
  double b(0.), val1, val2, val2a, val2b, val3, val3a;
  while (b<=8.) {
    val1  = diporig.FourierTransform(b);
    val2  = dipana.FourierTransform(b);
    val2a = dipana.AnalyticalFourierTransform(b);
    val2b = diphalf.FourierTransform(b);
    val3  = dipGauss.FourierTransform(b);
    val3a = dipGauss.AnalyticalFourierTransform(b);
    was<<" "<<b<<"   "<<val1
       <<"   "<<val2<<"   "<<val2a<<" ("<<(100.*(1.-val2/val2a))<<")    "
       <<"   "<<val2b
       <<"   "<<val3<<"   "<<val3a<<" ("<<(100.*(1.-val3/val3a))<<")"
       <<std::endl;
    if (b<=4.) b += m_deltab/10.;
    else b+= m_deltab/5.;
  }
  was.close();


  filename = std::string("Form_Factor_In_Bspace_kappa.")+tag+
    std::string(".dat");
  was.open(filename.c_str());
  was<<"#   Form factor in B space, dependence on kappa: \n"
     <<"# b   orig   kappa=0   kappa->-kappa\n";  
  b = 0.;
  while (b<=8.) {
    val1  = diporig.FourierTransform(b);
    val2  = dipzero.FourierTransform(b);
    val3  = dipmin.FourierTransform(b);
    was<<" "<<b<<"   "<<val1<<"   "<<val2<<"   "<<val3<<std::endl;
    if (b<=4.) b += m_deltab/10.;
    else b+= m_deltab/5.;
  }
  was.close();
}
*/




