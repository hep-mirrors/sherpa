#include "MODEL/Main/Running_AlphaS.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/MyStrStream.H"

namespace MODEL {
  Running_AlphaS * as =0;
  double Running_AlphaS::s_Nc    = 3.;
  double Running_AlphaS::s_CA    = s_Nc;
  double Running_AlphaS::s_CF    = (s_Nc*s_Nc-1.)/(2.*s_Nc);
  double Running_AlphaS::s_TF    = 1./2.;
  double Running_AlphaS::s_zeta2 = M_PI*M_PI/6.;
  double Running_AlphaS::s_zeta3 = 1.202056903;
}

using namespace MODEL;
using namespace ATOOLS;
using namespace std;

namespace MODEL {

  std::ostream &operator<<(std::ostream &str,const AS_Data &data)
  {
    str<<" AS_Data for scale->["<<data.mass2_low<<","<<data.mass2_up<<"]"
       <<" as->["<<data.as_low<<","<<data.as_up<<"]"<<std::endl
       <<"             nf->"<<data.nF<<" lam2->"<<data.lambda2<<std::endl
       <<"             bet0->"<<data.beta0<<" bet1->"<<data.beta1
       <<" bet2->"<<data.beta2<<" bet3->"<<data.beta3;
    return str;
  }

}

Running_AlphaS::Running_AlphaS(const double as_MZ,const double m2_Z,
			       const int order, const double fac) :
  m_asMZ(as_MZ), m_M2Z(m2_Z), m_order(2), 
  m_fac(fac), m_ascut(4.*M_PI), m_cutq2(1.)
{
  if(m_fac==1.0 && rpa.gen.Variable("RENORMALIZATION_SCALE_FACTOR")!="") {
    m_fac=ToType<double>(rpa.gen.Variable("RENORMALIZATION_SCALE_FACTOR"));
  }
  if (m_fac!=1.0) msg_Debugging()<<METHOD<<"(): Setting scale factor "<<m_fac<<"\n";
  FixASData();
}

Running_AlphaS::~Running_AlphaS()
{
  while (!m_asdata.empty()) {
    delete (*m_asdata.begin());
    m_asdata.erase(m_asdata.begin());
  }
}


void Running_AlphaS::FixASData() {
  AS_Data * asdata = new AS_Data;
  asdata->nF        = 3;
  asdata->mass2_low = 0.;
  asdata->mass2_up  = sqr(Flavour(kf_c).PSMass());
  asdata->beta0     = Beta0(asdata->nF);
  asdata->beta1     = Beta1(asdata->nF);
  asdata->beta2     = Beta2(asdata->nF);
  asdata->beta3     = Beta3(asdata->nF);
  m_asdata.insert(asdata);
  asdata = new AS_Data;
  asdata->nF        = 4;
  asdata->mass2_low = sqr(Flavour(kf_c).PSMass());
  asdata->mass2_up  = sqr(Flavour(kf_b).PSMass());
  asdata->beta0     = Beta0(asdata->nF);
  asdata->beta1     = Beta1(asdata->nF);
  asdata->beta2     = Beta2(asdata->nF);
  asdata->beta3     = Beta3(asdata->nF);
  m_asdata.insert(asdata);
  asdata = new AS_Data;
  asdata->nF        = 5;
  asdata->mass2_low = sqr(Flavour(kf_b).PSMass());
  asdata->mass2_up  = sqr(Flavour(kf_t).PSMass());
  asdata->beta0     = Beta0(asdata->nF);
  asdata->beta1     = Beta1(asdata->nF);
  asdata->beta2     = Beta2(asdata->nF);
  asdata->beta3     = Beta3(asdata->nF);
  m_asdata.insert(asdata);
  asdata = new AS_Data;
  asdata->nF        = 6;
  asdata->mass2_low = sqr(Flavour(kf_t).PSMass());
  asdata->mass2_up  = sqr(1.e20);
  asdata->beta0     = Beta0(asdata->nF);
  asdata->beta1     = Beta1(asdata->nF);
  asdata->beta2     = Beta2(asdata->nF);
  asdata->beta3     = Beta3(asdata->nF);
  m_asdata.insert(asdata);

  AS_Data * refregion(FindRegion(m_M2Z));
  FixLambda(refregion,m_M2Z,m_asMZ);
  asdata_iterator up,down;
  for (up=m_asdata.begin();up!=m_asdata.end();up++) {
    if ((*up)==refregion) break;
  }
  if (up!=m_asdata.end() && up!=m_asdata.begin()) {
    double asfactor(1.),ashigh;
    do {
      down     = up;
      down--;
      asfactor = MatchingFactorDown((*up),(*down));
      ashigh   = asfactor*(*up)->as_low;
      FixLambda((*down),(*down)->mass2_up,ashigh);
      up--;
      if (IsNan((*down)->as_low)) (*down)->as_low = 1.e30;
    } while (up!=m_asdata.begin());
  }
  FixSoftEnd();
  for (down=m_asdata.begin();down!=m_asdata.end();down++) 
    if ((*down)==refregion) break;
  up = down;
  up++;
  if (up!=m_asdata.end() && down!=m_asdata.begin()) {
    double asfactor(1.),aslow;
    do {
      asfactor = MatchingFactorUp((*down),(*up));
      aslow   = asfactor*(*down)->as_up;
      FixLambda((*up),(*up)->mass2_low,aslow);
      down = up;
      up++;
    } while (up!=m_asdata.end());
  }

  //std::cout<<"Found fixed region."<<std::endl<<(**down)<<std::endl;
  //std::cout<<"  Reference region."<<std::endl<<(*refregion)<<std::endl;
  //SelfTest();
}

void Running_AlphaS::FixLambda(AS_Data * region,
			       const double scale,const double asvalue) const
{
  //   std::cout<<"In "<<METHOD<<" for order = "<<m_order<<"."<<std::endl
  // 	   <<"   Fix lambda2 for as("<<sqrt(scale)<<") = "<<asvalue
  // 	   <<" in the region ("<<sqrt(region->mass2_low)<<" - "<<sqrt(region->mass2_up)<<" GeV)."
  // 	   <<std::endl
  // 	   <<"   nF = "<<region->nF<<" and Beta0 = "<<region->beta0<<", Beta1 = "<<region->beta1 
  // 	   <<", Beta2 = "<<region->beta2<<", Beta3 = "<<region->beta3<<"."<<std::endl;

  // Different conventions for the coefficients of the beta-function
  double beta0(region->beta0/4.), a(asvalue/M_PI);
  double exponent(1/a),C_MSbar(0.);
  if (m_order>0) {
    double b1(region->beta1/16./beta0);
    exponent += b1*log(a);
    C_MSbar   = b1/beta0*log(beta0);
    if (m_order>1) {
      double b2(region->beta2/64./beta0), b1_2(b1*b1);
      exponent += (b2-b1_2)*a;
      if (m_order>2) {
	double b3(region->beta3/256./beta0), b1_3(b1_2*b1);
	exponent += (b3/2.-b1*b2+b1_3/2.)*a*a;
      }
    }
  }
  exponent += C_MSbar;
  double lambda2 = scale*exp(-exponent/beta0);
  region->lambda2 = lambda2;
  
  double ascheck = LambdaForm(scale,region,m_order);
  //std::cout<<"Check this : as(fixed) = "<<asvalue<<" vs. as(lambda) = "<<ascheck<<"."<<std::endl
  //	   <<"   log(mu2/lambda2) = "<<log(scale/lambda2)<<" = "<<(-exponent/beta0)<<"."<<std::endl;
  while (dabs((ascheck-asvalue)/(ascheck+asvalue))>1.e-6) {
    double factor   = (1.-(ascheck-asvalue)/(ascheck+asvalue));
    region->lambda2 *= factor;
    ascheck         = LambdaForm(scale,region,m_order);
  }

  region->as_low = LambdaForm(region->mass2_low,region,m_order);
  region->as_up  = LambdaForm(region->mass2_up,region,m_order);
  //std::cout<<"Check this : as(fixed) = "<<asvalue<<" vs. as(lambda) = "<<ascheck<<"."<<std::endl
  //	   <<"   log(mu2/lambda2) = "<<log(scale/lambda2)<<" = "<<(-exponent/beta0)<<"."<<std::endl;
}

void Running_AlphaS::FixSoftEnd() {
  // make sure, three flavour running stops at m_ascut, update soft end of this regime.
  AS_Data * softend((*m_asdata.begin()));
  softend->mass2_low = softend->lambda2*(1.00000001);
  softend->as_low    = LambdaForm(softend->mass2_low,softend,m_order);

  m_cutq2 = Inverse(m_ascut);
  //std::cout<<"Check this : Cutoff m_cutq2 = "<<m_cutq2<<" --> as("<<sqrt(m_cutq2)<<") = "
  //	   <<LambdaForm(m_cutq2,(*m_asdata.begin()),m_order)<<"."<<std::endl;
  softend->mass2_low = m_cutq2;
  softend->as_low    = LambdaForm(m_cutq2,softend,m_order);
  // add a new nonperturbative regime, frozen.
  AS_Data * nonpert  = new AS_Data;
  nonpert->nF        = -1;
  nonpert->lambda2   = -1.;
  nonpert->as_low    = nonpert->as_up = softend->as_low;
  nonpert->mass2_low = 0.;
  nonpert->mass2_up  = softend->mass2_low;
  m_asdata.insert(nonpert);
}

AS_Data * Running_AlphaS::FindRegion(const double q2) const {
  double Q2(dabs(q2));
  asdata_iterator diter = m_asdata.begin();
  while (diter!=m_asdata.end()) {
    if (Q2>(*diter)->mass2_low && Q2<=(*diter)->mass2_up) return (*diter);
    diter++;
  }
  diter = m_asdata.end();
  diter--;
  return (*diter);
}

const int Running_AlphaS::Nf(const double q2) const
{
  return FindRegion(q2)->nF;
}

const double Running_AlphaS::Beta0(const double q2) const
{
  return Beta0(FindRegion(q2)->nF);
}

const double Running_AlphaS::Beta1(const double q2) const
{
  return Beta1(FindRegion(q2)->nF);
}

const double Running_AlphaS::Beta2(const double q2) const
{
  return Beta2(FindRegion(q2)->nF);
}

const double Running_AlphaS::Beta3(const double q2) const
{
  return Beta3(FindRegion(q2)->nF);
}

const double Running_AlphaS::CutQ2() const 
{ 
  return m_cutq2; 
}

const double Running_AlphaS::Beta0(const int nF) const
{
  return 11./3.*s_CA - 4./3.*s_TF*nF;
}

const double Running_AlphaS::Beta1(const int nF) const
{
  return 34./3.*sqr(s_CA) - 4.*s_CF*s_TF*nF - 20./3.*s_CA*s_TF*nF;
    
  // Check: 102.0-12.6667*nF
}

const double Running_AlphaS::Beta2(const int nF) const
{
  return 
    2857./54.*pow(s_CA,3.) + 2.*sqr(s_CF)*s_TF*nF - 
    205./9.*s_CF*s_CA*s_TF*nF - 1415./27.*sqr(s_CA)*s_TF*nF +
    44./9.*s_CF*sqr(s_TF*nF) + 158./27.*s_CA*sqr(s_TF*nF);

  // Check: 1428.50-279.611*nF+6.01852*sqr(nF)
}

const double Running_AlphaS::Beta3(const int nF) const
{
  return 
    (150653./486.-44./9.*s_zeta3)*pow(s_CA,4.) + 
    (-39143./81.+136./3.*s_zeta3)*pow(s_CA,3.)*s_TF*nF +
    (7073./243.-656./9.*s_zeta3)*pow(s_CA,2.)*s_CF*s_TF*nF +
    (-4204./27.-352./9.*s_zeta3)*pow(s_CF,2.)*s_CA*s_TF*nF +
    (7930./81.+224./9.*s_zeta3)*pow(s_CA,2.)*s_TF*nF +
    (1352./27.-704./9.*s_zeta3)*pow(s_CF,2.)*s_TF*nF +
    (17152./243.+448./9.*s_zeta3)*s_CF*s_CA*s_TF*nF +
    46.*pow(s_CF,3.)*s_TF*nF +
    424./243.*s_CA*pow(s_TF*nF,3.) +
    1232./243.*s_CF*pow(s_TF*nF,3.) +
    (-80./9.+704./3.*s_zeta3)*sqr(s_Nc)*(sqr(s_Nc)+36.)/24. +
    (512./9.-1664./3.*s_zeta3)*s_Nc*(sqr(s_Nc)+6.)/48.*nF +
    (-704./9.+512./3.*s_zeta3)*(pow(s_Nc,4.)-6.*sqr(s_Nc)+18.)/(96*sqr(s_Nc))*sqr(nF);
    
  // Check: 29243.0-6946.30*nF+405.089*sqr(nF)+1.49931*pow(nF,3.)
}

const double Running_AlphaS::MatchingFactorDown(const AS_Data * up,const AS_Data * down,
						const double refmass) const
{
  double factor(1.), a_up(up->as_low/M_PI), L(refmass<0.?0.:2.*log(down->mass2_up/refmass));
  if (m_order>0) {
    factor += -L/6.*a_up;
    if (m_order>1) {
      double L_2(L*L), a_up_2(a_up*a_up);
      double C2(-7./24.);
      factor += (L_2/36.-19./24.*L+C2)*a_up_2;
      if (m_order>2) {
	double L_3(L_2*L), a_up_3(a_up_2*a_up);
	int nL(down->nF);
	double C3(-80507./27648*s_zeta3-2./3.*s_zeta2*(log(2)/3.+1.)
		  -58933./124416.+double(nL)/9.*(s_zeta2+2479./3456.));
	factor += (-L_3/216.-131./576.*L_2+(-8521.+409.*double(nL))*L/1728.+C3)*a_up_3;
      }
    }
  }
  return factor;
}

const double Running_AlphaS::MatchingFactorUp(const AS_Data * down,const AS_Data * up,
					      const double refmass) const
{
  double factor(1.), a_down(down->as_up/M_PI), L(refmass<0.?0.:2.*log(up->mass2_low/refmass));
  if (m_order>0) {
    factor += +L/6.*a_down;
    if (m_order>1) {
      double L_2(L*L), a_down_2(a_down*a_down);
      double C2(+7./24.);
      factor += (L_2/36.+19./24.*L+C2)*a_down_2;
      if (m_order>2) {
	double L_3(L_2*L), a_down_3(a_down_2*a_down);
	int nL(up->nF);
	double C3(80507./27648*s_zeta3+2./3.*s_zeta2*(log(2)/3.+1.)+
		  58933./124416.-double(nL)/9.*(s_zeta2+2479./3456.));
	factor += (L_3/216.+511./576.*L_2+(8941.-409.*double(nL))*L/1728.+C3)*a_down_3;
      }
    }
  }
  return factor;
}

const double Running_AlphaS::operator()(const double q2) const {
  return LambdaForm(q2,FindRegion(q2),m_order);
}

const double  Running_AlphaS::AlphaS(const double q2) const {
  return operator()(q2);
}

const double Running_AlphaS::LambdaForm(const double q2,const AS_Data * asregion,
					const int order) const 
{
  if (asregion->nF<0) {
    return q2/(asregion->mass2_up-asregion->mass2_low)*
      (asregion->as_up-asregion->as_low)+asregion->as_low;
  }
  double Lambda2(asregion->lambda2), L(log(q2/Lambda2));
  double beta0(asregion->beta0/4.), beta0L(beta0*L);
  double a(1/beta0L);
  if (order>0) {
    double b1(asregion->beta1/16./beta0);
    double LogL(log(L)), beta0L_2(beta0L*beta0L);
    a += -b1*LogL/beta0L_2;
    if (order>1) {
      double b2(asregion->beta2/64./beta0);
      double b1_2(b1*b1), LogL_2(LogL*LogL), beta0L_3(beta0L_2*beta0L);
      a += (b1_2*(LogL_2-LogL-1.)+b2)/beta0L_3;
      if (order>2) {
	double b3(asregion->beta2/256./beta0); 
	double b1_3(b1_2*b1), LogL_3(LogL_2*LogL), beta0L_4(beta0L_3*beta0L);
	a += (b1_3*(-LogL_3+2.5*LogL_2+2.*LogL-0.5)+3.*b1*b2*LogL+0.5*b3)/beta0L_4;
      }
    }
  }
  return M_PI*a;
}

const double Running_AlphaS::Inverse(const double as) const {
  //std::cout<<"In "<<METHOD<<" for as = "<<as<<"."<<std::endl;
  asdata_iterator diter=m_asdata.begin();
  while (diter!=m_asdata.end()) {
    if ((*diter)->as_low>=as && (*diter)->as_up<=as) break;
    diter++;
  }
  if (diter==m_asdata.end()) {
    msg_Error()<<"Error in "<<METHOD<<":"<<std::endl
	       <<"   Inversion of alpha_S for as = "<<as<<" yielded no result."<<std::endl
	       <<"   Will return 0 and hope for the best."<<std::endl;
    return 0.;
  }
  double low((*diter)->mass2_low), up((*diter)->mass2_up), test((low+up)/2.);
  double astest(LambdaForm(test,(*diter),m_order));
  //std::cout<<"  ["<<low<<" ... "<<test<<" ... "<<up<<"] --> "<<astest<<std::endl;
  while (dabs((astest-as)/(astest+as))>1.e-6) {
    if (astest>as) low = test;
              else up = test;
    test   = (low+up)/2.;
    astest = LambdaForm(test,(*diter),m_order);
    //std::cout<<"  ["<<low<<" ... "<<test<<" ... "<<up<<"] --> "<<astest<<" vs. "<<as<<std::endl;
  }
  //std::cout<<"  ... out with "<<test<<"("<<sqrt(test)<<")."<<std::endl;
  return test;
}

void Running_AlphaS::SelfTest() const {
  int np=100001;

  double smin=0.0025;
  double smax=500*500;

  double mult=::exp(log(smax/smin)/double(np-2));
  double s= smin;

  double rs, as0, as1, as2, as3;
  AS_Data * region;
  std::ofstream was;
  was.open("as_test.dat");
  for (int i=1;i<np;++i) {
    rs =sqrt(s);
    region = FindRegion(s);
    as0 = LambdaForm(s,region,0);
    as1 = LambdaForm(s,region,1);
    as2 = LambdaForm(s,region,2);
    as3 = LambdaForm(s,region,3);
    was<<rs<<" "<<as0<<" "<<as1<<" "<<as2<<" "<<as3<<std::endl;
    //  msg_Out()<<"  Q = "<<rs<<" alphaS(Q) = "<<as<<std::endl;
    s*=mult;
  }
   was.close();
  //exit(1);
}

