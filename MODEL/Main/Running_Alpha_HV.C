#include "MODEL/Main/Running_Alpha_HV.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/MyStrStream.H"

#include <algorithm>

namespace MODEL {
  Running_Alpha_HV * as_HV =0;
}

using namespace MODEL;
using namespace ATOOLS;
using namespace std;

Running_Alpha_HV::Running_Alpha_HV(const double as_MZ,const double m2_MZ,const int order, const double fac, const std::string group, const double Nc) :   m_order(order), m_as_MZ(as_MZ), m_m2_MZ(m2_MZ), m_fac(fac), m_group(group), m_Nc(Nc)
{
  if(m_fac==1.0 && rpa.gen.Variable("RENORMALIZATION_SCALE_FACTOR")!="") {
    m_fac=ToType<double>(rpa.gen.Variable("RENORMALIZATION_SCALE_FACTOR"));
  }
  if (m_fac!=1.0) msg_Debugging()<<METHOD<<"(): Setting scale factor "<<m_fac<<"\n";
  m_type  = std::string("Running Coupling");
  p_thresh  = NULL;
  
  if (group==std::string("SU")) {
    m_CF = (m_Nc*m_Nc-1.)/(2.*m_Nc);        
    m_CA = m_Nc;           
    m_TR = 1./2.;
  }
  else if (group==std::string("SO")) {
    m_CF = 1./2.*(m_Nc-1);        
    m_CA = (m_Nc-2.);        
    m_TR=1;
    if (m_Nc==3) {
      m_CF *= 2.;
      m_CA *= 2.;
      m_TR *= 2.;
    }
  }
  else if (group==std::string("SP")) {
    m_CF = 1./4.*(m_Nc+1);
    m_CA = 0.5*(m_Nc+2.);     
    m_TR = 1./2.;
  }
  else {
    msg_Error()<<" Gauge Group not supported assume SU(3) instead! "<<std::endl; 
    m_Nc = 3.;
    m_CF = (m_Nc*m_Nc-1.)/(2.*m_Nc);        
    m_CA = m_Nc;           
    m_TR = 0.5;
  }

  //------------------------------------------------------------
  // HV thresholds for SU(N) interactions
  //------------------------------------------------------------
  m_nth = 0;

  for(KFCode_ParticleInfo_Map::const_iterator kfit(s_kftable.begin());
      kfit!=s_kftable.end();++kfit) {
    if (kfit->first>=9900001&&kfit->first<=9900021 && Flavour(kfit->first).Strong()) {
      m_nth++;
    }
  }

  p_thresh        = new AsDataSet[m_nth+1]; 
  double * masses = new double[m_nth];
  int count = 0;
  for(KFCode_ParticleInfo_Map::const_iterator kfit(s_kftable.begin());
      kfit!=s_kftable.end();++kfit) {
    if (kfit->first>=9900001 && kfit->first<=9900021) {
      Flavour flav(kfit->first);
      if (flav.Strong()) {
	masses[count] = sqr(flav.Mass(true));
	count++;
      }
    }
  }
  
  std::vector<double> sortmass(&masses[0],&masses[m_nth]);
  std::sort(sortmass.begin(),sortmass.end(),std::less<double>());
  for (int i(0);i<m_nth;++i) masses[i]=sortmass[i];

  int j   = 0; 
  m_mzset = 0;
  for (int i=0; i<m_nth; ++j) {
    if ((masses[i]>m_m2_MZ)&&(!m_mzset)) {
      //insert Z boson (starting point for any evaluation) 
      m_mzset               = j;
      p_thresh[j].low_scale = m_m2_MZ;
      p_thresh[j].as_low    = m_as_MZ;
      p_thresh[j].nf        = i-1;
    } 
    else {
      p_thresh[j].low_scale = masses[i];
      p_thresh[j].as_low    = 0.;
      p_thresh[j].nf        = i;
      ++i;
    }
    if (j>0) {
      p_thresh[j-1].high_scale = p_thresh[j].low_scale;
      p_thresh[j-1].as_high    = p_thresh[j].as_low;
    }
  }
  if (!m_mzset) {
    int j                    = m_nth;
    m_mzset                  = j;
    p_thresh[j].low_scale    = m_m2_MZ;
    p_thresh[j].as_low       = m_as_MZ;
    p_thresh[j-1].high_scale = m_m2_MZ;
    p_thresh[j-1].as_high    = m_as_MZ;
    p_thresh[j].nf           = p_thresh[j-1].nf;      
  }
  p_thresh[m_nth].high_scale   = 1.e20;  
  p_thresh[m_nth].as_high      = 0.;

  delete [] masses;

  for (int i=m_mzset;i<=m_nth;++i) {
    Lambda2(i);
    p_thresh[i].as_high       = Alpha_HV_Lam(p_thresh[i].high_scale,i);
    if (i<m_nth) {
      p_thresh[i+1].as_low    = p_thresh[i].as_high *
	InvZetaOS2(p_thresh[i].as_high,p_thresh[i].high_scale,p_thresh[i].high_scale,p_thresh[i].nf);
    }
  }

  for (int i=m_mzset-1;i>=0;--i) {
    double lam2               = Lambda2(i);
    p_thresh[i].as_low        = Alpha_HV_Lam(p_thresh[i].low_scale,i);
    if ((lam2>p_thresh[i].low_scale) || (p_thresh[i].as_low>1.)) ContinueAlpha_HV(i);
    else {
      if (i>0) {
	p_thresh[i-1].as_high = p_thresh[i].as_low *
	  ZetaOS2(p_thresh[i].as_low,p_thresh[i].low_scale,p_thresh[i].low_scale,p_thresh[i-1].nf);
      }
    }
  }
}


Running_Alpha_HV::~Running_Alpha_HV()
{
  if (p_thresh!=0) { delete [] p_thresh; p_thresh = NULL; }
}

double Running_Alpha_HV::Beta0(const int nf) {
  return 1./4. * (11./3.*m_CA - (4./3.)*m_TR*nf);
}

double Running_Alpha_HV::Beta1(const int nf) {
  return 1./16. * (34./3.*m_CA*m_CA - (10./3.*m_CA*m_TR+4.*m_CF*m_TR*nf));
}

double Running_Alpha_HV::Lambda2(const int nr) {
  double as  = p_thresh[nr].as_low; 
  double mu2 = p_thresh[nr].low_scale;
  if (as==0.) {
    as  = p_thresh[nr].as_high;
    mu2 = p_thresh[nr].high_scale;
  }

  const double a   = as/M_PI;

  int    & nf      = p_thresh[nr].nf;
  double & beta0   = p_thresh[nr].beta0;
  double * b       = p_thresh[nr].b;
  double & lambda2 = p_thresh[nr].lambda2;
  
  // calculate beta coefficients
  beta0 = Beta0(nf);
  b[1]  = Beta1(nf)/beta0;

  double betaL = 1./a;
  if (m_order>=1) {
    betaL     += b[1]*log(a);
  }

  lambda2         = ::exp(-betaL/beta0)*mu2;
  double tas1     = Alpha_HV_Lam(mu2,nr);
  double dlambda2 = 1.e-8;
  if (dabs(tas1-as)/as>1.e-11) {
    for (;(dabs(tas1-as)/as>1.e-11);) {
      lambda2     = lambda2+dlambda2;
      double tas2 = Alpha_HV_Lam(mu2,nr);
      dlambda2    = (as-tas2)/(tas2-tas1)*dlambda2;
      tas1        = tas2;
    }
  }

  return lambda2;
}

double Running_Alpha_HV::Alpha_HV_Lam(const double Q2,const int nr)
{
  // using shorter names
  double & beta0   = p_thresh[nr].beta0;
  double *  b      = p_thresh[nr].b;
  double & lambda2 = p_thresh[nr].lambda2;
  double L         = log(Q2/lambda2);
  double pref      = 1./(beta0*L);

  double a         = pref;
  if (m_order==0) return M_PI*a;

  double logL     = log(L);
  pref           *=1./(beta0*L);
  a              += -pref*(b[1] * logL);
  return M_PI*a;
}

double Running_Alpha_HV::ZetaOS2(const double as,const double mass2_os,
			       const double mu2,const int nl) {
  double zeta2g = 1.;

  // 0th order
  if (m_order==0) return zeta2g;

  // 1st order (one loop) corresponds to two loop lambda  
  double L      = log(mu2/mass2_os);
  double a      = as/M_PI;
  zeta2g       += - a*1./6.*L;
  return zeta2g; 
}

double Running_Alpha_HV::InvZetaOS2(const double as,const double mass2_os,
				  const double mu2,const int nl) {
  // might be simplified considerably when using mu2==mass2
  double zeta2g  = 1.;
  // 0th order   
  if (m_order==0) return zeta2g;

  // 1st order (one loop) corresponds to two loop lambda
  double L      = log(mu2/mass2_os);
  double a      = as/M_PI;
  zeta2g       += + a*1./6.*L;
  return zeta2g;
}



void Running_Alpha_HV::ContinueAlpha_HV(int & nr) {
  // shrink actual domain
  //  * to given t0        or
  //  * to alphaS=alphaCut
  double alpha_cut = 1.;
  double & beta0   = p_thresh[nr].beta0;
  double & lambda2 = p_thresh[nr].lambda2;
  double t0        = lambda2 * ::exp(M_PI/(alpha_cut*beta0));
  double as        = Alpha_HV_Lam(t0,nr);
  for (;dabs(as-alpha_cut)>1.e-8;) {
    double t1      = t0+0.00001;
    double as1     = Alpha_HV_Lam(t1,nr);
    double das     = (as -as1)/(t0-t1);
    t1             = (alpha_cut-as)/das + t0;
    t0             = t1;
    as             = Alpha_HV_Lam(t0,nr);    
  }

  m_cutq2 = t0;

  // modify lower domains
  p_thresh[nr].low_scale    = t0;
  p_thresh[nr-1].high_scale = t0;
  p_thresh[nr].as_low       = as;
  p_thresh[nr-1].as_high    = as;

  for (int i = nr-1; i>=0; --i) {
    p_thresh[i].nf          = -1;  // i.e. no ordinary running !!!
    p_thresh[i].lambda2     = 0.;
    p_thresh[i].as_low      = p_thresh[i].as_high/p_thresh[i].high_scale*p_thresh[i].low_scale;
    if (i>0) p_thresh[i-1].as_high=p_thresh[i].as_low;
  }
  nr =0;
}



double Running_Alpha_HV::operator()(double q2)
{
  double as;
  q2=q2*m_fac;
  if (q2<0.) q2=-q2;
  int i = m_mzset-1;
  if (q2<=m_m2_MZ) {
    for (;!((p_thresh[i].low_scale<q2)&&(q2<=p_thresh[i].high_scale));--i) {
      if (i<=0) break;
    }
    if (p_thresh[i].nf>=0) 
      as = Alpha_HV_Lam(q2,i);
    else 
      as = q2/p_thresh[i].high_scale * p_thresh[i].as_high;
  }
  else {
    ++i;
    for (;!((p_thresh[i].low_scale<q2)&&(q2<=p_thresh[i].high_scale));++i) {
      if (i>=m_nth) break;
    }
    as   = Alpha_HV_Lam(q2,i);
  }
  return as;
}  

double  Running_Alpha_HV::Alpha_HV(const double q2){
  return operator()(q2);
}

int Running_Alpha_HV::Nf(const double sc) 
{
  double q2(sc*m_fac);
  for (int i=0;i<=m_nth;++i) {
    if (q2<=p_thresh[i].high_scale && q2>p_thresh[i].low_scale )
      return p_thresh[i].nf;
  }
  return m_nth;
}

