#include "MODEL/Main/Running_AlphaS.H"
#include "PDF/Main/PDF_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/MyStrStream.H"

#include <algorithm>

namespace MODEL {
  Running_AlphaS * as =0;
}

using namespace MODEL;
using namespace ATOOLS;
using namespace PDF;
using namespace std;

namespace MODEL {

  std::ostream &operator<<(std::ostream &str,AsDataSet &set)
  {
    str<<"scale->["<<set.low_scale<<","<<set.high_scale<<"]";
    str<<" as->["<<set.as_low<<","<<set.as_high<<"]";
    str<<" nf->"<<set.nf<<" lam2->"<<set.lambda2;
    str<<" bet0->"<<set.beta0;
    return str;
  }

}

One_Running_AlphaS::One_Running_AlphaS(const double as_MZ,const double m2_MZ,
			       const int order, const int thmode,
			       PDF::PDF_Base *const aspdf,
				       One_Running_AlphaS *const mo) : 
  m_order(order), m_pdf(0),
  m_as_MZ(as_MZ), m_m2_MZ(m2_MZ),
  p_pdf(aspdf), p_sas(NULL)
{
  p_thresh  = NULL;

  m_CF    = 4./3.;        
  m_CA    = 3.;           

  //------------------------------------------------------------
  // SM thresholds for strong interactions, i.e. QCD
  //------------------------------------------------------------
  m_nth = 0;
  for(KFCode_ParticleInfo_Map::const_iterator kfit(s_kftable.begin());
      kfit!=s_kftable.end()&&kfit->first<=21;++kfit) {
    if (Flavour(kfit->first).Strong()) m_nth++;
  }

  p_thresh        = new AsDataSet[m_nth+1]; 
  double * masses = new double[m_nth];
  
  int count = 0;
  for(KFCode_ParticleInfo_Map::const_iterator kfit(s_kftable.begin());
      kfit!=s_kftable.end()&&kfit->first<=21;++kfit) {
    Flavour flav(kfit->first);
    if (flav.Strong()) {
      masses[count] = sqr(flav.Mass(thmode!=0));
      count++;
    }
  }

  if (p_pdf) {
    Data_Reader dataread(" ",";","!","=");
    dataread.AddComment("#");
    dataread.AddWordSeparator("\t");
    if (dataread.GetValue<int>("OVERRIDE_PDF_INFO",0)==1) {
      msg_Error()<<om::bold<<METHOD<<"(): "<<om::reset<<om::red
		 <<"Overriding \\alpha_s information from PDF. "
		 <<"Make sure you know what you are doing!"
		 <<om::reset<<std::endl;
    }
    else {
      const PDF::PDF_AS_Info &info(p_pdf->ASInfo());
      if (info.m_order>=0) {
      if (mo==NULL) {
      m_order=info.m_order;
      m_as_MZ=info.m_asmz;
      }
      if (dataread.GetValue<int>("USE_PDF_ALPHAS",0)==1) m_pdf=1;
      /*
      m_nth=info.m_flavs.size()+1;
      for (int i(0);i<m_nth;++i) {
	masses[i]=sqr(info.m_flavs[i].m_mass);
      }
      masses[m_nth-1]=0.0;
      */
      if (mo && m_as_MZ!=mo->m_as_MZ && m_pdf)
	THROW(fatal_error,"Cannot use PDF alphas to vary \\mu_R");
      if (mo==NULL || !IsEqual(m_as_MZ,mo->m_as_MZ)) {
      msg_Info()<<METHOD<<"() {\n  Setting \\alpha_s according to PDF\n"
		<<"  perturbative order "<<m_order
		<<"\n  \\alpha_s(M_Z) = "<<m_as_MZ;
      // msg_Info<<"\n  quark masses = { ";
      // for (int i(0);i<m_nth-1;++i) msg_Info()<<sqrt(masses[i])<<" ";
      msg_Info()<<"\n}"<<std::endl;
      }
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
    p_thresh[i].as_high       = AlphaSLam(p_thresh[i].high_scale,i);
    if (i<m_nth) {
      p_thresh[i+1].as_low    = p_thresh[i].as_high *
	InvZetaOS2(p_thresh[i].as_high,p_thresh[i].high_scale,p_thresh[i].high_scale,p_thresh[i].nf);
    }
  }
  for (int i=m_mzset-1;i>=0;--i) {
    double lam2               = Lambda2(i);
    p_thresh[i].as_low        = AlphaSLam(p_thresh[i].low_scale,i);
    if ((lam2>p_thresh[i].low_scale) || (p_thresh[i].as_low>1.)) ContinueAlphaS(i);
    else {
      if (i>0) {
	p_thresh[i-1].as_high = p_thresh[i].as_low *
	  ZetaOS2(p_thresh[i].as_low,p_thresh[i].low_scale,p_thresh[i].low_scale,p_thresh[i-1].nf);
      }
    }
  }

  if (mo==NULL) {
    double fac=ToType<double>(rpa->gen.Variable("RENORMALIZATION_SCALE_FACTOR"));
    if (fac!=1.0) msg_Info()<<METHOD<<"(): Setting scale factor "<<fac<<"\n";
    msg_Indent();
    p_sas = new One_Running_AlphaS((*this)(fac*m2_MZ),m2_MZ,order,thmode,aspdf,this);
  }
}


One_Running_AlphaS::~One_Running_AlphaS()
{
  if (p_thresh!=0) { delete [] p_thresh; p_thresh = NULL; }
  if (p_sas) delete p_sas;
}

double One_Running_AlphaS::Beta0(const int nf) {
  return 1./4. * (11. - (2./3.)*nf);
}

double One_Running_AlphaS::Beta1(const int nf) {
  return 1./16. * (102. - (38./3.)*nf);
}

double One_Running_AlphaS::Beta2(const int nf) {
  return 1./64. * (2857./2. - (5033./18.)*nf + (325./54.)*nf*nf);
}

double One_Running_AlphaS::Beta3(const int nf) {
  double zeta3 = 1.2020569031595942854;
  return 1./256. * ( (149753./6. + 3564.*zeta3) +
		     (-1078361./162. -6508./27.*zeta3)*nf +
		     (50065./162. +6472./81.*zeta3)*(nf*nf) +
		     (1093/729)*nf*nf*nf);
}

double One_Running_AlphaS::Lambda2(const int nr) {
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
  b[2]  = Beta2(nf)/beta0;
  b[3]  = Beta3(nf)/beta0;

  double betaL = 1./a;
  if (m_order>=1) {
    betaL     += b[1]*log(a);
    if (m_order>=2) {
      betaL   += (b[2]-b[1]*b[1])*a;
      if (m_order>=3) {
	betaL += (b[3]/2. - b[1] * b[2] + b[1]*b[1]*b[1]/2.)*a*a;
      }
    }
  }

  lambda2         = ::exp(-betaL/beta0)*mu2;
  double tas1     = AlphaSLam(mu2,nr);
  double dlambda2 = 1.e-8;
  if (dabs(tas1-as)/as>1.e-11) {
    for (;(dabs(tas1-as)/as>1.e-11);) {
      lambda2     = lambda2+dlambda2;
      double tas2 = AlphaSLam(mu2,nr);
      dlambda2    = (as-tas2)/(tas2-tas1)*dlambda2;
      tas1        = tas2;
    }
  }

  return lambda2;
}

double One_Running_AlphaS::AlphaSLam(const double Q2,const int nr)
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
  if (m_order==1) return M_PI*a;

  double log2L    = logL*logL;
  pref           *= 1./(beta0*L);
  a              += pref*(b[1]*b[1]*(log2L-logL-1.) + b[2]);
  if (m_order==2) return M_PI*a;

  // 3rd order (four loop) to be checked.
  double log3L    = logL*log2L;
  pref           *= 1./(beta0*L);
  a              += pref*(b[1]*b[1]*b[1]*(-log3L+2.5*log2L+2.*logL-0.5) 
			  - 3.*b[1]*b[2] + 0.5*b[3]);
  return M_PI*a;
}

double One_Running_AlphaS::ZetaOS2(const double as,const double mass2_os,
			       const double mu2,const int nl) {
  double zeta2g = 1.;

  // 0th order
  if (m_order==0) return zeta2g;

  // 1st order (one loop) corresponds to two loop lambda  
  double L      = log(mu2/mass2_os);
  double a      = as/M_PI;
  zeta2g       += - a*1./6.*L;
  if (m_order==1) return zeta2g; 

  // 2nd order
  double L2     = L*L;
  double a2     = a*a;
  zeta2g       += a2 *( 1./36.*L2 - 19./24.*L -7./24.);
  if (m_order==2) return zeta2g;

  // 3rd order : not yet checked ...
  double L3     = L2*L;
  double a3     = a2*a;
  double zeta2  = M_PI*M_PI/6.;
  double zeta3  = 1.2020569031595942854;
  zeta2g       += a3 * (-58933./124416. - 2./3.*zeta2*(1.+1./3.* log(2.)) 
			- 80507./27648.*zeta3 - 8521./1728.*L- 131./576. * L2 
			- 1./216.*L3 + nl*(2479./31104.+ zeta2/9. + 409./1728. * L ));
  return zeta2g;
}

double One_Running_AlphaS::InvZetaOS2(const double as,const double mass2_os,
				  const double mu2,const int nl) {
  // might be simplified considerably when using mu2==mass2
  double zeta2g  = 1.;
  // 0th order   
  if (m_order==0) return zeta2g;

  // 1st order (one loop) corresponds to two loop lambda
  double L      = log(mu2/mass2_os);
  double a      = as/M_PI;
  zeta2g       += + a*1./6.*L;
  if (m_order==1) return zeta2g;

  // 2nd order 
  double L2     = L*L;
  double a2     = a*a;
  zeta2g       += a2 *( 1./36.*L2 + 19./24.*L + 7./24.);
  if (m_order==2) return zeta2g;

  // 3rd order yet to be checked...
  double L3     = L2*L;
  double a3     = a2*a;
  double zeta2  = M_PI*M_PI/6.;
  double zeta3  = 1.2020569031595942854;
  zeta2g       += a3 * (58933./124416. + 2./3.*zeta2*(1.+1./3.* log(2.)) 
			+ 80507./27648.*zeta3 + 8941./1728.*L + 511./576. * L2 
			+ 1./216.*L3 + nl*(-2479./31104.- zeta2/9. - 409./1728. * L ));
  return zeta2g;
}



void One_Running_AlphaS::ContinueAlphaS(int & nr) {
  // shrink actual domain
  //  * to given t0        or
  //  * to alphaS=alphaCut
  double alpha_cut = 1.;
  double & beta0   = p_thresh[nr].beta0;
  double & lambda2 = p_thresh[nr].lambda2;
  double t0        = lambda2 * ::exp(M_PI/(alpha_cut*beta0));
  double as        = AlphaSLam(t0,nr);
  for (;dabs(as-alpha_cut)>1.e-8;) {
    double t1      = t0+0.00001;
    double as1     = AlphaSLam(t1,nr);
    double das     = (as -as1)/(t0-t1);
    t1             = (alpha_cut-as)/das + t0;
    t0             = t1;
    as             = AlphaSLam(t0,nr);    
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


double One_Running_AlphaS::ShowerCutQ2()
{
  if (p_sas==NULL) THROW(fatal_error,"Invalid call");
  return p_sas->m_cutq2;
}

double One_Running_AlphaS::operator[](double q2)
{
  if (p_sas==NULL) THROW(fatal_error,"Invalid call");
  return (*p_sas)(q2);
}

double One_Running_AlphaS::operator()(double q2)
{
  if (m_pdf) return p_pdf->AlphaSPDF(q2);
  double as;
  if (q2<0.) q2=-q2;
  int i = m_mzset-1;
  if (q2<=m_m2_MZ) {
    for (;!((p_thresh[i].low_scale<q2)&&(q2<=p_thresh[i].high_scale));--i) {
      if (i<=0) break;
    }
    if (p_thresh[i].nf>=0) 
      as = AlphaSLam(q2,i);
    else 
      as = q2/p_thresh[i].high_scale * p_thresh[i].as_high;
  }
  else {
    ++i;
    for (;!((p_thresh[i].low_scale<q2)&&(q2<=p_thresh[i].high_scale));++i) {
      if (i>=m_nth) break;
    }
    as   = AlphaSLam(q2,i);
  }
  return as;
}  

double  One_Running_AlphaS::AlphaS(const double q2){
  return operator()(q2);
}

int One_Running_AlphaS::Nf(const double sc)
{
  double q2(sc);
  for (int i=0;i<=m_nth;++i) {
    if (q2<=p_thresh[i].high_scale && q2>p_thresh[i].low_scale )
      return p_thresh[i].nf;
  }
  return m_nth;
}

void One_Running_AlphaS::SelfTest() {
  /*
  fastfunc ff_as;

  int np=1001;

  double smin=0.025;
  double smax=500*500;

  double mult=::exp(log(smax/smin)/double(np-2));

  ff_as.reset(np);
  double as0=operator()(0.);
  ff_as.insert(0.,as0);
  double s= smin;

  for (int i=1;i<np;++i) {
    double as =operator()(s);
    double rs =sqrt(s);
    ff_as.insert(rs,as);
    s*=mult;
  }
  ff_as.output(sqrt(smin),sqrt(smax),"alpha_s.test.dat");
  */
}


Running_AlphaS::Running_AlphaS(const double as_MZ,const double m2_MZ,
                               const int order, const int thmode,
                               const PDF::ISR_Handler_Map &isr)
{
  m_type  = std::string("Running Coupling");
  m_name  = "Alpha_QCD";
  for (ISR_Handler_Map::const_iterator it=isr.begin(); it!=isr.end(); ++it) {
    if (m_alphas.find(it->first)!=m_alphas.end())
      THROW(fatal_error, "Internal error.");
    PDF::PDF_Base *pdf(NULL);
    if (it->second->PDF(0)) pdf=it->second->PDF(0);
    if ((pdf==NULL||pdf->ASInfo().m_order<0) && it->second->PDF(1))
      pdf=it->second->PDF(1);
    m_alphas.insert(make_pair(it->first, new One_Running_AlphaS
                              (as_MZ,m2_MZ, order, thmode, pdf)));
  }
  SetActiveAs(PDF::isr::hard_process);
}

Running_AlphaS::~Running_AlphaS()
{
  for (AlphasMap::iterator it=m_alphas.begin(); it!=m_alphas.end(); ++it) {
    delete it->second;
  }
  m_alphas.clear();
}

void Running_AlphaS::SetActiveAs(PDF::isr::id id)
{
  AlphasMap::iterator it=m_alphas.find(id);
  if (it==m_alphas.end()) {
    THROW(fatal_error, "Internal Error");
  }
  else {
    p_active=it->second;
  }
}
