#include "Running_AlphaS.H"
#include "Sort.H"
#include "Run_Parameter.H"
#include "Message.H"
#include "MathTools.H"
// #include "Model.H"
// #include "Distribution_Functions.H"
#include <iostream>

namespace APHYTOOLS {
  Running_AlphaS * as =0;
}

using namespace APHYTOOLS;
using namespace AORGTOOLS;
using namespace AMATOOLS;

using namespace std;



double Running_AlphaS::Beta0(int nf) {
  return 1./4. * (11. - (2./3.)*nf);
}

double Running_AlphaS::Beta1(int nf){
  return 1./16. * (102. - (38./3.)*nf);
}

double Running_AlphaS::Beta2(int nf){
  return 1./64. * (2857./2. - (5033./18.)*nf + (325./54.)*nf*nf);
}

double Running_AlphaS::Beta3(int nf){
  double zeta3 = 1.2020569031595942854;
  return 1./256. * ( (149753./6. + 3564.*zeta3) +
		     (-1078361./162. -6508./27.*zeta3)*nf +
		     (50065./162. +6472./81.*zeta3)*(nf*nf) +
		     (1093/729)*nf*nf*nf);
}


// only called during initialisation (special order necessary!)
double Running_AlphaS::Lambda2(int nr) {
  double as  = thresh[nr].as_low; 
  double mu2 = thresh[nr].low_scale;
  if (as==0.) {
    as  = thresh[nr].as_high;
    mu2 = thresh[nr].high_scale;
  }

  const double a   = as/M_PI;

  // fill selected

  // using shorter names
  int &    nf     = thresh[nr].nf;
  double & beta0  = thresh[nr].beta0;
  double *b = thresh[nr].b;
  double & lambda2= thresh[nr].lambda2;
  
  // calculate beta coefficients
  beta0 = Beta0(nf);
  b[1]  = Beta1(nf)/beta0;
  b[2]  = Beta2(nf)/beta0;
  b[3]  = Beta3(nf)/beta0;

  double betaL;
  betaL = 1./a;
  if (order>=1) {
    betaL+= b[1]*log(a);
    if (order>=2) {
      betaL+= (b[2]-b[1]*b[1])*a;
      if (order>=3) {
	betaL+= (b[3]/2. - b[1] * b[2] + b[1]*b[1]*b[1]/2.)*a*a;
      }
    }
  }
  lambda2 = ::exp(-betaL/beta0)*mu2;

  // iterative improve lambda2
  double  tas1=AlphaSLam(mu2,nr);
  double dlambda2=1.e-8;
  if (dabs(tas1-as)/as>1.e-11) {
    for (;(dabs(tas1-as)/as>1.e-11);) {
      lambda2=lambda2+dlambda2;
      double tas2=AlphaSLam(mu2,nr);
      dlambda2 = (as-tas2)/(tas2-tas1)*dlambda2;
      tas1=tas2;
    }
  }

  return lambda2;
};

double Running_AlphaS::AlphaSLam(double Q2, int nr)
{
  // using shorter names
  double & beta0  = thresh[nr].beta0;
  double *  b     = thresh[nr].b;
  double & lambda2= thresh[nr].lambda2;

  double L    = log(Q2/lambda2);
  double pref = 1./(beta0*L);

  // 0th order (one loop)   Check: OK!
  double a;
  a  = pref;
  if (order==0) return M_PI*a;

  // 1st order (two loop)   Check: OK!
  double logL=log(L);
  pref*=1./(beta0*L);
  a += -pref*(b[1] * logL);
  if (order==1) return M_PI*a;

  // 2nd order (three loop)  Check: OK!
  double log2L=logL*logL;
  pref*=1./(beta0*L);
  a += pref*(b[1]*b[1]*(log2L-logL-1.) + b[2]);
  if (order==2) return M_PI*a;

  // 3rd order (four loop)
  double log3L=logL*log2L;
  pref*=1./(beta0*L);
  a += pref*(b[1]*b[1]*b[1]*(-log3L+2.5*log2L+2.*logL-0.5) 
	     - 3.*b[1]*b[2] + 0.5*b[3]);
  return M_PI*a;
}

// down:   alpha(nf-1) = zeta^2 * apha(nf)
double Running_AlphaS::ZetaOS2(double as, double mass2_os, double mu2, int nl) {
  // might be simplified considerably when using mu2==mass2
  double zeta2g;

  // 0th order        Check: OK! 
  zeta2g = 1.;
  if (order==0) return zeta2g;

  // 1st order (one loop) corresponds to two loop lambda  
  //                  Check: OK!
  double L=log(mu2/mass2_os);
  double a= as/M_PI;
  zeta2g+= - a*1./6.*L;
  if (order==1) {
    return zeta2g; 
  }
  // 2nd order        Check: OK!
  double L2=L*L;
  double a2=a*a;
  zeta2g+= a2 *( 1./36.*L2 - 19./24.*L -7./24.);
  if (order==2) return zeta2g;

  // 3rd order
  double L3=L2*L;
  double a3=a2*a;
  double zeta2=M_PI*M_PI/6.;
  double zeta3=1.2020569031595942854;
  zeta2g+= a3 * (-58933./124416. - 2./3.*zeta2*(1.+1./3.* log(2.)) 
              - 80507./27648.*zeta3 - 8521./1728.*L- 131./576. * L2 
              - 1./216.*L3 + nl*(2479./31104.+ zeta2/9. + 409./1728. * L ));
  return zeta2g;
}

double Running_AlphaS::InvZetaOS2(double as, double mass2_os, double mu2, int nl) {
  // might be simplified considerably when using mu2==mass2
  double zeta2g;

  // 0th order        Check: OK!
  zeta2g = 1.;
  if (order==0) return zeta2g;

  // 1st order (one loop) corresponds to two loop lambda
  //                  Check: OK!
  double L=log(mu2/mass2_os);
  double a= as/M_PI;
  zeta2g+= + a*1./6.*L;
  if (order==1) {
    return zeta2g;
  }

  // 2nd order        Check: OK!
  double L2=L*L;
  double a2=a*a;
  zeta2g+= a2 *( 1./36.*L2 + 19./24.*L + 7./24.);
  if (order==2) return zeta2g;

  // 3rd order
  double L3=L2*L;
  double a3=a2*a;
  double zeta2=M_PI*M_PI/6.;
  double zeta3=1.2020569031595942854;
  zeta2g+= a3 * (58933./124416. + 2./3.*zeta2*(1.+1./3.* log(2.)) 
              + 80507./27648.*zeta3 + 8941./1728.*L + 511./576. * L2 
              + 1./216.*L3 + nl*(-2479./31104.- zeta2/9. - 409./1728. * L ));
  return zeta2g;
}

Running_AlphaS::Running_AlphaS(){
  msg.Tracking()<<"Initialising Running_AlphaS"<<endl;
  //
  thresh  = 0;

  // further, in the moment the thresholds are assumed to coincide with 
  // the on-shell masses of the quark
  Init();
};


Running_AlphaS::~Running_AlphaS()
{
  if (thresh != 0) delete[] thresh;
}

void Running_AlphaS::Init()
{
  msg.Tracking()<<"in Running_AlphaS::Init"<<endl;

  CF    = 4./3.;        
  CA    = 3.;           

  if (rpa.consts.IsASRunning()==Switch::Off)
    mode =0;
  else
    mode =1;

  // read in order!!!!!!
  // working is: 0 == one loop,  1 == two loops (recom.), and  2 == three loops 
  //   attention! not working yet is:     3 == four loops
  order=0;

  Data_Read dr(rpa.GetPath()+std::string("/")+rpa.me.ModelFile());
  // Note: all alphas become inverted (see below) !!!
  as_MZ      = dr.GetValue<double>("alpha_S(MZ)");
  as_eff     = rpa.consts.FixedAlphaS();
  if (as_eff == NotDefined<double>()) {
    as_eff    = dr.GetValue<double>("alpha_S_fixed");
    rpa.consts.SetFixedAlphaS(as_eff);
  }
  if ( (mode==0) && (as_eff==NotDefined<double>())) {
    if (as_MZ==NotDefined<double>()) {
      msg.Error()<<" ERROR: alpha_S_fixed  has to be defined in "<< rpa.me.ModelFile()<<" ! "<<endl;
      msg.Error()<<"        if ""running coupling"" is not used"<<endl;
      abort();
    } 
    else {
      msg.Out()<<" WARNING: alpha_S_fixed  not defined in "<< rpa.me.ModelFile()<<endl;
      msg.Out()<<"          using alpha_S(MZ) instead"<<endl;
      as_eff=as_MZ;
      // or set mode = 2
      // and calculate as_eff after initializiation
    }
  }

  if ( (mode==1)  && (as_MZ==NotDefined<double>())) {
      msg.Error()<<" ERROR: alpha_S(MZ)  has to be defined in "<< rpa.me.ModelFile()<<" ! "<<endl;
      msg.Error()<<"        if ""running coupling"" is used"<<endl;
      abort();
  }
  if (as_MZ==NotDefined<double>()) {
      msg.Error()<<" WARNING: alpha_S(MZ)  has to be defined in "<< rpa.me.ModelFile()<<" ! "<<endl;
      // as_MZ = 0.118;
      as_MZ = 0.1172;
      // PDG 2002 value is : 0.1172(20)
  }


  m2_MZ = sqr(Flavour(kf::Z).PSMass());


  //------------------------------------------------------------
  // SM thresholds for strong interactions, i.e. QCD
  //------------------------------------------------------------
  short int count = 0;  // count strong flavours
  Fl_Iter fli;
  for (Flavour flav=fli.first();flav!=Flavour(kf::none);flav = fli.next()) 
    if (flav.Strong() && flav.IsOn() && flav.Kfcode()<30)  count++;

  nth = count;

  thresh = new AsDataSet[nth+1]; 
  double*  masses = new double[nth];
  
  // fill heavy quark thresholds
  count = 0;
  for (Flavour flav=fli.first();flav!=Flavour(kf::none);flav = fli.next()) {
    if (flav.Strong() && flav.IsOn() && flav.Kfcode()<30) {
      masses[count] = sqr(flav.PSMass());  // on shell masses
      count++;
    }
  }

  // sort thresholds with sort-tool 
  Bubble_Sort::Down<double>(masses,nth);

  // fill container (assuming one gluon and nth-1 quarks)
  int j=0; 
  mzset=0;
  for (int i=0; i<nth; ++j) {
    if ((masses[i]>m2_MZ)&&(!mzset)) {
      //insert Z boson (starting point for any evaluation) 
      mzset=j;
      thresh[j].low_scale=m2_MZ;
      thresh[j].as_low=as_MZ;
      thresh[j].nf=i-1;
    } 
    else {
      thresh[j].low_scale=masses[i];
      thresh[j].as_low=0.;
      thresh[j].nf=i;
      ++i;
    }


    if (j>0) {
      thresh[j-1].high_scale=thresh[j].low_scale;
      thresh[j-1].as_high=thresh[j].as_low;
    }
  }

  // set MZ point, in case all other masses are below
  if (!mzset) {
    int j=nth;
    mzset=j;
    thresh[j].low_scale=m2_MZ;
    thresh[j].as_low=as_MZ;
    thresh[j-1].high_scale=m2_MZ;
    thresh[j-1].as_high=as_MZ;
    thresh[j].nf=thresh[j-1].nf;      
  }
  thresh[nth].high_scale=1.e7;  // arbitrary upper limit
  thresh[nth].as_high=0.;

  // delete temp mass vector
  delete [] masses;

  // evolve up from as_MZ
  for (int i=mzset; i<=nth; ++i) {
    // calculate beta coefficients and lambda2
    Lambda2(i);

    // calculate alphaS at upper limit
    thresh[i].as_high=AlphaSLam(thresh[i].high_scale,i);

    // check matching condition
    if (i<nth) {
      // check arguments!!! especially "nf"
      thresh[i+1].as_low =thresh[i].as_high*
	InvZetaOS2(thresh[i].as_high,thresh[i].high_scale,thresh[i].high_scale,thresh[i].nf);
    }
  }

  // evolve down from as_MZ
  for (int i=mzset-1; i>=0; --i) {
    // calculate beta coefficients and lambda2
    double lam2= Lambda2(i);

    // calculate alphaS at lower limit
    thresh[i].as_low=AlphaSLam(thresh[i].low_scale,i);

    if ((lam2>thresh[i].low_scale)||(thresh[i].as_low>1.)) {
      // continuation of alpha_s over landau pole
      ContinueAlphaS(i);
    }
    else {
      // check matching condition
      if (i>0) {
      // check arguments!!! especially "nf"
	thresh[i-1].as_high =thresh[i].as_low*
	  ZetaOS2(thresh[i].as_low,thresh[i].low_scale,thresh[i].low_scale,thresh[i-1].nf);
      }
    }
  }

  // --- print status so far ---
  if (rpa.gen.Tracking()) {
    msg.Tracking()<<" Init (3) "<<endl;
    for (int i=0; i<=nth; ++i) {
      msg.Tracking()<<"  s_low ="<<thresh[i].low_scale
		    <<"  as_low="<<thresh[i].as_low
		    <<"  nf="<<thresh[i].nf<<endl;
      msg.Tracking()<<"  s_high ="<<thresh[i].high_scale
		    <<"  as_high="<<thresh[i].as_high
		    <<"  lambda2="<<thresh[i].lambda2<<endl;
    
    }
    msg.Tracking()<<endl;
  }

  msg.Tracking()<<" Initialisation of Alpha_S completed."<<endl;
}


void Running_AlphaS::ContinueAlphaS(int & nr){
  // shrink actual domain
  //  * to given t0        or
  //  * to alphaS=alphaCut
  double alpha_cut=1.0*M_PI;   // make parameter

  double & beta0   = thresh[nr].beta0;
  double & lambda2 = thresh[nr].lambda2;
  double t0,as;
  t0=lambda2 * ::exp(M_PI/(alpha_cut*beta0));
  as=AlphaSLam(t0,nr);
  for (;dabs(as-alpha_cut)>1.e-8;) {
    double t1=t0+0.00001;
    double as1=AlphaSLam(t1,nr);
    double das = (as -as1)/(t0-t1);
    t1 = (alpha_cut-as)/das + t0;
//     if (t1<1.e-6) {
//       t1= t0/2.;
//     }
    t0=t1;
    as=AlphaSLam(t0,nr);    
  }

  // modify lower domains
  thresh[nr].low_scale=t0;
  thresh[nr-1].high_scale=t0;
  thresh[nr].as_low=  as;
  thresh[nr-1].as_high=  as;

  for (int i = nr-1; i>=0; --i) {
    thresh[i].nf=-1;  // i.e. no ordinary running !!!
    thresh[i].lambda2=0.;
    thresh[i].as_low=thresh[i].as_high/thresh[i].high_scale*thresh[i].low_scale;
    if (i>0) thresh[i-1].as_high=thresh[i].as_low;
  }

  // return value:  
  nr =0;
  // finish initialisation.
}



double Running_AlphaS::operator()(double t)
{
  //  Check: OK! (up to 2nd order (3 loops))

  double as;
  if (t<0.) t=-t;
  int i=mzset-1;
  if (t<=m2_MZ) {
    // evolve down
    for (;!((thresh[i].low_scale<t)&&(t<=thresh[i].high_scale));--i) {
      if (i<=0) {
	break;
      }
    }
    if (thresh[i].nf>=0)
      as=AlphaSLam(t,i);
    else
      // linear (only used in conjuction with ContinueAlphaS())
      as=t/thresh[i].high_scale*thresh[i].as_high;
  }
  else {
    ++i;
    // evolve up
    for (;!((thresh[i].low_scale<t)&&(t<=thresh[i].high_scale));++i) {
      if (i>=nth) {
	break;
      }
    }
    as=AlphaSLam(t,i);
  }
  return as;
}  



void Running_AlphaS::SelfTest(){
  /*
  fastfunc ff_as;

  int np=1001;

  double smin=0.025;
  double smax=500*500;

  double mult=::exp(log(smax/smin)/double(np-2));

  ff_as.reset(np);

  // fill data
  msg.Tracking()<<endl;
  msg.Tracking()<<"============================================================"<<endl;
  msg.Tracking()<<" check selfconsistence at thresholds:"<<endl;
  for (int j=0;j<=nth;++j) {
    msg.Tracking()<<" scale="<<thresh[j].low_scale;
    msg.Tracking()<<"   ist="<<AlphaSLam(thresh[j].low_scale,j);
    msg.Tracking()<<"   soll="<<thresh[j].as_low<<endl;
    msg.Tracking()<<" scale="<<thresh[j].high_scale;
    msg.Tracking()<<"   ist="<<AlphaSLam(thresh[j].high_scale,j);
    msg.Tracking()<<"   soll="<<thresh[j].as_high<<endl;
  }

  msg.Tracking()<<" create data file .... "<<endl;
  double as0=operator()(0.);
  ff_as.insert(0.,as0);
  double s= smin;

  for (int i=1;i<np;++i) {
    double as =operator()(s);
    double rs =sqrt(s);
    ff_as.insert(rs,as);
    //    msg.Tracking()<<" rs="<<rs<<"  as="<<as<<endl;
    s*=mult;
  }

  // write file
  ff_as.output(sqrt(smin),sqrt(smax),"alpha_s.test.dat");
  */
}

double  Running_AlphaS::AlphaS(double t){
  if (mode) {
    // running coupling
    return operator()(t);
  }
  // fixed alpha_S!
  return as_eff;
}
