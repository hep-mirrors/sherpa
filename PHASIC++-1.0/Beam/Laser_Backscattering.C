#include "Laser_Backscattering.H"
#include "Run_Parameter.H"
#include "Message.H"
#include "Random.H"

using namespace BEAM;
using namespace AMATOOLS;
using namespace APHYTOOLS;
using namespace AORGTOOLS;
using namespace std;

Laser_Backscattering::Laser_Backscattering(Flavour _fl,double _pole)
{
  if ( (_fl != Flavour(kf::e)) && (_fl != Flavour(kf::e).Bar()) ) {
    msg.Error()<<"Tried to initialize Laser_Backscattering for flavour "<<_fl<<"."<<endl
	       <<"This option is not available. Terminate program."<<endl;
    abort();
  }
  type = std::string("Laser_Backscattering");


  E           = rpa.gen.Ecms()/2.;
  Ebounds[0]  = 50.;  
  Ebounds[1]  = 500.;
  pol = 0;

  // for testing purposes and to comply with 9302319
  mode        = 0;
  angles      = 0;
  Ebounds[0]  = 0.;  
  Ebounds[1]  = 500.;

  omegaL      = 1.17e-9;   // (in GeV !!!)  
  polE        = _pole;//0.85;       
  cout<<"LaserBackscattering polE :"<<polE<<endl;
  polL        = -1.;  
  //Setting pol flag if electrons are polarized
  if (polE!=0. || polL!=0.)   pol = 1;

  rho2        = 3.315865;   
  delta       = 1.387423/2.;
  nonlin1     = 0.06594662; nonlin2 = 0.7060851e-3;
  xe          = 4.*E*omegaL/sqr(_fl.PSMass());
  // to compare with 9302319

  polE        = 0.;         polL    = 0.;  
  rho2        = 0.;   
  delta       = 1.387423/2.;
  nonlin1     = 0.;         nonlin2 = 0.;
  xe          = 4.82;

  Emax   = E*xe/(1.+xe);

  xi     = nonlin1 + nonlin2 * E;
  xe    /= (1+xi);
  xmax   = xe/(1.+xe);
  xmax2  = 2.*xe/(1.+2.*xe);

  if (mode==0) upper = xmax;
          else upper = xmax2;
  peak   = xmax;
  //  AORGTOOLS::msg.Out()<<"*** mass :"<<_fl.PSMass()<<" , "<<omegaL<<" , "<<E<<" , "<<xe<<" , "<<peak<<endl;

  yfix   = 1./(1.+xe);
  yden   = log(1.+xe);
  ysteps = 50;

  totalC = 0.7115863 - 0.6776124e-3 * E + 0. * E * E; 
  total2 = totalC * 0.5540019 * (1.-exp(-37.38912 * xi * xi));
  totalE = totalC * (0.7257064 + 1.517959e-3 * E);

  totalC = total2 = totalE = 1.;


  msg.Tracking()<<"Initialised Laser-Backscattering ("<<mode<<") : "<<endl
		<<" xe,xmax = "<<xe<<", "<<xmax<<" for omegaL,mass ="
		<<omegaL<<", "<<_fl.PSMass()<<endl
		<<" with xi = "<<xi<<", norms   = "<<totalC<<"  /  "<<total2<<endl
		<<" with polE = "<<polE<<endl;
  //  PrintSpectra("spectrum.out");
}

Laser_Backscattering::~Laser_Backscattering() {}

void Laser_Backscattering::PrintSpectra(std::string filename) {
  bool flag = 0;
  ofstream ofile;
  if (filename != string("")) {
    ofile.open(filename.c_str());
    flag = 1;
  } 

  double z,res1,res2,res3,restot;
  double deg;
  for (int i=1;i<1510;i++) { 
    z   = xmax2*i*.0007;
    restot = deg  = 0.;
    restot += res1 = Compton(z,polE,polL,deg);
    restot += res2 = TwoPhotons(z,polE,polL,deg); 
    restot += res3 = Rescattering(z,polE,polL,deg);
    if (flag) ofile<<" "<<z<<"  "<<res1<<"  "<<res1+res2<<"  "<<res1+res2+res3;
    else  msg.Out()<<" "<<z<<"  "<<res1<<"  "<<res1+res2<<"  "<<res1+res2+res3;
    if (IsZero(restot)) {deg = 0.;restot = 1.e-17;}
    if (flag) ofile<<"  "<<deg/restot<<endl;
    else  msg.Out()<<"  "<<deg/restot<<endl;
  }

  if (flag) ofile.close();
}

bool Laser_Backscattering::CalculateWeight(double x,double scale) 
{
  if (!( (x*E >= Ebounds[0]) && (x*E <= Ebounds[1]))) {
      //cout<<"LSB : "<<x<<" : "<<Ebounds[0]<<" "<<x*E<<" "<<Ebounds[1]<<endl;
    weight = 0.;
    return 0;
  }

  polar = 0.;
  double spec;
  switch (mode) {
  case 1: 
    spec = Compton(x,polE,polL,polar);
    break;
  case 2: 
    spec = TwoPhotons(x,polE,polL,polar);
    break;
  case 3: 
    spec = Rescattering(x,polE,polL,polar);  
    break;
  default:
    spec = Compton(x,polE,polL,polar) + 
           TwoPhotons(x,polE,polL,polar) + 
           Rescattering(x,polE,polL,polar);  
    break;
  }
  polar  = polar/spec;
  weight = spec;
  //cout<<"weight(spec) = "<<spec<<endl;
   cout<<" ======================================== "<<endl;
   cout<<" x,scale ="<<x<<","<<scale<<endl;
   cout<<" polar= "<<polar<<endl;
   cout<<" weight="<<weight<<endl;

  return 1;
};

double Laser_Backscattering::Weight(Flavour flin)
{
  if (weight<=0.) return 0.;
  if (flin != Flavour(kf::photon)) return 0.;
  return weight;
}


double Laser_Backscattering::Compton(double x,double pole,double poll,double & deg)
{
  if ((x<0.) || (x>xmax) || (totalC < 0.) ) {
      return 0.;
  }

  double value  = SimpleCompton(x,xe,pole*poll);

  double g2    = xe/x - xe - 1;
  if (g2<0.) {
    if (pol) deg += value * Polarization(x,xe,pole,poll);
    return value;
  }

  double damp   = exp(-rho2 * g2/8.);
  if (pol) deg += damp * value * totalC * Polarization(x,xe,pole,poll);

  double wt = damp * totalC * value;
  return wt;
}

double Laser_Backscattering::TwoPhotons(double x,double pole,double poll,double & deg)
{
  if ((x<0.) || (x>xmax2) || (total2 < 0.)) return 0.;

  double value  = SimpleCompton(x,2.*xe,pole*poll);

  double g2    = 2.*xe/x - 2.*xe - 1;
  if (g2<0.) {
    if (pol) deg += value * total2 * Polarization(x,2.*xe,pole,poll);
    return value;
  }

  double damp   = exp(-rho2 * g2/8.) * pow(g2,delta);
  if (pol) deg += damp * value * total2 * Polarization(x,2.*xe,pole,poll);

  return damp * total2 * value;
}

double Laser_Backscattering::Rescattering(double x,double pole,double poll,double & deg)
{
  if ((x<0.) || (x>xmax) || (totalE < 0.)) return 0.;

  double yMin  = Max(yfix,0.5 * x * (1.+sqrt(4./(x*xe) + 1.)));
  if (yMin > 1.) return 0.;
  
  double y1, y2;
  double dy, dp, value, pvalue;
  double val1,val2,p1,p2;

  value    = pvalue  = 0.;
  y1       = y2      = yMin;
  y1      *= 1.000001;
  dy       = (1.-yMin)/ysteps;

  val1     = log(1.+y1*xe)/(y1 * yden) *
    SimpleCompton(x/y1,y1*xe,0.)*SimpleCompton(1-y1,xe,pole*poll);
  p1       = Polarization(x/y1,y1*xe,pole,poll);

  for (int i=0;i<ysteps;i++) {
    y2       += dy;
    val2      = log(1.+y2*xe)/(y2 * yden) *
      SimpleCompton(x/y2,y2*xe,0.)*SimpleCompton(1-y2,xe,pole*poll);
    value    += 0.5*(val1+val2)*dy;
    if (pol) {
      p2      = Polarization(x/y2,y2*xe,pole,poll);
      pvalue += 0.5*(val1*p1+val2*p2)*dy;
      p1      = p2;
    }
    val1    = val2;
  }

  if (pol) deg += pvalue*value*totalE;
  return totalE * value;
}



double Laser_Backscattering::SimpleCompton(double x,double z,double pol2) 
{
  double max   = z/(1.+z);
  //cout<<"SimpleCompton : "<<x<<", "<<z<<" : "<<max<<endl;
  if ((x<0.) || (x>max)) return 0.;

  double help  = x/(z*(1.-x));
  double value = 1.-x  + 1./(1.-x) - 4.*help + 4.*help*help; 
  value       -= pol2 * x*(2.-x)/(1.-x) * (2*help - 1.); 

  double norm  = (z*z*z+18.*z*z+32.*z+16.)/(2.*z*(z+1.)*(z+1.));
  norm        += (1.-4./z-8./(z*z)) * log(1.+z);
  norm        -= pol2 * (2. + z*z/((z+1.)*(z+1.)) - (1.+2./z) * log(z+1.));

  //cout<<"SimpleCompton : "<<x<<", "<<z/(1.+z)<<" : "<<value<<" : "<<norm<<endl;   
    return value/norm;
}

double Laser_Backscattering::Polarization(double x,double z,double pole,double poll)
{
  //cout<<" In Polarization."<<endl;
  double max   = z/(1.+z);
  if ((x<0.) || (x>max)) return 0.;

  double help1 = x/(z*(1.-x));
  double help2 = 1. - x + 1./(1.-x);
  double value = pole * help1 * z * (1.+(1.-x)*sqr(2.*help1-1.)) -
                 poll * help2 * (2.*help1-1.);
  double norm  = help2 + 4.*help1*(help1-1.) - pole*poll*help1*z*(2.-x)*(2.*help1-1.); 

  return value/norm;
}














