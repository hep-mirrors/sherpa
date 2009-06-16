#include "BEAM/Main/Laser_Backscattering.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Random.H"

#include <iomanip>

using namespace BEAM;
using namespace ATOOLS;
using namespace std;

Laser_Backscattering::Laser_Backscattering(const ATOOLS::Flavour _beam,
					   const double _energy,const double _polarisation,
					   const double _energyL,const double _polarisationL,
					   const int _mode,const int _angles,
					   const int _nonlin,const int _dir) :
  Beam_Base(string("Laser_Backscattering"),_beam,_energy,_polarisation,_dir),
  m_energyL(_energyL), m_polarisationL(_polarisationL), m_mode(_mode), m_angles(_angles)
{
  m_bunch        = Flavour(kf_photon);
  double disc    = 1.-sqr(m_bunch.Mass()/m_energy);
  m_vecout       = Vec4D(m_energy,0.,0.,_dir*m_energy*sqrt(disc));
  m_Ebounds[0]   = 0.;  
  m_Ebounds[1]   = 5.e10;

  if (m_energy>500. && m_mode!=-1 ) {
    msg_Out()<<" WARNING: The CompAZ spectrum is only valid for electron energies "<<endl
	     <<"          between 100 GeV and 400 GeV! "<<endl;

  } 

  if (m_angles!=0) {
    msg_Out()<<"WARNING:  Laser_Backscattering::Laser_Backscattering."<<endl
	     <<"   Angular distribution not implemented yet. Assume collinear beam."<<endl; 
    m_angles     = 0;
  }
  //if (m_angles==0) m_lab = Vec4D(m_energy,0.,0.,_dir*m_energy);


  // Setting m_pol flag if electrons ore Laser are polarized
  m_pol          = 0;
  if (m_polarisation!=0. || m_polarisationL!=0.) m_pol = 1;

  // Nonlinear corrections.
  m_rho2   = 3.315865;   
  m_delta  = 1.387423/2.;
  if (_nonlin==1 && m_mode!=-1) { m_nonlin1 = 0.06594662; m_nonlin2 = 0.7060851e-3; }
  else { m_nonlin1 = 0.;         m_nonlin2 = 0.;           }
  m_xe     = 4.*m_energy*m_energyL/sqr(ATOOLS::Flavour(kf_e).Mass(true));
  m_xi     = m_nonlin1 + m_nonlin2 * m_energy;
  m_xe    /= (1+m_xi);
  m_xmax   = m_xe/(1.+m_xe);
  m_xmax2  = 2.*m_xe/(1.+2.*m_xe);

  if (m_mode==1 || m_mode==-1) m_upper = m_xmax;
  else m_upper = m_xmax2;
  m_peak   = m_xmax;

  m_yfix   = 1./(1.+m_xe);
  m_yden   = log(1.+m_xe);
  m_ysteps = 100;

  if (m_mode==-1) {
    m_totalC=1.;
    m_total2=0.;
    m_totalE=0.;
  }
  else {
    m_totalC = 0.7115863 - 0.6776124e-3 * m_energy + 0. * m_energy * m_energy; 
    m_total2 = m_totalC * 0.5540019 * (1.-exp(-37.38912 * m_xi * m_xi));
    m_totalE = m_totalC * (0.7257064 + 1.517959e-3 * m_energy);
  }
}



Beam_Base * Laser_Backscattering::Copy() {
  if (m_nonlin1>0.) return new Laser_Backscattering(m_beam,m_energy,m_polarisation,
						    m_energyL,m_polarisationL,m_mode,m_angles,1,m_dir);
  return new Laser_Backscattering(m_beam,m_energy,m_polarisation,
				  m_energyL,m_polarisationL,m_mode,m_angles,0,m_dir);
}

Laser_Backscattering::~Laser_Backscattering() {}

void Laser_Backscattering::PrintSpectra(std::string filename,int mode) {
  
  if (mode==0) {

    bool flag = 0;
    ofstream ofile;
    if (filename != string("")) {
      ofile.open(filename.c_str());
      flag = 1;
    } 

    double z,res1,res2,res3,restot;
    double deg;
    for (int i=1;i<1510;i++) { 
      z   = m_xmax2*i*.0007;
      restot = deg  = 0.;
      restot += res1 = Compton(z,m_polarisation,m_polarisationL,deg);
      restot += res2 = TwoPhotons(z,m_polarisation,m_polarisationL,deg); 
      restot += res3 = Rescattering(z,m_polarisation,m_polarisationL,deg);
      if (flag) ofile<<" "<<z<<"  "<<res1<<"  "<<res1+res2<<"  "<<res1+res2+res3;
      else  msg_Out()<<" "<<z<<"  "<<res1<<"  "<<res1+res2<<"  "<<res1+res2+res3;
      if (IsZero(restot)) {deg = 0.;restot = 1.e-17;}
      if (flag) ofile<<"  "<<deg/restot<<endl;
      else  msg_Out()<<"  "<<deg/restot<<endl;
    }

    if (flag) ofile.close();
  }
  if (mode==1) {
    ofstream  bsp(filename.c_str());
    bsp.setf(ios::scientific);
    bsp.setf(ios::right,ios::adjustfield);
    double xmin=0;
    double xmax=1;
    int steps=50;
    for (int i=0;i<=steps;++i) {
      double z=xmin + double(i)/double(steps) *(xmax-xmin);
    
      double p1=0,p2=0,p3=0,pt=0;
      double r1,r2,r3,rt;
      r1 = Compton(z,m_polarisation,m_polarisationL,p1);
      r2 = TwoPhotons(z,m_polarisation,m_polarisationL,p2);
      r3 = Rescattering(z,m_polarisation,m_polarisationL,p3);  
      rt = r1 + r2 + r3;
      pt = p1 + p2 + p3;
      if (r1==0.) p1=0.; else p1=p1/r1;
      if (r2==0.) p2=0.; else p2=p2/r2;
      if (r3==0.) p3=0.; else p3=p3/r3;
      if (rt==0.) pt=0.; else pt=pt/rt;
    
      bsp<<z<<"  "<<r1<<"  "<<p1<<"  "<<r2<<"  "<<p2<<"  "
	 <<r3<<"  "<<p3<<"  "<<rt<<"  "<<pt<<endl;
    }
  }
}

bool Laser_Backscattering::CalculateWeight(double _x,double _scale) 
{
  m_x = _x; m_Q2 = _scale;
  
  m_polar = 0.;
  double spec;
  switch (m_mode) {
  case -1: 
  case 1: 
    spec = Compton(_x,m_polarisation,m_polarisationL,m_polar);
    break;
  case 2: 
    spec = TwoPhotons(_x,m_polarisation,m_polarisationL,m_polar);
    break;
  case 3: 
    spec = Rescattering(_x,m_polarisation,m_polarisationL,m_polar);  
    break;
  default:
    {
      spec = Compton(_x,m_polarisation,m_polarisationL,m_polar) + 
	TwoPhotons(_x,m_polarisation,m_polarisationL,m_polar) + 
	Rescattering(_x,m_polarisation,m_polarisationL,m_polar);  
      break;
    }
  }
  m_polar  = m_polar/spec;
  m_weight = spec;

  return 1;
}

double Laser_Backscattering::Weight(Flavour flin)
{
  if (m_weight<=0.) return 0.;
  //if (flin != Flavour(kf_photon)) return 0.;
  return m_weight;
}

ATOOLS::Vec4D Laser_Backscattering::OutMomentum() {
  if (m_angles==0) return m_x*m_vecout;
  msg_Error()<<"Error in Laser_Backscattering::OutMomentum()."<<endl
		     <<"    m_angles != 0 not implemented yet."<<endl;
  return m_x*m_vecout; 
}

double Laser_Backscattering::Compton(double x,double pole,double poll,double & deg)
{
  if ((x<=0.) || (x>m_xmax) || (m_totalC < 0.) ) return 0.;

  double value  = SimpleCompton(x,m_xe,pole*poll);

  double g2    = m_xe/x - m_xe - 1;

  if (g2<0. || m_mode==-1) {
    if (m_pol) deg += m_totalC * value * Polarisation(x,m_xe,pole,poll);
    return m_totalC * value;
  }

  double damp   = exp(-m_rho2 * g2/8.);
  if (m_pol) deg += damp * value * m_totalC * Polarisation(x,m_xe,pole,poll);

  double wt = damp * m_totalC * value;
  return wt;
}

double Laser_Backscattering::TwoPhotons(double x,double pole,double poll,double & deg)
{
  if ((x<=0.) || (x>m_xmax2) || (m_total2 < 0.) || m_mode==-1) return 0.;

  double value  = SimpleCompton(x,2.*m_xe,pole*poll);

  double g2    = 2.*m_xe/x - 2.*m_xe - 1;
  if (g2<0.) {
    if (m_pol) deg += value * m_total2 * Polarisation(x,2.*m_xe,pole,poll);
    return value;
  }

  double damp   = exp(-m_rho2 * g2/8.) * pow(g2,m_delta);
  if (m_pol) deg += damp * value * m_total2 * Polarisation(x,2.*m_xe,pole,poll);

  return damp * m_total2 * value;
}

double Laser_Backscattering::Rescattering(double x,double pole,double poll,double & deg)
{
  if ((x<=0.) || (x>m_xmax) || (m_totalE < 0.)|| m_mode==-1) return 0.;

  double yMin  = Max(m_yfix,0.5 * x * (1.+sqrt(4./(x*m_xe) + 1.)));
  if (yMin > 1.) return 0.;
  
  double y1, y2;
  double dy, value, pvalue;
  double val1,val2,p1,p2;

  value    = pvalue  = 0.;
  y1       = y2      = yMin;
  y1      *= 1.000001;
  dy       = (1.-yMin)/double(m_ysteps);

  val1     = log(1.+y1*m_xe)/(y1 * m_yden) *
    SimpleCompton(x/y1,y1*m_xe,0.)*SimpleCompton(1-y1,m_xe,pole*poll);
  p1       = Polarisation(x/y1,y1*m_xe,0.,poll);
  
  for (int i=0;i<m_ysteps;i++) {
    y2       += dy;
    val2      = log(1.+y2*m_xe)/(y2 * m_yden) *
      SimpleCompton(x/y2,y2*m_xe,0.)*SimpleCompton(1-y2,m_xe,pole*poll);
    value    += 0.5*(val1+val2)*dy;
    if (m_pol) {
      p2      = Polarisation(x/y2,y2*m_xe,0.,poll);
      pvalue += 0.5*(val1*p1+val2*p2)*dy;
      p1      = p2;
    }
    val1    = val2;
  }
  if (m_pol) deg += pvalue*m_totalE;
  return m_totalE * value;
}



double Laser_Backscattering::SimpleCompton(double x,double z,double pol2) 
{
  double max   = z/(1.+z);
  if ((x<0.) || (x>max)) return 0.;

  double help  = x/(z*(1.-x));
  double value = 1.-x  + 1./(1.-x) - 4.*help + 4.*help*help; 
  value       -= pol2 * x*(2.-x)/(1.-x) * (2*help - 1.);                   

  double norm  = (z*z*z+18.*z*z+32.*z+16.)/(2.*z*(z+1.)*(z+1.));           
  norm        += (1.-4./z-8./(z*z)) * log(1.+z);                           
  norm        -= pol2 * (2. + z*z/(2.*(z+1.)*(z+1.)) - (1.+2./z) * log(z+1.));

  return value/norm;
}

double Laser_Backscattering::Polarisation(double x,double z,double pole,double poll)
{
  double max   = z/(1.+z);
  if ((x<0.) || (x>max)) return 0.;

  double help1 = x/(z*(1.-x));
  double help2 = 1. - x + 1./(1.-x);
 
  double value = pole * help1 * z * (1.+(1.-x)*sqr(2.*help1-1.)) -
    poll * help2 * (2.*help1-1.);
  double norm  = help2 + 4.*help1*(help1-1.) - pole*poll*help1*z*(2.-x)*(2.*help1-1.); 
 
  return value/norm;
}


bool Laser_Backscattering::PolarisationOn() 
{
  if (m_polarisationL!=0. || m_polarisation!=0.) return 1;
  return 0;
}











