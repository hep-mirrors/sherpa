#include "Laser_Backscattering.H"
#include "Run_Parameter.H"
#include "Message.H"
#include "Random.H"

using namespace BEAM;
using namespace AMATOOLS;
using namespace APHYTOOLS;
using namespace AORGTOOLS;
using namespace std;

Laser_Backscattering::Laser_Backscattering(const APHYTOOLS::Flavour _beam,
					   const double _energy,const double _polarization,
					   const double _energyL,const double _polarizationL,
					   const int _mode,const int _angles,
					   const int _nonlin,const int _dir) :
  Beam_Base(string("Laser_Backscattering"),_beam,_energy,_polarization,_dir),
  m_energyL(_energyL), m_polarizationL(_polarizationL), m_mode(_mode), m_angles(_angles)
{
  m_Ebounds[0]   = 0.;  
  m_Ebounds[1]   = 500.;

  if (m_angles!=0) {
    msg.Error()<<"Warning in Laser_Backscattering::Laser_Backscattering."<<endl
	       <<"   Angular distribution not implemented yet. Assume collinear beam."<<endl; 
    m_angles     = 0;
  }
  if (m_angles==0) m_lab = m_vecout = Vec4D(m_energy,0.,0.,_dir*m_energy);


  if (m_energy>m_Ebounds[1] || m_energy<m_Ebounds[0]) {
    msg.Error()<<"Warning in Laser_Backscattering::Laser_Backscattering."<<endl
	       <<"   m_energy = "<<m_energy<<" out of bounds ... . Continue."<<endl;
  }
    
  // Setting m_pol flag if electrons ore Laser are polarized
  m_pol          = 0;
  if (m_polarization!=0. || m_polarizationL!=0.) m_pol = 1;

  // Nonlinear corrections.
  m_rho2   = 3.315865;   
  m_delta  = 1.387423/2.;
  if (_nonlin==1) { m_nonlin1 = 0.06594662; m_nonlin2 = 0.7060851e-3; }
             else { m_nonlin1 = 0.;         m_nonlin2 = 0.;           }
  m_xe     = 4.*m_energy*m_energyL/sqr(APHYTOOLS::Flavour(kf::e).PSMass());
  m_xi     = m_nonlin1 + m_nonlin2 * m_energy;
  m_xe    /= (1+m_xi);
  m_xmax   = m_xe/(1.+m_xe);
  m_xmax2  = 2.*m_xe/(1.+2.*m_xe);

  if (m_mode==0) m_upper = m_xmax;
            else m_upper = m_xmax2;
  m_peak   = m_xmax;

  m_yfix   = 1./(1.+m_xe);
  m_yden   = log(1.+m_xe);
  m_ysteps = 50;

  m_totalC = 0.7115863 - 0.6776124e-3 * m_energy + 0. * m_energy * m_energy; 
  m_total2 = m_totalC * 0.5540019 * (1.-exp(-37.38912 * m_xi * m_xi));
  m_totalE = m_totalC * (0.7257064 + 1.517959e-3 * m_energy);

  msg.Tracking()<<"Initialised Laser-Backscattering ("<<m_mode<<") : "<<endl
		<<" xe,xmax = "<<m_xe<<", "<<m_xmax
		<<" for energyL,mass ="<<m_energyL<<", "<<APHYTOOLS::Flavour(kf::e).PSMass()<<endl
		<<" with xi = "<<m_xi<<", norms   = "<<m_totalC<<"  /  "<<m_total2<<endl
		<<" with polarization = "<<m_polarization<<", polarizationL = "<<m_polarizationL<<endl;
}



Beam_Base * Laser_Backscattering::Copy() {
  if (m_nonlin1>0.) return new Laser_Backscattering(m_beam,m_energy,m_polarization,
						    m_energyL,m_polarizationL,m_mode,m_angles,1,m_dir);
  return new Laser_Backscattering(m_beam,m_energy,m_polarization,
				  m_energyL,m_polarizationL,m_mode,m_angles,0,m_dir);
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
    z   = m_xmax2*i*.0007;
    restot = deg  = 0.;
    restot += res1 = Compton(z,m_polarization,m_polarizationL,deg);
    restot += res2 = TwoPhotons(z,m_polarization,m_polarizationL,deg); 
    restot += res3 = Rescattering(z,m_polarization,m_polarizationL,deg);
    if (flag) ofile<<" "<<z<<"  "<<res1<<"  "<<res1+res2<<"  "<<res1+res2+res3;
    else  msg.Out()<<" "<<z<<"  "<<res1<<"  "<<res1+res2<<"  "<<res1+res2+res3;
    if (IsZero(restot)) {deg = 0.;restot = 1.e-17;}
    if (flag) ofile<<"  "<<deg/restot<<endl;
    else  msg.Out()<<"  "<<deg/restot<<endl;
  }

  if (flag) ofile.close();
}

bool Laser_Backscattering::CalculateWeight(double _x,double _scale) 
{
  m_x = _x; m_Q2 = _scale;
  if (!( (_x*m_energy>=m_Ebounds[0]) && (_x*m_energy<=m_Ebounds[1]))) {
    m_weight = 0.;
    return 0;
  }

  m_polar = 0.;
  double spec;
  switch (m_mode) {
  case 1: 
    spec = Compton(_x,m_polarization,m_polarizationL,m_polar);
    break;
  case 2: 
    spec = TwoPhotons(_x,m_polarization,m_polarizationL,m_polar);
    break;
  case 3: 
    spec = Rescattering(_x,m_polarization,m_polarizationL,m_polar);  
    break;
  default:
    spec = Compton(_x,m_polarization,m_polarizationL,m_polar) + 
           TwoPhotons(_x,m_polarization,m_polarizationL,m_polar) + 
           Rescattering(_x,m_polarization,m_polarizationL,m_polar);  
    break;
  }
  m_polar  = m_polar/spec;
  m_weight = spec;

  return 1;
};

double Laser_Backscattering::Weight(Flavour flin)
{
  if (m_weight<=0.) return 0.;
  if (flin != Flavour(kf::photon)) return 0.;
  return m_weight;
}

AMATOOLS::Vec4D Laser_Backscattering::OutMomentum() {
  if (m_angles==0) return m_x*m_vecout;
  AORGTOOLS::msg.Error()<<"Error in Laser_Backscattering::OutMomentum()."<<endl
			<<"    m_angles != 0 not implemented yet."<<endl;
  return m_x*m_vecout; 
}

double Laser_Backscattering::Compton(double x,double pole,double poll,double & deg)
{
  if ((x<0.) || (x>m_xmax) || (m_totalC < 0.) ) {
      return 0.;
  }

  double value  = SimpleCompton(x,m_xe,pole*poll);

  double g2    = m_xe/x - m_xe - 1;
  if (g2<0.) {
    if (m_pol) deg += value * Polarization(x,m_xe,pole,poll);
    return value;
  }

  double damp   = exp(-m_rho2 * g2/8.);
  if (m_pol) deg += damp * value * m_totalC * Polarization(x,m_xe,pole,poll);

  double wt = damp * m_totalC * value;
  return wt;
}

double Laser_Backscattering::TwoPhotons(double x,double pole,double poll,double & deg)
{
  if ((x<0.) || (x>m_xmax2) || (m_total2 < 0.)) return 0.;

  double value  = SimpleCompton(x,2.*m_xe,pole*poll);

  double g2    = 2.*m_xe/x - 2.*m_xe - 1;
  if (g2<0.) {
    if (m_pol) deg += value * m_total2 * Polarization(x,2.*m_xe,pole,poll);
    return value;
  }

  double damp   = exp(-m_rho2 * g2/8.) * pow(g2,m_delta);
  if (m_pol) deg += damp * value * m_total2 * Polarization(x,2.*m_xe,pole,poll);

  return damp * m_total2 * value;
}

double Laser_Backscattering::Rescattering(double x,double pole,double poll,double & deg)
{
  if ((x<0.) || (x>m_xmax) || (m_totalE < 0.)) return 0.;

  double yMin  = Max(m_yfix,0.5 * x * (1.+sqrt(4./(x*m_xe) + 1.)));
  if (yMin > 1.) return 0.;
  
  double y1, y2;
  double dy, dp, value, pvalue;
  double val1,val2,p1,p2;

  value    = pvalue  = 0.;
  y1       = y2      = yMin;
  y1      *= 1.000001;
  dy       = (1.-yMin)/m_ysteps;

  val1     = log(1.+y1*m_xe)/(y1 * m_yden) *
    SimpleCompton(x/y1,y1*m_xe,0.)*SimpleCompton(1-y1,m_xe,pole*poll);
  p1       = Polarization(x/y1,y1*m_xe,pole,poll);

  for (int i=0;i<m_ysteps;i++) {
    y2       += dy;
    val2      = log(1.+y2*m_xe)/(y2 * m_yden) *
      SimpleCompton(x/y2,y2*m_xe,0.)*SimpleCompton(1-y2,m_xe,pole*poll);
    value    += 0.5*(val1+val2)*dy;
    if (m_pol) {
      p2      = Polarization(x/y2,y2*m_xe,pole,poll);
      pvalue += 0.5*(val1*p1+val2*p2)*dy;
      p1      = p2;
    }
    val1    = val2;
  }

  if (m_pol) deg += pvalue*value*m_totalE;
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
  norm        -= pol2 * (2. + z*z/((z+1.)*(z+1.)) - (1.+2./z) * log(z+1.));

  return value/norm;
}

double Laser_Backscattering::Polarization(double x,double z,double pole,double poll)
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














