#include "BEAM/Main/EPA.H"

#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Math/Gauss_Integrator.H"

#include <string>
#include <fstream>

using namespace BEAM;
using namespace ATOOLS;

EPA::EPA(const Flavour _beam,const double _mass,
	 const double _charge,const double _energy,
	 const double _pol,const int _dir,Data_Reader *const read):
  Beam_Base("EPA",_beam,_energy,_pol,_dir),
  m_mass(_mass), m_charge(_charge)
{ 
  m_bunch=Flavour(kf_photon);
  m_vecout=Vec4D(m_energy,0.,0.,_dir*m_energy);
  std::string num(_dir>0?"1":"2");
  m_q2Max=read->GetValue<double>("EPA_q2Max_"+num,2.0);
  m_pt_min=read->GetValue<double>("EPA_ptMin_"+num,0.0);
  m_aqed=read->GetValue<double>("EPA_AlphaQED",0.0072992701);
  m_debug=read->GetValue<int>("EPA_Debug",0);
  m_formfactor=read->GetValue<int>("EPA_Form_Factor_"+num,m_beam.FormFactor());

  if (m_pt_min>1.0) {
    /* pt_min > 1 - according to approximation of 'qmi' calculation in CalculateWeight */
    THROW(critical_error,"Too big p_T cut ( "+ToString(m_pt_min)+")");
  }
  if (m_debug) {
    std::string filename = 
      read->GetValue<std::string>("EPA_Debug_File",
                                  "EPA_debugOutput");
    filename += num + ".log";
    this->selfTest(filename); 
  }
}

EPA::~EPA() 
{
}

double EPA::CosInt::GetCosInt(double X)
{
  if (X<0.) exit(1);
  ATOOLS::Gauss_Integrator integrator(this);
  return integrator.Integrate(X,100000.,1.e-4,1);
}

double EPA::phi(double x, double qq)
{
  if (abs(m_beam.Kfcode()) == kf_p_plus) {
    const double a = 7.16;
    const double b = -3.96;
    const double c = .028;
    double y,qq1,f;
    qq1=1+qq;
    y= x*x/(1-x);
    f=(1+a*y)*(-log(qq1/qq)+1/qq1+1/(2*qq1*qq1)+1/(3*qq1*qq1*qq1));
    f+=(1-b)*y/(4*qq*qq1*qq1*qq1);
    f+=c*(1+y/4)*(log((qq1-b)/qq1)+b/qq1+b*b/(2*qq1*qq1)+b*b*b/(3*qq1*qq1*qq1));
    return f;
  }
  if (m_beam.IsIon()) {
    // x := omega / omega0 is assumed in the following code!
    // ensure whether calls of phi for ions are done correctly 
    // x_omega=x*E/omega0=x*E*R/gamma
    double f = 0.;
    // needed for gaussian shaped nucleus
    const double q0 = 0.06;
    const int atomicNumber = m_beam.GetAtomicNumber();
    const double radius = 1.2/.197*pow(atomicNumber, 1./3.);
    CosInt Ci;
    // do form factor dependent calculation
    switch (m_formfactor) {
    //switch (2) {
    case 0: // point-like form factor
      f= log(1. + (1./(x*x)))/2. + 1./(1.+(1./(x*x)))/2.-1./2.;
      break;
    case 1: // homogeneously charged sphere
      f+= 3. / (16. * pow(x, 6.));
      f+= 3. / (8. * pow(x, 4.));
      f-= cos(2.*x)*3./(16*pow(x, 6.)) + cos(2.*x)*7./(40.*x*x);
      f-= cos(2.*x)*1./20.;
      f-= sin(2.*x)*3./(8*pow(x,5.)) +  sin(2.*x)*1./(10.*x*x*x);
      f+= sin(2.*x)*9./(20.*x) -  sin(2.*x)*x/10.;
      f-= Ci.GetCosInt(2.*x) * (1. + pow(x, 5.)/5.); // integral-cosine
      break;
    case 2: // gaussian shaped nucleus
      f=(1.+x*x / (q0*q0*radius*radius) );
      f*=ExpIntegral(1, x*x / (q0*q0*radius*radius));
      f-=exp(-x*x / (q0*q0*radius*radius));
      f/=2.;
      break;
    case 3: // homogeneously charged sphere (smooth function at low and high x)
      if (x < 0.003) { // make n(x) smooth at low x
        f=1.83698*pow(x,-0.00652101)*M_PI*m_energy; 
        //f=1.36549*pow(x,-0.059967)*M_PI*m_energy*atomicNumber;
        // prefactor*c*x^a with c and a from a fit to x_omega*n(x_omega)
        f/=(2*m_aqed*m_charge*m_charge*radius*m_beam.Mass());
      }
      else if (x > 1.33086) { // cut off oscillating parts at high x
        f=0.;
      }
      else { // normal homogenously charged sphere
	f+= 3. / (16. * pow(x, 6.));
	f+= 3. / (8. * pow(x, 4.));
	f-= cos(2.*x)*3./(16*pow(x, 6.)) + cos(2.*x)*7./(40.*x*x);
	f-= cos(2.*x)*1./20.;
	f-= sin(2.*x)*3./(8*pow(x,5.)) +  sin(2.*x)*1./(10.*x*x*x);
	f+= sin(2.*x)*9./(20.*x) -  sin(2.*x)*x/10.;
	f-= Ci.GetCosInt(2.*x) * (1. + pow(x, 5.)/5.); // integral-cosine
      }
      break;
    default:
      THROW(fatal_error,"Unknown ion form factor chosen");
    }
    return (double)f;
  }
  return 0.;
}

void EPA::selfTest(std::string filename)
{
  std::ofstream debugOutput;
  debugOutput.open(filename.c_str());

  debugOutput << "# EPA::selfTest() starting ..." << std::endl;

  // select output format
  debugOutput.setf(std::ios::scientific, std::ios::floatfield);
  debugOutput.precision(10);

  double x_omega=.1e-2;
  const int atomicNumber = m_beam.GetAtomicNumber();
  const double radius = 1.2/.197*pow(atomicNumber, 1./3.);
  double omega0, gamma;
  gamma = m_energy / m_beam.Mass(); 
  //gamma = m_energy * atomicNumber / m_beam.Mass(); 
  // energy is defined as sqrt[s_NN], N=nucleon
  // but recalculated already in the Beam_Spectra_Handler
  omega0 = gamma / radius;

  // write parameters
  debugOutput << "# Form Factor: " << m_formfactor << std::endl;
  debugOutput << "# A= " << atomicNumber << std::endl;
  debugOutput << "# R= " << radius << std::endl;
  debugOutput << "# E= " << m_energy << std::endl;
  debugOutput << "# Z= " << m_charge << std::endl;
  debugOutput << "# M_Ion=" << m_beam.Mass() << std::endl;
  debugOutput << "# gamma= " << gamma << std::endl;
  debugOutput << "# omega0= " << omega0 << std::endl;

  // write spectrum
  while (x_omega < 5) {
    x_omega*=1.005;
    CalculateWeight(x_omega*omega0/m_energy, 0); // m_weight = n(x)
    debugOutput << x_omega << "\t" << x_omega * m_weight << std::endl;
  }

  debugOutput << "# EPA::selfTest() finished" << std::endl << std::endl;
  debugOutput.close();
  return;
}


bool EPA::CalculateWeight(double x,double q2)
{
  // x = omega/E = (E-E')/E  ; E,E' - incoming and outgoing protons energy
  //                           omega = E-E' - energy of emitted photon

  if (abs(m_beam.Kfcode()) == kf_p_plus) {
    const double alpha = m_aqed;
    const double qz = 0.71;
    m_x = x; m_Q2 = q2;
    if (x==1) {
      m_weight=0.0;
      return 1;
    }
    double f, qmi, qma;
    qma=m_q2Max/qz;
    // x = omega/E = (E-E')/E  ; E,E' - incoming and outgoing protons energy
    //                           omega = E-E' - energy of emitted photon
    qmi= m_mass*m_mass * x*x /(1-x)/qz;
    qmi+=m_pt_min*m_pt_min /(1-x)/qz;

    f = alpha/M_PI*(phi(x,qma)-phi(x,qmi))*(1-x)/x;
    f *= m_charge*m_charge;
    if (f < 0) f = 0.;
    m_weight = f;
    return 1;
  }
  else if (m_beam.IsIon()) { // n(x)
    const double alpha = m_aqed;
    m_x = x; m_Q2 = q2;
    const int atomicNumber = m_beam.GetAtomicNumber();
    const double radius = 1.2/.197*pow(atomicNumber, 1./3.);
    double f, omega0, gamma;
    gamma = m_energy / m_beam.Mass();
    //gamma = m_energy * atomicNumber / m_beam.Mass();
    // energy is defined as sqrt[s_NN], N=nucleon
    // but recalculated already in the Beam_Spectra_Handler
    omega0 = gamma / radius;
    /*
    std::cout << "radius=" << radius << std::endl;
    std::cout << "omega0=" << omega0 << std::endl;
    std::cout << "gamma=" << gamma << std::endl;
    */
    /*
    f = 2 * alpha * m_charge * m_charge / M_PI / (m_x * omega0);
    f *= phi(m_x, m_Q2);
    */
    f = 2 * alpha * m_charge * m_charge / M_PI / m_x;
    // since CalculateWeight() is dn=N(x)*dx/x and not dn=N(omega)*domega/omega
    //f = 2 * alpha * m_charge * m_charge / M_PI / (m_x * m_energy);
    f *= phi(m_x*m_energy/omega0, m_Q2); // phi(x_omega, m_Q2)
    // x_omega=m_x*m_energy/omega0 
 
    m_weight = f;
    return 1;
  }
  return 0;
}

Beam_Base *EPA::Copy() 
{
  return new EPA(*this);
}

double EPA::Weight(Flavour fl)                
{ 
  return m_weight; 
}

ATOOLS::Vec4D EPA::OutMomentum()              
{ 
  return m_x*m_vecout; 
}

ATOOLS::Flavour EPA::Remnant()                
{ 
  return m_beam;
}



