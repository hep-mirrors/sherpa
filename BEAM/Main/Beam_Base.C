#include "BEAM/Main/Beam_Base.H"

#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/My_MPI.H"

using namespace BEAM;

Beam_Base::Beam_Base(std::string _type,const ATOOLS::Flavour _beam,
		     const double _energy,const double _polarisation,
		     const int _dir) :
  m_type(_type), m_beam(_beam), m_bunch(m_beam), m_dir(_dir), 
  m_energy(_energy), m_polarisation(_polarisation),
  m_x(1.), m_Q2(0.), m_weight(1.)
{
  Init();
}

bool Beam_Base::Init(int mode)
{
  double disc      =  mode?1.0:1.0-ATOOLS::sqr(m_beam.Mass()/m_energy);
  if (disc<0) {
    msg_Error()<<"Error in Beam_Base :"<<m_type<<std::endl
		       <<"   Mismatch of energy and mass of beam particle : "
		       <<m_beam<<" / "<<m_energy<<std::endl
		       <<"   Will lead to termination of program."<<std::endl;
    ATOOLS::Abort();
  }
  m_lab    = ATOOLS::Vec4D(m_energy,0.,0.,m_dir*m_energy*sqrt(disc));
  m_vecout = ATOOLS::Vec4D(m_energy,0.,0.,m_dir*m_energy*sqrt(disc));
  return true;
}

Beam_Base::~Beam_Base() 
{
}

ATOOLS::Flavour Beam_Base::Beam()         
{
  return m_beam; 
}

ATOOLS::Flavour Beam_Base::Bunch()        
{
  return m_bunch; 
}

ATOOLS::Flavour Beam_Base::Remnant()      
{
  return m_beam; 
}

ATOOLS::Vec4D Beam_Base::OutMomentum()  
{
  return m_vecout; 
}

ATOOLS::Vec4D Beam_Base::InMomentum()   
{
  return m_lab; 
}

double Beam_Base::Energy()       
{
  return m_energy; 
}

double Beam_Base::Polarisation() 
{
  return m_polarisation; 
}

double Beam_Base::Weight(ATOOLS::Flavour fl) 
{ 
  return m_weight; 
}

bool Beam_Base::On()           
{
  return 0; 
}

bool Beam_Base::PolarisationOn() 
{
  return (m_polarisation!=0.); 
}

std::string Beam_Base::Type()         
{
  return m_type; 
}

double Beam_Base::Exponent()     
{
  return 0.; 
}

double Beam_Base::Xmax()         
{ 
  return 1.; 
}

double Beam_Base::Peak()         
{
  return 1.; 
}

