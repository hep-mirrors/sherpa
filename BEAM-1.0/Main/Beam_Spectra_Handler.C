#include "Beam_Spectra_Handler.H"
#include "Monochromatic.H"
#include "Spectrum_Reader.H"
#include "Laser_Backscattering.H"
#include "Run_Parameter.H" 
#include "Message.H"
#include <stdio.h>


using namespace ATOOLS;
using namespace BEAM;
using namespace std;


Beam_Spectra_Handler::Beam_Spectra_Handler(Data_Read * dataread) : 
  p_BeamBase(NULL) 
{
  p_BeamBase = new Beam_Base*[2];
  for (short int i=0;i<2;i++) p_BeamBase[i] = NULL;

  if (!(SpecifySpectra(dataread) && InitKinematics(dataread))) {
    msg.Error()<<"Error in Beam_Spectra_Handler::Beam_Spectra_Handler :"<<endl
	       <<"    Could not init spectra or kinematics. Abort program."<<endl;
    abort();
  }

  m_mode = 0;
  m_polarisation = 0;
  for (short int i=0;i<2;i++) {
    if (p_BeamBase[i]->On()) m_mode += i+1;
    if (p_BeamBase[i]->PolarisationOn()) m_polarisation += i+1;
  }
  ATOOLS::rpa.gen.SetBeam1(p_BeamBase[0]->Beam());
  ATOOLS::rpa.gen.SetBeam2(p_BeamBase[1]->Beam());
}

Beam_Spectra_Handler::~Beam_Spectra_Handler() { 
  for (short int i=0;i<2;i++) {
    if (p_BeamBase[i]) { delete p_BeamBase[i]; p_BeamBase[i] = NULL; }
  }
  if (p_BeamBase) { delete [] p_BeamBase; p_BeamBase = NULL; }
}


bool Beam_Spectra_Handler::SpecifySpectra(Data_Read * dataread)
{
  bool okay   = 1;
  char help[20];
  int beam_spectrum,spectrum_generator;
  for (short int num=0;num<2;num++) {
    sprintf(help,"%i",num+1);
    std::string number   = string(help); 

    beam_spectrum        = dataread->GetValue<Beam_Type::code>("BEAM_SPECTRUM_"+number);
    if ((beam_spectrum!=Beam_Type::Monochromatic) && (beam_spectrum!=Beam_Type::Gaussian)) 
      spectrum_generator = dataread->GetValue<Beam_Generator::code>("SPECTRUM_"+number); 
    switch (beam_spectrum) {
    case Beam_Type::Monochromatic :
      okay = okay&&InitializeMonochromatic(dataread,num);
      break;
    case Beam_Type::Gaussian :
      msg.Error()<<"Error in Beam_Initialization::SpecifySpectra :"<<endl
		 <<"   Gaussian beam spectra still have to be implemented."<<endl 
		 <<"   Will read in parameters, check the procedure and abort later."<<endl;
      okay = 0;
      break;
    case Beam_Type::Simple_Compton :
      dataread->SetValue("LASER_MODE","-1");
    case Beam_Type::Laser_Back :
      okay = okay&&InitializeLaserBackscattering(dataread,num);
      break;
    case Beam_Type::Spec_Read :
      okay = okay&&InitializeSpectrumReader(dataread,num);
      break;
    default :
      msg.Error()<<"Warning in Beam_Initialization::SpecifySpectra :"<<endl
		 <<"   No beam sprectum specified for beam "<<num+1<<endl
		 <<"   Will initialize monochromatic beam."<<endl;
      okay = okay&&InitializeMonochromatic(dataread,num);
      break;
    }
  }
  return okay;
}

bool Beam_Spectra_Handler::InitializeLaserBackscattering(Data_Read * dataread,int num) {
  char help[20];
  sprintf(help,"%i",num+1);
  std::string number        = string(help); 
  int     flav              = dataread->GetValue<int>("BEAM_"+number);  
  Flavour beam_particle     = Flavour(kf::code(abs(flav)));
  if (flav<0) beam_particle = beam_particle.Bar();
  double  beam_energy       = dataread->GetValue<double>("BEAM_ENERGY_"+number);
  double  beam_polarization = dataread->GetValue<double>("BEAM_POL_"+number);

  if ( (beam_particle!=Flavour(kf::e)) && (beam_particle!=Flavour(kf::e).Bar()) ) {
    msg.Error()<<"Error in Beam_Initialization::SpecifySpectra :"<<endl
	       <<"   Tried to initialize Laser_Backscattering for "<<beam_particle<<"."<<endl
	       <<"   This option is not available. Result will be to terminate program."<<endl;
    return 0;
  }      
  double Laser_energy       = dataread->GetValue<double>("E_LASER_"+number);
  double Laser_polarization = dataread->GetValue<double>("P_LASER_"+number);
  int mode                  = dataread->GetValue<int>("LASER_MODE");
  int angles                = dataread->GetValue<Switch::code>("LASER_ANGLES");
  int nonlin                = dataread->GetValue<Switch::code>("LASER_NONLINEARITY");

  p_BeamBase[num]          = new Laser_Backscattering(beam_particle,beam_energy,beam_polarization,
						      Laser_energy,Laser_polarization,
						      mode,angles,nonlin,1-2*num);
  return 1;
}

bool Beam_Spectra_Handler::InitializeSpectrumReader(Data_Read * dataread,int num) {
  char help[20];
  sprintf(help,"%i",num+1);
  std::string number        = string(help); 
  int     flav              = dataread->GetValue<int>("BEAM_"+number);  
  Flavour beam_particle     = Flavour(kf::code(abs(flav)));
  if (flav<0) beam_particle = beam_particle.Bar();
  double beam_energy        = dataread->GetValue<double>("BEAM_ENERGY_"+number);
  double beam_polarization  = dataread->GetValue<double>("BEAM_POL_"+number);
  double laser_energy       = dataread->GetValue<double>("E_LASER_"+number);
  double laser_polarization = dataread->GetValue<double>("P_LASER_"+number);

  std::string spectrumfile  = dataread->GetValue<std::string>("SPECTRUM_FILE_"+number);

  p_BeamBase[num] = new Spectrum_Reader(beam_particle,beam_energy,beam_polarization,
					laser_energy, laser_polarization,
					spectrumfile,1-2*num);
  return 1;
}

bool Beam_Spectra_Handler::InitializeMonochromatic(Data_Read * dataread,int num) {
  char help[20];
  sprintf(help,"%i",num+1);
  std::string number = string(help); 
  int     flav              = dataread->GetValue<int>("BEAM_"+number);  
  Flavour beam_particle     = Flavour(kf::code(abs(flav)));
  if (flav<0) beam_particle = beam_particle.Bar();
  double  beam_energy       = dataread->GetValue<double>("BEAM_ENERGY_"+number);
  double  beam_polarization = dataread->GetValue<double>("BEAM_POL_"+number);
  p_BeamBase[num]           = new Monochromatic(beam_particle,beam_energy,beam_polarization,1-2*num);
  return 1;
}


bool Beam_Spectra_Handler::InitKinematics(Data_Read * dataread) {
 
  // cms system from beam momenta - this is for potential assymmetric collisions.
  Vec4D  P      = p_BeamBase[0]->InMomentum()+p_BeamBase[1]->InMomentum();
  double s      = P.Abs2();
  double E      = sqrt(s);
  rpa.gen.SetEcms(E);

  m_splimits[0] = s*dataread->GetValue<double>("BEAM_SMIN");
  m_splimits[1] = s*ATOOLS::Min(dataread->GetValue<double>("BEAM_SMAX"),Upper1()*Upper2());
  m_splimits[2] = s;
  m_ylimits[0]  = -10.;
  m_ylimits[1]  = 10.;
  m_exponent[0] = .5;
  m_exponent[1] = .98 * ( p_BeamBase[0]->Exponent() + p_BeamBase[1]->Exponent());
  m_mass12      = sqr(p_BeamBase[0]->Bunch().PSMass());
  m_mass22      = sqr(p_BeamBase[1]->Bunch().PSMass());
  double x      = 1./2.+(m_mass12-m_mass22)/(2.*E*E);
  double E1     = x*E;
  double E2     = E-E1;
  m_fiXVECs[0]  = Vec4D(E1,0.,0., sqrt(sqr(E1)-m_mass12));
  m_fiXVECs[1]  = Vec4D(E2,0.,0.,-sqrt(sqr(E1)-m_mass12));

  m_asymmetric  = 0;
  if ((dabs((m_fiXVECs[0]+(-1.)*p_BeamBase[0]->InMomentum()).Abs2())>0.0000001) ||
      (dabs((m_fiXVECs[1]+(-1.)*p_BeamBase[1]->InMomentum()).Abs2())>0.0000001) ) m_asymmetric = 1;


  m_type = p_BeamBase[0]->Type() + std::string("*") + p_BeamBase[1]->Type();
  return 1;
}


void Beam_Spectra_Handler::Output() {
  msg.Out()<<"Beam_Spectra_Handler : "<<endl
	   <<"   type = "<<m_type<<endl
	   <<"   for    "<<p_BeamBase[0]->Beam()<<"  ("<<p_BeamBase[0]->InMomentum()<<")"<<endl
	   <<"   and    "<<p_BeamBase[1]->Beam()<<"  ("<<p_BeamBase[1]->InMomentum()<<")"<<endl;
}


bool Beam_Spectra_Handler::CheckConsistency(ATOOLS::Flavour * _beams,
					    ATOOLS::Flavour * _bunches) {
  bool fit = 1;
  for (int i=0;i<2;i++) {
    if ((_beams[i]!=GetBeam(i)->Beam()) || (_bunches[i]!=GetBeam(i)->Bunch())) {
      fit = 0;
      break;
    }
    /*
      if (p_BeamBase[i]->Type() == string("Laser_Backscattering")) {
      if (! ( ((_beams[i]==Flavour(kf::e)) || (_beams[i]==Flavour(kf::e).Bar())) &&
      (_bunches[i]==Flavour(kf::photon))         ) ) {
      fit = 0;
      break;
      }
      }
      if (p_BeamBase[i]->Type() == string("Beam_Strahlung")) {
      if (! ( ((_beams[i] == Flavour(kf::e)) || (_beams[i] == Flavour(kf::e).Bar())) &&
      (_beams[i] == _bunches[i])         ) ) {
      fit = 0;
      break;
      }
      }
      if (p_BeamBase[i]->Type() == string("Monochromatic") ||
      p_BeamBase[i]->Type() == string("Gaussian") ) {
      if (_bunches[i]!=_beams[i]) {
      fit = 0;
      break;
      }
      }
    */
  }
  return fit;
}

bool Beam_Spectra_Handler::CheckConsistency(ATOOLS::Flavour * _bunches) {
  bool fit = 1;
  for (int i=0;i<2;i++) {
    if (_bunches[i]!=GetBeam(i)->Bunch()) {
      fit = 0;
      break;
    }
    /*
      if (p_BeamBase[i]->Type() == string("Laser_Backscattering")) {
      if (_bunches[i]!=Flavour(kf::photon)) {
      fit = 0;
      break;
      }
      }
      if (p_BeamBase[i]->Type() == string("Beam_Strahlung")) {
      if ((_bunches[i]!=Flavour(kf::e) && _bunches[i]!=Flavour(kf::e).Bar()) ||
      (_bunches[i]!=GetBeam(i)->Bunch()) ){
      fit = 0;
      break;
      }
      }
      if (p_BeamBase[i]->Type() == string("Monochromatic") ||
      p_BeamBase[i]->Type() == string("Gaussian") ) {
      if (_bunches[i]!=GetBeam(i)->Bunch()) {
      fit = 0;
      break;
      }
      }
    */
  }
  return fit;
}


bool Beam_Spectra_Handler::MakeBeams(Vec4D * p,double sprime,double y) 
{
  if (m_mode==0) {
    m_x1 = m_x2 = 1.;
    p[0] = m_fiXVECs[0];
    p[1] = m_fiXVECs[1];
    return 1;
  }
  else {
    if ( (sprime<m_splimits[0]) || (sprime>m_splimits[1]) || m_splimits[0]==m_splimits[1] ) {
      return 0;
    }

    double E      = sqrt(m_splimits[2]);
    double Eprime = sqrt(sprime);
    double x      = 1./2.+(m_mass12-m_mass22)/(2.*sprime);
    double E1     = x*Eprime;
    double E2     = Eprime-E1;
    
    p[0]          = Vec4D(E1,0.,0.,sqrt(sqr(E1)-m_mass12));
    p[1]          = Vec4D(E2,(-1.)*Vec3D(p[0]));
    E1            = exp(y);  
    E2            = exp(-y);  


    m_CMSBoost    = Poincare(Vec4D(E1+E2,0.,0.,E1-E2));
    
    Vec4D p1      = p[0];
    Vec4D p2      = p[1];
    m_CMSBoost.BoostBack(p1);
    m_CMSBoost.BoostBack(p2);
    m_x1          = 2.*p1[0]/E;
    m_x2          = 2.*p2[0]/E;

    if (m_mode==1) m_x2 = 1.;
    if (m_mode==2) m_x1 = 1.;

    return 1;
  }
}

/* ----------------------------------------------------------------

   Weight calculation 

   ---------------------------------------------------------------- */


bool Beam_Spectra_Handler::CalculateWeight(double scale) 
{
  switch (m_mode) {
  case 3 :
    if ( (p_BeamBase[0]->CalculateWeight(m_x1,scale)) && 
	 (p_BeamBase[1]->CalculateWeight(m_x2,scale)) ) return 1;
    break;
  case 2 :
    if (p_BeamBase[1]->CalculateWeight(m_x2,scale))     return 1;
    break;
  case 1 :
    if (p_BeamBase[0]->CalculateWeight(m_x1,scale))     return 1;
    break;
  }
  return 0;
};


double Beam_Spectra_Handler::Weight(Flavour * flin)
{
  if (flin==NULL) return (p_BeamBase[0]->Weight() * p_BeamBase[1]->Weight());
  return (p_BeamBase[0]->Weight(flin[0]) * p_BeamBase[1]->Weight(flin[1]));
}


double Beam_Spectra_Handler::Weight(int * pol_types, double *dofs)
{
  double weight = 1.;
  for (int i=0;i<2;++i) {
    if (p_BeamBase[i]->PolarisationOn()) {
      if (pol_types[i]!=99) {
	double hel=(double)pol_types[i];
	double pol=p_BeamBase[i]->Polarisation();
	double dof=dofs[i];
	if (hel*pol>0.) 
	  weight*=(1.+dabs(pol)*(dof-1.))/dof;
	else
	  weight*=(1.-dabs(pol))/dof;

	//assuming 2 degrees of freedom
	//	weight*=dabs(hel+pol)/2.;
      }
      else {
	msg.Out()<<"ERROR: unpolarised cross section for polarised beam!! "<<endl;
      } 
    }
  }
  return weight; 
}

/* ----------------------------------------------------------------

   Boosts

   ---------------------------------------------------------------- */


void  Beam_Spectra_Handler::BoostInCMS(Vec4D* p,int n) {
  for (int i=0; i<n; ++i) m_CMSBoost.Boost(p[i]);
}

void  Beam_Spectra_Handler::BoostInLab(Vec4D* p,int n) {
  for (int i=0; i<n; ++i) m_CMSBoost.BoostBack(p[i]);
}

void   Beam_Spectra_Handler::SetSprimeMin(double _spl)      
{ 
  m_splimits[0]  = Max(m_splimits[0],_spl);
  if (m_splimits[0]>m_splimits[1])  m_splimits[0]=m_splimits[1];
}

void   Beam_Spectra_Handler::SetSprimeMax(double _spl)      { m_splimits[1]  = Min(m_splimits[1],_spl); }

void Beam_Spectra_Handler::AssignKeys(ATOOLS::Integration_Info *const info)
{
  ATOOLS::msg.Tracking()<<"Beam_Spectra_Handler::AssignKeys(..):"
			<<"Creating initial mapping keys ...\n";
  m_spkey.Assign("s' beam",4,0,info);
  m_ykey.Assign("y beam",3,0,info);
  m_xkey.Assign("x beam",5,0,info);
  ATOOLS::msg.Tracking()<<"... done."<<std::endl;
}

void Beam_Spectra_Handler::SetLimits() 
{
  for (short int i=0;i<2;++i) {
    m_spkey[i]=m_splimits[i];
    m_ykey[i]=m_ylimits[i];
  }
  m_spkey[2]=ATOOLS::sqr(ATOOLS::rpa.gen.Ecms());
  m_xkey[0]=-std::numeric_limits<double>::max();
  m_xkey[2]=-std::numeric_limits<double>::max();
  m_xkey[1]=log(Upper1());
  m_xkey[3]=log(Upper2());
}

