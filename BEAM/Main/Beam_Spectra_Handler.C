#include "BEAM/Main/Beam_Spectra_Handler.H"
#include "BEAM/Main/Monochromatic.H"
#include "BEAM/Main/Spectrum_Reader.H"
#include "BEAM/Main/Laser_Backscattering.H"
#include "BEAM/Main/EPA.H"
#include "ATOOLS/Org/Run_Parameter.H" 
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Message.H"
#include <stdio.h>
#include "ATOOLS/Org/Exception.H"


using namespace ATOOLS;
using namespace BEAM;
using namespace std;


Beam_Spectra_Handler::Beam_Spectra_Handler(Data_Reader * dataread) : 
  p_BeamBase(NULL) 
{
  p_BeamBase = new Beam_Base*[2];
  for (short int i=0;i<2;i++) p_BeamBase[i] = NULL;

  if (!(SpecifySpectra(dataread) && InitKinematics(dataread))) {
    msg_Error()<<"Error in Beam_Spectra_Handler::Beam_Spectra_Handler :"<<endl
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
  ATOOLS::rpa.gen.SetPBeam(0,p_BeamBase[0]->InMomentum());
  ATOOLS::rpa.gen.SetPBeam(1,p_BeamBase[1]->InMomentum());
}

Beam_Spectra_Handler::~Beam_Spectra_Handler() { 
  for (short int i=0;i<2;i++) {
    if (p_BeamBase[i]) { delete p_BeamBase[i]; p_BeamBase[i] = NULL; }
  }
  if (p_BeamBase) { delete [] p_BeamBase; p_BeamBase = NULL; }
}


bool Beam_Spectra_Handler::SpecifySpectra(Data_Reader * dataread)
{
  bool okay   = 1;
  char help[20];
  int beam_spectrum,spectrum_generator;
  for (short int num=0;num<2;num++) {
    sprintf(help,"%i",num+1);
    std::string number   = string(help); 

    std::string  bs(dataread->GetValue<std::string>("BEAM_SPECTRUM_"+number,"Monochromatic"));
    if (bs=="Monochromatic") beam_spectrum=0;
    else if (bs=="Gaussian") beam_spectrum=1;
    else if (bs=="Laser_Backscattering") beam_spectrum=2;
    else if (bs=="Simple_Compton") beam_spectrum=3;
    else if (bs=="Spectrum_Reader") beam_spectrum=4;
    else if (bs=="EPA") beam_spectrum=5;
    else beam_spectrum = 99;
    if ((beam_spectrum!=Beam_Type::Monochromatic) && (beam_spectrum!=Beam_Type::Gaussian)) 
      spectrum_generator = dataread->GetValue<int>("SPECTRUM_"+number,(int)Beam_Generator::Unknown); 
    switch (beam_spectrum) {
    case Beam_Type::Monochromatic :
      okay = okay&&InitializeMonochromatic(dataread,num);
      break;
    case Beam_Type::Gaussian :
      msg_Error()<<"Error in Beam_Initialization::SpecifySpectra :"<<endl
		 <<"   Gaussian beam spectra still have to be implemented."<<endl 
		 <<"   Will read in parameters, check the procedure and abort later."<<endl;
      okay = 0;
      break;
    case Beam_Type::Simple_Compton :
      dataread->AddCommandLine("LASER_MODE = -1");
    case Beam_Type::Laser_Back :
      okay = okay&&InitializeLaserBackscattering(dataread,num);
      break;
    case Beam_Type::EPA :
      okay = okay&&InitializeEPA(dataread,num);
      break;
    case Beam_Type::Spec_Read :
      okay = okay&&InitializeSpectrumReader(dataread,num);
      break;
    default :
      msg_Error()<<"Warning in Beam_Initialization::SpecifySpectra :"<<endl
		 <<"   No beam spectrum specified for beam "<<num+1<<endl
		 <<"   Will initialize monochromatic beam."<<endl;
      okay = okay&&InitializeMonochromatic(dataread,num);
      break;
    }
  }
  return okay;
}

bool Beam_Spectra_Handler::InitializeLaserBackscattering(Data_Reader * dataread,int num) {
  char help[20];
  sprintf(help,"%i",num+1);
  std::string number        = string(help); 
  int     flav              = dataread->GetValue<int>("BEAM_"+number,0);  
  InitializeFlav((kf_code)abs(flav));
  Flavour beam_particle     = Flavour((kf_code)abs(flav));
  if (flav<0) beam_particle = beam_particle.Bar();
  double  beam_energy       = dataread->GetValue<double>("BEAM_ENERGY_"+number,0.0);
  double  beam_polarization = dataread->GetValue<double>("BEAM_POL_"+number,0.0);

  if ( (beam_particle!=Flavour(kf_e)) && (beam_particle!=Flavour(kf_e).Bar()) ) {
    msg_Error()<<"Error in Beam_Initialization::SpecifySpectra :"<<endl
	       <<"   Tried to initialize Laser_Backscattering for "<<beam_particle<<"."<<endl
	       <<"   This option is not available. Result will be to terminate program."<<endl;
    return 0;
  }      
  double Laser_energy       = dataread->GetValue<double>("E_LASER_"+number,0.0);
  double Laser_polarization = dataread->GetValue<double>("P_LASER_"+number,0.0);
  int mode                  = dataread->GetValue<int>("LASER_MODE",0);
  std::string anglesx       = dataread->GetValue<std::string>("LASER_ANGLES","Off");
  std::string nonlinx       = dataread->GetValue<std::string>("LASER_NONLINEARITY","Off");
  int angles(anglesx=="On"?1:0), nonlin(nonlinx=="On"?1:0);

  p_BeamBase[num]          = new Laser_Backscattering(beam_particle,beam_energy,beam_polarization,
						      Laser_energy,Laser_polarization,
						      mode,angles,nonlin,1-2*num);
  return 1;
}

bool Beam_Spectra_Handler::InitializeSpectrumReader(Data_Reader * dataread,int num) {
  char help[20];
  sprintf(help,"%i",num+1);
  std::string number        = string(help); 
  int     flav              = dataread->GetValue<int>("BEAM_"+number,0);  
  InitializeFlav((kf_code)abs(flav));
  Flavour beam_particle     = Flavour((kf_code)abs(flav));
  if (flav<0) beam_particle = beam_particle.Bar();
  double beam_energy        = dataread->GetValue<double>("BEAM_ENERGY_"+number,0.0);
  double beam_polarization  = dataread->GetValue<double>("BEAM_POL_"+number,0.0);
  double laser_energy       = dataread->GetValue<double>("E_LASER_"+number,0.0);
  double laser_polarization = dataread->GetValue<double>("P_LASER_"+number,0.0);

  std::string spectrumfile  = dataread->GetValue<std::string>("SPECTRUM_FILE_"+number,"");

  p_BeamBase[num] = new Spectrum_Reader(beam_particle,beam_energy,beam_polarization,
					laser_energy, laser_polarization,
					spectrumfile,1-2*num);
  return 1;
}

bool Beam_Spectra_Handler::InitializeMonochromatic(Data_Reader * dataread,int num) {
  char help[20];
  sprintf(help,"%i",num+1);
  std::string number = string(help); 
  int     flav              = dataread->GetValue<int>("BEAM_"+number,0);  
  InitializeFlav((kf_code)abs(flav));
  Flavour beam_particle     = Flavour((kf_code)abs(flav));
  if (flav<0) beam_particle = beam_particle.Bar();
  double  beam_energy       = dataread->GetValue<double>("BEAM_ENERGY_"+number,0.0);
  double  beam_polarization = dataread->GetValue<double>("BEAM_POL_"+number,0.0);
  p_BeamBase[num]           = new Monochromatic(beam_particle,beam_energy,beam_polarization,1-2*num);
  return 1;
}


bool Beam_Spectra_Handler::InitializeEPA(Data_Reader * dataread,int num) 
{
  std::string number(ToString(num+1));
  int     flav              = dataread->GetValue<int>("BEAM_"+number,num+1);  
  InitializeFlav((kf_code)abs(flav));
  Flavour beam_particle     = Flavour((kf_code)(abs(flav)));
  if (flav<0) beam_particle = beam_particle.Bar();

  //  msg.Out()<<" Flavour(kf::p_plus)="<<Flavour(kf::p_plus)<<endl;
    if ( (beam_particle!=Flavour(kf_p_plus)) && (beam_particle!=Flavour(kf_p_plus).Bar()) ) {
    msg_Error()<<"Error in Beam_Initialization::SpecifySpectra :"<<endl
	       <<"   Tried to initialize EPA for "<<beam_particle<<"."<<endl
	       <<"   This option is not available. Result will be to terminate program."<<endl;
    return 0;
    }
    //  msg.Out()<<"InitializeEPA; flav "<<flav<<" kf::code="<< kf::code(flav) << endl; // 
  
  double  beam_energy       = dataread->GetValue<double>("BEAM_ENERGY_"+number,0.0);
  double  beam_polarization = dataread->GetValue<double>("BEAM_POL_"+number,0.0);
  double  beam_mass         = beam_particle.Mass(true);
  double  beam_charge       = beam_particle.Charge();
  p_BeamBase[num]           = new EPA(beam_particle,beam_mass,beam_charge,
				      beam_energy,beam_polarization,1-2*num,dataread);
  return 1;
}

bool Beam_Spectra_Handler::InitKinematics(Data_Reader * dataread) {
 
  // cms system from beam momenta - this is for potential asymmetric collisions.
  Vec4D  P      = p_BeamBase[0]->InMomentum()+p_BeamBase[1]->InMomentum();
  double s      = P.Abs2();
  double E      = sqrt(s);
  rpa.gen.SetEcms(E);
  Read_Write_Base::AddGlobalTag("E_CMS",ToString(rpa.gen.Ecms()));

  m_splimits[0] = s*dataread->GetValue<double>("BEAM_SMIN",1e-10);
  m_splimits[1] = s*ATOOLS::Min(dataread->GetValue<double>("BEAM_SMAX",1.0),Upper1()*Upper2());
  m_splimits[2] = s;
  m_ylimits[0]  = -10.;
  m_ylimits[1]  = 10.;
  m_exponent[0] = .5;
  m_exponent[1] = .98 * ( p_BeamBase[0]->Exponent() + p_BeamBase[1]->Exponent());
  m_mass12      = sqr(p_BeamBase[0]->Bunch().Mass());
  m_mass22      = sqr(p_BeamBase[1]->Bunch().Mass());
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
  msg_Out()<<"Beam_Spectra_Handler : "<<endl
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
      if (! ( ((_beams[i]==Flavour(kf_e)) || (_beams[i]==Flavour(kf_e).Bar())) &&
      (_bunches[i]==Flavour(kf_photon))         ) ) {
      fit = 0;
      break;
      }
      }
      if (p_BeamBase[i]->Type() == string("Beam_Strahlung")) {
      if (! ( ((_beams[i] == Flavour(kf_e)) || (_beams[i] == Flavour(kf_e).Bar())) &&
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
      if (_bunches[i]!=Flavour(kf_photon)) {
      fit = 0;
      break;
      }
      }
      if (p_BeamBase[i]->Type() == string("Beam_Strahlung")) {
      if ((_bunches[i]!=Flavour(kf_e) && _bunches[i]!=Flavour(kf_e).Bar()) ||
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


bool Beam_Spectra_Handler::MakeBeams(Vec4D * p) 
{
  double sprime=m_spkey[3], y=m_ykey[2];
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

void Beam_Spectra_Handler::InitializeFlav(kf_code flav)
{
  if (s_kftable.find(flav)==s_kftable.end()) {
    bool initialize_diquarks(false);
    if (flav==kf_p_plus) {
      s_kftable[flav]=
        new Particle_Info(kf_p_plus,0.938272,0,3,1,1,1,1,"P+","P+");
      initialize_diquarks=true;
    }
    else if (flav==kf_n) {
      s_kftable[flav]=
        new Particle_Info(kf_n,0.939566,7.424e-28,0,0,1,1,1,"n","n");
      initialize_diquarks=true;
    }
    else if (flav==kf_e) {
      s_kftable[flav]=
        new Particle_Info(kf_e,0.000511,.0,-3,-1,0,1,0,1,1,0,"e-","e^-");
    }
    else if (flav==kf_photon) {
      s_kftable[flav]=
        new Particle_Info(22,.0,.0,0,0,0,2,-1,1,1,0,"P","\\gamma");
    }
    else {
      THROW(fatal_error,"You specified a beam particle "+ToString(flav)+
            "which is not contained in your chosen model. Will abort.");
    }
      

    if (initialize_diquarks) {
      s_kftable[1103]=new Particle_Info(1103,0.77133,0,-2,0,-3,2,0,0,1,0,"dd_1","dd_1");
      s_kftable[2101]=new Particle_Info(2101,0.57933,0,1,0,-3,0,0,0,1,0,"ud_0","ud_0");
      s_kftable[2103]=new Particle_Info(2103,0.77133,0,1,0,-3,2,0,0,1,0,"ud_1","ud_1");
      s_kftable[2203]=new Particle_Info(2203,0.77133,0,4,0,-3,2,0,0,1,0,"uu_1","uu_1");
      s_kftable[3101]=new Particle_Info(3101,0.80473,0,-2,0,-3,0,0,0,1,0,"sd_0","sd_0");
      s_kftable[3103]=new Particle_Info(3103,0.92953,0,-2,0,-3,2,0,0,1,0,"sd_1","sd_1");
      s_kftable[3201]=new Particle_Info(3201,0.80473,0,1,0,-3,0,0,0,1,0,"su_0","su_0");
      s_kftable[3203]=new Particle_Info(3203,0.92953,0,1,0,-3,2,0,0,1,0,"su_1","su_1");
      s_kftable[3303]=new Particle_Info(3303,1.09361,0,-2,0,-3,2,0,0,1,0,"ss_1","ss_1");
      s_kftable[4101]=new Particle_Info(4101,1.96908,0,1,0,-3,0,0,0,1,0,"cd_0","cd_0");
      s_kftable[4103]=new Particle_Info(4103,2.00808,0,1,0,-3,2,0,0,1,0,"cd_1","cd_1");
      s_kftable[4201]=new Particle_Info(4201,1.96908,0,4,0,-3,0,0,0,1,0,"cu_0","cu_0");
      s_kftable[4203]=new Particle_Info(4203,2.00808,0,4,0,-3,2,0,0,1,0,"cu_1","cu_1");
      s_kftable[4301]=new Particle_Info(4301,2.15432,0,1,0,-3,0,0,0,1,0,"cs_0","cs_0");
      s_kftable[4303]=new Particle_Info(4303,2.17967,0,1,0,-3,2,0,0,1,0,"cs_1","cs_1");
      s_kftable[4403]=new Particle_Info(4403,3.27531,0,4,0,-3,2,0,0,1,0,"cc_1","cc_1");
      s_kftable[5101]=new Particle_Info(5101,5.38897,0,-2,0,-3,0,0,0,1,0,"bd_0","bd_0");
      s_kftable[5103]=new Particle_Info(5103,5.40145,0,-2,0,-3,2,0,0,1,0,"bd_1","bd_1");
      s_kftable[5201]=new Particle_Info(5201,5.38897,0,1,0,-3,0,0,0,1,0,"bu_0","bu_0");
      s_kftable[5203]=new Particle_Info(5203,5.40145,0,1,0,-3,2,0,0,1,0,"bu_1","bu_1");
      s_kftable[5301]=new Particle_Info(5301,5.56725,0,-2,0,-3,0,0,0,1,0,"bs_0","bs_0");
      s_kftable[5303]=new Particle_Info(5303,5.57536,0,-2,0,-3,2,0,0,1,0,"bs_1","bs_1");
      s_kftable[5401]=new Particle_Info(5401,6.67143,0,1,0,-3,0,0,0,1,0,"bc_0","bc_0");
      s_kftable[5403]=new Particle_Info(5403,6.67397,0,1,0,-3,2,0,0,1,0,"bc_1","bc_1");
      s_kftable[5503]=new Particle_Info(5503,10.07354,0,-2,0,-3,2,0,0,1,0,"bb_1","bb_1");
    }
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
}


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
	msg_Out()<<"ERROR: unpolarised cross section for polarised beam!! "<<endl;
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
  m_spkey.Assign("s' beam",4,0,info);
  m_ykey.Assign("y beam",3,0,info);
  m_xkey.Assign("x beam",5,0,info);
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

