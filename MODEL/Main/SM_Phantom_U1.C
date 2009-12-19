#include "MODEL/Main/SM_Phantom_U1.H"
#include "MODEL/Main/Standard_Model.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "MODEL/Main/Effective_Higgs_Coupling.H"
#include "ATOOLS/Org/Message.H"
#include <iomanip>

using namespace MODEL;
using namespace ATOOLS;
using namespace std;

DECLARE_GETTER(SM_Phantom_U1_Getter,"SM+Phantom_U1",Model_Base,Model_Arguments);

Model_Base *SM_Phantom_U1_Getter::operator()(const Model_Arguments &args) const
{
  return new SM_Phantom_U1(args.m_path,args.m_file,args.m_elementary);
}

void SM_Phantom_U1_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"The Standard Model + U(1) phantom Higgs\n"
     <<std::setw(width+4)<<" "<<"{\n"
     <<std::setw(width+7)<<" "<<"parameter specification [keyword=value]\n"
     <<std::setw(width+7)<<" "<<"- all the SM parameters\n"
     <<std::setw(width+7)<<" "<<"- ...\n"
     <<std::setw(width+4)<<" "<<"}\n";
}


SM_Phantom_U1::SM_Phantom_U1(std::string _dir,std::string _file,bool _elementary) :
  Model_Base(_dir,_file,_elementary)
{
  p_sm = new Standard_Model(m_dir,m_file,false);
  ParticleInit();
  if (m_elementary) {
    ATOOLS::OutputParticles(msg->Info());
    ATOOLS::OutputContainers(msg->Info());
  }
}

bool SM_Phantom_U1::ModelInit()
{
  if (m_elementary)
    msg_Info()<<"Initialize the Standard Model plus U(1) phantom Higgs from "
	      <<m_dir<<" / "<<m_file<<std::endl;
  m_name      = std::string("SM+Phantom_U1");

  p_numbers          = p_sm->ExtractScalarNumbers();
  p_constants        = p_sm->ExtractScalarConstants();
  p_complexconstants = p_sm->ExtractComplexConstants();
  p_functions        = p_sm->ExtractScalarFunctions();
  p_matrices         = p_sm->ExtractComplexMatrices();

  delete p_sm;

  FillSpectrum();

  if (!SanityChecks()) {
    msg_Error()<<"Potential Error in "<<METHOD<<":"<<endl
	       <<"   Sanity checks not passed."<<endl
	       <<"   Continue and hope for the best."<<endl;
  }
  return true;
}

SM_Phantom_U1::~SM_Phantom_U1() 
{ }

void SM_Phantom_U1::ParticleInit() {
  //add Higgs particles and new gauge boson
  //kf_code,mass,width,3*charge,icharge,strong,2*spin,majorana,take,stable,massive,idname,tex_name
  s_kftable[kf_H0]   = new Particle_Info(kf_H0,1000.,10.0,0,0,0,0,-1,1,1,1,"H0","H_0");
  s_kftable[kf_A0]   = new Particle_Info(kf_A0,1000.,10.0,0,0,0,0,-1,1,1,1,"A0","A_0");
  s_kftable[kf_Z0_2] = new Particle_Info(kf_Z0_2,1000.,10.0,0,0,0,2,-1,0,1,1,"Z'","Z'");
  //add pseudo gluon
  s_kftable[kf_shgluon] = new Particle_Info(kf_shgluon,0.,0.0,0,0,8,2,-1,1,1,0,"shgluon","shgluon",1);

  ReadParticleData();
}

void SM_Phantom_U1::FillSpectrum() {
  p_dataread = new Data_Reader(" ",";","!","=");
  p_dataread->AddComment("#");
  p_dataread->AddWordSeparator("\t");
  p_dataread->SetInputPath(m_dir);
  p_dataread->SetInputFile(m_file);
  p_constants->insert(make_pair(string("Tan(Beta)"),    
				p_dataread->GetValue<double>("Tan(Beta)",1.)));
  p_constants->insert(make_pair(string("Tan(Theta)"),    
				p_dataread->GetValue<double>("Tan(Theta)",0.)));
  p_constants->insert(make_pair(string("M_H1"),    
				p_dataread->GetValue<double>("M_H1",-1.)));
  p_constants->insert(make_pair(string("M_H2"),    
				p_dataread->GetValue<double>("M_H2",-1.)));
  p_constants->insert(make_pair(string("M_Z'"),    
				p_dataread->GetValue<double>("M_Z'",-1.)));
  p_constants->insert(make_pair(string("g'_1"),    
				p_dataread->GetValue<double>("g'_1",0.)));

  CMatrix HiggsMix(2);
  HiggsMix[0][0] = HiggsMix[1][1] = sqrt(1./(1.+sqr(ScalarConstant(string("Tan(Theta)")))));
  HiggsMix[0][1] = sqrt(1.-sqr(abs(HiggsMix[1][1])));
  HiggsMix[1][0] = -HiggsMix[0][1];

  p_matrices->insert(std::make_pair(std::string("HiggsMix"),HiggsMix));
  
  Flavour flav;
  flav = Flavour(kf_h0);
  flav.SetMass(ScalarConstant(string("M_H1")));
  flav.SetHadMass(ScalarConstant(string("M_H1")));
  flav.SetMassOn(true);
  flav.SetStable(false);
  flav.SetWidth(-1.);
  flav = Flavour(kf_H0);
  flav.SetMass(ScalarConstant(string("M_H2")));
  flav.SetHadMass(ScalarConstant(string("M_H2")));
  flav.SetMassOn(true);
  flav.SetStable(false);
  flav.SetWidth(-1.);
  flav = Flavour(kf_A0);
  flav.SetMass(0.);
  flav.SetMassOn(true);
  flav.SetWidth(0.);
  flav = Flavour(kf_Z0_2);
  flav.SetMass(ScalarConstant(string("M_Z'")));
  flav.SetHadMass(ScalarConstant(string("M_Z'")));
  flav.SetMassOn(true);
  flav.SetWidth(0.);


  //Effective coupling for Higgs-Gluon-Gluon / Higgs-3 Gluon /Higgs-4 Gluon vertices 
  Complex ehf(0.,0.),eHf(0.,0.),ehW(0.,0.),eHW(0.,0.);
  //top
  double mt(Flavour(kf_t).Mass()), mW(Flavour(kf_Wplus).Mass()),
    mh(Flavour(kf_h0).Mass()), mH(Flavour(kf_H0).Mass());
  Effective_Higgs_Coupling ehc(mh);
  ehf = ehc.GetFermionContribution(mt,false);
  ehW = ehc.GetVectorContribution(mW);
  Effective_Higgs_Coupling eHc(mH);
  eHf = eHc.GetFermionContribution(mt,false);
  eHW = ehc.GetVectorContribution(mW);

  //std::cout<<"Check this.  Photons :"
  //	   <<(3.*ehf*sqr(Flavour(kf_t).Charge())+ehW*sqr(Flavour(kf_Wplus).Charge()))<<" / "
  //	   <<(3.*eHf*sqr(Flavour(kf_t).Charge())+eHW*sqr(Flavour(kf_Wplus).Charge()))
  //	   <<" and gluons : "
  //	   <<ehf<<" / "<<eHf<<std::endl;
  p_complexconstants->insert(std::make_pair(std::string("h0_gg_fac"),ehf));
  p_complexconstants->insert(std::make_pair(std::string("H0_gg_fac"),eHf));
  p_complexconstants->insert(std::make_pair(std::string("h0_pp_fac"),
  					    (3.*ehf*sqr(Flavour(kf_t).Charge())+
  					     ehW*sqr(Flavour(kf_Wplus).Charge()))));
  p_complexconstants->insert(std::make_pair(std::string("H0_pp_fac"),
  					    (3.*eHf*sqr(Flavour(kf_t).Charge())+
  					     eHW*sqr(Flavour(kf_Wplus).Charge()))));
}

bool SM_Phantom_U1::SanityChecks() {
  if (ScalarConstant(string("Tan(Beta)"))<1.e-6 ||
      ScalarConstant(string("M_H1"))<0.       || 
      ScalarConstant(string("M_H2"))<0.) return false;
  return true;
}



