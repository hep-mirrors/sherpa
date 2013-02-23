#include "MODEL/Main/SM_MSSM_EHC.H"
#include "MODEL/Main/Standard_Model.H"
#include "MODEL/Main/THDM.H"
#include "MODEL/Main/MSSM.H"
#include "MODEL/Main/Effective_Higgs_Coupling.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Data_Reader.H"
#include <iomanip>

using namespace MODEL;
using namespace ATOOLS;

DECLARE_GETTER(SM_EHC,"SM+EHC",Model_Base,Model_Arguments);

Model_Base *Getter<Model_Base,Model_Arguments,SM_EHC>::
operator()(const Model_Arguments &args) const
{
  return new SM_EHC(args.m_path,args.m_file,args.m_elementary);
}

void Getter<Model_Base,Model_Arguments,SM_EHC>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"The Standard Model + Effective Higgs Couplings to Photons/Gluons\n"
     <<std::setw(width+4)<<" "<<"{\n"
     <<std::setw(width+7)<<" "<<"parameter specification [keyword=value]\n"
     <<std::setw(width+7)<<" "<<"- all the SM parameters\n"
     <<std::setw(width+7)<<" "<<"- FINITE_TOP_MASS (0(1) neglect(consider) top mass effects in loop)\n"
     <<std::setw(width+7)<<" "<<"- FINITE_W_MASS (0(1) neglect(consider) W mass effects in loop)\n"
     <<std::setw(width+7)<<" "<<"- DEACTIVATE_GGH (0(1) deactivate(activate) ggH coupling)\n"
     <<std::setw(width+7)<<" "<<"- DEACTIVATE_PPH (0(1) deactivate(activate) ppH coupling)\n"
     <<std::setw(width+4)<<" "<<"}\n";
}

SM_EHC::SM_EHC(std::string _dir,std::string _file,bool _elementary) :
  Model_Base(_dir,_file,_elementary)
{
  p_sm = new Standard_Model(m_dir,m_file,false);
  ParticleInit();
  if (m_elementary) {
    ATOOLS::OutputParticles(msg->Info());
    ATOOLS::OutputContainers(msg->Info());
  }
}

bool SM_EHC::ModelInit(const PDF::ISR_Handler_Map& isr)
{
  if (m_elementary)
    msg_Info()<<"Initialize the Standard Model \\w EHC from "<<m_dir<<" / "<<m_file<<std::endl;
  m_name      = std::string("SM+EHC");
  p_numbers   = new ScalarNumbersMap();
  p_constants = new ScalarConstantsMap();
  p_complexconstants = new ComplexConstantsMap();
  p_functions = new ScalarFunctionsMap();
  p_matrices  = new ComplexMatricesMap();

  p_sm->ModelInit(isr);
  p_numbers   = p_sm->ExtractScalarNumbers();
  p_constants = p_sm->ExtractScalarConstants();
  p_complexconstants = p_sm->ExtractComplexConstants();
  p_functions = p_sm->ExtractScalarFunctions();
  p_matrices  = p_sm->ExtractComplexMatrices();

  delete p_sm;
  
  FillSpectrum(isr);

  return true;
}

void SM_EHC::ParticleInit() {
  //add pseudo gluon
  //kf_code,mass,width,3*charge,icharge,strong,2*spin,majorana,take,stable,massive,idname,tex_name
  s_kftable[kf_shgluon] = new Particle_Info(kf_shgluon,0.,0.0,0,0,8,2,-1,1,1,0,"shgluon","shgluon",1);
  
  //particle data read by SM
}

void SM_EHC::FillSpectrum(const PDF::ISR_Handler_Map& isr) {
  p_dataread->RereadInFile();
  
  //Effective coupling for Higgs-Gluon-Gluon / Higgs-3 Gluon /Higgs-4 Gluon vertices 
  Complex eh(2./3.,0.);
  double ph(2.*47./18.); double pfac(1.);
  if (p_dataread->GetValue<int>("FINITE_TOP_MASS",0)==1) {
    double hm=Flavour(kf_h0).Mass();
    double tm=Flavour(kf_t).Mass();
    Effective_Higgs_Coupling ehc(hm);
    eh = ehc.GetFermionContribution(tm);
    double taut(sqr(hm/2./tm));
    pfac += -56./705.*taut - 32./987.*sqr(taut);
  }
  if (p_dataread->GetValue<int>("FINITE_W_MASS",0)==1) {
    double hm=Flavour(kf_h0).Mass();
    double Wm=Flavour(kf_Wplus).Mass();
    double tauW(sqr(hm/2./Wm));
    pfac += 66./235.*tauW + 228./1645.*sqr(tauW) + 696./8225.*tauW*sqr(tauW)
            + 5248./90475.*sqr(sqr(tauW)) + 1280./29939.*tauW*sqr(sqr(tauW))
            + 54528./1646645.*sqr(tauW*sqr(tauW));
  }
  if (p_dataread->GetValue<int>("DEACTIVATE_GGH",0)==1) eh=Complex(0.,0.);
  if (p_dataread->GetValue<int>("DEACTIVATE_PPH",0)==1) ph=0.;
  p_constants->insert(std::make_pair(std::string("h0_pp_fac"),ph*pfac));
  p_constants->insert(std::make_pair(std::string("h0_gg_fac"),real(eh)));
}


DECLARE_GETTER(MSSM_EHC,"MSSM+EHC",Model_Base,Model_Arguments);

Model_Base *Getter<Model_Base,Model_Arguments,MSSM_EHC>::
operator()(const Model_Arguments &args) const
{
  return new MSSM_EHC(args.m_path,args.m_file,args.m_elementary);
}

void Getter<Model_Base,Model_Arguments,MSSM_EHC>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"MSSM + Effective Higgs Coupling to Gluons\n"
    <<std::setw(width+4)<<" "<<"{\n"
     <<std::setw(width+7)<<" "<<"parameter specification [keyword=value]\n"
     <<std::setw(width+7)<<" "<<"- all the MSSM parameters\n"
     <<std::setw(width+7)<<" "<<"- FINITE_TOP_MASS (0(1) neglect(consider) top mass effects in loop)\n"
     <<std::setw(width+4)<<" "<<"}\n";
}

MSSM_EHC::MSSM_EHC(std::string _dir,std::string _file,bool _elementary) :
  Model_Base(_dir,_file,_elementary)
{
  p_mssm = new MSSM(m_dir,m_file,false);
  ParticleInit();
}

bool MSSM_EHC::ModelInit(const PDF::ISR_Handler_Map& isr)
{
  if (m_elementary)
    msg_Info()<<"Initialize the MSSM \\w EHCs from "<<m_dir<<" / "<<m_file<<std::endl;
  m_name      = std::string("MSSM+EHC");
  p_numbers   = new ScalarNumbersMap();
  p_constants = new ScalarConstantsMap();
  p_complexconstants = new ComplexConstantsMap();
  p_functions = new ScalarFunctionsMap();
  p_matrices  = new ComplexMatricesMap();
 
  p_mssm->ModelInit(isr);
  p_numbers   = p_mssm->ExtractScalarNumbers();
  p_constants = p_mssm->ExtractScalarConstants();
  p_complexconstants = p_mssm->ExtractComplexConstants();
  p_functions = p_mssm->ExtractScalarFunctions();
  p_matrices  = p_mssm->ExtractComplexMatrices();

  delete p_mssm;
  
  FillSpectrum(isr);

  if (m_elementary) {
    ATOOLS::OutputParticles(msg->Info());
    ATOOLS::OutputContainers(msg->Info());
  }

  return true;
}

void MSSM_EHC::ParticleInit() {
  //add pseudo gluon
  //kf_code,mass,width,3*charge,icharge,strong,2*spin,majorana,take,stable,massive,idname,tex_name
  s_kftable[kf_shgluon] = new Particle_Info(kf_shgluon,0.,0.0,0,0,8,2,-1,1,1,0,"shgluon","shgluon",1);
  //particle data read by SM
}


void MSSM_EHC::FillSpectrum(const PDF::ISR_Handler_Map& isr) {
  p_dataread->RereadInFile();
  
  //Effective coupling for Higgs-Gluon-Gluon / Higgs-3 Gluon /Higgs-4 Gluon vertices 

  double sina = real(ComplexMatrixElement(std::string("Z_R"),1,1)); 
  double cosa = real(ComplexMatrixElement(std::string("Z_R"),0,1)); 
  double sinb = real(ComplexMatrixElement(std::string("Z_H"),1,1)); 
  double cosb = real(ComplexMatrixElement(std::string("Z_H"),1,0)); 

  int fm = p_dataread->GetValue<int>("FINITE_TOP_MASS",0);

  //h0
  Complex eh0(0.,0.);          
  { //top
    double hm=Flavour(kf_h0).Mass();
    Effective_Higgs_Coupling ehc(hm);
    double mass = Flavour(kf_t).Mass();
    if (!fm) mass=-1.;
    eh0 += (cosa/sinb)*ehc.GetFermionContribution(mass);
  }
  //to be added: squarks!

  Complex eH0(0.,0.);
  { //top
    double hm=Flavour(kf_H0).Mass();
    Effective_Higgs_Coupling ehc(hm);
    double mass = Flavour(kf_t).Mass();
    if (!fm) mass=-1.;
    eH0 += (sina/sinb)*ehc.GetFermionContribution(mass);
  }
  //to be added: squarks!

  Complex eA0(0.,0.);
  { //top
    double hm=Flavour(kf_A0).Mass();
    Effective_Higgs_Coupling ehc(hm);
    double mass = Flavour(kf_t).Mass();
    if (!fm) mass=-1.;
    eA0 += (cosb/sinb)*ehc.GetFermionContribution(mass,1);
  }

  p_constants->insert(std::make_pair(std::string("h0_gg_fac"),real(eh0)));
  p_constants->insert(std::make_pair(std::string("H0_gg_fac"),real(eH0)));
  p_constants->insert(std::make_pair(std::string("A0_gg_fac"),real(eA0)));
}


