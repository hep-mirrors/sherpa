#include "MODEL/Main/FeynRules_Model.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "MODEL/Main/Running_Fermion_Mass.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "MODEL/Main/FeynRules_Spectrum.H"
#include "ATOOLS/Org/Message.H"
#include "MODEL/Main/Spectrum_Generator_Base.H"
#include <iomanip>

using namespace MODEL;
using namespace ATOOLS;
using namespace std;

DECLARE_GETTER(FeynRules_Model_Getter,"FeynRules",Model_Base,Model_Arguments);

Model_Base *FeynRules_Model_Getter::operator()(const Model_Arguments &args) const
{
  return new FeynRules_Model(args.m_path,args.m_file,args.m_elementary);
}

void FeynRules_Model_Getter::PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"Physical Model defined through FeynRules inputs\n";
  str<<std::setw(width+4)<<" "<<"{\n"
     <<std::setw(width+7)<<" "<<"parameter specification [keyword=value]\n"
     <<std::setw(width+7)<<" "<<"- FR_PARTICLES (particle data input file) \n"
     <<std::setw(width+7)<<" "<<"- FR_IDENTFILE ()\n"
     <<std::setw(width+7)<<" "<<"- FR_PARAMFILE ()\n"
     <<std::setw(width+7)<<" "<<"- FR_PARAMCARD ()\n"
     <<std::setw(width+7)<<" "<<"- FR_INTERACTIONS ()\n"
     <<std::setw(width+4)<<" "<<"}";
}


FeynRules_Model::FeynRules_Model(std::string _dir,std::string _file,bool _elementary) :
  Model_Base(_dir,_file,_elementary)
{
  ParticleInit();
  CustomContainerInit();
  if (m_elementary) {
    msg_Info()<<"\n";
    ATOOLS::OutputParticles(msg->Info());
    ATOOLS::OutputContainers(msg->Info());
  }
}

bool FeynRules_Model::ModelInit(const PDF::ISR_Handler_Map& isr)
{
  p_dataread->RereadInFile();
  if (m_elementary)
    msg_Info()<<"Initialize a FeynRules Model from "<<m_dir<<" / "<<m_file<<std::endl;
  m_name             = std::string("FeynRules");
  p_numbers          = new ScalarNumbersMap();
  p_constants        = new ScalarConstantsMap();
  p_complexconstants = new ComplexConstantsMap();
  p_functions        = new ScalarFunctionsMap();
  p_matrices         = new ComplexMatricesMap();

  FillSpectrum(isr);

  rpa->gen.AddCitation
    (1,"The Sherpa interface to FeynRules is published under \\cite{Christensen:2009jx}.");
  return true;
}

void FeynRules_Model::ParticleInit() {
  m_partfile = p_dataread->GetValue<string>("FR_PARTICLES",std::string("Particle.dat"));
  Data_Reader reader = Data_Reader(" ",";","!","=");
  p_dataread->AddComment("#");
  reader.AddWordSeparator("\t");
  reader.SetAddCommandLine(false);
  reader.SetInputPath(m_dir);
  reader.SetInputFile(m_partfile);
  
  kf_code kfc;
  int     charge,icharge,spin,strong,Majorana;
  bool    Take,stable,massive;
  double  mass,width;
  std::string idname, texname;
  
  if (!reader.OpenInFile()) {
    msg_Error()<<METHOD<<"(): FeynRules particle data file '"
	       <<m_partfile<<"' not found."<<std::endl;
    return;
  }
  std::vector<std::vector<std::string> > helpsvv;
  reader.MatrixFromFile(helpsvv);
  for(size_t i(1);i<helpsvv.size();++i) {
    if (helpsvv[i].size()!=13) {
      msg_Error()<<METHOD<<"(): Inconsistent entry in line "<<i
  		 <<" of '"<<reader.InputFile()<<"'."<<std::endl;
      continue;
    }

    if (helpsvv[i][0]==std::string("kf")) continue;

    kfc=ToType<int>(helpsvv[i][0]); 
    mass=ToType<double>(helpsvv[i][1]); width=ToType<double>(helpsvv[i][2]);
    charge=ToType<int>(helpsvv[i][3]); icharge=ToType<int>(helpsvv[i][4]);
    strong=ToType<int>(helpsvv[i][5]); spin=ToType<int>(helpsvv[i][6]);
    Majorana=ToType<int>(helpsvv[i][7]); Take=ToType<int>(helpsvv[i][8]);
    stable=ToType<int>(helpsvv[i][9]); massive=ToType<int>(helpsvv[i][10]);
    idname=helpsvv[i][11]; texname=helpsvv[i][12];
    s_kftable[kfc] = new Particle_Info(kfc,mass,width,charge,icharge,strong,spin,
				       Majorana,Take,stable,massive,idname,texname);
  }
  ReadParticleData();
  
  // add containers
  s_kftable[kf_none] = new
    Particle_Info(kf_none,-1,0,0,0,0,0,-1,0,1,0,"no_particle","no particle");
  s_kftable[kf_resummed] = new
    Particle_Info(kf_resummed,0.,0.,0,0,1,2,0,1,1,0,"r","resummed");
  s_kftable[kf_bjet] = new
    Particle_Info(kf_bjet,0.,0.,0,0,1,2,0,1,1,0,"bj","bjet");

  s_kftable[kf_fermion] = new
    Particle_Info(kf_fermion,0.,0., 0,0,0,1,0,1,1,0,"fermion","fermion");
  s_kftable[kf_jet] = new
    Particle_Info(kf_jet,0.,0.,0,0,1, 2,0,1,1,0,"j","jet");
  s_kftable[kf_quark] = new
    Particle_Info(kf_quark,0.,0.,0, 0,1,1,0,1,1,0,"Q","Quark");
  s_kftable[kf_lepton] = new
    Particle_Info(kf_lepton,0.,0.,-3,-1,0,1,0,1,1,0,"lepton","lepton");
  s_kftable[kf_neutrino] = new
    Particle_Info(kf_neutrino,0.,0.,0,1,0, 1,0,1,1,0,"neutrino","neutrino");
  s_kftable[kf_fermion]->Clear();
  s_kftable[kf_jet]->Clear();
  s_kftable[kf_resummed]->Clear();
  s_kftable[kf_resummed]->m_resummed=true;
  s_kftable[kf_quark]->Clear();
  s_kftable[kf_lepton]->Clear();
  s_kftable[kf_neutrino]->Clear();
  for (int i=1;i<7;i++) {
    Flavour addit((kf_code)i);
    if (addit.Mass()==0.0 || (!addit.IsMassive() && addit.IsOn())) {
      s_kftable[kf_jet]->Add(addit);
      s_kftable[kf_jet]->Add(addit.Bar());
      s_kftable[kf_quark]->Add(addit);
      s_kftable[kf_quark]->Add(addit.Bar());
      s_kftable[kf_fermion]->Add(addit);
      s_kftable[kf_fermion]->Add(addit.Bar());
    }
  }
  s_kftable[kf_jet]->Add(Flavour(kf_gluon));
  s_kftable[kf_jet]->SetResummed();
  for (int i=11;i<17;i+=2) {
    Flavour addit((kf_code)i);
    if (addit.Mass()==0.0 || (!addit.IsMassive() && addit.IsOn())) {
      s_kftable[kf_lepton]->Add(addit);
      s_kftable[kf_lepton]->Add(addit.Bar());
      s_kftable[kf_fermion]->Add(addit);
      s_kftable[kf_fermion]->Add(addit.Bar());
    }
  }
  for (int i=12;i<18;i+=2) {
    Flavour addit((kf_code)i);
    if (addit.Mass()==0.0 && addit.IsOn()) {
      s_kftable[kf_neutrino]->Add(addit);
      s_kftable[kf_neutrino]->Add(addit.Bar());
      s_kftable[kf_fermion]->Add(addit);
      s_kftable[kf_fermion]->Add(addit.Bar());
    }
  }
}


void FeynRules_Model::FillSpectrum(const PDF::ISR_Handler_Map& isr) {
  RunSpectrumGenerator(isr);
}

void FeynRules_Model::RunSpectrumGenerator(const PDF::ISR_Handler_Map& isr) {
  p_spectrumgenerator = new FeynRules_Spectrum(p_dataread,this,m_dir);
  p_spectrumgenerator->Run(isr);
}  

bool FeynRules_Model::CheckFlavours(int nin, int nout, Flavour* flavs)
{
  return true;
}
