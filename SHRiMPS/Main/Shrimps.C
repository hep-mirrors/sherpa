#include "SHRiMPS/Main/Shrimps.H"
#include "SHRiMPS/Main/Hadron_Init.H"
#include "SHRiMPS/Eikonals/Eikonal_Creator.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Phys/Flavour.H"
#include "MODEL/Main/Strong_Coupling.H"
#include "MODEL/Main/Model_Base.H"
#include <string>
#include <vector>

using namespace SHRIMPS;
using namespace std;

Shrimps::Shrimps(ATOOLS::Data_Reader * dr,
		 BEAM::Beam_Spectra_Handler *const beam,
		 PDF::ISR_Handler *const isr) :  
  p_beamremnants(NULL), p_generator(NULL)
{
  ATOOLS::rpa->gen.AddCitation(1,"SHRiMPS is not published yet.");
  MBpars.Init(dr);
  if (MBpars.RunMode()==run_mode::unknown) {
    msg_Error()<<"Error in "<<METHOD<<":\n   unknown runmode.  Will exit.\n";
    exit(0);
  }
  else if (MBpars.RunMode()==run_mode::xsecs_only) {
    msg_Events()<<METHOD<<": run Shrimps to generate cross sections only.\n"
		<<"   Will write out results and exit afterwards.\n";
    GenerateXsecs();
    exit(1);
  }
  InitialiseTheRun(beam,isr);
  if (MBpars.RunMode()==run_mode::test) {
    TestShrimps();
    exit(1);
  }
}

Shrimps::~Shrimps() 
{
  if (p_beamremnants) delete p_beamremnants;
  if (p_generator)    delete p_generator;
}

void Shrimps::InitialiseTheRun(BEAM::Beam_Spectra_Handler *const beam,
			       PDF::ISR_Handler *const isr) {
  ResetTheFramework();
  InitialiseFormFactors();
  InitialiseSingleChannelEikonals();
  InitialiseBeamRemnants(beam,isr);
  InitialiseTheEventGenerator();  
  //Hadron_Init hadroninit;
  //hadroninit.Init();
}

void Shrimps::ResetTheFramework() {
  m_pdfs.clear();
}

void Shrimps::InitialiseFormFactors() {
  for (size_t i=0;i<MBpars.NGWStates();i++) {
    FormFactor_Parameters params(MBpars.GetFFParameters());
    params.number = i;
    if (i==1) params.kappa *= -1.;
    Form_Factor * ff = new Form_Factor(params);
    ff->Initialise();
    MBpars.AddFormFactor(ff);
  }
  msg_Info()<<METHOD<<" done.\n";
}

void Shrimps::InitialiseSingleChannelEikonals() 
{
  msg_Info()<<METHOD<<" for "<<MBpars.GetFormFactors()->size()
	    <<" form factors.\n";
  Eikonal_Creator creator;
  list<Form_Factor *> * ffs(MBpars.GetFormFactors());
  for (list<Form_Factor *>::iterator ff1=ffs->begin();ff1!=ffs->end();ff1++) {
    for (list<Form_Factor *>::iterator ff2=ffs->begin();ff2!=ffs->end();ff2++) {
      creator.SetFormFactors((*ff1),(*ff2));
      Omega_ik * eikonal(creator.InitialiseEikonal());
      MBpars.AddEikonal(eikonal);
    }
  }
  msg_Info()<<METHOD<<" done.\n";
}

void Shrimps::InitialiseBeamRemnants(BEAM::Beam_Spectra_Handler * const beam,
				     PDF::ISR_Handler *const isr) { 
  for (size_t i=0;i<2;i++) 
    m_pdfs.push_back(Continued_PDF(isr->PDF(i),isr->Flav(i)));

  p_beamremnants = new Beam_Remnant_Handler(beam,m_pdfs);
}

void Shrimps::InitialiseTheEventGenerator() {
  Cross_Sections xsecs;
  xsecs.CalculateCrossSections();
  p_generator = new Event_Generator();
  p_generator->Initialise(p_beamremnants);
}

int Shrimps::GenerateEvent(ATOOLS::Blob_List * blobs) {
  return p_generator->MinimumBiasEvent(blobs);
}

ATOOLS::Return_Value::code Shrimps::FillBeamBlobs(ATOOLS::Blob_List * blobs) {
  return p_beamremnants->FillBeamBlobs(blobs,p_generator->GetEikonal(),
				       p_generator->Smin());
}

void Shrimps::CleanUp(const size_t & mode) {
  p_generator->Reset();
  p_beamremnants->Reset(mode);
}

void Shrimps::GenerateXsecs() {
  msg_Out()<<METHOD<<".\n";  
  InitialiseFormFactors();
  std::vector<double> energies;
  ReadEnergiesFromFile(energies);

  std::vector<double> xsectot, xsecinel, xsecelas;
  Cross_Sections xsecs;
  for (size_t i;i<energies.size();i++) {
    double energy = energies[i];
    MBpars.UpdateForNewEnergy(energy);
    InitialiseSingleChannelEikonals();
    xsecs.CalculateCrossSections();
    xsectot.push_back(xsecs.SigmaTot()/1.e9);
    xsecinel.push_back(xsecs.SigmaInel()/1.e9);
    xsecelas.push_back(xsecs.SigmaEl()/1.e9);
  }
  WriteOutXSecsYodaFile(energies, xsectot, xsecinel, xsecelas);
}

// if (Elastics.find((*energy))!=Elastics.end()) {
//       Elastic_Event_Generator elastic(m_cross.GetSigmaElastic(),NULL,-1);
//       m_cross.GetSigmaElastic()->PrintDifferentialelasticXsec(false,tuning,
// 								dirname);
//       m_cross.GetSigmaSD()->PrintDifferentialElasticAndSDXsec(false,dirname);
//       m_cross.GetSigmaDD()->PrintDifferentialElasticAndDiffXsec(false,
// 								  dirname);
//     }
//   }
//   was.close();
// }
// else {
//     std::vector<double> Elastics;
//     Elastics.push_back(62.5);
//     Elastics.push_back(546.);
//     Elastics.push_back(1800.);
//     Elastics.push_back(7000.);
//     std::vector<double> xsectot, xsecinel,xsecelas;
//     for (int i=0; i<Energies.size(); i++) {
//       InitialiseSingleChannelEikonals((Energies[i]));
//       InitialiseCrossSections((Energies[i]));
//     }
//     for (int i=0; i<Elastics.size(); i++) {
//       InitialiseSingleChannelEikonals((Elastics[i]));
//       InitialiseCrossSections((Elastics[i]));
//       m_cross.GetSigmaElastic()->PrintDifferentialelasticXsec(false,tuning,
// 							      dirname);
//     }
//   }
//   */
// }

void Shrimps::WriteOutXSecsYodaFile(const std::vector<double> & energies,
				    const std::vector<double> & xsectot,
				    const std::vector<double> & xsecinel,
				    const std::vector<double> & xsecelas) {
  std::string dirname = std::string("InclusiveQuantities");
  ATOOLS::MakeDir(dirname);
  std::string filename(dirname+std::string("/xsecs.dat"));
  std::ofstream was;
  was.open(filename.c_str());
  was<<"# BEGIN HISTOGRAM /XSECS/d01-x01-y01\n";
  was<<"AidaPath=/XSECS/d01-x01-y01"<<std::endl;
  for (size_t i=0; i<energies.size(); i++){
    was<<energies[i]<<"   "<<energies[i]<<"   "<<xsectot[i]<<"   0.0\n";
  }
  was<<"# END HISTOGRAM\n"<<std::endl;
  was<<"# BEGIN HISTOGRAM /XSECS/d02-x01-y01\n";
  was<<"AidaPath=/XSECS/d02-x01-y01"<<std::endl;
  for (size_t i=0; i<energies.size(); i++){
    was<<energies[i]<<"   "<<energies[i]<<"   "<<xsecinel[i]<<"   0.0\n";
  }
  was<<"# END HISTOGRAM\n"<<std::endl;
  was<<"# BEGIN HISTOGRAM /XSECS/d03-x01-y01\n";
  was<<"AidaPath=/XSECS/d03-x01-y01"<<std::endl;
  for (size_t i=0; i<energies.size(); i++){
    was<<energies[i]<<"   "<<energies[i]<<"   "<<xsecelas[i]<<"   0.0\n";
  }
  was<<"# END HISTOGRAM"<<std::endl;
  was.close();
}
  
void Shrimps::ReadEnergiesFromFile(std::vector<double> & energies) {
  std::string infile("energies_xsecs.dat");
  std::ifstream input;
  input.open(infile.c_str());
  if(!input){
    msg_Error()<<"File "<<infile<<" does not exist, will exit now.\n";
    exit(1);
  }
  std::string test;
  while (!input.eof()) {
    input>>test;
    energies.push_back(std::atof(test.c_str()));
  }
  input.close();
}



void Shrimps::TestShrimps() {
  msg_Info()<<"Start testing SHRiMPS.\n";
  std::string dirname = std::string("Tests");
  ATOOLS::MakeDir(dirname);

  PrintAlphaS(dirname);
  PrintPDFs(dirname);
  MBpars.GetFormFactors()->front()->Test(dirname); 
  TestEikonalGrids(dirname);
  TestCrossSections(dirname);
  
  msg_Info()<<"Tests done.  Results to be found in "<<dirname<<".\n";
}

void Shrimps::PrintPDFs(const std::string & dirname) {
  int nxval(100);
  double xmin(1.e-5),x;
  for (int i=0; i<5; i++){
    double Q2 = double(i)/2.;
    std::ostringstream ostr;ostr<<Q2;std::string Q2str = ostr.str();    
    std::string filename(dirname+"/pdfs_"+Q2str+".dat");
    std::ofstream was;
    was.open(filename.c_str());
    was<<"# x   u   ubar   d   dbar  s   g"<<std::endl;
    was<<"# Q^2 = "<<Q2<<" GeV^2"<<std::endl;
    for (int j=0;j<=nxval; j++){
      x = pow(10.,double(j)/double(nxval)*log10(xmin));
      m_pdfs[0].Calculate(x,Q2);
      was<<x<<"   "
	 <<m_pdfs[0].XPDF(ATOOLS::Flavour(kf_u))<<"   "
	 <<m_pdfs[0].XPDF(ATOOLS::Flavour(kf_u).Bar())<<"   "
	 <<m_pdfs[0].XPDF(ATOOLS::Flavour(kf_d))<<"   "
	 <<m_pdfs[0].XPDF(ATOOLS::Flavour(kf_d).Bar())<<"   "
	 <<m_pdfs[0].XPDF(ATOOLS::Flavour(kf_s))<<"   "
	 <<m_pdfs[0].XPDF(ATOOLS::Flavour(kf_gluon))<<"\n";
    }
    was.close();
  }
}

void Shrimps::PrintAlphaS(const std::string & dirname) {
  int    nQ2val(1000);
  double Q2max(ATOOLS::sqr(100.)),Q2min(ATOOLS::sqr(1e-3)),Q2;
  double logstepsize((log(Q2max)-log(Q2min))/nQ2val);
  MODEL::Strong_Coupling * alphaS(static_cast<MODEL::Strong_Coupling *>
	   (MODEL::s_model->GetScalarFunction(std::string("strong_cpl"))));

  std::string filename(dirname+"/alphas.dat");
  std::ofstream was;
  was.open(filename.c_str());
  was<<"# Q [GeV]    alpha_s(Q^2)"<<"\n";
  for (int i=0; i<nQ2val; i++){
    Q2 = exp(log(Q2min) + i*logstepsize);
    was<<sqrt(Q2)<<"    "<<(*alphaS)(Q2)<<std::endl;
  }
  was.close();
}

void Shrimps::TestEikonalGrids(const std::string & dirname) {
  Form_Factor * ff(MBpars.GetFormFactors()->front());
  double Delta(MBpars.GetEikonalParameters().Delta);
  double Ymax(MBpars.GetEikonalParameters().Ymax);
  Analytic_Contributor ana12(ff,Delta,Ymax,+1);
  Analytic_Contributor ana21(ff,Delta,Ymax,-1);  
  Omega_ik * eikonal(MBpars.GetEikonals()->front());
  eikonal->TestIndividualGrids(&ana12,&ana21,Ymax,dirname);

  Analytic_Eikonal anaeik;
  eikonal->TestEikonal(&anaeik,dirname);
}

void Shrimps::TestCrossSections(const std::string & dirname) {
  Cross_Sections cross;
  cross.CalculateCrossSections();
  cross.Test(dirname);
}



