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

Shrimps::Shrimps(ATOOLS::Data_Reader * dr,
		 BEAM::Beam_Spectra_Handler *const beam,
		 PDF::ISR_Handler *const isr) :  
  p_generator(NULL)
{
  ATOOLS::rpa->gen.AddCitation(1,"SHRiMPS is not published yet.");
  MBpars.Init(dr);
  if (MBpars.RunMode()==run_mode::unknown) {
    msg_Error()<<"Error in "<<METHOD<<":\n   unknown runmode.  Will exit.\n";
    exit(0);
  }
  if (MBpars.RunMode()==run_mode::xsecs_only) {
    msg_Events()<<METHOD<<": run Shrimps to generate cross sections only.\n"
		<<"   Will write out results and exit afterwards.\n";
    GenerateXsecs();
    exit(1);
  }
  InitialiseTheRun(beam,isr);
  if (MBpars("TestShrimps")!=0) {
    TestShrimps();
    exit(1);
  }
}

Shrimps::~Shrimps() 
{
  while (!m_eikonals.empty()) {
    delete m_eikonals.back();
    m_eikonals.pop_back();
  }
  m_eikonals.clear();
  if (p_beamremnants) delete p_beamremnants;
  if (p_generator) delete p_generator;
}

void Shrimps::InitialiseTheRun(BEAM::Beam_Spectra_Handler *const beam,
			       PDF::ISR_Handler *const isr) {
  ResetTheFramework();
  InitialiseFormFactors();
  InitialiseSingleChannelEikonals(ATOOLS::rpa->gen.Ecms());
  //InitialiseCrossSections(ATOOLS::rpa->gen.Ecms());
  InitialiseBeamRemnants(beam,isr);
  //InitialiseTheEventGenerator();  
  //Hadron_Init hadroninit;
  //hadroninit.Init();
}

void Shrimps::ResetTheFramework() {
  m_pdfs.clear();
  while (!m_eikonals.empty()) {
    delete m_eikonals.back();
    m_eikonals.pop_back();
  }
  m_eikonals.clear();
  m_ffs.clear();
}

void Shrimps::InitialiseFormFactors() {
  msg_Info()<<METHOD<<" for "<<MBpars("NGWstates")<<" states.\n";
  for (size_t i=0;i<MBpars("NGWstates");i++) {
    FormFactor_Parameters params(MBpars.GetFFParameters());
    params.number = i;
    if (i==1) params.kappa *= -1.;
    m_ffs.push_back(Form_Factor(params));
    m_ffs.back().Initialise();
  }
  msg_Info()<<METHOD<<" done.\n";
}

void Shrimps::InitialiseSingleChannelEikonals(const double & Ecms) 
{
  msg_Info()<<METHOD<<"(Y = "<<MBpars("originalY")<<", "
	   <<"deltaY = "<<MBpars("deltaY")<<").\n";
  Eikonal_Creator creator(MBpars.GetEikonalParameters());
  size_t NGWstates(MBpars("NGWstates"));
  for (int i=0;i<NGWstates;i++) {
    for (int j=0;j<NGWstates;j++) {
      msg_Info()<<"   *** create eikonal for channel ["<<i<<" "<<j<<"].\n";
      creator.SetFormFactors(&m_ffs[i],&m_ffs[j]);
      Omega_ik * eikonal(creator.InitialiseEikonal());
      m_eikonals.push_back(eikonal);
    }
  }
  msg_Info()<<METHOD<<" done.\n";
}


void Shrimps::InitialiseCrossSections(const double & energy) {
  m_cross = Cross_Sections(&m_eikonals,energy,m_test);
  m_cross.CalculateTotalCrossSections();
}

void Shrimps::InitialiseBeamRemnants(BEAM::Beam_Spectra_Handler * const beam,
				     PDF::ISR_Handler *const isr) { 
  for (size_t i=0;i<2;i++) 
    m_pdfs.push_back(Continued_PDF(isr->PDF(i),isr->Flav(i)));

  p_beamremnants = new Beam_Remnant_Handler(beam,m_pdfs);
}

void Shrimps::InitialiseTheEventGenerator() {
  //p_generator = new Event_Generator(m_runmode,m_weightmode);
  //p_generator->Initialise(&m_cross,p_beamremnants,m_test);
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

void Shrimps::TestShrimps() {
  msg_Info()<<"Start testing SHRiMPS.\n";
  std::string dirname = std::string("Tests");
  ATOOLS::MakeDir(dirname);

  PrintAlphaS(dirname);
  PrintPDFs(dirname);

  m_ffs[0].Test(dirname); 

  double Delta(MBpars.GetEikonalParameters().Delta);
  double Ymax(MBpars.GetEikonalParameters().Ymax);
  Analytic_Contributor ana12(&m_ffs[0],Delta,Ymax,+1);
  Analytic_Contributor ana21(&m_ffs[0],Delta,Ymax,-1);  
  m_eikonals.front()->TestIndividualGrids(&ana12,&ana21,Ymax,dirname);

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








void Shrimps::GenerateXsecs() {
  InitialiseFormFactors();
  std::string dirname = std::string("InclusiveQuantities");
  ATOOLS::MakeDir(dirname);

  bool tuning(false);
  
  if(!tuning){
    std::list<double> Energies;
    Energies.push_back(50.);
    Energies.push_back(62.5);
    Energies.push_back(100.);
    Energies.push_back(546.);
    Energies.push_back(630.);
    Energies.push_back(900.);
    Energies.push_back(1000.);
    Energies.push_back(1800.);
    Energies.push_back(1960.);
    Energies.push_back(2360.);
    Energies.push_back(7000.);
    Energies.push_back(8000.);
    Energies.push_back(14000.);
    Energies.push_back(55000.);
    Energies.push_back(100000.);
    std::set<double> Elastics;
    Elastics.insert(62.5);
    Elastics.insert(546.);
    Elastics.insert(900.);
    Elastics.insert(1800.);
    Elastics.insert(7000.);


    std::string filename(dirname+std::string("/xsecs_total.dat"));
    std::ofstream was;
    was.open(filename.c_str());
    for (std::list<double>::iterator energy=Energies.begin();
       energy!=Energies.end();energy++) {
      InitialiseSingleChannelEikonals((*energy));
      InitialiseCrossSections((*energy));
      msg_Events()<<"E = "<<ATOOLS::om::red<<(*energy)<<ATOOLS::om::reset;
      msg_Events()<<" sigma_tot = "<<m_cross.SigmaTot()/1.e9
		  <<" sigma_inel = "<<m_cross.SigmaInel()/1.e9
		  <<" sigma_SD = "<<m_cross.SigmaSD()/1.e9
		  <<" sigma_DD = "<<m_cross.SigmaDD()/1.e9
		  <<" sigma_el = "<<m_cross.SigmaEl()/1.e9
		  <<" el.slope = "<<m_cross.ElasticSlope()
		  <<std::endl;
      was<<(*energy)<<"  "
         <<m_cross.SigmaTot()/1.e9<<"  "
         <<m_cross.SigmaInel()/1.e9<<"  "
         <<m_cross.SigmaSD()/1.e9<<"  "
         <<m_cross.SigmaDD()/1.e9<<"  "
         <<m_cross.SigmaEl()/1.e9<<"  "
         <<m_cross.ElasticSlope()<<"  "
         <<std::endl;
      if (Elastics.find((*energy))!=Elastics.end()) {
        Elastic_Event_Generator elastic(m_cross.GetSigmaElastic(),NULL,-1);
        m_cross.GetSigmaElastic()->PrintDifferentialelasticXsec(false,tuning,
								dirname);
        m_cross.GetSigmaSD()->PrintDifferentialElasticAndSDXsec(false,dirname);
        m_cross.GetSigmaDD()->PrintDifferentialElasticAndDiffXsec(false,
								  dirname);
      }
    }
    was.close();
  }
  else {
    std::vector<double> Energies;
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
      Energies.push_back(std::atof(test.c_str()));
    }
    input.close();
    std::vector<double> Elastics;
    Elastics.push_back(62.5);
    Elastics.push_back(546.);
    Elastics.push_back(1800.);
    Elastics.push_back(7000.);

    std::vector<double> xsectot, xsecinel,xsecelas;
    
    for (int i=0; i<Energies.size(); i++) {
      InitialiseSingleChannelEikonals((Energies[i]));
      InitialiseCrossSections((Energies[i]));
      xsectot.push_back(m_cross.SigmaTot()/1.e9);
      xsecinel.push_back(m_cross.SigmaInel()/1.e9);
      xsecelas.push_back(m_cross.SigmaEl()/1.e9);
    }
    std::string filename(dirname+std::string("/xsecs_tuning.dat"));
    std::ofstream was;
    was.open(filename.c_str());
    was<<"# BEGIN HISTOGRAM /XSECS/d01-x01-y01\n";
    was<<"AidaPath=/XSECS/d01-x01-y01"<<std::endl;
    for (int i=0; i<Energies.size(); i++){
      was<<Energies[i]<<"   "<<Energies[i]<<"   "<<xsectot[i]<<"   0.0\n";
    }
    was<<"# END HISTOGRAM\n"<<std::endl;
    was<<"# BEGIN HISTOGRAM /XSECS/d02-x01-y01\n";
    was<<"AidaPath=/XSECS/d02-x01-y01"<<std::endl;
    for (int i=0; i<Energies.size(); i++){
      was<<Energies[i]<<"   "<<Energies[i]<<"   "<<xsecinel[i]<<"   0.0\n";
    }
    was<<"# END HISTOGRAM\n"<<std::endl;
    was<<"# BEGIN HISTOGRAM /XSECS/d03-x01-y01\n";
    was<<"AidaPath=/XSECS/d03-x01-y01"<<std::endl;
    for (int i=0; i<Energies.size(); i++){
      was<<Energies[i]<<"   "<<Energies[i]<<"   "<<xsecelas[i]<<"   0.0\n";
    }
    was<<"# END HISTOGRAM"<<std::endl;
    was.close();

    for (int i=0; i<Elastics.size(); i++) {
      InitialiseSingleChannelEikonals((Elastics[i]));
      InitialiseCrossSections((Elastics[i]));
      m_cross.GetSigmaElastic()->PrintDifferentialelasticXsec(false,tuning,
							      dirname);
    }
  }
}

