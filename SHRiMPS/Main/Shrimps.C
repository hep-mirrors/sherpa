#include "SHRiMPS/Main/Shrimps.H"
#include "SHRiMPS/Main/Hadron_Init.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Phys/Flavour.H"
#include "MODEL/Main/Strong_Coupling.H"
#include "MODEL/Main/Model_Base.H"
#include <string>

using namespace SHRIMPS;

Shrimps::Shrimps(ATOOLS::Data_Reader * dr,
	     BEAM::Beam_Spectra_Handler *const beam,
	     PDF::ISR_Handler *const isr) :
  m_runmode(run_mode::unknown), 
  m_weightmode(weight_mode::unknown),
  m_test(dr->GetValue<int>("TestShrimps",0)),
  p_generator(NULL)
{
  ATOOLS::rpa->gen.AddCitation
    (1,"SHRiMPS is not published yet. \
Please contact the authors if you are using it.");
  double Ecms(ATOOLS::rpa->gen.Ecms());
  MBpars.Init(dr);
  MBpars.Set(std::string("originalY"),
	     log(Ecms/ATOOLS::Flavour(kf_p_plus).HadMass()));

  m_NGWstates  = m_test?1:int(MBpars("NGWstates"));	     
  m_runmode    = MBpars.RunMode();
  m_weightmode = MBpars.WeightMode();
  switch (m_runmode) {
  case run_mode::xsecs_only:
    msg_Events()<<"Run Shrimps to generate cross sections only."<<std::endl;
    GenerateXsecs();
    exit(1);
  case run_mode::elastic_events:
  case run_mode::single_diffractive_events:
  case run_mode::double_diffractive_events:
  case run_mode::quasi_elastic_events:
  case run_mode::inelastic_events:
  case run_mode::all_min_bias:  
  case run_mode::underlying_event:
    InitialiseFormFactors();
    InitialiseSingleChannelEikonals(Ecms);
    InitialiseCrossSections(Ecms);
    InitialiseBeamRemnants(beam,isr);
    InitialiseEventGenerator();
    break;
  case run_mode::unknown: 
  default:
    abort();
  }

  if(m_test){
    PrintPDFandAlphaS();
    exit(1);
  }

  Hadron_Init hadroninit;
  hadroninit.Init();
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

void Shrimps::InitialiseFormFactors() {
  if (m_NGWstates!=2 && m_NGWstates!=1) abort();

  std::vector<double> params;
  params.push_back(m_test?1.:MBpars("FFpref"));
  params.push_back(MBpars("Lambda2"));
  params.push_back(MBpars("beta0"));
  params.push_back(MBpars("kappa"));
  params.push_back(MBpars("xi"));
  params.push_back(MBpars("bmax"));
  params.push_back(MBpars("accu"));

  m_ffs.clear();
  for (int i=0;i<m_NGWstates;i++) {
    if (i==1) params[3] *= -1;
    m_ffs.push_back(Form_Factor(i,m_test));
    m_ffs.back().Initialise(params);
  }
}

void Shrimps::InitialiseSingleChannelEikonals(const double & Ecms) 
{
  MBpars.Set(std::string("originalY"),
	     log(Ecms/ATOOLS::Flavour(kf_p_plus).HadMass()));
  msg_Tracking()<<METHOD<<"(Y = "<<MBpars("originalY")<<", "
		<<"deltaY = "<<MBpars("deltaY")<<").\n";
  while (!m_eikonals.empty()) {
    delete m_eikonals.back();
    m_eikonals.pop_back();
  }
  m_eikonals.clear();
  
  Eikonal_Creator creator(m_test);

  for (int i=0;i<m_NGWstates;i++) {
    for (int j=0;j<m_NGWstates;j++) {
      msg_Tracking()
	<<"Initialise and produce grids for single channel eikonal: "
	<<i<<" "<<j<<" --> "<<&m_ffs[i]<<" "<<&m_ffs[j]<<" ("
	<<ATOOLS::sqr(m_ffs[i].Prefactor()*m_ffs[j].Prefactor())<<").\n";
      m_eikonals.push_back(creator.CreateEikonal(&m_ffs[i],&m_ffs[j]));
    }
  }
  msg_Info()
    <<"Initialised eikonal for Y = "<<MBpars("originalY")<<", "
    <<"deltaY = "<<MBpars("deltaY")<<"."<<std::endl; 
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

void Shrimps::InitialiseEventGenerator() {
  p_generator = new Event_Generator(m_runmode,m_weightmode);
  p_generator->Initialise(&m_cross,p_beamremnants,m_test);
}

int Shrimps::GenerateEvent(ATOOLS::Blob_List * blobs) {
  msg_Tracking()<<"   -->"<<METHOD<<"("<<blobs->size()<<" blobs)\n";
  /*
    if (m_runmode==run_mode::underlying_event) {
    ATOOLS::Blob * blob;
    for (size_t i=0;i<blobs->size();++i) {
    blob = (*blobs)[i];
    msg_Out()<<"   * deal with blob: "<<std::endl<<(*blob)<<std::endl;
    if (blob->Has(ATOOLS::blob_status::needs_beams) &&
    !p_beamremnants->FindInTreatedBlobs(blob)) {
    for (int i=0;i<blob->NOutP();i++) {
    msg_Out()<<"   check ("<<blob->Type()<<"["<<i<<"]: "
    <<blob->OutParticle(i)->Number()<<": ";
    if (blob->OutParticle(i)->DecayBlob()) {
    if (blob->OutParticle(i)->DecayBlob()->Type()==
    ATOOLS::btp::Signal_Process) {
    msg_Out()<<blob->OutParticle(i)->DecayBlob()->Type()<<".\n";
    p_generator->DressShowerBlob(blob);
    msg_Out()<<".\n";
    exit(1);
    }
    msg_Out()<<".\n";
    }
    }
    }
    }
    }
  */
  return p_generator->MinimumBiasEvent(blobs);
}

ATOOLS::Return_Value::code Shrimps::FillBeamBlobs(ATOOLS::Blob_List * blobs) {
  msg_Tracking()
    <<"===================== "<<METHOD<<"("<<blobs->size()<<" blobs) "
    <<"=====================\n";
  return p_beamremnants->FillBeamBlobs(blobs,p_generator->GetEikonal(),
				       p_generator->Smin());
}

void Shrimps::CleanUp(const size_t & mode) {
  p_beamremnants->Reset(mode);
}

void Shrimps::GenerateXsecs() {
  InitialiseFormFactors();
  std::string dirname = std::string("InclusiveQuantities");
  ATOOLS::MakeDir(dirname);

  std::list<double> Energies;
  Energies.push_back(50.);
  Energies.push_back(62.5);
  Energies.push_back(100.);
  //Energies.push_back(200.);
  //Energies.push_back(300.);
  //Energies.push_back(400.);
  Energies.push_back(546.);
  Energies.push_back(630.);
  //Energies.push_back(800.);
  //Energies.push_back(900.);
  Energies.push_back(1000.);
  //Energies.push_back(1200.);
  //Energies.push_back(1500.);
  Energies.push_back(1800.);
  Energies.push_back(1960.);
  Energies.push_back(2360.);
  //Energies.push_back(5000.);
  Energies.push_back(7000.);
  Energies.push_back(8000.);
  Energies.push_back(14000.);
  //Energies.push_back(100000.);
  std::set<double> Elastics;
  Elastics.insert(62.5);
  Elastics.insert(546.);
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
      m_cross.GetSigmaElastic()->PrintDifferentialelasticXsec(false,dirname);
      m_cross.GetSigmaSD()->PrintDifferentialElasticAndSDXsec(false,dirname);
      m_cross.GetSigmaDD()->PrintDifferentialElasticAndDiffXsec(false,dirname);
    }
  }
  was.close();
}

void Shrimps::PrintPDFandAlphaS()
{
  std::string filename("InclusiveQuantities/pdfs.dat");
  std::ofstream was;
  was.open(filename.c_str());
  int nxval(100);
  double x(1.),xmin(1.e-5),Q2(0.),updf,ubarpdf,dpdf,spdf,gpdf;
  was<<"# x   u   ubar   d   s   g"<<std::endl;
  for (int i=0; i<=2; i++){
    Q2 = double(i);
    was<<"# Q^2 = "<<Q2<<" GeV^2"<<std::endl;
    for (int j=0;j<=nxval; j++){
      x = j*(1.-xmin)/(nxval)+xmin;
      m_pdfs[0].Calculate(x,Q2);
      updf    = m_pdfs[0].XPDF(ATOOLS::Flavour(kf_u));
      ubarpdf = m_pdfs[0].XPDF(ATOOLS::Flavour(kf_u).Bar());
      dpdf    = m_pdfs[0].XPDF(ATOOLS::Flavour(kf_d));
      spdf    = m_pdfs[0].XPDF(ATOOLS::Flavour(kf_s));
      gpdf    = m_pdfs[0].XPDF(ATOOLS::Flavour(kf_gluon));
      was<<x<<"   "<<updf<<"   "<<ubarpdf<<"   "<<dpdf<<"   "<<spdf<<"   "<<gpdf<<std::endl;
    }
    was<<std::endl<<std::endl;
  }
  was.close();

  filename="InclusiveQuantities/alphas.dat";
  was.open(filename.c_str());
  int nQ2val(1000);
  double Q2max(ATOOLS::sqr(100.)),Q2min(ATOOLS::sqr(1e-3));
  MODEL::Strong_Coupling * alphaS(static_cast<MODEL::Strong_Coupling *>
	   (MODEL::s_model->GetScalarFunction(std::string("strong_cpl"))));
  was<<"# Q^2 [GeV^2]    alpha_s(Q^2)"<<std::endl;
  for (int i=0; i<nQ2val; i++){
    Q2 = exp(i*(log(Q2max)-log(Q2min))/nQ2val+log(Q2min));
    was<<Q2<<"    "<<(*alphaS)(Q2)<<std::endl;
  }
  was.close();
  return;
}

