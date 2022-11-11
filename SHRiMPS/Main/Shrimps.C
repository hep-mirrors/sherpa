#include "SHRiMPS/Main/Shrimps.H"
#include "SHRiMPS/Main/Hadron_Init.H"
#include "SHRiMPS/Eikonals/Eikonal_Creator.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Phys/Flavour.H"
#include "MODEL/Main/Strong_Coupling.H"
#include "MODEL/Main/Model_Base.H"
#include "SHRiMPS/Glauber/Glauber.H"
#include "SHRiMPS/Glauber/Nucleus.H"
#include <string>
#include <vector>

#include "SHRiMPS/Cross_Sections/Sigma_Total.H"

using namespace SHRIMPS;
using namespace std;


Shrimps::Shrimps(PDF::ISR_Handler *const isr) :
  p_remnants(NULL), p_generator(NULL),
  m_ana(true)
{
  ATOOLS::rpa->gen.AddCitation(1,"SHRiMPS is not published yet.");
  MBpars.Init();
  if (MBpars.RunMode()==run_mode::unknown) {
    msg_Error()<<"Error in "<<METHOD<<":\n   unknown runmode.  Will exit.\n";
    exit(0);
  }
  else if (MBpars.RunMode()==run_mode::xsecs_only) {
    msg_Events()<<METHOD<<": run Shrimps to generate cross sections only.\n"
		<<"   Will write out results and exit afterwards.\n";
    GenerateXsecs();
    exit(0);
  }
  else if (MBpars.RunMode()==run_mode::test) {
    msg_Events()<<METHOD<<": run tests.\n";
    TestShrimps(isr);
  }
  InitialiseTheRun(isr);
  if (m_ana) {
    m_histos[std::string("Yasym_core")]    = new ATOOLS::Histogram(0,  0.0,  8.0, 32);
    m_histos[std::string("Yasym_hard")]    = new ATOOLS::Histogram(0,  0.0,  8.0, 32);
    m_histos[std::string("Yasym_shower")]  = new ATOOLS::Histogram(0,  0.0,  8.0, 32);
    m_histos[std::string("Yasym_soft")]    = new ATOOLS::Histogram(0,  0.0,  8.0, 32);
    m_histos[std::string("Yasym_frag")]    = new ATOOLS::Histogram(0,  0.0,  8.0, 32);
    m_histos[std::string("Yasym_frag_in")] = new ATOOLS::Histogram(0,  0.0,  8.0, 32);
  }
}

Shrimps::~Shrimps() 
{
  if (p_xsecs)     delete p_xsecs;
  if (p_remnants)  delete p_remnants;
  if (p_generator) delete p_generator;
  if (m_ana) {
    std::string name  = std::string("Ladder_Analysis/");
    for (std::map<std::string, ATOOLS::Histogram * >::iterator hit=m_histos.begin();
	 hit!=m_histos.end();hit++) {
      hit->second->Finalize();
      hit->second->Output(name+hit->first);
      delete hit->second;
    }
  }
}

void Shrimps::InitialiseTheRun(PDF::ISR_Handler *const isr) {
  Hadron_Init().Init();
  InitialiseFormFactors();
  InitialiseSingleChannelEikonals();
  InitialiseRemnants(isr);
  InitialiseTheEventGenerator();  
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
}

void Shrimps::InitialiseSingleChannelEikonals() 
{
  msg_Info()<<METHOD<<" for "<<MBpars.GetFormFactors()->size()
  	    <<" form factors.\n";
  Eikonal_Creator creator;
  vector<Form_Factor *> * ffs(MBpars.GetFormFactors());
  MBpars.ResetEikonals(ffs->size());
  for (size_t i=0;i<ffs->size();i++) {
    for (size_t j=0;j<ffs->size();j++) {
      creator.SetFormFactors((*ffs)[i],(*ffs)[j]);
      Omega_ik * eikonal(creator.InitialiseEikonal());
      MBpars.AddEikonal(i,j,eikonal);
      //msg_Out()<<"   * ["<<i<<j<<"]: "<<eikonal<<"\n";
    }
  }
}

void Shrimps::InitialiseRemnants(PDF::ISR_Handler *const isr) { 
  p_remnants = new Remnant_Handler(isr);
}

void Shrimps::InitialiseTheEventGenerator() {
  p_xsecs = new Cross_Sections();
  p_xsecs->CalculateCrossSections();
  p_generator = new Event_Generator(p_xsecs,false);
  p_generator->Initialise(p_remnants,&m_cluster);
  p_remnants->SetColourGenerator(p_generator->GetColourGenerator());
  m_cluster.SetYmax(p_generator->Ymax());
  m_cluster.SetMinKT2(p_generator->MinKT2());
  m_cluster.SetShowerParams(ShowerMode(),ShowerMinKT2());
  m_cluster.SetShowerFac(ShowerFac());
  p_generator->Reset();
  p_remnants->Reset();
}

int Shrimps::InitMinBiasEvent(ATOOLS::Blob_List * blobs) {
  if (blobs->FindFirst(ATOOLS::btp::Fragmentation)!=NULL &&
      blobs->FindFirst(ATOOLS::btp::Hadron_Decay)==NULL) Analyse(blobs);
  return p_generator->InitMinimumBiasEvent(blobs);
}

void Shrimps::Analyse(ATOOLS::Blob_List * blobs) {
  msg_Out()<<"  * "<<METHOD<<"("<<blobs->size()<<" blobs).\n";
  msg_Out()<<"   - "<<METHOD<<"(yhat = "<<p_generator->Yhat()<<")\n";
  m_histos[std::string("Yasym_core")]->Insert(p_generator->Yhat());
  ATOOLS::Blob * hard = blobs->FindFirst(ATOOLS::btp::Hard_Collision);
  Analyse(hard,std::string("Yasym_hard"));  
  ATOOLS::Blob * shower = blobs->FindFirst(ATOOLS::btp::Shower);
  Analyse(shower,std::string("Yasym_shower"));  
  ATOOLS::Blob * soft = blobs->FindFirst(ATOOLS::btp::Soft_Collision);
  Analyse(soft,std::string("Yasym_soft"));  
  ATOOLS::Blob * frag = blobs->FindFirst(ATOOLS::btp::Fragmentation);
  Analyse(frag,std::string("Yasym_frag"));  
  Analyse(frag,std::string("Yasym_frag_in"));  
}

void Shrimps::Analyse(ATOOLS::Blob * blob,std::string tag) {
  msg_Out()<<"   - "<<METHOD<<"("<<blob->Type()<<", "<<blob->NOutP()<<" outgoing particles.)\n";
  if (tag==std::string("Yasym_frag_in")) {
    for (size_t i=0;i<blob->NInP();i++) {
      double y = blob->InParticle(i)->Momentum().Y();
      m_histos[tag]->Insert(ATOOLS::dabs(y),(y>0?1.:-1.));
    }
  }
  else {
    for (size_t i=0;i<blob->NOutP();i++) {
      double y = blob->OutParticle(i)->Momentum().Y();
      m_histos[tag]->Insert(ATOOLS::dabs(y),(y>0?1.:-1.));
    }
  }
}


ATOOLS::Blob * Shrimps::GenerateEvent() {
  msg_Out()<<"  * "<<METHOD<<".\n";
  //Glauber glauberTest(p_xsecs,1,208);
  return p_generator->GenerateEvent();
}

ATOOLS::Cluster_Amplitude * Shrimps::ClusterConfiguration(ATOOLS::Blob *const blob) {
  //m_cluster.SetMinKT2(p_shrimps->ShowerMinKT2());
  //m_cluster.SetRescatt(p_shrimps->IsLastRescatter());
  //m_cluster.SetTMax(p_shrimps->LadderTMax());
  //m_cluster.SetNLad(p_shrimps->NLadders());
  if (!m_cluster.Cluster(blob)) {
    msg_Error()<<"Error in "<<METHOD<<": could not cluster blob.\n"
	       <<(*blob)<<"\n";
    return NULL;
  }
  return m_cluster.Amplitude();
}

ATOOLS::Return_Value::code Shrimps::MakeBeamBlobs(ATOOLS::Blob_List * blobs) {
  p_remnants->SetFormFactors(p_generator->GetEikonal()->FF1(),
			     p_generator->GetEikonal()->FF2());
  return p_remnants->FillBeamBlobs(blobs, p_generator->B());
}

void Shrimps::CleanUp(const size_t & mode) {
  p_generator->Reset();
  p_remnants->Reset();
}

void Shrimps::GenerateXsecs() {
  std::string dirname = std::string("InclusiveQuantities");
  ATOOLS::MakeDir(dirname);

  InitialiseFormFactors();
  std::set<double> energies, energies_tot, energies_inel, energies_el, energies_sd, energies_dd, elastics;
  ReadEnergiesFromFile(energies_tot,"energies_xsecs_tot.dat");
  ReadEnergiesFromFile(energies_inel,"energies_xsecs_inel.dat");
  ReadEnergiesFromFile(energies_el,"energies_xsecs_el.dat");
  ReadEnergiesFromFile(energies_sd,"energies_xsecs_sd.dat");
  ReadEnergiesFromFile(energies_dd,"energies_xsecs_dd.dat");
  ReadEnergiesFromFile(elastics,"energies_elastics.dat");
  energies = energies_tot;
  for (std::set<double>::iterator siter = energies_inel.begin(); siter != energies_inel.end(); ++siter) {
      if (energies.find(*siter) == energies.end()) energies.insert(*siter);
  }
  for (std::set<double>::iterator siter = energies_el.begin(); siter != energies_el.end(); ++siter) {
      if (energies.find(*siter) == energies.end()) energies.insert(*siter);
  }
  for (std::set<double>::iterator siter = energies_sd.begin(); siter != energies_sd.end(); ++siter) {
      if (energies.find(*siter) == energies.end()) energies.insert(*siter);
  }
  for (std::set<double>::iterator siter = energies_dd.begin(); siter != energies_dd.end(); ++siter) {
      if (energies.find(*siter) == energies.end()) energies.insert(*siter);
  }

  std::vector<double> xsectot, xsecinel, xsecelas, xsecsd, xsecdd;

  p_xsecs = new Cross_Sections();
  for (std::set<double>::iterator energy_iter=energies.begin();
       energy_iter!=energies.end();energy_iter++) {
    double energy = (*energy_iter);
    MBpars.UpdateForNewEnergy(energy);
    InitialiseSingleChannelEikonals();
    msg_Info()<<"Calculate cross sections for c.m. energy E = "<<energy<<"\n";
    p_xsecs->CalculateCrossSections();
    if (energies_tot.find(energy) != energies_tot.end()) xsectot.push_back(p_xsecs->SigmaTot()/1.e9);
    if (energies_inel.find(energy) != energies_inel.end()) xsecinel.push_back(p_xsecs->SigmaInel()/1.e9);
    if (energies_el.find(energy) != energies_el.end()) xsecelas.push_back(p_xsecs->SigmaEl()/1.e9);
    if (energies_sd.find(energy) != energies_sd.end()) xsecsd.push_back(p_xsecs->SigmaSD(0)/1.e9);
    if (energies_dd.find(energy) != energies_dd.end()) xsecdd.push_back(p_xsecs->SigmaDD()/1.e9);
    msg_Out()<<"** "<<energy<<" -> "
             <<"xstot = "<<p_xsecs->SigmaTot()<<", "
             <<"xsinel = "<<p_xsecs->SigmaInel()<<", "
             <<"xsel = "<<p_xsecs->SigmaEl()<<", "
             <<"xssd = "<<p_xsecs->SigmaSD(0)<<", "
             <<"xsdd = "<<p_xsecs->SigmaDD()<<".\n";
    if (elastics.find(energy)!=elastics.end()) {
      WriteOutElasticsYodaFile(energy,dirname);
      WriteOutDiffXSecs(energy,dirname);
    }
  }
  WriteOutXSecsYodaFile(energies_tot, energies_inel, energies_el, energies_sd, energies_dd, xsectot, xsecinel, xsecelas, xsecsd, xsecdd, dirname);
}

void Shrimps::WriteOutElasticsYodaFile(const double & energy,
				       std::string dirname) {
    std::string Estring(ATOOLS::ToString(energy));
    std::set<double> tvals;
    std::string infile(std::string("tvals_dsigma_el_dt_"+Estring+".dat"));
    ReadEnergiesFromFile(tvals,infile);

    std::string filename(dirname+std::string("/Dsigma_el_by_dt_"+Estring+".dat"));
    std::ofstream was;
    was.open(filename.c_str());
    was<<"# BEGIN HISTO1D /DSIGMA_EL_BY_DT_"+Estring+"/d01-x01-y01\n";
    was<<"Path=/DSIGMA_EL_BY_DT_"+Estring+"/d01-x01-y01"<<std::endl;

    double value(0.),vallow,valhigh,a,b;
    double t,tlow, thigh;
    const double & tmin = p_xsecs->GetSigmaElastic()->Tmin();
    const double & tmax = p_xsecs->GetSigmaElastic()->Tmax();
    const size_t & steps = p_xsecs->GetSigmaElastic()->Steps();
    const std::vector<double> & eldiffgrid = p_xsecs->GetSigmaElastic()->GetDiffGrid();
    unsigned int ilow, ihigh;

    //     msg_Out()<<"Calculating differential elastic cross sections for tuning."<<std::endl;
    for (std::set<double>::iterator iter=tvals.begin(); iter != tvals.end(); iter++) {
//       msg_Out()<<"calculating for t = "<<tvals[i]<<" GeV^2"<<std::endl;
        t = *iter;
      ilow=int((t-tmin)/(tmax-tmin)*steps);
      if(ilow>steps) ilow=steps-1;
      ihigh=int((t-tmin)/(tmax-tmin)*steps)+1;
      if(ihigh>steps) ilow=steps;
      tlow=tmin+(tmax-tmin)*ilow/steps;
      thigh=tmin+(tmax-tmin)*ihigh/steps;
      vallow=eldiffgrid[ilow];
      valhigh=eldiffgrid[ihigh];
      a=(valhigh-vallow)/(thigh-tlow);
      b=vallow-a*tlow;
      value=a*t+b;

      was<<t<<"   "<<t<<"   "<<value<<"   0.0   0.0\n";
    }
    was<<"# END HISTO1D\n"<<std::endl;
    was.close();

}

void Shrimps::WriteOutDiffXSecs(double energy, std::string dirname){
    std::string Estring(ATOOLS::ToString(energy));
    std::string filename(dirname+std::string("/Differential_XSecs_"+Estring+".dat"));
    std::ofstream was;
    was.open(filename.c_str());
    was<<"# total cross section "+Estring+"\n";
    was<<"# B     dsigma/(2 pi B dB) \n";
    double Bmax(MBpars.GetEikonalParameters().bmax);
    int Nstep(25);

    Sigma_Tot sigmatot;
    sigmatot.Calculate();
    was<<"# 0 - sigma_tot\n";
    for (int i = 0; i < Nstep; ++i) {
        double B(i*Bmax/Nstep);
        double val(sigmatot.GetCombinedValue(B));
        was<<B<<"     "<<val<<endl;
    }
    was<<"\n\n";

    was<<"# 1 - sigma_inel\n";
    for (int i = 0; i < Nstep; ++i) {
        double B(i*Bmax/Nstep);
        double val(p_xsecs->GetSigmaInelastic()->GetCombinedValue(B));
        was<<B<<"     "<<val<<endl;
    }
    was<<"\n\n";

    was<<"# 2 - sigma_el\n";
    for (int i = 0; i < Nstep; ++i) {
        double B(i*Bmax/Nstep);
        double val(p_xsecs->GetSigmaElastic()->GetCombinedValue(B));
        was<<B<<"     "<<val<<endl;
    }
    was<<"\n\n";
    was<<"# 3 - sigma_qe\n";
    for (int i = 0; i < Nstep; ++i) {
        double B(i*Bmax/Nstep);
        double val(p_xsecs->GetSigmaSD()->GetCombinedValue(B));
        was<<B<<"     "<<val<<endl;
    }
    was<<"\n\n";
    was<<"# 4 - sigma_sd0\n";
    for (int i = 0; i < Nstep; ++i) {
        double B(i*Bmax/Nstep);
        double val(p_xsecs->GetSigmaSD()->GetCombinedValueSD0(B));
        was<<B<<"     "<<val<<endl;
    }
    was<<"\n\n";
    was<<"# 5 - sigma_sd1\n";
    for (int i = 0; i < Nstep; ++i) {
        double B(i*Bmax/Nstep);
        double val(p_xsecs->GetSigmaSD()->GetCombinedValueSD1(B));
        was<<B<<"     "<<val<<endl;
    }
    was<<"\n\n";
    was<<"# 6 - sigma_dd\n";
    for (int i = 0; i < Nstep; ++i) {
        double B(i*Bmax/Nstep);
        double val(p_xsecs->GetSigmaSD()->GetCombinedValueDD(B));
        was<<B<<"     "<<val<<endl;
    }
    was.close();
}


void Shrimps::WriteOutXSecsYodaFile(const std::set<double> & energies_tot,
                    const std::set<double> & energies_inel,
                    const std::set<double> & energies_el,
                    const std::set<double> & energies_sd,
                    const std::set<double> & energies_dd,
                    const std::vector<double> & xsectot,
				    const std::vector<double> & xsecinel,
				    const std::vector<double> & xsecelas,
                    const std::vector<double> & xsecsd,
                    const std::vector<double> & xsecdd,
                    std::string dirname) {
  std::string filename(dirname+std::string("/xsecs.dat"));
  std::ofstream was;
  was.open(filename.c_str());
  was<<"# BEGIN HISTO1D /XSECS/total\n";
  was<<"Path=/XSECS/total"<<std::endl;
  size_t i(0);
  for (std::set<double>::iterator energy_iter=energies_tot.begin();
       energy_iter!=energies_tot.end();energy_iter++) {
    was<<(*energy_iter)<<"   "<<(*energy_iter)<<"   "
       <<xsectot[i++]<<"   0.0   0.0\n";
  }
  was<<"# END HISTO1D\n"<<std::endl;
  was<<"# BEGIN HISTO1D /XSECS/inel\n";
  was<<"Path=/XSECS/inel"<<std::endl;
  i = 0;
  for (std::set<double>::iterator energy_iter=energies_inel.begin();
       energy_iter!=energies_inel.end();energy_iter++) {
    was<<(*energy_iter)<<"   "<<(*energy_iter)<<"   "
       <<xsecinel[i++]<<"   0.0   0.0\n";
  }
  was<<"# END HISTO1D\n"<<std::endl;
  was<<"# BEGIN HISTO1D /XSECS/el\n";
  was<<"Path=/XSECS/el"<<std::endl;
  i = 0;
  for (std::set<double>::iterator energy_iter=energies_el.begin();
       energy_iter!=energies_el.end();energy_iter++) {
    was<<(*energy_iter)<<"   "<<(*energy_iter)<<"   "
       <<xsecelas[i++]<<"   0.0   0.0\n";
  }
  was<<"# END HISTO1D"<<std::endl;
  was<<"# BEGIN HISTO1D /XSECS/sd\n";
  was<<"Path=/XSECS/sd"<<std::endl;
  i = 0;
  for (std::set<double>::iterator energy_iter=energies_sd.begin();
       energy_iter!=energies_sd.end();energy_iter++) {
    was<<(*energy_iter)<<"   "<<(*energy_iter)<<"   "
       <<xsecsd[i++]<<"   0.0   0.0\n";
  }
  was<<"# END HISTO1D"<<std::endl;
  was<<"# BEGIN HISTO1D /XSECS/dd\n";
  was<<"Path=/XSECS/dd"<<std::endl;
  i = 0;
  for (std::set<double>::iterator energy_iter=energies_dd.begin();
       energy_iter!=energies_dd.end();energy_iter++) {
    was<<(*energy_iter)<<"   "<<(*energy_iter)<<"   "
       <<xsecdd[i++]<<"   0.0   0.0\n";
  }
  was<<"# END HISTO1D"<<std::endl;
  was.close();
}
  
void Shrimps::ReadEnergiesFromFile(std::set<double> & energies,
				   std::string infile) {
  std::ifstream input;
  input.open(infile.c_str());
  if(!input){
    msg_Error()<<"File "<<infile<<" does not exist, will exit now.\n";
    exit(1);
  }
  std::string test;
  while (!input.eof()) {
    input>>test;
    energies.insert(std::atof(test.c_str()));
  }
  input.close();
}



void Shrimps::TestShrimps(PDF::ISR_Handler *const isr) {
  msg_Info()<<"Start testing SHRiMPS.\n";
  std::string dirname = std::string("Tests");
  ATOOLS::MakeDir(dirname);
  InitialiseFormFactors();
  InitialiseRemnants(isr);
  InitialiseSingleChannelEikonals();

  PrintAlphaS(dirname);
  PrintPDFs(dirname);
  MBpars.GetFormFactors()->front()->Test(dirname); 
  TestEikonalGrids(dirname);
  TestCrossSections(dirname);
  TestEventGeneration(dirname);
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
    Continued_PDF * pdf = p_remnants->GetHadronDissociation(0)->GetPDF();
    for (int j=0;j<=nxval; j++){
      x = pow(10.,double(j)/double(nxval)*log10(xmin));
      pdf->Calculate(x,Q2);
      was<<x<<"   "
	 <<pdf->XPDF(ATOOLS::Flavour(kf_u))<<"   "
	 <<pdf->XPDF(ATOOLS::Flavour(kf_u).Bar())<<"   "
	 <<pdf->XPDF(ATOOLS::Flavour(kf_d))<<"   "
	 <<pdf->XPDF(ATOOLS::Flavour(kf_d).Bar())<<"   "
	 <<pdf->XPDF(ATOOLS::Flavour(kf_s))<<"   "
	 <<pdf->XPDF(ATOOLS::Flavour(kf_gluon))<<"\n";
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
  Omega_ik * eikonal((*MBpars.GetEikonals())[0][0]);
  eikonal->TestIndividualGrids(&ana12,&ana21,Ymax,dirname);

  Analytic_Eikonal anaeik;
  eikonal->TestEikonal(&anaeik,dirname);
}

void Shrimps::TestCrossSections(const std::string & dirname) {
  Cross_Sections cross;
  cross.CalculateCrossSections();
  cross.Test(dirname);
}

void Shrimps::TestEventGeneration(const std::string & dirname) {
  Event_Generator generator(p_xsecs,true);
  generator.Test(dirname);
}



