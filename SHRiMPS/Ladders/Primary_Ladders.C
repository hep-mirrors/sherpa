#include "SHRiMPS/Ladders/Primary_Ladders.H"
#include "SHRiMPS/Ladders/Ladder_Generator_QT.H"
#include "SHRiMPS/Ladders/Ladder_Generator_KT.H"
#include "SHRiMPS/Ladders/Ladder_Generator_Eik.H"
#include "SHRiMPS/Tools/MinBias_Parameters.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace SHRIMPS;
using namespace ATOOLS;

Primary_Ladders::Primary_Ladders() :
  p_laddergenerator(new Ladder_Generator_Eik()),
  m_Ecms(rpa->gen.Ecms()/2.),
  m_test(false),
  n_calls(0), n_start(0), n_gen(0)
{
  if (m_test) {
    m_histos[std::string("N_start")]             = new Histogram(0, -0.5, 19.5, 20);
    m_histos[std::string("N_gen")]               = new Histogram(0, -0.5, 19.5, 20);
    m_histos[std::string("n_trial")]             = new Histogram(0, -0.5, 19.5, 20);
    m_histos[std::string("n_accept")]            = new Histogram(0, -0.5, 19.5, 20);
    m_histos[std::string("n_trial_highpt")]      = new Histogram(0, -0.5,  9.5, 10);
    m_histos[std::string("n_accept_highpt")]     = new Histogram(0, -0.5,  9.5, 10);
    m_histos[std::string("Yasym_trial")]         = new Histogram(0,  0.0,  8.0, 32);
    m_histos[std::string("Yasym_trial_pt")]      = new Histogram(0,  0.0,  8.0, 32);
    m_histos[std::string("KTasym_trial")]        = new Histogram(0,  0.0,  8.0, 32);
    m_histos[std::string("Yasym_accept")]        = new Histogram(0,  0.0,  8.0, 32);
    m_histos[std::string("Yasym_accept_pt")]     = new Histogram(0,  0.0,  8.0, 32);
    m_histos[std::string("KTasym_accept")]       = new Histogram(0,  0.0,  8.0, 32);
    m_histos[std::string("Yasym_trial_highpt")]  = new Histogram(0,  0.0,  8.0, 32);
    m_histos[std::string("Yasym_accept_highpt")] = new Histogram(0,  0.0,  8.0, 32);
  }
}

Primary_Ladders::~Primary_Ladders() {
  Reset();
  if (p_laddergenerator) delete p_laddergenerator;
  if (m_test) {
    std::string name  = std::string("Ladder_Analysis/");
    for (std::map<std::string, Histogram * >::iterator hit=m_histos.begin();
	 hit!=m_histos.end();hit++) {
      hit->second->Finalize();
      hit->second->Output(name+hit->first);
      delete hit->second;
    }
    msg_Out()<<METHOD<<" produced "<<n_gen<<" ladders (start = "<<n_start<<") in "
	     <<n_calls<<" events.\n";
  }
}

void Primary_Ladders::Initialise(Remnant_Handler * remnants) {
  p_laddergenerator->Initialise(remnants);
}

void Primary_Ladders::Test() { return; if (m_test) p_laddergenerator->Test(); }

bool Primary_Ladders::operator()(Omega_ik * eikonal,const double & B,const size_t & N) {
  Reset();
  //msg_Out()<<"     -------------------------------------------------------------\n"
  //	   <<"     --- Make "<<N<<" new ladders at B = "<<B<<"\n";
  p_laddergenerator->InitCollision(eikonal,B);
  size_t Ngen = 0, trials = 0;
  double b1, b2;
  bool   contains_one_inelastic = false;
  msg_Out()<<"--------------------------------------------------------------\n";
  while (Ngen<N) {
    
    Vec4D position = eikonal->SelectB1B2(b1,b2,B);
    p_laddergenerator->SetImpactParameters(b1,b2);
    p_laddergenerator->SetMaximalScale(m_E[0],m_E[1]);
    Ladder * ladder = (*p_laddergenerator)(position);
    if (m_test && ladder) FillAnalysis(ladder,"trial");
    if (IsAllowed(ladder) && m_colourgenerator(ladder)) {	
      p_laddergenerator->QuarkReplace();
      p_laddergenerator->FixLadderType();
      if (ladder->Type()==ladder_type::inelastic) contains_one_inelastic = true;
      Add(ladder);
      Ngen++;
      trials = 0;
      if (m_test) FillAnalysis(ladder,"accept");
    }
    else {
      if (ladder) delete ladder;
      if (Ngen>0 && (trials++)>100) break;
    }
  }
  if (m_test) {
    n_gen += Ngen; n_start += N; n_calls++;
    m_histos[std::string("N_start")]->Insert(N);
    m_histos[std::string("N_gen")]->Insert(Ngen);
  }
  msg_Out()<<"contains one inelastic ladder = "<<contains_one_inelastic<<"\n"
	   <<"--------------------------------------------------------------\n";
  return contains_one_inelastic;
}
 
void Primary_Ladders::Reset() {
  m_E[0] = m_E[1] = m_Ecms;
  while (!m_ladders.empty()) {
    delete (m_ladders.back());
    m_ladders.pop_back();
  }
  m_ladders.clear();
  m_colourgenerator.Reset();
}

bool Primary_Ladders::IsAllowed(Ladder * ladder) {
  if (ladder==NULL) return false;
  for (size_t i=0;i<2;i++) {
    if (m_E[i]-ladder->InPart(i)->Momentum()[0] < 5.) return false;
  }
  return true;
}

void Primary_Ladders::Add(Ladder * ladder) {
  for (size_t i=0;i<2;i++) m_E[i] -= ladder->InPart(i)->Momentum()[0];
  m_ladders.push_back(ladder);
}
 
void Primary_Ladders::FillAnalysis(Ladder * ladder,const std::string & tag) {
  m_histos[std::string("n_")+tag]->Insert(ladder->GetEmissions()->size());
  for (LadderMap::iterator pit=ladder->GetEmissions()->begin(); pit!=ladder->GetEmissions()->end();pit++) {
    double kt = pit->second.Momentum().PPerp();
    m_histos[std::string("KTasym_")+tag]->Insert(dabs(pit->first), pit->first>0.?kt:-kt);
    m_histos[std::string("Yasym_")+tag]->Insert(dabs(pit->first), pit->first>0.?1.:-1.);
    m_histos[std::string("Yasym_")+tag]->Insert(dabs(pit->first), pit->first>0.?1.:-1.);
    m_histos[std::string("Yasym_")+tag+std::string("_pt")]->Insert(dabs(pit->first),pit->first>0.?1.:-1.);
    if (pit->second.Momentum().PPerp()>2.5)
      m_histos[std::string("Yasym_")+tag+std::string("_highpt")]->
	Insert(dabs(pit->first),pit->first>0.?1.:-1.);
  }
}
