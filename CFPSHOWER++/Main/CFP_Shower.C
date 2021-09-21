#include "CFPSHOWER++/Main/CFP_Shower.H"
#include "CFPSHOWER++/Tools/CFP_Parameters.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"

using namespace CFPSHOWER;
using namespace PDF;
using namespace ATOOLS;
using namespace std;

CFP_Shower::CFP_Shower(const Shower_Key &key) :
  Shower_Base("CFP_Shower"), p_massselector(NULL), m_test(false)
{
  cfp_pars = new CFP_Parameters();
  cfp_pars->Init(key.p_reader);
  m_shower.Init(key.p_model, key.p_isr);
  Parton::Reset(0);
  p_cluster = m_shower.GetClusterDefinitions();
  if (m_test) TestShower();
}
  
CFP_Shower::~CFP_Shower() {
  delete cfp_pars;
}

bool CFP_Shower::
PrepareShower(Cluster_Amplitude * const ampl,const bool & soft) {
  p_massselector = ampl->MS();
  m_shower.SetMassSelector(p_massselector);
  //Configuration * config = MakeSimpleFFConfiguration();
  //m_configs.push_back(config);
  //return true;
  std::map<Cluster_Leg *,Parton *> lmap;
  Cluster_Amplitude * current = InitConfigs(ampl);
  // Starting from the hardest part of the history, work backwards to the full
  // parton configuration.  In each step find the splitter-spectator pair
  // and the three resulting partons and fill an ME splitting - their
  // sequence represents the parton shower history of the ME event
  while (current) {
    Configuration * config = new Configuration(current,lmap,p_cluster);
    m_configs.push_back(config);
    current=current->Prev();
  }
  m_configs.front()->SetT(m_muQ2);
  //msg_Out()<<METHOD<<":\n"<<(*ampl)<<"\n"<<(*m_configs.front())<<"\n";
  return true;
}

Cluster_Amplitude * CFP_Shower::InitConfigs(Cluster_Amplitude * const ampl) {
  // Generic settings for the parton masses in the parton shower -
  // TODO: maybe we should add a security measure here?
  p_massselector = ampl->MS();
  m_shower.SetMassSelector(p_massselector);
  // Find the innermost, hardest cluster amplitude and extract the overall
  // parton shower starting scale
  Cluster_Amplitude * current = ampl;
  while (current->Next()) current = current->Next();
  m_muQ2 = current->MuQ2();
  return current;
}

Configuration * CFP_Shower::Convert(Cluster_Amplitude * const ampl,
				    map<Cluster_Leg*,Parton*> & lmap)
{
  Configuration * config = new Configuration(ampl,lmap,p_cluster);
  return config;
}

int CFP_Shower::PerformShowers() {
  //msg_Out()<<"#########################################################\n"
  //	   <<"#########################################################\n"
  //	   <<"#########################################################\n"
  //	   <<METHOD<<":\n";
  m_weight = 1.;
  for (size_t i=0;i<m_configs.size();i++) {
    if (m_shower.Evolve(m_configs[i])) m_weight *= m_shower.Weight();
  }
  m_shower.Reset();
  return 1;
}

int CFP_Shower::PerformDecayShowers() {
  return PerformShowers();
}

bool CFP_Shower::ExtractPartons(Blob_List *const bl) {
  //msg_Out()<<"############################################################\n"
  //	   <<"############################################################\n"
  //	   <<METHOD<<"\n"
    //	   <<(*bl)
  //	   <<"############################################################\n";
  Blob * blob(bl->FindLast(btp::Shower));
  if (blob==NULL) THROW(fatal_error,"No Shower blob");
  blob->SetTypeSpec("CFP_SHOWER");
  for (int i=0;i<blob->NInP();++i)
    blob->InParticle(i)->SetStatus(part_status::decayed);
  for (int i=0;i<blob->NOutP();++i)
    blob->OutParticle(i)->SetStatus(part_status::decayed);
  blob->SetStatus(blob_status::needs_beams |
		  blob_status::needs_hadronization);
  bool nois(blob->NOutP()==0);
  //if (nois) msg_Out()<<(*blob)<<"\n";
  while (!m_configs.empty()) {
    Configuration * config = m_configs.back();
    //msg_Out()<<"Extract partons from configuration "
    //	     <<"("<<m_configs.size()<<") with "
    //	     <<config->size()<<"("<<Parton::Count()<<") partons:\n"
    //	     <<(*config)<<"\n";
    while (!config->empty()) {
      Parton * parton = config->back();
      //msg_Out()<<config->size()<<" ("<<Parton::Count()<<") left: "
      //       <<"["<<parton->On()<<"]"<<(*parton);
      if (!(parton->Beam()>0 && nois) && parton->On())
	ExtractParton(blob,parton);
      delete parton;
      config->pop_back();
      //msg_Out()<<config->size()<<" ("<<Parton::Count()<<") left.\n";
    }
    delete config;
    m_configs.pop_back();
  }
  //msg_Out()<<"############################################################\n"
  //	   <<(*blob)<<"\n"
  //	   <<"finished "<<METHOD<<" with "<<Parton::Count()<<" partons "
  //	   <<"and "<<Splitting::Count()<<" splittings undeleted.\n"
  //	   <<"############################################################\n"
  //	   <<"############################################################\n";
  return true;
}

void CFP_Shower::ExtractParton(Blob *const blob,Parton *const parton) {
  //msg_Out()<<" --> "<<parton->Flav()<<": "<<parton->Mom()<<" "
  //	   <<"("<<parton->Beam()<<")\n";
  Particle * part = (parton->Beam()>0?
		     new Particle(-1,parton->Flav().Bar(),-parton->Mom(),'I'):
		     new Particle(-1,parton->Flav(),parton->Mom(),'F'));
  //msg_Out()<<" --> "<<part<<"\n";
  part->SetNumber(0);
  part->SetFinalMass(p_massselector->Mass(parton->Flav()));
  if (parton->Beam()==0) {
    //msg_Out()<<"Colour:\n";
    //msg_Out()<<"   ["<<parton->GetColor()[0]<<", "
    //<<parton->GetColor()[1]<<"]\n";
    part->SetFlow(1,parton->GetColor()[0]);
    part->SetFlow(2,parton->GetColor()[1]);
    blob->AddToOutParticles(part);
  }
  else {
    part->SetFlow(1,parton->GetColor()[1]);
    part->SetFlow(2,parton->GetColor()[0]);
    part->SetBeam(parton->Beam()-1);
    blob->AddToInParticles(part);
  }
  //msg_Out()<<" --> "<<(*part)<<"\n";
}

void CFP_Shower::CleanUp() {
  Parton::Reset(0);
  m_shower.Reset();
}

PDF::Cluster_Definitions_Base * CFP_Shower::GetClusterDefinitions() {}

void CFP_Shower::TestShower(const double & E,const Flavour & flav) {
  msg_Info()<<"In "<<METHOD
	    <<" will start trivial testing for FF configuration:\n";
  for (size_t i=0;i<1000;i++) {
    Configuration * config = MakeSimpleFFConfiguration(E,flav);
    m_shower.Evolve(config);
    m_shower.Reset();
    delete config;
  }
  msg_Info()<<"############################################################\n"
	    <<"############################################################\n"
	    <<"############################################################\n"
	    <<"finished "<<METHOD<<" with "<<Parton::Count()<<" partons "
	    <<"and "<<Splitting::Count()<<" splittings undeleted.\n"
	    <<"############################################################\n"
	    <<"############################################################\n"
	    <<"############################################################\n";
  exit(1);
}

Configuration * CFP_Shower::MakeSimpleFFConfiguration(const double & E,
						      const Flavour & flav) {
  double costheta  = 1.-2.*ran->Get(), sintheta = sqrt(1.-sqr(costheta));
  double phi       = 2.*M_PI*ran->Get();
  Vec3D  mom       = Vec3D(sintheta*cos(phi), sintheta*sin(phi), costheta);
  Configuration * config = new Configuration(E*E,1.);
  Flavour electron       = Flavour(kf_e);
  Parton * partonep      = new Parton(electron.Bar(),-E*Vec4D(1.,0.,0.,1.));
  partonep->SetBeam(1);
  Parton * partonem      = new Parton(electron,-E*Vec4D(1.,0.,0.,-1.));
  partonem->SetBeam(2);
  config->AddParton(partonem);
  config->AddParton(partonep);
  unsigned int color     = Flow::Counter();
  //ATOOLS::Vec4D mom1 = Vec4D(45.6,-9.86354, 34.7153,-27.8732);
  //ATOOLS::Vec4D mom2 = Vec4D(45.6, 9.86354,-34.7153, 27.8732);
  Parton * parton1       = new Parton(flav,      E*Vec4D(1,mom));
  Parton * parton2       = new Parton(flav.Bar(),E*Vec4D(1,-mom));
  parton1->SetColor(Color(color,0));
  parton1->AddSpectator(parton2);
  parton2->SetColor(Color(0,color));
  parton2->AddSpectator(parton1);
  config->AddParton(parton1);
  config->AddParton(parton2);
  config->SetT(sqr(2.*E));
  return config;
}


DECLARE_GETTER(CFP_Shower,"CFP",Shower_Base,Shower_Key);

Shower_Base * Getter<Shower_Base,Shower_Key,CFP_Shower>::
operator()(const Shower_Key & key) const
{
  return new CFP_Shower(key);
}

void Getter<Shower_Base,Shower_Key,CFP_Shower>::
PrintInfo(std::ostream & str,const size_t width) const
{ 
  str<<"The CFP Parton Shower"; 
}
