#include "SHRiMPS/Main/Cluster_Algorithm.H"
#include "SHRiMPS/Event_Generation/Inelastic_Event_Generator.H"
#include "SHRiMPS/Cross_Sections/Sigma_Inelastic.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/Histogram.H"

using namespace SHRIMPS;
using namespace ATOOLS;
using namespace std;

Inelastic_Event_Generator::
Inelastic_Event_Generator(Sigma_Inelastic * sigma,const int & test) :
  Event_Generator_Base(sigma), p_sigma(sigma), 
  m_primaries(Primary_Ladders()),
  m_mustinit(true),
  p_collemgen(new Collinear_Emission_Generator())
{}

Inelastic_Event_Generator::~Inelastic_Event_Generator() {
  Reset();
  while (!m_Bgrids.empty()) {
    delete m_Bgrids.begin()->second;
    m_Bgrids.erase(m_Bgrids.begin());
  }
  m_Bgrids.clear();
  delete p_collemgen;
}

void Inelastic_Event_Generator::
Initialise(Remnant_Handler * remnants,Cluster_Algorithm * cluster) {
  m_sigma = 0.;
  Sigma_Inelastic sigma;
  vector<vector<Omega_ik *> > * eikonals(MBpars.GetEikonals());
  for (size_t i=0;i<eikonals->size();i++) {
    for (size_t j=0;j<(*eikonals)[i].size();j++) {
      Omega_ik * eikonal = (*eikonals)[i][j];
      m_Bgrids[eikonal] = sigma.FillBGrid(eikonal);
      m_sigma += m_xsecs[eikonal] = m_Bgrids[eikonal]->back() * rpa->Picobarn();
    }
  }
  msg_Info()<<METHOD<<" yields effective inelastic cross section "
	    <<"sigma = "<<m_sigma/1.e9<<" mbarn.\n";
  p_cluster  = cluster;
  m_primaries.Initialise(remnants);
  //m_primaries.Test();
  Reset();
}

void Inelastic_Event_Generator::Reset() {
  m_mustinit = true;
  m_primaries.Reset();
}

int Inelastic_Event_Generator::GenerateEvent(Blob_List * blobs,const bool & flag) {
  //msg_Out()<<"######################################################################\n"
  //	   <<"######################################################################\n"
  //	   <<"######################################################################\n"
  //	   <<METHOD<<"\n";
  if (m_mustinit) {
    Blob * blob(blobs->FindFirst(btp::Soft_Collision));
    if (!blob ||
	!blob->Has(blob_status::needs_minBias)) return 0;
    if (!InitEvent(blobs) ||
	!m_primaries(p_eikonal,m_B,m_Nladders)) return -1;
  }
  return MakePrimaryScatterBlobs(blobs);
}

bool Inelastic_Event_Generator::MakePrimaryScatterBlobs(ATOOLS::Blob_List * blobs) {
  if (m_primaries.GetLadders()->empty()) return 0;
  Ladder * ladder = m_primaries.GetLadders()->front();
  Blob * blob     = NULL;
  if (m_mustinit) {
    blob = blobs->FindFirst(btp::Soft_Collision);
    blob->AddData("Weight",new Blob_Data<double>(m_sigma));
    blob->AddData("Factorisation_Scale",new Blob_Data<double>(1.));
    blob->AddData("Renormalization_Scale",new Blob_Data<double>(1.));
    m_mustinit = false;
  }
  else {
    blob = new Blob();
    blob->SetId();
    blobs->push_back(blob);
    blob->AddData("Weight",new Blob_Data<double>(1.));
  }
  blob->SetPosition(ladder->Position());
  blob->SetType(btp::Hard_Collision);
  blob->SetTypeSpec("MinBias");
  blob->UnsetStatus(blob_status::needs_minBias);
  blob->SetStatus(blob_status::needs_showers);
  for (size_t i=0;i<2;i++) blob->AddToInParticles(ladder->InPart(i)->GetParticle());
  for (LadderMap::iterator lmit=ladder->GetEmissions()->begin();
       lmit!=ladder->GetEmissions()->end();lmit++) {
    Particle * part = lmit->second.GetParticle();
    blob->AddToOutParticles(part);
    if (dabs(lmit->first)>m_primaries.Ymax()) part->SetInfo('B');
  }
  delete ladder;
  m_primaries.GetLadders()->pop_front();
  return p_collemgen->GenerateEmissions(blobs);
  //PRINT_VAR(*blob);
  //return 1;
}

bool Inelastic_Event_Generator::InitEvent(Blob_List * blobs) {
  p_eikonal  = 0; m_B = -1;
  for (size_t trials=0;trials<1000;trials++) {
    if (SelectEikonal() && SelectB()) {
      m_Nladders = int(1./p_eikonal->Prefactor() *
		       ran->Poissonian((*p_eikonal)(m_B)));
      if (m_Nladders>0) return true;
    }
  }
  return false;
}

bool Inelastic_Event_Generator::SelectEikonal() {
  p_eikonal = 0;
  while (p_eikonal==NULL) {
    double disc = ran->Get()*m_sigma;
    for (std::map<Omega_ik *,double>::iterator eikiter=m_xsecs.begin();
	 eikiter!=m_xsecs.end();eikiter++) {
      disc-=eikiter->second;
      if (disc<=1.e-12) { p_eikonal = eikiter->first; break; }
    }
  }
  return (p_eikonal!=0);
}

bool Inelastic_Event_Generator::SelectB() {
  if (p_eikonal==0) {
    msg_Error()<<"Error in "<<METHOD<<": no eikonal selected.\n";
    return false;
  }
  std::vector<double> * grid = m_Bgrids[p_eikonal];  
  double deltaB(p_eikonal->DeltaB());
  m_B = -1.;
  do {
    double random = ran->Get()*(*grid)[grid->size()-1];
    size_t bin(0);
    while (bin<grid->size()-1 && (random-(*grid)[bin]>=0)) bin++;
    if (bin>=grid->size()) continue;
    double inthigh((*grid)[bin]), intlow((*grid)[bin-1]);
    double Bhigh(bin*deltaB), Blow((bin-1)*deltaB);
    m_B  = (Blow*(random-intlow)+Bhigh*(inthigh-random))/(inthigh-intlow);
  } while (m_B<0.);
  return (m_B>=0.);
}


void Inelastic_Event_Generator::Test(const std::string & dirname) {}
