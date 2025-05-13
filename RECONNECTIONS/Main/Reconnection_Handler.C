#include "RECONNECTIONS/Main/Reconnection_Handler.H"
#include "RECONNECTIONS/Main/Reconnect_By_Singlet.H"
#include "RECONNECTIONS/Main/Reconnect_Statistical.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Scoped_Settings.H"

using namespace RECONNECTIONS;
using namespace ATOOLS;
using namespace std;

Reconnection_Handler::Reconnection_Handler(const bool & on) :
  m_on(on),
  p_reconnector(new Reconnect_Statistical()),
  m_nfails(0)
{
  const auto wgts {p_reconnector->get_variationweights()};

  // TODO fix
  wgt_sums.resize(10,0);
}

Reconnection_Handler::~Reconnection_Handler() {
  if (m_on) {
    msg_Info()<<METHOD<<": reconnection handler winds down with "
	      <<m_nfails<<" errors overall.\n";
  }
  delete p_reconnector;
}

void Reconnection_Handler::Initialize() {
  if (m_on) p_reconnector->Initialize();
}

void Reconnection_Handler::Reset() {
  if (m_on) p_reconnector->Reset();
}

Return_Value::code Reconnection_Handler::operator()(Blob_List *const blobs,
						    Particle_List *const parts) {
  if (!m_on) return Return_Value::Nothing;
  switch ((*p_reconnector)(blobs)) {
  case -1:
    // things went wrong, try new event and hope it works better
    msg_Tracking()<<"Error in "<<METHOD<<": reconnections didn't work out.\n"
		  <<"   Ask for new event and hope for the best.\n";
    p_reconnector->Reset();
    m_nfails++;
    return Return_Value::New_Event;
  case 1:
    // added colour reconnections, but produce a reconnection blob.
    AddReconnectionBlob(blobs);
  case 0:
    // didn't find any blob that needed reconnections
    break;
  }
  p_reconnector->Reset();

  // add weights
  const auto wgts_cr  = p_reconnector->get_variationweights();
  Blob *blob(blobs->FindFirst(btp::Signal_Process));
  if (blob == NULL)
    blob = blobs->FindFirst(btp::Hard_Collision);
  auto & wgtmap = (*blob)["WeightsMap"]->Get<Weights_Map>();
  //const bool found {wgtmap.find("RECONNECTIONS") == wgtmap.end() ? false : true};
  //std::cout << " ---> " << wgts_cr << std::endl;
  m_nreconnections++;
  for(int i{0}; i<wgts_cr.size(); i++) {
    const std::string name {"v"+std::to_string(i)};
    const auto wgt {wgts_cr[i]};
    wgt_sums[i] += wgt;

    double norm {1.};
    if(m_nreconnections > 100) {
      norm = wgt_sums[i] / m_nreconnections;
    }
    
    const auto clipped_wgt {std::min(wgt,pow(10,2))};
    wgtmap["RECONNECTIONS"][name] = clipped_wgt / norm;
  }
  //std::cout << wgts_cr << std::endl;

  const auto n_vars {wgts_cr.size()};
  p_reconnector->reset_variationweights(n_vars);
  return Return_Value::Success; 
}

void Reconnection_Handler::AddReconnectionBlob(Blob_List *const blobs) {
  Blob * blob = new Blob();
  blob->AddStatus(blob_status::needs_hadronization);
  blob->SetType(btp::Fragmentation);
  blob->SetId();
  Part_List * particles = p_reconnector->GetParticles();
  while (!particles->empty()) {
    Particle * part = particles->front();
    part->DecayBlob()->AddToOutParticles(part);
    part->SetDecayBlob(NULL);
    blob->AddToInParticles(part);
    particles->pop_front();
  }
  blobs->push_back(blob);
}

