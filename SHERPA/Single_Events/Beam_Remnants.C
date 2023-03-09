#include "SHERPA/Single_Events/Beam_Remnants.H"
#include "ATOOLS/Org/Message.H"

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

Beam_Remnants::Beam_Remnants(Beam_Remnant_Handler * _beamremnant) :
  p_beamremnanthandler(_beamremnant),
  m_ana(false)
{
  m_name = "Beam_Remnants:"+(p_beamremnanthandler->Fill()==1?
			     p_beamremnanthandler->Name():string("None"));
  m_type = eph::Hadronization;
  if (m_ana) InitHistos();
}

Beam_Remnants::~Beam_Remnants() {
  if (m_ana && !m_histos.empty()) {
    for (map<string, Histogram * >::iterator hit=m_histos.begin();
	 hit!=m_histos.end();hit++) {
      string name  = string("ImpactParameter_Analysis/");
      hit->second->Finalize();
      hit->second->Output(name+hit->first+string(".dat"));
      delete hit->second;
    }   
  }
}

Return_Value::code Beam_Remnants::Treat(ATOOLS::Blob_List* bloblist)
{
  if (bloblist->empty()) {
    msg_Error()<<"Beam_Remnants::Treat("<<bloblist<<"): "<<endl
	       <<"   Blob list contains "<<bloblist->size()<<" entries."<<endl
	       <<"   Continue and hope for the best."<<endl;
    return Return_Value::Error;
  }
  Blob *signal(bloblist->FindFirst(btp::Signal_Process));
  if (signal && signal->NInP()<2) return Return_Value::Nothing;
  bool onlyBunch = false;
  if (!signal || signal->Has(blob_status::needs_signal)) {
    Blob * hard  = bloblist->FindFirst(btp::Hard_Collision);
    Blob * qelas = bloblist->FindFirst(btp::Elastic_Collision);
    if (!qelas) qelas = bloblist->FindFirst(btp::Soft_Diffractive_Collision);        
    if (!qelas) qelas = bloblist->FindFirst(btp::Quasi_Elastic_Collision);        
    if (!hard && !qelas) return Return_Value::Nothing;
    if (qelas) onlyBunch = true;
  }
  else {
    btp::code signal_type = signal->Type();
    if (signal_type==btp::Elastic_Collision ||
	signal_type==btp::Soft_Diffractive_Collision ||
	signal_type==btp::Quasi_Elastic_Collision)
      onlyBunch = true;
  }
  Blob * beam(bloblist->FindFirst(btp::Beam));
  if (beam && !beam->Has(blob_status::needs_beams)) return Return_Value::Nothing;
  if (m_ana) Analyse(bloblist);
  return p_beamremnanthandler->FillBeamAndBunchBlobs(bloblist,onlyBunch);
}


void Beam_Remnants::CleanUp(const size_t & mode) 
{
  p_beamremnanthandler->CleanUp(mode);
}

void Beam_Remnants::Finish(const std::string &) {}

void Beam_Remnants::InitHistos() {
  m_histos["phi_1"] = new Histogram(0,0.,360.,36);
  m_histos["phi_2"] = new Histogram(0,0.,360.,36);
  m_histos["b_1"]   = new Histogram(0,0.,10.,100);
  m_histos["b_2"]   = new Histogram(0,0.,10.,100);
  m_histos["B"]     = new Histogram(0,0.,20.,200);
}

void Beam_Remnants::Analyse(ATOOLS::Blob_List * blobs) {
  Vec4D B1 = p_beamremnanthandler->GetRemnants()->GetRemnant(0)->Position(), B2 = -B1;
  m_histos["B"]->Insert(2.*B1.PPerp());
  for (deque<Blob *>::iterator bit=blobs->begin();bit!=blobs->end();bit++) {
    Blob * blob = *bit;
    if (blob->Type()==btp::Hard_Collision) {
      Vec4D pos   = blob->Position();
      double b1   = (pos-B1).PPerp(), b2 = (pos-B2).PPerp();
      double phi1 = acos((pos[1]-B1[1])/b1) * 360./(2.*M_PI);
      double phi2 = acos((pos[1]-B2[1])/b2) * 360./(2.*M_PI);
      m_histos["b_1"]->Insert(b1);
      m_histos["b_2"]->Insert(b2);
      m_histos["phi_1"]->Insert(phi1);
      m_histos["phi_2"]->Insert(phi2);
    }
  }
}

