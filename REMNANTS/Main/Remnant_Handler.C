#include "REMNANTS/Main/Remnant_Handler.H"
#include "REMNANTS/Main/Hadron_Remnant.H"
#include "REMNANTS/Main/Electron_Remnant.H"
#include "REMNANTS/Main/Photon_Remnant.H"
#include "REMNANTS/Main/No_Remnant.H"
#include "PDF/Main/ISR_Handler.H"
#include "BEAM/Main/Beam_Spectra_Handler.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"

using namespace REMNANTS;
using namespace ATOOLS;
using namespace std;

Remnant_Handler::
Remnant_Handler(PDF::ISR_Handler * isr,BEAM::Beam_Spectra_Handler * beam,
		const std::string & path,const std::string & file) :
  m_check(true), m_output(true)
{
  InitializeRemnants(isr,beam);
  DefineRemnantStrategy();
  InitializeKinematicsAndColours(path,file);
}

Remnant_Handler::~Remnant_Handler() {
  for (size_t i(0);i<2;++i) { if (p_remnants[i]!=NULL) delete p_remnants[i]; }
}

void Remnant_Handler::
InitializeRemnants(PDF::ISR_Handler * isr,BEAM::Beam_Spectra_Handler * beam) {
  for (size_t i=0;i<2;++i) {
    p_remnants[i] = NULL;
    Flavour flav  = isr->Flav(i);
    if (isr->PDF(i)!=0) {
      if (flav.IsHadron())
	p_remnants[i] = new Hadron_Remnant(isr->PDF(i),i);
      else if (flav.IsLepton())
	p_remnants[i] = new Electron_Remnant(isr->PDF(i),i);
      else if (flav.IsPhoton()) {
	// TODO: This is a bit more tricky once we assume a photon with hadronic structure.
	// Without a structure we do not need a remnant and will get away with No_Remnant.
	msg_Error()<<METHOD<<": Photon remnants not implemented yet.\n"
		   <<"   Will continue and assume point-particles.\n";
      }
    }
    if (p_remnants[i]==NULL) p_remnants[i] = new No_Remnant(i);
  }
  // Finish the initialisation of the Remnant_Bases: make sure they know each other,
  // their beam, the Colour_Generator, and hand them also to the ISR_Handler.
  // TODO: this latter part may become obsolete - I will have to check this.
  for (size_t i=0;i<2;++i) {
    p_remnants[i]->SetPartner(p_remnants[1-i]);
    p_remnants[i]->SetBeam(beam->GetBeam(i));
    p_remnants[i]->SetColours(&m_colours);
    p_remnants[i]->Reset();
    isr->SetRemnant(p_remnants[i],i);
  }
}

void Remnant_Handler::DefineRemnantStrategy() {
  // Pretty self-explanatory.  This covers the way of how we make sure that potential intrinsic
  // transverse momenta in the beam breakups do not lead to violation of four-momentum
  // conservation.  We have pretty much 4 options:
  // - simple, where the remnants do not break up;
  // - ll where the breakup is collinear only and we only need to boost the system in the end;
  // - DIS, where one breakup is collinear and the other involves transverse momenta; and
  //   (TODO: this is the one I still need to implement)
  // - hh, where both breakups generate transverse momenta.
  // For the latter two we will realize four-momentum conservation through insertion of a
  // "soft" blob, mainly a garbage collection where we collect partiles and shuffle them
  // in a pretty minimal fashion.
  if (p_remnants[0]->Type()==rtp::intact && p_remnants[1]->Type()==rtp::intact)
    m_type = strat::simple;
  else if (p_remnants[0]->Type()==rtp::lepton && p_remnants[1]->Type()==rtp::lepton)
    m_type = strat::ll;
  else if (p_remnants[0]->Type()==rtp::hadron && p_remnants[1]->Type()==rtp::lepton)
    m_type = strat::DIS1;
  else if (p_remnants[0]->Type()==rtp::lepton && p_remnants[1]->Type()==rtp::hadron)
    m_type = strat::DIS2;
  else if (p_remnants[0]->Type()==rtp::hadron && p_remnants[1]->Type()==rtp::hadron)
    m_type = strat::hh;
  else {
    msg_Error()<<METHOD<<" throws error: no strategy found for remnants "
	       <<int(p_remnants[0]->Type())<<" & "<<int(p_remnants[1]->Type())<<"\n"
	       <<"   Will exit the run.\n";
    exit(1);
  }
}

void Remnant_Handler::
InitializeKinematicsAndColours(const std::string & path,const std::string & file) {
  m_kinematics.Initialize(this,path,file);
  m_colours.Initialize(this);
}

bool Remnant_Handler::ExtractShowerInitiators(Blob *const showerblob) {
  // This method is called after each successful parton shower off a hard scatter.
  // It extracts the two initial particles from the shower and extracts them
  // from the corresponding beam remnant.
  
  // Make sure only shower blobs with exactly two initiators are treated, and only once.
  if (!(showerblob->Type()==btp::Shower) ||
      m_treatedshowerblobs.find(showerblob)!=m_treatedshowerblobs.end()) return true;
  size_t countIn = 0;
  for (size_t i=0;i<showerblob->NInP();++i) {
    if (!showerblob->InParticle(i)->ProductionBlob()) countIn++;
  }
  if (countIn!=2) return true;
  // Now extract the shower initiators from the remnants - they will get added
  // to the lists of extracted particles for each remnant and their colour will
  // be added to the Colour_Generator in each beam.
  for (size_t i=0;i<showerblob->NInP();++i) {
    Particle * part = showerblob->InParticle(i);
    if (part->ProductionBlob()!=NULL) continue;
    // Make sure extraction works out - mainly subject to energy conservation
    if (!Extract(part,part->Beam())) return false;
  }
  m_treatedshowerblobs.insert(showerblob);
  return true;
}

void Remnant_Handler::ConnectColours(ATOOLS::Blob *const showerblob) {
  // After each showering step, we try to compensate some of the colours.  In the absence
  // of multiple parton interactions this will not involve anything complicated with colours.
  // In each step, the shower initiators are checked for being a valence quark.  In case they
  // are, remnants (equivalent to recoilers) are generated, usually diquarks for protons,
  // if they are seaquarks, suitable spectators are generated.  The colours of the shower
  // initiators and, possibly, spectators will be added to a stack which will in turn partially
  // replace the new colours.  This is handled in the Colour_Generator.

  //bool printit = false;
  //Vec4D check(0.,0.,0.,0.);
  //for (size_t i=0;i<showerblob->NOutP();i++) {
  //  Particle * part = showerblob->OutParticle(i);
  //  if (part->Info()=='F') check += part->Momentum();
  //  if (part->Momentum()[0]<0.) printit = true;
  //}
  //for (size_t i=0;i<showerblob->NInP();i++) {
  //  Particle * part = showerblob->InParticle(i);
  //  if (part->Momentum()[0]<0.) printit = true;
  //}
  //if (printit || check.PPerp()>1.) {
  //  msg_Out()<<"#########################################################################\n"
  //  	     <<"#########################################################################\n"
  //  	     <<METHOD<<" for showerblob ("<<showerblob->Id()<<") with transverse momentum "
  //	     <<"kt = "<<check.PPerp()<<"\n"<<(*showerblob)<<"\n";
  //x}
  m_colours.ConnectColours(showerblob);
  //msg_Out()<<METHOD<<" done.\n"
  //	   <<"#########################################################################\n";
}

Return_Value::code Remnant_Handler::MakeBeamBlobs(Blob_List *const bloblist,
						  Particle_List *const particlelist)
{
  //msg_Out()<<"#########################################################################\n"
  //	   <<"Enter "<<METHOD<<".\n";
    
  // Adding the blobs related to the breakup of incident beams: one for each beam,
  // plus, potentially a third one to balance transverse momenta. 
  InitBeamAndSoftBlobs(bloblist);
  // Fill in the transverse momenta though the Kinematics_Generator.
  if (!m_kinematics.FillBlobs(bloblist) || !CheckBeamBreakup(bloblist)) {
    Reset();
    return Return_Value::Retry_Event;
  }
  Reset();
  return Return_Value::Success;
}

void Remnant_Handler::InitBeamAndSoftBlobs(Blob_List *const bloblist) {
  // Making a new blob (softblob) to locally compensate 4 momentum.  Ultimately,
  // it will reflect different strategies of how to compensate intrinsic kperp:
  // hadron colliders vs. DIS (DIS still needs to be implemented).
  if (m_type!=strat::simple && m_type!=strat::ll) {
    p_softblob = m_kinematics.MakeSoftBlob();
    bloblist->push_front(p_softblob);
  }  
  // Look for shower blobs that need beams and unset the flag
  for (Blob_List::iterator bit=bloblist->begin();bit!=bloblist->end();++bit) {
    if ((*bit)->Has(blob_status::needs_beams) && (*bit)->Type()==btp::Shower) {
      (*bit)->UnsetStatus(blob_status::needs_beams);
    }
  }
  // Remnant bases will generate their beam blobs, reset the incoming four-momenta
  m_colours.ResetFlags();
  for (size_t beam=0;beam<2;beam++) {
    bloblist->push_front(p_remnants[beam]->MakeBlob());
  }
}

bool Remnant_Handler::CheckBeamBreakup(Blob_List * bloblist) {
  // Final checks on beam breakup: four-momentum and colour conservation
  if (m_type==strat::simple) return 1;
  if (!m_check || (bloblist->FourMomentumConservation() &&
		   bloblist->ColorConservation())) return true;
  if (m_output) msg_Error()<<"Error in "<<METHOD<<": "
			   <<"colour or four-momentum not conserved.\n";
  bool mom(false), col(false);
  for (Blob_List::iterator bit=bloblist->begin();
       bit!=bloblist->end();++bit) {
    Vec4D checkmom = (*bit)->CheckMomentumConservation();
    if (dabs(checkmom.Abs2())>0.01 || dabs(checkmom[0])>0.1) {
      if (m_output) msg_Error()<<"   momentum non-conservation ("<<checkmom<<") in\n"
			       <<(**bit)<<"\n";
      mom = true;
    }
  }
  for (Blob_List::iterator bit=bloblist->begin();
       bit!=bloblist->end();++bit) {
    bool transient = ((*bit)->Type()==btp::Signal_Process ||
		      ((*bit)->Type()==btp::Shower &&
		       (*bit)->InParticle(0)->ProductionBlob()->Type()==btp::Signal_Process));
    if (!(*bit)->CheckColour(transient)) {
      if (m_output) msg_Error()<<"   colour non-conservation in\n"
			       <<(**bit)<<"\n";
      col = true;
    }
  }
  if (col) exit(1);
  return false;
}

bool Remnant_Handler::Extract(ATOOLS::Particle * part,const unsigned int beam) {
  // Extracting a particle from a remnant only works for positive energies.
  if (part->Momentum()[0]<0.) {
    msg_Error()<<METHOD<<" yields shower with negative incoming energies.\n"
	       <<(*part->DecayBlob())<<"\n";
    return false;
  }
  return p_remnants[beam]->Extract(part);
}

void Remnant_Handler::Reset() {
  for (size_t beam=0;beam<2;beam++) p_remnants[beam]->Reset();
  m_treatedshowerblobs.clear();
  m_kinematics.Reset();
  m_colours.Reset();
}
