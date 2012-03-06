#include "SINIC++/Beam_Remnants/Beam_Remnant_Handler.H"
#include "SINIC++/Tools/MinBias_Parameters.H"
#include "BEAM/Main/Beam_Spectra_Handler.H"
#include "ATOOLS/Phys/Momenta_Stretcher.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Message.H"

using namespace SINIC;
using namespace ATOOLS;
using namespace std;

Beam_Remnant_Handler::
Beam_Remnant_Handler(BEAM::Beam_Spectra_Handler * beamspectra,
		     vector<Continued_PDF> & pdfs) :
  p_blob(NULL), p_pdfs(&pdfs), m_paircounter(0)
{
  m_checkmom.push_back(Vec4D(0,0,0,0)); 
  m_checkmom.push_back(Vec4D(0,0,0,0));
  for (int beam=0;beam<2;beam++) {
    m_beams.push_back(beamspectra->GetBeam(beam));
    double E(beamspectra->GetBeam(beam)->OutMomentum()[0]);
    int dir(beamspectra->GetBeam(beam)->OutMomentum()[3]>0?1:-1);
    m_beamvecs.push_back(Vec4D(E,0,0,dir*E));
    m_hadrons.push_back(new Hadron_Dissociation(&pdfs[beam]));
  }
  p_colour        = new Colour_Generator(&m_hadrons);
  p_reconnections = new Colour_Reconnections();
}

Beam_Remnant_Handler::~Beam_Remnant_Handler() {
  while (!m_hadrons.empty()) {
    delete m_hadrons.back();
    m_hadrons.pop_back();
  }
  delete p_colour;
  delete p_reconnections;
}

bool Beam_Remnant_Handler::
InitialiseCollision(const int & N,Omega_ik * eikonal) {
  Reset();
  if (eikonal==NULL && N==0) {
    for (size_t beam=0;beam<2;beam++) {
      m_hadrons[beam]->DeleteParticles();
    }
    return true;
  }
  double eta(eikonal->EffectiveIntercept());
  if (!m_hadrons[0]->DefineDissociation(N,0.0001,eta,eikonal->FF1()) ||
      !m_hadrons[1]->DefineDissociation(N,0.0001,eta,eikonal->FF2())) {
    for (size_t beam=0;beam<2;beam++) m_hadrons[beam]->DeleteParticles();
    return false;
  }
  bool ok;
  Flavour flav[2];
  do {  
    ok = true;
    for (int i=0;i<N;i++) {
      for (size_t beam=0;beam<2;beam++) {
	flav[beam] = m_hadrons[beam]->GetParticle(i)->Flav();
      }
      if (flav[0].IsQuark() && flav[1].IsQuark()) {
	if ((flav[0].IsAnti() && flav[1].IsAnti()) ||
	    (!flav[0].IsAnti() && !flav[1].IsAnti())) {
	  ok = false;
	  m_hadrons[int(ran->Get()>=0.5)]->Reshuffle(i);
	  break;
	}
      }
    }
  } while (!ok);

  p_blob = new Blob();
  p_blob->SetType(btp::Soft_Collision);
  p_blob->SetTypeSpec("Four_Momentum_Compensation");
  p_blob->SetId();
  p_blob->SetStatus(blob_status::needs_hadronization |
		    blob_status::needs_beams);
  for (size_t beam=0;beam<2;beam++) 
    m_hadrons[beam]->AddParticlesToBlob(p_blob,beam);
  p_colour->SetSoftBlob(p_blob);
  m_paircounter = 0;
  return true;
}

void Beam_Remnant_Handler::AddBeamBlobs(ATOOLS::Blob_List * blobs) {
  for (size_t beam=0;beam<2;beam++) {
    m_hadrons[beam]->FillBeamBlob();
    blobs->push_front(m_hadrons[beam]->GetBeamBlob());
  }
}

Return_Value::code Beam_Remnant_Handler::
FillBeamBlobs(Blob_List * blobs,Omega_ik * eikonal,const double & smin) {
  AddBeamBlobs(blobs);
  Blob * blob;
  for (Blob_List::iterator biter=blobs->begin();biter!=blobs->end();biter++) {
    blob = (*biter);
    if (blob->Has(blob_status::needs_beams)) {
      if (blob->Type()==btp::Beam || blob->Type()==btp::Shower) {
	blob->UnsetStatus(blob_status::needs_beams);
	LinkShowerInitiators(blob);
      }
      else if (blob->Type()==btp::QElastic_Collision) {
	if (m_hadrons[0]->Elastic() && m_hadrons[1]->Elastic()) {
	  if (blob->NInP()==2) {
	    Particle * bp0(blob->InParticle(0));
	    Particle * bp1(blob->InParticle(1));
	    Particle * hp0(m_hadrons[0]->GetBeamBlob()->InParticle(0));
	    Particle * hp1(m_hadrons[1]->GetBeamBlob()->InParticle(0));
	    if ((bp0->Momentum()-hp0->Momentum()).Abs2()<1.e-3 &&
		(bp1->Momentum()-hp1->Momentum()).Abs2()<1.e-3) {
	      m_hadrons[0]->GetBeamBlob()->AddToOutParticles(bp0);
	      m_hadrons[1]->GetBeamBlob()->AddToOutParticles(bp1);
	      return Return_Value::Success;
	    }
	    if ((bp0->Momentum()-hp1->Momentum()).Abs2()<1.e-3 &&
		(bp1->Momentum()-hp0->Momentum()).Abs2()<1.e-3) {
	      m_hadrons[1]->GetBeamBlob()->AddToOutParticles(bp0);
	      m_hadrons[0]->GetBeamBlob()->AddToOutParticles(bp1);
	      return Return_Value::Success;
	    }
	  }
	  msg_Error()<<"Problem in "<<METHOD<<":\n"
		     <<"  Could not map \n"
		     <<(*blob)<<"\n"
		     <<(*m_hadrons[0]->GetBeamBlob())<<"\n"
		     <<(*m_hadrons[1]->GetBeamBlob())<<"\n"
		     <<"will continue with new event and hope for the best.\n";
	  return Return_Value::Retry_Event;
	}
      }
    }
  }

  AddSpectators();

  for (int beam=0;beam<2;beam++) {
    if (m_hadrons[beam]->GetBeamBlob()->
	CheckMomentumConservation().Abs2()>1.e-6) {
      msg_Error()<<"Problem in "<<METHOD<<":\n"
		 <<"   Beam blob ("<<m_hadrons[beam]->GetBeamBlob()->Id()
		 <<") seems fishy.\n"
		 <<(*m_hadrons[beam]->GetBeamBlob())<<"\n";
    }
  }
  if (p_blob->CheckMomentumConservation().Abs2()>1.) {
    msg_Error()<<"Problem in "<<METHOD<<":\n"
	       <<"   Soft blob ("<<p_blob->Id()<<") seems fishy:\n"
	       <<"   Delta = "<<p_blob->CheckMomentumConservation()<<"\n"
	       <<(*p_blob)<<"\n";
  }
  if (!p_blob->CheckColour()) {
    msg_Error()<<"Problem in "<<METHOD<<":\n   Extra blob "
	       <<"("<<p_blob->Id()<<") seems fishy: "
	       <<"Bad colour configuration.\n"<<(*p_blob)<<"\n";
    return Return_Value::Retry_Event;
  }
  for (size_t beam=0;beam<2;beam++) {
    if (!m_hadrons[beam]->GetBeamBlob()->CheckColour()) {
      msg_Error()<<"Problem in "<<METHOD<<":\n   Extra blob "
		 <<"("<<m_hadrons[beam]->GetBeamBlob()->Id()<<") seems fishy: "
		 <<"Bad colour configuration.\n"
		 <<(*m_hadrons[beam]->GetBeamBlob())<<"\n";
      return Return_Value::Retry_Event;
    }
  }
  if (p_reconnections->FinishConfiguration(blobs,smin)) 
    return Return_Value::Success;
  return Return_Value::Retry_Event;
}


void Beam_Remnant_Handler::LinkShowerInitiators(Blob * blob) {
  Particle_Vector inps(blob->GetInParticles());
  Particle_Vector outps(p_blob->GetOutParticles());
  for (size_t i=0;i<inps.size();i++) {
    if (inps[i]->ProductionBlob()) continue;
    for (size_t j=0;j<outps.size();j++) {
      if (inps[i]->Flav()==outps[j]->Flav() &&
	  inps[i]->Momentum()==outps[j]->Momentum()) {
	p_blob->DeleteOutParticle(outps[j]);
	p_blob->AddToOutParticles(inps[i]);
      }
    }
  }
}

void Beam_Remnant_Handler::AddSpectators() {
  p_colour->FinalColours();

  size_t size(m_hadrons[0]->Size());
  Particle * part1, * part2;
  while (m_paircounter<size && NextIS(part1,part2)) {};
  msg_Tracking()<<"After "<<METHOD<<":\n"<<(*p_blob)<<"\n";
}

bool Beam_Remnant_Handler::
NextIS(Particle *& part1,Particle *& part2) {
  for (int beam=0;beam<2;beam++) {
    p_part[beam] = new Particle(*m_hadrons[beam]->GetParticle(m_paircounter));
    if (!p_part[beam]) return false;
    p_part[beam]->SetNumber();
  }
  double xp[2], xm[2], xt2[2];
  for (int beam=0;beam<2;beam++) 
    m_hadrons[beam]->GetXs(m_paircounter,xp[beam],xm[beam],xt2[beam],m_shat);
  double Xp = xp[0] + xp[1],        Xm = xm[0] + xm[1];
  double tp = (xt2[0] + xt2[1])/Xm, tm = (xt2[0] - xt2[1])/Xm;
  double Px = 0.,                   Mx = 0.;

  Px += xp[0] = 1./2. * (Xp*(1.+sqrt(1.-2.*tp-tm*tm)) - tm);  
  Mx += xm[1] = 1./2. * (Xm*(1.+sqrt(1.-2.*tp-tm*tm)) - tm);  
  Mx += xm[0] = xt2[0]/xp[0];
  Px += xp[1] = xt2[1]/xm[1];

  for (int beam=0;beam<2;beam++) {
    p_part[beam]->SetMomentum(Vec4D(xp[beam]*m_beamvecs[0]+
				    xm[beam]*m_beamvecs[1]+
				    m_hadrons[beam]->Kperp(m_paircounter)));
    p_blob->AddToOutParticles(p_part[beam]);
    m_checkmom[beam] += p_part[beam]->Momentum();
    if (m_paircounter+1 == m_hadrons[beam]->Size()
        && !ATOOLS::IsEqual(m_checkmom[beam],m_beamvecs[beam],1e-3))
      msg_Out()<<METHOD<<" Four Momentum not conserved in intial state: "
                       <<(m_checkmom[beam]-m_beamvecs[beam])/
	m_beamvecs[beam][0]<<std::endl;
  }	   
  m_paircounter++;
  part1 = p_part[0]; part2 = p_part[1];

  return true;
}

bool Beam_Remnant_Handler::UpdateColours(Ladder * ladder,const bool & last) {
  msg_Tracking()<<METHOD<<"(last = "<<last<<"):\n"<<(*ladder);
  if (!(*p_colour)(ladder,p_part,m_paircounter-1)) {
    msg_Error()<<"Error in "<<METHOD<<"(last = "<<last<<"):\n"<<(*ladder);
    return false;
  }
  return true;
}

void Beam_Remnant_Handler::Reset(const size_t & mode) {
  for (int beam=0;beam<2;beam++) 
    m_hadrons[beam]->Reset(m_beamvecs[beam]);
  m_shat        = (m_beams[0]->OutMomentum()+m_beams[1]->OutMomentum()).Abs2();
  m_paircounter = 0;
  m_checkmom[0] = ATOOLS::Vec4D(0,0,0,0);
  m_checkmom[1] = ATOOLS::Vec4D(0,0,0,0);
  p_colour->Reset();

  if (mode>0) {
    if (p_blob && (p_blob->NInP()>0 || p_blob->NOutP()>0)) { 
      delete p_blob; p_blob = NULL; 
    }
  }
}
