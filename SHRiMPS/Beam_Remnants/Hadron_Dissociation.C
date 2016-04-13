#include "SHRiMPS/Beam_Remnants/Hadron_Dissociation.H"
#include "BEAM/Main/Beam_Base.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"

using namespace SHRIMPS;
using namespace BEAM;
using namespace ATOOLS;

Hadron_Dissociation::
Hadron_Dissociation(Beam_Base * beambase,Continued_PDF * pdf) :
  p_pdf(pdf),
  m_beamvec(beambase->OutMomentum()), m_outmom(Vec4D(0.,0.,0.,0.)),
  m_beamflav(beambase->Bunch()),
  m_dir(m_beamvec[3]>0.?1:-1), m_xmin(2./m_beamvec[0]), m_QT2max(4.),
  p_blob(NULL)
{ }

void Hadron_Dissociation::Reset() {
  m_outmom = m_beamvec;
  for (size_t i=0;i<2;i++) m_cols[i].clear();
}

bool Hadron_Dissociation::FillBeamBlob(Blob_List * blobs) {
  SpecifyBeamBlob();
  AddPartonsFromCollision(blobs);
  IdentifyAndFillSoftBlob(blobs);
  //msg_Out()<<METHOD<<": "<<m_cols[0].size()<<" "<<m_cols[1].size()<<"\n";
  return true;
}

void Hadron_Dissociation::SpecifyBeamBlob() {
  if (!p_blob) p_blob = new Blob();
  p_blob->SetType(btp::Beam);
  p_blob->SetTypeSpec("Shrimps");
  p_blob->SetStatus(blob_status::inactive);
  p_blob->SetId();
}
  
void Hadron_Dissociation::AddPartonsFromCollision(Blob_List * blobs) {
  for (Blob_List::iterator biter=blobs->begin();biter!=blobs->end();biter++) {
    if ((*biter)->Has(blob_status::needs_beams) &&
	(*biter)->Type()==btp::Shower) {
      HarvestBlob((*biter));
    }
  }
}

void Hadron_Dissociation::HarvestBlob(Blob * blob) {
  for (size_t in=0;in<blob->NInP();in++) {
    Particle * part(blob->InParticle(in));
    if (!part->ProductionBlob() && m_dir*part->Momentum()[3]>0) {
      p_blob->AddToOutParticles(part);
      m_outmom -= part->Momentum();
      for (size_t i=0;i<2;i++) {
	std::set<int>::iterator cit(m_cols[i].find(part->GetFlow(1+i)));
	if (cit==m_cols[i].end()) m_cols[1-i].insert(part->GetFlow(1+i));
	else m_cols[i].erase(cit);
      }
    }
  }
}

void Hadron_Dissociation::IdentifyAndFillSoftBlob(Blob_List * blobs) {
  Blob * softblob = blobs->FindFirst(btp::Soft_Collision);
  softblob->SetTypeSpec("Four_Momentum_Compensation");
  softblob->UnsetStatus(blob_status::needs_minBias);
  softblob->SetStatus(blob_status::needs_hadronization);
  AddSpectatorPartons(softblob);
}

void Hadron_Dissociation::AddSpectatorPartons(Blob * softblob) {
  FixConstituentFlavours();
  Vec4D qmom,dimom;
  CalculateParallelMomenta(qmom,dimom);
  Particle * quark(new Particle(0,m_quark,qmom,'B'));
  quark->SetNumber(-1);
  quark->SetFlow(1,(*m_cols[0].begin()));
  p_blob->AddToOutParticles(quark);
  softblob->AddToInParticles(quark);
  Particle * diquark(new Particle(0,m_diquark,dimom,'B'));
  diquark->SetNumber(-1);
  diquark->SetFlow(2,(*m_cols[1].begin()));
  p_blob->AddToOutParticles(diquark);
  softblob->AddToInParticles(diquark);

  Particle * outquark(new Particle(*quark));
  outquark->SetNumber(-1);
  outquark->SetInfo('F');
  softblob->AddToOutParticles(outquark);
  Particle * outdiquark(new Particle(*diquark));
  outdiquark->SetNumber(-1);
  outdiquark->SetInfo('F');
  softblob->AddToOutParticles(outdiquark);

  (*m_qtmap)[outquark]   = Vec4D(0.,0.,0.,0.);
  (*m_qtmap)[outdiquark] = Vec4D(0.,0.,0.,0.);
}


void Hadron_Dissociation::
CalculateParallelMomenta(Vec4D & qmom,Vec4D & dimom) {
  // assume that diquark has at least 2 GeV energy
  double xmax((m_outmom[0]-2.)/m_beamvec[0]),x(-1.);
  int trials(0);
  while (trials<1000) {
    x = m_xmin+ran->Get()*(xmax-m_xmin);
    p_pdf->Calculate(x,0.);
    if (p_pdf->XPDF(m_quark)/p_pdf->XPDFMax(m_quark)>ran->Get()) break;
  }
  qmom  = x*m_beamvec;
  dimom = m_outmom-qmom;
}

void Hadron_Dissociation::SelectTrialTransverseMomenta() {
  for (std::map<Particle *,Vec4D>::iterator pvit=m_qtmap->begin();
       pvit!=m_qtmap->end();pvit++) {
    double qt  = sqrt(p_ff->SelectQT2(m_QT2max,0.));
    double phi = ran->Get()*2.*M_PI;
    pvit->second = Vec4D(0.,cos(phi),sin(phi),0.);
  }
}

void Hadron_Dissociation::FixConstituentFlavours() {
  double random(ran->Get());
  if (m_beamflav==Flavour(kf_p_plus)) {
    if (random<1./3.) {
      m_quark   = Flavour(kf_d);
      m_diquark = Flavour(kf_uu_1);
    }     
    else if (random<1./2.) {
      m_quark   = Flavour(kf_u);
      m_diquark = Flavour(kf_ud_1);
    }
    else {
      m_quark   = Flavour(kf_u);
      m_diquark = Flavour(kf_ud_0);
    }
  }
  else if (m_beamflav==Flavour(kf_p_plus).Bar()) {
    if (random<1./3.) {
      m_quark   = Flavour(kf_d).Bar();
      m_diquark = Flavour(kf_uu_1).Bar();
    }     
    else if (random<1./2.) {
      m_quark   = Flavour(kf_u).Bar();
      m_diquark = Flavour(kf_ud_1).Bar();
    }
    else {
      m_quark   = Flavour(kf_u).Bar();
      m_diquark = Flavour(kf_ud_0).Bar();
    }
  }
  else {
    msg_Error()<<"Error in "<<METHOD<<"(bunch = "<<m_beamflav<<"):\n"
	       <<"   No parton dissociation found.  Will exit.\n";
    exit(1);
  }
}
