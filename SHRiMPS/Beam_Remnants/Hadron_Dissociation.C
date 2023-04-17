#include "SHRiMPS/Beam_Remnants/Hadron_Dissociation.H"
#include "REMNANTS/Main/Remnant_Base.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"

using namespace SHRIMPS;
using namespace BEAM;
using namespace ATOOLS;

Hadron_Dissociation::
Hadron_Dissociation(const int & beam,
		    const ATOOLS::Vec4D & inmom,const ATOOLS::Flavour & inflav,
		    Continued_PDF * pdf) :
  m_beam(beam),p_pdf(pdf),
  m_beamvec(inmom), m_outmom(Vec4D(0.,0.,0.,0.)),m_beamflav(inflav),
  m_dir(m_beamvec[3]>0.?1:-1), m_xmin(2./m_beamvec[0]), m_QT2max(4.), m_expo(2.),
  p_blob(NULL)
{ }

Hadron_Dissociation::~Hadron_Dissociation() {
  if (p_pdf) delete p_pdf;
}

void Hadron_Dissociation::Reset() {
  m_outmom = m_beamvec;
  for (size_t i=0;i<2;i++) m_beamcols[i].clear();
  m_beamflavs.clear();
  m_qtmap.clear();
  p_blob = NULL;
}

bool Hadron_Dissociation::FillBeamBlob(Blob_List * blobs, const double & B) {
  AddBeamBlob(blobs, B);
  HarvestCollision(blobs);
  if (m_outmom[0] < 2.) {
    msg_Tracking()<<METHOD<<" arrives at residual mom = "<<m_outmom<<"\n"; 
    return false;
  }
  if (!CompensateFlavours() || !CompensateColours()) {
    THROW(fatal_error,"Error compensating flavours or colours");
  }
  AddSpectatorPartons();
  if (!CheckResiduals()) {
    THROW(fatal_error,"Error with residual flavours");
  }
  return true;
}

void Hadron_Dissociation::AddBeamBlob(Blob_List * blobs,const double & B) {
  if (!p_blob) p_blob = new Blob();
  p_blob->SetType(btp::Beam);
  p_blob->SetTypeSpec("Shrimps");
  p_blob->SetStatus(blob_status::inactive);
  p_blob->SetId();
  Particle * inpart = new Particle(-1,m_beamflav,m_beamvec);
  inpart->SetNumber();
  inpart->SetBeam(m_beam);
  p_blob->AddToInParticles(inpart);
  blobs->push_front(p_blob);
  if (m_beam==1) p_blob->SetPosition(Vec4D(0.,B,0.,0.));
}
  
void Hadron_Dissociation::HarvestCollision(Blob_List * blobs) {
  for (Blob_List::iterator biter=blobs->begin();biter!=blobs->end();biter++) {
    Blob * blob = (*biter);
    if (blob->Has(blob_status::needs_beams)) {
      for (size_t in=0;in<blob->NInP();in++) {
	Particle * part(blob->InParticle(in));
	if (!part->ProductionBlob() && m_dir*part->Momentum()[3]>0) {
	  p_blob->AddToOutParticles(part);
	  m_outmom -= part->Momentum();
	  //msg_Out()<<m_beam<<" subtracts: "<<part->Momentum()<<" --> "<<m_outmom[0]<<"\n";
	  if (!part->Flav().IsGluon()) m_beamflavs.push_back(part->Flav());
	  for (size_t i=0;i<2;i++) {
	    if (part->GetFlow(1+i)==0) continue;
	    std::set<int>::iterator cit(m_beamcols[i].find(part->GetFlow(1+i)));
	    if (cit==m_beamcols[i].end()) m_beamcols[1-i].insert(part->GetFlow(1+i));
	    else m_beamcols[i].erase(cit);
	  }
	}
      }
    }
  }
}

bool Hadron_Dissociation::CompensateFlavours() {
  if (m_beamflavs.empty()) return true;
  while (!m_beamflavs.empty()) {
    //msg_Out()<<"Add compensator for "<<m_beamflavs.front()<<"\n";
    if (!AddFlavourCompensator(m_beamflavs.front())) return false;
    m_beamflavs.pop_front();
  }
  return true;
}

bool Hadron_Dissociation::CompensateColours() {
  CleanColours();
  if (m_beamcols[0].size()==1 &&
      m_beamcols[0].size()==m_beamcols[1].size()) return true;
  while (m_beamcols[0].size()>1) {
    Particle * compensator = new Particle(0,Flavour(kf_gluon),Vec4D(0.,0.,0.,0.));
    if (!SelectCompensatorMomentum(compensator)) {
      delete compensator;
      return false;
    }
    for (size_t i=0;i<2;i++) {
      compensator->SetFlow(i+1,(*m_beamcols[i].begin()));
      m_beamcols[i].erase(compensator->GetFlow(i+1));
    }
    compensator->SetNumber();
    compensator->SetBeam(m_beam);
    compensator->SetInfo('B');
    compensator->SetPosition(p_blob->Position());
    p_blob->AddToOutParticles(compensator);
    p_softblob->AddToInParticles(compensator);
    //msg_Out()<<"Added compensator to blob:\n"<<(*compensator)<<"\n";
    
    Particle * out = new Particle(*compensator);
    out->SetPosition(p_blob->Position());
    p_softblob->AddToOutParticles(out);
    m_qtmap[out] = Vec4D(0.,0.,0.,0.);
  }
  return true;
}

void Hadron_Dissociation::CleanColours() {
  std::set<int>::iterator cit[2];
  for (size_t i=0;i<2;i++) cit[i] = m_beamcols[i].begin();
  while (cit[0]!=m_beamcols[0].end()) {
    cit[1] = m_beamcols[1].begin();
    do {
      if ((*cit[0])==(*cit[1])) {
	msg_Out()<<METHOD<<" deletes identical colour in triplet and anti-triplet:"<<(*cit[0])<<"\n";
	for (size_t i=0;i<2;i++) cit[i] = m_beamcols[i].erase(cit[i]);
	break;
      }
      else cit[1]++;
    } while (cit[1]!=m_beamcols[1].end());
    cit[0]++;
  }
}
  
bool Hadron_Dissociation::AddFlavourCompensator(const Flavour & flav) {
  Particle * compensator = new Particle(0,flav.Bar(),Vec4D(0.,0.,0.,0.));
  if (!SelectCompensatorMomentum(compensator) ||
      !SelectCompensatorTripletColours(compensator)) {
    msg_Error()<<"Error in "<<METHOD<<": could not compensate "<<flav<<".\n"
	       <<"   Return false and hope for the best.\n";
    delete compensator;
    return false;
  }
  compensator->SetNumber();
  compensator->SetBeam(m_beam);
  compensator->SetInfo('B');
  p_blob->AddToOutParticles(compensator);
  p_softblob->AddToInParticles(compensator);
  //msg_Out()<<"Added compensator to blob:\n"<<(*compensator)<<"\n";

  Particle * out = new Particle(*compensator);
  p_softblob->AddToOutParticles(out);
  m_qtmap[out] = Vec4D(0.,0.,0.,0.);
  return true;
}

bool Hadron_Dissociation::SelectCompensatorMomentum(Particle * part) {
  double xmax = m_outmom[0]/m_beamvec[0] - 2.*m_xmin;
  if (xmax<m_xmin) return false;
  double rand = ran->Get();
  double x    = pow(1./rand * (pow(xmax,1.-m_expo) + (1.-rand) * pow(m_xmin,1.-m_expo)),
		    1./(1.-m_expo));
  part->SetMomentum(x * m_beamvec);
  m_outmom -= part->Momentum();
  //msg_Out()<<METHOD<<" selected x = "<<x<<" in ["<<m_xmin<<", "<<xmax<<"]: "
  //	   <<"residual momentum = "<<m_outmom<<"\n";
  return true;
}

bool Hadron_Dissociation::SelectCompensatorTripletColours(Particle * part) {
  if (part->Flav().IsQuark()) {
    size_t index = part->Flav().IsAnti()?1:0; 
    if (m_beamcols[index].size()>1) {
      part->SetFlow(index+1,(*m_beamcols[index].begin()));
      m_beamcols[index].erase(part->GetFlow(index+1));
      return true;
    }
    else {
      part->SetFlow(index+1,-1);
      m_beamcols[1-index].insert(part->GetFlow(index+1));
      //msg_Out()<<METHOD<<" adds "<<part->GetFlow(index+1)<<" to beamcols["<<(1-index)<<"]\n";
      return true;
    }
  }
  return false;
}

void Hadron_Dissociation::AddSpectatorPartons() {
  FixConstituentFlavours();
  Vec4D qmom,dimom;
  CalculateParallelMomenta(qmom,dimom);
  Particle * quark(new Particle(0,m_quark,qmom,'B'));
  quark->SetNumber();
  quark->SetBeam(m_beam);
  quark->SetFlow(1,(*m_beamcols[0].begin()));
  quark->SetPosition(p_blob->Position());
  p_blob->AddToOutParticles(quark);
  p_softblob->AddToInParticles(quark);
  Particle * diquark(new Particle(0,m_diquark,dimom,'B'));
  diquark->SetNumber();
  diquark->SetFlow(2,(*m_beamcols[1].begin()));
  diquark->SetBeam(m_beam);
  diquark->SetPosition(p_blob->Position());
  p_blob->AddToOutParticles(diquark);
  p_softblob->AddToInParticles(diquark);

  Particle * outquark(new Particle(*quark));
  outquark->SetNumber();
  outquark->SetInfo('B');
  outquark->SetBeam(m_beam);
  outquark->SetPosition(p_blob->Position());
  p_softblob->AddToOutParticles(outquark);
  Particle * outdiquark(new Particle(*diquark));
  outdiquark->SetNumber();
  outdiquark->SetInfo('B');
  outdiquark->SetBeam(m_beam);
  outdiquark->SetPosition(p_blob->Position());
  p_softblob->AddToOutParticles(outdiquark);

  m_qtmap[outquark]   = Vec4D(0.,0.,0.,0.);
  m_qtmap[outdiquark] = Vec4D(0.,0.,0.,0.);
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
    trials++;
  }
  qmom  = x*m_beamvec;
  dimom = m_outmom-qmom;
}

void Hadron_Dissociation::SelectTrialTransverseMomenta() {
  for (std::map<Particle *,Vec4D>::iterator pvit=m_qtmap.begin();
       pvit!=m_qtmap.end();pvit++) {
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
    THROW(fatal_error,"No parton dissociation found.")
  }
}


bool Hadron_Dissociation::CheckResiduals() {
  if (m_beamcols[0].size()>1 || m_beamcols[1].size()>1 || !m_beamflavs.empty()) {
    msg_Out()<<METHOD<<": "
	     <<"colours = ("<<m_beamcols[0].size()<<" "<<m_beamcols[1].size()<<"), "
	     <<"flavour = "<<m_beamflavs.size()<<".\n";
    for (size_t i=0;i<2;i++) {
      msg_Out()<<"   colours["<<i<<"]: ";
      for (std::set<int>::iterator cit=m_beamcols[i].begin();
	   cit!=m_beamcols[i].end();cit++) msg_Out()<<(*cit)<<" ";
      msg_Out()<<"\n";
    }
    return false;
  }
  return true;
}
