#include "REMNANTS/Tools/Kinematics_Generator.H"
#include "REMNANTS/Main/Remnant_Handler.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Message.H"

using namespace REMNANTS;
using namespace ATOOLS;
using namespace std;


Kinematics_Generator::Kinematics_Generator() {}

Kinematics_Generator::~Kinematics_Generator() {}

void Kinematics_Generator::
Initialize(Remnant_Handler * const rhandler,const string & path, const string & file) {
  p_rhandler = rhandler;
  for (size_t beam=0;beam<2;beam++) {
    p_remnants[beam]   = p_rhandler->GetRemnant(beam);
    p_extracted[beam]  = p_remnants[beam]->GetExtracted();
    p_spectators[beam] = p_remnants[beam]->GetSpectators();
  }
  if (p_rhandler->Type()!=strat::simple) m_kperpGenerator.Initialize(path,file);
  SetKinType(rhandler);
}

void Kinematics_Generator::SetKinType(Remnant_Handler * const rhandler) {
  switch (size_t(p_rhandler->Type())) {
  case size_t(strat::simple):
    m_kintype = kin_type::intact;
    break;
  case size_t(strat::ll):
    m_kintype = kin_type::coll;
    break;
  case size_t(strat::DIS1):
    m_kintype = kin_type::DIS1;
    break;
 case size_t(strat::DIS2):
    m_kintype = kin_type::DIS2;
    break;
  case size_t(strat::hh):
    m_kintype = kin_type::hh;
    break;
  }
}

void Kinematics_Generator::Reset() {
  if (m_kintype==kin_type::intact) return;
  for (size_t beam=0;beam<2;beam++) {
    m_ktmap.clear();
    m_shuffledmap.clear();
  }
  m_boostedblobs.clear();
}

Blob * Kinematics_Generator::MakeSoftBlob() {
  // Create blob for local four-momentum conservation during kinematics
  // reshuffling due to construction of transverse momenta
  p_softblob = new Blob();
  p_softblob->SetType(btp::Soft_Collision);
  p_softblob->SetStatus(blob_status::needs_hadronization);
  p_softblob->SetId();
  return p_softblob;
}
  
bool Kinematics_Generator::FillBlobs(Blob_List * blobs) {
  // Collinear kinematics only - the remnants will only fill extracted particles
  // and possible remnant particles without any need of kinematic reshuffling
  if ((m_kintype==kin_type::intact || m_kintype==kin_type::coll)) {
    return CollinearKinematics();
  }
  // Full kinematics including transverse degrees of freedom - as a consequence
  // this is more involved.
  return TransverseKinematics();
}

bool Kinematics_Generator::CollinearKinematics() {
  for (size_t beam=0;beam<2;beam++) {
    // if beam blobs cannot be filled return false and trigger retrial
    // By far and large here we have a fixed spectator, if necessary, and assign
    // the four-momentum difference between incoming beam particle and outgoing
    // shower initiator to it.
    if (!p_remnants[beam]->FillBlob()) return false;
    m_inmom[beam] = p_remnants[beam]->InMomentum();
  }
  return true;
}

bool Kinematics_Generator::TransverseKinematics() {
  // Remnants fill the beam blobs with spectators and copies of the shower initiators
  // (they will become outging particles of the soft blob, as well as copies of the spectators).
  for (size_t beam=0;beam<2;beam++) {
    // if beam blobs cannot be filled return false and trigger retrial
    if (!p_remnants[beam]->FillBlob(&m_ktmap)) return false;
    m_inmom[beam] = p_remnants[beam]->InMomentum();
  }
  // Initialise particle-momentum maps to track the transverse momenta
  InitKTMaps();
  // Distribute transverse momenta - this will involve mostly minor reshuffling of longitudinal
  // momenta in the remnant break-up.  The check is to make sure this does not violate momentum
  // conservation, as a by-product is already produces the momenta used in the boosting
  // of the connected blobs.  If we produce too large transverse momenta we start by scaling them
  // down by factors of 10 after 100 mistrials.  If we have to scale them down by a factor of
  // 1000 we force all of them to be exactly zero.
  // I still have to think about an error treatment in case this goes wrong.
  size_t maxnum = 100;
  double scale  = 1.;
  do {
    for (short unsigned int beam=0;beam<2;++beam)
      m_kperpGenerator.CreateBreakupKinematics(beam,&m_ktmap,scale);
    maxnum--;
    if (maxnum<=0)   {
      maxnum = 100; scale *= 0.1;
      msg_Error()<<"Warning: "<<METHOD<<" reduces overall prescale for kt to scale = "<<scale<<"\n";
    }
    if (scale<1.e-3) scale = 0.;
  } while (!Check() && scale>0. && maxnum>=0);
  // Fill particles from remnant break-up into soft blob, unless we have simple
  // collinear kinematics with no momentum shuffling
  for (size_t beam=0;beam<2;beam++) {
    Blob * beamblob = p_remnants[beam]->GetBlob();
    for (size_t i=0;i<beamblob->NOutP();i++) {
      Particle * part = beamblob->OutParticle(i);
      p_softblob->AddToInParticles(part);
    }
  }
  // Adjust the kinematics, proceed in two steps:
  // - first go through pairs of shower initiators, reconstruct their momenta,
  //   and boost the connected blobs into the new frame (with kT).
  // - then adjust the remnant partons that are not initiating a shower/collision.
  // First, reset momenta of beam particles - successively we will subtract momenta
  // of extracted particles and use this to monitor longitudinal momentum conservation.
  if ((scale<1.e-4 && maxnum<=0) ||
      !AdjustShowerInitiators() || !AdjustRemnants()) return false;
  return true;
}

void Kinematics_Generator::InitKTMaps() {
  // Put shower initiator (extracted) and spectator partons from both beam remnants into one
  // map of particles into produced transverse momenta, to be filled by the KPerp_Generator.
  for (short unsigned int beam=0;beam<2;++beam) {
    for (Part_Iterator pit=p_extracted[beam]->begin();
	 pit!=p_extracted[beam]->end();pit++) {
      m_ktmap[(*pit)] = Vec4D();
    }
    for (Part_Iterator pit=p_spectators[beam]->begin();
	 pit!=p_spectators[beam]->end();pit++) {
      m_ktmap[(*pit)] = Vec4D();
    }
  }
}

bool Kinematics_Generator::Check() {
  bool success  = true;
  Part_Iterator plit[2];
  for (size_t beam=0;beam<2;beam++) {
    plit[beam]       = p_extracted[beam]->begin();
    m_checkmom[beam] = m_inmom[beam];
  }
  Particle * part[2];
  while (plit[0]!=p_extracted[0]->end()) {
    for (size_t beam=0;beam<2;beam++) part[beam] = *plit[beam];    
    if (CheckScatter(part)) { for (size_t beam=0;beam<2;beam++) plit[beam]++; }
    else { success =false; break; }
  }
  if (success) { success = success && CheckRemnants(); }
  return success;
}

bool Kinematics_Generator::CheckScatter(Particle * part[2]) {
  Vec4D labmom[2], kperp[2], oldLab(0.,0.,0.,0.), Kperp(0.,0.,0.,0.);
  double mt2[2];
  // Harvest momenta etc. from both shower/scatter initiators: will cehck, scatter-by-scatter,
  // if new momenta can be constructed and if there is enough energy in the system.
  for (size_t beam=0;beam<2;beam++) {
    oldLab   += labmom[beam] = part[beam]->Momentum();
    Kperp    += kperp[beam]  = m_ktmap[part[beam]];
    mt2[beam] = sqr(part[beam]->Flav().Mass())-kperp[beam].Abs2();
  }
  double M2  = oldLab.Abs2(), M = sqrt(M2), Y = oldLab.Y();
  double KT2 = -Kperp.Abs2(), MT2 = M2+KT2, MT = sqrt(MT2);
  // Make sure the transverse momentum is kinematically viable
  if (sqr(MT2-mt2[0]-mt2[1])<4.*mt2[0]*mt2[1]) {
    msg_Error()<<METHOD<<" throws error: transverse momentum not viable.\n";
    return false;
  }
  // Reconstruct momenta, assuming the rapidity of the produced system is conserved
  // and check here that energies and beam kinematics still work out.
  double pz    = sqrt(sqr(MT2-mt2[0]-mt2[1])-4.*mt2[0]*mt2[1])/(2.*MT), e[2];
  double coshy = cosh(Y), sinhy = sinh(Y);
  for (size_t beam=0;beam<2;beam++) {
    double e  = (MT2+mt2[beam]-mt2[1-beam])/(2.*MT);
    double E  = e*coshy + (beam==0?1.:-1.)*pz*sinhy;
    double PZ = e*sinhy + (beam==0?1.:-1.)*pz*coshy; 
    m_shuffledmap[part[beam]] = labmom[beam] = Vec4D(E,0.,0.,PZ) + kperp[beam];
    m_checkmom[beam] -= Vec4D(E,0.,0.,PZ);
    if (m_checkmom[beam][0]<0.) {
      msg_Error()<<METHOD<<" throws error: no momentum left in beam "
		 <<beam<<", "<<m_checkmom[beam]<<" for M = "<<M<<" and Y = "<<Y<<".\n";
      return false;
    }
  }
  return true;
}

bool Kinematics_Generator::CheckRemnants() {
  // we use the Momenta_Stretcher to adjust the momenta of the remnants:
  // - add all off-shell momenta (longitudinal and transverse momenta) of both beams
  // - boost them in their c.m. system
  // - stretch them to include the particle masses and restore on-shell condition of
  //   momenta with the Momenta_Stretcher
  // - boost back into the lab system
  // - special role of "recoiler" - in hadrons usually the di-quark:
  //   its momentum is reconstructed from the remaining longitudinal momentum of the
  //   beam remnant, after all other particles have been subtracted, this ensures
  //   longitudinal momentum conservation (the transverse one is taken care of
  //   in the KPerp_Generator through its internal compensation, remnant-by-remnant)
  vector<Vec4D> moms;
  vector<double> masses;
  vector<Particle *>parts;
  Vec4D tot(0.,0.,0.,0.), mom, kt[2];
  for (size_t beam=0;beam<2;beam++) {
    Particle * recoiler = p_remnants[beam]->GetRecoiler();
    for (Part_Iterator plit=p_spectators[beam]->begin();
	 plit!=p_spectators[beam]->end();plit++) {
      Particle * part = (*plit);
      if (part==recoiler) continue;
      tot += mom = part->Momentum()+m_ktmap[part];
      parts.push_back(part);
      moms.push_back(mom);
      masses.push_back(part->Flav().Mass());
      m_checkmom[beam] -= part->Momentum();
    }
    // Check if energies still positive - otherwise we are in deep truoble
    if (m_checkmom[beam][0]<0.) {
      msg_Debugging()<<METHOD<<" throws error: no momentum left in beam "
		     <<beam<<", "<<m_checkmom[beam]<<"\n";
      return false;
    }
    tot += mom = m_checkmom[beam] + m_ktmap[recoiler];
    parts.push_back(recoiler);
    moms.push_back(mom);
    masses.push_back(recoiler->Flav().Mass());
  }
  Poincare residualcms(tot);
  // After boosting into their c.m. frame, use the Momenta_Stretcher to rescale
  // particles onto their mass shells and to account for their transvers momenta.  
  for (size_t i=0;i<moms.size();i++) residualcms.Boost(moms[i]);
  if (!m_stretcher.ZeroThem(0,moms) ||
      !m_stretcher.MassThem(0,moms,masses)) {
    if (msg->LevelIsDebugging()) {
      msg_Out()<<METHOD<<" throws error: rescaling impossible.\n";
      for (size_t beam=0;beam<2;beam++) {
	Particle * recoiler = p_remnants[beam]->GetRecoiler();
	for (Part_Iterator plit=p_spectators[beam]->begin();
	     plit!=p_spectators[beam]->end();plit++) {
	  mom = (*plit)->Momentum()+m_ktmap[(*plit)];
	  msg_Out()<<"  "<<(*plit)->Number()<<": "<<mom
		   <<" --> "<<(*plit)->Flav().Mass()<<"\n";
	}
      }
    }
    return false;
  }
  for (size_t i=0;i<moms.size();i++) {
    residualcms.BoostBack(moms[i]);
    m_shuffledmap[parts[i]] = moms[i];
  }
  return true;
}

bool Kinematics_Generator::AdjustShowerInitiators() {
  // iterate pairwise over both sets of extracted particles, i.e. the shower initiators,
  // obtain constructed momenta with kperp from the Kinematocs_Generator and
  // boost the blobs downstream.
  Part_Iterator plit[2];
  for (size_t beam=0;beam<2;beam++) plit[beam]=p_remnants[beam]->GetExtracted()->begin();
  bool runit = true;
  Particle * part[2];
  while (runit) {
    Vec4D oldLab(0.,0.,0.,0.), newLab(0.,0.,0.,0.);
    for (size_t beam=0;beam<2;beam++) {
      part[beam] = (*plit[beam]);
      oldLab    += part[beam]->Momentum();
      newLab    += ShuffledMomentum(part[beam]);
    }
    Blob * blob   = part[0]->DecayBlob();
    m_oldcmsboost = Poincare(oldLab); 
    m_newcmsboost = Poincare(newLab);
    size_t catchit(0);
    if (blob!=part[1]->DecayBlob() || !BoostConnectedBlob(blob,catchit)) {
      msg_Error()<<"Error in "<<METHOD<<": catchit = "<<catchit<<" for \n"<<(*blob)<<"\n";
      exit(1);
    }
    for (size_t beam=0;beam<2;beam++) {
      plit[beam]++;
      part[beam]->SetMomentum(ShuffledMomentum(part[beam]));
      p_softblob->AddToOutParticles(part[beam]);
      if (plit[beam]==p_remnants[beam]->GetExtracted()->end()) runit = false;
    }
  }
  return true;
}


bool Kinematics_Generator::BoostConnectedBlob(ATOOLS::Blob * blob,size_t & catchit) {
  // Iterate recursively through blobs and boost them into their new systems.
  if (++catchit>100) {
    msg_Error()<<METHOD<<": Error\n"<<"   Blob nesting is too deep.\n";
    return false;
  }
  if (blob==NULL || m_boostedblobs.find(blob)!=m_boostedblobs.end()) return true;
  m_boostedblobs.insert(blob);
  btp::code btype   = blob->Type();
  for (size_t i=0;i<blob->NOutP();++i) {
    Particle * part = blob->OutParticle(i);
    Blob * decblob  = part->DecayBlob();
    // only boost blobs that are NOT signal hard processes or hard collisions -
    // i.e. leave the fully perturbative, unshowered parton level untouched.
    btp::code dtype = decblob==NULL?btp::Unspecified:decblob->Type();
    if (btype!=btp::Signal_Process && btype!=btp::Hard_Decay && btype!=btp::Hard_Collision &&
	dtype!=btp::Signal_Process && dtype!=btp::Hard_Decay && dtype!=btp::Hard_Collision) {
      Vec4D mom = part->Momentum();
      m_oldcmsboost.Boost(mom);
      m_newcmsboost.BoostBack(mom);
      part->SetMomentum(mom);
    }
    // boost blobs downstream
    if (!BoostConnectedBlob(part->DecayBlob(),catchit)) return false;
  }
  return true;
}

bool Kinematics_Generator::AdjustRemnants() {
  // iterate over all spectators and give them the new momenta from the Kinematics_Generator
  for (size_t beam=0;beam<2;beam++) {
    for (Part_Iterator plit=p_remnants[beam]->GetSpectators()->begin();
	 plit!=p_remnants[beam]->GetSpectators()->end();plit++) {
      Particle * part = (*plit);
      part->SetMomentum(ShuffledMomentum(part));
      p_softblob->AddToOutParticles(part);
    }
  }
  return true;
}

const Vec4D & Kinematics_Generator::ShuffledMomentum(Particle *const part) {
  if (m_shuffledmap.find(part)!=m_shuffledmap.end())
    return m_shuffledmap.find(part)->second;
  msg_Error()<<"Error in "<<METHOD<<": did not find\n"
	     <<(*part)<<"\n"
	     <<"   will return original momentum.\n";
  abort();
  return part->Momentum();
}



