#include "PHASIC++/Decays/Decay_Clustering.H"

#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "ATOOLS/Org/Return_Value.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Math/Random.H"

#include "PHASIC++/Decays/Decay_Handler.H"

using namespace PHASIC;
using namespace ATOOLS;
using namespace std;

void Decay_Clustering::DefineInitialShowerConditions(Blob* initialblob, Blob* showerblob)
{
  DEBUG_FUNC(showerblob->Id());

  Blob_Data_Base * bdb((*showerblob)["ClusterAmplitude"]);
  if (!bdb) THROW(fatal_error, "Internal Error");
  Cluster_Amplitude* ampl=bdb->Get<Cluster_Amplitude*>();

  DEBUG_VAR(*ampl);
  for (int i=0; i<initialblob->NOutP(); ++i) {
    ampl->Leg(initialblob->NInP()+i)->SetMom
      (initialblob->OutParticle(i)->Momentum());
  }
  if (ampl->NIn()==2) {
    for (Cluster_Amplitude *campl(ampl);
	 campl;campl=campl->Next()) {
      if (-campl->Leg(0)->Mom()[0]>rpa->gen.PBeam(0)[0] ||
	  -campl->Leg(1)->Mom()[0]>rpa->gen.PBeam(1)[0]) {
        PRINT_INFO("Decay clustering failed. Retry event.");
	throw Return_Value::Retry_Event;
      }
    }
  }
  size_t imax=ampl->Legs().size()-1;
  for (int i=0; i<initialblob->NOutP(); ++i) {
    if (!initialblob->OutParticle(i)->Flav().IsStable()) {
      AddDecayClustering(ampl, initialblob->OutParticle(i)->DecayBlob(),
                         imax, 1<<(initialblob->NInP()+i));
    }
  }
}

class ParticlePairFirstEnergySort {
public:
  bool operator()(const ParticlePair& a,const ParticlePair& b)
  { return (a.first->Momentum()[0]<b.first->Momentum()[0]); }
};

void Decay_Clustering::AddDecayClustering(ATOOLS::Cluster_Amplitude*& ampl,
                                       Blob* blob,
                                       size_t& imax,
                                       size_t idmother)
{
  DEBUG_FUNC("blob->Id()="<<blob->Id()<<" idmother="<<ID(idmother));
  DEBUG_VAR(*blob);
  Particle_Vector daughters;
  ParticlePair_Vector photons;
  for (auto p: blob->GetOutParticles()) {
    if (p->Info()=='S') photons.push_back(make_pair(p,p));
    else                daughters.push_back(p);
  }
  std::sort(photons.begin(),photons.end(),ParticlePairFirstEnergySort());
  msg_Debugging()<<"daughters: ";
  for (size_t i(0);i<daughters.size();++i)
    msg_Debugging()<<daughters[i]->Flav().IDName()<<" ";
  msg_Debugging()<<" +  "<<photons.size()<<" soft photon(s)"<<std::endl;
  AssignPhotons(daughters,photons);
  if (daughters.size()==2) {
    msg_Debugging()<<"1 to 2 case"<<std::endl;
    Cluster_Amplitude* copy=ampl->InitPrev();
    copy->CopyFrom(ampl);
    copy->SetNLO(0);
    copy->SetFlag(1);
    copy->SetMS(ampl->MS());
    Cluster_Leg *lij(ampl->IdLeg(idmother));
    copy->SetKT2(lij->Mom().Abs2());
    for (size_t i=0; i<ampl->Legs().size(); ++i)
      ampl->Leg(i)->SetStat(ampl->Leg(i)->Stat()|1);
    lij->SetStat(1|2|4);
    size_t idk(0);
    for (size_t i=0; i<copy->Legs().size(); ++i) {
      copy->Leg(i)->SetK(0);
      if (copy->Leg(i)->Id()!=idmother &&
          (copy->Leg(i)->Col().m_i==lij->Col().m_j ||
           copy->Leg(i)->Col().m_j==lij->Col().m_i))
        idk=copy->Leg(i)->Id();
    }
    if (lij->Col().m_i==0 && lij->Col().m_j==0) {
      // Ad hoc EW partner
      size_t ampl_nout=ampl->Legs().size()-ampl->NIn();
      if (ampl_nout==1) idk=ampl->Leg(0)->Id();
      else {
        size_t select(0);
        do {
          select=ampl->NIn()+floor(ran->Get()*ampl_nout);
        } while (ampl->Leg(select)->Id()&idmother ||
                 select>ampl->Legs().size()-1);
        idk=ampl->Leg(select)->Id();
      }
    }
    if (idk==0) THROW(fatal_error,"Colour partner not found");
    lij->SetK(idk);
    Cluster_Leg *d1(copy->IdLeg(idmother));
    size_t stat1(0), stat2(0);
    d1->SetMom(RecombinedMomentum(daughters[0],photons,stat1));
    d1->SetStat(stat1);
    d1->SetFlav(daughters[0]->Flav());
    copy->CreateLeg(RecombinedMomentum(daughters[1],photons,stat2),
                    daughters[1]->RefFlav());
    size_t idnew=1<<(++imax);
    copy->Legs().back()->SetId(idnew);
    copy->Legs().back()->SetStat(stat2);
    Cluster_Amplitude::SetColours(ampl->IdLeg(idmother),
                                  copy->IdLeg(idmother),
                                  copy->Legs().back());
    copy->SetIdNew(idnew);
    DEBUG_VAR(*copy);
    Cluster_Amplitude* tmp=copy;
    while (tmp->Next()) {
      tmp=tmp->Next();
      if (tmp->IdNew()&idmother) tmp->SetIdNew(tmp->IdNew()|idnew);
      for (size_t i=0; i<tmp->Legs().size(); ++i) {
        if (tmp->Leg(i)->Id()&idmother) {
          tmp->Leg(i)->SetId(tmp->Leg(i)->Id()|idnew);
	}
        if (tmp->Leg(i)->K()&idmother) {
          tmp->Leg(i)->SetK(tmp->Leg(i)->K()|idnew);
        }
      }
      DEBUG_VAR(*tmp);
    }
    std::vector<size_t> ids;
    ids.push_back(idmother);
    ids.push_back(idnew);
    while (photons.size())
      AddPhotonsClustering(copy, daughters, photons, imax, ids);
    if (!daughters[0]->Flav().IsStable())
      AddDecayClustering(copy, daughters[0]->DecayBlob(), imax, idmother);
    if (!daughters[1]->Flav().IsStable())
      AddDecayClustering(copy, daughters[1]->DecayBlob(), imax, idnew);
    ampl=copy;
  }
  else if (daughters.size()==3) {
    PRINT_INFO("Error: showering of 1->3 decays not implemented yet.");
    throw Return_Value::Retry_Event;
    /*
    msg_Debugging()<<"1 to 3 case"<<std::endl;
    // structure m -> 0 P[->1 2]
    // propagator always combines daughters 1+2
    Cluster_Amplitude* step1=ampl->InitPrev();
    step1->CopyFrom(ampl);
    step1->SetNLO(0);
    step1->SetFlag(1);
    step1->SetMS(ampl->MS());
    Cluster_Leg *lij(ampl->IdLeg(idmother));
    step1->SetKT2(lij->Mom().Abs2());
    if (!lij) THROW(fatal_error,"Cluster leg of id "+ToString(idmother)
                                +" not found.");
    for (size_t i=0; i<ampl->Legs().size(); ++i)
      ampl->Leg(i)->SetStat(ampl->Leg(i)->Stat()|1);
    lij->SetStat(1|2|4);
    size_t idk(0);
    for (size_t i=0; i<step1->Legs().size(); ++i) {
      step1->Leg(i)->SetK(0);
      if (step1->Leg(i)->Id()!=idmother)
	if (step1->Leg(i)->Col().m_i==lij->Col().m_j ||
	    step1->Leg(i)->Col().m_j==lij->Col().m_i) 
	  idk=step1->Leg(i)->Id();
    }
    if (lij->Col().m_i==0 && lij->Col().m_j==0) {
      // Ad hoc EW partner
      size_t ampl_nout=ampl->Legs().size()-ampl->NIn();
      if (ampl_nout==1) idk=ampl->Leg(0)->Id();
      else {
        size_t select(0);
        do {
          select=ampl->NIn()+floor(ran->Get()*ampl_nout);
        } while (ampl->Leg(select)->Id()&idmother || select>ampl->Legs().size()-1);
        idk=ampl->Leg(select)->Id();
      }
    }
    if (idk==0) THROW(fatal_error,"Colour partner not found");
    lij->SetK(idk);
    Cluster_Leg *d1(step1->IdLeg(idmother));
    size_t stat1(0),stat2(0),stat3(0);
    d1->SetMom(RecombinedMomentum(daughters[0],photons,stat1));
    d1->SetStat(stat1);
    d1->SetFlav(daughters[0]->Flav());
    // todo: 1->2 qcd shower with ew fs recoil partner
    // d1->SetK(idmother);// not that simple: w->qq' has color connection in fs
    Decay_Channel* dc(NULL);
    Blob_Data_Base* data = (*blob)["dc"];
    if (data) {
      dc=data->Get<Decay_Channel*>();
      DEBUG_VAR(*dc);
    }
    else THROW(fatal_error, "Internal error.");
    Comix1to3* amp=dynamic_cast<Comix1to3*>(dc->GetDiagrams()[0]);
    if (!amp) THROW(fatal_error, "Internal error.");
    Flavour prop_flav=amp->Prop();
    Vec4D momd2=RecombinedMomentum(daughters[1],photons,stat2);
    Vec4D momd3=RecombinedMomentum(daughters[2],photons,stat3);
    Vec4D prop_mom=momd2+momd3;
    step1->CreateLeg(prop_mom, prop_flav);
    size_t idnew1=1<<(++imax);
    step1->Legs().back()->SetId(idnew1);
    step1->Legs().back()->SetStat(0);
    Cluster_Amplitude::SetColours(ampl->IdLeg(idmother),
                                  step1->IdLeg(idmother),
                                  step1->Legs().back());
    step1->SetIdNew(idnew1);
    DEBUG_VAR(*step1);
    Cluster_Amplitude* tmp=step1;
    while (tmp->Next()) {
      tmp=tmp->Next();
      if (tmp->IdNew()&idmother) tmp->SetIdNew(tmp->IdNew()|idnew1);
      for (size_t i=0; i<tmp->Legs().size(); ++i) {
        if (tmp->Leg(i)->Id()&idmother) {
          tmp->Leg(i)->SetId(tmp->Leg(i)->Id()|idnew1);
	}
        if (tmp->Leg(i)->K()&idmother) {
          tmp->Leg(i)->SetK(tmp->Leg(i)->K()|idnew1);
        }
      }
      DEBUG_VAR(*tmp);
    }

    
    Cluster_Amplitude* step2=step1->InitPrev();
    step2->CopyFrom(step1);
    step2->SetNLO(0);
    step2->SetFlag(1);
    step2->SetMS(step1->MS());
    for (size_t i=0; i<step1->Legs().size(); ++i)
      step1->Leg(i)->SetStat(step1->Leg(i)->Stat()|1);
    step1->IdLeg(idnew1)->SetStat(1|4);
    step1->IdLeg(idnew1)->SetK(idk);
    for (size_t i=0; i<step2->Legs().size(); ++i) step2->Leg(i)->SetK(0);
    Cluster_Leg *d2(step2->IdLeg(idnew1));
    d2->SetMom(momd2);
    d2->SetStat(stat2);
    d2->SetFlav(daughters[1]->Flav());
    step2->CreateLeg(momd3, daughters[2]->Flav());
    size_t idnew2=1<<(++imax);
    step2->Legs().back()->SetId(idnew2);
    step2->Legs().back()->SetStat(stat3);
    Cluster_Amplitude::SetColours(step1->IdLeg(idnew1),
                                  step2->IdLeg(idnew1),
                                  step2->Legs().back());
    step2->SetIdNew(idnew2);
    DEBUG_VAR(*step2);
    tmp=step2;
    while (tmp->Next()) {
      tmp=tmp->Next();
      if (tmp->IdNew()&(idmother|idnew1))
	tmp->SetIdNew(tmp->IdNew()|idnew2);
      for (size_t i=0; i<tmp->Legs().size(); ++i) {
        if (tmp->Leg(i)->Id()&idnew1) {
          tmp->Leg(i)->SetId(tmp->Leg(i)->Id()|idnew2);
	}
        if (tmp->Leg(i)->K()&idnew1) {
          tmp->Leg(i)->SetK(tmp->Leg(i)->K()|idnew2);
        }
      }
      DEBUG_VAR(*tmp);
    }

    std::vector<size_t> ids;
    ids.push_back(idmother);
    ids.push_back(idnew1);
    ids.push_back(idnew2);
    while (photons.size())
      AddPhotonsClustering(step2,daughters,photons,imax,ids);
    if (daughters[0]->DecayBlob())
      AddDecayClustering(step2, daughters[0]->DecayBlob(), imax, idmother);
    if (daughters[1]->DecayBlob())
      AddDecayClustering(step2, daughters[1]->DecayBlob(), imax, idnew1);
    if (daughters[2]->DecayBlob())
      AddDecayClustering(step2, daughters[2]->DecayBlob(), imax, idnew2);
    ampl=step2;
    */
  }
  else {
    PRINT_VAR(*blob);
    THROW(fatal_error, "1 -> n not implemented yet.");
  }
  DEBUG_VAR(*ampl);
}

void Decay_Clustering::AddPhotonsClustering(Cluster_Amplitude*& ampl,
                                         const Particle_Vector daughters,
                                         ParticlePair_Vector& photons,
                                         size_t& imax,
                                         const std::vector<size_t>& ids)
{
  DEBUG_FUNC(photons.size()<<" photons to be clustered");
  Particle * photon(photons.back().first);
  Particle * daughter(photons.back().second);
  photons.pop_back();
  size_t idmother(0);
  if      (daughter==daughters[0]) idmother=ids[0];
  else if (daughter==daughters[1]) idmother=ids[1];
  else if (daughter==daughters[2]) idmother=ids[2];
  else THROW(fatal_error,"Did not find id for "+daughter->Flav().IDName());
  msg_Debugging()<<"Cluster "<<photon->Flav()<<" "<<photon->Momentum()
                <<" with "<<daughter->Flav()<<" "<<ID(idmother)<<std::endl;
  Cluster_Amplitude* copy=ampl->InitPrev();
  copy->CopyFrom(ampl);
  copy->SetNLO(0);
  copy->SetFlag(1);
  copy->SetMS(ampl->MS());
  Cluster_Leg *lij(ampl->IdLeg(idmother));
  for (size_t i=0; i<ampl->Legs().size(); ++i)
    ampl->Leg(i)->SetStat(ampl->Leg(i)->Stat()|1);
  lij->SetStat(1|2|4);
  size_t idk(0);
  for (size_t i=0; i<copy->Legs().size(); ++i) {
    copy->Leg(i)->SetK(0);
  }
  if (lij->Col().m_i!=0 || lij->Col().m_j!=0) {
    THROW(fatal_error,"Adding QED to coloured particle.");
  }
  // Ad hoc QED partner, must not be another soft photon
  size_t ampl_nout=ampl->Legs().size()-ampl->NIn();
  if (ampl_nout==1) idk=ampl->Leg(0)->Id();
  else {
    size_t select(0);
    size_t nvalid(0);
    for (size_t i(ampl->NIn());i<ampl->Legs().size();++i) {
      if (!(ampl->Leg(i)->Id()&idmother || i>ampl->Legs().size()-1 ||
            ampl->Leg(i)->Flav().Kfcode()==kf_photon)) {
        nvalid++;
      }
    }
    if (nvalid==0) select=0;
    else {
      do {
        select=ampl->NIn()+floor(ran->Get()*ampl_nout);
      } while (ampl->Leg(select)->Id()&idmother ||
          select>ampl->Legs().size()-1 ||
          ampl->Leg(select)->Flav().Kfcode()==kf_photon);
    }
    msg_Debugging()<<"choose ("<<ID(ampl->Leg(select)->Id())<<") "
      <<ampl->Leg(select)->Flav()<<std::endl;
    idk=ampl->Leg(select)->Id();
  }
  if (idk==0) THROW(fatal_error,"Colour partner not found");
  lij->SetK(idk);
  Cluster_Leg *d1(copy->IdLeg(idmother));
  size_t stat1(0), stat2(0);
  d1->SetMom(RecombinedMomentum(daughter,photons,stat1));
  d1->SetStat(stat1);
  d1->SetFlav(daughter->Flav());
  copy->CreateLeg(photon->Momentum(),photon->RefFlav());
  size_t idnew=1<<(++imax);
  copy->Legs().back()->SetId(idnew);
  copy->Legs().back()->SetStat(stat2);
  Cluster_Amplitude::SetColours(ampl->IdLeg(idmother),
                                copy->IdLeg(idmother),
                                copy->Legs().back());
  copy->SetIdNew(idnew);
  DEBUG_VAR(*copy);
  Cluster_Amplitude* tmp=copy;
  while (tmp->Next()) {
    tmp=tmp->Next();
    if (tmp->IdNew()&idmother) tmp->SetIdNew(tmp->IdNew()|idnew);
    for (size_t i=0; i<tmp->Legs().size(); ++i) {
      if (tmp->Leg(i)->Id()&idmother) {
        tmp->Leg(i)->SetId(tmp->Leg(i)->Id()|idnew);
      }
      if (tmp->Leg(i)->K()&idmother) {
        tmp->Leg(i)->SetK(tmp->Leg(i)->K()|idnew);
      }
    }
    DEBUG_VAR(*tmp);
  }
  ampl=copy;
}

void Decay_Clustering::AssignPhotons(const Particle_Vector& daughters,
                                  ParticlePair_Vector& photons)
{
  // for every photon, find charged particle that's closest
  // ignore radiation off charged resonance for now
  if (photons.size()) {
    Particle_Vector cdaughters;
    for (size_t i(0);i<daughters.size();++i)
      if (daughters[i]->Flav().Charge()) cdaughters.push_back(daughters[i]);
    if (cdaughters.size()==1) {
      for (size_t i(0);i<photons.size();++i)
        photons[i].second=cdaughters[0];
    }
    else {
      Vec4D cmom(0.,0.,0.,0.);
      Vec4D_Vector cmoms;
      for (size_t i(0);i<cdaughters.size();++i) {
        cmoms.push_back(cdaughters[i]->Momentum());
        cmom+=cmoms[i];
      }
      Poincare ccms(cmom);
      for (size_t i(0);i<cdaughters.size();++i) ccms.Boost(cmoms[i]);
      for (size_t i(0);i<photons.size();++i){
        Vec4D pmom(photons[i].first->Momentum());
        ccms.Boost(pmom);
        size_t id(0);
        double dR(pmom.DR(cmoms[0]));
        for (size_t j(1);j<cmoms.size();++j) {
          double dRj(pmom.DR(cmoms[j]));
          if (dRj<dR) { id=j; dR=dRj; }
        }
        photons[i].second=cdaughters[id];
      }
    }
    for (size_t i(0);i<photons.size();++i) {
      if (photons[i].first==photons[i].second)
        THROW(fatal_error,"Photon has not been assigned.");
      msg_Debugging()<<photons[i].first->Flav()<<" "
                     <<photons[i].first->Momentum()
                     <<" assigned to "<<photons[i].second->Flav()<<std::endl;
    }
  }
}

Vec4D Decay_Clustering::RecombinedMomentum(const Particle * daughter,
                                        const ParticlePair_Vector& photons,
                                        size_t& stat)
{
  Vec4D mom(0.,0.,0.,0.);
  for (size_t i(0);i<photons.size();++i) {
    if (photons[i].second==daughter) {
      mom+=photons[i].first->Momentum();
      stat|=2|4;
    }
  }
  msg_Debugging()<<daughter->Flav()<<": "<<mom<<" "<<stat<<std::endl;
  return mom+daughter->Momentum();
}
