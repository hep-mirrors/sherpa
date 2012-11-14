#include "SHERPA/Single_Events/Decay_Handler_Base.H"

#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Phys/Momenta_Stretcher.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Return_Value.H"

#include "PHASIC++/Decays/Decay_Channel.H"
#include "PHASIC++/Decays/Decay_Table.H"
#include "PHASIC++/Decays/Decay_Map.H"

#include "METOOLS/SpinCorrelations/Spin_Density.H"
#include "METOOLS/SpinCorrelations/Decay_Matrix.H"
#include "METOOLS/SpinCorrelations/Amplitude2_Tensor.H"

using namespace SHERPA;
using namespace PHASIC;
using namespace ATOOLS;
using namespace METOOLS;
using namespace std;

Decay_Handler_Base::Decay_Handler_Base() :
  p_decaymap(NULL), p_bloblist(NULL), p_ampl(NULL),
  m_extraqed(false), m_spincorr(false), m_mass_smearing(1)
{
}

Decay_Handler_Base::~Decay_Handler_Base()
{
  if (p_decaymap) delete p_decaymap; p_decaymap=NULL;
}

class Decay_Width_Sorter {
  PHASIC::Decay_Map* p_decaymap;
public:
  Decay_Width_Sorter(PHASIC::Decay_Map* decaymap) : p_decaymap(decaymap) {}

  bool operator()(Particle* p1, Particle* p2) {
    Decay_Table* table1=p_decaymap->FindDecay(p1->Flav());
    Decay_Table* table2=p_decaymap->FindDecay(p2->Flav());
    double width1(0.0), width2(0.0);
    if (table1) width1=table1->TotalWidth();
    if (table2) width2=table2->TotalWidth();
    return width1 < width2;
  }
};

void Decay_Handler_Base::SetMasses(ATOOLS::Blob* blob, bool usefinalmass)
{
  DEBUG_FUNC(blob->GetOutParticles().size());
  Particle_Vector daughters = blob->GetOutParticles();
  if (m_mass_smearing==0 || daughters.size()==1) return;
  sort(daughters.begin(), daughters.end(), Decay_Width_Sorter(p_decaymap));

  Vec4D total(0.0,0.0,0.0,0.0);
  size_t nr_daughters(0);
  for(size_t i=0; i<daughters.size(); i++) {
    if (!daughters[i]->DecayBlob() || daughters[i]->DecayBlob()->NOutP()==0) {
      total+=daughters[i]->Momentum();
      ++nr_daughters;
    }
  }
  double max_mass;
  if (usefinalmass) max_mass=blob->InParticle(0)->FinalMass();
  else max_mass=total.Mass();
  if (nr_daughters<2) return;
  
  bool success=true; size_t cnt=0;
  do {
    success=true;
    double max = max_mass;
    for(PVIt it=daughters.begin();it!=daughters.end();++it) {
      if(m_mass_smearing==2 && !Decays((*it)->Flav())) continue;
      if ((*it)->DecayBlob() && (*it)->DecayBlob()->NOutP()>0) continue;
      else if ((*it)->DecayBlob()) {
        if (!DiceMass(*it,max)) {
          success=false; ++cnt;
	  if (cnt>22) {
            msg_Error()<<METHOD<<" failed to set masses, retrying event."<<endl;
            throw Return_Value::Retry_Event;
          }
          break;
        }
      }
      else {
        double mass = (*it)->RefFlav().RelBWMass(0.0, max, this->Mass((*it)->RefFlav()));
        (*it)->SetFinalMass(mass);
        DEBUG_INFO(max<<" > "<<"m["<<(*it)->RefFlav()<<"]="<<mass);
      }
      max-=(*it)->FinalMass();
    }
  } while(success==false);
}

bool Decay_Handler_Base::DiceMass(ATOOLS::Particle* p, double max)
{
  Blob* decayblob=p->DecayBlob();
  Blob_Data_Base* data = (*decayblob)["dc"];
  if (data) {
    Decay_Channel* dc = data->Get<Decay_Channel*>();
    double width = p_decaymap->FindDecay(p->Flav())->TotalWidth();
    double mass=dc->GenerateMass(max, width);
    if (mass>0.0) p->SetFinalMass(mass);
    else return false;
  }
  return true;
}

void Decay_Handler_Base::BoostAndStretch(Blob* blob, const Vec4D& labmom)
{
  DEBUG_FUNC("");
  DEBUG_VAR(blob->MomentumConserved());
  // 1.
  Particle* inpart = blob->InParticle(0);
  Vec4D mom(inpart->Momentum());
  double m02=sqr(inpart->FinalMass());
  double p02=mom.PSpat2();
  double E02=sqr(mom[0]);
  double factor=sqrt((m02+p02)/E02);
  DEBUG_VAR(factor);
  mom[0]*=factor;
  inpart->SetMomentum(mom);
  Particle_Vector daughters = blob->GetOutParticles();
  for(PVIt it=daughters.begin();it!=daughters.end();++it) {
    mom=(*it)->Momentum();
    mom[0]*=factor;
    (*it)->SetMomentum(mom);
  }
  DEBUG_VAR(blob->MomentumConserved());

  // 2.
  Poincare twiddle2rest(inpart->Momentum());
  Poincare labboost(labmom);
  labboost.Invert();

  blob->Boost(twiddle2rest);
  blob->Boost(labboost);
  DEBUG_VAR(blob->MomentumConserved());

  // 3.
  Momenta_Stretcher stretch;
  if (!stretch.StretchBlob(blob)) {
    msg_Error()<<METHOD<<" failed to stretch blob, retrying event."<<endl
               <<*blob<<endl;
    throw Return_Value::Retry_Event;
  }
  DEBUG_VAR(blob->MomentumConserved());
}

void Decay_Handler_Base::TreatInitialBlob(ATOOLS::Blob* blob,
                                          METOOLS::Amplitude2_Tensor* amps,
                                          const Particle_Vector& origparts)
{
  DEBUG_FUNC("");
  DEBUG_VAR(*blob);
  // random shuffle, against bias in spin correlations and mixing
  Particle_Vector daughters = blob->GetOutParticles();
  std::vector<size_t> shuffled(daughters.size());
  for (size_t i=0; i<daughters.size(); ++i) shuffled[i]=i;
  for (size_t i=0; i<daughters.size(); ++i) {
    if (abs(daughters[i]->Momentum().Abs2()-sqr(daughters[i]->FinalMass()))>1e-6) {
      PRINT_INFO("Initial particle "<<daughters[i]->Flav()<<" not onshell: "
                 <<"p^2="<<daughters[i]->Momentum().Mass()
                 <<" vs. m^2="<<daughters[i]->FinalMass());
      //      throw Return_Value::Retry_Event;
    }
  }
  random_shuffle(shuffled.begin(), shuffled.end(), *ran);
  
  // initial blobs still contain on-shell particles, stretch them off-shell
  for (size_t i=0; i<daughters.size(); ++i) {
    if (!Decays(daughters[shuffled[i]]->Flav()) ||
        daughters[shuffled[i]]->DecayBlob()) {
      continue;
    }
    CreateDecayBlob(daughters[shuffled[i]]);
  }
  SetMasses(blob, false);
  Momenta_Stretcher stretch;
  if (!stretch.StretchBlob(blob)) {
    msg_Error()<<METHOD<<" failed to stretch blob, retrying event."<<endl
               <<*blob<<endl;
    throw Return_Value::Retry_Event;
  }
  DEBUG_VAR(*blob);

  for (size_t ii(0); ii<daughters.size(); ++ii) {
    size_t i=shuffled[ii];
    DEBUG_INFO("treating "<<*daughters[i]);
    if (m_spincorr) {
      if (amps && origparts[i]) {
        DEBUG_VAR(*amps);
        Decay_Matrix* D(NULL);
        if (!Decays(daughters[i]->Flav()) || !daughters[i]->DecayBlob() ||
            daughters[i]->DecayBlob()->NOutP()>0) {
          D=new Decay_Matrix(origparts[i]); // delta
        }
        else {
          DEBUG_INFO("treating: "<<origparts[i]->Flav());
          Spin_Density sigma(origparts[i],amps);
          sigma.SetParticle(daughters[i]);
          DEBUG_VAR(sigma);
          D=FillDecayTree(daughters[i]->DecayBlob(), &sigma);
          D->SetParticle(origparts[i]);
        }
        if (amps->Contains(origparts[i])) {
          DEBUG_INFO("contracting with D["<<D->Particle()<<"]");
          amps->Contract(D);
        }
        delete D;
      }
      else {
        Spin_Density sigma(daughters[i]);
        if (Decays(daughters[i]->Flav())) {
          Decay_Matrix* D=FillDecayTree(daughters[i]->DecayBlob(), &sigma);
          delete D;
        }
      }
    }
    else {
      if (Decays(daughters[i]->Flav())) {
        FillDecayTree(daughters[i]->DecayBlob(), NULL);
      }
    }
  }
}

Decay_Matrix* Decay_Handler_Base::FillDecayTree(Blob * blob, Spin_Density* s0)
{
  Particle* inpart = blob->InParticle(0);
  DEBUG_FUNC(inpart->RefFlav()<<" "<<inpart->Number());
  if (s0) DEBUG_VAR(*s0);
  Vec4D labmom = inpart->Momentum();
  
  // fill decay blob all on-shell
  Blob_Data_Base* data = (*blob)["p_onshell"];
  if (data) inpart->SetMomentum(data->Get<Vec4D>());
  else {
    msg_Error()<<METHOD<<" could not find p_onshell, retrying event."<<endl
               <<*blob<<endl;
    throw Return_Value::Retry_Event;
  }
  Amplitude2_Tensor* amps=FillOnshellDecay(blob, s0);
  inpart->SetStatus(part_status::decayed);
  inpart->SetInfo('D');

  Particle_Vector daughters = blob->GetOutParticles();
  random_shuffle(daughters.begin(), daughters.end(), *ran);
  for (PVIt it=daughters.begin();it!=daughters.end();++it) {
    if (Decays((*it)->Flav())) {
      if (!CanDecay(inpart->Flav())) {
        msg_Error()<<METHOD<<" Particle '"<<inpart->Flav()<<"' set unstable, "
                   <<"but decay handler doesn't know how to deal with it.";
        throw Return_Value::Retry_Event;
      }
      CreateDecayBlob(*it);
    }
  }

  SetMasses(blob, true);
  BoostAndStretch(blob, labmom);
  DEBUG_VAR(*blob);
  if (m_extraqed) AttachExtraQED(blob);

  DEBUG_INFO("recursively treating the created daughter decay blobs:");
  if (m_spincorr) DEBUG_VAR(*amps);
  for (size_t i(0); i<daughters.size();++i) {
    // have to ignore photons from soft photon handler
    if (daughters[i]->Info()=='S') continue;
    DEBUG_VAR(daughters[i]->Flav());

    if (!Decays(daughters[i]->Flav())) {
      DEBUG_INFO("is stable.");
      if (m_spincorr) {
        Decay_Matrix* D=new Decay_Matrix(daughters[i]);
        amps->Contract(D);
        delete D;
      }
      if (m_spincorr) DEBUG_VAR(*amps);
    }
    else {
      if (!CanDecay(inpart->Flav())) {
        msg_Error()<<METHOD<<" Particle '"<<inpart->Flav()<<"' set unstable, "
                   <<"but decay handler doesn't know how to deal with it."
                   <<endl<<*blob<<endl;
      }
      Spin_Density* si=m_spincorr?new Spin_Density(daughters[i],s0,amps):NULL;
      DEBUG_INFO("decaying with spin density "<<si);
      Decay_Matrix* Di=FillDecayTree(daughters[i]->DecayBlob(), si);
      if (si) delete si;
      if (Di) {
        amps->Contract(Di);
        delete Di;
      }
    }
  }
  DEBUG_INFO("finished daughters of "<<inpart->RefFlav()<<" "<<inpart->Number());
  if (m_spincorr) DEBUG_VAR(*amps);
  return m_spincorr?new Decay_Matrix(inpart,amps):NULL;
}

Amplitude2_Tensor* Decay_Handler_Base::FillOnshellDecay(Blob *blob,
                                                        Spin_Density* sigma)
{
  DEBUG_FUNC("");
  Decay_Channel* dc(NULL);
  Blob_Data_Base* data = (*blob)["dc"];
  if(data) {
    dc=data->Get<Decay_Channel*>();
    DEBUG_VAR(*dc);
  }
  else {
    Decay_Table* table=p_decaymap->FindDecay(blob->InParticle(0)->Flav());
    if (table==NULL) {
      msg_Error()<<METHOD<<"Error: Did not find "<<blob->InParticle(0)->Flav()
                 <<" in decay tables."<<endl;
      throw Return_Value::Retry_Event;
    }
    dc=table->Select();
    blob->AddData("dc",new Blob_Data<Decay_Channel*>(dc));
  }
  
  Particle* inpart=blob->InParticle(0);
  inpart->SetStatus(part_status::decayed);
  Flavour flav; Particle* particle=NULL;
  for (size_t i=1; i<dc->Flavs().size(); ++i) {
    flav=dc->Flavs()[i];
    if (inpart->Flav().IsAnti()) flav = flav.Bar();
    particle = new Particle(0, flav);
    particle->SetFinalMass(Mass(flav));
    particle->SetStatus(part_status::active);
    particle->SetNumber();
    particle->SetInfo('D');
    blob->AddToOutParticles( particle );
  }
  
  size_t n=1+blob->NOutP();
  vector<Vec4D> moms(n);
  moms[0]=inpart->Momentum();

  Amplitude2_Tensor* amps(NULL);
  if (sigma) {
    std::deque<Particle*> parts;
    parts.insert(parts.end(), blob->InParticle(0));
    Particle_Vector outparts=blob->GetOutParticles();
    parts.insert(parts.end(), outparts.begin(), outparts.end());
    dc->GenerateKinematics(moms,inpart->Flav().IsAnti(),sigma,parts);
    amps=dc->Amps();
  }
  else {
    dc->GenerateKinematics(moms,inpart->Flav().IsAnti());
  }
  for (size_t i=1; i<n; i++) blob->GetOutParticles()[i-1]->SetMomentum(moms[i]);
  DEBUG_VAR(*blob);
  return amps;
}

Cluster_Amplitude* Decay_Handler_Base::ClusterConfiguration(Blob *const bl)
{
  msg_Debugging()<<METHOD<<"() {\n";
  msg_Indent();
  p_ampl = Cluster_Amplitude::New();
  p_ampl->SetMS(this);
  for (int i(0);i<bl->NInP();++i) {
    Particle *p(bl->InParticle(i));
    ColorID col(p->GetFlow(2),p->GetFlow(1));
    p_ampl->CreateLeg(-p->Momentum(),p->Flav().Bar(),col,1<<i);
  }
  p_ampl->SetNIn(bl->NInP());
  for (int i(0);i<bl->NOutP();++i) {
    Particle *p(bl->OutParticle(i));
    ColorID col(p->GetFlow(1),p->GetFlow(2));
    p_ampl->CreateLeg(p->Momentum(),p->Flav(),col,1<<(i+p_ampl->NIn()));
  }
  while (p_ampl->Legs().size()>p_ampl->NIn()+2) {
    msg_Debugging()<<*p_ampl<<"\n";
    Cluster_Amplitude *ampl(p_ampl);
    p_ampl = p_ampl->InitNext();
    p_ampl->SetMS(this);
    for (size_t i(0);i<ampl->NIn();++i) {
      Cluster_Leg *cl(ampl->Leg(i));
      p_ampl->CreateLeg(cl->Mom(),cl->Flav(),cl->Col(),cl->Id());
    }
    p_ampl->SetNIn(ampl->NIn());
    Cluster_Leg *lij(NULL);
    for (size_t i(ampl->NIn());i<ampl->Legs().size()-1;++i) {
      Cluster_Leg *li(ampl->Leg(i));
      for (size_t j(i+1);j<ampl->Legs().size();++j) {
        Cluster_Leg *lj(ampl->Leg(j));
        ColorID nc;
        if (li->Col().m_i==0 && li->Col().m_j==0) {
          nc=lj->Col();
        }
        else if (lj->Col().m_i==0 && lj->Col().m_j==0) {
          nc=li->Col();
        }
        else if (li->Col().m_i && li->Col().m_i==lj->Col().m_j) {
          nc.m_i=lj->Col().m_i;
          nc.m_j=li->Col().m_j;
        }
        else if (li->Col().m_j && li->Col().m_j==lj->Col().m_i) {
          nc.m_i=li->Col().m_i;
          nc.m_j=lj->Col().m_j;
        }
        if (nc.m_i>=0 && nc.m_j>=0) {
          Flavour fl(kf_photon);
          if (nc.m_i && nc.m_j) fl=Flavour(kf_gluon);
          else if (nc.m_i) fl=Flavour(kf_d);
          else if (nc.m_j) fl=Flavour(kf_d).Bar();
          p_ampl->CreateLeg(li->Mom()+lj->Mom(),fl,nc,li->Id()+lj->Id());
          lij=p_ampl->Legs().back();
          break;
        }
      }
      if (lij) break;
    }
    if (lij==NULL) THROW(fatal_error,"Internal eror");
    for (size_t i(ampl->NIn());i<ampl->Legs().size();++i) {
      Cluster_Leg *cl(ampl->Leg(i));
      if (cl->Id()&lij->Id()) continue;
      p_ampl->CreateLeg(cl->Mom(),cl->Flav(),cl->Col(),cl->Id());
    }
  }
  double mu2=p_ampl->Leg(0)->Mom().Abs2();
  p_ampl->SetMuF2(mu2);
  p_ampl->SetKT2(mu2);
  msg_Debugging()<<*p_ampl<<"\n";
  while (p_ampl->Prev()) {
    p_ampl=p_ampl->Prev();
    p_ampl->SetMuF2(mu2);
    p_ampl->SetKT2(mu2);
  }
  msg_Debugging()<<"}\n";
  return p_ampl;
}

void Decay_Handler_Base::CleanUp()
{
  if (p_decaymap) p_decaymap->ResetCounters();
}

bool Decay_Handler_Base::Decays(const ATOOLS::Flavour& flav)
{
  if (!flav.IsOn()) return false;
  if (flav.IsStable()) return false;
  return true;
}

bool Decay_Handler_Base::CanDecay(const ATOOLS::Flavour& flav)
{
  if (p_decaymap) return p_decaymap->Knows(flav);
  else return false;
}
