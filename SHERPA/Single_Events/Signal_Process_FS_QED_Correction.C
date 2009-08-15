#include "SHERPA/Single_Events/Signal_Process_FS_QED_Correction.H"

#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Phys/Momenta_Stretcher.H"
#include "ATOOLS/Phys/Particle.H"


using namespace SHERPA;
using namespace ATOOLS;
using namespace PHASIC;
using namespace std;

////////////////////////////////////////////////////////////////////////////////
////                                                                        ////
////     in the documentation LEPTON is synonymous                          ////
////     to EVERYTHING NOT STRONGLY CHARGED                                 ////
////                                                                        ////
////////////////////////////////////////////////////////////////////////////////


Signal_Process_FS_QED_Correction::Signal_Process_FS_QED_Correction
(MEHandlersMap *_mehandlers, Soft_Photon_Handler *_sphotons) :
  p_mehandlers(_mehandlers), p_sphotons(_sphotons)
{
  m_name      = string("Lepton_FS_QED_Corrections:");
  m_type      = eph::Perturbative;
  // general switch
  Data_Reader reader(" ",";","!","=");
  std::string on = reader.GetValue<std::string>("ME_QED","On");
  m_on = (on=="On")?true:false;

  if (p_mehandlers->size()>1) {
    m_on = false;
    msg_Error()<<METHOD<<"() {\n"
               <<"  there are "<<p_mehandlers->size()<<" ME Handlers present ...\n"
               <<"  do not know how to continue ...\n"
               <<"  event generation will commence without QED corrections to leptons from hard scatter ...\n"
               <<"}\n";
  }

  Process_Vector pvec(p_mehandlers->begin()->second->AllProcesses());
  if (m_on) {
    m_name += p_sphotons->SoftQEDGenerator();
    msg_Debugging()<<METHOD<<"(){\n";
    for (size_t i=0;i<pvec.size();++i) {
      for (size_t j=0;j<pvec[i]->Size();++j) {
        SubInfoVector siv;
        FindSubProcessInfosContainingLeptons((*pvec[i])[j]->Info(),siv);
        msg_Debugging()<<"  found process: "<<(*pvec[i])[j]->Name()<<" with "
                       <<siv.size()<<" defined lepton production blobs...\n";
        for (size_t k=0;k<siv.size();++k)
          msg_Debugging()<<*siv[k]<<endl;
        m_proc_lep_map.insert(make_pair((*pvec[i])[j]->Name(),siv));
      }
    }
    msg_Debugging()<<"}\n";
  }
  else
    m_name += "None";
}

Signal_Process_FS_QED_Correction::~Signal_Process_FS_QED_Correction() {}


Return_Value::code Signal_Process_FS_QED_Correction::Treat
(Blob_List * bloblist, double & weight)
{
  if (bloblist->empty()) {
    msg_Error()<<"Signal_Process_FS_QED_Correction::Treat("<<bloblist<<","<<weight<<"): "<<endl
               <<"   Blob list contains "<<bloblist->size()<<" entries."<<endl
               <<"   Continue and hope for the best."<<endl;
    return Return_Value::Error;
  }
  // look for QCD corrected hard process in need for QED
  Blob * sigblob(bloblist->FindLast(btp::Shower));
  if (!sigblob) return Return_Value::Nothing;
  // if already treated -> nothing to do
  if (sigblob->TypeSpec()=="YFS-type QED Corrections to ME")
    return Return_Value::Nothing;
  if (sigblob->TypeSpec()=="setting leptons on-shell")
    return Return_Value::Nothing;
  // extract FS leptons
  // two vectors -> the ones from the blob and the ones to be massive
  Particle_Vector fslep(sigblob->GetOutParticles());
  Particle_Vector mfslep;
  for (Particle_Vector::iterator it=fslep.begin();it!=fslep.end();) {
    if ((*it)->Flav().Strong() || (*it)->DecayBlob()!=NULL) {
      fslep.erase(it);
    }
    else {
      mfslep.push_back(new Particle(-1,(*it)->Flav(),(*it)->Momentum(),'F'));
      (*mfslep.rbegin())->SetNumber(0);
      ++it;
    }
  }
  // if no leptons, nothing to do
  // if only one lepton, cannot do anything
  if (fslep.size()<2) {
    sigblob->UnsetStatus(blob_status::needs_extraQED);
    for (Particle_Vector::iterator it=mfslep.begin();it!=mfslep.end();++it)
      delete *it;
    return Return_Value::Nothing;
  }
  // put them on-shell (spoils consistency of pertubative calculation,
  // but necessary for YFS)
  if (!PutOnMassShell(mfslep)) {
    msg_Error()<<"Signal_Process_FS_QED_Correction::Treat("<<bloblist<<","<<weight<<"): "<<endl
               <<"  Leptons could not be put on their mass shell."<<endl
               <<"  Trying new event."<<endl;
    for (Particle_Vector::iterator it=mfslep.begin();it!=mfslep.end();++it)
      delete *it;
    return Return_Value::New_Event;
  }
  // if switched off or no need for QED stop here and build a blob
  if (!m_on || !sigblob->Has(blob_status::needs_extraQED)) {
    Blob * onshellblob = bloblist->AddBlob(btp::QED_Radiation);
    onshellblob->SetTypeSpec("setting leptons on-shell");
    if (sigblob->Has(blob_status::needs_extraQED))
      sigblob->UnsetStatus(blob_status::needs_extraQED);
    for (Particle_Vector::iterator it=fslep.begin();it!=fslep.end();++it) {
      (*it)->SetInfo('H');
      (*it)->SetStatus(part_status::decayed);
      onshellblob->AddToInParticles(*it);
    }
    for (Particle_Vector::iterator it=mfslep.begin();it!=mfslep.end();++it) {
      onshellblob->AddToOutParticles(*it);
    }
    onshellblob->SetStatus(blob_status::needs_hadronization);
    return Return_Value::Success;
  }
  // build effective verteces for resonant production
  // use subprocess infos if possible
  Blob_Vector blobs = BuildResonantBlobs(mfslep);
  // add radiation
  for (Blob_Vector::iterator it=blobs.begin();it!=blobs.end();++it) {
    // do nothing if no resonance determined
    if ((*it)->InParticle(0)->Flav().Kfcode()!=kf_none) {
      (*it)->SetStatus(blob_status::needs_extraQED);
      if (!p_sphotons->AddRadiation(*it)) {
        msg_Error()<<"Signal_Process_FS_QED_Correction::Treat("<<bloblist<<","<<weight<<"): "<<endl
                   <<"  Higher order QED corrections failed."<<endl
                   <<"  Retrying event."<<endl;
        for (Particle_Vector::iterator it=mfslep.begin();it!=mfslep.end();++it)
          delete *it;
        for (Blob_Vector::iterator it=blobs.begin();it!=blobs.end();++it)
          delete *it;
        return Return_Value::Retry_Event;
      }
    }
  }
  for (Blob_Vector::iterator it=blobs.begin();it!=blobs.end();++it) {
    msg_Debugging()<<**it<<endl;
    (*it)->DeleteInParticles();
  }
  sigblob->UnsetStatus(blob_status::needs_extraQED);
  // build new QED radiation blob
  Blob * QEDblob = bloblist->AddBlob(btp::QED_Radiation);
  QEDblob->SetTypeSpec("YFS-type QED Corrections to ME");
  for (Particle_Vector::iterator it=fslep.begin();it!=fslep.end();++it) {
    // set info back to hard process, otherwise
    // check for momentum conservation does not work
    (*it)->SetInfo('H');
    (*it)->SetStatus(part_status::decayed);
    QEDblob->AddToInParticles(*it);
  }
  for (Blob_Vector::iterator it=blobs.begin();it!=blobs.end();++it) {
    while ((*it)->NOutP()) {
      Particle * part((*it)->RemoveOutParticle((*it)->NOutP()-1,true));
      QEDblob->AddToOutParticles(part);
    }
  }
  QEDblob->SetStatus(blob_status::needs_hadronization);
  for (size_t i=0;i<blobs.size();++i) {
    delete blobs[i];
    blobs[i]=NULL;
  }
  return Return_Value::Success;
}

bool Signal_Process_FS_QED_Correction::PutOnMassShell(const Particle_Vector& partvec)
{
  // if massless in ME put on mass shell for YFS
  bool allonshell(true);
  std::vector<double>masses(partvec.size(),0.);
  for (size_t i=0;i<partvec.size();++i) {
    masses[i]=partvec[i]->Flav().Mass(1);
    if (!IsEqual(partvec[i]->Momentum().Abs2(),sqr(masses[i]),1E-4))
      allonshell=false;
  }
  if (allonshell) return true;
  Momenta_Stretcher momstretch;
  return momstretch.StretchMomenta(partvec,masses);
}

Flavour Signal_Process_FS_QED_Correction::DetermineResonanceFlavour
(const Particle_Vector& partvec)
{
  msg_Debugging()<<"determining resonance flavour on the basis of:"<<endl;
  for (size_t j=0;j<partvec.size();++j) msg_Debugging()<<*partvec[j]<<endl;
  SubInfoVector siv;
  FindSubProcessInfosContainingLeptons
  (p_mehandlers->begin()->second->Process()->Info(),siv);
  // if two leptons of same flavour, take Z for now
  if ((partvec.size()==2) &&
      (partvec[0]->Flav()==partvec[1]->Flav().Bar())) {
    return Flavour(kf_Z);
  }
  // if lepton and corresponding neutrino, take W for now
  else if ((partvec.size()==2) &&
          (partvec[0]->Flav().LeptonFamily()==partvec[1]->Flav().LeptonFamily())) {
    if ((partvec[0]->Flav().Charge()+partvec[1]->Flav().Charge())==1.)
      return Flavour(kf_Wplus);
    else
      return Flavour(kf_Wplus).Bar();
  }
  // guess on basis of sum of charges
  else {
    int chargesum(0);
    for (size_t i=0;i<partvec.size();++i)
      chargesum+=partvec[i]->Flav().IntCharge();
    if (chargesum==0)  return Flavour(kf_Z);
    if (chargesum==3)  return Flavour(kf_Wplus);
    if (chargesum==-3) return Flavour(kf_Wplus).Bar();
  }
  // i got no clue what this might be
  return Flavour(kf_none);
}

Vec4D Signal_Process_FS_QED_Correction::MomentumSum
(const Particle_Vector& partvec)
{
  Vec4D sum(0.,0.,0.,0.);
  for (size_t i=0;i<partvec.size();++i) {
    sum += partvec[i]->Momentum();
  }
  return sum;
}

void Signal_Process_FS_QED_Correction::FindSubProcessInfosContainingLeptons
(const Process_Info& pi, SubInfoVector& siv)
{
  // loop over FS of process -> find subprocs, that are subsequent decays
  for (size_t i=0;i<pi.m_fi.m_ps.size();++i) {
    if (pi.m_fi.m_ps[i].m_ps.size()>1)
      FindSubProcessInfosContainingLeptons(pi.m_fi.m_ps[i],siv);
  }
}

void Signal_Process_FS_QED_Correction::FindSubProcessInfosContainingLeptons
(const Subprocess_Info& spi, SubInfoVector& siv)
{
  // assume connected leptons are produced in same subprocess info
  size_t count(0), leps(0);
  for (size_t i=0;i<spi.m_ps.size();++i) {
    if (spi.m_ps[i].m_ps.size()==0) {
      count++;
      if (!spi.m_ps[i].m_fl.Strong()) leps++;
    }
    else {
      FindSubProcessInfosContainingLeptons(spi.m_ps[i],siv);
    }
  }
  if (count==spi.m_ps.size() && leps!=0) {
    // -> final subprocess info
    // if contains leptons, add to siv
    siv.push_back(&spi);
  }
}

Blob_Vector Signal_Process_FS_QED_Correction::BuildResonantBlobs
(Particle_Vector& pv)
{
  // get production subprocesses for the active process
  std::string name(p_mehandlers->begin()->second->Process()->Name());
  SubInfoVector siv(m_proc_lep_map[name]);
  // create blobs accordingly (only if lepton list is unambiguous)
  Blob_Vector blobs;
  DEBUG_INFO("pv unambiguous: "<<ContainsNoAmbiguities(pv));
  DEBUG_INFO(siv.size()<<" subprocess infos for process "<<name);
  if (siv.size() && ContainsNoAmbiguities(pv)) {
    for (size_t i=0;i<siv.size();++i) {
      blobs.push_back(new Blob(Vec4D(0.,0.,0.,0.)));
      FillBlob(blobs[i],*siv[i],pv);
      msg_Debugging()<<"built blob:"<<endl;
      msg_Debugging()<<*blobs[i]<<endl;
    }
  }
  // otherwise create global resonant blob
  // if there are leptons not contained in defined resonant blobs
  if (pv.size()) {
    blobs.push_back(new Blob(Vec4D(0.,0.,0.,0.)));
    FillBlob(*blobs.rbegin(),DetermineResonanceFlavour(pv),pv);
    msg_Debugging()<<"built generic blob:"<<endl;
    msg_Debugging()<<**blobs.rbegin()<<endl;
  }
  return blobs;
}

bool Signal_Process_FS_QED_Correction::ContainsNoAmbiguities
(const Particle_Vector& pv)
{
  // look whether any particle occurs more than once
  std::set<std::string> checklist;
  for (size_t i=0;i<pv.size();++i) {
    if (checklist.find(pv[i]->Flav().IDName()) == checklist.end())
      checklist.insert(pv[i]->Flav().IDName());
    else return false;
  }
  return true;
}

void Signal_Process_FS_QED_Correction::FillBlob
(Blob * blob, const Subprocess_Info& spi, Particle_Vector& pv)
{
  // find the leptons owned by this subprocess
  Particle_Vector localpv;
  bool onlyleptons(true);
  for (size_t i=0;i<spi.m_ps.size();++i) {
    if (spi.m_ps[i].m_fl.Strong()) onlyleptons=false;
    for (Particle_Vector::iterator it=pv.begin();it!=pv.end();) {
      if ((*it)->Flav()==spi.m_ps[i].m_fl) {
        localpv.push_back(*it);
        pv.erase(it);
      }
      else ++it;
    }
  }
  if (onlyleptons) FillBlob(blob,spi.m_fl,localpv);
  else FillBlob(blob,DetermineResonanceFlavour(localpv),localpv);
}

void Signal_Process_FS_QED_Correction::FillBlob
(Blob * blob, const Flavour& resflav, Particle_Vector& pv)
{
  Vec4D sum(MomentumSum(pv));
  for (Particle_Vector::iterator it=pv.begin();it!=pv.end();) {
    blob->AddToOutParticles(*it);
    pv.erase(it);
  }
  blob->AddToInParticles(new Particle(-1,resflav,sum,'R'));
  msg_Debugging()<<"resonance is:\n"<<*blob->InParticle(0)<<endl;
  blob->InParticle(0)->SetFinalMass(blob->InParticle(0)->Momentum().Mass());
}

void Signal_Process_FS_QED_Correction::CleanUp() {}

void Signal_Process_FS_QED_Correction::Finish(const std::string &) {}

