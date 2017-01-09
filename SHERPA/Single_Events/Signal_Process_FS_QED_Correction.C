#include "SHERPA/Single_Events/Signal_Process_FS_QED_Correction.H"

#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Org/Default_Reader.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Phys/Momenta_Stretcher.H"
#include "ATOOLS/Phys/Particle.H"
#include "MODEL/Main/Model_Base.H"
#include "MODEL/Main/Single_Vertex.H"


using namespace SHERPA;
using namespace ATOOLS;
using namespace PHASIC;
using namespace MODEL;
using namespace std;

////////////////////////////////////////////////////////////////////////////////
////                                                                        ////
////     in the documentation LEPTON is synonymous                          ////
////     to EVERYTHING NOT STRONGLY CHARGED                                 ////
////                                                                        ////
////////////////////////////////////////////////////////////////////////////////


Signal_Process_FS_QED_Correction::Signal_Process_FS_QED_Correction
(Matrix_Element_Handler *_mehandler, Soft_Photon_Handler *_sphotons) :
  m_on(true), m_qed(true), m_forceqedonnloqcd(false),
  p_mehandler(_mehandler), p_sphotons(_sphotons), p_newsublist(NULL)
{
  DEBUG_FUNC("");
  m_name      = string("Lepton_FS_QED_Corrections:");
  m_type      = eph::Perturbative;
  // general switch
  Default_Reader reader;
  reader.SetInputFile(rpa->gen.Variable("ME_DATA_FILE"));
  bool impliciteon = (reader.Get<std::string>("ME_QED","On")=="On");
  bool expliciteon = (reader.Get<std::string>("ME_QED","")=="On");
  // look whether there is any hadronisation following
  // if not, do not even put them on-shell -> switch everthing off
  msg_Debugging()<<"impl="<<impliciteon<<", expl="<<expliciteon<<std::endl;
  if (!impliciteon) {
    m_qed = false;
    Default_Reader fragmentation_reader;
    fragmentation_reader.SetInputFile(rpa->gen.Variable("FRAGMENTATION_DATA_FILE"));
    string fragmentation_model(fragmentation_reader.GetStringNormalisingNoneLikeValues("FRAGMENTATION", ""));
    m_on = (fragmentation_model != "None");
  }
  // switch off if there are hard decays, have their own QED corrections,
  // cannot tell here what has been corrected and what not -- OR --
  // if NLO_Mode Fixed_Order, switch off completely, unless explicitely stated
  std::string hd(reader.GetStringNormalisingNoneLikeValues("HARD_DECAYS","None"));
  if (!expliciteon &&
      (p_mehandler->HasNLO() == 1 || !(hd == "None"))) {
    m_on = false; m_qed = false;
  }
  if (expliciteon && p_mehandler->HasNLO()==1) {
    m_forceqedonnloqcd = reader.Get<size_t>("ME_QED_ON_NLO_QCD",0);
  }
  msg_Debugging()<<"on="<<m_on<<" ,  qed="<<m_qed<<std::endl;
  msg_Debugging()<<"force on on nlo qcd: "<<m_forceqedonnloqcd<<std::endl;

  if (m_on && m_qed) m_name += p_sphotons->SoftQEDGenerator();
  else               m_name += "None";
}

Signal_Process_FS_QED_Correction::~Signal_Process_FS_QED_Correction()
{
  if (p_newsublist) { DeleteNewSubList(); delete p_newsublist; }

}


Return_Value::code Signal_Process_FS_QED_Correction::Treat
(Blob_List * bloblist, double & weight)
{
  if (!m_on) return Return_Value::Nothing;
  if (bloblist->empty()) {
    msg_Error()<<"Signal_Process_FS_QED_Correction::Treat"
	       <<"("<<bloblist<<","<<weight<<"): "<<endl
               <<"   Blob list contains "<<bloblist->size()<<" entries."<<endl
               <<"   Continue and hope for the best."<<endl;
    return Return_Value::Error;
  }
  // look for QCD corrected hard process in need for QED
  Blob * sigblob(bloblist->FindLast(btp::Shower));
  // if NLO QCD, no shower blob, take ME blob instead, status not set yet
  if (m_forceqedonnloqcd) {
    sigblob = bloblist->FindFirst(btp::Signal_Process);
    if (!sigblob) return Return_Value::Nothing;
    msg_Debugging()<<sigblob->TypeSpec()<<std::endl;
    sigblob->AddStatus(blob_status::needs_extraQED);
  }
  if (!sigblob) return Return_Value::Nothing;
  // if already treated -> nothing to do
  if (sigblob->TypeSpec()=="YFS-type_QED_Corrections_to_ME")
    return Return_Value::Nothing;
  if (sigblob->TypeSpec()=="setting_leptons_on-shell")
    return Return_Value::Nothing;
  // extract FS leptons
  // two vectors -> the ones from the blob and the ones to be massive
  DEBUG_FUNC(m_qed);
  Particle_Vector fslep(sigblob->GetOutParticles());
  Particle_Vector mfslep;
  for (Particle_Vector::iterator it=fslep.begin();it!=fslep.end();) {
    if ((*it)->Flav().Strong() || (*it)->Flav().IsDiQuark() || 
	(*it)->DecayBlob()!=NULL) {
      fslep.erase(it);
    }
    else {
      mfslep.push_back(new Particle(-1,(*it)->Flav(),(*it)->Momentum(),'F'));
      (*mfslep.rbegin())->SetNumber(0);
      (*mfslep.rbegin())->SetOriginalPart(*it);
      (*mfslep.rbegin())->SetFinalMass((*it)->FinalMass());
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
  // if switched off or no need for QED stop here and build a blob
  if (!m_qed || !sigblob->Has(blob_status::needs_extraQED)) {
    Blob * onshellblob = bloblist->AddBlob(btp::QED_Radiation);
    onshellblob->SetTypeSpec("setting_leptons_on-shell");
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
  // put them on-shell (spoils consistency of pertubative calculation,
  // but necessary for YFS)
  if (!PutOnMassShell(mfslep)) {
    msg_Error()<<"Signal_Process_FS_QED_Correction::Treat("
	       <<bloblist<<","<<weight<<"): \n"
               <<"  Leptons could not be put on their mass shell.\n"
               <<"  Trying new event.\n"
               <<"  The event contained a ";
    for (Particle_Vector::iterator it=mfslep.begin();it!=mfslep.end();++it)
       msg_Error()<<(*it)->Flav().ShellName()<<"-";
    if (mfslep.size()==2) msg_Error()<<"pair";
    else                  msg_Error()<<"set";
    msg_Error()<<" of too little invariant mass to be put\n"
	       <<"  on their mass shell. If you are sensitive to this specific"
	       <<" signature consider\n  to set the respective particles"
	       <<" massive in the perturbative calculation using\n"
	       <<"  'MASSIVE[<id>]=1' to avoid this problem.\n";
    for (Particle_Vector::iterator it=mfslep.begin();it!=mfslep.end();++it)
      delete *it;
    return Return_Value::New_Event;
  }
  // add radiation
  Blob_Vector blobs;
  if (!p_sphotons->AddRadiation(mfslep,blobs)) {
    msg_Error()<<"Signal_Process_FS_QED_Correction::Treat("<<bloblist
               <<","<<weight<<"): "<<endl
               <<"  Higher order QED corrections failed."<<endl
               <<"  Retrying event."<<endl;
    for (Particle_Vector::iterator it=mfslep.begin();it!=mfslep.end();++it)
      delete *it;
    for (Blob_Vector::iterator it=blobs.begin();it!=blobs.end();++it)
      delete *it;
    return Return_Value::Retry_Event;
  }
  sigblob->UnsetStatus(blob_status::needs_extraQED);
  // build new QED radiation blob
  Blob * QEDblob = bloblist->AddBlob(btp::QED_Radiation);
  QEDblob->SetTypeSpec("YFS-type_QED_Corrections_to_ME");
  for (Particle_Vector::iterator it=fslep.begin();it!=fslep.end();++it) {
    // set info back to hard process, otherwise
    // check for momentum conservation does not work
    (*it)->SetInfo('H');
    (*it)->SetStatus(part_status::decayed);
    QEDblob->AddToInParticles(*it);
  }
  // if fixed-order NLO RS and on, copy QED radiation also to subevtlist
  if (m_forceqedonnloqcd) ModifySubEvtList(sigblob,fslep,blobs);
  // first fill in all LO particles
  for (Blob_Vector::iterator it=blobs.begin();it!=blobs.end();++it) {
    while ((*it)->NOutP() && (*it)->OutParticle(0)->Info()!='S') {
      Particle * part((*it)->RemoveOutParticle(0,true));
      QEDblob->AddToOutParticles(part);
    }
  }
  // then append all photons
  for (Blob_Vector::iterator it=blobs.begin();it!=blobs.end();++it) {
    while ((*it)->NOutP()) {
      Particle * part((*it)->RemoveOutParticle(0,true));
      QEDblob->AddToOutParticles(part);
    }
  }
  QEDblob->SetStatus(blob_status::needs_hadronization);
  // clean up
  for (size_t i=0;i<blobs.size();++i) {
    delete blobs[i];
    blobs[i]=NULL;
  }
  return Return_Value::Success;
}

bool Signal_Process_FS_QED_Correction::PutOnMassShell
(const Particle_Vector& partvec)
{
  // if massless in ME put on mass shell for YFS
  bool allonshell(true); kf_code kfc;
  std::vector<double>masses(partvec.size(),0.);
  for (size_t i=0;i<partvec.size();++i) {
    kfc=partvec[i]->Flav().Kfcode();
    if(kfc==kf_graviton || kfc==kf_gscalar)
      masses[i]=sqrt(fabs(partvec[i]->Momentum().Abs2()));
    else masses[i]=partvec[i]->Flav().Mass(1);
    //If one of the two squared masses is zero, IsEqual always returns 0.
    if (!IsEqual(partvec[i]->Momentum().Abs2(),sqr(masses[i]),1E-4))
      allonshell=false;
  }
  if (allonshell) return true;
  Momenta_Stretcher momstretch;
  return momstretch.StretchMomenta(partvec,masses);
}

bool Signal_Process_FS_QED_Correction::ModifySubEvtList
(ATOOLS::Blob *sigblob, const ATOOLS::Particle_Vector& fslep,
 const ATOOLS::Blob_Vector& blobs)
{
  DEBUG_FUNC("");
  NLO_subevtlist* sublist(NULL);
  Blob_Data_Base * bdb((*sigblob)["NLO_subeventlist"]);
  if (bdb) sublist=bdb->Get<NLO_subevtlist*>();
  if (!sublist) return true;

  // QED corrections in blobs have to be transferred onto each subevt
  // en lieu of the fsleps

  if (p_newsublist) DeleteNewSubList();
  else p_newsublist=new NLO_subevtlist();

  // determine new size of final state for S (R has one more)
  size_t newnreal(2+sigblob->NOutP()-fslep.size());
  for (size_t i(0);i<blobs.size();++i) newnreal+=blobs[i]->NOutP();
  size_t newnsub(newnreal-1);
  msg_Debugging()<<"new n(sub)  = "<<newnsub<<std::endl;
  msg_Debugging()<<"new n(real) = "<<newnreal<<std::endl;

  // determine boost for S events
  Vec4D fslepmomR(0.,0.,0.,0.);
  for (size_t i(0);i<fslep.size();++i) fslepmomR+=fslep[i]->Momentum();
  msg_Debugging()<<"fslepmomR = "<<fslepmomR<<std::endl;
  Poincare oldframe(fslepmomR);

  // book-keep which post-rad fsleps have been filled, and photons
  Particle_Vector fslepQED,fsphotons;
  for (size_t i(0);i<blobs.size();++i)
    for (size_t j(0);j<blobs[i]->NOutP();++j)
      if (!blobs[i]->OutParticle(j)->Flav().IsPhoton())
        fslepQED.push_back(blobs[i]->OutParticle(j));
      else
        fsphotons.push_back(blobs[i]->OutParticle(j));
  if (fslep.size()!=fslepQED.size()) THROW(fatal_error,"Internal error.");

  // order fslepQED same as fslep -> can fill them in order
  for (size_t i(0);i<fslep.size();) {
    if (fslep[i]->Flav()==fslepQED[i]->Flav()) ++i;
    else {
      for (size_t j(0);j<fslep.size()-1;++j) {
        std::swap(fslepQED[j],fslepQED[j+1]);
      }
    }
  }
  msg_Debugging()<<fslep.size()<<" fsleps(QED): "<<std::endl;
  for (size_t i(0);i<fslep.size();++i)
    msg_Debugging()<<i<<": "<<fslep[i]->Flav()<<" -> "
                            <<fslepQED[i]->Flav()<<std::endl;

  for (size_t i(0);i<sublist->size();++i) {
    // iterate over sub events and replace radiating particles
    NLO_subevt* sub((*sublist)[i]);
    msg_Debugging()<<i<<": "<<*sub<<std::endl;
    for (size_t j(0);j<sub->m_n;++j)
      msg_Debugging()<<sub->p_fl[j]<<": "<<sub->p_mom[j]<<std::endl;

    // structure will be: first n entries as before, (newn-n) photons appended
    size_t newn(newnsub);
    if (sub->IsReal()) newn=newnreal;
    Flavour* newfls = new Flavour[newn];
    Vec4D* newmoms = new Vec4D[newn];
    size_t* newids = new size_t[newn];
    for (size_t n(0);n<newn;++n) newids[n]=0;
    NLO_subevt* newsub=new NLO_subevt(*sub);
    newsub->m_n=newn;
    newsub->m_delete=true;

    Vec4D fslepmomS(0.,0.,0.,0.);
    if (sub->IsReal()) fslepmomS=fslepmomR;
    else for (size_t j(0);j<sub->m_n;++j)
      if (!sub->p_fl[j].Strong()) fslepmomS+=sub->p_mom[j];
    msg_Debugging()<<"fslepmomS = "<<fslepmomS<<std::endl;
    Poincare newframe(fslepmomS);
    newframe.Invert();

    // replace leptons
    size_t idx(0);
    for (size_t n(0);n<sub->m_n;++n) {
      newids[n]=sub->p_id[n];
      newfls[n]=sub->p_fl[n];
      if (n<sub->Proc<PHASIC::Process_Base>()->NIn() || sub->p_fl[n].Strong()) {
        newmoms[n]=sub->p_mom[n];
      }
      else {
        if (sub->IsReal()) newmoms[n]=fslepQED[idx]->Momentum();
        else newmoms[n]=newframe*(oldframe*(fslepQED[idx]->Momentum()));
        ++idx;
      }
    }
    if (idx!=fslepQED.size()) THROW(fatal_error,"Not all leptons replaced.");

    //add photons
    for (size_t n(0);n<fsphotons.size();++n) {
      newids[sub->m_n+n]=0;
      newfls[sub->m_n+n]=Flavour(kf_photon);
      if (sub->IsReal()) newmoms[sub->m_n+n]=fsphotons[n]->Momentum();
      else newmoms[sub->m_n+n]=newframe*(oldframe*(fsphotons[n]->Momentum()));
    }

    newsub->p_id=newids;
    newsub->p_fl=newfls;
    newsub->p_mom=newmoms;
    p_newsublist->push_back(newsub);

    msg_Debugging()<<i<<": "<<*newsub<<std::endl;
    for (size_t j(0);j<newsub->m_n;++j)
      msg_Debugging()<<newsub->p_fl[j]<<": "<<newsub->p_mom[j]<<std::endl;
  }

  // set new sublist as evt-sublist
  bdb->Set<NLO_subevtlist*>(p_newsublist);

  return true;
}

void Signal_Process_FS_QED_Correction::DeleteNewSubList()
{
  for (size_t i=0; i<p_newsublist->size(); ++i) delete (*p_newsublist)[i];
  p_newsublist->clear();
}

void Signal_Process_FS_QED_Correction::CleanUp(const size_t & mode) {}

void Signal_Process_FS_QED_Correction::Finish(const std::string &) {}

