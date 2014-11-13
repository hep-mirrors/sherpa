#include "SHERPA/Single_Events/Signal_Process_FS_QED_Correction.H"

#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Phys/Momenta_Stretcher.H"
#include "ATOOLS/Phys/Particle.H"
#include "MODEL/Main/Model_Base.H"
#include "MODEL/Interaction_Models/Single_Vertex.H"


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
  m_on(true), m_qed(true),
  p_mehandler(_mehandler), p_sphotons(_sphotons)
{
  DEBUG_FUNC("");
  m_name      = string("Lepton_FS_QED_Corrections:");
  m_type      = eph::Perturbative;
  // general switch
  Data_Reader reader(" ",";","!","=");
  reader.AddComment("#");
  reader.AddWordSeparator("\t");
  reader.SetInputFile(rpa->gen.Variable("ME_DATA_FILE"));
  bool impliciteon = (reader.GetValue<std::string>("ME_QED","On")=="On");
  bool expliciteon = (reader.GetValue<std::string>("ME_QED","")=="On");
  // look whether there is any hadronisation following
  // if not, do not even put them on-shell -> switch everthing off
  if (!impliciteon) {
    m_qed = false;
    Data_Reader reader1(" ",";","!","=");
    reader1.AddComment("#");
    reader1.AddWordSeparator("\t");
    reader1.SetInputFile(rpa->gen.Variable("FRAGMENTATION_DATA_FILE"));
    m_on = (reader1.GetValue<std::string>("FRAGMENTATION","")!="Off");
  }
  // switch off if there are hard decays, have their own QED corrections,
  // cannot tell here what has been corrected and what not -- OR --
  // if NLO_Mode Fixed_Order, switch off completely, unless explicitely stated
  if (!expliciteon &&
      (p_mehandler->HasNLO()==1 ||
       !(reader.GetValue<std::string>("HARD_DECAYS","Off")=="Off"))) {
    m_on = false; m_qed = false;
  }
  msg_Debugging()<<"on="<<m_on<<" ,  qed="<<m_qed<<std::endl;

  if (m_on && m_qed) m_name += p_sphotons->SoftQEDGenerator();
  else               m_name += "None";
}

Signal_Process_FS_QED_Correction::~Signal_Process_FS_QED_Correction() {}


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

void Signal_Process_FS_QED_Correction::CleanUp(const size_t & mode) {}

void Signal_Process_FS_QED_Correction::Finish(const std::string &) {}

