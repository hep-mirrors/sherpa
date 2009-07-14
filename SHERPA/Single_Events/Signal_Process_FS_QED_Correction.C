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
using namespace std;

Signal_Process_FS_QED_Correction::Signal_Process_FS_QED_Correction
(MEHandlersMap *_mehandlers, Soft_Photon_Handler *_sphotons) :
  p_mehandlers(_mehandlers), p_sphotons(_sphotons)
{
  m_name      = string("Lepton_FS_QED_Corrections:")
                  +p_sphotons->SoftQEDGenerator();
  m_type      = eph::Perturbative;
  // TODO: extract the resonance information from process infos known to MEHs

  // general switch
  Data_Reader reader(" ",";","!","=");
  std::string on = reader.GetValue<std::string>("ME_QED","off");
  m_on = (on=="on")?true:false;
  PRINT_VAR(m_on);
}

Signal_Process_FS_QED_Correction::~Signal_Process_FS_QED_Correction() {}


Return_Value::code Signal_Process_FS_QED_Correction::Treat
(Blob_List * bloblist, double & weight)
{
  if (!m_on) return Return_Value::Nothing;
  if (bloblist->empty()) {
    msg_Error()<<"Signal_Process_FS_QED_Correction::Treat("<<bloblist<<","<<weight<<"): "<<endl
               <<"   Blob list contains "<<bloblist->size()<<" entries."<<endl
               <<"   Continue and hope for the best."<<endl;
    return Return_Value::Error;
  }
  // look for hard process
  Blob_List sigblobs(bloblist->Find(btp::Signal_Process));
  Blob * sigblob(NULL);
  if (sigblobs.size()==1) sigblob=sigblobs[0];
  else {
    msg_Error()<<"Signal_Process_FS_QED_Correction::Treat("<<bloblist<<","<<weight<<"): "<<endl
               <<"  Blob list contains "<<sigblobs.size()<<" signal blobs."<<endl
               <<"  Continue and hope for the best."<<endl;
    return Return_Value::Error;
  }
  // extract FS leptons
  // two vectors -> the ones from the blob and the ones to be massive
  Particle_Vector fslep(sigblob->GetOutParticles());
  Particle_Vector mfslep;
  for (Particle_Vector::iterator it=fslep.begin();it!=fslep.end();++it) {
    if (!((*it)->Flav().IsLepton())) fslep.erase(it);
    else mfslep.push_back(new Particle(**it));
  }
  // put them on-shell (spoils consistency of pertubative calculation,
  // but necessary for YFS)
  if (!PutOnMassShell(mfslep)) {
    msg_Error()<<"Signal_Process_FS_QED_Correction::Treat("<<bloblist<<","<<weight<<"): "<<endl
               <<"  Leptons could not be put on their mass shell."<<endl
               <<"  Continue and hope for the best."<<endl;
    return Return_Value::Error;
  }
  // build effective vertex for resonant production
  Blob * effblob = new Blob();
  Particle * resonance = new Particle(0,DetermineResonanceFlavour(mfslep),
                                      MomentumSum(mfslep),'R');
  effblob->AddToInParticles(resonance);
  for (Particle_Vector::iterator it=mfslep.begin();it!=mfslep.end();++it) {
    effblob->AddToOutParticles(new Particle(**it));
  }
  effblob->SetStatus(blob_status::needs_extraQED);
  // add radiation
  if (!p_sphotons->AddRadiation(effblob)) {
    msg_Error()<<"Signal_Process_FS_QED_Correction::Treat("<<bloblist<<","<<weight<<"): "<<endl
               <<"  Higher order QED corrections failed."<<endl
               <<"  Continue and hope for the best."<<endl;
    return Return_Value::Error;
  }
  // build new QED radiation blob
  Blob * QEDblob = bloblist->AddBlob(btp::Shower);
  QEDblob->SetTypeSpec("YFS-type QED Corrections to ME");
  for (Particle_Vector::iterator it=fslep.begin();it!=fslep.end();++it) {
    QEDblob->AddToInParticles(*it);
  }
  for (size_t i=0;i<effblob->GetOutParticles().size();++i) {
    QEDblob->AddToOutParticles(effblob->GetOutParticles()[i]);
  }
//   msg_Out()<<*bloblist<<endl;
  delete resonance;
  return Return_Value::Nothing;
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
  // if two leptons of same flavour, take Z for now
  if ((partvec.size()==2) &&
      (partvec[0]->Flav()==partvec[1]->Flav().Bar())) {
    return Flavour(kf_Z);
  }
  // if lepton and corresponding neutrino, take W for now
  else if ((partvec.size()==2) &&
          (partvec[0]->Flav().LeptonFamily()==partvec[1]->Flav().LeptonFamily())) {
    if (partvec[0]->Flav().Charge()+partvec[1]->Flav().Charge())
      return Flavour(kf_Wplus);
    else
      return Flavour(kf_Wplus).Bar();
  }
  // TODO: take resonance from resonantly specified processes
  //       -> use process infos from MEHs
  // TODO: invent some placeholder otherwise
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


void Signal_Process_FS_QED_Correction::CleanUp() {}

void Signal_Process_FS_QED_Correction::Finish(const std::string &) {}
