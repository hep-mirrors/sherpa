#include "SHERPA/Single_Events/Soft_Photons.H"

#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Scoped_Settings.H"
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


Soft_Photons::Soft_Photons(
  Matrix_Element_Handler *_mehandler) :
  p_mehandler(_mehandler), p_yfshandler(p_mehandler->GetYFS())
{

  DEBUG_FUNC("");
  m_on = p_yfshandler->GetMode();
  m_name      = string("YFS_Lepton_QED_Corrections:");
  m_type      = eph::Perturbative;
  if (m_on) m_name += "YFS";
  else               m_name += "None";
}

Soft_Photons::~Soft_Photons(){

}


ATOOLS::Return_Value::code  Soft_Photons::Treat(Blob_List* bloblist) {
  if (!m_on) return Return_Value::Nothing;
  if (bloblist->empty()) {
    msg_Error() << METHOD
                << "(" << bloblist << "): " << endl
                << "   Blob list contains " << bloblist->size() << " entries." << endl
                << "   Continue and hope for the best." << endl;
    return Return_Value::Error;
  }
  // Blob * QEDblob;

  // if(bloblist->FindFirst(btp::QED_Radiation)){
  //   // QEDblob = bloblist->FindFirst(btp::QED_Radiation);
  //   bloblist->Delete(bloblist->FindFirst(btp::QED_Radiation));
  // }
  // QEDblob = bloblist->AddBlob(btp::QED_Radiation);
  // QEDblob->Reset();
  // QEDblob->SetTypeSpec("YFS-type_QED_Corrections_for_Lepton_collision");
  Blob * sigblob(bloblist->FindFirst(btp::Signal_Process));
  if (!sigblob) return Return_Value::Nothing;
  sigblob->AddStatus(blob_status::needs_yfs);
  // sigblob->SetStatus(blob_status::needs_reconnections | blob_status::needs_hadronization);

  // sigblob->AddStatus(blob_status::needs_showers);
  // sigblob->AddStatus(blob_status::needs_extraQED);
  // if(!p_yfshandler->MakeYFS()) return Return_Value::New_Event;
  // Particle_Vector islep(sigblob->GetInParticles());
  // Particle_Vector fslep(sigblob->GetOutParticles());
  // Particle_Vector mfslep, mislep;

  // for (Particle_Vector::iterator it=fslep.begin();it!=fslep.end();) {
  //     mfslep.push_back(new Particle(-1,(*it)->Flav(),(*it)->Momentum(),'F'));
  //     (*mfslep.rbegin())->SetNumber(0);
  //     (*mfslep.rbegin())->SetOriginalPart(*it);
  //     // (*mfslep.rbegin())->SetFinalMass((*it)->FinalMass());
  //     ++it;
  // }

  // for (Particle_Vector::iterator it=islep.begin();it!=islep.end();) {
  //     mislep.push_back(new Particle(-1,(*it)->Flav(),(*it)->Momentum(),'I'));
  //     (*mislep.rbegin())->SetNumber(0);
  //     (*mislep.rbegin())->SetOriginalPart(*it);
  //     // (*mislep.rbegin())->SetFinalMass((*it)->FinalMass());
  //     ++it;
  // }

  int Nin=2;
  int Nout = p_mehandler->Process()->NOut();
  int N=Nin+Nout;
  // if(QEDblob->GetOutParticles().size()>Nout){
  //   for (int i = Nout; i < QEDblob->GetOutParticles().size(); ++i)
  //   {
  //    QEDblob->RemoveOutParticle(QEDblob->GetOutParticles()[i]);
  //   }
  // }
  // if(QEDblob->GetInParticles().size()>Nin){
  //   for (int i = Nin; i < QEDblob->GetInParticles().size(); ++i)
  //   {
  //     QEDblob->RemoveInParticle(QEDblob->GetInParticles()[i]);
  //   }
  // }

  // for (Particle_Vector::iterator it=mislep.begin();it!=mislep.end();++it) {
  //   // check for momentum conservation does not work
  //   (*it)->SetInfo('H');
  //   (*it)->SetStatus(part_status::active);
  //   QEDblob->AddToInParticles(*it);
  // }

  // for (Particle_Vector::iterator it=mfslep.begin();it!=mfslep.end();++it) {
  // //   // check for momentum conservation does not work
  //   (*it)->SetInfo('H');
  //   (*it)->SetStatus(part_status::active);
  //   QEDblob->AddToOutParticles(*it);
  // }
  
  if(p_yfshandler->GetFSRMode()!=0){
    // Add the fsr corrected final states
      Particle_Vector out = sigblob->GetOutParticles();
      Particle_Vector yfsout = p_yfshandler->m_particles;
      ATOOLS::ParticleMomMap yfsoutMap = p_yfshandler->m_outparticles;
      if(out.size()!=(yfsout.size()-2)){
        msg_Error()<<METHOD<<" Missmatch in outparitcles for YFS"<<std::endl
                            <<"Born Out size = "<< out.size()<<std::endl
                            <<"YFS Out size = "<< yfsout.size()<<std::endl;
      }
      for(int i=0; i<out.size(); i++){
        sigblob->OutParticle(i)->SetMomentum(yfsoutMap[yfsout[i+2]]); // remove born momenta
      }
    }
  if (p_yfshandler->GetMode() != 0) {
    // blob->SetStatus(blob_status::needs_yfs);
    ATOOLS::Vec4D_Vector isrphotons = p_yfshandler->GetISRPhotons();
    ATOOLS::Vec4D_Vector fsrphotons;
    Particle *particle;
    if (p_yfshandler->GetFSRMode() != 0) {
      fsrphotons = p_yfshandler->GetFSRPhotons();
    }
    if (p_yfshandler->FillBlob()) {
      for (int i = 0; i < isrphotons.size(); ++i)
      {
        particle = new Particle(-1, Flavour(22),
                                isrphotons[i]);
        particle->SetNumber(0);
        particle->SetInfo('S');
        sigblob->AddToOutParticles(particle);
      }
      for (int i = 0; i < fsrphotons.size(); ++i)
      {
        particle = new Particle(-1, Flavour(22),
                                fsrphotons[i]);
        particle->SetNumber(0);
        particle->SetInfo('S');
        sigblob->AddToOutParticles(particle);
      }
    }
  }
  sigblob->CheckMomentumConservation();
  // QEDblob->SetStatus(blob_status::needs_reconnections | blob_status::needs_hadronization);
  return Return_Value::Nothing;
}



void Soft_Photons::CleanUp(const size_t & mode) {}

void Soft_Photons::Finish(const std::string &) {}

