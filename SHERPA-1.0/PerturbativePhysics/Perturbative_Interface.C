#include "Perturbative_Interface.H"

using namespace SHERPA;

Perturbative_Interface::Perturbative_Interface(Matrix_Element_Handler *mehandler,
					       Shower_Handler *shower):
  p_mehandler(mehandler), 
  p_shower(shower), 
  m_ini(shower->ISROn()), 
  m_fin(shower->FSROn()), 
  m_maxjetnumber(mehandler->MaxJets()), 
  m_weight(1.),
  p_fl(NULL), 
  p_moms(NULL) {}

Perturbative_Interface::~Perturbative_Interface() 
{
  if (p_fl)   { delete p_fl;   p_fl   = NULL; }
  if (p_moms) { delete p_moms; p_moms = NULL; }
}

void Perturbative_Interface::DeleteMergingBlobs(ATOOLS::Blob *const blob,
						const bool forward)
{ 
  for (size_t i=0;!forward && i<(size_t)blob->NInP();++i) {
    ATOOLS::Particle *cur=blob->InParticle(i);
    if (cur->ProductionBlob()!=NULL) DeleteMergingBlobs(cur->ProductionBlob(),false);
  } 
  for (size_t i=0;forward && i<(size_t)blob->NOutP();++i) {
    ATOOLS::Particle *cur=blob->OutParticle(i);
    if (cur->DecayBlob()!=NULL) DeleteMergingBlobs(cur->DecayBlob(),true);
  } 
}

void Perturbative_Interface::CleanBlobList(ATOOLS::Blob_List *const bloblist,
					   const ATOOLS::btp::code type)
{
  for (Blob_Iterator blit=bloblist->begin();blit!=bloblist->end();++blit) {
    if ((*blit)->Type()==type) {
      DeleteMergingBlobs(*blit,true);
      DeleteMergingBlobs(*blit,false);
      blit=bloblist->begin();
    }
  }
}

void Perturbative_Interface::Reset()
{
}

void Perturbative_Interface::UpdateEnergy(const double energy,const size_t i)
{
}
