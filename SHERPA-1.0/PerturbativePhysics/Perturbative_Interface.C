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

void Perturbative_Interface::RemoveBackward(ATOOLS::Blob *const blob,
					    ATOOLS::Blob_List *const bloblist)
{ 
  while (blob->NInP()>0) {
    ATOOLS::Particle *cur=blob->InParticle(0);
    if (cur->ProductionBlob()!=NULL) {
      ATOOLS::Blob *prod=cur->ProductionBlob();
      RemoveBackward(prod,bloblist);
      for (ATOOLS::Blob_List::iterator bit=bloblist->begin();
	   bit!=bloblist->end();++bit) {
	if (*bit==prod) {
	  bloblist->erase(bit);
	  break;
	}
      }
      delete prod;
    }
    else blob->DeleteInParticle(cur);
  } 
}

void Perturbative_Interface::RemoveForward(ATOOLS::Blob *const blob,
					    ATOOLS::Blob_List *const bloblist)
{
  while (blob->NOutP()>0) {
    ATOOLS::Particle *cur=blob->OutParticle(0);
    if (cur->DecayBlob()!=NULL) {
      ATOOLS::Blob *dec=cur->DecayBlob();
      RemoveForward(dec,bloblist);
      for (ATOOLS::Blob_List::iterator bit=bloblist->begin();
	   bit!=bloblist->end();++bit) {
	if (*bit==dec) {
	  bloblist->erase(bit);
	  break;
	}
      }
      delete dec;
    }
    else blob->DeleteOutParticle(cur);
  } 
}

void Perturbative_Interface::CleanBlobList(ATOOLS::Blob_List *const bloblist,
					   const ATOOLS::btp::code type)
{
  for (ATOOLS::Blob_Iterator blit=bloblist->begin();blit!=bloblist->end();++blit) {
    if ((*blit)->Type()==type) {
      ATOOLS::Blob *blob=*blit;
      RemoveForward(blob,bloblist);
      RemoveBackward(blob,bloblist);
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
