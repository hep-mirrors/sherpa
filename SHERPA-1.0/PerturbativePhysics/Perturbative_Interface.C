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

void Perturbative_Interface::CleanBlobs(ATOOLS::Blob_List * bl)
{
  for (Blob_Iterator blit=bl->begin();blit!=bl->end();) {
    if ((*blit)->Type()==btp::ME_PS_Interface_FS || 
	(*blit)->Type()==btp::ME_PS_Interface_IS) blit=bl->erase(blit);
    else ++blit;
  }
}

