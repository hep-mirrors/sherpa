#include "CS_Shower.H"

using namespace CS_SHOWER;
using namespace ATOOLS;
using namespace std;

CS_Shower::CS_Shower() : p_allsinglets(new All_Singlets) 
{
}

CS_Shower::~CS_Shower() 
{
  if (p_allsinglets) { delete p_allsinglets; p_allsinglets = NULL; }
}



bool CS_Shower::PerformShowers(All_Singlets * allsinglets) {
  if (allsinglets!=NULL) {
    if (p_allsinglets->size()>0) {
      for (ASiter asit=p_allsinglets->begin();asit!=p_allsinglets->end();asit++) {
	Singlet * sing = (*asit);
	if (sing) { delete sing; sing = NULL; }
      }
      p_allsinglets->clear();
    }
    p_allsinglets = allsinglets;
  }
  return m_shower.EvolveShower(p_allsinglets);
}

bool CS_Shower::ExtractPartons(ATOOLS::Blob_List * bloblist) {
  for (deque<Blob*>::iterator blit=bloblist->begin();blit!=bloblist->end();blit++) {
    if ((*blit)->Type()==btp::Shower) {
      ExtractPartons((*blit));
      return true;
    }
  }
  return false;
}

void CS_Shower::ExtractPartons(ATOOLS::Blob * blob) {
  for (ASiter asit=p_allsinglets->begin();asit!=p_allsinglets->end();asit++) {
    (*asit)->ExtractFSPartons(blob);
  }
}

void CS_Shower::PrepareAllSinglets()
{
  for (ASiter asit=p_allsinglets->begin();asit!=p_allsinglets->end();asit++) {
    Singlet * sing = (*asit);
    if (sing) { delete sing; sing = NULL; }
  }
  p_allsinglets->clear();
}
