#include "AddOns/Analysis/Main/Analysis_Object.H"

#include "AddOns/Analysis/Tools/Particle_Qualifier.H"
#include "ATOOLS/Org/Message.H"

namespace ANALYSIS {

  class List_Creator: public Analysis_Object {
  private:

    std::string  m_outlist;

  public:

    List_Creator(const std::string &outlist);

    void CreatePrimordialHadronsList(const ATOOLS::Blob_List *bl);
    void CreateIntermediateHadronsList(const ATOOLS::Blob_List *bl);
    void CreateChargedParticleList(const ATOOLS::Blob_List *bl);
    void CreateUEParticleList(const ATOOLS::Blob_List *bl);

    void Evaluate(const ATOOLS::Blob_List & ,double weight, double ncount);

    Analysis_Object *GetCopy() const;

  };

} // namespace ANALYSIS

using namespace ANALYSIS;
using namespace ATOOLS;

DECLARE_GETTER(List_Creator_Getter,"CreateList",
 	       Analysis_Object,Argument_Matrix);

void List_Creator_Getter::PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"list";
}

Analysis_Object *
List_Creator_Getter::operator()(const Argument_Matrix &parameters) const
{
  std::string outlist("");
  if (parameters.size()!=1) return NULL;
  if (parameters[0].size()!=1) return NULL;
  if (parameters[0][0]!="PrimordialHadrons" &&
      parameters[0][0]!="IntermediateHadrons" &&
      parameters[0][0]!="ChargedParticle" &&
      parameters[0][0]!="UEPartons") return NULL;
  return new List_Creator(parameters[0][0]);
}

#include "AddOns/Analysis/Main/Primitive_Analysis.H"

List_Creator::List_Creator(const std::string &outlist): 
  m_outlist(outlist) {}

void List_Creator::Evaluate
(const ATOOLS::Blob_List &bl,double weight,double ncount)
{
  Particle_List *list(p_ana->GetParticleList(m_outlist,true));
  if (list!=NULL) return;
  if (m_outlist=="PrimordialHadrons") 
    CreatePrimordialHadronsList(&bl);
  else if (m_outlist=="IntermediateHadrons") 
    CreateIntermediateHadronsList(&bl);
  else if (m_outlist=="ChargedParticle") 
    CreateChargedParticleList(&bl);
  else if (m_outlist=="UEPartons") 
    CreateUEParticleList(&bl);
}

void List_Creator::CreatePrimordialHadronsList(const ATOOLS::Blob_List *bl)
{
  Particle_List *pl(new Particle_List);
  for (Blob_List::const_iterator blit(bl->begin());
       blit!=bl->end();++blit) {
    if ((*blit)->Type()==btp::Fragmentation) {
      for (int i=0;i<(*blit)->NOutP();++i) {
	Particle * p = (*blit)->OutParticle(i);
	if (p->Flav().IsHadron()) pl->push_back(p);
      }
    }
  }
  p_ana->AddParticleList("PrimordialHadrons",pl);
}

void List_Creator::CreateIntermediateHadronsList(const ATOOLS::Blob_List *bl)
{
  Particle_List *pl(new Particle_List);
  for (Blob_List::const_iterator blit(bl->begin());
       blit!=bl->end();++blit) {
    if (((*blit)->Type()==btp::Hadron_Decay || 
	 (*blit)->Type()==btp::Fragmentation) &&
	(*blit)->NOutP()>1) {
      for (int i=0;i<(*blit)->NOutP();++i) {
	Particle * p = (*blit)->OutParticle(i);
	if (p->Flav().IsHadron()) pl->push_back(p);
      }
    }
  }
  p_ana->AddParticleList("IntermediateHadrons",pl);
}

void List_Creator::CreateChargedParticleList(const ATOOLS::Blob_List *bl)
{
  Particle_List *pl_fs(p_ana->GetParticleList("FinalState"));
  if (pl_fs==NULL) {
    msg_Error()<<METHOD<<"(): Final state particle list not found."<<std::endl;
    return;
  }
  Particle_List *pl(new Particle_List);
  copy_if(pl_fs->begin(),pl_fs->end(),
	  back_inserter(*pl),Is_Charged());
  p_ana->AddParticleList("ChargedParticle",pl);
}


void List_Creator::CreateUEParticleList(const ATOOLS::Blob_List *bl)
{
  Particle_List *pl(new Particle_List);
  std::map<int,Blob*> bmap;
  for (Blob_List::const_iterator blit(bl->begin());
       blit!=bl->end();++blit) {
    if ((*blit)->Type()==btp::Hard_Collision) {
      bmap[(*blit)->Id()]=(*blit);
    }
  }
  int cnt;
  do {
    cnt=0;
    for (std::map<int,Blob*>::iterator mit(bmap.begin());
	 mit!=bmap.end();++mit) {
      Blob *blb(mit->second);
      if (blb->Type()!=btp::Shower) {
	for (int i=0;i<blb->NOutP();++i) {
	  std::map<int,Blob*>::iterator src
	    (bmap.find(blb->OutParticle(i)->DecayBlob()->Id()));
	  if (src==bmap.end()) {
	    cnt++;
	    bmap[blb->OutParticle(i)->DecayBlob()->Id()]=
	      blb->OutParticle(i)->DecayBlob();
	  }
	}
	bmap.erase(blb->Id());
      }
      if (cnt>0) break;
    }
  } while(cnt>0);   
  for (std::map<int,Blob*>::iterator mit(bmap.begin());
       mit!=bmap.end();++mit) {
    Blob *blb(mit->second);
    for (int i(0);i<blb->NOutP();++i) {
      Particle *p(blb->OutParticle(i));
      pl->push_back(p);	
    }
  }
  p_ana->AddParticleList("UEPartons",pl);
}


Analysis_Object *List_Creator::GetCopy() const
{
  return new List_Creator(m_outlist);
}
