#include "Trigger_Base.H"

#include "Primitive_Analysis.H"
#include "Message.H"

using namespace ANALYSIS;
using namespace ATOOLS;

Trigger_Base::Trigger_Base(const std::string &inlist,
			   const std::string &outlist):
  m_inlist(inlist), m_outlist(outlist) {}

void Trigger_Base::Evaluate(const ATOOLS::Blob_List &bl, 
			    double weight, int ncount)
{
  Particle_List *inlist(p_ana->GetParticleList(m_inlist));
  if (inlist==NULL) {
    msg.Error()<<METHOD<<"(): List '"<<m_inlist
		       <<"' not found."<<std::endl;
    return;
  }
  Particle_List *outlist(new Particle_List());
  Evaluate(*inlist,*outlist,weight,ncount);
  p_ana->AddParticleList(m_outlist,outlist);
}

Two_List_Trigger_Base::Two_List_Trigger_Base
(const std::string &inlist,const std::string &reflist,
 const std::string &outlist):
  m_inlist(inlist), m_reflist(reflist), m_outlist(outlist) {}

void Two_List_Trigger_Base::Evaluate(const ATOOLS::Blob_List &bl, 
				     double weight, int ncount)
{
  Particle_List *inlist(p_ana->GetParticleList(m_inlist));
  if (inlist==NULL) {
    msg.Error()<<METHOD<<"(): List '"<<m_inlist
		       <<"' not found."<<std::endl;
    return;
  }
  Particle_List *reflist(p_ana->GetParticleList(m_reflist));
  if (reflist==NULL) {
    msg.Error()<<METHOD<<"(): List '"<<m_reflist
	       <<"' not found."<<std::endl;
    return;
  }
  Particle_List *outlist(new Particle_List());
  Evaluate(*inlist,*reflist,*outlist,weight,ncount);
  p_ana->AddParticleList(m_outlist,outlist);
}


