#include "AddOns/Analysis/Triggers/Trigger_Base.H"

#include "AddOns/Analysis/Main/Primitive_Analysis.H"
#include "ATOOLS/Org/Message.H"

using namespace ANALYSIS;
using namespace ATOOLS;

Trigger_Base::Trigger_Base(const std::string &inlist,
			   const std::string &outlist):
  m_inlist(inlist), m_outlist(outlist) {}

void Trigger_Base::Evaluate(const ATOOLS::Blob_List &bl, 
			    double weight, double ncount)
{
  Particle_List *outlist(new Particle_List());
  Particle_List *inlist(p_ana->GetParticleList(m_inlist));
  if (inlist==NULL) {
    msg_Error()<<METHOD<<"["<<typeid(*this).name()
	       <<"](): List '"<<m_inlist<<"' not found."<<std::endl;
    p_ana->AddParticleList(m_outlist,outlist);
    return;
  }
  Evaluate(*inlist,*outlist,weight,ncount);
  p_ana->AddParticleList(m_outlist,outlist);
}

Two_List_Trigger_Base::Two_List_Trigger_Base
(const std::string &inlist,const std::string &reflist,
 const std::string &outlist):
  m_inlist(inlist), m_reflist(reflist), m_outlist(outlist) {}

void Two_List_Trigger_Base::Evaluate(const ATOOLS::Blob_List &bl, 
				     double weight, double ncount)
{
  Particle_List *outlist(new Particle_List());
  Particle_List *inlist(p_ana->GetParticleList(m_inlist));
  if (inlist==NULL) {
    msg_Error()<<METHOD<<"["<<typeid(*this).name()
	       <<"](): List '"<<m_inlist<<"' not found."<<std::endl;
    p_ana->AddParticleList(m_outlist,outlist);
    return;
  }
  Particle_List *reflist(p_ana->GetParticleList(m_reflist));
  if (reflist==NULL) {
    msg_Error()<<METHOD<<"["<<typeid(*this).name()
	       <<"](): List '"<<m_reflist<<"' not found."<<std::endl;
    p_ana->AddParticleList(m_outlist,outlist);
    return;
  }
  Evaluate(*inlist,*reflist,*outlist,weight,ncount);
  p_ana->AddParticleList(m_outlist,outlist);
}

N_List_Trigger_Base::N_List_Trigger_Base
(const std::vector<std::string> &inlists,const std::string &outlist):
  m_inlists(inlists), m_outlist(outlist) {}

void N_List_Trigger_Base::Evaluate(const ATOOLS::Blob_List &bl, 
				   double weight, double ncount)
{
  Particle_List *outlist(new Particle_List());
  std::vector<const Particle_List*> inlists(m_inlists.size());
  for (size_t i(0);i<m_inlists.size();++i) {
    inlists[i]=p_ana->GetParticleList(m_inlists[i]);
    if (inlists[i]==NULL) {
      msg_Error()<<METHOD<<"["<<typeid(*this).name()
		 <<"](): List '"<<m_inlists[i]<<"' not found."<<std::endl;
      p_ana->AddParticleList(m_outlist,outlist);
      return;
    }
  }
  Evaluate(inlists,*outlist,weight,ncount);
  p_ana->AddParticleList(m_outlist,outlist);
}


