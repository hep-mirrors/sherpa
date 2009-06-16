#ifndef Analysis_Triggers_Trigger_Base_H
#define Analysis_Triggers_Trigger_Base_H

#include "AddOns/Analysis/Main/Analysis_Object.H"

#include "ATOOLS/Phys/Particle_List.H"

namespace ANALYSIS {

  class Trigger_Base: public Analysis_Object {
  protected:

    std::string m_inlist, m_outlist;

  public:

    Trigger_Base(const std::string &inlist,const std::string &outlist);

    virtual void Evaluate(const ATOOLS::Particle_List &inlist,
			  ATOOLS::Particle_List &outlist,
			  double weight, double ncount) = 0;

    void Evaluate(const ATOOLS::Blob_List &bl, 
		  double weight, double ncount);

  };// end of class Trigger_Base

  class Two_List_Trigger_Base: public Analysis_Object {
  protected:

    std::string m_inlist, m_reflist, m_outlist;

  public:
    
    Two_List_Trigger_Base(const std::string &inlist,
			  const std::string &reflist,
			  const std::string &outlist);

    virtual void Evaluate(const ATOOLS::Particle_List &inlist,
			  const ATOOLS::Particle_List &reflist,
			  ATOOLS::Particle_List &outlist,
			  double weight, double ncount) = 0;

    void Evaluate(const ATOOLS::Blob_List &bl, 
		  double weight, double ncount);

  };// end of class Two_List_Trigger_Base

  class N_List_Trigger_Base: public Analysis_Object {
  protected:

    std::vector<std::string> m_inlists;

    std::string m_outlist;

  public:
    
    N_List_Trigger_Base(const std::vector<std::string> &inlists,
			const std::string &outlist);

    virtual void Evaluate
    (const std::vector<const ATOOLS::Particle_List*> &inlists,
     ATOOLS::Particle_List &outlist,double weight, double ncount) = 0;

    void Evaluate(const ATOOLS::Blob_List &bl, 
		  double weight, double ncount);

  };// end of class N_List_Trigger_Base

}// end of namespace ANALYSIS

#endif
