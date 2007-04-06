#ifndef Trigger_Base_H
#define Trigger_Base_H

#include "Analysis_Object.H"

#include "Particle_List.H"

namespace ANALYSIS {

  class Trigger_Base: public Analysis_Object {
  protected:

    std::string m_inlist, m_outlist;

  public:

    Trigger_Base(const std::string &inlist,const std::string &outlist);

    virtual void Evaluate(const ATOOLS::Particle_List &inlist,
			  ATOOLS::Particle_List &outlist,
			  double weight, int ncount) = 0;

    void Evaluate(const ATOOLS::Blob_List &bl, 
		  double weight, int ncount);

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
			  double weight, int ncount) = 0;

    void Evaluate(const ATOOLS::Blob_List &bl, 
		  double weight, int ncount);

  };// end of class Two_List_Trigger_Base

}// end of namespace ANALYSIS

#endif
