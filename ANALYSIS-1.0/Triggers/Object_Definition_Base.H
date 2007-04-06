#ifndef Object_Definition_Base_h
#define Object_Definition_Base_h

#include "Analysis_Object_Base.H"
#include "Ordering.H"
#include "Particle_List.H"
#include "Getter_Function.H"

namespace ANALYSIS {
  class Object_Definition_Data {
  private:
    ATOOLS::Order_Base  * p_order;
    ATOOLS::Particle_List m_particles;
  public:
    Object_Definition_Data(std::string);
    ~Object_Definition_Data();
    void ResetPList();
    void AddPToPList(ATOOLS::Particle * const);
    void SortPList();

    ATOOLS::Particle_List & GetParticleList() { return m_particles; }
  };

  class Object_Definition_Base : public Analysis_Object_Base {
  protected:
    Object_Definition_Data * p_data;
  public:
    Object_Definition_Base() {}
    virtual ~Object_Definition_Base();
    virtual void Evaluate(const ATOOLS::Particle_List & plist,double weight, int ncount)=0;
    virtual void Evaluate(const ATOOLS::Blob_List & blobs, double weight, int ncount)=0;
  };

  typedef ATOOLS::Getter_Function<Object_Definition_Base,std::string> Object_Definition_Getter;
}
#endif
