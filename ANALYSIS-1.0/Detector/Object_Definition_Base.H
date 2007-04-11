#ifndef Object_Definition_Base_h
#define Object_Definition_Base_h

#include "Analysis_Object.H"
#include "Ordering.H"
#include "Particle_List.H"
#include "Getter_Function.H"



namespace ANALYSIS {
  class Primitive_Detector;

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

  class Object_Definition_Base : public Analysis_Object {
  protected:
    std::string              m_name;
    ATOOLS::kf::code         m_code;
    Object_Definition_Data * p_data;
  public:
    Object_Definition_Base() {}
    Object_Definition_Base(const std::string,const ATOOLS::kf::code,const std::string);
    virtual ~Object_Definition_Base();
    virtual Analysis_Object * GetCopy() { return NULL; }

    virtual void Evaluate(const ATOOLS::Particle_List &,double, int) {};
    virtual void Evaluate(const ATOOLS::Blob_List &, double, int) {};

    virtual void FillMCTruthList(ATOOLS::Particle_List *,Primitive_Detector *)=0;
    virtual void FillSimpleDetectorList(ATOOLS::Particle_List *,Primitive_Detector *)=0;

    std::string Name() { return m_name; }
  };

  typedef ATOOLS::Getter_Function<Object_Definition_Base,std::string> Object_Definition_Getter;

  typedef std::map<std::string,Object_Definition_Base *> String_ObjectDefinition_Map;
  typedef String_ObjectDefinition_Map::iterator          String_ObjectDefinition_Iter;
}
#endif
