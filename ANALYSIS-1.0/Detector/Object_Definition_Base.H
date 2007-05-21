#ifndef Object_Definition_Base_h
#define Object_Definition_Base_h

#include "Tracker.H"
#include "ElMag_Calorimeter.H"
#include "Had_Calorimeter.H"
#include "Muon_Chambers.H"
#include "Reconstructed_Object.H"
#include "Analysis_Object.H"
#include "Particle_List.H"
#include "Particle.H"
#include "Ordering.H"
#include <map>

namespace ANALYSIS {
  class Primitive_Analysis;

  class Object_Definition_Base : public Analysis_Object {
  protected:
    bool   m_elements;
    std::string m_name;

    ATOOLS::kf::code m_kfcode;

    ATOOLS::Particle_List * p_myparticles;
    ATOOLS::Order_Base    * p_order;

    Tracker               * p_tracker;
    ElMag_Calorimeter     * p_ECal;
    Had_Calorimeter       * p_HCal;
    Muon_Chambers         * p_chambers;

    ObjectList              m_objects;
    virtual void GetElements();
    virtual void DropUsedCells();
    virtual void DropUsedTracks();
  public:
    Object_Definition_Base(Primitive_Analysis *,
			   const std::string name,const std::string="ET_UP");
    virtual ~Object_Definition_Base();

    virtual void Evaluate(const ATOOLS::Blob_List &,double,int) {}
    virtual Analysis_Object *GetCopy() const { return NULL; }

    virtual void ReconstructObjects(ATOOLS::Particle_List *,ATOOLS::Vec4D & METvector)=0;
    virtual void Reset();

    virtual ATOOLS::Particle_List * GetParticles() const { return p_myparticles; }
    virtual std::string Name() const { return m_name; }

    virtual void SetOrdering(const std::string);
  };


  class odb_sorter {
  public:
    bool operator() (const Object_Definition_Base * a,
		     const Object_Definition_Base * b) {
      if (a->Name()=="ElectronMaker") return true;
      if (b->Name()=="ElectronMaker") return false;
      if (a->Name()=="PhotonMaker")   return true;
      if (b->Name()=="PhotonMaker")   return false;
      if (a->Name()=="MuonMaker")     return true;
      if (b->Name()=="MuonMaker")     return false;
      if (a->Name()=="JetMaker")      return true;
      if (b->Name()=="JetMaker")      return false;
      if (a->Name()=="METMaker")      return true;
      return false;
    };
  };

}

#endif
