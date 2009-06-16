#ifndef Analysis_Detector_Object_Definition_Base_H
#define Analysis_Detector_Object_Definition_Base_H

#include "AddOns/Analysis/Detector/Tracker.H"
#include "AddOns/Analysis/Detector/ElMag_Calorimeter.H"
#include "AddOns/Analysis/Detector/Had_Calorimeter.H"
#include "AddOns/Analysis/Detector/Muon_Chambers.H"
#include "AddOns/Analysis/Detector/Reconstructed_Object.H"
#include "AddOns/Analysis/Main/Analysis_Object.H"
#include "ATOOLS/Phys/Particle_List.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Phys/Ordering.H"
#include <map>

namespace ANALYSIS {
  class Primitive_Analysis;

  class Object_Definition_Base : public Analysis_Object {
  protected:
    bool   m_elements;
    std::string m_name;

    kf_code m_kfcode;

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

    virtual void Evaluate(const ATOOLS::Blob_List &,double,double) {}
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
