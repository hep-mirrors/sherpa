#ifndef Primitive_Detector_H
#define Primitive_Detector_H

#include "Primitive_Observable_Base.H"
#include "Primitive_Detector_Element.H"
#include "Object_Definition_Base.H"
#include "Particle_List.H"

namespace ANALYSIS {
  struct det_mode {
    enum code {
      MC_truth    = 0,
      simple_det  = 10,
      object_reco = 20
    };
  };

  class Primitive_Detector: public Primitive_Observable_Base {
    det_mode::code              m_mode;
    std::string                 m_inlistname;
    String_DetectorElement_Map  m_elements;
    String_ObjectDefinition_Map m_objectdefs;
    ATOOLS::Particle_List     * m_internalparticles;

    void ProduceMCTruthObjectLists();
    void ProduceSimpleDetectorObjectLists();
    void ProduceReconstructedObjectLists();
    void FillDetector();
  public:
    Primitive_Detector(const det_mode::code,const std::string &inlist="FinalState");
    ~Primitive_Detector();

    void                         Print();
    void                         SetName(std::string);
    std::string                  Name() const;
    void                         AddDetectorElement(Primitive_Detector_Element *);
    Primitive_Detector_Element * GetElement(std::string);
    void                         AddObjectDefinition(Object_Definition_Base *);
    Object_Definition_Base     * GetObjectDefinition(std::string);


    Primitive_Observable_Base* Copy() const;
    
    void Evaluate(const ATOOLS::Blob_List &,double=1.,int=1);
    
    void Reset();
  };
}

#endif
