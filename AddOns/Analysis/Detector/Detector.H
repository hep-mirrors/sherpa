#ifndef Analysis_Detector_Detector_H
#define Analysis_Detector_Detector_H

#include "AddOns/Analysis/Main/Analysis_Object.H"
#include "AddOns/Analysis/Main/Primitive_Analysis.H"
#include "AddOns/Analysis/Detector/Detector_Element.H"
#include "AddOns/Analysis/Detector/Particle_Smearer_Base.H"
#include "AddOns/Analysis/Detector/Object_Definition_Base.H"
#include "ATOOLS/Phys/Particle_List.H"

namespace ANALYSIS {
  class Detector : public Analysis_Object {
  private:
    Primitive_Analysis * p_ana;
    std::string          m_inlist, m_outlist;
    ATOOLS::Vec4D        m_METvector;

    std::map<std::string,Detector_Element *>        m_elements;
    std::map<std::string,Particle_Smearer_Base *>   m_smearers;
    std::set<Object_Definition_Base *,odb_sorter>   m_definitions;

    bool m_test;
    std::map<std::string,ATOOLS::Histogram *> m_histograms;

    void                    TestRandomIsotropicEvent(const int);
    void                    TestIsotropicEvent(const int,const int,const ATOOLS::Flavour);
    void                    TestReconstructionCodes(const int,const int,
						    const ATOOLS::Flavour,const bool);
    ATOOLS::Particle_List * ProduceParticleList(const std::string,int,const double,ATOOLS::Flavour);
    void                    InitHistograms(const int mode);
    void                    FillHistograms(const int j,const ATOOLS::Flavour);
    double                  Energy(const int);
  public:
    Detector(Primitive_Analysis * =NULL);
    ~Detector();

    void Evaluate(const ATOOLS::Blob_List &,double,double);
    Analysis_Object * GetCopy() const { return NULL; }

    void Fill(ATOOLS::Particle_List *);
    void ReconstructObjects(ATOOLS::Particle_List *&,ATOOLS::Vec4D & METvector);
    void Reset();

    bool IsActive() const { return (m_elements.size()>0); }

    void                     AddDetectorElement(Detector_Element *);
    Detector_Element       * GetElement(std::string);
    void                     AddParticleSmearer(Particle_Smearer_Base *);
    Particle_Smearer_Base  * GetParticleSmearer(std::string);
    void                     AddObjectDefinition(Object_Definition_Base *);
    Object_Definition_Base * GetObjectDefinition(std::string);

    void Test(const int mode=0);
    void Print();
  };
}

#endif
