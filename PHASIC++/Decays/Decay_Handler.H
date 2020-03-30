#ifndef Decay_Handler_h
#define Decay_Handler_h

#include "ATOOLS/Phys/Flavour.H"

namespace ATOOLS {
  class Blob_List;
  class Blob;
  class Particle;
  class Cluster_Amplitude;
  template<typename Scalar>
  class Lorentz_Ten2;
  typedef Lorentz_Ten2<double>  Lorentz_Ten2D;
}

namespace METOOLS {
  class Spin_Density;
  class Decay_Matrix;
  class Amplitude2_Tensor;
}

namespace PHASIC {

  class Decay_Map;

  /**
   * The Mass_Selector knows which type of mass is to be used for all decays
   * associated with this decay map. This can be the HadMass() if the decays
   * belong to the non-perturbative regime of the event, or the Mass() if
   * they belong to the perturbative part. 
   */
  class Decay_Handler : public ATOOLS::Mass_Selector {
  protected:
    PHASIC::Decay_Map* p_decaymap;

  public :
    /**
     * Constructor
     */
    Decay_Handler();
    
    /**
     * Destructor
     */
    virtual ~Decay_Handler();

    virtual void CreateDecayBlob(ATOOLS::Blob_List* bloblist, ATOOLS::Particle* inpart) = 0;
    virtual METOOLS::Amplitude2_Tensor* FillOnshellDecay(ATOOLS::Blob* blob, METOOLS::Spin_Density* sigma);
    virtual bool VetoDecayAndPrepForNew(ATOOLS::Blob* blob) { return false; }
    virtual void AfterTreatInitialBlob(ATOOLS::Blob* blob, ATOOLS::Blob_List* bloblist) {}
    virtual void BeforeFillDecayTree(ATOOLS::Blob * blob) {}

    virtual double DecayWidth(const ATOOLS::Flavour& flav);
    virtual bool DiceMass(ATOOLS::Particle* part, double max);

    virtual ATOOLS::Flavour PropFlav(ATOOLS::Blob* blob);
    
    virtual void CleanUp();
    virtual std::string Name() = 0;
  };
}

#endif
