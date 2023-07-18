#ifndef SHERPA_PerturbativePhysics_Perturbative_Interface_H
#define SHERPA_PerturbativePhysics_Perturbative_Interface_H

#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Phys/Weights.H"
#include "ATOOLS/Org/Return_Value.H"

namespace ATOOLS {

  class Cluster_Amplitude;
}

namespace REMNANTS {
  class Remnant_Handler;
}

namespace SHERPA {
  class Matrix_Element_Handler;
  class Hard_Decay_Handler;
  class Decay_Handler_Base;
  class MI_Handler;
  class Soft_Collision_Handler;
  class Shower_Handler;

  struct mets_bbar_mode {
    enum code {
      none         = 0,
      enabled      = 1,
      lowestmulti  = 2,
      exclcluster  = 4
    };
  };

  inline mets_bbar_mode::code operator|(const mets_bbar_mode::code mm1,
                                        const mets_bbar_mode::code mm2)
  { return (mets_bbar_mode::code)((int)mm1|(int)mm2); }
  inline const mets_bbar_mode::code &operator|=(mets_bbar_mode::code &mm1,
                                                const mets_bbar_mode::code mm2)
  { return mm1=(mets_bbar_mode::code)((int)mm1|(int)mm2); }
  inline mets_bbar_mode::code operator&(const mets_bbar_mode::code mm1,
                                        const mets_bbar_mode::code mm2)
  { return (mets_bbar_mode::code)((int)mm1&(int)mm2); }
  inline const mets_bbar_mode::code &operator&=(mets_bbar_mode::code &mm1,
                                                const mets_bbar_mode::code mm2)
  { return mm1=(mets_bbar_mode::code)((int)mm1&(int)mm2); }
  inline mets_bbar_mode::code operator^(const mets_bbar_mode::code mm1,
                                        const mets_bbar_mode::code mm2)
  { return (mets_bbar_mode::code)((int)mm1^(int)mm2); }
  inline const mets_bbar_mode::code &operator^=(mets_bbar_mode::code &mm1,
                                                const mets_bbar_mode::code mm2)
  { return mm1=(mets_bbar_mode::code)((int)mm1^(int)mm2); }

  std::ostream & operator<<(std::ostream & s,const mets_bbar_mode::code & mm);
  std::istream & operator>>(std::istream & s,mets_bbar_mode::code &mm);

  class Perturbative_Interface {
  protected :

    Matrix_Element_Handler    * p_me;
    Hard_Decay_Handler        * p_dec;
    MI_Handler                * p_mi;
    Decay_Handler_Base        * p_hd;
    Soft_Collision_Handler    * p_sc;
    Shower_Handler            * p_shower;        
    REMNANTS::Remnant_Handler * p_remnants;

    ATOOLS::Blob *p_hard;
    ATOOLS::Blob_List *p_bloblist;

    ATOOLS::Poincare m_cms;

    ATOOLS::Cluster_Amplitude *p_ampl;

    ATOOLS::Weights_Map m_weightsmap, m_lkfweightsmap;
    double m_globalkfac, m_maxkfac;
    mets_bbar_mode::code m_bbarmode;

    long int m_fails_Moms, m_fails_Ampls, m_fails_Masses;
    bool LocalKFactor(ATOOLS::Cluster_Amplitude* ampl);

  public :

    Perturbative_Interface(Matrix_Element_Handler* const,
                           Hard_Decay_Handler* const,
                           Shower_Handler* const);
    Perturbative_Interface(Decay_Handler_Base *const hdh,
			   Shower_Handler *const psh);
    Perturbative_Interface(MI_Handler *const hdh,
			   Shower_Handler *const psh);
    Perturbative_Interface(Soft_Collision_Handler *const sch,
			   Shower_Handler *const psh);

    ~Perturbative_Interface();

    ATOOLS::Return_Value::code DefineInitialConditions(ATOOLS::Blob*,
                                                       ATOOLS::Blob_List*);

    bool FillBlobs();

    int  PerformShowers();
    int  PerformDecayShowers();

    void CleanUp();

    inline void SetRemnantHandler(REMNANTS::Remnant_Handler * remnants) { p_remnants = remnants; }

    inline ATOOLS::Weights_Map WeightsMap() const { return m_weightsmap; }

    inline Matrix_Element_Handler    * MEHandler()      const { return p_me; }
    inline ATOOLS::Cluster_Amplitude * Amplitude()      const { return p_ampl; }
    inline Shower_Handler            * Shower()               { return p_shower; }
    inline MI_Handler                * MIHandler()      const { return p_mi; }
    inline REMNANTS::Remnant_Handler * RemnantHandler() const { return p_remnants; }
  };// end of class Perturbative_Interface

  typedef std::map<std::string,Perturbative_Interface *> PertInterfaceMap;
  typedef PertInterfaceMap::iterator                     PertInterfaceIter;

}// end of namespace SHERPA

#endif
