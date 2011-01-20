#ifndef MODEL_Main_Spectrum_Generator_Base_H
#define MODEL_Main_Spectrum_Generator_Base_H

#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Data_Reader.H"

namespace ATOOLS { class Decay_Table; }

namespace MODEL {
  class Spectrum_Generator_Base {
  protected:
    ATOOLS::Data_Reader * p_dataread;
    Model_Base          * p_model;

    std::vector<double>          m_masses;
    std::vector<double>          m_widths;
    std::vector<ATOOLS::Flavour> m_particles;
    std::vector<ATOOLS::Decay_Table *> m_decays;
  public:
    Spectrum_Generator_Base(ATOOLS::Data_Reader *,Model_Base *);
    virtual ~Spectrum_Generator_Base();
    virtual void Run(PDF::ISR_Handler *const isr=NULL) = 0;
  };

  /*!
    \class Spectrum_Generator_Base
    \brief is the base class for all spectrum generators to be linked.

    Spectrum_Generator_Base is the base class used for wrappers around various
    spectrum generators, like the LesHouches interface.
    If these codes allow for it, the base class is not only useful to 
    store particles and their masses, but also widths and branching ratios can
    be stored in corresponding vectors. However, at the moment this option is not 
    fully exploited.
  */
  /*!    
    \var Model_Base * Spectrum_Generator_Base::p_model
    A pointer to the model class that owns the spectrum geenrator - this allows for
    simple transfer of the data already contained in the model class.
  */
  /*!
    ATOOLS::Data_Read * Spectrum_Generator_Base::p_dataread
    A pointer to the model's Data_Reader object. It is handed over through the constructor.
    This might become important, if more information is needed.
  */
  /*!
    std::vector<double> Spectrum_Generator_Base::m_widths
    A vector containing all widths that have been calculated, this works only in conjunction 
    with the particle vector.
  */
  /*!
    std::vector<double> Spectrum_Generator_Base::m_masses
    A vector containing all masses that have been calculated, this works only in conjunction 
    with the particle vector.
  */
  /*!
    std::vector<ATOOLS::Flavour> Spectrum_Generator_Base::m_particles
    A vector containing all particles that have been dealt with.
  */
  /*!
    std::vector<Decay_Table *> Spectrum_Generator_Base::m_decays
    A vector with all information about caluclated decays.
    \todo Work on the decays, they have to be filled into the vector of decayproducts,
          couplings have to be extracted, branching ratios transmitted to the rest of the code.
  */
  /*!
    Spectrum_Generator_Base::Spectrum_Generator_Base(ATOOLS::Data_Read *,Model_Base *)
    This constructor fills the pointers to the ATOOLS::Data_Read object and to the model
    that owns the spectrum generator.
  */
  /*!
    Spectrum_Generator_Base::Spectrum_Generator_Base::~Spectrum_Generator_Base()
    In the destructor, the list of Decay_Channels and its components is deleted
  */
  /*!
    virtual void Spectrum_Generator_Base::Run(std::string)
    Empty method, to be filled by the corresponding explicit wrapper.
  */
  /*!
    virtual void Spectrum_Generator_Base::FillMasses()
    Empty method, to be filled by the corresponding explicit wrapper.
  */
}

#endif


