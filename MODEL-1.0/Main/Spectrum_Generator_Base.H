#ifndef Spectrum_Generator_Base_h
#define Spectrum_Generator_Base_h

#include "Model_Base.H"
#include "MathTools.H"
#include "Data_Read.H"

namespace MODEL {
  class DecayChannel {
  private :
    APHYTOOLS::Flavour              m_flin;
    std::vector<APHYTOOLS::Flavour> m_flouts;
    double                          m_width;
  public :
    DecayChannel(const APHYTOOLS::Flavour);
    void AddDecayProduct(const APHYTOOLS::Flavour);
    void SetDecayWidth(const double);
    void Output();
    int                Nout();
    APHYTOOLS::Flavour Flin();
    APHYTOOLS::Flavour Flout(const int i);
    double             DecayWidth();
  };


  class Spectrum_Generator_Base {
  protected:
    Model_Base           * p_model;
    AORGTOOLS::Data_Read * p_dataread;

    std::vector<double>             widths;
    std::vector<double>             masses;
    std::vector<APHYTOOLS::Flavour> particles;

    std::vector<DecayChannel *>     decayproducts;
  public:
    Spectrum_Generator_Base(AORGTOOLS::Data_Read *,Model_Base *);
    ~Spectrum_Generator_Base();
    virtual void Run(std::string) = 0;
    virtual void FillMasses()     = 0;
    virtual void FillDecays()     = 0;
  };

  /*!
    \class Spectrum_Generator_Base
    \brief is the base class for all spectrum geenrators to be linked.

    Spectrum_Generator_Base is the base class used for wrappers around various
    spectrum generators, so far the two Fortran codes Hdecay and Isasusy have been 
    linked. If these codes allow for it, the base class is not only useful to 
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
    AORGTOOLS::Data_Read * Spectrum_Generator_Base::p_dataread
    A pointer to the model's dataread object. It is handed over through the constructor.
    This might become important, if moer information is needed.
  */
  /*!
    std::vector<double> Spectrum_Generator_Base::widths
    A vector containing all widths that have been calculated, this works only in conjunction 
    with the particle vector.
  */
  /*!
    std::vector<double> Spectrum_Generator_Base::masses
    A vector containing all masses that have been calculated, this works only in conjunction 
    with the particle vector.
  */
  /*!
    std::vector<APHYTOOLS::Flavour> Spectrum_Generator_Base::particles
    A vector containing all particles that have been dealt with.
  */
  /*!
    std::vector<DecayChannel *> Spectrum_Generator_Base::decayproducts
    A vector with all information about caluclated decays. Presumably, more work has to here.
    \todo Work on the decays, they have to be filled into the vector of decayproducts,
          couplnigs have to be extracted, branching ratios transmitted to the rest of the code.
  */
  /*!
    Spectrum_Generator_Base::Spectrum_Generator_Base(AORGTOOLS::Data_Read *,Model_Base *)
    This constructor fills the pointers to the AORGTOOLS::Data_Read object and to the model
    that owns the spectrum generator.
  */
  /*!
    Spectrum_Generator_Base::Spectrum_Generator_Base::~Spectrum_Generator_Base()
    In the destructor, the list of DecayChannels and its components is deleted
  */
  /*!
    virtual void Spectrum_Generator_Base::Run(std::string)
    Empty method, to be filled by the corresponding explicit wrapper.
  */
  /*!
    virtual void Spectrum_Generator_Base::FillMasses()
    Empty method, to be filled by the corresponding explicit wrapper.
  */
  /*!
    virtual void Spectrum_Generator_Base::FillDecays()
    Empty method, to be filled by the corresponding explicit wrapper.
  */
}

#endif


