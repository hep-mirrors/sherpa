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
}

#endif


