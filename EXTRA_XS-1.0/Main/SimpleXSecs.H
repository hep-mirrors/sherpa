#ifndef SimpleXSecs_h
#define SimpleXSecs_h

#include "XS_Group.H"
#include "XS_Base.H"

#include "Model_Base.H"
#include "Beam_Spectra_Handler.H"
#include "ISR_Handler.H"

#include "Selector.H"
#include "Flavour.H"

#include "Data_Read.H"


namespace EXTRAXS {
  class SimpleXSecs : public XS_Group {
  private :
    std::string                  m_path,m_file;
    std::string                  m_respath,m_resfile;

    int                          m_nmax,m_count,m_maxjet;

    AORGTOOLS::Data_Read       * p_dataread;
    APHYTOOLS::Selector_Data   * p_seldata;
    BEAM::Beam_Spectra_Handler * p_beam;
    PDF::ISR_Handler           * p_isr;
  public :
    SimpleXSecs(std::string,std::string,MODEL::Model_Base * _model);
    ~SimpleXSecs();
    bool InitializeProcesses(BEAM::Beam_Spectra_Handler *,PDF::ISR_Handler *);
    bool CalculateTotalXSec();

    void        SelectOne();
    bool        UnweightedEvent(int mode=0);
    double      WeightedEvent(int mode=0);
    bool        SameEvent()    ;
    double      SameWeightedEvent();

    std::string          ProcessName() { return p_selected->Selected()->Name();    }
    int                  Nin()         { return p_selected->Selected()->Nin();     }
    int                  Nout()        { return p_selected->Selected()->Nout();    }
    int               ** Colours()     { return p_selected->Selected()->Colours(); }
    APHYTOOLS::Flavour * Flavs()       { return p_selected->Selected()->Flavs();   }
    AMATOOLS::Vec4D    * Momenta()     { return p_selected->Selected()->Momenta(); }
    int                  MaxJets()     { return m_maxjet;                          }

    void SingleEvents();
    bool PrepareXSecTables();
    bool LookUpXSec(double,bool,std::string);

  };
}

#endif
