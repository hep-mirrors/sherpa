#ifndef Simple_XS_H
#define Simple_XS_H

#include "XS_Group.H"

namespace ATOOLS { class Data_Read;    }
namespace MODEL  { class Model_Base;   }
namespace PDF    { class Remnant_Base; }

namespace EXTRAXS {

  class Simple_XS: public XS_Group {
  private :

    std::string m_path, m_file, m_respath, m_resfile;

    size_t m_nmax, m_count, m_minqcdjet, m_maxqcdjet, m_maxjet;

    ATOOLS::Data_Read *p_dataread;

    std::string         m_regulator;
    std::vector<double> m_regulation;

    PDF::Remnant_Base *p_remnants[2];

    int ExtractFlavours(std::vector<ATOOLS::Flavour> &flavours,
			std::string buffer);

    void InitializeProcess(ATOOLS::Flavour *flavs,std::string &efunc,
			   std::string &printgraphs,bool psmc,
			   double &inisum,double &finsum,
			   size_t order_ew,size_t order_strong,
			   size_t nin,size_t nout,
			   const PHASIC::cls::scheme &clsc,
			   const PHASIC::hls::scheme &hlsc);

  public :

    // constructor
    Simple_XS(const std::string &path,const std::string &file,
	      MODEL::Model_Base *const model);

    // destructor
    ~Simple_XS();

    // member functions
    bool InitializeProcesses(BEAM::Beam_Spectra_Handler *const beamhandler,
			     PDF::ISR_Handler *const isrhandler,
			     const bool construct=true);
    bool CalculateTotalXSec(const std::string &resultpath);

    bool SelectOne();
    ATOOLS::Blob_Data_Base *OneEvent(const int mode=0);

    void SingleEvents();
    bool PrepareXSecTables();
    bool LookUpXSec(double dummy1,bool dummy2,std::string dummy3);

    // inline functions
    inline const std::string &ProcessName() const 
    { return p_selected->Selected()->Name(); }

    inline size_t NIn() const  { return p_selected->Selected()->NIn();  }
    inline size_t NOut() const { return p_selected->Selected()->NOut(); }

    inline const ATOOLS::Flavour *Flavours() const 
    { return p_selected->Selected()->Flavours(); }
    inline const ATOOLS::Vec4D   *Momenta() const  
    { return p_selected->Selected()->Momenta();  }
    
    inline int **const Colours() const 
    { return dynamic_cast<XS_Base*>(p_selected->Selected())->Colours();  }

    inline bool InSwaped() const { return p_selected->Selected()->InSwaped(); }

    inline size_t MinQCDJets() const { return m_minqcdjet; }
    inline size_t MaxQCDJets() const { return m_maxqcdjet; }

    inline size_t MaxJets() const { return m_maxjet; }

    void ResetSelector(ATOOLS::Selector_Data *const selectordata);

  }; // end of class Simple_XS

} // end of namespace EXTRAXS

#endif
