#ifndef XS_Group_H
#define XS_Group_H

#include "XS_Base.H"

#include <vector>

namespace EXTRAXS {

  class XS_Selector; 

  class XS_Group: public XS_Base {
  protected:

    std::vector<XS_Base *> m_xsecs;

    bool m_atoms, m_channels, m_foundown;

    XS_Selector *p_xsselector;

    ATOOLS::Selector_Data *p_selectordata;

    XS_Base *Matching(const std::string &name);
    void     PrepareTerminate();

  public :


    // constructors
    XS_Group(); 
    XS_Group(const size_t nin,const size_t nout,const std::string &name);
    XS_Group(const size_t nin,const size_t nout,const ATOOLS::Flavour *flavours);
    XS_Group(const size_t nin,const size_t nout,const ATOOLS::Flavour *flavours,
	     const int scalescheme,const int kfactorscheme,
	     BEAM::Beam_Spectra_Handler *const beamhandler,
	     PDF::ISR_Handler *const isrhandler,
	     ATOOLS::Selector_Data *const selectordata);

    // destructor
    ~XS_Group();

    // member functions
    void Add(XS_Base *const xsec);
    bool Remove(XS_Base *const xsec);
    bool Delete(XS_Base *const xsec);
    void Clear();

    virtual void SetAtoms(const bool atoms);
    virtual void CreateFSRChannels();
    virtual void CreateISRChannels();

    size_t Size() const;

    XS_Base *const operator[](const size_t i) const;

    double Differential2();

    void SetISR(PDF::ISR_Handler *const isrhandler);

    void SetMax(const double max,const int flag);
    void SetTotal();
    void Reset();

    double Differential(const double s,const double t,const double u);
    bool   SetColours(const double s,const double t,const double u);
    double operator()(const double s,const double t,const double u);
    void   AddPoint(const double xs);
    void   ResetMax(int);
    void   OptimizeResult();

    bool CalculateTotalXSec(const std::string &resultpath,
			    const bool create=false);
    void WriteOutXSecs(std::ofstream &outfile);
    bool SelectOne();
    bool SelectOneFromList();
    void DeSelect();
    bool ReSelect(int);
    void SetEvents(const double number);

    void SetISRThreshold(const double threshold);
    void AddEvent(const double xs,const double validxs,const int ncounts);
    void ResetEvents();

    bool OneEvent();

    ATOOLS::Blob_Data_Base *WeightedEvent(const int mode=0);

    void SetPSHandler(PHASIC::Phase_Space_Handler *const pshandler); 
    void ResetSelector(ATOOLS::Selector_Data *const selectordata);

    void Print();

    // inline functions 
    inline void SetSelectorData(ATOOLS::Selector_Data *const data) 
    { p_selectordata=data; }

    inline void PushBack(XS_Base *const xsec) { m_xsecs.push_back(xsec); }

    inline ATOOLS::Selector_Data *const SelectorData() const
    { return p_selectordata; }

    inline XS_Selector *const XSSelector() const { return p_xsselector; }

  };// end of class XS_Group

}// end of namespace EXTRAXS

#endif


