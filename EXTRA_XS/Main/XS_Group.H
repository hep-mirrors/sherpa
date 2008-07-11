#ifndef XS_Group_H
#define XS_Group_H

#include "XS_Base.H"

#include <vector>

namespace EXTRAXS {

  class XS_Selector; 

  class XS_Group: public XS_Base {
  protected:

    std::vector<XS_Base *> m_xsecs;

    bool m_atoms, m_foundown;

    XS_Selector *p_xsselector;

    ATOOLS::Selector_Data *p_selectordata;

    bool CheckFlavours(std::vector<ATOOLS::Flavour> &cfl);
    bool ConstructProcesses(const size_t &oew,const size_t &os,
			    std::vector<ATOOLS::Flavour> cfl,
			    const size_t &ci,const double &fixscale);

  public :


    // constructors
    XS_Group(); 
    XS_Group(const size_t nin,const size_t nout,const std::string &name);
    XS_Group(const size_t nin,const size_t nout,
	     const ATOOLS::Flavour *flavours,XS_Model_Base *const model);
    XS_Group(const size_t nin,const size_t nout,
	     const ATOOLS::Flavour *flavours,
	     const PHASIC::scl::scheme scalescheme,const int kfactorscheme,
	     BEAM::Beam_Spectra_Handler *const beamhandler,
	     PDF::ISR_Handler *const isrhandler,
	     ATOOLS::Selector_Data *const selectordata,
	     XS_Model_Base *const model);

    // destructor
    ~XS_Group();

    // member functions
    bool ConstructProcesses(const size_t &oew,const size_t &os,
			    const double &fixscale);
    
    void Add(XS_Base *const xsec);
    bool Remove(XS_Base *const xsec);
    bool Delete(XS_Base *const xsec);
    void Clear();

    virtual void SetAtoms(const bool atoms);

    void SetFactorizationScale(const std::string &muf2);
    void SetRenormalizationScale(const std::string &mur2);
    void SetCoreMaxJetNumber(const int &n);

    size_t Size() const;

    XS_Base *const operator[](const size_t i) const;

    double Differential2();

    void SetISR(PDF::ISR_Handler *const isrhandler);

    void SetMax(const double max,const int flag);
    void SetTotal(const int mode=1);
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
    void DeSelect();
    bool ReSelect(int);

    void SetISRThreshold(const double threshold);

    ATOOLS::Blob_Data_Base *OneEvent();
    ATOOLS::Blob_Data_Base *WeightedEvent(const int mode=0);

    void SetPSHandler(PHASIC::Phase_Space_Handler *const pshandler); 
    void ResetSelector(ATOOLS::Selector_Data *const selectordata);

    void Print();

    bool SetPDFFlavours();
    bool Tests();

    void SetScales(const double &scale);

    // inline functions 
    inline void SetSelectorData(ATOOLS::Selector_Data *const data) 
    { p_selectordata=data; }

    inline void PushBack(XS_Base *const xsec) { m_xsecs.push_back(xsec); }

    inline ATOOLS::Selector_Data *const SelectorData() const
    { return p_selectordata; }

    inline XS_Selector *const XSSelector() const { return p_xsselector; }

    XS_Base *Matching(const std::string &name);

  };// end of class XS_Group

}// end of namespace EXTRAXS

#endif


