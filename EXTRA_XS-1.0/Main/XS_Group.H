#ifndef XS_Group_H
#define XS_Group_H

#include "XS_Base.H"

#include <vector>

namespace EXTRAXS {

  class XS_Selector; 

  class XS_Group: public XS_Base {
  protected:

    std::vector<XS_Base *> m_xsecs;

    bool m_atoms, m_channels;

    XS_Selector *p_xsselector;

    void PrepareTerminate();

  public :


    // constructors
    XS_Group(); 
    XS_Group(const size_t nin,const size_t nout,const std::string &name);
    XS_Group(const size_t nin,const size_t nout,const ATOOLS::Flavour *flavours);
    XS_Group(const size_t nin,const size_t nout,const ATOOLS::Flavour *flavours,
	     const int scalescheme,const int kfactorscheme,const double scalefactor,
	     BEAM::Beam_Spectra_Handler *const beamhandler,PDF::ISR_Handler *const isrhandler,
	     ATOOLS::Selector_Data *const selectordata);

    // destructor
    ~XS_Group();

    void Add(XS_Base *const xsec);
    bool Delete(XS_Base *const xsec);
    void Clear();

    virtual void SetAtoms(const bool atoms);
    virtual void CreateFSRChannels();

    size_t Size();

    XS_Base *const operator[](const size_t i) const;

    double Differential2();

    void SetISR(PDF::ISR_Handler *const isrhandler);

    void SetMax(const double max,const int flag);
    void SetTotalXS();

    double Differential(const double s,const double t,const double u);
    bool   SetColours(const double s,const double t,const double u);
    double operator()(const double s,const double t,const double u);
    void   AddPoint(const double xs);

    bool CalculateTotalXSec(const std::string &resultpath);
    void WriteOutXSecs(std::ofstream &outfile);
    void SelectOne();
    void DeSelect();

    bool   OneEvent();
    ATOOLS::Blob_Data_Base *WeightedEvent();

    // inline functions 
    inline void PushBack(XS_Base *const xsec) { m_xsecs.push_back(xsec); }

    inline XS_Selector *const XSSelector() { return p_xsselector; }

  };// end of class XS_Group

}// end of namespace EXTRAXS

#endif


