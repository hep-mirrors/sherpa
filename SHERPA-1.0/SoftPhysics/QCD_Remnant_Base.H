#ifndef QCD_Remnant_Base_H
#define QCD_Remnant_Base_H

#include "Remnant_Base.H"
#include "Color_Dipole.H"
#include "Parton_Finder.H"
#include "ISR_Handler.H"

namespace SHERPA {

  class QCD_Remnant_Base: public Remnant_Base {
  public:

    typedef std::vector<Color_Dipole*> Dipole_Vector;

  protected:

    std::vector<ATOOLS::Flavour> m_constit;

    Color_Dipole  *p_start;
    Dipole_Vector  m_connected, m_attached;

    PDF::PDF_Base *p_pdfbase;

    double       m_deltax, m_scale, m_ecms;
    unsigned int m_xscheme, m_maxtrials, m_set;
    size_t       m_initial, m_constrained;

    bool m_dupdf;

    Color_Dipole *FindClosest(const Color_Dipole *info,
			      const qri::type type);

    bool Connect(const bool sorted); 

    void AssignRemnants();
    bool ConnectRemnants();
    void FillRemnants();
 
  public:
    
    // constructor
    QCD_Remnant_Base(PDF::ISR_Handler *isrhandler,
		     const unsigned int _m_beam,
		     const double _m_scale,const rtp::code type);

    // member functions
    void Clear();

  };// end of class QCD_Remnant_Base

}// end of namespace SHERPA

#endif
