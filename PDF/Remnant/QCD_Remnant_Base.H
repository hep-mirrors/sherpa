#ifndef PDF_Remnant_QCD_Remnant_Base_H
#define PDF_Remnant_QCD_Remnant_Base_H

#include "PDF/Remnant/Remnant_Base.H"
#include "PDF/Remnant/Color_Dipole.H"
#include "ATOOLS/Phys/Parton_Finder.H"
#include "PDF/Main/ISR_Handler.H"

namespace PDF {

  class QCD_Remnant_Base: public Remnant_Base {
  public:

    typedef std::vector<Color_Dipole*>   Dipole_Vector;
    typedef std::set<ATOOLS::Particle*>  Particle_Set;
    typedef std::vector<ATOOLS::Flavour> Flavour_Vector;

  protected:

    Flavour_Vector m_constit;

    Color_Dipole  *p_start;
    Dipole_Vector  m_connected;

    PDF::PDF_Base *p_pdfbase;

    double       m_deltax;
    unsigned int m_xscheme, m_maxtrials, m_set;
    size_t       m_initial, m_constrained;

    double *p_string;

    Color_Dipole *FindClosest(const Color_Dipole *info,
			      const qri::type type);
    Color_Dipole *FindRandom(const Color_Dipole *info,
			     const qri::type type);
    Color_Dipole *Find(const Color_Dipole *info,
		       const qri::type type);

    bool Connect(const bool sorted); 

    void AssignRemnants();
    void FillRemnants();
 
  public:
    
    // constructor
    QCD_Remnant_Base(PDF::ISR_Handler *isrhandler,
		     const unsigned int _m_beam,const rtp::code type);

    // destructor
    ~QCD_Remnant_Base();
    
    // member functions
    void Clear();

    bool AdjustColors();

    inline void SetStringDrawing(const double opt,const size_t i)
    { p_string[i]=opt; }

  };// end of class QCD_Remnant_Base

}// end of namespace PDF

#endif
