#ifndef Ladder_H
#define Ladder_H

#include "BFKL_Sudakov.H"
#include "ISR_Handler.H"
#include "XS_Group.H"

namespace EXTRAXS {

  class Ladder: public XS_Group {
  public:

    typedef std::vector<ATOOLS::Vec4D>   Vector_Vector;
    typedef std::vector<ATOOLS::Flavour> Flavour_Vector;

  private:

    BFKL_Sudakov *p_sudakov;

    PDF::PDF_Base *p_pdfs[2];

    Vector_Vector  m_moms;
    Flavour_Vector m_flavs;

    ATOOLS::Vec4D m_pa, m_pb, m_k1, m_k2, m_q;

    double m_weight, m_ya, m_yb, m_kt2min, m_q2;
    double m_y1, m_kt21, m_phi1, m_z1;
    double m_a1, m_a2, m_z2, m_kt22, m_mu21, m_mu22;

    bool GeneratePSJet();
    bool ConstructRung();
    bool GenerateLadder();

    double Jacobian();
    bool   SetScales();

  public:

    // constructor
    Ladder(const size_t nin,const size_t nout,
	   const ATOOLS::Flavour *flavours,
	   const int scalescheme,const int kfactorscheme,
	   BEAM::Beam_Spectra_Handler *const beamhandler,
	   PDF::ISR_Handler *const isrhandler,
	   ATOOLS::Selector_Data *const selectordata);

    // destructor
    ~Ladder();

    // member functions
    bool Initialize();

    bool SetColours(const ATOOLS::Vec4D *momenta);
    bool SetColours(const double s,const double t,const double u);

    double Differential(const ATOOLS::Vec4D *momenta);
    double Differential2();
    
    void CreateFSRChannels();
    void CreateISRChannels();

    bool CalculateTotalXSec(const std::string &rpath,
			    const bool create=false);

    bool SelectOne();
    void DeSelect();

  };// end of class Ladder

}// end of namespace EXTRAXS

#endif
