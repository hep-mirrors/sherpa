#ifndef Ladder_H
#define Ladder_H

#include "BFKL_Sudakov.H"
#include "ISR_Handler.H"
#include "XS_Group.H"

namespace PDF { class Doubly_Unintegrated_PDF; }

namespace EXTRAXS {

  class Ladder: public XS_Group {
  public:

    typedef std::vector<ATOOLS::Vec4D>   Vector_Vector;
    typedef std::vector<ATOOLS::Flavour> Flavour_Vector;

  private:

    BFKL_Sudakov *p_sudakov;

    PDF::Doubly_Unintegrated_PDF *p_pdfs[2];

    Vector_Vector  m_moms;
    Flavour_Vector m_flavs;

    ATOOLS::Vec4D m_pa, m_pb, m_k1, m_kn, m_q, m_oldq, m_p1, m_p2;

    double m_weight, m_ya, m_yb, m_kt2min, m_q2, m_y1, m_yn;
    double m_a1, m_an, m_z1, m_zn, m_kt21, m_kt2n, m_mu21, m_mu2n;
    double m_y, m_kt2, m_phi, m_lasty;

    size_t m_ncols, m_sudmode, m_nfixed, m_multimode, m_splitmode;

    bool GeneratePDFJet();
    bool ConstructRung();
    bool ConstructIncoming();
    bool GenerateLadder();

    double Flux() const;

    bool TestEmission();
    bool SetScales();
    void InitSudakov();

    bool DiceOneEmission();

  public:

    // constructor
    Ladder(const size_t nin,const size_t nout,
	   const ATOOLS::Flavour *flavours,
	   const PHASIC::scl::scheme scalescheme,const int kfactorscheme,
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

    void AddEvent(const double xs,const double validxs,
		  const int ncounts);

  };// end of class Ladder

}// end of namespace EXTRAXS

#endif
