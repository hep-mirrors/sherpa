#ifndef Ladder_H
#define Ladder_H

#include "BFKL_Sudakov.H"
#include "ISR_Handler.H"
#include "XS_Group.H"
#include "Jet_Finder.H"

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
    Flavour_Vector m_flavs, m_props;

    ATOOLS::Vec4D   m_pa, m_pb, m_k1, m_kn, m_q, m_oldq, m_p1, m_p2;
    ATOOLS::Flavour m_fl1, m_fln;

    double m_weight, m_ya, m_yb, m_kt2min, m_q2, m_y1, m_yn, m_sab;
    double m_a1, m_an, m_z1, m_zn, m_kt21, m_kt2n, m_mu21, m_mu2n;
    double m_y, m_kt2, m_phi, m_lasty, m_dxs, m_bfkldxs, m_splitweight;
    double m_asecms;

    size_t m_ncols, m_sudmode, m_nfixed, m_multimode, m_splitmode;
    size_t m_memode, m_pnmode, m_nmaxme, m_hfp;

    std::map<std::string,XS_Base*> m_xsmes;

    ATOOLS::Poincare m_cms, m_zaxis;

    ATOOLS::Jet_Finder *p_jf;

    double YMax(const double &kt2,const double &sab,
		const double &s1=0.0,const double &sn=0.0) const;
    double KT2Max(const double &y,const double &sab,
		  const double &s1=0.0,const double &sn=0.0) const;

    bool GeneratePDFJet();
    bool BoostToCMS();
    bool BoostToLab();
    bool ConstructRung(const bool in=false);
    bool ConstructIncoming();
    bool GenerateLadder();

    double Flux() const;

    bool TestEmission();
    bool SetScales();

    bool DiceFixedMulti();
    int  DiceVariableMulti();

    void SetColours();
    void SetColour(const size_t &i,const bool fw,int *const lc);

    std::string MapProcess(const std::vector<ATOOLS::Flavour> &fl);

    double MECorrection();

  public:

    // constructor
    Ladder(const size_t nin,const size_t nout,
	   const ATOOLS::Flavour *flavours,
	   const PHASIC::scl::scheme scalescheme,const int kfactorscheme,
	   BEAM::Beam_Spectra_Handler *const beamhandler,
	   PDF::ISR_Handler *const isrhandler,
	   ATOOLS::Selector_Data *const selectordata,
	   XS_Model_Base *const model);

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
