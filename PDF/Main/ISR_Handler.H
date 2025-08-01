#ifndef PDF_Main_ISR_Handler_H
#define PDF_Main_ISR_Handler_H

#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Org/Info_Key.H"
#include "PDF/Main/ISR_Base.H"

namespace ATOOLS   { class Blob_Data_Base; }
namespace BEAM     { class Beam_Base;      }
namespace REMNANTS { class Remnant_Base;   }

namespace PDF {
  struct isr {
    enum id {
      none            =  0,
      hard_process    =  1,
      hard_subprocess =  2,
      bunch_rescatter =  3,
      unknown         = 99
    };
  };

  struct isrmode {
    enum code {
      none           = 0,
      hadron_hadron  = 1,
      lepton_hadron  = 2,
      hadron_lepton  = 3,
      lepton_lepton  = 4,
      unknown        = 99
    };
  };
  std::ostream& operator<<(std::ostream&,isrmode::code mode);


  class ISR_Handler;
  typedef std::map<isr::id,PDF::ISR_Handler*> ISR_Handler_Map;
  typedef std::vector<double> Double_Vector;

  class ISR_Handler {

  protected:
    std::array<ISR_Base*, 2> p_isrbase;

    isr::id       m_id;
    std::string   m_name;
    isrmode::code m_type;

    int m_mode, m_rmode, m_swap;

    std::array<double, 2> m_mass2, m_exponent, m_x, m_mu2;
    ATOOLS::Info_Key m_sprimekey, m_ykey, m_xkey;
    double m_splimits[3],m_ylimits[2];
    double m_fixed_smin, m_fixed_smax;

    ATOOLS::Poincare m_cmsboost;

    std::vector<double> m_info_lab, m_info_cms;

    std::array<BEAM::Beam_Base*, 2>                        p_beam;
    std::array<std::shared_ptr<REMNANTS::Remnant_Base>, 2> p_remnants;

    double m_xf1, m_xf2;

    bool m_freezePDFforLowQ;

    void FixType();

  public:
    ISR_Handler(std::array<ISR_Base *, 2> isrbase,const isr::id & id=isr::hard_process);
    ~ISR_Handler();

    void Init();

    bool CheckConsistency(ATOOLS::Flavour *bunches,ATOOLS::Flavour *partons);
    bool CheckConsistency(ATOOLS::Flavour *partons);
    void SetPartonMasses(const ATOOLS::Flavour_Vector &fl);
    void SetMasses(const ATOOLS::Flavour_Vector &fl);
    bool CheckMasses();

    bool MakeISR(const double &sp,const double &y,
		 ATOOLS::Vec4D_Vector& p,const ATOOLS::Flavour_Vector &flavs);
    bool GenerateSwap(const ATOOLS::Flavour &f1,const ATOOLS::Flavour &f2);
    bool AllowSwap(const ATOOLS::Flavour &f1,const ATOOLS::Flavour &f2) const;

    double PDFWeight(int mode,
                     ATOOLS::Vec4D p1,ATOOLS::Vec4D p2,double Q12,double Q22,
                     ATOOLS::Flavour fl1,ATOOLS::Flavour fl2,int warn=1);
    static double Flux(const ATOOLS::Vec4D& p1);
    static double Flux(const ATOOLS::Vec4D& p1, const ATOOLS::Vec4D& p2);
    double CalcX(ATOOLS::Vec4D p);

    bool BoostInCMS(ATOOLS::Vec4D *p,size_t n);
    bool BoostInLab(ATOOLS::Vec4D *p,size_t n);
    ATOOLS::Poincare* GetCMSBoost() { return &m_cmsboost; }

    void AssignKeys(ATOOLS::Integration_Info *info);
    void SetLimits(double beamy = 0.);
    void Reset();
    void Reset(size_t i) const;

    std::array<std::shared_ptr<REMNANTS::Remnant_Base>, 2> GetRemnants() const
    {
      return p_remnants;
    }
    std::shared_ptr<REMNANTS::Remnant_Base> GetRemnant(const size_t& beam) const
    {
      return beam < 2 ? p_remnants[beam] : nullptr;
    }
    ATOOLS::Blob_Data_Base* Info(int frame) const;

    inline void SetSprimeMin(const double spmin) { m_splimits[0]=ATOOLS::Max(m_fixed_smin,spmin); }
    inline void SetSprimeMax(const double spmax) { m_splimits[1]=ATOOLS::Min(m_fixed_smax,spmax); }
    inline void SetFixedSprimeMin(const double spmin) { m_fixed_smin = m_splimits[0] = spmin; }
    inline void SetFixedSprimeMax(const double spmax) { m_fixed_smax = m_splimits[1] = spmax; }
    inline void SetPole(const double pole)      { m_splimits[2] = pole; }
    inline void SetYMin(const double ymin)      { m_ylimits[0]  = ymin; }
    inline void SetYMax(const double ymax)      { m_ylimits[1]  = ymax; }
    inline void SetName(const std::string& name) { m_name = name;        }
    inline std::string Name() const             { return m_name; }
    inline isrtype::code Type(const int & i=-1) const { return p_isrbase[i]->Type(); }
    inline isr::id Id() const                   { return m_id; }
    inline const isrmode::code & Mode() const   { return m_type; }
    inline int    On() const                    { return m_mode;    }
    inline double * SprimeRange()               { return m_splimits;    }
    inline double * YRange()                    { return m_ylimits;     }
    inline double  Exponent(const int i) const  { return m_exponent[i]; }
    inline double SprimeMin() const             { return m_splimits[0]; }
    inline double SprimeMax() const             { return m_splimits[1]; }
    inline double Pole() const                  { return m_splimits[2]; }
    inline double YMin() const                  { return m_ylimits[0];  }
    inline double YMax() const                  { return m_ylimits[1];  }
    inline double Upper1() const                { return p_isrbase[0]->XMax(); }
    inline double Upper2() const                { return p_isrbase[1]->XMax(); }
    inline double XMax(const size_t beam) const { return p_isrbase[beam]->XMax(); }
    inline double XMin(const size_t beam) const { return p_isrbase[beam]->XMin(); }
    inline double X1() const                    { return m_x[0]; }
    inline double X2() const                    { return m_x[1]; }
    inline double MuF2(int beam) const          { return m_mu2[beam]; }
    inline void   SetMuF2(double mu2, int beam) { m_mu2[beam]=mu2; }
    inline double XF1() const                   { return m_xf1; }
    inline double XF2() const                   { return m_xf2; }
    inline void   SetXF1(double xf)             { m_xf1 = xf; }
    inline void   SetXF2(double xf)             { m_xf2 = xf; }
    inline int Swap() const                     { return m_swap; }

    inline void SetPDF(PDF_Base *pdf) {
      SetPDF(pdf, 0); SetPDF(pdf, 1);
    }
    inline void SetPDF(PDF_Base *pdf, const size_t beam) {
      if (beam<2) p_isrbase[beam]->SetPDF(pdf);
    }
    inline PDF_Base * PDF(const size_t beam) {
      return (beam<2?p_isrbase[beam]->PDF():nullptr);
    }
    inline void SetPDFMember() const {
      for (auto isr : p_isrbase)
	if (isr->On()) isr->PDF()->SetPDFMember();
    }

    inline ATOOLS::Flavour  Flav(const size_t beam) {
      return p_isrbase[beam]->Flavour();
    }
    inline void SetRemnants(
            std::array<std::shared_ptr<REMNANTS::Remnant_Base>, 2> remnants)
    {
      p_remnants = remnants;
    }
    inline void SetRescaleFactor(const double & rescale,const size_t beam) {
      p_isrbase[beam]->SetRescaleFactor(rescale);
    }
    inline void ResetRescaleFactor(const size_t beam) { SetRescaleFactor(1.,beam); }
    inline void SetBeam(BEAM::Beam_Base *const beambase,const size_t beam) {
      if (beam<2) p_beam[beam]=beambase;
    }
    inline BEAM::Beam_Base * GetBeam(const size_t beam) const {
      return beam<2?p_beam[beam]:nullptr;
    }
    inline void SetRunMode(const int &rmode) { m_rmode=rmode; }

    void Output();
  };// end of class ISR_Handler

  /*!
    \namespace PDF
    The namespace PDF houses all classes that are employed to generate
    parton spectra. In the framework of both the SHERPA package and of the
    program AMEGIC the following nomenclature is assumed :
    - There are incoming beams at a certain energy, the nominal energy of the
      beam in the collider, which then result in bunches of interacting
      particles which have an energy distribution, and, maybe, a \f$k_\perp\f$
      distribution of transverse momenta w.r.t.the beam axis.
    - The bunch particles then can have a substructure, i.e. they might consist of
      partons, which then interact in a hard subprocess.

    As an illustrative example, consider the case of DIS of an electron on a photon.
    The incoming electron beam emits bunches of photons that might or might not
    resolved into partons during the interaction with the proton. In the PDF namespace,
    the parton distribution funcitons of both the photon and the proton are handled.
  */
  /*!
    \class ISR_Handler
    \brief Manager of all Initial State Radiation that can be identifeid as parton
           distributions.
    This class manages all initial state radiation (ISR) according to the parton
    distribution functions (PDFs) that are handed over. The ISR_Handler is initialized from
    the SHERPA package or from Amegic. Before coming into full effect during integration
    or event generation, it initalises suitable ISR treatment through ISR_Bases that will
    contain the PDFs for each of the bunches.
  */
  /*!
    \var ISR_Base ** ISR_Handler::p_isrbase
    Pointers to the two ISR bases, one for each bunch.

    \sa ISR_Base
  */
  /*!
     \var int ISR_Handler::m_mode
     The m_mode flag indicates what kind of ISR treatment is to be considered:
     - 0 no ISR for both bunchs
     - 1 only bunch 1 experiences ISR
     - 2 only bunch 2 experiences ISR
     - 3 both bunches experience ISR.
  */
  /*!
    \var ATOOLS::Poincare ISR_Handler::m_cmsboost
    A boost from the c.m. system of the incoming bunches to the c.m. system of the
    outgoing partons, which form the initial state of the hard interaction.
  */
  /*!
    \var double ISR_Handler::m_exponent[2]
    Characteristic exponents used for the integration.
  */
  /*!
    \var double ISR_Handler::m_splimits[3]
    \f$s'\f$-limits and characteristics:
    m_splimits[0,1] = \f$s'_{\rm min, max}\f$
    and
    m_splimits[2] = \f$s_{\rm beams}\f$.
  */
  /*!
    \var ISR_Handler::m_ylimits[2]
    The rapidity region covered. It is per default set to the range \f$y \in [-10,10]\f$. In fact
    this range should be calculated from the range of the BeamBases.
    \todo Rapidity range from BeamBases.
  */
  /*!
    \var double ISR_Handler::m_mass2
    Squares of the masses of the incoming particles.
  */
  /*!
    \var double ISR_Handler::m_x
    The energy fractions each outgoing bunch has w.r.t. the corresponding
    incoming beams.
  */
  /*!
    var std::string ISR_Handler::m_name
    Name of the ISR_Handler.
  */
  /*!
    var std::string ISR_Handler::m_type
    Type of the ISR_Handler, it consists of the types of the BeamBases.
  */
  /*!
    \fn bool ISR_Handler::CheckConsistency(ATOOLS::Flavour *,ATOOLS::Flavour *)
    This checks whether the two sets of flavours match the flavours of the incoming
    and outgonig particles of the two p_isrbases. If this is the case, true is returned.
    This method is largely similar to the corresponding one in the BEAM::Beam_Spectra_Handler.
  */
  /*!
    \fn bool ISR_Handler::CheckConsistency(ATOOLS::Flavour *)
    This checks whether the flavours are allowed to be used as outgonig flavours in the
    two p_isrbases. If this is the case, true is returned. This method is largely similar to
    the corresponding one in the BEAM::Beam_Spectra_Handler.
  */
  /*!
    \fn void ISR_Handler::SetPartonMasses(ATOOLS::Flavour * _fl)
    This sets the masses squared and the vectors such that they fit to the masses of the
    flavours.
  */
  /*!
    \fn bool ISR_Handler::MakeISR(ATOOLS::Vec4D *,double,double);
    Depending on the \f$s'\f$-value handed over as arguments, two matching vectors for the
    outgoing partons in their c.m. frame (out) are constructed. Then the energy fractions in the
    c.m. system (in) of the incoming bunches are determined with help of the other argument, the
    rapidity \f$y\f$ according to
    \f[
    \hat E^{(in)}_{1,2} = \exp\left(\pm y\right)
    \f]
    and the boost linking the two frames, CMBoost is initialized. This boost is then used
    to bring the c.m. vectors into the correct frame, i.e. the c.m. frame
    of the bunches, i.e.
    \f[
    p^{(out)}_{1,2} \Longrightarrow p^{(in)}_{1,2}\,.
    \f]
    This method is largely similar to the corresponding one in the
    BEAM::Beam_Spectra_Handler.
  */
}

#endif
