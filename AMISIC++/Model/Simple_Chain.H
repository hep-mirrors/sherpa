#ifndef AMISIC_Model_Simple_Chain_H
#define AMISIC_Model_Simple_Chain_H
/*!
  \file Simple_Chain.H
  \brief Declares the class Simple_Chain
*/

#include "AMISIC++/Main/MI_Base.H"
#include "AMISIC++/Model/Grid_Creator.H"
#include "AMISIC++/Model/Profile_Function_Base.H"
#include "AMISIC++/Tools/Fixed_Variable_Channel.H"
#include "AMISIC++/Tools/Semihard_QCD.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "MODEL/Main/Model_Base.H"

namespace ATOOLS { class Data_Reader;  }
namespace PDF    { class Remnant_Base; }

namespace AMISIC {

  class Simple_Chain: public MI_Base {
  public:

    typedef Amisic_Histogram<double> Amisic_Histogram_Type;
    typedef std::map<std::string,Amisic_Histogram_Type*>  Amisic_Histogram_Map;

    typedef PHASIC::Fixed_Variable_Channel FSR_Channel;

    typedef std::map<std::string,PHASIC::Process_Base*> Process_Map;
    typedef std::multimap<double,std::string,
			  std::greater<double> >    Sort_Map;

  private:

    Amisic_Histogram_Map   m_differentials;
    Amisic_Histogram_Type *p_differential, *p_total;
    Grid_Creator          *p_gridcreator;

    Process_Map m_processmap;

    double m_norm, m_enhance, m_maxreduction, m_sigma_nd_fac;

    std::string m_xsextension, m_pathextra, m_selectorfile;
    std::string m_resdir, m_ressuffix;

    std::vector<Semihard_QCD*> p_processes;
    std::vector<FSR_Channel*>  p_fsrinterface;

    MODEL::Model_Base          *p_model;
    BEAM::Beam_Spectra_Handler *p_beam;
    PDF::ISR_Handler           *p_isr;

    Profile_Function_Base *p_profile;

    std::string m_selected;

    unsigned int m_maxtrials, m_check, m_vegas, m_pi;
    double       m_sigmahard, m_error, m_ecms;
    bool         m_external, m_regulate, m_heavy_flavour;

    std::string         m_regulator; 
    std::vector<double> m_regulation;

    std::string m_xsec_output;
    
    PDF::Remnant_Base *p_remnants[2];

    ATOOLS::Info_Key m_spkey, m_ykey, m_xkey, m_isrspkey, m_isrykey;

    ATOOLS::Data_Reader *p_read;

    void Init();
    void CleanUp();

    bool GeneratePathName();
    bool ReadInData();

    bool CreateGrid();
    void InitializeProcessList(const ATOOLS::Flavour& in1,
                               const ATOOLS::Flavour& in2,
                               const ATOOLS::Flavour& out1,
                               const ATOOLS::Flavour& out2);
    void CalculateSigmaND();
    /*!
      \fn void CalculateSigmaND()
      \brief Calculates the non-diffractive inelastic cross section
      according to the model by G. A. Schuler and T. Sjostrand 
      (see Phys. Rev. D <b>49</b> (1994) 2257).
    */
    int CalculateTotal();
    /*!
      \fn bool CalculateTotal()
      \brief Determines the total inelastic cross section 
      according to the parton model from the initialized histograms.
    */
    bool SetUpInterface();
    /*!
      \fn bool SetUpInterface()
      \brief Sets up the phasespace integrators for event generation.
    */

    void SetISRRange();
    /*!
      \fn void SetISRRange()
      \brief Sets the isr range for event generation, which is
      \f$4p_\perp^2\,\le\,\hat{s}\,\le\,4E_{1, max}E_{2, max}\f$,
      \f$\log(p_\perp/E_{2, max})\,\le\,y\,\le\,\log(E_{1, max}/p_\perp)\f$
    */
    void ResetISRRange();
    /*!
      \fn void ResetISRRange()
      \brief Resets the ISR range to the ISR handler default.
    */

    bool CreateMomenta();
    /*!
      \fn bool CreateMomenta()
      \brief Employs the selected process to determine the kinematics
      of the interaction. The outgoing transverse momentum is kept fixed. 
      It is generated in advance by GenerateOrderingParameter.      
    */

    bool ReadInStatus(const std::string &path);
    /*!
      \fn void PrepareTerminate()
      \brief Restores grid calculations in case of debugging.
    */

    void PrepareTerminate();
    /*!
      \fn void PrepareTerminate()
      \brief Stores grid calculations in case of program abort.
    */

  public:

    // constructor
    Simple_Chain();
    /*!
      \fn Simple_Chain()
      \brief The default constructor.
      Creates an ISR handler, a beam spectra handler
      and a model handler to be used within Simple_Chain.
    */
    Simple_Chain(MODEL::Model_Base *_p_model,
		 BEAM::Beam_Spectra_Handler *_p_beam,
		 PDF::ISR_Handler *_p_isr);
    /*!
      \fn Simple_Chain(MODEL::Model_Base *_p_model,
      BEAM::Beam_Spectra_Handler *_p_beam,
      PDF::ISR_Handler *_p_isr)
      \brief The standard constructor.
      Initializes Simple_Chain with external ISR, beam and 
      model handlers.
    */

    // destructor
    ~Simple_Chain();

    // member functions
    bool Initialize();
    void Update(const MI_Base *mibase);
    void Reset();

    bool GenerateEnhanceFactor();
    /*!
      \fn bool GenerateEnhanceFactor()
      \brief Generates an enhance factor according to the 
      hadron profile in impact parameter space. Employs
      Profile_Function_Base.
    */
    bool GenerateOrderingParameter();
    /*!
      \fn bool GenerateOrderingParameter()
      \brief Determines the maximum outgoing transverse momentum 
      of the next hard interaction.
    */

    bool VetoProcess(ATOOLS::Blob *blob);
    /*!
      \fn bool VetoProcess(ATOOLS::Blob *blob)
      \brief Vetoes the given hard process according to 
      \f$p(p_\perp)\,\propto\;
      {\rm d}\sigma\,/\,{\rm d}p_\perp\f$.
     */
    bool GenerateProcess();
    /*!
      \fn bool GenerateProcess()
      \brief Selects a process type for the interaction according to 
      the differential probabilities. Sets and resets the ISR ranges.
    */

    // inline functions
    inline double Norm() { return m_norm; }

    inline Profile_Function_Base *ProfileFunction() { return p_profile; }

    inline void SetNorm(double norm) { m_norm=norm; }

  }; // end of class Simple_Chain
  /*!
    \class Simple_Chain
    \brief Simulates multiple interactions according to the model
    by Sjostrand and Zijl (Phys. Rev. D <b>36</b> (1987) 2019)

    This class handles the multiple interaction generation according 
    to the model by T. Sjostrand and M. van Zijl (see Phys. Rev. D 
    <b>36</b> (1987) 2019). 
    To determine the differential interaction probability 
    \f[
    p(p_{\perp})\,=\;\frac{1}{\sigma_{\rm ND}}\,
    \frac{{\rm d}\sigma}{{\rm d}p_{\perp}}
    \f] 
    the corresponding processes of the internal process library 
    are employed. The class Grid_Creator is used to fill the 
    corresponding histograms which are then added and integrated 
    numerically.

    During the event generation Simple_Chain is steered by the 
    class Amisic. From it it receives information about the last hard 
    scale \f$p_{\perp i-1}\f$ and the maximally available beam energies.
    The next interaction is then found by solving for \f$p_{\perp i}\f$
    the relation
    \f[
    p_{\perp\,i}\,=\;I^{-1}\left[\;I(p_{\perp\,i-1})-\log R\;\right]\;,
    \f]
    where the integrated probability \f$I(p_\perp)\f$ is defined by
    \f[
    I(p_\perp)\,=\;\int_{p_\perp^2}^{\,s/4}{\rm d} p_\perp'^2\,
    p(p_\perp')\;.
    \f]
    If the interacting hadrons are considered extended objects in
    impact parameter space then this relation reads
    \f[
    p_{\perp\,i}\,=\;I^{-1}\left[\;I(p_{\perp\,i-1})-
    \frac{\log R}{e(b)}\;\right]\;,
    \f]
    where \f$e(b)\f$ is the enhance factor according to
    \f[
    p(p_\perp,b)\,=\;e(b)\,\frac{1}{\sigma_{\rm ND}}\,
    \frac{{\rm d} \sigma}{{\rm d} p_\perp^2}\;.
    \f]
    It is found through a suitable veto algorithm.<br>
    The non-diffractive hadronic cross section is determined 
    according to a model by G. A. Schuler and T. Sjostrand 
    (see Phys. Rev. D <b>49</b> (1994) 2257).
  */

} // end of namespace AMISIC

#endif
