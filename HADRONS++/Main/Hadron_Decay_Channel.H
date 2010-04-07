#ifndef HADRONS_Main_Hadron_Decay_Channel_H
#define HADRONS_Main_Hadron_Decay_Channel_H

#include <string>
#include <utility>
#include "ATOOLS/Org/Data_Reader.H"
#include "HADRONS++/Main/Tools.H"
#include "ATOOLS/Phys/Decay_Table.H"
#include "ATOOLS/Math/Vector.H"

namespace ATOOLS {
  class Particle;
}
namespace METOOLS {
  class Spin_Amplitudes;
  class Amplitude_Tensor;
}

namespace HADRONS {
  class Current_Base;
  class HD_ME_Base;
  class HD_PS_Base;

  typedef std::pair<Complex,std::pair<Current_Base*,Current_Base*> >  CurrentsPair;
  typedef std::pair<Complex,HD_ME_Base*> MEPair;

  class Hadron_Decay_Channel : public ATOOLS::Decay_Channel {
  private:
    ATOOLS::Flavour       * p_flavours;
    ATOOLS::Flavour       * p_physicalflavours;
    ATOOLS::Vec4D         * p_momenta;
    std::vector<MEPair>     m_mes;    // MEs with their relative factors
    std::vector<CurrentsPair> m_currents;
    HD_PS_Base            * p_ps;
    std::string             m_path;
    std::string             m_filename;
    std::string             m_origin;
    bool                    m_always_integrate;
    double                  m_symmetry;
    // CP asymmetry coefficients as used in PDG B0 listings
    double                  m_cp_asymmetry_C;
    double                  m_cp_asymmetry_S;
    Complex                 m_cp_asymmetry_lambda;
    double                  m_original_width;
    GeneralModel            m_startmd;
    METOOLS::Spin_Amplitudes * p_amplitudes;
    
    void   ProcessOptions(std::vector<std::vector<std::string> > helpsvv);
    void   ProcessPhasespace(std::vector<std::vector<std::string> > ps_svv,
                             ATOOLS::Data_Reader                  & reader,
                             GeneralModel                   const & model_for_ps );
    void   ProcessME( std::vector<std::vector<std::string> > me_svv,
                      ATOOLS::Data_Reader                  & reader,
                      GeneralModel                         & model_for_ps );
    void   ProcessResult(std::vector<std::vector<std::string> > result_svv);
    GeneralModel Parameters2Model(std::vector<std::vector<std::string> > helpsvv,
                                  GeneralModel * md = NULL);
    bool   WriteOut( std::vector<double> results, bool newfile=false );
    std::vector<double> CalculateResults();
    
  public:
    Hadron_Decay_Channel( ATOOLS::Flavour flin, std::string _path );
    ~Hadron_Decay_Channel();
    void   Initialise(GeneralModel startmd);
    double Differential();
    double Differential( ATOOLS::Vec4D * mom, bool anti );

    void   SetAlwaysIntegrate( bool a ) { m_always_integrate = a; }
    double CPAsymmetryC() { return m_cp_asymmetry_C; }
    double CPAsymmetryS() { return m_cp_asymmetry_S; }
    Complex CPAsymmetryLambda() { return m_cp_asymmetry_lambda; }
    double OriginalWidth() { return m_original_width; }
    void   SetOriginalWidth(double origwidth) { m_original_width = origwidth; }
    std::string FileName() { return m_filename; }
    void   SetFileName(std::string filename="");
    inline const std::string& Origin() { return m_origin; }
    inline void SetOrigin(const std::string& origin) { m_origin=origin; }
    
    HD_ME_Base            * SelectME(std::string me_string);
    Current_Base          * SelectCurrent(std::string current_string);

    std::vector<MEPair>     GetMEs() { return m_mes; }
    std::vector<CurrentsPair> GetCurrents() { return m_currents; }
    HD_PS_Base            * GetPS() { return p_ps; }
    ATOOLS::Flavour       * Flavours() { return p_flavours; }
    int                     GetN() { return 1+NOut(); }

    void   CalculateAmplitudes(ATOOLS::Vec4D* moms, METOOLS::Spin_Amplitudes* amps, bool anti);
    void   CalculateAmplitudes(ATOOLS::Vec4D* moms, METOOLS::Amplitude_Tensor* amps, bool anti);
    bool   SetColorFlow(std::vector<ATOOLS::Particle*> outparts,int n_q, int n_g);
    
  }; // end of class
   
  /*!
    \file Hadron_Decay_Channel.H
    \brief Declares the class HADRONS::Hadron_Decay_Channel

    This code can be found in the directory \c Main.
  */    

  /*!
    \fn Hadron_Decay_Channel::Hadron_Decay_Channel( ATOOLS::Decay_Channel * _dc, std::string _path )
    \brief Constructor

    This is the constructor for a hadron decay channel. It uses the information 
    from the pointer <var> _dc </var>
    for initialisation. It calls the HD_ME_Selector::GetME method in order
    to set the pointer Hadron_Decay_Channel::p_me on the right ME class. 
    */
  /*!
    \fn Hadron_Decay_Channel::InitialisePhaseSpace(std::vector<std::string> & PStype, GeneralModel startmd) 
    \brief Reads the decay table

    This method actually reads the decay table, e.g. \c TauDecays.dat.
    The information as a string matrix of the form
    <TABLE>
    <TR> <TD>1st column</TD><TD>2nd column</TD><TD>3rd column</TD> </TR>
    <TR> <TD>decay products</TD><TD>branching ratio</TD><TD>name of DC file</TD> </TR>
    </TABLE>
    If a name of the DC file is \b not given, this method creates a new DC file with a generic name.
    This DC file can be seen as template with the most basic entries (the intregrator is Rambo and 
    no parameters are given).

    It calls the HD_PS_Base constructor for reading the Decay Channel File.
    At the end the parameters that were written in the DC file are passed to the corresponding
    ME class (that's where the pointer <var>p_me</var> points to) by calling 
    the virtual method HD_ME_Base::SetModelParameters.
    */
  /*!
    \fn Hadron_Decay_Channel::CreateTrivial( ATOOLS::Spin_Density_Matrix * sigma )
    \brief Create trivial amplitude tensor

    Creates a trivial amplitude tensor that is of correct size (in case the ME is unable to produce
    such a tensor. It gets stored in Hadron_Decay_Channel::p_ampls. The index list
    is stored in Hadron_Decay_Channel::p_indices.
    */
    

  /*!
    \class Hadron_Decay_Channel
    \brief Tool to manage a certain hadronic decay channel

    This class is a subclass of PHASIC::Integrable_Base.
  */    
  
  /*!
    \fn Hadron_Decay_Channel::Hadron_Decay_Channel( ATOOLS::Decay_Channel * _dc )
    \brief Constructor for a hadronic decay channel

    Sets initial values such as result path and file name, decay name, flavours of involved particles
    and uses HD_ME_Selector::GetME to map the corresponding decay's ME.
  */
  /*!
    \fn Hadron_Decay_Channel::InitialisePhaseSpace(std::vector<std::string> & PStype )
    \brief Method to have the phase space initialised.

    Uses HD_PS_Base::HD_PS_Base in order to set up the phase space correctly and calls
    HD_PS_Base::Initialise for initialisation if the PS generator has to be initialised,
    i.e. no settings from a last run is written in file. If there are, these integration channel settings
    are used for PS integration.
  */    
  /*!
    \fn Hadron_Decay_Channel::Differential()
    \brief Returns differential rate \f$d\Gamma \times dlips\f$

    This method makes the phase space generator choose a phase space point (the
    mother particle decays from rest) and its weight factor.
    Hereby, the pointer on the phase space generator Hadron_Decay_Channel::p_ps plays a key role.
    Afterwards, it calls the ME class (via the pointer Hadron_Decay_Channel::p_me) 
    which has to fill the amplitude tensor Hadron_Decay_Channel::p_ampls
    and indix list Hadron_Decay_Channel::p_indices. Both are mere C++ vectors.
    The amplitude tensor contains the value of the matrix element (\b not squared) for
    each helicity combination. The index tensor is necessary for the correct book keeping
    (see HD_ME_Base class for more details). However, since this Differential method
    does not care about spin correlations it merely sums up all absolute squared entries of 
    the amplitude tensor Hadron_Decay_Channel::p_ampls in order to obtain a value for
    the distribution.
  */
  /*!
    \fn double Hadron_Decay_Channel::Differential( Vec4D * mom, Spin_Density_Matrix * sigma )
    \brief Returns differential rate \f$d\Gamma \times dlips\f$ with spin correlations

    If a a pointer on a spin density matrix <var>sigma</var> is given, this methods
    returns the "matrix element squared" value with spin correlations.
    It makes the phase space generator choose a phase space point (the
    mother particle decays from rest) and its weight factor.
    Afterwards, all momenta are boosted into the lab system for correct spin correlations.
    It then calls the ME class (via the pointer Hadron_Decay_Channel::p_me) 
    which has to fill the amplitude tensor Hadron_Decay_Channel::p_ampls
    and indix list Hadron_Decay_Channel::p_indices. Both are mere C++ vectors.
    The amplitude tensor contains the value of the matrix element (\b not squared) for
    each helicity combination. The index tensor is necessary for the correct book keeping
    (see HD_ME_Base class for more details). 
    It sums up the ME squares for each helicity combination weighting with the entries
    of the spin density matrix. For that it uses the Spin_Correlation_Tensor::Contract and
    Spin_Correlation_Tensor::Trace method.
     
    However, if the ME class does not give an index list things are traeted without
    spin correlation, that is treating the spin density matrix as if it were a unit matrix.

    There are also effective ME's that do not have information for each helicity combination.
    Their amplitude tensor only contains one entry. Then a trivial amplitude tenser is created
    in order to ensure correct behaviour of Hadrons::PerformDecay method (no crash), because
    the Hadrons::PerformDecay method expects an amplitude tensor of correct size.
   */
  /*!
    \var Hadron_Decay_Channel::p_me
    Pointer on ME.
  */    
  /*!
    \var Hadron_Decay_Channel::p_ps
    Pointer on PS generator.
  */
  /*!
    \var Hadron_Decay_Channel::m_metype
    Type of ME gained from ATOOLS::Decay_Channel class.
  */ 
  /*!
    \var Hadron_Decay_Channel::m_fulldecay
    Has information about decay: syntax "(binary number) = decimal number"
    - (00) = 0 ... no full decay, anti particle does not decay
    - (01) = 1 ... no full decay, anti particle does decay
    - (10) = 2 ... full decay, anti particle does not decay
    - (11) = 3 ... full decay, anti particle does decay
    .
    "no full decay" means that unstable daughters do not decay.
    "anti particle does not decay" means that the anti-particle is
    seen to be stable.
  */
  /*!
    \var Hadron_Decay_Channel::p_ampls
    Pointer on the amplitude tensor.
    */
  /*!
    \var Hadron_Decay_Channel::p_indices
    Pointer on the index list.
    */

    
    

  
  
} // end of namespace


#endif
