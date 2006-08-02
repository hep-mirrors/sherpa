#ifndef Hadron_Decay_Channel_H
#define Hadron_Decay_Channel_H

#include "Integrable_Base.H"
#include "HD_ME_Base.H"
#include "Current_Base.H"
#include "HD_PS_Base.H"
#include "Decay_Table.H"
#include <string>
#include "Vector.H"
#include "Spin_Density_Matrix.H"
#include <utility>
#include "Data_Reader.H"


#define HBAR_C 197.326968e-15    // in [GeV mm] used for decay position

namespace HADRONS {

  typedef std::pair<double,std::pair<Current_Base*,Current_Base*> >  CurrentsPair;
  typedef std::pair<double,HD_ME_Base*> MEPair;

  class Hadron_Decay_Channel : public PHASIC::Integrable_Base {
  private:
    ATOOLS::Decay_Channel * p_dc;
    std::vector<MEPair>     m_mes;    // MEs with their relative factors
    std::vector<CurrentsPair> m_currents;
    HD_PS_Base            * p_ps;
//     std::string             m_metype;
    std::string             m_path;
    std::string             m_filename;
    short                   m_fulldecay;
    bool                    m_mass_smearing;
    bool                    m_always_integrate;
    bool                    m_createbooklet;
    std::string             m_channelname, m_chnamenumbers;
    std::vector<std::pair<int,int> > * p_indices;
    std::vector<Complex>  * p_ampls;
    double                  m_symmetry;
    GeneralModel            m_startmd;
    
    void   ProcessOptions(std::vector<std::vector<std::string> > helpsvv);
    void   ProcessPhasespace(std::vector<std::vector<std::string> > ps_svv,
                             Data_Reader                          & reader,
                             GeneralModel                   const & model_for_ps );
    void   ProcessME( std::vector<std::vector<std::string> > me_svv,
                      Data_Reader           & reader,
                      GeneralModel          & model_for_ps );
    void   ProcessResult(std::vector<std::vector<std::string> > result_svv);
    GeneralModel Parameters2Model(std::vector<std::vector<std::string> > helpsvv,
                                  GeneralModel * md = NULL);
    bool   WriteOut( std::vector<double> results, bool newfile=false );
    void   AddAmpls( std::vector<Complex>* newampls, std::vector<std::pair<int,int> >* newindices );
    
  public:
    Hadron_Decay_Channel( ATOOLS::Decay_Channel * _dc, std::string _path );
    ~Hadron_Decay_Channel();
    bool   Initialise(std::vector<std::string> & PStype, GeneralModel startmd);
//     bool   InitialisePhaseSpace(std::vector<std::string> & PStype, GeneralModel startmd );
    void   WriteModelOnScreen( GeneralModel _locmd );
    void   SwapInOrder()    {}
    void   RestoreInOrder() {}
    void   DeSelect()       {}
    bool   SelectOne()      { return 1; }
    bool   ReSelect( int )  { return 1; }
    bool   SelectOneFromList() { return 1; }
    double Differential();
    double Differential( ATOOLS::Vec4D * mom, ATOOLS::Spin_Density_Matrix * );
//    double Differential(double&);
    double Differential(const ATOOLS::Vec4D *momenta) { return 0.; }
    double Differential2()                            { return 0.; }

    void   SetFullDecay( short a ) { m_fulldecay = a; }
    short  FullDecay() { return m_fulldecay; }
    void   SetMassSmearing( bool a ) { m_mass_smearing = a; }
    void   SetAlwaysIntegrate( bool a ) { m_always_integrate = a; }
    bool   MassSmearing() { return m_mass_smearing; }
    void   SetCreateBooklet( bool a ) { m_createbooklet = a; }
    bool   CreateBooklet() { return m_createbooklet; }
    void   CreateTrivial( ATOOLS::Spin_Density_Matrix * sigma );
    
    HD_ME_Base            * SelectME(std::string me_string);
    Current_Base          * SelectCurrent(std::string current_string);

//     HD_ME_Base            * GetME() { return p_me; }
//     Current_Base         ** GetCurrents() { return p_currents; }
    std::vector<MEPair>     GetMEs() { return m_mes; }
    std::vector<CurrentsPair> GetCurrents() { return m_currents; }
    HD_PS_Base            * GetPS() { return p_ps; }
    ATOOLS::Vec4D         * Momenta() { return p_momenta; }        // not as const
    ATOOLS::Vec4D           Momentum( int i ) { return p_momenta[i]; }
    ATOOLS::Decay_Channel * DecayChannel() { return p_dc; }

    std::string                ChannelName()     { return m_channelname; }
    std::string                ChannelNameNumbers() { return m_chnamenumbers; }
    double                  MeanLifetime()    { return HBAR_C/p_ps->Result(); } // in [mm]
    double                  InMass()          { return (p_dc->GetDecaying()).Mass(); }

    std::vector<std::pair<int,int> > * GetIndexList() { return p_indices; }
    std::vector<Complex>  * GetAmplitudeTensor()      { return p_ampls; }
    
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
    \fn Hadron_Decay_Channel::WriteModelOnScreen( GeneralModel _locmd )
    \brief Output for all parameters written in DC file

    This methods writes the values of the parameters that are written in the Decay Channel File
    in the screen.
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
    \var Hadron_Decay_Channel::p_dc
    Pointer on decay channel.
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
    \var Hadron_Decay_Channel::m_mass_smearing
      Stores whether MassSmearing was set to 0 or 1 in the decay channel file.
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
