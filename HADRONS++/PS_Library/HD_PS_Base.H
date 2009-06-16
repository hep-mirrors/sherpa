#ifndef HADRONS_PS_Library_HD_PS_Base_H
#define HADRONS_PS_Library_HD_PS_Base_H

#include "PHASIC++/Channels/Single_Channel.H"
#include "PHASIC++/Channels/Multi_Channel.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/Vector.H"
#include <string>
#include "PHASIC++/Channels/Vegas.H"
#include "HADRONS++/Main/Tools.H"

namespace HADRONS {



  class Hadron_Decay_Channel;

  struct ChannelInformation
  {
    std::string name;                // name of integrator
    int nRes;                        // number of resonances
    std::string res1, res2, res3;    // resonances
    int a,b,c,d;                    // index for outgoing particles
  };

  class HD_Channel_Selector {
    bool DecomposeChannel( std::string name, ChannelInformation & ci );
  public:
    PHASIC::Single_Channel * GetChannel(
        int nin, int nout, const ATOOLS::Flavour * flavs, std::string name, GeneralModel const &);
  }; // end of class
  
  class HD_PS_Base : public PHASIC::Multi_Channel {
  private:
    HD_Channel_Selector  * p_channelselector;
    Hadron_Decay_Channel * p_hdc;
    double                 m_res,m_error,m_max,m_flux;
  public:
    HD_PS_Base( Hadron_Decay_Channel * hdc );
    std::vector<double> CalculateNormalisedWidth();
    virtual ~HD_PS_Base();
    bool           AddChannel(std::string,double, GeneralModel const &);
    bool           IsChannel( std::string );
//     void           Initialise();
    inline double  Result() { return m_res; }
    inline double  Error()  { return m_error; }
    inline double  RelError() { return m_error/m_res*100.; }    // in %
    inline double  Maximum() { return m_max; }
    inline void    SetResult(double res)  { m_res=res; }
    inline void    SetError(double error) { m_error=error; }
    inline void    SetMaximum(double max) { m_max=max; }
  }; // end of class

  /*!
    \file HD_PS_Base.H
    \brief Declares the classes HADRONS::HD_Channel_Selector, HADRONS::ChannelInformation, and HADRONS::HD_PS_Base.

    This file can be found in the directory \c PS_Library.
  */

  /*!
    \class HD_Channel_Selector
    \brief Tool to select an integration channel
  */
  /*!
    \fn HD_Channel_Selector::GetChannel(int nin, int nout, const ATOOLS::Flavour * flavs, std::string name, GeneralModel & md )
    \brief Returns pointer on a PS generator

    #in <var>nin</var> must be 1 and the #out <var>nout</var>greater equal than 2.
    The argument <var>name</var> corresponds to the name of the integrator such as
    - Isotropic (Rambo)
    - Dalitz_<var>prop</var>_<var>ij</var>
    - TwoResonances_<var>prop1</var>_<var>k</var>_<var>prop2</var>_<var>ij</var>
    - ...
    .
    At first the name is decomposed by HD_Channel_Selector::DecomposeChannel and stored as
    ChannelInformation. Depending in the information in ChannelInformation it returns
    a pointer on the corresponding PS integrator class such as
    - Rambo
    - Dalitz
    - TwoResonances
    - ...
    .
  */    
  /*!
    \fn HD_Channel_Selector::DecomposeChannel( string name, ChannelInformation & ci )
    \brief Reads name and extract the information

    This method extracts the required information from the name <var>name</var> and stores
    them in the ChannelInformation class variable <var>ci</var>.
    */

  /*!
    \class ChannelInformation
    \brief Auxilliary class for setting the correct integrator

    This class is a mere bookkeeping class for all information that is extracted from the name
    of the integrator:
    - which are the propagators
    - which particles are external legs
    - ...
    .
    It can cope with the following integrators
    <TABLE>
    <TR> <TD>full name</TD> <TD></TD> <TD>Isotropic</TD>
        <TD>Dalitz_<var>prop</var>_<var>ij</var></TD>
        <TD>TwoResonances_<var>prop1</var>_<var>k</var>_<var>prop2</var>_<var>ij</var></TD> </TR>
    <TR> <TD>diagram</TD> <TD></TD> <TD></TD> <TD>\f$prop \to i+j\f$</TD> 
        <TD>\f$prop_1 \to k+ prop_2;\: prop_2 \to i+j\f$</TD> </TR>
    <TR> <TD>type</TD> <TD>name</TD> <TD>Isotropic</TD> <TD>Dalitz</TD> <TD>TwoResonances</TD> </TR>
    <TR> <TD>number of resonances</TD> <TD>nRes</TD> <TD>0</TD> <TD>1</TD> <TD>2</TD> </TR>
    <TR> <TD>resonance 1</TD> <TD>res1</TD> <TD>-</TD> <TD><var>prop</var></TD> <TD><var>prop1</var></TD> </TR>
    <TR> <TD>resonance 2</TD> <TD>res2</TD> <TD>-</TD> <TD>-</TD> <TD><var>prop2</var></TD> </TR>
    <TR> <TD>resonance 3</TD> <TD>res3</TD> <TD>-</TD> <TD>-</TD> <TD>-</TD> </TR>
    <TR> <TD>1st leg</TD> <TD>a</TD> <TD>-</TD> <TD><var>i</var></TD> <TD><var>k</var></TD> </TR>
    <TR> <TD>2nd leg</TD> <TD>b</TD> <TD>-</TD> <TD><var>j</var></TD> <TD><var>i</var></TD> </TR>
    <TR> <TD>3rd leg</TD> <TD>c</TD> <TD>-</TD> <TD>-</TD> <TD><var>j</var></TD> </TR>
    <TR> <TD>4th leg</TD> <TD>d</TD> <TD>-</TD> <TD>-</TD> <TD>-</TD> </TR>
    </TABLE>
    */
  /*!
    \var ChannelInformation::name
    Type of the integrator:
    - \c Isotripic
    - \c Dalitz
    - \c TwoResonances
    .
    */
  /*!
    \var ChannelInformation::nRes
    Number of the resonances (internal propagators)
    */
  /*!
    \var ChannelInformation::res1
    Name of the resonance 1
    */
  /*!
    \var ChannelInformation::res2
    Name of the resonance 2
    */
  /*!
    \var ChannelInformation::res3
    Name of the resonance 3
    */
  /*!
    \var ChannelInformation::a
    Number of the particle in the 1st leg
    */
  /*!
    \var ChannelInformation::b
    Number of the particle in the 2nd leg
    */
  /*!
    \var ChannelInformation::c
    Number of the particle in the 3rd leg
    */
  /*!
    \var ChannelInformation::d
    Number of the particle in the 4th leg
    */
    

  /*!
    \class HD_PS_Base
    \brief Base class for the phase space generator

    This class is a subclass of PHASIC::MultiChannel
  */    
  /*!
    \fn HD_PS_Base::HD_PS_Base(Hadron_Decay_Channel * hdc, std::vector<std::string> & _pst, bool & mustinit, struct GeneralModel &_locmd, bool read_dc )
    \brief Constructor and initialisation of private attributes

    The constructor calls HD_PS_Base::Construct to set up the PS.
    The argument <var>mustinit</var> becomes \c true if the PS of the considered decay channel
    needs to be initialised (HD_PS_Base::Initialise, i.e. integration and finding the maximums), 
    otherwise it is \c false.
  */    
  /*!
    \fn HD_PS_Base::Construct ( GeneralModel & _md )
    \brief Reads the decay channel file

    This method reads the DC file, such as \f$\tau \to \pi \nu_\tau\f$ and
    extracts the stored information about:
    - \b Channels: 
        - The switch \c AlwaysIntegrate: 1 means that this channel will be reintegrated; 0 means that
          this channel will only be reintegrated if no results from a previous run are given.
        - Reads the channel name(s) and its (their) weight factor(s) 
          (if no weight factor is given it is set to 1).  
        - These channels are piped into HD_PS_Base::AddChannel which is tasked with setting pointers right. 
        .
    - \b Resonances: which resonance does decay into which particles.
    - \b Parameters: given parameters for the decay amplitude
      These values are stored in the map <var>_md</var> (HD_ME_Base::GeneralModel)
      If the map contains a certain entry already (e.g. because this parameter was given in
      \c HadronConstants.dat) the value of this parameter will be overwritten (only for
      this channel of course).
      Parameters may also be complex. 
    - \b Result: result, error, maximum. These three values are stored in the private attributes 
      and reused for the run. That is, if these values are given, there is no need to reintegrate.
      - result: the total decay width
      - error: its error
      - maximum: the maximal value of \f$|M|^2 \times dlips\f$ in the complete phase space
      .
    .
    If results are written in the DC file this method return \c false otherwise \c true. This information
    is necessery to whether the integrator must be initialised or not. Moreover, if the
    switch \c AlwaysIntegrate is switched-on this method returns \c true as well (then a channel
    becomes reintegrated no matter if some results are written out or not)
  */    
  /*!
    \fn HD_PS_Base::AddChannel(std::string,double=1.)
    \brief Set pointers right

    This method calls HD_Channel_Selector::GetChannel and creates a pointer which after setting the
    integration channel weight according to the weight factor given is used to finally add this
    (integration!) channel to the vector
    of channels in PHASIC::Multi_Channel class (PHASIC::Multi_Channel::Add).
  */
  /*!
    \fn void HD_PS_Base::Initialise()
    \brief Initialisation of the PS integrator

    This method uses HD_PS_Base::CalculateNormalisedWidth and HD_PS_Base::WriteOut to
    do that job.
  */    
  /*!
    \fn void HD_PS_Base::CalculateNormalisedWidth()
    \brief Integrates the decay channel

    This method calculates the total decay width of the decay channel. It uses Vegas Optimisation 
    if it is switched-on. Furthermore, it determines the maximum value of the channel. 
    Basically, it calls the method Hadron_Decay_Channel::Differential for getting the
    value and the corresponding phase space point.
  */    
  /*!
    \fn HD_PS_Base::WriteOut()
    \brief Writes out a complete DC file

    This method genererates the the DC file corresponding to the considered decay channel.
    The file is written out and the original file is stored with the ending \c .dat. This
    becomes necessary when a channel has  to be reintegrated so that the result, error,
    and maximum might have been changed.
  */    
  
  
} // end of namespace


#endif
