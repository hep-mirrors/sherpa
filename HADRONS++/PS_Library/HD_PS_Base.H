#ifndef HADRONS_PS_Library_HD_PS_Base_H
#define HADRONS_PS_Library_HD_PS_Base_H

#include "PHASIC++/Channels/Single_Channel.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/Vector.H"
#include <string>
#include "HADRONS++/Main/Tools.H"

namespace HADRONS {

  struct ChannelInformation
  {
    std::string name;                // name of integrator
    int nRes;                        // number of resonances
    std::string res1, res2, res3;    // resonances
    int a,b,c,d;                    // index for outgoing particles
  };

  class HD_Channel_Selector {
    static bool DecomposeChannel( std::string name, ChannelInformation & ci );
  public:
    static PHASIC::Single_Channel*
    GetChannel(int nin, int nout, const ATOOLS::Flavour * flavs,
               std::string name, GeneralModel const &,
               const ATOOLS::Mass_Selector* ms);
  };
  
  /*!
    \file HD_PS_Base.H
    \brief Declares the classes HADRONS::HD_Channel_Selector and HADRONS::ChannelInformation.

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
    
} // end of namespace


#endif
