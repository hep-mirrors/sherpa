#ifndef HADRONS_ME_Library_HD_ME_Base_H
#define HADRONS_ME_Library_HD_ME_Base_H

#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/Vector.H"
#include <string>
#include <utility>
#include "HELICITIES/Main/Spin_Structure.H"
#include "HADRONS++/Main/Tools.H"
#include "ATOOLS/Org/Getter_Function.H"

#define DEFINE_ME_GETTER(CLASS,NAME,TAG)                              \
DECLARE_GETTER(NAME,TAG,HD_ME_Base,Flavour_Info);                     \
HD_ME_Base* NAME::operator()(const Flavour_Info &parameters) const    \
{ return new CLASS(parameters.flavs, parameters.n, parameters.indices,TAG); }

namespace HADRONS {

  class HD_ME_Base {
  protected:
    std::string       m_name;
    int               m_n;
    ATOOLS::Flavour * m_flavs;
    double          * p_masses,* p_masses2;
    std::string       m_path;
    bool              m_anti;
    int             * p_i;  // "external" indices

    virtual double lambdaNorm(const double,const double,const double);
  public:
    HD_ME_Base(ATOOLS::Flavour* flavs, int n, int* decayindices, std::string name);
    virtual ~HD_ME_Base();
    void SetPath( std::string s ) { m_path=s; }
    virtual std::string Name() { return m_name; }
    void SetAnti(bool anti) { m_anti = anti; }

    virtual void   operator()(const ATOOLS::Vec4D* p,
                              HELICITIES::Spin_Amplitudes* amps) = 0;
    virtual void   AddVector(double,double) {}
    virtual void   SetModelParameters( GeneralModel _md ) {};
    virtual bool   SetColorFlow(std::vector<ATOOLS::Particle*> outparts,int n_q, int n_g);
  }; // end of class HD_ME_Base

  typedef ATOOLS::Getter_Function<HD_ME_Base,Flavour_Info> HD_ME_Getter_Function;

  /*!
  \file HD_ME_Base.H
  \brief Declares the classes HADRONS::HD_ME_Base as well as the structure
    HADRONS::GeneralModel.

  This file can be found in the directory \c ME_Library.

  */

  /*!
    \class GeneralModel
    \brief A map with all relevant parameters for a certain decay amplitude

    It is a map that both contains the name of the parameter (as it is written in the
    DC file) and its value:
    <TABLE>
    <TR> <TD>1st column</TD><TD>2nd column</TD> </TR>
    <TR> <TD>name</TD><TD>value</TD> </TR>
    </TABLE>
    */
  /*!
    \fn GeneralModel::operator()(const std::string &tag,const double &def) const 
    \brief Returns the value of the parameter with the name <var>name</var>

    Returns the value of the parameter with the name <var>tag</var>. If such a parameter
    does not exist, it returns the default value given by <var>def</var>.
    */
    

  /*!
  \class HD_ME_Base
  \brief Base class for ME classes
     
    \section spin The Implemention of Spin Correlations

    At this stage, the implementation of the spin correlation algorithm will be explained a bit.
    Each ME class (Tau_Lepton, Tau_Pseudo, ...) has to set up an HELICITIES::Amplitude_Tensor.
    This is a structure where an entry stands for a value of the ME (\e not squared) with
    a certain helicity combination. That is, each helicity combination has its own entry in 
    the Amplitude_Tensor. To make it easy to insert an amplitude for a specific helicity
    combination, it provides a number of methods. Have a look at the documentation in
    HELICITIES::Amplitude_Tensor for them.
  */

  /*!
  \fn HD_ME_Base::HD_ME_Base( int,ATOOLS::Flavour *)
  \brief Constructor and initialisation of private attributes
  */

  /*!
    \fn HD_ME_Base::operator()(const Vec4D* p, Amplitude_Tensor * amps, int k0_n )
    \brief Main operator returning the Amplitude_Tensor of the decay amplitude

    This operator is the most important one. Therefore, each ME class <b>has to have</b>
    such an operator. It fills the given Amplitude_Tensor with one amplitude for each
    possible spin combination.
  */

  /*!
    \fn HD_ME_Base::SetModelParameters( GeneralModel _md )
    \brief This method sets the parameters for the decay amplitude

    Every subclass of HD_ME_Base (the ME class) \b must have this method if it needs
    parameters that are written in the DC file.
    Basically, this method gets called after reading the corresponding DC file and it is
    its job to set the parameters in the ME class to the values written in the DC file.
    This is done by passing the variable <var>_md</var> (see GeneralModel) 
    which contains the values and names
    of all relevant parameters.
    */
  /*!
  \var HD_ME_Base::p_masses
  Pointer of masses.
  */
  /*!
  \var HD_ME_Base::p_masses2
  Pointer of masses squared.
  */
  /*!
  \var HD_ME_Base::p_flavs
  Pointer on flavours.
  */
  /*!
  \var HD_ME_Base::m_nout
  Number of outgoing particles.
  */
  /*!
  \var HD_ME_Base::m_metype
  Type of ME.
  */

} // end of namespace


#endif
