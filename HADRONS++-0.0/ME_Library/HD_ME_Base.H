#ifndef HD_ME_Base_H
#define HD_ME_Base_H

#include "Flavour.H"
#include "Vector.H"
#include <string>
#include <utility>
#include "Message.H"
#include "Tools.H"
#include "Getter_Function.H"

namespace HADRONS {

  class HD_ME_Base {
  protected:
    std::string       m_metype;
    int               m_nout;
    ATOOLS::Flavour * p_flavs;
    double          * p_masses,* p_masses2;
	std::string       m_path;
  public:
    HD_ME_Base(int,ATOOLS::Flavour *);
    virtual ~HD_ME_Base();
    void SetPath( std::string s ) { m_path=s; }
	std::string METype() { return m_metype; }

    virtual void   operator()( 
        const ATOOLS::Vec4D      * _p, 
        std::vector<Complex>     * _ampls_tensor, 
        std::vector<std::pair<int,int> > * _indices,
        int k0_n )
    {
      _ampls_tensor->clear();
      _indices->clear();
    }
    virtual void   AddVector(double,double) {}
	virtual void   SetModelParameters( struct GeneralModel _md ) {};
  }; // end of class HD_ME_Base

  typedef ATOOLS::Getter_Function<HD_ME_Base,Flavour_Info> HD_ME_Getter_Function;

  class Isotropic : public HD_ME_Base {
  public:
    Isotropic(int,ATOOLS::Flavour *,
	      std::string _met=std::string(""));
  }; // end of class Isotropic
  
  /*!
	\file HD_ME_Base.H
	\brief Declares the classes HADRONS::HD_ME_Base and HADRONS::Isotropic as well as the structure
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
    Each ME class (Tau_Lepton, Tau_Pseudo, ...) has to set up a so-called \b amplitude-tensor.
    This is a C++ vector where an entry stands for a value of the ME (\e not squared) with 
    a certain helicity combination. That is, each helicity combination has its own entry in 
    the amplitude tensor. Since the order plays an important role because these values are to
    be weighted with a spin density matrix, it is essential to create an <b>index list</b> that
    tells where which index in the tensor is located. This will be explained in the following example:
    Imagine a 3-body decay of the form
    \f[ M(0) \to 1 + 2 + 3, \f] where each particle has an helicity \f$\alpha_i\f$ with
    \f$0 \leq \alpha_i \leq m_i-1\f$. Particles with spin 0 do not have an \f$\alpha_i\f$ and therefore
    no index in the amplitude tensor. Particles with spin \f$s\f$ usually have \f$m_i=2s+1\f$.
    Let's say that the values if the ME's are written into the amplitude tensor in the following order
    <TABLE>
    <TR> <TD>entry</TD> 
         <TD>\f$\alpha_3\f$</TD><TD>\f$\alpha_0\f$</TD><TD>\f$\alpha_2\f$</TD><TD>\f$\alpha_1\f$</TD> </TR>
    <TR> <TD>0</TD> <TD>0</TD><TD>0</TD><TD>0</TD><TD>0</TD> </TR>
    <TR> <TD>1</TD> <TD>0</TD><TD>0</TD><TD>0</TD><TD>1</TD> </TR>
    <TR> <TD>2</TD> <TD>0</TD><TD>0</TD><TD>0</TD><TD>2</TD> </TR>
    <TR> <TD>\f$\vdots\f$</TD> <TD>0</TD><TD>0</TD><TD>0</TD><TD>\f$\vdots\f$</TD> </TR>
    <TR> <TD>\f$m_1-1\f$</TD> <TD>0</TD><TD>0</TD><TD>0</TD><TD>\f$m_1-1\f$</TD> </TR>
    <TR> <TD>\f$m_1\f$</TD> <TD>0</TD><TD>0</TD><TD>1</TD><TD>0</TD> </TR>
    <TR> <TD>\f$\vdots\f$</TD> <TD>0</TD><TD>0</TD><TD>1</TD><TD>\f$\vdots\f$</TD> </TR>
    <TR> <TD>\f$2m_1-1\f$</TD> <TD>0</TD><TD>0</TD><TD>1</TD><TD>\f$m_1-1\f$</TD> </TR>
    <TR> <TD>\f$2m_1\f$</TD> <TD>0</TD><TD>0</TD><TD>2</TD><TD>0</TD> </TR>
    <TR> <TD>\f$\vdots\f$</TD> <TD>0</TD><TD>0</TD><TD>\f$\vdots\f$</TD><TD>\f$\vdots\f$</TD> </TR>
    <TR> <TD>\f$m_2m_1-1\f$</TD> <TD>0</TD><TD>0</TD><TD>\f$m_2-1\f$</TD><TD>\f$m_1-1\f$</TD> </TR>
    <TR> <TD>\f$m_2m_1\f$</TD> <TD>0</TD><TD>1</TD><TD>0</TD><TD>0</TD> </TR>
    <TR> <TD>\f$\vdots\f$</TD> <TD>0</TD><TD>\f$\vdots\f$</TD><TD>\f$\vdots\f$</TD><TD>\f$\vdots\f$</TD> </TR>
    <TR> <TD>\f$m_0m_2m_1-1\f$</TD> <TD>0</TD><TD>\f$m_0-1\f$</TD><TD>\f$m_2-1\f$</TD><TD>\f$m_1-1\f$</TD> </TR>
    <TR> <TD>\f$m_0m_2m_1\f$</TD> <TD>1</TD><TD>0</TD><TD>0</TD><TD>0</TD> </TR>
    <TR> <TD>\f$\vdots\f$</TD> <TD>\f$\vdots\f$</TD><TD>\f$\vdots\f$</TD><TD>\f$\vdots\f$</TD><TD>\f$\vdots\f$</TD> </TR>
    <TR> <TD>\f$m_3m_0m_2m_1-1\f$</TD> <TD>\f$m_3-1\f$</TD><TD>\f$m_0-1\f$</TD><TD>\f$m_2-1\f$</TD><TD>\f$m_1-1\f$</TD> </TR>
    </TABLE>
    The index list corresponds to the list with "two columns" (it is a C++ vector of std::pair). The
    first column corresponds to the particle and the second to its spin times two. For this example it must
    have four entries (because there are four indices in the amplitude tensor \f$\alpha_0 \dots \alpha_3\f$):
    <TABLE>
    <TR> <TD>entry</TD> <TD>particle</TD> <TD>spin</TD> </TR>
    <TR> <TD>0</TD> <TD>1</TD> <TD>\f$2s_1 = m_1-1\f$</TD> </TR>
    <TR> <TD>1</TD> <TD>2</TD> <TD>\f$2s_2 = m_2-1\f$</TD> </TR>
    <TR> <TD>2</TD> <TD>0</TD> <TD>\f$2s_0 = m_0-1\f$</TD> </TR>
    <TR> <TD>3</TD> <TD>3</TD> <TD>\f$2s_3 = m_3-1\f$</TD> </TR>
    </TABLE>
    It totally depends on the programming in which order the amplitude tensor is filled. Hence, 
    it the duty of the programmer to tell which index is to be found at which position.
    If the ME class does return an \b empty index list, the decay is treated without 
    spin correlations.
  */	
  /*!
	\fn HD_ME_Base::HD_ME_Base( int,ATOOLS::Flavour *)
	\brief Constructor and initialisation of private attributes
  */	
  /*!
    \fn HD_ME_Base::operator()( const ATOOLS::Vec4D * _p, std::vector<Complex> * _ampls_tensor, std::vector<std::pair<int,int> > * _indices, int k0_n )
    \brief Main operator returning the amplitude tensor of the decay amplitude

    This operator is the most important one. Therefore, each ME class <b>has to have</b> such an operator,
    otherwise Hadron_Decay_Channel::Differential calls an operator which returns basically nothing.
    It needs a pointer on the amplitude tensor and index list to be filled.
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

  /*!
	\class Isotropic
	\brief Not worth mentioning

	This class is a subclass of HADRONS::HD_ME_Base.
  */	

} // end of namespace


#endif
