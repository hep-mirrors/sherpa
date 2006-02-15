#ifndef HD_ME_Base_H
#define HD_ME_Base_H

#include "Flavour.H"
#include "Vector.H"
#include <string>
#include <utility>

namespace HADRONS {

  typedef std::map<std::string,double> String_Map;

  struct GeneralModel: public String_Map {
  public:
    inline double operator()(const std::string &tag,const double &def) const 
    {
      String_Map::const_iterator fit(find(tag));
      return fit!=end()?fit->second:def;
    }
  }; 

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
        int                        k0_n )
    {
      _ampls_tensor->clear();
      _indices->clear();
    }
    virtual void   AddVector(double,double) {}
	virtual void   SetModelParameters( struct GeneralModel _md ) {};
  }; // end of class HD_ME_Base

  class Isotropic : public HD_ME_Base {
  public:
    Isotropic(int,ATOOLS::Flavour *,
	      std::string _met=std::string(""));
  }; // end of class Isotropic
  
  /*!
	\file HD_ME_Base.H
	\brief Declares the classes HADRONS::HD_ME_Base and HADRONS::Isotropic

	This file can be found in directory \c ME_Library.
  */

  /*!
	\class HD_ME_Base
	\brief Tool to set up the connection to the corresponding ME method
  */	
  /*!
	\fn HD_ME_Base::HD_ME_Base
	\brief Constructor and initialisation of private attributes
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
