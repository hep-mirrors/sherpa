#include "METOOLS/Explicit/Current.H"
#include "METOOLS/Currents/C_Pseudo.H"
/*!
  @file P_C.C
  @brief Implements the class CP.
*/

namespace METOOLS {

  template <typename SType>
  class CP: public Current,
	    public Current_Contractor<SType> {
  public:

    typedef std::complex<SType>   SComplex;
    typedef std::vector<SComplex> SComplex_Vector;

    typedef CAsT4<SType> CAsT4Type;
    typedef std::vector<CAsT4Type*> CAsT4Type_Vector;

  protected:

    /// @brief propagator factor
    SComplex m_prop;

    /// @brief returns the label for the graph output
    std::string CLabel() const;

  public:

    /// @brief constructor using Current_Key key
    CP(const Current_Key &key);

    /// @brief constructs current
    void ConstructJ(const ATOOLS::Vec4D &p,const int ch,
		    const int cr,const int ca,const int mode);

    /// @brief sets the gauge vector
    void SetGauge(const ATOOLS::Vec4D &k);

    /// @brief adds the propagator to the current
    void AddPropagator();

    /// @brief contracts a current
    void SContract
    (const Current &c,const Int_Vector &pols,
     SComplex_Vector &ress,const size_t &offset) const;

    /// @brief formats the current into a string
    std::string Format(const CObject *c) const;

    /// @brief returns the current type
    char Type() const;

  };// end of class CP
  /*!
    @class CP
    @brief This class represents a pseudo-scalar current.

    This class represents a pseudo-scalar current.
  */

}// end of namespace METOOLS

#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/STL_Tools.H"
#include "ATOOLS/Org/MyStrStream.H"

#define M_I SComplex(0.0,1.0)

using namespace METOOLS;
using namespace ATOOLS;

template <typename SType>
CP<SType>::CP(const Current_Key &key): 
  Current(key), m_prop(-M_I) 
{
}

template <typename SType>
void CP<SType>::ConstructJ(const ATOOLS::Vec4D &p,const int ch,
			   const int cr,const int ca,const int mode)
{
}

template <typename SType>
void CP<SType>::SetGauge(const ATOOLS::Vec4D &k)
{
}

template <typename SType>
void CP<SType>::AddPropagator()
{
  // add propagator for off-shell leg
  for (size_t i(0);i<m_j.size();++i) {
  CAsT4Type_Vector *j(m_j[i].template Get<CAsT4Type>());
  for (typename CAsT4Type_Vector::iterator 
	 jit(j->begin());jit!=j->end();++jit) **jit*=m_prop;
  }
}

template <typename SType> void CP<SType>::SContract
(const Current &c,const Int_Vector &pols,
 SComplex_Vector &ress,const size_t &offset) const
{
  THROW(fatal_error,"Multiplication of tensor particles not allowed");
}

template <typename SType>
std::string CP<SType>::Format(const CObject *c) const
{
  return ToString(*(CAsT4Type*)c,6);
}

template <typename SType>
char CP<SType>::Type() const
{
  return 'P';
}

template <typename SType>
std::string CP<SType>::CLabel() const
{
  return "double,label.side=right,label=$"+
    (this->m_out.empty()?this->m_fl.Bar():this->m_fl).TexName()+"$";
}

DECLARE_GETTER(CP<double>,"DP",Current,Current_Key);

Current *ATOOLS::Getter<Current,Current_Key,CP<double> >::
operator()(const Current_Key &key) const
{
  if (key.m_fl.IsTensor() && key.m_fl.IsDummy()) 
    return new CP<double>(key);
  return NULL;
}

void ATOOLS::Getter<Current,Current_Key,CP<double> >::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"tensor current (double)";
}
