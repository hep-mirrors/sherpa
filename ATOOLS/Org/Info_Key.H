#ifndef ATOOLS_Org_Info_Key_H
#define ATOOLS_Org_Info_Key_H

#include "ATOOLS/Org/Smart_Pointer.H"
#include "ATOOLS/Math/Vector.H"
#include <string>
#include <vector>
#include <iostream>

namespace ATOOLS {

  typedef std::vector<double>        Double_Container;
  typedef std::vector<ATOOLS::Vec4D> Vector_Container;

  namespace si {

    enum code { 
      idle      =  0,
      reset     =  1,
      dicing    =  2,
      diced     =  3,
      weighting =  4,
      weighted  =  5,
      error     =  6  
    };

  }// end of namespace si

  class Integration_Info;

  class Info_Key {
  private:

    Integration_Info *p_info;
    std::string       m_name, m_info;
    
    size_t m_valuekey, m_weightkey;

    friend class Integration_Info;

  public:

    // constructors
    Info_Key();

    // destructor
    ~Info_Key();

    // member functions
    void SetInfo(const std::string info);

    // setters
    inline si::code SetStatus(const si::code sicode) const;

    // getters
    inline const std::string &Name() const;
    inline const std::string &Info() const;

    si::code Status() const;

    // member functions
    void Assign(const std::string name,const size_t doubles,
		const size_t vectors,Integration_Info *const info);
    void Assign(const std::string name,const size_t doubles,
		const size_t vectors,const SP(Integration_Info) &info);

    inline Double_Container &Doubles() const;
    inline Vector_Container &Vectors() const;

    inline double        Double(const size_t i) const;
    inline ATOOLS::Vec4D Vector(const size_t i) const;
    inline double        Weight() const;

    inline bool Assigned() const;

    inline double        &operator[](const size_t i);
    inline ATOOLS::Vec4D &operator()(const size_t i);

    inline void operator<<(const double weight);

    inline Info_Key &operator=(const Info_Key &key);
    inline bool      operator==(const Info_Key &key);

    friend std::ostream &operator<<(std::ostream &str,const Info_Key &key);

  };// end of class Info_Key

  /*!
    \class Info_Key
    \brief Handles access to Integration_Info

    This class provides access to Integration_Info which handles value and 
    weight storage. Each Info_Key is assigned a reference point in Integration_Info 
    during Assign(<name>,<double_size>,<vector_size>,<info>).
    In the following the key calls Integration_Info to obtain and to set values
    and weights and to test for the status of the corresponding variable (see 
    the example below).
    \code
       class Example_ISR_Channel: public Single_Channel {
       private:
         double   m_exponent
         Info_Key m_spkey, ...;
       public:
         ...
       };

       void Example_ISR_Channel::GeneratePoint(const double *rans)
       {
         // m_spkey owns three values where m_spkey[0] -> lower limit
	 //                                 m_spkey[1] -> upper limit
	 //                                 m_spkey[2] -> value
         m_spkey[2]=CE::MasslessPropMomenta(m_exponent,m_spkey[0],m_spkey[1],ran[0]);
	 ...
       }

       void Example_ISR_Channel::GenerateWeight()
       {
         // operator<<(..) sets the unique weight
         m_spkey<<1./CE::MasslessPropWeight(m_exponent,m_spkey[0],m_spkey[1],m_spkey[2]);
	 ...
       }
    \endcode
  */

}// end of namespace ATOOLS

#include "ATOOLS/Org/Info_Key.inl.H"

#endif
