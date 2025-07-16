#ifndef ATOOLS_Phys_Simple_Polarisation_Info_H
#define ATOOLS_Phys_Simple_Polarisation_Info_H

#include <iostream>

namespace ATOOLS {

  class Simple_Polarisation_Info {
  protected:
    char   m_pol_type;
    double m_angle;
    int    m_info;
  public:
    Simple_Polarisation_Info() :
      m_pol_type(' '), m_angle(0.), m_info(0) {};
    Simple_Polarisation_Info(int info, double angle) :
      m_pol_type('l'), m_angle(angle), m_info(info) 
    {
      if (m_info==0) m_pol_type=' ';
    };
    double Angle()  const { return m_angle;}
    int    Info()   const { return m_info;}
    char   GetPol() const { return m_pol_type;}
  };

  std::ostream& operator<< (std::ostream& s, const Simple_Polarisation_Info & polinfo);

}

#endif
