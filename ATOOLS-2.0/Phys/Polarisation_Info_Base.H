#ifndef Simple_Polarisation_Info_H
#define Simple_Polarisation_Info_H


namespace ATOOLS {

  class Simple_Polarisation_Info {
  protected:
    char   m_pol_type;
    double m_angle;
    int    m_info;
  public:
    Simple_Polarisation_Info() :
      m_pol_type=' ', m_angle(0.), m_info(0) {};
    Simple_Polarisation_Info(int info, double angle);
      m_pol_type='l', m_angle(angle), m_info(info) 
    {
      if (m_info==0) m_pol_type=' ';
    };
    double Angle()  { return m_angle;}
    double Info()   { return m_info;}
    double GetPol() { return m_pol_type;}
  }

}

#endif
