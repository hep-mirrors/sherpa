#ifndef PI_Interface_H
#define PI_Interface_H

#include "Primitive_Integrator.H"
#include "Phase_Space_Handler.H"

namespace PHASIC {

  class PI_Interface: public ATOOLS::Primitive_Integrand {
  private:

    ATOOLS::Primitive_Integrator *p_integrator;

    Phase_Space_Handler *p_pshandler;

    std::string m_key;
    bool        m_initialize, m_stop;
    psm::code   m_mode, m_cmode;

    std::vector<double> m_point;
    double              m_integral, m_value;

    size_t m_nopt, m_nmax, m_ncells;

  public:

    // constructor 
    PI_Interface(Phase_Space_Handler *const pshandler,
		 const std::string &key,const size_t dim);

    // destructor
    ~PI_Interface();

    // member functions
    bool Initialize();

    void   GeneratePoint(const psm::code &mode=psm::normal_call);
    double GenerateWeight();

    double operator()(const std::vector<double> &x);
    void   AddPoint(const double weight);
    void   FinishConstruction(const double apweight);

    bool WriteOut(const std::string &path) const;
    bool ReadIn(const std::string &path);

    // inline functions
    inline const std::string         &Key() const   { return m_key;   }
    inline const std::vector<double> &Point() const { return m_point; }

    inline void SetMode(const psm::code mode) { m_mode=mode;   }
    inline psm::code Mode() const             { return m_mode; }

    inline double Integral() const { return m_integral; }
    inline double Value() const    { return m_value;    }

    inline bool Stop() const { return m_stop; }

    inline ATOOLS::Primitive_Integrator *operator->() const 
    { return p_integrator; }

  };// end of class PI_Interface

}// end of namespace PHASIC

#endif
