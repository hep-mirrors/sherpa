#ifndef COMIX_Phasespace_PS_Current_H
#define COMIX_Phasespace_PS_Current_H

#include "COMIX/Phasespace/PS_Vertex.H"

namespace COMIX {

  class PS_Current: public Current {
  public:

    typedef std::vector<Complex> Complex_Vector;

  protected:

    Current *p_scc;

  public:

    // constructor
    inline PS_Current(const Current_Key &key):
      Current(key), p_scc(0) { m_j.resize(1); }

    // member functions
    void ConstructJ(const ATOOLS::Vec4D &p,const int ch,
		    const int cr,const int ca,const int mode);
    void SetGauge(const ATOOLS::Vec4D &k);

    void AddPropagator();

    std::string Format(const CObject *c) const;

    char Type() const;    

    void SetSCC(Current *const scc);

    // inline functions
    inline void SetMass(const double & mass)   { m_mass=mass;   }
    inline void SetWidth(const double & width) { m_width=width; }

    inline Current *SCC() const { return p_scc; }

  };// end of class PS_Current

}// end of namespace COMIX

#endif
