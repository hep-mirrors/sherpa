#ifndef COMIX_Phasespace_PS_Current_H
#define COMIX_Phasespace_PS_Current_H

#include "COMIX/Phasespace/PS_Vertex.H"

namespace COMIX {

  class PS_Current: public Current_Base {
  public:

    typedef std::vector<Complex> Complex_Vector;

  protected:

    std::vector<PS_Info> m_j;

    Current_Base *p_scc;

  public:

    // constructor
    inline PS_Current(const Current_Key &key):
      Current_Base(key), p_scc(0) {}

    // member functions
    void ResetJ();
    void AddJ(const PS_Info &j);

    const std::vector<PS_Info> &J() const;
    
    void ConstructJ(const ATOOLS::Vec4D &p,const int ch,
		    const int cr,const int ca);
    void SetGauge(const ATOOLS::Vec4D &k);

    void AddPropagator();

    Complex Contract(const Current_Base &c,
		     const long unsigned int &hm,
		     const long unsigned int &hp) const;
    void Contract(const Current_Base &c,
		  const LongUIntMap_Matrix &hmp,
		  Complex_Vector &ress) const;

    char Type() const;    

    void Evaluate();

    void Print() const;

    void ResetZero();

    void SetSCC(Current_Base *const scc);

    // inline functions
    inline void SetMass(const double & mass)   { m_mass=mass;   }
    inline void SetWidth(const double & width) { m_width=width; }

    inline Current_Base *SCC() const { return p_scc; }

  };// end of class PS_Current

}// end of namespace COMIX

#endif
