#ifndef DIRE__Tools__Parton_H
#define DIRE__Tools__Parton_H

#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/Vector.H"
#include "DIRE/Tools/Color.H"

namespace DIRE {

  class Amplitude;

  class Parton {
  public:

    struct Weight {
      double m_t, m_w;
      inline Weight(const double &t=-1.0,
		    const double &w=0.0):
	m_t(t), m_w(w) {}
    };// end of struct Weight

    typedef std::vector<Weight> Weight_Vector;

    typedef std::map<const Parton*,Weight_Vector> Weight_Map;

  private:

    Amplitude *p_ampl;

    ATOOLS::Flavour m_f;
    ATOOLS::Vec4D   m_p;

    Color m_c;
    int   m_h, m_b;

    Parton *p_out[2], *p_in;

    size_t m_id;

    Weight_Map m_ws;

    static size_t s_cnt;

  public:

    Parton(Amplitude *const ampl,const ATOOLS::Flavour &f,
	   const ATOOLS::Vec4D &p=ATOOLS::Vec4D(),
	   const Color &c=Color(),const int h=0);

    ~Parton();

    double GetXB() const;

    void AddWeight(const Parton *s,const double &t,const double &w);

    double GetWeight(const double &t) const;

    void SetColor(const Color &c);

    inline const ATOOLS::Flavour &Flav() const { return m_f; }
    inline const ATOOLS::Vec4D   &Mom() const  { return m_p; }
    inline const Color           &Col() const  { return m_c; }

    inline void SetFlav(const ATOOLS::Flavour &f) { m_f=f; }
    inline void SetMom(const ATOOLS::Vec4D &p)    { m_p=p; }
    inline void SetCol(const Color &c)            { m_c=c; }

    inline int  Hel() const   { return m_h; }
    inline void SetHel(int h) { m_h=h;      }

    inline int  Beam() const   { return m_b; }
    inline void SetBeam(int b) { m_b=b;      }

    inline size_t Id() const              { return m_id; }
    inline void   SetId(const size_t &id) { m_id=id;     }

    inline Parton *Out(int i) const              { return p_out[i]; }
    inline void    SetOut(int i,Parton *const p) { p_out[i]=p;      }

    inline Parton *In() const             { return p_in; }
    inline void    SetIn(Parton *const p) { p_in=p;      }

    inline Amplitude *Ampl() const { return p_ampl; }

    inline static size_t Counter() { return s_cnt; }

    inline const Weight_Map &Weights() const { return m_ws; }

    inline void ClearWeights() { m_ws.clear(); }

  };// end of class Parton

  std::ostream &operator<<(std::ostream &s,const Parton &p);

  typedef std::vector<Parton*> Parton_Vector;

}// end of namespace DIRE

#endif
