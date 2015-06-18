#ifndef DIRE__Tools__Weight_H
#define DIRE__Tools__Weight_H

namespace DIRE {

  struct Weight {
  public:

    double m_f, m_g, m_h;

    inline Weight(const double &f,
		  const double &g,const double &h):
      m_f(f), m_g(g), m_h(h) {}

    inline double MC() const { return m_f/m_g; }
    inline double Accept() const { return m_g/m_h; }
    inline double Reject() const
    { return m_g/m_h*(m_h-m_f)/(m_g-m_f); }

  };// end of struct Weight

  inline std::ostream &operator<<(std::ostream &s,const Weight &w)
  { return s<<"W(f="<<w.m_f<<",g="<<w.m_g<<",h="<<w.m_h<<")"; }

}// end of namespace DIRE

#endif
