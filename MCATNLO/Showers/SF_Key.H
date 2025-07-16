#ifndef MCATNLO_Showers_SF_Key_H
#define MCATNLO_Showers_SF_Key_H

#include <string>

namespace ATOOLS {
  class Mass_Selector;
}

namespace MODEL {
  class Model_Base;
  class Single_Vertex;
}

namespace MCATNLO {

  class SF_Coupling;

  struct cstp {
    enum code {
      FF   = 11,
      FI   = 12,
      IF   = 21,
      II   = 22,
      none = 0
    };
  };

  struct SF_Key {
  public:
    const ATOOLS::Mass_Selector *&p_ms;
    const MODEL::Single_Vertex *p_v;
    int m_mode, m_kfmode, m_ewmode, m_col;
    cstp::code m_type;
    SF_Coupling *p_cf;
    std::pair<double, double> m_pdfmin;
  public:
    inline SF_Key(const ATOOLS::Mass_Selector *&ms,
		  const MODEL::Single_Vertex *v,
		  const int mode,const cstp::code type,
		  const int kfmode,const int ewmode,const int col,
                  std::pair<double, double>& pdfmin):
      p_ms(ms), p_v(v), m_mode(mode), m_kfmode(kfmode),
      m_ewmode(ewmode), m_col(col), m_type(type),
      p_cf(NULL), m_pdfmin{ pdfmin } {}
    std::string ID(const int mode) const;
  };

  std::ostream &operator<<(std::ostream &str,const SF_Key &k);

}

#endif
