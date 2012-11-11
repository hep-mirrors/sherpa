#include "METOOLS/Main/Spin_Structure.H"

namespace METOOLS {
  class Current;
  class Vertex;
}

namespace EXTRAXS {
  class Comix1to3 : public METOOLS::Spin_Amplitudes {
    std::vector<METOOLS::Current*> m_cur;
    METOOLS::Current *m_scur, *m_fcur;
    METOOLS::Vertex *m_v1, *m_v2;

    std::vector<METOOLS::Current*> m_anticur;
    METOOLS::Current *m_antiscur, *m_antifcur;
    METOOLS::Vertex *m_antiv1, *m_antiv2;

    std::vector<size_t> m_nhel;
    
    ATOOLS::Flavour m_prop;

    size_t NHel(const ATOOLS::Flavour& fl);
    METOOLS::Vertex* GetVertex(METOOLS::Current* cur1,
                                  METOOLS::Current* cur2,
                                  METOOLS::Current* prop);      
  public:
    Comix1to3(const std::vector<ATOOLS::Flavour>& flavs,
              const ATOOLS::Flavour& prop,
              size_t nonprop, size_t propi, size_t propj);
    ~Comix1to3();
    void Calculate(const ATOOLS::Vec4D_Vector& momenta, bool anti=false);
    ATOOLS::Flavour Prop() const { return m_prop; }
  };
}

