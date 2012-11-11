#include "METOOLS/Main/Spin_Structure.H"

namespace METOOLS {
  class Current;
  class Vertex;
}

namespace EXTRAXS {

  class Comix1to2 : public METOOLS::Spin_Amplitudes {
    std::vector<METOOLS::Current*> m_cur;
    METOOLS::Current *m_fcur;
    METOOLS::Vertex *m_v1;

    std::vector<METOOLS::Current*> m_anticur;
    METOOLS::Current *m_antifcur;
    METOOLS::Vertex *m_antiv1;

    std::vector<size_t> m_nhel;

    size_t NHel(const ATOOLS::Flavour& fl);
    METOOLS::Vertex* GetVertex(METOOLS::Current* cur1,
                                  METOOLS::Current* cur2,
                                  METOOLS::Current* prop);      
  public:
    Comix1to2(const std::vector<ATOOLS::Flavour>& flavs);
    ~Comix1to2();
    void Calculate(const ATOOLS::Vec4D_Vector& momenta, bool anti=false);
  }; 

}

