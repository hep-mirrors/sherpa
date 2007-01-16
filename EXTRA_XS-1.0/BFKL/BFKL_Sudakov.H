#ifndef BFKL_Sudakov_H
#define BFKL_Sudakov_H

#include "LEV_Group.H"

namespace EXTRAXS {

  class BFKL_Sudakov {
  public:

    typedef std::vector<ATOOLS::Vec4D> Vector_Vector;

  private:

    LEV_Group *p_levs;

    double m_ya, m_yb, m_kt2min, m_kt2max;
    double m_asmax, m_integral, m_gamma, m_weight;
    double m_y, m_kt2, m_phi;

    ATOOLS::Vec4D   m_q;
    ATOOLS::Flavour m_flq;

    int m_ktscheme, m_splitmode;

    static size_t s_nf;

  public:

    // constructor
    BFKL_Sudakov();

    // destructor
    ~BFKL_Sudakov();

    // member functions
    bool Initialize();

    bool Init();
    bool Dice();
    bool Approve(const ATOOLS::Vec4D &k1,const ATOOLS::Vec4D &q1,
		 const ATOOLS::Vec4D &k2,const ATOOLS::Vec4D &q2);

    // inline functions
    inline void SetYA(const double &ya) { m_ya=ya; }
    inline void SetYB(const double &yb) { m_yb=yb; }

    inline void SetKT2Min(const double &kt2min) { m_kt2min=kt2min; }
    inline void SetKT2Max(const double &kt2max) { m_kt2max=kt2max; }

    inline void SetInMomentum(const ATOOLS::Vec4D &q)  { m_q=q;   }
    inline void SetInFlavour(const ATOOLS::Flavour &f) { m_flq=f; }

    inline void SetKTScheme(const int &scheme) { m_ktscheme=scheme; }
    inline void SetSplitMode(const int &mode)  { m_splitmode=mode;  }

    inline double GetY() const     { return m_y;    }
    inline double GetKT2() const   { return m_kt2;  }
    inline double GetPhi() const   { return m_phi;  }

    inline double GetWeight() const { return m_weight;  }

    inline LEV_Base *Selected() const { return p_levs->Selected(); }

    inline int KTScheme() const  { return m_ktscheme;  }
    inline int SplitMode() const { return m_splitmode; }

    inline size_t Nf() const { return s_nf; }

  };// end of class BFKL_Sudakov

}// end of namespace EXTRAXS

#endif
