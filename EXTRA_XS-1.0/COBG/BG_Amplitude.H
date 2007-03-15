#ifndef BG_Amplitude_H
#define BG_Amplitude_H

#include "C_Vector.H"
#include "Spinor.H"
#include "Flavour.H"

namespace EXTRAXS {

  typedef std::vector<int>        Int_Vector;
  typedef std::vector<Int_Vector> Int_Matrix;
  typedef std::vector<Complex>    Complex_Vector;

  typedef std::vector<ATOOLS::Vec4D>   Vec4D_Vector;
  typedef std::vector<ATOOLS::Flavour> Flavour_Vector;

  class BG_Amplitude {
  private:

    ATOOLS::Vec4D  *p_p, **p_ps, *p_ssp;
    ATOOLS::CVec4D *p_ep, *p_em, *p_ssep, *p_ssem, ***p_j;

    int *p_ch, *p_id;

    ATOOLS::Spinor m_kp, m_km;

    size_t m_n, m_nh, m_maxid;

    Int_Matrix     m_chirs, m_schirs;
    Complex_Vector m_ress;
    Int_Vector     m_hmap, m_hamps;

    ATOOLS::CVec4D EP(const ATOOLS::Vec4D &p);
    ATOOLS::CVec4D EM(const ATOOLS::Vec4D &p);

    ATOOLS::CVec4D V3L(const ATOOLS::Vec4D &p1,const ATOOLS::Vec4D &p2,
		       const ATOOLS::CVec4D &j1,const ATOOLS::CVec4D &j2);
    ATOOLS::CVec4D V4L(const ATOOLS::CVec4D &j1,const ATOOLS::CVec4D &j2,
		       const ATOOLS::CVec4D &j3);
  
    ATOOLS::CVec4D VL(const size_t &i,const size_t &j,const size_t &ih);
    ATOOLS::CVec4D JL(const size_t &i,const size_t &j,const size_t &ih);
    
    bool CheckChirs(const Int_Vector &chirs);
    bool MapChirs(const Int_Vector &chirs);
    bool Construct(Int_Vector chirs,const size_t &i);

    void CalcJL();

    void CleanUp();
    bool Construct(const size_t &n);

  public:

    BG_Amplitude();

    ~BG_Amplitude();

    bool EvaluateAll(const Int_Vector &perm);

    bool Construct(Flavour_Vector flavs,Int_Vector chirs);
    bool GaugeTest(const Vec4D_Vector &moms);

    void SetMomenta(const Vec4D_Vector &moms);

    size_t     MakeId(const Int_Vector &ids) const;
    Int_Vector MakeId(const size_t &id) const;

    inline const Int_Matrix &Chiralities() const { return m_chirs; }
    inline const Complex_Vector &Results() const { return m_ress; }

  };// end of class BG_Amplitude

}// end of namespace EXTRAXS

#endif
