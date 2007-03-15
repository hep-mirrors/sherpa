#ifndef CDBG_Amplitude_H
#define CDBG_Amplitude_H

#include "Spinor.H"
#include "C_Vector.H"
#include "C_Spinor.H"
#include "Flavour.H"
#include "Vertex.H"

namespace EXTRAXS {

  typedef std::vector<Complex> Complex_Vector;

  typedef std::vector<Current*>       Current_Vector;
  typedef std::vector<Current_Vector> Current_Matrix;

  typedef std::vector<int>        Int_Vector;
  typedef std::vector<Int_Vector> Int_Matrix;

  class CDBG_Amplitude {
  private:

    ATOOLS::Vec4D   *p_p, *p_ssp;
    ATOOLS::Flavour *p_fl;
    ATOOLS::CSpinor *p_sp, *p_sm;
    ATOOLS::CVec4D  *p_ep, *p_em, *p_ssep, *p_ssem;

    int *p_ch, **p_cl;

    ATOOLS::Spinor m_kp, m_km;

    size_t m_n, m_maxid;
    int    m_mode, m_nf;

    Int_Vector     m_hmap, m_hamps;
    Int_Matrix     m_chirs;
    Complex_Vector m_ress;
    Current_Matrix m_cur;

    ATOOLS::CVec4D m_jv;

    ATOOLS::CVec4D EP(const ATOOLS::Vec4D &p);
    ATOOLS::CVec4D EM(const ATOOLS::Vec4D &p);

    void CalcJL(const size_t &n,const size_t &id);

    bool CheckChirs(const Int_Vector &chirs);
    bool MapChirs(const Int_Vector &chirs);
    bool ConstructChirs(Int_Vector chirs,const size_t &i);

    void CalcJL();

    void CleanUp();
    void Construct(Int_Vector ids,const size_t &n);
    bool Construct(const std::vector<ATOOLS::Flavour> &flavs);

    bool MatchIndices(const Int_Vector &ids,const size_t &n,
		      const size_t &i,const size_t &j,const size_t &k);

    void AddCurrent(const Int_Vector &ids,const size_t &n,
		    const ATOOLS::Flavour &fl);

  public:

    CDBG_Amplitude();

    ~CDBG_Amplitude();

    Complex Evaluate(const Int_Vector &chirs);
    bool    EvaluateAll();

    bool    EvaluateAll(const Int_Vector &perm);

    bool Construct(std::vector<ATOOLS::Flavour> flavs,
		   Int_Vector chirs);
    bool GaugeTest(const std::vector<ATOOLS::Vec4D> &moms);

    void SetMomenta(const std::vector<ATOOLS::Vec4D> &moms);
    void SetColors(const Int_Vector &rc,
		   const Int_Vector &ac);

    bool DWITest();
    bool CyclicityTest();
    bool ReflectionTest();

    size_t     MakeId(const Int_Vector &ids) const;
    Int_Vector MakeId(const size_t &id) const;

    inline void SetMode(const int mode) { m_mode=mode; }

    inline const ATOOLS::CVec4D &JV() const { return m_jv; }

    inline int Mode() const { return m_mode; }

    inline const Int_Matrix &Chiralities() const { return m_chirs; }
    inline const Complex_Vector &Results() const { return m_ress;  }

  };// end of class CDBG_Amplitude

}// end of namespace EXTRAXS

#endif
