#ifndef CDBG_Amplitude_H
#define CDBG_Amplitude_H

#include "Spinor.H"
#include "C_Vector.H"
#include "C_Spinor.H"
#include "Flavour.H"
#include "Vertex.H"

namespace EXTRAXS {

  typedef std::vector<Complex>     Complex_Vector;

  typedef std::vector<Current_Base*>  Current_Vector;
  typedef std::vector<Current_Vector> Current_Matrix;

  typedef std::vector<ATOOLS::Vec4D>   Vec4D_Vector;
  typedef std::vector<ATOOLS::Flavour> Flavour_Vector;

  class CDBG_Amplitude {
  private:

    String_Vector  m_model;
    Flavour_Vector m_fl;

    Vec4D_Vector m_p;
    Int_Vector   m_ch; 
    Int_Matrix   m_cl;

    size_t m_n, m_maxid, m_nf;

    Int_Vector     m_hmap, m_hamps;
    Int_Matrix     m_chirs;
    Complex_Vector m_ress;
    Current_Matrix m_cur;

    void CleanUp();

    bool MatchIndices(const Int_Vector &ids,const size_t &n,
		      const size_t &i,const size_t &j,const size_t &k);
    void AddCurrent(const Int_Vector &ids,const size_t &n,
		    const ATOOLS::Flavour &fl);

    void Construct(Int_Vector ids,const size_t &n);
    bool Construct(const Flavour_Vector &flavs);

    bool CheckChirs(const Int_Vector &chirs);
    bool MapChirs(const Int_Vector &chirs);
    bool ConstructChirs(Int_Vector chirs,const size_t &i);

    void SetGauge(const size_t &n);

    void CalcJL();

    void WriteOutGraph(std::ostream &str,Graph_Node *graph,
		       size_t &ng) const;

  public:

    CDBG_Amplitude();

    ~CDBG_Amplitude();

    Complex Evaluate(const Int_Vector &chirs);
    bool    EvaluateAll(const Int_Vector &perm);
    bool    EvaluateAll();

    bool Construct(const Flavour_Vector &flavs,
		   const String_Vector &models);
    bool GaugeTest(const Vec4D_Vector &moms);

    void SetMomenta(const Vec4D_Vector &moms);
    void SetColors(const Int_Vector &rc,
		   const Int_Vector &ac);

    void WriteOutGraphs(const std::string &file) const;

    inline const Int_Matrix &Chiralities() const { return m_chirs; }
    inline const Complex_Vector &Results() const { return m_ress;  }

    inline const Flavour_Vector &Flavours() const { return m_fl; }

  };// end of class CDBG_Amplitude

}// end of namespace EXTRAXS

#endif
