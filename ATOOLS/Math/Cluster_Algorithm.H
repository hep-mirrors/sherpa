#ifndef ATOOLS_Math_Cluster_Algorithm_H
#define ATOOLS_Math_Cluster_Algorithm_H

#include <vector>

namespace ATOOLS {

  struct cs {

    enum code {
      num  = 0,
      dist = 1
    };

  };// end of struct cs

  template <class PointType,class MeasureType,class RecombinationType>
  class Cluster_Algorithm {
  public:

    typedef PointType         Point_Type;
    typedef MeasureType       Measure_Type;
    typedef RecombinationType Recombination_Type;

    typedef std::vector<Point_Type> Point_Vector;
    typedef std::vector<size_t>     Size_Vector;
    typedef std::vector<int>        Int_Vector;

    typedef std::vector<double>        Double_Vector;
    typedef std::vector<Double_Vector> Double_Matrix;

  private:

    Point_Vector m_p, m_lp, m_sp;
    Size_Vector  m_i;
    Int_Vector   m_j;
    size_t       m_n, m_imin, m_jmin;

    Double_Matrix m_d;
    Double_Vector m_r;
    double        m_dmin;

    Measure_Type       m_measure;
    Recombination_Type m_recom;

    bool m_recalc;

    inline void SetDMin(const size_t &i, const size_t &j,const double &dij) 
    {
      if (dij<m_dmin) {
	m_dmin=dij;
	m_imin=i;
	m_jmin=j;
      }
    }

    bool ArrangePoints();

  public:

    // constructor
    Cluster_Algorithm();

    // member functions
    bool Cluster(const double &crit,const cs::code &code=cs::dist);

    // inline functions
    inline void SetPoints(const Point_Vector &p)   { m_sp=m_p=p; }
    inline void SetClusters(const Point_Vector &p) { m_p=p;      }

    inline void SetRecalc(const bool &recalc) { m_recalc=recalc; }

    inline bool Recalc() const { return m_recalc; }

    inline Measure_Type       &Measure()       { return m_measure; }
    inline Recombination_Type &Recombination() { return m_recom;   }

    inline const Point_Vector  &Points() const    { return m_sp; }
    inline const Point_Vector  &Clusters() const  { return m_p;  }
    inline const Point_Vector  &Losts() const     { return m_lp; }
    inline const Int_Vector    &Positions() const { return m_j;  }
    inline const Double_Vector &DMins() const     { return m_r;  }

    inline size_t NClusters() const { return m_n;    }
    inline double DMin() const      { return m_dmin; }

    inline const Size_Vector &CPos() const { return m_i; }

  };// end of class Cluster

}// end of namespace ATOOLS

#endif
