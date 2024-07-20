#ifndef YFS_MAIN_Weight_H
#define YFS_MAIN_Weight_H

#include "ATOOLS/Math/Vector.H"

namespace YFS {
  class Weight {
    protected:
      double            m_weight;
      double            m_maxweight;

      virtual void CalculateWeight() = 0;
      virtual void CalculateMax() = 0;

    public:
      Weight();
     ~Weight();
      double m_jacW, m_massW;
      // void JacobianWeight(ATOOLS::Vec4D_Vector k);
      void TotalWeight();
      inline double SetJacobianWeight(double jacW)  { m_jacW = jacW; }
      inline double SetMassWeight(double massW)  { m_massW = massW; }
      inline double GetMax()    { return m_maxweight; }
  };


}

#endif
