#ifndef AMISIC_Perturbative_MI_Integrator_H
#define AMISIC_Perturbative_MI_Integrator_H

namespace PDF {
  class ISR_Handler;
}

namespace AMISIC {
  class MI_Processes;
  class Matter_Overlap;
  class axis;
}

namespace AMISIC {
  class MI_Integrator {
  private:
    MI_Processes * p_procs;
    double         m_pt2min, m_xmin[2], m_xmax[2], m_Emin;
    unsigned int   m_MCpoints;
    double         m_xt, m_ymax, m_yvol, m_y[2], m_x[2];
    double         m_pt2, m_shat, m_that, m_uhat, m_cost, m_b;
    double         m_xsec, m_uncert, m_xsmax;
  public:
    MI_Integrator(MI_Processes * procs);
    void   Initialize(PDF::ISR_Handler * isr);
    double operator()(const double & s,Matter_Overlap * mo,const double & b);
    bool   TrialEvent(const double & s,Matter_Overlap * mo);
    bool   MakeKinematics(const double & pt2, const double & s);
    
    inline const double & Integral()    const { return m_xsec; } 
    inline const double & Uncertainty() const { return m_uncert; } 

    inline const double & SHat() const { return m_shat; }
    inline const double & THat() const { return m_that; }
    inline const double & UHat() const { return m_uhat; }
    inline const double & Yvol() const { return m_yvol; }
    inline const double & PT2()  const { return m_pt2;  }
    inline const double & B()    const { return m_b;    }
    inline const double & X(unsigned short int i) const { return m_x[i]; }
    inline const double & Y(unsigned short int i) const { return m_y[i]; }
  };
}
#endif
