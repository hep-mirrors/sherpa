#ifndef Reggeon_Trajectory_H
#define Reggeon_Trajectory_H

namespace AMISIC {

  class Reggeon_Trajectory {
  private:
    
    double m_alpha_0, m_alpha_prime, m_s, m_s_0;
    
  public:
    
    // constructor
    Reggeon_Trajectory(const double alpha_0,const double alpha_prime);
    
    // destructor
    ~Reggeon_Trajectory();
    
    // member functions
    void Fit(const double t,const double dsigma);
    
    double DSigma(const double t) const;
    double GetT(const double tmin,const double tmax,const double ran) const;
    
    // inline functions
    inline double operator()(const double t=0.0) const { return DSigma(t); }

    inline void SetS(const double s) { m_s=s; }
    
  };// end of class Reggeon_Trajectory

}// end of namespace AMISIC

#endif
