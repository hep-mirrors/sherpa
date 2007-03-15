#ifndef Helicity_Integrator_H
#define Helicity_Integrator_H

#include <vector>

namespace PHASIC {

  typedef std::vector<int>    Int_Vector;
  typedef std::vector<double> Double_Vector;

  class Helicity_Integrator {
  private:

    Int_Vector    m_chirs, m_n;
    Double_Vector m_weights, m_asum, m_sum, m_sum2;

    size_t m_valid, m_id, m_iter;
    double m_weight;
    bool   m_new;

    void Construct(Int_Vector chirs,const size_t i);

  public:

    // constructor
    Helicity_Integrator();

    // destructor
    ~Helicity_Integrator();

    // member functions
    bool Construct(const size_t &n);

    bool   GeneratePoint();
    double Weight();
    
    void AddPoint(const double &weight);

    void Optimize();
    void WriteOut(const std::string &pid);
    void ReadIn(const std::string &pid);

    size_t     MakeId(const Int_Vector &ids) const;
    Int_Vector MakeId(const size_t &id) const;

    // inline functions
    inline const std::vector<int> &Chiralities() const 
    { return m_chirs; }

  };// end of class Helicity_Integrator

}// end of namespace PHASIC

#endif
