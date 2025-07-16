#ifndef PHASIC_Main_Helicity_Integrator_H
#define PHASIC_Main_Helicity_Integrator_H

#include "ATOOLS/Phys/Flavour.H"
#include <vector>

namespace PHASIC {

  struct hls {

    enum scheme {
      unknown   = 0,
      sum       = 1,
      sample    = 2
    };

    static std::map<std::string,std::string> HelicitySchemeTags();

  };// end of struct hls

  std::ostream &operator<<(std::ostream &str,const hls::scheme &s);

  typedef std::vector<int>    Int_Vector;
  typedef std::vector<double> Double_Vector;

  class Helicity_Integrator {
  private:

    Int_Vector     m_chirs, m_n;
    Double_Vector  m_weights, m_asum, m_sum, m_sum2;

    ATOOLS::Flavour_Vector m_flavs;

    size_t m_valid, m_id, m_iter;
    double m_weight;
    bool   m_new, m_on;

    bool CheckChirs(const Int_Vector &chirs);
    void Construct(Int_Vector chirs,const size_t i);

  public:

    // constructor
    Helicity_Integrator();

    // destructor
    ~Helicity_Integrator();

    // member functions
    bool Construct(const ATOOLS::Flavour_Vector &flavs);

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

    inline void SetOn(const bool on) { m_on=on; }

    inline bool On() const { return m_on; }

  };// end of class Helicity_Integrator

}// end of namespace PHASIC

#endif
