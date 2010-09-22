#ifndef PHASIC_Main_RB_Data_H
#define PHASIC_Main_RB_Data_H

#include "ATOOLS/Phys/Flavour.H"

namespace ATOOLS {

  class NLO_subevt;
  class Histogram;

  struct RB_Key {
    int m_type;
    ATOOLS::Flavour m_flij, m_fli, m_flj;
  public:
    RB_Key(const int type,const ATOOLS::Flavour &flij,
           const ATOOLS::Flavour &fli,const ATOOLS::Flavour &flj);
    RB_Key(const NLO_subevt *sub);
    bool operator<(const RB_Key &rb) const;
    inline RB_Key SwapIJ() const
    { return RB_Key(m_type,m_flij,m_flj,m_fli); }
  };// end of struct RB_Key

  std::ostream &operator<<(std::ostream &str,const RB_Key &rbk);
  std::istream &operator>>(std::istream &str,RB_Key &rbk);

  struct RB_Data {
    double m_sum, m_max, m_n, m_rbint, m_xssum;
    double m_ktres, m_bmax;
    Histogram *p_wh;
  public:
    RB_Data();
    ~RB_Data();
    void AddPoint(const double &rb,const double &xs,const double &b);
    void InitWeightHisto();
    void ReadInHisto(const std::string &fname);
    void WriteOutHisto(const std::string &fname) const;
    double GetMaxEps(const double &eps);
    RB_Data &operator=(const RB_Data &rb);
    inline double Mean() const { return m_rbint/m_xssum; }
  };// end of struct RB_Data

  std::ostream &operator<<(std::ostream &str,const RB_Data &rbd);
  std::istream &operator>>(std::istream &str,RB_Data &rbd);

  typedef std::map<RB_Key,RB_Data*> RB_Map;

}// end of namespace ATOOLS

#endif
