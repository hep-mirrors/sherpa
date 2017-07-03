#ifndef PHASIC_Process_Subprocess_Info_H
#define PHASIC_Process_Subprocess_Info_H

#include "ATOOLS/Phys/Decay_Info.H"
#include "ATOOLS/Phys/NLO_Types.H"

namespace PHASIC {

  std::string PSId(const size_t &id);

  struct Subprocess_Info {
  public:

    ATOOLS::Flavour m_fl;
    std::string     m_id, m_pol, m_sv, m_mpl;

    std::vector<Subprocess_Info> m_ps;
    ATOOLS::DecayInfo_Vector m_decins;
    
    size_t m_nmax, m_nmin;
    int    m_tag, m_osf;

    ATOOLS::nlo_type::code m_nlotype;
    std::vector<double> m_nlocpl;

    void GetExternal(std::vector<ATOOLS::Flavour> &fl) const;
    void SetExternal(const std::vector<ATOOLS::Flavour> &fl,size_t &n);

    bool GetExternal(ATOOLS::Flavour &fl,const size_t &i,size_t &n) const;
    bool SetExternal(const ATOOLS::Flavour &fl,const size_t &i,size_t &n);

    size_t GetDecayInfos(ATOOLS::DecayInfo_Vector &ids,size_t &n,bool init);
    void DeleteDecayInfos();

    void SetTags(const std::vector<int>& tags,int &n);

  public:

    Subprocess_Info(const ATOOLS::Flavour &fl=ATOOLS::Flavour(kf_none),
		    const std::string &id="",const std::string &pol="",
		    const std::string &mpl="");
    ~Subprocess_Info();

    std::string MultiplicityTag() const;

    double Factorial(const double &n) const;

    double ISSymmetryFactor() const;
    double FSSymmetryFactor() const;

    size_t NExternal() const;
    size_t NTotalExternal() const;
    size_t NMinExternal() const;
    size_t NMaxExternal() const;
    bool   IsGroup() const;

    void SetExternal(const std::vector<ATOOLS::Flavour> &fl);
    void SetExternal(const ATOOLS::Flavour &fl,const size_t &i);

    std::vector<ATOOLS::Flavour> GetExternal() const;
    ATOOLS::Flavour              GetExternal(const size_t &i) const;

    void Add(const Subprocess_Info &info);
    bool AddDecay(const Subprocess_Info &ii,const Subprocess_Info &fi,int osf=0);

    inline ATOOLS::DecayInfo_Vector GetDecayInfos() const { return m_decins; }
    void BuildDecayInfos(size_t nin=2);

    void Print(std::ostream &ostr,const size_t &ni=0) const;

    void SetNMax(const Subprocess_Info &ref);
    void GetNMax(const Subprocess_Info &ref);

    inline ATOOLS::nlo_type::code NLOType() const { return m_nlotype; }
    inline void SetNLOType(ATOOLS::nlo_type::code nt) { m_nlotype=nt; }

    void SetTags(int& start);
    void SetTags(const std::vector<int>& tags);
    void GetTags(std::vector<int>& tags) const;

    int Combine(const size_t &i,const size_t &j,
		const ATOOLS::Flavour &flij,int &cnt);

    void ExtractMPL(std::vector<ATOOLS::Flavour_Vector> &fl) const;

    bool operator<(const Subprocess_Info &pi) const;
    bool operator==(const Subprocess_Info &pi) const;

  };// end of struct Subprocess_Info

  std::ostream &operator<<(std::ostream &ostr,const Subprocess_Info &info);

}// end of namespace PHASIC

#endif

