#ifndef PHASIC_Process_Subprocess_Info_H
#define PHASIC_Process_Subprocess_Info_H

#include "ATOOLS/Phys/Flavour.H"
#include "PHASIC++/Process/NLO_Helpers.H"

#include <vector>

namespace PHASIC {

  std::string PSId(const size_t &id);

  struct Decay_Info {
    long unsigned int m_id;
    ATOOLS::Flavour m_fl;
    size_t m_nmax;
    int m_osd;
    Decay_Info(const long unsigned int &id,
	       const ATOOLS::Flavour &fl,
	       const size_t &nmax,const int osd):
      m_id(id), m_fl(fl), m_nmax(nmax), m_osd(osd) {}
  };// end of struct Decay_Info

  typedef std::vector<Decay_Info> Decay_Info_Vector;

  struct Subprocess_Info {
  public:

    ATOOLS::Flavour m_fl;
    std::string     m_id, m_pol;

    std::vector<Subprocess_Info> m_ps;

    size_t m_nmax;
    int    m_tag,m_osf;

    nlo_type::code m_nloqcdtype, m_nloewtype;

    void GetExternal(std::vector<ATOOLS::Flavour> &fl) const;
    void SetExternal(const std::vector<ATOOLS::Flavour> &fl,size_t &n);

    bool GetExternal(ATOOLS::Flavour &fl,const size_t &i,size_t &n) const;
    bool SetExternal(const ATOOLS::Flavour &fl,const size_t &i,size_t &n);

    size_t GetDecayInfos(Decay_Info_Vector &ids,size_t &n) const;

  public:

    Subprocess_Info(const ATOOLS::Flavour &fl=ATOOLS::Flavour(kf_none),
		    const std::string &id="",const std::string &pol="");

    std::string MultiplicityTag() const;

    size_t NExternal() const;
    size_t NTotalExternal() const;
    size_t NMaxExternal() const;

    void SetExternal(const std::vector<ATOOLS::Flavour> &fl);
    void SetExternal(const ATOOLS::Flavour &fl,const size_t &i);

    std::vector<ATOOLS::Flavour> GetExternal() const;
    ATOOLS::Flavour              GetExternal(const size_t &i) const;

    void Add(const Subprocess_Info &info);
    bool AddDecay(const Subprocess_Info &ii,const Subprocess_Info &fi,int osf=0);

    Decay_Info_Vector GetDecayInfos() const;

    void Print(std::ostream &ostr,const size_t &ni=0) const;

    void SetNMax(const Subprocess_Info &ref);
    void GetNMax(const Subprocess_Info &ref);

    nlo_type::code NLOType() const;
    void SetNLOType(nlo_type::code nlotype);

  };// end of struct Subprocess_Info

  std::ostream &operator<<(std::ostream &ostr,const Subprocess_Info &info);

}// end of namespace PHASIC

#endif

