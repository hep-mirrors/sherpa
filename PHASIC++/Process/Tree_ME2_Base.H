#ifndef PHASIC_Process_Tree_ME2_Base_H
#define PHASIC_Process_Tree_ME2_Base_H

#include "PHASIC++/Process/Process_Info.H"
#include "ATOOLS/Math/Vector.H"
#include "MODEL/Main/Model_Base.H"

namespace MODEL { struct Coupling_Map; }

namespace PHASIC {

  class External_ME_Args;

  class Tree_ME2_Base {
  public:

    struct Map_Info {
      std::vector<int> m_perm, m_hels;
      size_t m_id;
    };

  protected:

    const ATOOLS::Flavour_Vector m_flavs;
    MODEL::Coupling_Data* p_aqcd, * p_aqed;

    size_t m_namps;
    double m_norm;

    void CompleteCombinations
    (std::set<std::pair<size_t,size_t> > &combs,
     std::map<size_t,ATOOLS::Flavour_Vector> &fls) const;

  public:

    Tree_ME2_Base(const External_ME_Args& args);

    virtual ~Tree_ME2_Base();

    virtual double Calc(const ATOOLS::Vec4D_Vector &p) = 0;

    virtual void SetCouplings(const MODEL::Coupling_Map& cpls);
    virtual double AlphaQCD() const;
    virtual double AlphaQED() const;

    virtual std::vector<Complex> GetAmplitudes(const size_t &id);
    virtual Complex GetPhase(const size_t &id);

    virtual Complex GetHelicityPhase(const ATOOLS::Vec4D &pijt,
				     const ATOOLS::Vec4D &eps1);

    virtual std::vector<Map_Info> GetFlavourHelicityMap();

    virtual void FillCombinations
    (std::set<std::pair<size_t,size_t> > &combs,
     std::map<size_t,ATOOLS::Flavour_Vector> &fls);

    virtual int OrderQCD(const int &id=-1) const;
    virtual int OrderEW(const int &id=-1) const;

    virtual double TR() const;

    static Tree_ME2_Base *GetME2(const PHASIC::External_ME_Args& pi);
    static Tree_ME2_Base *GetME2(const std::string& tag,
				 const PHASIC::External_ME_Args& pi);

    // Legacy method, implemented as wrapper around new versions
    static Tree_ME2_Base *GetME2(const PHASIC::Process_Info& pi);

    inline size_t NAmps() const { return m_namps; }

    inline void SetNorm(const double &norm) { m_norm=norm; }

  };

  std::ostream &operator<<(std::ostream &str,
			   const Tree_ME2_Base::Map_Info &mi);

  class Trivial_Tree : public Tree_ME2_Base {
  public:

    Trivial_Tree(const PHASIC::External_ME_Args& args) :
      Tree_ME2_Base(args) {}

    double Calc(const ATOOLS::Vec4D_Vector &p);

  };
}

#define DECLARE_TREEME2_GETTER(NAME,TAG)			   \
  DECLARE_GETTER(NAME,TAG,PHASIC::Tree_ME2_Base,PHASIC::External_ME_Args);	   \
  void ATOOLS::Getter<PHASIC::Tree_ME2_Base,PHASIC::External_ME_Args,NAME>::	   \
  PrintInfo(std::ostream &str,const std::size_t width) const		   \
  {                                                                \
    str<<#TAG;                                                     \
  }

#endif
