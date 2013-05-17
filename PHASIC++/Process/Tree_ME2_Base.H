#ifndef PHASIC_Process_Tree_ME2_Base_H
#define PHASIC_Process_Tree_ME2_Base_H

#include "PHASIC++/Process/Process_Info.H"
#include "ATOOLS/Math/Vector.H"
#include "MODEL/Main/Model_Base.H"

namespace MODEL { struct Coupling_Map; }

namespace PHASIC {

  class Tree_ME2_Base {
  public:

    struct Map_Info {
      std::vector<int> m_perm, m_hels;
      size_t m_id;
    };

  protected:

    const Process_Info &m_pinfo;
    const ATOOLS::Flavour_Vector m_flavs;
    MODEL::Coupling_Data* p_aqcd, * p_aqed;

  public:

    Tree_ME2_Base(const PHASIC::Process_Info &pi,
                  const ATOOLS::Flavour_Vector &flavs);

    virtual ~Tree_ME2_Base();

    virtual double Calc(const ATOOLS::Vec4D_Vector &p) = 0;

    virtual void SetCouplings(const MODEL::Coupling_Map& cpls);
    virtual double AlphaQCD() const;
    virtual double AlphaQED() const;

    virtual std::vector<Complex> GetAmplitudes();

    virtual Complex GetHelicityPhase(const ATOOLS::Vec4D &pijt,
				     const ATOOLS::Vec4D &eps1);

    virtual std::vector<Map_Info> GetFlavourHelicityMap();

    virtual void FillCombinations
    (std::set<std::pair<size_t,size_t> > &combs,
     std::map<size_t,ATOOLS::Flavour_Vector> &fls);

    virtual int OrderQCD(const int &id=-1);
    virtual int OrderEW(const int &id=-1);

    virtual double TR() const;

    static Tree_ME2_Base *GetME2(const PHASIC::Process_Info& pi);
    static Tree_ME2_Base *GetME2(const std::string& tag,
                                    const PHASIC::Process_Info& pi);
  };

  std::ostream &operator<<(std::ostream &str,
			   const Tree_ME2_Base::Map_Info &mi);

  class Trivial_Tree : public Tree_ME2_Base {
  public:

    Trivial_Tree(const PHASIC::Process_Info &pi,
                 const ATOOLS::Flavour_Vector &flavs) :
      Tree_ME2_Base(pi,flavs) {}

    double Calc(const ATOOLS::Vec4D_Vector &p);

  };
}

#define DECLARE_TREEME2_GETTER(NAME,TAG)			   \
  DECLARE_GETTER(NAME,TAG,Tree_ME2_Base,Process_Info);		   \
  void ATOOLS::Getter<Tree_ME2_Base,Process_Info,NAME>::	   \
  PrintInfo(std::ostream &str,const size_t width) const		   \
  {                                                                \
    str<<#TAG;                                                     \
  }

#endif
