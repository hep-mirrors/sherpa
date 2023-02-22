#ifndef PHASIC_Process_Process_Info_H
#define PHASIC_Process_Process_Info_H

#include "ATOOLS/Phys/NLO_Types.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "PHASIC++/Process/Subprocess_Info.H"
#include "PHASIC++/Main/Color_Integrator.H"
#include "PHASIC++/Main/Helicity_Integrator.H"

namespace PHASIC {
  
  class ME_Generators;

  struct Process_Info {
  public:

    Subprocess_Info m_ii, m_fi;

    cls::scheme m_cls;
    hls::scheme m_hls;

    std::vector<double> m_maxcpl, m_mincpl, m_maxacpl, m_minacpl;

    size_t m_nminq, m_nmaxq, m_ckkw, m_itmin, m_rsitmin;
    size_t m_amegicmhv, m_ntchan, m_mtchan, m_sort;

    ATOOLS::nlo_mode::code m_nlomode;

    std::string m_scale, m_coupling, m_kfactor, m_gpath, m_addname;
    std::string m_megenerator, m_rsmegenerator, m_loopgenerator, m_special;
    std::string m_integrator, m_rsintegrator;

    ATOOLS::Scoped_Settings m_selectors;

    ME_Generators* p_gens;

    ATOOLS::Flavour_Vector m_nodecs;

  public:

    inline Process_Info(const Subprocess_Info &ii=Subprocess_Info(),
			const Subprocess_Info &fi=Subprocess_Info()):
      m_ii(ii), m_fi(fi), m_cls(cls::unknown), m_hls(hls::sum),
      m_maxcpl(2,99), m_mincpl(2,0), m_maxacpl(2,99), m_minacpl(2,0),
      m_nminq(0), m_nmaxq(99), m_ckkw(0), m_itmin(1000), m_rsitmin(1000),
      m_amegicmhv(1), m_ntchan(0), m_mtchan(99), m_sort(3),
      m_nlomode(ATOOLS::nlo_mode::none),
      m_scale("VAR"), m_kfactor("None"), m_loopgenerator("Internal"),
      p_gens(NULL) {}

    ATOOLS::Flavour_Vector ExtractFlavours() const;
    bool Has(ATOOLS::nlo_type::code nlotype) const;

    int Combine(const size_t &i,const size_t &j,
		const ATOOLS::Flavour &flij);

    std::vector<ATOOLS::Flavour_Vector> ExtractMPL() const;

    bool operator<(const Process_Info &pi) const;

  };// end of struct Process_Info

  std::ostream &operator<<(std::ostream &ostr,const Process_Info &info);

}// end of namespace PHASIC

#endif

