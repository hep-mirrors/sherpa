#ifndef PHASIC_Process_Process_Info_H
#define PHASIC_Process_Process_Info_H

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

    size_t m_oew, m_oqcd, m_maxoew, m_maxoqcd;
    size_t m_nminq, m_nmaxq, m_ckkw, m_nmax, m_itmin;
    size_t m_nlomode, m_amegicmhv, m_ntchan, m_mtchan;

    std::string m_scale, m_coupling, m_kfactor, m_gpath;
    std::string m_megenerator, m_rsmegenerator, m_loopgenerator;

    std::string m_selectorfile, m_integrator, m_rsintegrator;

    bool m_mpiprocess;

    ME_Generators* p_gens;

    ATOOLS::Flavour_Vector m_nodecs;

  public:

    inline Process_Info(const Subprocess_Info &ii=Subprocess_Info(),
			const Subprocess_Info &fi=Subprocess_Info()):
      m_ii(ii), m_fi(fi), m_cls(cls::sum), m_hls(hls::sum),
      m_oew(99), m_oqcd(99), m_maxoew(99), m_maxoqcd(99),
      m_nminq(0), m_nmaxq(99), m_ckkw(0), m_nmax(0), m_itmin(5000),
      m_nlomode(0), m_amegicmhv(1), m_ntchan(0), m_mtchan(99),
      m_scale("VAR"), m_kfactor("NO"), m_loopgenerator("Internal"),
      m_selectorfile(""), m_mpiprocess(false), p_gens(NULL) {}

    ATOOLS::Flavour_Vector ExtractFlavours() const;
    bool Has(nlo_type::code nlotype) const;

    int Combine(const size_t &i,const size_t &j,
		const ATOOLS::Flavour &flij);

  };// end of struct Process_Info

  std::ostream &operator<<(std::ostream &ostr,const Process_Info &info);

}// end of namespace PHASIC

#endif

