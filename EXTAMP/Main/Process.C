#include "MODEL/Main/Model_Base.H"
#include "EXTAMP/Main/External_ME_Interface.H"
#include "EXTAMP/Main/Process.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"

#include <assert.h>

namespace EXTAMP {

  Process::Process(const PHASIC::Process_Info& pi) {

    /* Also is done in the base class method
       PHASIC::Process_Base::Init but we don't call that method
       because it takes beam information and is thus only applicable
       to more high-level processes */
    m_flavs = pi.ExtractFlavours();
    m_nin   = pi.m_ii.NExternal();
    m_nout  = pi.m_fi.NExternal();
    m_pinfo = pi;

    /* Not done in PHASIC base classes */
    m_mincpl = pi.m_mincpl;
    m_maxcpl = pi.m_maxcpl;

    /* Factor due to spin and color averaging in the initial state,
       i.e. number of unobserved degrees of freedom in initial state
       (polarizations, colors) */
    m_norm =pi.m_fi.FSSymmetryFactor();

    /* Count the number of occurrences of each flavour n_i in the
       final state and calculate \prod_i factorial(n_i) */
    m_norm*=pi.m_ii.ISSymmetryFactor();

    /* m_cpls is member of PHASIC::Process_Base, but does not get
       populated for some reason in PHASIC::Process_Base. Instead,
       it is initialised in the various 'Single_Process' types as
       done below. */
    MODEL::s_model->GetCouplings(m_cpls);

    FillPartonIndices();

    m_cluster_flav_map
      = External_ME_Interface::ConstructCombinableMap(Flavours(),
						      Info(),
						      NIn());
  }


  void Process::FillPartonIndices()
  {
    m_parton_indices.clear();
    const ATOOLS::Flavour& jet(kf_jet);
    for(size_t i(0); i<Flavours().size(); i++)
      if(jet.Includes(Flavours()[i]))
	m_parton_indices.push_back(i);
  }


  bool Process::Combinable(const size_t &idi,const size_t &idj)
  {
    bool ret = m_cluster_flav_map.find(idi | idj)!=m_cluster_flav_map.end();
    return ret;
  }


  const ATOOLS::Flavour_Vector &Process::CombinedFlavour(const size_t &idij)
  {
    std::map<size_t, ATOOLS::Flavour_Vector>::const_iterator it = m_cluster_flav_map.find(idij);
    if(it==m_cluster_flav_map.end())
      THROW(fatal_error, "Internal error");
    return it->second;
  }

  int Process::PerformTests()
  {
    ATOOLS::Vec4D_Vector p(m_flavs.size(), ATOOLS::Vec4D());
    PHASIC::Phase_Space_Handler::TestPoint(&p.front(),&Info(),Generator(),1);
    return 1;
  }

}
