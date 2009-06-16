#include "EXTRA_XS/NLO/Virtual_ME2_Base.H"
#include "EXTRA_XS/NLO/Loop_ME_Base.H"
#include "EXTRA_XS/Main/ME_Base.H"
#include "ATOOLS/Org/Exception.H"

using namespace EXTRAXS;
using namespace PHASIC;
using namespace ATOOLS;
using namespace HELICITIES;

namespace EXTRAXS {
  class TreeLoop_Virtual : public Virtual_ME2_Base {
    ME_Base* p_tree;
    Loop_ME_Base* p_loop;
  public:
    TreeLoop_Virtual(const Process_Info& pi, const Flavour_Vector& flavs,
                     ME_Base* tree, Loop_ME_Base* loop) :
      Virtual_ME2_Base(pi, flavs), p_tree(tree), p_loop(loop)
    {
    }

    ~TreeLoop_Virtual() {
      if (p_tree) { delete p_tree; }
      if (p_loop) { delete p_loop; }
    }

    void Calc(const ATOOLS::Vec4D_Vector& mom);

    // eps scheme factor?
  };
}


void TreeLoop_Virtual::Calc(const Vec4D_Vector& mom) {
  p_tree->Calc(mom);
  p_loop->Calc(mom);
  
  if (p_tree->Result().size()!=p_loop->Result().size())
    THROW(fatal_error, "n_hel(p_tree) != n_hel(p_loop)");

  DivArrC sum(vector<Complex>(5, Complex(0.0, 0.0)));
  for (size_t ihel=0; ihel<p_tree->Result().size(); ++ihel) {
    sum+=conj(p_tree->Result()[ihel])*p_loop->Result()[ihel];
  }
  if (!IsZero(imag(sum))) THROW(fatal_error, "imag(sum tree*loop)!=0");
  m_res=2.0*real(sum);
}


DECLARE_VIRTUALME2_GETTER(TreeLoop_Virtual_Getter,"TreeLoop_Virtual")
Virtual_ME2_Base *TreeLoop_Virtual_Getter::operator()(const Process_Info &pi) const
{
  if (pi.m_loopgenerator!="Internal") return NULL;
  if ((pi.m_fi.m_nloqcdtype==nlo_type::loop || pi.m_fi.m_nloewtype==nlo_type::loop) &&
      !(pi.m_fi.m_nloqcdtype==nlo_type::loop && pi.m_fi.m_nloewtype==nlo_type::loop)) {
    Loop_ME_Base* loop_me=Loop_ME_Base::GetME(pi);
    if (!loop_me) return NULL;
    Process_Info tree_pi(pi);
    tree_pi.m_fi.m_nloqcdtype=nlo_type::lo;
    tree_pi.m_fi.m_nloewtype=nlo_type::lo;
    ME_Base* tree_me=ME_Base::GetME(tree_pi);
    if (!tree_me) return NULL;

    Flavour_Vector fl=pi.ExtractFlavours();
    return new TreeLoop_Virtual(pi, fl, tree_me, loop_me);
  }
  return NULL;
}
