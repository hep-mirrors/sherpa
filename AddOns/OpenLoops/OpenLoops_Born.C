#include "AddOns/OpenLoops/OpenLoops_Born.H"

#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Library_Loader.H"

using namespace PHASIC;
using namespace MODEL;
using namespace ATOOLS;
using namespace std;

namespace OpenLoops {

OpenLoops_Interface* OpenLoops_Born::s_interface=NULL;

OpenLoops_Born::OpenLoops_Born(const Process_Info& pi,
                               const Flavour_Vector& flavs,
                               Amp2Func amp2,
                               PermutationFunc permutationfunc,
                               std::vector<int> permutation,
                               std::string functag) :
  Tree_ME2_Base(pi, flavs),
  m_amp2(amp2), m_permutationfunc(permutationfunc),
  m_permutation(permutation)
{
  m_symfac=pi.m_fi.FSSymmetryFactor();
  m_symfac*=pi.m_ii.ISSymmetryFactor();
  m_oew=pi.m_oew;
  m_oqcd=pi.m_oqcd;
}

OpenLoops_Born::~OpenLoops_Born()
{
}

double OpenLoops_Born::Calc(const Vec4D_Vector& momenta)
{
  Vec4D_Vector m_moms(momenta);

  s_interface->SetParameter("alpha", AlphaQED());
  s_interface->SetParameter("alphas", AlphaQCD());

  double M2L0;
  std::vector<double> M2L1(3), M2L2(5), IRL1(3), IRL2(5);

  m_permutationfunc(&m_permutation[0]);
  m_amp2(&m_moms[0][0], &M2L0, &M2L1[0], &IRL1[0], &M2L2[0], &IRL2[0]);

  if (IsZero(M2L1[1]) && IsZero(M2L1[2]) &&
      IsZero(IRL1[1]) && IsZero(IRL1[2]) &&
      IsZero(M2L2[1]) && IsZero(M2L2[2]) && IsZero(M2L2[3]) && IsZero(M2L2[4]) &&
      IsZero(IRL2[1]) && IsZero(IRL2[2]) && IsZero(IRL2[3]) && IsZero(IRL2[4])) {
    // OL returns ME2 including 1/symfac, but Calc is supposed to return it
    // without 1/symfac, thus multiplying with symfac here
    if (IsZero(M2L0)) return m_symfac*M2L2[0];
    else return m_symfac*M2L0;
  }
  else {
    PRINT_INFO("Poles non-zero. Returning 0.");
    return 0.0;
  }
}

}

using namespace OpenLoops;

DECLARE_TREEME2_GETTER(OpenLoops_Born,"OpenLoops_Born")
Tree_ME2_Base *ATOOLS::Getter<Tree_ME2_Base,Process_Info,OpenLoops_Born>::
operator()(const Process_Info &pi) const
{
  DEBUG_FUNC(pi);
  if (pi.m_loopgenerator!="OpenLoops") return NULL;
  if (pi.m_fi.m_nloewtype!=nlo_type::lo) return NULL;
  if (pi.m_fi.m_nloqcdtype!=nlo_type::lo &&
      pi.m_fi.m_nloqcdtype!=nlo_type::born &&
      pi.m_fi.m_nloqcdtype!=nlo_type::real) return NULL;
  if (MODEL::s_model->Name()!="SM") return NULL;

  Flavour_Vector flavs=pi.ExtractFlavours();
  Flavour_Vector map_flavs=OpenLoops_Interface::MapFlavours(flavs);
  msg_Tracking()<<endl<<flavs<<" --> "<<map_flavs<<endl;

  vector<int> permutation;
  string process=OpenLoops_Interface::GetProcessPermutation(map_flavs, permutation);
  pair<string, string> groupsub=OpenLoops_Interface::ScanFiles(process, pi.m_oew, pi.m_oqcd, 0);
  string grouptag=groupsub.first;
  string subid=groupsub.second;
  if (grouptag!="") {
    // symbols in fortran are always defined as lower case
    string lc_functag(grouptag+"_"+process+"_"+subid+"_");
    for (size_t i(0);i<lc_functag.length();++i)
      lc_functag[i]=tolower(lc_functag[i]);
    vector<string> suffixes;
    suffixes.push_back("1sL"); // deprecated
    suffixes.push_back("ls");
    suffixes.push_back("lst");
    suffixes.push_back("lps");
    suffixes.push_back("lpst");
    void *ampfunc, *permfunc;
    for (size_t i=0; i<suffixes.size(); ++i) {
      string libraryfile="openloops_"+grouptag+"_"+suffixes[i];
      ampfunc=s_loader->GetLibraryFunction(libraryfile,"vamp2_"+lc_functag);
      permfunc=s_loader->GetLibraryFunction(libraryfile,"set_permutation_"+lc_functag);
      if (ampfunc!=NULL && permfunc!=NULL) break;
    }
    if (ampfunc==NULL || permfunc==NULL) {
      PRINT_INFO("Didn't find functions");
      return NULL;
    }

    msg_Info()<<endl;
    PRINT_INFO("Initialising OpenLoops Born for "<<flavs<<": "<<lc_functag);
    return new OpenLoops_Born(pi, flavs, (Amp2Func) ampfunc,
                              (PermutationFunc) permfunc, permutation, lc_functag);
  }
  else {
    return NULL;
  }
}


