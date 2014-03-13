#include "OpenLoops_Virtual.H"

#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Library_Loader.H"
#include "ATOOLS/Math/Poincare.H"


using namespace PHASIC;
using namespace ATOOLS;
using namespace std;

namespace OpenLoops {

OpenLoops_Interface* OpenLoops_Virtual::s_interface=NULL;

OpenLoops_Virtual::OpenLoops_Virtual(const Process_Info& pi,
                                     const Flavour_Vector& flavs,
                                     Amp2Func amp2,
                                     PermutationFunc permutationfunc,
                                     std::vector<int> permutation,
                                     std::string functag) :
  Virtual_ME2_Base(pi, flavs),
  m_amp2(amp2), m_permutationfunc(permutationfunc),
  m_permutation(permutation), m_analyse(true)
{
  m_oew=pi.m_oew;
  m_oqcd=pi.m_oqcd;
}

OpenLoops_Virtual::~OpenLoops_Virtual()
{
}


void OpenLoops_Virtual::Calc(const Vec4D_Vector& momenta) {
  Vec4D_Vector m_moms(momenta);

  s_interface->SetParameter("alpha", AlphaQED());
  s_interface->SetParameter("alphas", AlphaQCD());
  s_interface->SetParameter("mu", sqrt(m_mur2));
  s_interface->FlushParameters();
  
  double M2L0;
  vector<double> M2L1(3), M2L2(5), IRL1(3), IRL2(5);

  MyTiming* timing;
  if (msg_LevelIsDebugging()) {
    timing = new MyTiming();
    timing->Start();
  }
  m_permutationfunc(&m_permutation[0]);
  m_amp2(&m_moms[0][0], &M2L0, &M2L1[0], &IRL1[0], &M2L2[0], &IRL2[0]);
  if (msg_LevelIsDebugging()) {
    timing->Stop();
    PRINT_INFO(momenta[2][0]<<" "<<m_flavs<<" user="<<timing->UserTime()
               <<" real="<<timing->RealTime()<<" sys="<<timing->SystemTime());
  }

  double B(M2L0), V_finite(M2L1[0]), V_eps(M2L1[1]), V_eps2(M2L1[2]);

  // factor which by Sherpa convention has to be divided out at this stage
  double factor=B*AlphaQCD()/2.0/M_PI;

  m_born=B;
  m_res.Finite()=(V_finite/factor);
  m_res.IR()=(V_eps/factor);
  m_res.IR2()=(V_eps2/factor);
}

}


using namespace OpenLoops;

void dummyamp2func(double* moms, double* M2L0, double* M2L1, double* IRL1, double* M2L2, double* IRL2)
{ THROW(normal_exit, "Shopping list generated."); }
void dummypermfunc(int* permutation)
{ THROW(normal_exit, "Shopping list generated."); }

DECLARE_VIRTUALME2_GETTER(OpenLoops_Virtual,"OpenLoops_Virtual")
Virtual_ME2_Base *ATOOLS::Getter<Virtual_ME2_Base,Process_Info,OpenLoops_Virtual>::
operator()(const Process_Info &pi) const
{
  DEBUG_FUNC(pi);
  if (pi.m_loopgenerator!="OpenLoops") return NULL;
  if (pi.m_fi.m_nloewtype!=nlo_type::lo) return NULL;
  if (pi.m_fi.m_nloqcdtype!=nlo_type::loop) return NULL;
  if (MODEL::s_model->Name()!="SM") return NULL;

  Flavour_Vector flavs=pi.ExtractFlavours();
  Flavour_Vector map_flavs=OpenLoops_Interface::MapFlavours(flavs);
  msg_Tracking()<<endl<<flavs<<" --> "<<map_flavs<<endl;

  vector<int> permutation;
  string process=OpenLoops_Interface::GetProcessPermutation(map_flavs, permutation);
  pair<string, string> groupsub=OpenLoops_Interface::ScanFiles(process, pi.m_oew, pi.m_oqcd-1, 1);
  string grouptag=groupsub.first;
  string subid=groupsub.second;
  if (grouptag!="") {
    // symbols in fortran are always defined as lower case
    string lc_functag(grouptag+"_"+process+"_"+subid+"_");
    for (size_t i(0);i<lc_functag.length();++i)
      lc_functag[i]=tolower(lc_functag[i]);
    vector<string> suffixes;
    suffixes.push_back("1tL"); // deprecated
    suffixes.push_back("1L"); // deprecated
    suffixes.push_back("1ptL"); // deprecated
    suffixes.push_back("0L"); // deprecated
    suffixes.push_back("l");
    suffixes.push_back("lt");
    suffixes.push_back("lpt");
    suffixes.push_back("ls");
    suffixes.push_back("lp");
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
    msg_Info()<<"Trying OpenLoops for "<<flavs;
    msg_Debugging()<<" -> "<<lc_functag;
    msg_Info()<<" ... "<<std::endl;
    return new OpenLoops_Virtual(pi, flavs, (Amp2Func) ampfunc,
                                 (PermutationFunc) permfunc, permutation, lc_functag);
  }
  else {
    if (OpenLoops_Interface::s_generate_list) {
      if (OpenLoops_Interface::s_shoppinglist.find(process)==OpenLoops_Interface::s_shoppinglist.end()) {
        ofstream list("OL_list.m", ios::app);
        list<<"(* "<<process<<" *)"<<endl;
        list<<"SubProcess["<<OpenLoops_Interface::s_shoppinglist.size()+1<<"] = {\n"
            <<" FeynArtsProcess -> {"
            <<OpenLoops_Interface::s_particles[map_flavs[0].HepEvt()].m_faname<<", "
            <<OpenLoops_Interface::s_particles[map_flavs[1].HepEvt()].m_faname<<"} -> {";
        for (size_t i=2; i<map_flavs.size()-1; ++i) {
          list<<OpenLoops_Interface::s_particles[map_flavs[i].HepEvt()].m_faname<<", ";
        }
        list<<OpenLoops_Interface::s_particles[map_flavs[map_flavs.size()-1].HepEvt()].m_faname<<"},\n"
            <<" SelectCoupling -> (Exponent[#, gQCD] == "<<pi.m_oqcd-1<<"+2*#2 &),\n"
            <<" SortExternal -> True,\n"
            <<" InsertFieldsOptions -> {Restrictions -> {ExcludeParticles -> {S[2|3], SV}, NoQuarkMixing}}";
        list<<"\n};\n"<<endl;
        list.close();
        PRINT_INFO("Generated list entry for "<<process);
        OpenLoops_Interface::s_shoppinglist.insert(process);
      }
      return new OpenLoops_Virtual(pi, flavs, dummyamp2func,
                                   dummypermfunc, permutation, "dummy");
    }
  }

  return NULL;
}
