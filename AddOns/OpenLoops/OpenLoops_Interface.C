#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include "MODEL/Main/Model_Base.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Library_Loader.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include <algorithm>
#include <sys/stat.h>

#include "OpenLoops_Interface.H"
#include "OpenLoops_Virtual.H"
#include "OpenLoops_Born.H"

using namespace PHASIC;
using namespace MODEL;
using namespace ATOOLS;
using namespace std;

namespace OpenLoops {

  std::map<int, OpenLoops_Particle> OpenLoops_Interface::s_particles;
  std::string OpenLoops_Interface::s_olprefix;
  std::vector<std::string> OpenLoops_Interface::s_allowed_libs;
  int OpenLoops_Interface::s_pole_mode;
  int OpenLoops_Interface::s_amp_switch;
  int OpenLoops_Interface::s_amp_switch_rescue;
  int OpenLoops_Interface::s_nf;
  int OpenLoops_Interface::s_nq_nondecoupled;
  bool OpenLoops_Interface::s_generate_list;
  double OpenLoops_Interface::s_overwrite_mb;
  set<string> OpenLoops_Interface::s_shoppinglist;

  bool OpenLoops_Interface::Initialize(const string &path,const string &file,
                                       MODEL::Model_Base *const model,
                                       BEAM::Beam_Spectra_Handler *const beam,
                                       PDF::ISR_Handler *const isr)
  {
    InitializeParticles();

    struct stat st;
    Data_Reader reader(" ",";","#","=");
    s_olprefix = rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/OpenLoops";
    if(stat(s_olprefix.c_str(),&st) != 0) s_olprefix = OPENLOOPS_PREFIX;
    s_olprefix = reader.GetValue<string>("OL_PREFIX", s_olprefix);

    Data_Reader liblist(" ",";","#","=");
    liblist.SetAddCommandLine(false);
    liblist.SetInputFile(s_olprefix+"/lib/librarylist.info");
    vector<vector<string> > ollibs;
    liblist.MatrixFromFile(ollibs);
    if (ollibs.size()==0) {
      THROW(fatal_error, "Did not find "+s_olprefix+"/lib/librarylist.info");
    }
    s_loader->AddPath(s_olprefix+"/lib");
    s_loader->AddPath(s_olprefix+"/proclib");
    for (size_t i=0; i<ollibs.size(); ++i) {
      if (!s_loader->LoadLibrary(ollibs[i][0])) {
        THROW(fatal_error, "Failed to load OpenLoops library "+ollibs[i][0]);
      }
    }

    reader.VectorFromFile(s_allowed_libs, "OL_ALLOWED_LIBS");
    s_pole_mode=reader.GetValue<int>("OL_POLE_MODE", 0);
    s_amp_switch=reader.GetValue<int>("OL_AMP_SWITCH", 1);
    s_amp_switch_rescue=reader.GetValue<int>("OL_AMP_SWITCH_RESCUE", 7);
    s_nf=reader.GetValue<int>("OL_NF", 6);
    s_nq_nondecoupled=reader.GetValue<int>("OL_RUNNING_FLAVOURS",Flavour(kf_quark).Size()/2);
    s_generate_list=reader.GetValue<size_t>("OL_GENERATE_LIST", false);
    s_overwrite_mb=reader.GetValue<double>("OL_OVERWRITE_MB", -1.0);


    OpenLoops_Virtual::SetInterface(this);
    OpenLoops_Born::SetInterface(this);

    int lenws=1200;
    char welcomestr[lenws];
    openloops_welcome_(welcomestr, &lenws);
    msg_Info()<<std::string(welcomestr,lenws)<<std::endl;
    PRINT_INFO("Initialised OpenLoops generator:");
    PRINT_VAR(s_olprefix);
    PRINT_VAR(s_amp_switch);
    PRINT_VAR(s_amp_switch_rescue);
    PRINT_VAR(s_pole_mode);
    PRINT_VAR(s_nf);
    PRINT_VAR(s_nq_nondecoupled);
    PRINT_VAR(s_allowed_libs);
    PRINT_VAR(s_overwrite_mb);

    MyStrStream cite;
    cite<<"The OpenLoops library~\\cite{Cascioli:2011va} of virtual"<<endl
        <<"matrix elements has been used. "<<endl;
    if (s_amp_switch==1 || s_amp_switch_rescue==1 ||
        s_amp_switch==7 || s_amp_switch_rescue==7) {
      cite<<"It is partly based on the tensor integral reduction described "<<endl
          <<"in~\\cite{Denner:2002ii,Denner:2005nn,Denner:2010tr}."<<endl;
    }
    rpa->gen.AddCitation(1,cite.str());

    return true;
  }

  void OpenLoops_Interface::InitializeParticles()
  {
    s_particles.clear();
    s_particles[kf_nue] = OpenLoops_Particle("ne", "F[1,{1}]", 0);
    s_particles[-kf_nue] = OpenLoops_Particle("nex", "-F[1,{1}]", 1);
    s_particles[kf_numu] = OpenLoops_Particle("nm", "F[1,{2}]", 2);
    s_particles[-kf_numu] = OpenLoops_Particle("nmx", "-F[1,{2}]", 3);
    s_particles[kf_nutau] = OpenLoops_Particle("nl", "F[1,{3}]", 4);
    s_particles[-kf_nutau] = OpenLoops_Particle("nlx", "-F[1,{3}]", 5);
    s_particles[kf_e] = OpenLoops_Particle("e", "F[2,{1}]", 6);
    s_particles[-kf_e] = OpenLoops_Particle("ex", "-F[2,{1}]", 7);
    s_particles[kf_mu] = OpenLoops_Particle("m", "F[2,{2}]", 8);
    s_particles[-kf_mu] = OpenLoops_Particle("mx", "-F[2,{2}]", 9);
    s_particles[kf_tau] = OpenLoops_Particle("l", "F[2,{3}]", 10);
    s_particles[-kf_tau] = OpenLoops_Particle("lx", "-F[2,{3}]", 11);
    s_particles[kf_u] = OpenLoops_Particle("u", "F[3,{1}]", 12);
    s_particles[-kf_u] = OpenLoops_Particle("ux", "-F[3,{1}]", 13);
    s_particles[kf_c] = OpenLoops_Particle("c", "F[3,{2}]", 14);
    s_particles[-kf_c] = OpenLoops_Particle("cx", "-F[3,{2}]", 15);
    s_particles[kf_t] = OpenLoops_Particle("t", "F[3,{3}]", 16);
    s_particles[-kf_t] = OpenLoops_Particle("tx", "-F[3,{3}]", 17);
    s_particles[kf_d] = OpenLoops_Particle("d", "F[4,{1}]", 18);
    s_particles[-kf_d] = OpenLoops_Particle("dx", "-F[4,{1}]", 19);
    s_particles[kf_s] = OpenLoops_Particle("s", "F[4,{2}]", 20);
    s_particles[-kf_s] = OpenLoops_Particle("sx", "-F[4,{2}]", 21);
    s_particles[kf_b] = OpenLoops_Particle("b", "F[4,{3}]", 22);
    s_particles[-kf_b] = OpenLoops_Particle("bx", "-F[4,{3}]", 23);
    s_particles[kf_h0] = OpenLoops_Particle("h", "S[1]", 24);
    s_particles[kf_photon] = OpenLoops_Particle("a", "V[1]", 25);
    s_particles[kf_Z] = OpenLoops_Particle("z", "V[2]", 26);
    s_particles[-kf_Wplus] = OpenLoops_Particle("w", "V[3]", 27);
    s_particles[kf_Wplus] = OpenLoops_Particle("wx", "-V[3]", 28);
    s_particles[kf_gluon] = OpenLoops_Particle("g", "V[5]", 29);
  }

  void OpenLoops_Interface::OpenLoopsInit(double mur2,
                                          double alpha_QED, double alpha_S)
  {
    double Mass_E=Flavour(kf_e).Mass();
    double Mass_M=Flavour(kf_mu).Mass();
    double Mass_L=Flavour(kf_tau).Mass();
    double Mass_U=Flavour(kf_u).Mass();
    double Mass_D=Flavour(kf_d).Mass();
    double Mass_S=Flavour(kf_s).Mass();
    double Mass_C=Flavour(kf_c).Mass();
    double Mass_B=s_overwrite_mb<0.0 ? Flavour(kf_b).Mass() : s_overwrite_mb;
    double Mass_T=Flavour(kf_t).Mass();
    double Mass_W=Flavour(kf_Wplus).Mass();
    double Mass_Z=Flavour(kf_Z).Mass();
    double Mass_H=Flavour(kf_h0).Mass();
    double Width_C=Flavour(kf_c).Width();
    double Width_B=Flavour(kf_b).Width();
    double Width_T=Flavour(kf_t).Width();
    double Width_W=Flavour(kf_Wplus).Width();
    double Width_Z=Flavour(kf_Z).Width();
    double Width_H=Flavour(kf_h0).Width();
    int last_switch=1;
    int amp_switch=s_amp_switch;
    int amp_switch_rescue=s_amp_switch_rescue;
    int use_coli_cache(true);
    int check_Ward_tree=false;
    int check_Ward_loop=false;
    int out_symmetry=true;
    int leading_colour=false;
    parameters_init_(&Mass_E, &Mass_M, &Mass_L, &Mass_U, &Mass_D, &Mass_S, &Mass_C, &Width_C, &Mass_B, &Width_B, &Mass_T, &Width_T,
                     &Mass_W, &Width_W, &Mass_Z, &Width_Z, &Mass_H, &Width_H, &alpha_QED, &alpha_S,
                     &last_switch, &amp_switch, &amp_switch_rescue, &use_coli_cache,
                     &check_Ward_tree, &check_Ward_loop, &out_symmetry, &leading_colour);


    double renscale=sqrt(mur2);
    double fact_UV=1.0;
    double fact_IR=1.0;
    double pole1_UV=0.0;
    double pole1_IR=0.0;
    double pole2_IR=0.0;
    int polenorm_swi=0;
    double opp_rootsvalue=1000.0;
    double opp_limitvalue=0.01;
    double opp_thrs=0.000001;
    int opp_idig=0;
    int opp_scaloop=2;
    int sam_isca=2;
    int sam_verbosity=0;
    int sam_itest=0;
    int fermion_loops=1;
    int nonfermion_loops=1;
    int CT_on=1;
    int R2_on=1;
    int IR_on=1;
    int polemode=(s_pole_mode>0?1:0);
    double set_C_PV_threshold=1.E-10;
    double set_D_PV_threshold=1.E-10;
    int set_DD_red_mode=2;
    loop_parameters_init_(&renscale, &fact_UV, &fact_IR, &pole1_UV, &pole1_IR, &pole2_IR, &polenorm_swi, &s_nf, &s_nq_nondecoupled,
                          &opp_rootsvalue, &opp_limitvalue, &opp_thrs, &opp_idig, &opp_scaloop,
                          &sam_isca, &sam_verbosity, &sam_itest, &fermion_loops, &nonfermion_loops, &CT_on, &R2_on, &IR_on, &polemode,
                          &set_C_PV_threshold, &set_D_PV_threshold, &set_DD_red_mode);
  }

  std::string OpenLoops_Interface::MatchOptions(vector<string> options, int oew, int oqcd, int nloop) {
    for (size_t i=2; i<options.size(); ++i) {
      string option=options[i].substr(0, options[i].find("="));
      string value=options[i].substr(options[i].find("=")+1);

      if (option=="EW" && value!=ToString(oew)+",0") return "0";
      if (option=="QCD" && value!=ToString(oqcd)+","+ToString(nloop)) return "0";
      if (option=="CKMORDER") {
        int ckmorder=ToType<int>(value);
        if (ckmorder<3) {
          if (s_model->ComplexMatrixElement("CKM", 0,2)!=Complex(0.0,0.0) ||
              s_model->ComplexMatrixElement("CKM", 2,0)!=Complex(0.0,0.0)) {
            return "0";
          }
        }
        if (ckmorder<2) {
          if (s_model->ComplexMatrixElement("CKM", 1,2)!=Complex(0.0,0.0) ||
              s_model->ComplexMatrixElement("CKM", 2,1)!=Complex(0.0,0.0)) {
            return "0";
          }
        }
        if (ckmorder<1) {
          if (s_model->ComplexMatrixElement("CKM", 0,1)!=Complex(0.0,0.0) ||
              s_model->ComplexMatrixElement("CKM", 1,0)!=Complex(0.0,0.0)) {
            return "0";
          }
        }
      }
      if (option=="nf" && ToType<int>(value)!=s_nf) return "0";
      if (option=="MD" && Flavour(kf_d).Mass()>0.0) return "0";
      if (option=="MU" && Flavour(kf_u).Mass()>0.0) return "0";
      if (option=="MS" && Flavour(kf_s).Mass()>0.0) return "0";
      if (option=="MC" && Flavour(kf_c).Mass()>0.0) return "0";
      if (option=="MB" && Flavour(kf_b).Mass()>0.0 && s_overwrite_mb<0.0) return "0";
      if (option=="MT" && Flavour(kf_t).Mass()>0.0) return "0";
      if (option=="ME" && Flavour(kf_e).Mass()>0.0) return "0";
      if (option=="MM" && Flavour(kf_mu).Mass()>0.0) return "0";
      if (option=="MT" && Flavour(kf_tau).Mass()>0.0) return "0";

      if (option=="map") return value;
    }
    return "1";
  }


  int OpenLoops_Interface::SelectInfo(const dirent *entry)
  {
    string fname(entry->d_name);
    if (fname.find(".info")!=string::npos) {
      return true;
    }
    else return false;
  }

  pair<string, string> OpenLoops_Interface::ScanFiles(string& process, int oew, int oqcd, int nloop)
  {
    struct dirent **entries;
    string procdatapath=s_olprefix+"/proclib";
    int n(scandir(procdatapath.c_str(),&entries,&SelectInfo,alphasort));
    if (n<0) THROW(fatal_error, "OpenLoops process dir "+procdatapath+" not found.");
    vector<string> files;
    for (int ifile=0; ifile<n; ++ifile) {
      files.push_back(string(entries[ifile]->d_name));
      free(entries[ifile]);
    }
    free(entries);

    for (int ifile=0; ifile<files.size(); ++ifile) {
      Data_Reader reader(" ",";","#","");
      reader.SetAddCommandLine(false);
      reader.SetInputFile(procdatapath+"/"+files[ifile]);
      reader.SetMatrixType(mtc::transposed);
      vector<vector<string> > content;
      reader.MatrixFromFile(content);
      for (size_t i=0; i<content.size(); ++i) {
        if (content[i][1]==process) {
          string match=MatchOptions(content[i], oew, oqcd, nloop);
          if (match=="0") {
            PRINT_INFO("Ignoring process with incompatible options.");
            continue;
          }
          else if (match!="1") {
            PRINT_INFO("Mapping "<<process<<" to "<<match<<".");
            process=match;
            return ScanFiles(process, oew, oqcd, nloop);
          }
          string grouptag=content[i][0];
          string process_subid=content[i][2];
          if (s_allowed_libs.size()>0) {
            bool allowed=false;
            for (size_t i=0; i<s_allowed_libs.size(); ++i) {
              if (grouptag==s_allowed_libs[i]) allowed=true;
            }
            if (!allowed) continue;
          }
          return make_pair(grouptag, process_subid);
        }
      }
    }
    PRINT_INFO("Didn't find info file matching process "<<process);
    return make_pair("", "0");
  }


  bool OpenLoops_Interface::Order(const pair<size_t, Flavour> &a,const pair<size_t, Flavour> &b)
  {
    if (s_particles.find(a.second.HepEvt())==s_particles.end() ||
        s_particles.find(b.second.HepEvt())==s_particles.end()) {
      THROW(fatal_error, "Unknown particle.");
    }
    return s_particles[a.second.HepEvt()].m_order<s_particles[b.second.HepEvt()].m_order;
  }


  string OpenLoops_Interface::GetProcessPermutation(const Flavour_Vector& flavs_orig,
                                                    vector<int>& permutation) {
    vector<pair<size_t, Flavour> > flavs_ol(flavs_orig.size());
    for (size_t i=0; i<2; ++i) {
      flavs_ol[i]=make_pair(i,flavs_orig[i]);
    }
    for (size_t i=2; i<flavs_orig.size(); ++i) {
      flavs_ol[i]=make_pair(i,flavs_orig[i].Bar());
    }

    sort(flavs_ol.begin(), flavs_ol.end(), Order);

    permutation.clear();
    permutation.resize(flavs_ol.size());
    for (size_t i=0; i<flavs_ol.size(); ++i) {
      for (size_t j=0; j<flavs_ol.size(); ++j) {
        if (i==flavs_ol[j].first) {
          permutation[i]=int(j+1); // +1 for fortran openloops
          break;
        }
      }
    }

    string process;
    for (size_t i=0; i<flavs_ol.size(); ++i) {
      process += s_particles[flavs_ol[i].second.HepEvt()].m_name;
    }
    return process;
  }

  bool SortByFirstDec(pair<size_t, size_t> p1, pair<size_t, size_t> p2) {
    return p1.first>p2.first;
  }

  Flavour_Vector OpenLoops_Interface::MapFlavours(const Flavour_Vector& orig)
  {
    for (size_t i=2; i<orig.size(); ++i) {
      if (orig[i].Width()!=0.0) {
        THROW(fatal_error, "Non-zero width of final state particle.");
      }
    }

    /* Concept:
      For each family i=0,...,2:
      (1) Given a final state, determine the four (anti)lepton/neutrino
          multiplicities in the given process:
          a_nubar, a_nu, a_lbar, a_l

      (2) Compute the discriminant N[i] as
          N[i] = Ngen - (i+1) + Ngen*(a_nubar + a_nu*Nmax + a_lbar*Nmax^2 + a_l*Nmax^3)
          where Nmax should be chosen such that a_...<=Nmax.
          In practice one can safely set Nmax=10 and it will work for any
          process with <= 20 final-state leptons.
          It is also convenient to set Ngen=10, although i runs only from 1 to 3.

      (3) Reassign the lepton generations with a permutation
          p1 -> 1, p2 -> 2, p3 -> 3  such that  N[p1] > N[p2] > N[p3] */

    multiset<int> hepevt;
    for (size_t i=0; i<orig.size(); ++i) { hepevt.insert(orig[i].HepEvt()); }

    size_t Ngen(10), Nmax(10);
    vector<pair<size_t, size_t> > N(3);
    for (size_t i=0; i<3; ++i) {
      int nu_gen=10+2*(i+1);
      int l_gen=9+2*(i+1);

      size_t a_nu=hepevt.count(nu_gen);
      size_t a_nubar=hepevt.count(-nu_gen);
      size_t a_l=hepevt.count(l_gen);
      size_t a_lbar=hepevt.count(-l_gen);

      N[i]=make_pair(Ngen-(i+1) + Ngen*
                     (a_nubar+a_nu*Nmax+a_lbar*Nmax*Nmax+a_l*Nmax*Nmax*Nmax),
                     i);
    }

    sort(N.begin(), N.end(), SortByFirstDec);

    Flavour_Vector ret(orig);
    for (size_t i=0; i<3; ++i) {
      int nu_gen=10+2*(N[i].second+1);
      int nu_gen_new=10+2*(i+1);
      int l_gen=9+2*(N[i].second+1);
      int l_gen_new=9+2*(i+1);

      for (size_t j=0; j<orig.size(); ++j) {
        if (orig[j].Kfcode()==nu_gen) {
          ret[j]=Flavour(nu_gen_new, orig[j].IsAnti());
        }
        if (orig[j].Kfcode()==l_gen) {
          ret[j]=Flavour(l_gen_new, orig[j].IsAnti());
        }
      }
    }
    return ret;
  }

  void dummyamp2func(double* moms, double* M2L0, double* M2L1, double* IRL1,
                     double* M2L2, double* IRL2)
  {
    THROW(normal_exit, "Shopping list generated.");
  }
  
  void dummypermfunc(int* permutation)
  {
    THROW(normal_exit, "Shopping list generated.");
  }

}

using namespace OpenLoops;

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
      suffixes.push_back("1t");
      suffixes.push_back("1");
      suffixes.push_back("1pt");
      suffixes.push_back("0");
      void *ampfunc, *permfunc;
      for (size_t i=0; i<suffixes.size(); ++i) {
        string libraryfile="openloops_"+grouptag+"_"+suffixes[i]+"L";
        ampfunc=s_loader->GetLibraryFunction(libraryfile,"vamp2_"+lc_functag);
        permfunc=s_loader->GetLibraryFunction(libraryfile,"set_permutation_"+lc_functag);
        if (ampfunc!=NULL && permfunc!=NULL) break;
      }
      if (ampfunc==NULL || permfunc==NULL) {
        PRINT_INFO("Didn't find functions");
        return NULL;
      }

      msg_Info()<<endl;
      PRINT_INFO("Initialising OpenLoops Virtual for "<<flavs<<": "<<lc_functag);
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
              <<" SelectCoupling -> (Exponent[#, gQCD] >= "<<pi.m_oqcd-1<<"+2*#2 &),\n"
              <<" SortExternal -> True,\n"
              <<" InsertFieldsOptions -> {Restrictions -> {ExcludeParticles -> {S[2|3], SV}, NoQuarkMixing}}";
          if (OpenLoops_Interface::s_nf!=6) {
            list<<",\n SetParameters -> JoinOptions[{nf -> "<<OpenLoops_Interface::s_nf<<"}]";
          }
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
      suffixes.push_back("1s");
      void *ampfunc, *permfunc;
      for (size_t i=0; i<suffixes.size(); ++i) {
        string libraryfile="openloops_"+grouptag+"_"+suffixes[i]+"L";
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

  DECLARE_GETTER(OpenLoops_Interface,"OpenLoops",ME_Generator_Base,ME_Generator_Key);

  ME_Generator_Base *ATOOLS::Getter<ME_Generator_Base,ME_Generator_Key,
				    OpenLoops_Interface>::
  operator()(const ME_Generator_Key &key) const
  {
    return new OpenLoops::OpenLoops_Interface();
  }

  void ATOOLS::Getter<ME_Generator_Base,ME_Generator_Key,OpenLoops_Interface>::
  PrintInfo(ostream &str,const size_t width) const
  { 
    str<<"Interface to the OpenLoops loop ME generator"; 
  }

