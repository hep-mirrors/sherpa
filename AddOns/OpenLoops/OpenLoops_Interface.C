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


extern "C" {

  void openloops_welcome_(char* str, int* len_str);

  void loop_parameters_write_();

  void ol_getparameter_double_c_(const char* key, double* val, int* err);
  void ol_getparameter_int_c_(const char* key, int* val, int* err);
  void ol_setparameter_double_c_(const char* key, double* val, int* err);
  void ol_setparameter_int_c_(const char* key, int* val, int* err);
  void ol_setparameter_string_c_(const char* key, const char* val, int* err);
  void ol_parameters_flush_();

}


namespace OpenLoops {

  std::map<int, OpenLoops_Particle> OpenLoops_Interface::s_particles;
  std::string OpenLoops_Interface::s_olprefix;
  std::vector<std::string> OpenLoops_Interface::s_allowed_libs;
  bool OpenLoops_Interface::s_generate_list;
  set<string> OpenLoops_Interface::s_shoppinglist;

  bool OpenLoops_Interface::Initialize(const string &path,const string &file,
                                       MODEL::Model_Base *const model,
                                       BEAM::Beam_Spectra_Handler *const beam,
                                       PDF::ISR_Handler *const isr)
  {
    InitializeParticles();

    // find OL installation prefix with several overwrite options
    struct stat st;
    Data_Reader reader(" ",";","#","=");
    s_olprefix = rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/OpenLoops";
    if(stat(s_olprefix.c_str(),&st) != 0) s_olprefix = OPENLOOPS_PREFIX;
    s_olprefix = reader.GetValue<string>("OL_PREFIX", s_olprefix);
    msg_Info()<<"Initialising OpenLoops generator from "<<s_olprefix<<endl;

    // load library dynamically
    s_loader->AddPath(s_olprefix+"/lib");
    s_loader->AddPath(s_olprefix+"/proclib");
    if (!s_loader->LoadLibrary("openloops")) THROW(fatal_error, "Failed to load libopenloops.");

    // set particle masses/widths
    int tmparr[] = {kf_e, kf_mu, kf_tau, kf_u, kf_d, kf_s, kf_c, kf_b, kf_t, kf_Wplus, kf_Z, kf_h0};
    vector<int> pdgids (tmparr, tmparr + sizeof(tmparr) / sizeof(tmparr[0]) );
    int dummy;
    for (size_t i=0; i<pdgids.size(); ++i) {
      if (Flavour(pdgids[i]).Mass()>0.0) SetParameter("mass("+ToString(pdgids[i])+")", Flavour(pdgids[i]).Mass());
      if (Flavour(pdgids[i]).Width()>0.0) SetParameter("width("+ToString(pdgids[i])+")", Flavour(pdgids[i]).Width());
    }

    // set remaining OL parameters specified by user
    vector<string> parameters;
    reader.VectorFromFile(parameters,"OL_PARAMETERS");
    for (size_t i=1; i<parameters.size(); i=i+2) SetParameter(parameters[i-1], parameters[i]);
    ol_parameters_flush_();

    reader.VectorFromFile(s_allowed_libs, "OL_ALLOWED_LIBS");
    s_generate_list=reader.GetValue<size_t>("OL_GENERATE_LIST", false);
    if (s_allowed_libs.size()) PRINT_VAR(s_allowed_libs);

    int lenws=1200;
    char welcomestr[lenws];
    openloops_welcome_(welcomestr, &lenws);
    msg_Info()<<std::string(welcomestr,lenws)<<std::endl;
    loop_parameters_write_();

    MyStrStream cite;
    cite<<"The OpenLoops library~\\cite{Cascioli:2011va} of virtual"<<endl
        <<"matrix elements has been used. "<<endl
        <<"It is partly based on the tensor integral reduction described "<<endl
        <<"in~\\cite{Denner:2002ii,Denner:2005nn,Denner:2010tr}."<<endl;
    rpa->gen.AddCitation(1,cite.str());

    OpenLoops_Virtual::SetInterface(this);
    OpenLoops_Born::SetInterface(this);
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
      const string nfs("nf");
      double nfsv=GetIntParameter(nfs);
      if (option=="nf" && ToType<int>(value)!=nfsv) return "0";
      if (option=="MD" && (double)GetDoubleParameter("mass(1)")>0.0) return "0";
      if (option=="MU" && (double)GetDoubleParameter("mass(2)")>0.0) return "0";
      if (option=="MS" && (double)GetDoubleParameter("mass(3)")>0.0) return "0";
      if (option=="MC" && (double)GetDoubleParameter("mass(4)")>0.0) return "0";
      if (option=="MB" && (double)GetDoubleParameter("mass(5)")>0.0) return "0";
      if (option=="MT" && (double)GetDoubleParameter("mass(6)")>0.0) return "0";
      if (option=="ME" && (double)GetDoubleParameter("mass(11)")>0.0) return "0";
      if (option=="MM" && (double)GetDoubleParameter("mass(13)")>0.0) return "0";
      if (option=="MT" && (double)GetDoubleParameter("mass(15)")>0.0) return "0";

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

        bool autoload=true;
        if (content[i].size()>0 && content[i][0]=="options") {
          for (size_t j=0; j<content[i].size(); ++j) {
            if (content[i][j]=="autoload=0") {
              autoload=false;
            }
          }
        }

        if (content[i].size()>2 && content[i][1]==process) {
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
          if (s_allowed_libs.size()>0 || autoload==false) {
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

  double OpenLoops_Interface::GetDoubleParameter(const std::string & key) {
    int err(0);
    double value;
    ol_getparameter_double_c_(key.c_str(), &value, &err);
    return value;
  }
  int OpenLoops_Interface::GetIntParameter(const std::string & key) {
    int err(0);
    int value;
    ol_getparameter_int_c_(key.c_str(), &value, &err);
    return value;
  }
  template <class ValueType>
  void HandleParameterStatus(int err, const std::string & key, ValueType value) {
    if (err==0) {
      msg_Tracking()<<"Setting OpenLoops parameter: "<<key<<" = "<<value<<endl;
    }
    else if (err==1) {
      THROW(fatal_error, "Unknown OpenLoops parameter: "+key+" = "+ToString(value));
    }
    else if (err==2) {
      THROW(fatal_error, "Error setting OpenLoops parameter: "+key+" = "+ToString(value));
    }
  }
  void OpenLoops_Interface::SetParameter(const std::string & key, double value) {
    int err(0);
    ol_setparameter_double_c_(key.c_str(), &value, &err);
    HandleParameterStatus(err, key, value);
  }
  void OpenLoops_Interface::SetParameter(const std::string & key, int value) {
    int err(0);
    ol_setparameter_int_c_(key.c_str(), &value, &err);
    HandleParameterStatus(err, key, value);
  }
  void OpenLoops_Interface::SetParameter(const std::string & key, std::string value) {
    int err(0);
    ol_setparameter_string_c_(key.c_str(), value.c_str(), &err);
    HandleParameterStatus(err, key, value);
  }
  void OpenLoops_Interface::FlushParameters() {
    ol_parameters_flush_();
  }


}

using namespace OpenLoops;

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

