#include "PHASIC++/EWSudakov/Clustered_EWSudakov_Calculator.H"

#include "ATOOLS/Org/Run_Parameter.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PHASIC++/Process/ME_Generators.H"
#include "SHERPA/PerturbativePhysics/Matrix_Element_Handler.H"
#include "SHERPA/SoftPhysics/Resonance_Finder.H"

using namespace ATOOLS;
using namespace PHASIC;
using namespace SHERPA;

Clustered_EWSudakov_Calculator::Clustered_EWSudakov_Calculator(Process_Base* _proc)
  : proc{_proc}
{
  Scoped_Settings meqedsettings{
    Settings::GetMainSettings()["ME_QED"] };
  m_resdist =
    meqedsettings["CLUSTERING_THRESHOLD"].SetDefault(10.0).Get<double>();
  m_disabled =
    Settings::GetMainSettings()["EWSUDAKOV_CLUSTERING_DISABLED"].SetDefault(false).Get<bool>();
  m_zzhack =
    Settings::GetMainSettings()["EWSUDAKOV_ZZHACK_ENABLED"].SetDefault(false).Get<bool>();

  //auto level = msg->Level();
  //msg->SetLevel(15);

  const Flavour_Vector& flavs = proc->Flavours();

  if (m_disabled) {
    calculators.emplace(
        std::make_pair(flavs, new EWSudakov_Calculator{proc}));
  } else {
    AddCalculators(flavs, 0);
    for (size_t i {0}; i < flavs.size(); ++i)
      msg_Debugging() << flavs[i] << ' ';
    msg_Debugging() << "added " << calculators.size() << " calculators\n";
  }

  //msg->SetLevel(level);
}

void Clustered_EWSudakov_Calculator::AddCalculators(const Flavour_Vector& flavs, size_t clusterings)
{
  DEBUG_FUNC(clusterings);
  size_t nflavs {flavs.size()};
  for (size_t i {0}; i < nflavs; ++i)
    msg_Debugging() << flavs[i];
  msg_Debugging() << '\n';
  for (size_t i {0}; i < nflavs; ++i) {
    if (flavs[i].IsLepton()) {
      for (size_t j {i + 1}; j < nflavs; ++j) {
        if (flavs[j] == flavs[i].Bar()) {
          Flavour_Vector newflavs = flavs;
          newflavs[i] = Flavour{kf_Z};
          newflavs.erase(newflavs.begin() + j);
          AddCalculators(newflavs, clusterings + 1);
        } else if (flavs[j] == flavs[i].IsoWeakPartner()) {
          Flavour_Vector newflavs = flavs;
          newflavs[i] = Flavour{kf_Wplus, (flavs[i].Charge() + flavs[j].Charge() < 0)};
          newflavs.erase(newflavs.begin() + j);
          AddCalculators(newflavs, clusterings + 1);
        }
      }
    }
  }
  AddCalculator(flavs, clusterings);
}

Clustered_EWSudakov_Calculator::~Clustered_EWSudakov_Calculator()
{
}

void Clustered_EWSudakov_Calculator::AddCalculator(const Flavour_Vector& flavs, size_t clusterings)
{
  auto it = calculators.find(flavs);
  if (it != calculators.end())
    return;

  // build process info for clustered process
  Process_Info pi;
  pi.m_addname = "__Sudakov";
  pi.m_megenerator = "Comix";
  for (size_t i{0}; i < proc->NIn(); ++i) {
    pi.m_ii.m_ps.push_back(Subprocess_Info(flavs[i], "", ""));
  }
  for (size_t i{proc->NIn()}; i < proc->NIn() + proc->NOut() - clusterings; ++i) {
    pi.m_fi.m_ps.push_back(Subprocess_Info(flavs[i], "", ""));
  }
  pi.m_maxcpl = proc->Info().m_maxcpl;
  pi.m_mincpl = proc->Info().m_mincpl;
  pi.m_maxacpl = proc->Info().m_maxacpl;
  pi.m_minacpl = proc->Info().m_minacpl;
  pi.m_mincpl[1] -= clusterings;
  pi.m_maxcpl[1] -= clusterings;

  // initialize process
  auto clustered_proc =
    proc->Generator()->Generators()->InitializeProcess(pi, false);
  if (clustered_proc == NULL) {
    msg_Debugging()
      << "WARNING: Clustered_EWSudakov_Calculator::AddCalculator can not"
      << "initialize process for process info: " << pi << '\n';
    return;
  }
  clustered_proc->SetSelector(Selector_Key{});
  clustered_proc->SetScale(Scale_Setter_Arguments(
        MODEL::s_model, "VAR{" + ToString(sqr(rpa->gen.Ecms())) + "}",
        "Alpha_QCD 1"));
  clustered_proc->SetKFactor(KFactor_Setter_Arguments("None"));
  msg_Debugging()
    << "Clustered_EWSudakov_Calculator::AddCalculator initialized "
    << clustered_proc->Name() << '\n';

  // add calculator
  calculators.emplace(
      std::make_pair(flavs, new EWSudakov_Calculator{clustered_proc}));
}

EWSudakov_Log_Corrections_Map
Clustered_EWSudakov_Calculator::CorrectionsMap(Vec4D_Vector mom)
{
  auto level = msg->Level();
  //msg->SetLevel(15);

  if (m_disabled) {
    return calculators.begin()->second->CorrectionsMap(mom);
  }

  Flavour_Vector flavs = proc->Flavours();

  if (m_zzhack) {
    // TODO: delete this hack and the corresponding switch later

    assert(flavs[2].Kfcode() == kf_e);
    assert(flavs[3].Kfcode() == kf_mu);
    assert(flavs[4].Kfcode() == kf_e);
    assert(flavs[5].Kfcode() == kf_mu);

    flavs[3] = Flavour{kf_Z};
    mom[3] += mom[5];
    flavs.erase(flavs.begin() + 5);
    mom.erase(mom.begin() + 5);

    flavs[2] = Flavour{kf_Z};
    mom[2] += mom[4];
    flavs.erase(flavs.begin() + 4);
    mom.erase(mom.begin() + 4);

  } else {

    msg_Debugging() << "Will use input process for EWSudakov: ";
    for (const auto& flav : flavs) {
      msg_Debugging() << flav.ShellName() << " ";
    }
    msg_Debugging() << '\n';

    size_t nflavs {flavs.size()};
    std::map<double, std::vector<long int>> restab;
    for (long int i {0}; i < nflavs; ++i) {
      if (flavs[i].IsLepton()) {
        for (long int j {i + 1}; j < nflavs; ++j) {
          long int kf {kf_none};
          if (flavs[j] == flavs[i].Bar()) {
            kf = kf_Z;
          } else if (flavs[j] == flavs[i].IsoWeakPartner()) {
            kf = kf_Wplus;
          }
          if (kf != kf_none) {
            Flavour resonance{kf};
            double invariant_mass {(mom[i]+mom[j]).Mass()};
            double mdist{
              std::abs(invariant_mass - resonance.Mass()) / resonance.Width()};
            msg_Debugging()<<"found resonance candidate " << i << ", " << j << " -> " << kf << " (" << mdist << ")";
            if (mdist < m_resdist) {
              msg_Debugging()<<"-> accept\n";
              long int ida[3]={i,j,kf};
              restab[mdist]=std::vector<long int>(ida,ida+3);
            } else {
              msg_Debugging()<<"-> reject\n";
            }
          }
        }
      }
    }

    // replace resonances starting with the least off-shell one
    std::vector<double> clusterings;
    std::unordered_set<size_t> clustered_indizes;
    for (const auto& mdist_ida_pair : restab) {
      if (clustered_indizes.find(mdist_ida_pair.second[0]) == clustered_indizes.end()
          && clustered_indizes.find(mdist_ida_pair.second[1]) == clustered_indizes.end()) {
        clusterings.push_back(mdist_ida_pair.first);
        clustered_indizes.insert(mdist_ida_pair.second[0]);
        clustered_indizes.insert(mdist_ida_pair.second[1]);
      }
    }
    std::set<size_t> removelist;
    for (const auto& mdist : clusterings) {
      if (restab[mdist][2] == kf_Wplus) {
        flavs[restab[mdist][0]] =
          Flavour{kf_Wplus, (flavs[restab[mdist][1]].Charge() + flavs[restab[mdist][2]].Charge() < 0)};
      } else {
        flavs[restab[mdist][0]] =
          Flavour{kf_Z};
      }
      mom[restab[mdist][0]] += mom[restab[mdist][1]];
      removelist.insert(restab[mdist][1]);
    }
    for (auto i = removelist.rbegin(); i != removelist.rend(); ++ i) {
      flavs.erase(flavs.begin() + *i);
      mom.erase(mom.begin() + *i);
    }
    msg_Debugging() << "Will use clustered process for EWSudakov: ";
    for (const auto& flav : flavs) {
      msg_Debugging() << flav.ShellName() << " ";
    }
    msg_Debugging() << '\n';
  }

  msg->SetLevel(level);

  assert(calculators.find(flavs) != calculators.end());
  return calculators[flavs]->CorrectionsMap(mom);
}
