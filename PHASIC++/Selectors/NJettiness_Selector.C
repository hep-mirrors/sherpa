#include "PHASIC++/Selectors/Selector.H"
#include "ATOOLS/Phys/Fastjet_Helpers.H"

#include "PHASIC++/Main/Process_Integrator.H"
#include "ATOOLS/Org/Run_Parameter.H"

#include <numeric>

/*--------------------------------------------------------------------

 NJettiness arXiv:1004.2489

 --------------------------------------------------------------------*/

using namespace ATOOLS;

namespace PHASIC {
  class NJettiness_Selector: public Selector_Base {
  private:
    size_t m_N;
    double m_taucut;
  public:
    NJettiness_Selector(Process_Base* const proc,
                          ATOOLS::Scoped_Settings s);
    ~NJettiness_Selector() = default;

  protected:
  public:
    static void PrintCommonInfoLines(std::ostream& str, size_t width);

    bool Trigger(ATOOLS::Selector_List &sl) override;
    void BuildCuts(Cut_Data *) override;
  };

  NJettiness_Selector::NJettiness_Selector(Process_Base* const proc,
                                           ATOOLS::Scoped_Settings s)
    : Selector_Base("",proc)
    {
      DEBUG_INFO("Hello from n-jettiness");

      /// Operates on a fixed number of hard directions, N
      m_N = s["N"].SetDefault(0).Get<size_t>();

      /// Passes if the value of N-jettiness is above a given cut
      m_taucut = s["TauN"].SetDefault(1.).Get<double>();

      /// additional parameters might be frame, which underlying algo to
      /// determine the hard-directions... for the moment use anti-kt in
      /// in the lab frame
    }

  bool NJettiness_Selector::Trigger(ATOOLS::Selector_List &sl)
  {
    DEBUG_FUNC(p_proc->Flavours());

    /// fastjet clustering of fs particles
    std::vector<fjcore::PseudoJet> input, jets;

    /// prepare input partons
    for(size_t i{m_nin}; i< sl.size();++i){
      if(Flavour(kf_jet).Includes(sl[i].Flavour())){
        input.push_back(MakePseudoJet(sl[i].Flavour(), sl[i].Momentum()));
      }
    }
    if(m_N > input.size()){
      THROW(fatal_error, METHOD + " Cannot have more hard directions than partons");
    }

    if(m_N > 0){
      /// create cluster sequence and get the hard fs directions (these need to be added to beam)
      /// in a hadron-collider setup
      fjcore::ClusterSequence cs(input,
                                 fjcore::JetDefinition(fjcore::antikt_algorithm,1));
      jets = cs.exclusive_jets(static_cast<int>(m_N));
    }

    /// add beams if necessary
    for(size_t i{0}; i< m_nin; ++i){
      if(Flavour(kf_jet).Includes(sl[i].Flavour())) {
        jets.push_back(MakePseudoJet(sl[i].Flavour(), sl[i].Momentum()));
      }
    }

    /// compute TauN, i.e. the sum of the minimum for each FS partons w.r.t to the hard directions
    std::vector<double> TauN;
    for(auto& i: input){
      std::vector<double> TauN_i;
      for(auto& j : jets){
        TauN_i.push_back(fjcore::dot_product(j, i)/j.E()); /// there are of course different choices here one can make
      }
      TauN.push_back(*std::min_element(TauN_i.begin(),TauN_i.end()));
    }
    const double TauN_Val = std::accumulate(TauN.begin(), TauN.end(), 0.0);

    const bool trigger = (TauN_Val > m_taucut);

    if(!trigger){
      msg_Debugging() << "Point discarded by " << METHOD << std::endl;
    } else {
      msg_Debugging() << "Point passed: " << METHOD << std::endl;
    }

    return (1 - m_sel_log->Hit(1 - trigger));
  }



  void NJettiness_Selector::BuildCuts(Cut_Data *) {return;}
}

using namespace PHASIC;

DECLARE_GETTER(NJettiness_Selector,"NJettiness",Selector_Base,Selector_Key);

Selector_Base *ATOOLS::Getter<Selector_Base,Selector_Key,NJettiness_Selector>::
operator()(const Selector_Key &key) const
{
  auto s = key.m_settings["NJettiness"];

  return new NJettiness_Selector(key.p_proc, s);
}

void ATOOLS::Getter<Selector_Base,Selector_Key,NJettiness_Selector>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"NJettiness:\n";
  str<<width<<"  N: number of hard directions\n"
     <<width<<"  TauN: NJettiness-cut\n";
  str<<"cite: arXiv:1004.2489";
}
