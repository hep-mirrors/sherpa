#include "SHERPA/Tools/Userhook_Base.H"
#include "ATOOLS/Org/Message.H"
#include "SHERPA/Single_Events/Event_Handler.H"

#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Phys/Weight_Info.H"
#include <algorithm>
#include "MODEL/Main/Running_AlphaS.H"
#include "SHERPA/Main/Sherpa.H"
#include "SHERPA/Initialization/Initialization_Handler.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "CSSHOWER++/Main/CS_Shower.H"
#include "PDF/Main/Shower_Base.H"
#include "PHASIC++/Process/MCatNLO_Process.H"

/*
 * counter-terms applied to the X+bb events in a fused sample.
 */

using namespace ATOOLS;
using namespace SHERPA;

class Fusing_Direct_Hook : public Userhook_Base {

private:
  MODEL::Running_AlphaS *p_as;
  Sherpa* p_sherpa;
  double m_factor;
public:

  Fusing_Direct_Hook(const Userhook_Arguments args) :
    Userhook_Base("Fusing_Direct"), p_as(MODEL::as), p_sherpa(args.p_sherpa)
  {
    msg_Debugging()<<"Fusing_Direct Hook active."<<std::endl;
    Settings& s = Settings::GetMainSettings();
    m_factor = s["FUSING_DIRECT_FACTOR"].SetDefault(2.).Get<double>();
  }

  ~Fusing_Direct_Hook() {}

  ATOOLS::Return_Value::code Run(ATOOLS::Blob_List* blobs) {
    DEBUG_FUNC(p_sherpa->GetInitHandler()->GetMatrixElementHandler()->Process()->Parent()->Name());

    auto me_w_info = (*blobs->FindFirst(btp::Signal_Process))["MEWeightInfo"]->Get<ME_Weight_Info*>();
    double new_weight_factor = 1.0;

    //apply only to S-Events
    if (me_w_info->m_type != (ATOOLS::mewgttype::METS |ATOOLS::mewgttype::H)){

      PHASIC::MCatNLO_Process * mcproc = dynamic_cast<PHASIC::MCatNLO_Process* >
        (p_sherpa->GetInitHandler()->GetMatrixElementHandler()->Process()->Parent());
      if (mcproc==NULL){
        THROW(fatal_error,"no MC@NLO process found! For use with separate LO-Process, use K-Factor!");
      }
      // TODO: make more efficient, without whole amplitude copying!
      Cluster_Amplitude * ampl = mcproc->GetAmplitude();
      double muf2 = ampl->KT2();
      ampl->Delete();

      double mur2 = me_w_info->m_mur2;
      ATOOLS::Flavour bquark(kf_b);
      double mb2 = bquark.Mass()*bquark.Mass();
      double alphas((*p_as)(mur2));
      double TR = 0.5;
      double correction = 0.0;

      //  **********    gg initial state    ****************
      if ((me_w_info->m_fl1 ==21) && (me_w_info->m_fl2 ==21)){
        double l = log(mur2/muf2);
        correction = alphas * 2.*TR/(3.*M_PI) * l ;
        DEBUG_INFO("gg initial state. correction(" << correction <<").");
      }

      //  **********    qq initial state    ****************
      if ((me_w_info->m_fl1 !=21) && (me_w_info->m_fl2 !=21)){
        double l = log(mur2/mb2);
        correction = alphas * 2.*TR/(3.*M_PI) * l ;
        DEBUG_INFO("qq initial state. correction(" << correction <<").");
      }

      correction *= m_factor;
      
      double born_weight = me_w_info->m_B;
      double sum_meweight = me_w_info->m_B + me_w_info->m_K + me_w_info->m_KP + me_w_info->m_VI;
      new_weight_factor = (1. - correction* born_weight/sum_meweight);
    }
    else {
      msg_Debugging() << "H-Event, skip alpha_s correction." << std::endl;
      new_weight_factor = 1.0;
      // cannot return here, because each event needs to have the same weight structure!
    }

    DEBUG_VAR(new_weight_factor);
    Weights_Map& wmap = (*blobs->FindFirst(btp::Signal_Process))["WeightsMap"]->Get<Weights_Map>();
    wmap["Fusing"] = new_weight_factor;
    wmap["Fusing"]["Nominal"] = new_weight_factor;
    wmap["Fusing"]["NoDirectCorrection"] = 1.0;
    // TODO: calculate counter-terms based on the muR variations. not done yet, since numerical impact is small.

    return Return_Value::Nothing;
  }

  void Finish() {
  }

};

DECLARE_GETTER(Fusing_Direct_Hook,"Fusing_Direct",
               Userhook_Base,Userhook_Arguments);

Userhook_Base *ATOOLS::Getter<Userhook_Base,Userhook_Arguments,Fusing_Direct_Hook>::
operator()(const Userhook_Arguments &args) const
{
  return new Fusing_Direct_Hook(args);
}

void ATOOLS::Getter<Userhook_Base,Userhook_Arguments,Fusing_Direct_Hook>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"Fusing_Direct userhook";
}
