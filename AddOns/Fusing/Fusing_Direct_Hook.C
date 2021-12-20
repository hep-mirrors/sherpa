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
    Userhook_Base("Example"), p_as(MODEL::as), p_sherpa(args.p_sherpa)
  {
    msg_Debugging()<<"Fusing_Direct Hook active."<<std::endl;
    ATOOLS::Data_Reader *reader =  p_sherpa->GetInitHandler()->DataReader();
    m_factor = reader->GetValue<double>("FUSING_DIRECT_FACTOR", 2.);
  }

  ~Fusing_Direct_Hook() {}

  ATOOLS::Return_Value::code Run(ATOOLS::Blob_List* blobs, double &weight) {
    String_BlobDataBase_Map  bdmap = blobs->FindFirst(btp::Signal_Process)->GetData();
    String_BlobDataBase_Map::iterator search_mewinfo(bdmap.find("MEWeightInfo")), search_weight(bdmap.find("Weight"));
    if (search_mewinfo==bdmap.end() || search_weight==bdmap.end()) {
      THROW(fatal_error,"Internal error: No weight info found in signal blob!");
    }

    ME_Weight_Info  * me_w_info = search_mewinfo->second->Get<ME_Weight_Info *>();
    double mur2 = me_w_info->m_mur2;
    ATOOLS::Flavour bquark = ATOOLS::Flavour(5);
    double mb2 = bquark.Mass()*bquark.Mass();
    double alphas((*p_as)(mur2));
    double TR = 0.5;
    double born_weight = me_w_info->m_B;
    double correction(0);
    double sum_meweight = me_w_info->m_B + me_w_info->m_K + me_w_info->m_KP + me_w_info->m_VI;

    //apply only to S-Events
    if ( me_w_info->m_type ==(ATOOLS::mewgttype::METS |ATOOLS::mewgttype::H)  ){
      msg_Debugging() << "H-Event, skip alpha_s correction." << std::endl;
      return Return_Value::Nothing;
    }

    PHASIC::MCatNLO_Process * mcproc = dynamic_cast<PHASIC::MCatNLO_Process* >
      (p_sherpa->GetInitHandler()->GetMatrixElementHandler()->Process()->Parent());
    if (mcproc==NULL){
      THROW(fatal_error,"no MC@NLO process found! For use with separate LO-Process, use K-Factor!");
    }
    // TODO: make more efficient, without whole amplitude copying!
    Cluster_Amplitude * ampl = mcproc->GetAmplitude();
    double muf2 = ampl->KT2();
    ampl->Delete();


    //  **********    gg initial state    ****************
    if ((me_w_info->m_fl1 ==21) && (me_w_info->m_fl2 ==21)){
      double l = log(mur2/muf2);
      correction = alphas * 2.*TR/(3.*M_PI) * l ;
      msg_Debugging() << "gg initial state. correction(" << correction <<")." << std::endl;
    }

    //  **********    qq initial state    ****************
    if ((me_w_info->m_fl1 !=21) && (me_w_info->m_fl2 !=21)){
      double l = log(mur2/mb2);
      correction = alphas * 2.*TR/(3.*M_PI) * l ;
      msg_Debugging() << "qq initial state. correction(" << correction <<")." << std::endl;
    }


    correction *= m_factor;

    double new_weight = search_weight->second->Get<double>() * (1. - correction* born_weight/sum_meweight);
    (*blobs->FindFirst(btp::Signal_Process))["Weight"]->Set(new_weight);
    weight = new_weight; // obsolete?

    // TODO: calculate counter-terms based on the muR variations. not done yet, since numerical impact is small.
    String_BlobDataBase_Map::iterator search_varweights = bdmap.find("Variation_Weights");
    if (search_varweights==bdmap.end()) {
      THROW(fatal_error,"No VarWeight found in signal blob!");
    }
    Variation_Weights  var_weights= search_varweights->second->Get<Variation_Weights >();
    var_weights*=(1. - correction* born_weight/sum_meweight);
    (*blobs->FindFirst(btp::Signal_Process))["Variation_Weights"]->Set(var_weights);

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
