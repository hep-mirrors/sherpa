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
  Blob_List * p_bl;
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

    p_bl = blobs;
    String_BlobDataBase_Map  bdmap = p_bl->FindFirst(btp::Signal_Process)->GetData();

    auto search_mewinfo = bdmap.find("MEWeightInfo");
    if (search_mewinfo==bdmap.end()) {
        THROW(fatal_error,"No MEWeightinfo found in singnal blob!");
      }
    auto search_weight = bdmap.find("Weight");
    if (search_weight==bdmap.end()) {
        THROW(fatal_error,"No Weight found in singnal blob!");
      }

    ME_Weight_Info  * me_w_info = search_mewinfo->second->Get<ME_Weight_Info *>();
    double weight_bl = search_weight->second->Get<double>();
    double mur2 = me_w_info->m_mur2;
    ATOOLS::Flavour bquark = ATOOLS::Flavour(5);
    double mb2 = bquark.Mass()*bquark.Mass();
    double alphas((*p_as)(mur2));
    double TR = 0.5;


    double born_weight = me_w_info->m_B;


    double new_weight(weight_bl);
    double correction(0);
    double sum_meweight = me_w_info->m_B + me_w_info->m_K + me_w_info->m_KP + me_w_info->m_VI;


    //apply only to S-Events
    if ( me_w_info->m_type ==(ATOOLS::mewgttype::code::METS |ATOOLS::mewgttype::code::H)  ){
            msg_Debugging() << "H-Event, skip alpha_s correction." << std::endl;
            return Return_Value::Nothing;

        }

    PHASIC::Process_Vector procs = p_sherpa->GetInitHandler()->GetMatrixElementHandler()->AllProcesses();
    if (procs.size()!=1)
        THROW(fatal_error,"Too many procs!");
    PHASIC::MCatNLO_Process  * mcproc = dynamic_cast<PHASIC::MCatNLO_Process* >(procs.at(0));
    if (mcproc==NULL){
            THROW(fatal_error,"no MC@MLO process found! For use with separate LO-Process, use K-Factor!");

        }
    Cluster_Amplitude * ampl = mcproc->GetAmplitude();
    double muf2 = ampl->KT2();


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
    new_weight = weight_bl* (1. - correction* born_weight/sum_meweight);

    weight = new_weight;
    p_bl->FindFirst(btp::Signal_Process)->AddData("Weight",new Blob_Data<double>(new_weight));


	// TODO: calculate counter-terms based on the muR variations. not done yet, since numerical impact is small.
    auto search_varweights = bdmap.find("Variation_Weights");
    if (search_varweights==bdmap.end()) {
            THROW(fatal_error,"No VarWeight found in singnal blob!");
        }
    Variation_Weights  var_weights= search_varweights->second->Get<Variation_Weights >();
    var_weights*=(1. - correction* born_weight/sum_meweight);
    p_bl->FindFirst(btp::Signal_Process)->AddData("Variation_Weights",new Blob_Data<Variation_Weights>(var_weights));

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
