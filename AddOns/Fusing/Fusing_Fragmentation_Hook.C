#include "SHERPA/Tools/Userhook_Base.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"

#include "SHERPA/Single_Events/Event_Handler.H"

#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Phys/Weight_Info.H"
#include <algorithm>
#include "MODEL/Main/Running_AlphaS.H"
#include "SHERPA/Main/Sherpa.H"
#include "SHERPA/Initialization/Initialization_Handler.H"
#include "ATOOLS/Org/Data_Reader.H"


/*
 * implements the fragmentation component of the fusing algorithm based
 * on a combined cluster amplitude.
*/


using namespace ATOOLS;
using namespace SHERPA;



class Fusing_Fragmentation_Hook : public Userhook_Base {

private:
  MODEL::Running_AlphaS *p_as;
  Blob_List * p_bl;
  Sherpa* p_sherpa;
  size_t m_tops_in_proc;
  bool m_nlo, m_check_ho;
  bool m_store;



public:

  Fusing_Fragmentation_Hook(const Userhook_Arguments args) :
    Userhook_Base("Fusing_Fragmentation"), p_as(MODEL::as), p_sherpa(args.p_sherpa)
  {
       msg_Debugging()<<"Fusing_Fragmentation Hook active."<<std::endl;
       ATOOLS::Data_Reader *reader =  p_sherpa->GetInitHandler()->DataReader();

       m_nlo = reader->GetValue<int>("FUSING_FRAGMENTATION_NLO_Mode", 1);
       m_store = reader->GetValue<int>("FUSING_FRAGMENTATION_STORE_AS_WEIGHT", 0);
       m_check_ho = reader->GetValue<int>("FUSING_FRAGMENTATION_CHECK_HO", 1);

       msg_Debugging()<<METHOD<< "NLO(" << m_nlo << ") \n }" << std::endl;

	   // count number of tops for ttbb in case of rare g->tt cluster configurations	
       m_tops_in_proc=0; 
       PHASIC::Process_Base  * proc = p_sherpa->GetInitHandler()->GetMatrixElementHandler()->AllProcesses().at(0);
       ATOOLS::Flavour_Vector flavs = proc->Flavours();
       for (int i=2; i< flavs.size(); i++){
           ATOOLS::Flavour flav = flavs.at(i);
           if  (flav.Kfcode()==6) m_tops_in_proc ++;
         }


 }

  ~Fusing_Fragmentation_Hook() {}

  ATOOLS::Return_Value::code Run(ATOOLS::Blob_List* blobs, double &weight) {
    p_bl = blobs;

       ATOOLS::String_BlobDataBase_Map  bdmap = p_bl->FindFirst(ATOOLS::btp::Shower)->GetData();

       auto search = bdmap.find("AllAmplitudes");
       if (search==bdmap.end()) {
           THROW(fatal_error,"No matching amplitude found in blob. This algorithm works only with rel-2-2-7 or later!");
         }

       ATOOLS::Cluster_Amplitude * ampl = search->second->Get<ATOOLS::Cluster_Amplitude*>();
       bool veto   = true;

       //find first amplitude with a bb-pair
       size_t pos=0;
       while (true){
           if (ampl==NULL) break;
           pos = CheckBPair(ampl);
           if (pos==0){
                   ampl = ampl->Next();
           } else break;
         }

       if (pos==0) veto=false; // no bb-pair found
       else veto = !CheckFinalDir(ampl);

       // redo veto for certain configurations (see code below)
       if (veto && m_check_ho) {
            if (CheckHigherOrder(ampl)){
                msg_Debugging() << "skip Fusing Veto, configuration is of higher order.";
                veto=false;
            }
       }

       return ReturnFunc(!veto, weight);




  }

  Return_Value::code ReturnFunc(bool retval, double & event_weight){
      //retval=false: veto Event
      //retval=true:  keep Event

      String_BlobDataBase_Map  bdmap = p_bl->FindFirst(btp::Signal_Process)->GetData();
      auto search_weight = bdmap.find("Weight");
      if (search_weight==bdmap.end()) {
              THROW(fatal_error,"No Weight found in singnal blob!");
          }
      double weight_bl = search_weight->second->Get<double>();




      if (retval==false){
              msg_Debugging() << "false, veto configuration.  \n";
          }

      // if not storing weights:  do either Veto or not
      if (!m_store) {
              if (retval) return Return_Value::Nothing;
              else return Return_Value::New_Event;
          }
      // if m_store: store additional weight(zero) in bloblist and slip veto
      double new_weight = weight_bl;
      if (retval==false) new_weight=0;

      p_bl->FindFirst(ATOOLS::btp::Signal_Process)->AddData("UserHook",new ATOOLS::Blob_Data<double>(new_weight));
      return Return_Value::Nothing;


  }

  ///////////////////////   helper functions Veto


      bool CheckFinalDir(ATOOLS::Cluster_Amplitude *ampl){
        /*return only true, if a sufficient high number of light partons is present in this amplitude's final state
          if the tops are clustered to a gluon before, there is one more light parton allowed in this step */

        //  true:   keep amplitude
        //  false:  do veto
        size_t num_q(0);
        size_t num_g(0);
        size_t num_tops(0);
        for (int i=2; i<ampl->Legs().size(); i++){
            ATOOLS::Cluster_Leg *leg = ampl->Legs().at(i);
            if (leg->Flav().IsGluon() && !leg->FromDec()) num_g++;
            if (abs(leg->Flav().Kfcode())<5 && !leg->FromDec()) num_q++;
            if (abs(leg->Flav().Kfcode())==6 && !leg->FromDec()) num_tops++;
        }


        if (m_nlo){
            if (num_tops < m_tops_in_proc){
                if ((num_q + num_g)>2) return true;
                else return false;
              }
            else{
                if ( (num_q + num_g)>1) return true;
                else return false;
             }
         }
        else{
            if (num_tops < m_tops_in_proc){
                if ((num_q >0) || (num_g>1)) return true;
                else return false;
              }
            else{
              if ((num_q >0) || (num_g>0)) return true;
              else return false;
            }
         }
      }



      size_t CheckBPair(ATOOLS::Cluster_Amplitude *ampl){
        if (ampl==NULL) return 0;
        bool bfound(false), bbfound(false);
        for (int i=2; i<ampl->Legs().size(); i++){
            ATOOLS::Cluster_Leg *leg = ampl->Legs().at(i);
            if ( !leg->FromDec() && leg->Flav().Kfcode()==5){
              if (leg->Flav().IsAnti()) {
                bbfound=true;
              }else{
                 bfound=true;
                }
              }
          }
        if (bfound && bbfound) return 1;
        return 0;
      }


      bool CheckHigherOrder (ATOOLS::Cluster_Amplitude *ampl){
        /*check for configurations as bb->ttbb or bb->bb which should not be vetoed,
         * since this is not part of 4F ttbb.
         implemantation: redo veto, if bb-> x bb (j) and in x is no particle of the 93 container.

         returns true, if such a config is found
        */
        if (!CheckBPair(ampl)) {
            THROW(fatal_error,"amplitude is missing b-quarks!");
          }
        msg_Debugging() << "Fusing, check higher order: " << ampl;
        for (int i =0; i < 2; i++ ){
          ATOOLS::Flavour flav = ampl->Leg(i)->Flav();
          if (flav.Kfcode()!=5) return false;
        }
        return true;
      }




  void Finish() {
  }



};

DECLARE_GETTER(Fusing_Fragmentation_Hook,"Fusing_Fragmentation",
               Userhook_Base,Userhook_Arguments);

Userhook_Base *ATOOLS::Getter<Userhook_Base,Userhook_Arguments,Fusing_Fragmentation_Hook>::
operator()(const Userhook_Arguments &args) const
{
  return new Fusing_Fragmentation_Hook(args);
}

void ATOOLS::Getter<Userhook_Base,Userhook_Arguments,Fusing_Fragmentation_Hook>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"Fusing_Fragmentation_UserHook";
}
