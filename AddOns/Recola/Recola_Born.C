#include "AddOns/Recola/Recola_Born.H"
#include "AddOns/Recola/Recola_Interface.H"
#include "PHASIC++/Process/External_ME_Args.H"
#include "ATOOLS/Org/Run_Parameter.H" 
#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Library_Loader.H"


using namespace PHASIC;

namespace Recola {


  Recola_Born::Recola_Born(const External_ME_Args& args,
			   unsigned int recola_id, int amptype) :

    Tree_ME2_Base(args), 
    m_recola_id(recola_id), 
    m_amptype(amptype)

  {
    m_symfac =Flavour::FSSymmetryFactor(args.m_outflavs);
    m_symfac*=Flavour::ISSymmetryFactor(args.m_inflavs);
  }


  double Recola_Born::Calc(const Vec4D_Vector& momenta) 
  {
    if (!Recola_Interface::checkProcGeneration()){
        Recola_Interface::GenerateProcesses(AlphaQED(),AlphaQCD(),
                                            100.,100.,100.);
    }
     
    MyTiming* timing;
    if (msg_LevelIsDebugging()) {
      timing = new MyTiming();
      timing->Start();
    }
    double aqcd=AlphaQCD();
    // TODO: where to get this from?
    double mur=100.;
    int defflav=Recola_Interface::GetDefaultFlav();
    set_alphas_rcl(aqcd,mur,defflav);

    
    double res(0.0);
    if (m_amptype==12||m_amptype==1) 
    {
	    double born;
	    Recola_Interface::EvaluateBorn(m_recola_id, momenta, born, m_amptype);
	    res = born;
    }
    else
    {
	    THROW(not_implemented, "Unknown amplitude type");
    }
   
    // Recola returns ME2 including 1/symfac, but Calc is supposed to return it
    // without 1/symfac, thus multiplying with symfac here
    return res*m_symfac;
  } 
  
  bool Recola_Born::IsMappableTo(const PHASIC::Process_Info& pi){
    return false;
  }
}

using namespace Recola;

DECLARE_TREEME2_GETTER(Recola_Born,
		       "Recola_Born")

Tree_ME2_Base *ATOOLS::Getter<Tree_ME2_Base,
			      External_ME_Args,
			      Recola_Born>::
operator()(const External_ME_Args& args) const
{
  if(args.m_source.length() &&
     args.m_source != "Recola") return NULL;
  PRINT_VAR("IN Recola born getter");
  int amptype = Recola_Interface::GetAmpType();
  int id = Recola_Interface::RegisterProcess(args, amptype);
  if (id<=0) return NULL;
 
  return new Recola_Born(args,id,amptype);
}

