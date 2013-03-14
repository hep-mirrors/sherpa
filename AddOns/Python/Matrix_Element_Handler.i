//%module Matrix_Element_Handler
%{
#include <MODEL/Main/Model_Base.H>
#include <BEAM/Main/Beam_Spectra_Handler.H>
#include <PDF/Main/ISR_Handler.H>
#include <PHASIC++/Main/Phase_Space_Handler.H>
#include <PHASIC++/Process/Process_Base.H>
#include <PHASIC++/Process/ME_Generators.H>
#include <PHASIC++/Process/MCatNLO_Process.H>
#include <ATOOLS/Org/CXXFLAGS_PACKAGES.H>
#include <SHERPA/PerturbativePhysics/Matrix_Element_Handler.H>
#include <ATOOLS/Phys/Cluster_Amplitude.H>
#include <PHASIC++/Main/Process_Integrator.H>

#include <map> 

  typedef std::vector<PHASIC::NLOTypeStringProcessMap_Map*> ProcMap_Vector;
%}



namespace BEAM { class Beam_Spectra_Handler; }
namespace PDF {
  class ISR_Handler;
  class NLOMC_Base;
}

namespace PHASIC { 
  class Single_Process;
  class Selector_Key;
  struct nlo_type;
}

namespace SHERPA {

  class Shower_Handler;

  class Matrix_Element_Handler {
  public:

    typedef std::vector<PHASIC::NLOTypeStringProcessMap_Map*> ProcMap_Vector;

    typedef std::map<std::string,std::pair<int,int> >         MPIV_Map;
    typedef std::map<std::string,std::pair<int,double> >      MPDV_Map;
    typedef std::map<std::string,std::pair<int,std::string> > MPSV_Map;


  public :

    Matrix_Element_Handler(const std::string &path,const std::string &file,
                           const std::string &processfile,
                           const std::string &selectorfile);

    ~Matrix_Element_Handler();

    bool GenerateOneEvent(); 

    inline PHASIC::Process_Vector AllProcesses() const { return m_procs; }

    %extend{
       PHASIC::Process_Base* GetProcess(const ATOOLS::Cluster_Amplitude &ampl) const
      {   
	double weight(0.);
	std::string name(PHASIC::Process_Base::GenerateName(&ampl));
	for (int i(0); i<$self->ProcMaps().size(); i++){
	  PHASIC::StringProcess_Map::const_iterator pit($self->ProcMaps()[i]->find(PHASIC::nlo_type::lo)->second->find(name));
	  if(pit == $self->ProcMaps()[i]->find(PHASIC::nlo_type::lo)->second->end())
	    continue;
	  else
	    return pit->second;	
	}
	std::string error("Could not find process ");
	error.append(name);
	throw  error.c_str();
      }

       const double GetWeight(const ATOOLS::Cluster_Amplitude &ampl) const
      {   
	double weight(0.);
	std::string name(PHASIC::Process_Base::GenerateName(&ampl));
	for (int i(0); i<$self->ProcMaps().size(); i++){
	  PHASIC::StringProcess_Map::const_iterator pit($self->ProcMaps()[i]->find(PHASIC::nlo_type::lo)->second->find(name));
	  if(pit == $self->ProcMaps()[i]->find(PHASIC::nlo_type::lo)->second->end())
	    continue;
	  else{
	    SP(PHASIC::Color_Integrator) CI = (pit->second->Integrator()->ColorIntegrator());
	    if(CI!=0)
	      {
                CI->GeneratePoint();
                for (size_t i=0; i<ampl.Legs().size(); ++i)
                ampl.Leg(i)->SetCol(ATOOLS::ColorID(CI->I()[i],CI->J()[i]));
               }
	    weight+=pit->second->Differential(ampl);
	    return weight; 
	  }
	}
	std::string error("Could not find process ");
	error.append(name);
	throw  error.c_str();
      }
    }// end of %extend

  };// end of class Matrix_Element_Handler

}// end of namespace SHERPA
