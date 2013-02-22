//%module Process_Base
%include "std_string.i"
%{
#include <PHASIC++/Process/Process_Info.H>
#include <PHASIC++/Selectors/Selector.H>
#include <PHASIC++/Selectors/Cut_Data.H>
#include <PHASIC++/Scales/Scale_Setter_Base.H>
#include <PHASIC++/Scales/KFactor_Setter_Base.H>
#include <ATOOLS/Phys/NLO_Subevt.H>
#include <PHASIC++/Process/Process_Base.H>
#include <PHASIC++/Main/Process_Integrator.H>
#include <ATOOLS/Phys/Cluster_Leg.H>

  using namespace PHASIC;
  using namespace ATOOLS;
  %}

namespace ATOOLS { 
  class Cluster_Leg;
  class Cluster_Amplitude; 
  class Histogram;
  struct Decay_Info;
  typedef std::vector<Decay_Info* > DecayInfo_Vector;
}

namespace BEAM { class Beam_Spectra_Handler; }

namespace PDF { 
  class ISR_Handler;
  class Shower_Base;
}

namespace METOOLS { class Spin_Amplitudes; }

namespace PHASIC {

  class Process_Integrator;
  class Phase_Space_Handler;
  class Combined_Selector;
  class ME_Generator_Base;

  struct Weight_Info;

  class Process_Base;

  class Process_Base {
  protected:
    
  public:

    Process_Base();
    virtual ~Process_Base();

    virtual double Differential(const ATOOLS::Vec4D_Vector &p) = 0;
    virtual double Differential(const ATOOLS::Cluster_Amplitude &ampl,
                                int mode=0);

    static void SortFlavours(ATOOLS::Cluster_Amplitude *const ampl);

    static std::string GenerateName(const ATOOLS::Cluster_Amplitude *ampl);
    static std::string GenerateName(const ATOOLS::NLO_subevt *sub,
				    const size_t &nin);

    static std::string GenerateName(const Subprocess_Info &info);
    static std::string GenerateName(const Subprocess_Info &ii,
				    const Subprocess_Info &fi);

    inline size_t NIn() const  { return m_nin;  }
    inline size_t NOut() const { return m_nout; }

    inline const ATOOLS::Flavour_Vector &Flavours() const { return m_flavs; }

    inline const std::string &Name() const { return m_name; }

    %extend{
      double GeneratePoint(ATOOLS::Cluster_Amplitude &amp)
      {
	PHASIC::Color_Integrator *CI(&*$self->Integrator()->ColorIntegrator());
	if (!CI)
	  THROW(fatal_error, "No color integrator. Make sure Comix is used when calling this function")
	CI->GeneratePoint();
	for (size_t i=0; i<amp.Legs().size(); ++i)
	  amp.Leg(i)->SetCol(ATOOLS::ColorID(CI->I()[i],CI->J()[i]));
	return CI->GlobalWeight();
      }
    }

    %extend{
      void SetColors(ATOOLS::Cluster_Amplitude &amp)
      { 
	PHASIC::Int_Vector ci(amp.Legs().size());
	PHASIC::Int_Vector cj(amp.Legs().size());
	PHASIC::Color_Integrator *CI(&*$self->Integrator()->ColorIntegrator());
	if (!CI)
	  THROW(fatal_error, "No color integrator. Make sure Comix is used when calling this function")
	CI->GeneratePoint();
	for (size_t i=0; i<amp.Legs().size(); ++i)
	  {
	    ci[i] = amp.Leg(i)->Col().m_i;
	    cj[i] = amp.Leg(i)->Col().m_j;
	  }
	CI->SetI(ci);
	CI->SetJ(cj);
      }
    }

  };// end of class Process_Base

}// end of namespace ATOOLS
