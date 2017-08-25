#include "PDF/Main/Cluster_Definitions_Base.H"
#include "PHASIC++/Process/Tree_ME2_Base.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PHASIC++/Process/Process_Group.H"
#include "PHASIC++/Process/Subprocess_Info.H"
#include "PHASIC++/Process/CS_Dipole.H"
#include "EXTAMP/Main/External_ME_Interface.H"
#include "EXTAMP/Main/Process_Group.H"
#include "EXTAMP/Main/Process.H"
#include "EXTAMP/Main/Born_Process.H"
#include "EXTAMP/Main/RS_Process.H"
#include "EXTAMP/Main/BVI_Process.H"

#include <assert.h>

namespace EXTAMP{


  External_ME_Interface::External_ME_Interface() : PHASIC::ME_Generator_Base("External") {}


  bool External_ME_Interface::Initialize(const std::string &path,const std::string &file,
					 MODEL::Model_Base *const model,
					 BEAM::Beam_Spectra_Handler *const beam,
					 PDF::ISR_Handler *const isr)
  {

    /* Pure virtual in PHASIC::ME_Generator_Base. Store beam, isr, and
       model since they need to be passed down to any
       PHASIC::Process_Base when calling its init method
       PHASIC::Process_Base::Init(...) */

    p_beam = beam;
    p_isr  = isr;
    return true;
  }

  External_ME_Interface::~External_ME_Interface() {}

  PHASIC::Process_Base* External_ME_Interface::InitializeProcess(const PHASIC::Process_Info &pi, bool add)
  {

    /* This is pure virtual in ME_Generator_Base and supposed to
       instantiate a high-level group-type process deriving from
       PHASIC::Process_Group. This inherits from PHASIC::Process_Base,
       which requires PHASIC::Process_Base::Init(...) to be called.

       Then, PHASIC::Process_Group::ConstructProcesses() is used in
       Sherpa to (recursively) initialize the PARTONIC sub-processes
       of the group. Within
       PHASIC::Process_Group::ConstructProcesses(), several pure
       virtuals of PHASIC::Process_Group are used for that purpose:
       
       - Process_Base* PHASIC::Process_Group::GetProcess(const PHASIC::Process_Info &pi)
       - bool          PHASIC::Process_Group::Initialize(Process_Base *const proc)
       
       The former must be implemented by EXTAMP::Process_Group in such
       a way as to instantiate the PARTONIC (i.e. non-group type)
       process corresponding to the Process_Info, and return a valid
       pointer to it. Does not make any sense structurally, so just
       implemented it as a wrapper around
       External_ME_Interface::InstantiatePartonicProcess(PHASIC::Process_Info). */
    
    EXTAMP::Process_Group *newproc = new EXTAMP::Process_Group();
    newproc->Init(pi,p_beam,p_isr);
    newproc->ConstructProcesses();
    newproc->SetGenerator(this);

    /* In case no valid partonic channels are found, return NULL so the
       Matrix_Element_Hanlder knows and throws a 'No hard process found' */
    return  (newproc->Size()>0) ? newproc : NULL;
  };


  /* Determine if a PARTONIC born Process as specified by pi exists */
  bool External_ME_Interface::PartonicProcessExists(PHASIC::Process_Info &pi)
  {
    if (PHASIC::Tree_ME2_Base::GetME2(pi)!=NULL) return true;

    /* Could be a real correction, increase order and try again. */
    if(pi.m_maxcpl[0]>pi.m_mincpl[0])
      pi.m_mincpl[0]+= 1;
    
    return PHASIC::Tree_ME2_Base::GetME2(pi)!=NULL;
  };


  PHASIC::Tree_ME2_Base* External_ME_Interface::GetExternalBornME(const PHASIC::Process_Info& pi)
  {
    /* Use the getter to find an implementation of born matrix element */
    PHASIC::Process_Info cpi(pi);
    
    PHASIC::Tree_ME2_Base* me2 = PHASIC::Tree_ME2_Base::GetME2(cpi);
    if(me2) return me2;

    /* Could be a real correction! */
    if(cpi.m_maxcpl[0]>cpi.m_mincpl[0])
      cpi.m_mincpl[0]+= 1;
    
    me2 = PHASIC::Tree_ME2_Base::GetME2(cpi);  
    if(me2) return me2;
    
    std::stringstream msg; msg << "No ME found for this process:\n" << pi;
    THROW(fatal_error, msg.str());
  }

  
  PHASIC::CS_Dipole* External_ME_Interface::GetExternalDipole(const PHASIC::Dipole_Info& di)
  {
    /* Use the getter to find implementation of correlated matrix element */
    PHASIC::CS_Dipole* me = PHASIC::CS_Dipole::Get(di);
    if(!me)
      THROW(fatal_error, "No dipole found");
    return me;
  }

  
  /* Instantiate a partonic (i.e. non-goup type) EXTAMP::Single_Process. */
  PHASIC::Process_Base* External_ME_Interface::InstantiatePartonicProcess
  (const PHASIC::Process_Info &pi)
  {
    ATOOLS::nlo_type::code nlotype=pi.m_fi.m_nlotype;
    
    /* Need to set orders according to our conventions before passing
       on: through previous call to PartonicProcessExists, it is
       ensured that a born-like ME with orders specified by mincpl
       exists. Can therefore set mincpl=maxcpl for non-loop
       contributions (born and real). For BVI, need to increase
       maxcpl[0] by one though. */
    PHASIC::Process_Info cpi(pi);
    cpi.m_maxcpl = cpi.m_mincpl;

    if( nlotype==ATOOLS::nlo_type::lo )
      return new Born_Process(cpi);
    
    if ( nlotype&ATOOLS::nlo_type::vsub )
      {
	cpi.m_maxcpl[0] += 1;
	return new BVI_Process(cpi);
      }
    
    if ( nlotype&ATOOLS::nlo_type::rsub )
      return new RS_Process(cpi);

    THROW(fatal_error, "Internal error");
    return NULL;
  };


  External_ME_Interface::Combinable_Map
  External_ME_Interface::ConstructCombinableMap(const ATOOLS::Flavour_Vector& flavs,
						const PHASIC::Process_Info&  pinfo,
						const size_t& nin)
  {
    /* Simple ad-hoc prescription good enough for reconstructing one
       splitting */
    DEBUG_FUNC(pinfo); assert(flavs.size()>0); Combinable_Map ret; 
    for(size_t i(0); i<flavs.size()-1; i++){
      for(size_t j=std::max(nin,i+1); j<flavs.size(); j++){

	/* Convert incoming flavour to outgoing flavour of opposite
	   charge (convention for clustering flavours in Sherpa) */
	const ATOOLS::Flavour& fl_i  = i<nin ? flavs[i].Bar() : flavs[i];
	const ATOOLS::Flavour& fl_j  = flavs[j];

	/* Combined flavour */
	ATOOLS::Flavour fl_ij;
	
	/* Combined id used to lookup this combination */
	size_t ij  = (1<<i)|(1<<j);
	
	/* g->gg or g->qqbar splitting */
	if(fl_i.Bar()==fl_j && (fl_i.IsQuark() || fl_i.IsGluon()))
	  fl_ij = ATOOLS::Flavour(kf_gluon);
	/* q->qg or q->gq splitting */
	if(fl_i.IsGluon() && fl_j.IsQuark())
	  fl_ij = fl_j;
	if(fl_j.IsGluon() && fl_i.IsQuark())
	  fl_ij = fl_i;

	/* No valid SM QCD splitting ij -> i,j */
	if(fl_ij.Kfcode()==kf_none)
	  continue;
	
	/* Check if this clustering yields a valid partonic process */
	PHASIC::Process_Info cpi(pinfo);
	/* Need to reverse the flavour of incoming particles yet again
	   because convention for PHASIC::Process_Info is different */
	cpi.Combine(i,j, i<nin ? fl_ij.Bar() : fl_ij); cpi.m_mincpl[0] -= 1;
	if(!External_ME_Interface::PartonicProcessExists(cpi)){
	  msg_Debugging() << "Discarding clustering "
			  << fl_i << "[" <<  i << "] + " << fl_j << "[" << j << "] "
			  << "since no born ME found for process\n " << cpi << "\n";
	  continue;
	}

	msg_Debugging() << "Inserging clustering " << fl_i << "[" <<  i << "] + " << fl_j << "[" << j << "]"<<"\n";
	ret.insert(std::make_pair(ij,ATOOLS::Flavour_Vector(1,fl_ij)));
      }
    }
    return ret;
  }

}

using namespace EXTAMP;

DECLARE_GETTER(External_ME_Interface,"External",
	       PHASIC::ME_Generator_Base,
	       PHASIC::ME_Generator_Key);

PHASIC::ME_Generator_Base *ATOOLS::Getter
<PHASIC::ME_Generator_Base,
 PHASIC::ME_Generator_Key,
 External_ME_Interface>::
operator()(const PHASIC::ME_Generator_Key &key) const
{
  return new External_ME_Interface();
}

void ATOOLS::Getter<PHASIC::ME_Generator_Base,
		    PHASIC::ME_Generator_Key,
		    External_ME_Interface>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"Generic ME interface"; 
}

