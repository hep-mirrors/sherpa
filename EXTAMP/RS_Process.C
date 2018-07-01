#include "EXTAMP/RS_Process.H"
#include "EXTAMP/CS_Dipoles.H"
#include "EXTAMP/Dipole_Wrapper_Process.H"
#include "EXTAMP/External_ME_Interface.H"

#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Phys/NLO_Subevt.H"

#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Process/Tree_ME2_Base.H"
#include "PHASIC++/Process/External_ME_Args.H"
#include "PHASIC++/Process/Spin_Color_Correlated_ME2.H"
#include "PHASIC++/Selectors/Combined_Selector.H"
#include "YODA/WriterYODA.h"

#include <assert.h>

namespace EXTAMP {

  RS_Process::RS_Process(const PHASIC::Process_Info& pi) : Process(pi) {
    DEBUG_FUNC(pi);

    std::vector<double> orders = pi.m_borncpl; orders[0] += 1;
    PHASIC::External_ME_Args args(pi.m_ii.GetExternal(),
				  pi.m_fi.GetExternal(),
				  orders);
    p_real_me = External_ME_Interface::GetExternalBornME(args);
    p_real_me->SetCouplings(m_cpls);

    /* Construct dipoles, these do not require any input beyond flavour config */
    m_dipoles   = ConstructDipoles();

    /* Create subevents. Not all properties set here, will be done
       later by Dipole_Wrapper_Processes */
    m_subevents = ConstructSubevents(m_dipoles.size());

    /* Create Coupling_Data instances for each dipole and associate
       with corresponding NLO_subevts */
    ConstructRunningCouplings(m_cpls,m_subevents,m_dipoles);

    /* Summary of constructed dipoles to debugging output */
    msg_Debugging() << "Constructed " << m_dipoles.size() << " dipoles:\n";
    for(std::vector<CS_Dipole*>::const_iterator it=m_dipoles.begin();
	it!=m_dipoles.end(); ++it)
      msg_Debugging() << (*it)->Info() << "\n";

    /* Get parameters for subevent smearing */
    ATOOLS::Default_Reader reader;
    m_alpha_0     = reader.Get<double>("NLO_SMEAR_THRESHOLD", 0.0);
    m_smear_power = reader.Get<double>("NLO_SMEAR_POWER",     0.5);

    /* In AMEGIC, m_alpha_0>0.0 is used for a functional form different from alpha. */
    if(m_alpha_0>0.0) THROW(not_implemented, "Smearing only implemented for alpha parameter");
    
//    int n_bin = 100;
//        std::vector<double> bins;
//        for(int i=0; i<=n_bin;i++) { bins.push_back(pow(10,10.*(1.*i/n_bin - 1))); }
//    m_hist_mreal_y46 = YODA::Histo1D(bins, "RSterm/y46", "Mreal_y46");
//    m_hist_dipole_y46 = YODA::Histo1D(bins, "RSterm/y46", "dipole_y46");
//    m_hist_alpha_min = YODA::Histo1D(bins, "alpha_min", "alpha_min");
    m_myfile.open ("RSdata.dat");
  }


  void RS_Process::Init(const PHASIC::Process_Info &pi,
			BEAM::Beam_Spectra_Handler *const beamhandler,
			PDF::ISR_Handler *const isrhandler,const int mode)
  {
    PHASIC::Process_Base::Init(pi, beamhandler, isrhandler, mode);

    /* Create Dipole wrapper procs, PHASIC::Process_Base::Init() needs
       to be called first */
    for(size_t i(0); i<m_dipoles.size(); i++)
      m_dipole_wrappers.
	push_back(new Dipole_Wrapper_Process(*this,m_dipoles[i],
					     Integrator()->Beam(),
					     Integrator()->ISR()));

    /* Now set properties of NLO_subevts properly. This is done by the
       dipole wrapper procs since their flavour ordering must be used
       in the NLO_subevts */
    assert(m_subevents.size() == m_dipole_wrappers.size()+1);
    for(size_t i(0); i<m_dipoles.size(); i++)
      {
	m_dipole_wrappers[i]->SetSubEventProperties(*m_subevents[i]);
	m_dipole_wrappers[i]->AssignSubEvent(m_subevents[i]);
      }


    /* Can set name of subevent only after Process_Base::Init has been
       called */
    ATOOLS::NLO_subevt* evt = m_subevents.back();
    evt->m_pname  =  Name();
    evt->m_pname  =  evt->m_pname.substr(0,evt->m_pname.rfind("__"));

    /* With the integrator available, can now set momentum pointer of
       subevt */
    evt->p_mom    = &(Integrator()->Momenta()[0]);

  }


  double RS_Process::Partonic(const ATOOLS::Vec4D_Vector &p,
			      const int mode)
  {
    /* Calculate dipole kinematics and update subevents accordingly */
    CalculateKinematics(p);
//DEBUG_VAR(p);

    /* for ID */
//    const ATOOLS::Vec4D n = p[0]+p[1]-p[3]-p[5]-p[6];
//    m_hist_alpha_min.fill( p[4]*p[5] / (p[5]*n), 1.);  // vi_tilde for b-quark
    /* for CS */
//    m_hist_alpha_min.fill( p[4]*p[5] / (p[4]*p[5]+p[4]*p[6]+p[5]*p[6]), 1.);  // y for b-quark
    if(!PassesAlphaMin(m_dipoles, p))
      {
	SetSubEventsToZero(m_subevents);
	return m_lastxs = 0.0;
      }

    /* Check trigger and set m_trig of subevents accordingly */
    Selector()->RSTrigger(&m_subevents);
    
    /* Trigger the calculation scales (including those of dipoles) */
    ScaleSetter()->CalculateScale(p);

    /* Calculate sum of dipole subtraction terms */
    double S(0.0);
    for(size_t i=0; i<m_dipoles.size(); i++)
      {
	/* Check if dipole kinematic passes cuts before evaluating.
	   Dipoles approximate +real emission term. Sign convention in
	   Dipole_Wrapper_Processes is such that ADDING their
	   contribution cancels divergencies. */
	bool sub_trig   = m_subevents[i]->m_trig;
//DEBUG_VAR(sub_trig);
	double sub_dxs = (sub_trig ? m_dipole_wrappers[i]->Calc(m_subevents[i]) : 0.0);
	S += sub_dxs;
      if (msg_LevelIsDebugging()) {
        std::cout << "dipole[" << i << "] = " << sub_dxs << "\n"; 
      }
      }

    /* Check if kinematics passes trigger before calculating real
       emission matrix element */
    bool trig = m_subevents.back()->m_trig;
    double R  = (trig ? p_real_me->Calc(p)/NormFac() : 0.0);

    /* Update real emission subevent */
    m_subevents.back()->m_trig   = trig;
    m_subevents.back()->m_me     = R;
    m_subevents.back()->m_mewgt  = R;
    m_subevents.back()->m_result = R;

    /* Apply smearing to reduce binning fluctuations */
    if(m_alpha_0!=0.0) SmearSubEvents(m_dipoles, m_subevents, ATOOLS::dabs(m_alpha_0), m_smear_power);
    DEBUG_VAR(R);
    DEBUG_VAR(S);
//    m_hist_dipole_y46.fill(p[4]*p[6]/(p[1]*p[0]), S);
//    m_hist_mreal_y46.fill( p[4]*p[6]/(p[1]*p[0]), R);

    if( R!=0 && S!=0){
    m_myfile << std::setprecision(9) << p[4]*p[5]/(p[1]*p[0]) << " " 
             << p[4][0] << " "
             << R << " " << S << std::endl;
    }

    return m_lastxs = R + S;
  }


  void RS_Process::CalculateKinematics(const ATOOLS::Vec4D_Vector& p)
  {
    /* Calculate Born-like kinematics for all dipoles and set the
       momenta of corresponding NLO_subevents (used by scale setter,
       selector, etc ... ) */
    for(size_t i(0); i<m_dipoles.size(); i++)
	m_dipole_wrappers[i]->CalcKinematics(p);
  }


  bool RS_Process::PassesAlphaMin(const Dipole_Vector& dv, const ATOOLS::Vec4D_Vector& p) const
  {
    for(Dipole_Vector::const_iterator it=dv.begin(); it!=dv.end(); ++it)
      if(! ((*it)->PassesAlphaMin(p)))
	 return false;
    return true;
  }


  void RS_Process::SetSubEventsToZero(ATOOLS::NLO_subevtlist subs) const
  {
    for(ATOOLS::NLO_subevtlist::iterator it=subs.begin();it!=subs.end();++it)
      (*it)->m_trig = (*it)->m_me = (*it)->m_result = (*it)->m_mewgt = 0;
  }


  bool RS_Process::ApplyIDdipole(const Dipole_Info &di) const
  {
    /* return true, if both real config. and Born-config. have (at least) one b- and 
       one bbar-quark */
//    /* AND if there is no g->gg splitting*/
//
//    /* check whether there is g->gg splitting */
//    if (di.m_flav_type == FlavourType::gtogg) return false;
    
    /* count b&bbar quarks in final state of real config,
       require, to have at least one b-quark and
       one bbar-quark*/
    int cnt_b=0, cnt_bbar=0;
    for(size_t cnt=2; cnt<di.m_real_flavs.size();cnt++){
      if(di.m_real_flavs[cnt].IsbQuark())    cnt_b++;
      if(di.m_real_flavs[cnt].IsbbarQuark()) cnt_bbar++;
    }
    if(cnt_b<1 || cnt_bbar<1) return false;
  
    /* initialisation of auxiliary dipoles, necessary to have access to Born flavours */
    CS_Dipole *dipole;
    switch(di.m_split_type){
    case SplittingType::FF:
      dipole = new FF_Dipole(di); break;
    case SplittingType::FI:
      dipole = new FI_Dipole(di); break;
    case SplittingType::IF:
      dipole = new IF_Dipole(di); break;
    case SplittingType::II:
      dipole = new II_Dipole(di); break;
    default:
      THROW(fatal_error, "Internal error");
    }

    /* count b&bbar quarks in final state of Born-config,
       require, to have at least one b-quark and
       one bbar-quark*/
    cnt_b=0, cnt_bbar=0;
    for(size_t cnt=2; cnt<dipole->Flavours().size();cnt++){
      if(dipole->Flavours()[cnt].IsbQuark())    cnt_b++;
      if(dipole->Flavours()[cnt].IsbbarQuark()) cnt_bbar++;
    }
    if(cnt_b<1 || cnt_bbar<1) return false;

    return true;
  }
  

  RS_Process::Dipole_Vector RS_Process::ConstructDipoles()
  {
    /* Get subtraction parameters */
    ATOOLS::Default_Reader reader;
    int subtraction_type = ATOOLS::ToType<int>(ATOOLS::rpa->gen.Variable("NLO_SUBTRACTION_SCHEME"));
    double alphamin = reader.Get<double>("DIPOLE_AMIN" , ATOOLS::Max(ATOOLS::Accu(),1.0e-8));
    double alphamax = reader.Get<double>("DIPOLE_ALPHA", 1.0);

    /* Build dipoles */
    Dipole_Vector ret;
    for(std::vector<size_t>::const_iterator j=PartonIndices().begin(); j!=PartonIndices().end(); ++j)
      for(std::vector<size_t>::const_iterator i=j+1; i!=PartonIndices().end(); ++i)
	for(std::vector<size_t>::const_iterator k=PartonIndices().begin(); k!=PartonIndices().end(); ++k)
	  {
	    if(i==k || j==k) continue;
	    if(!Combinable((1<<(*i)), (1<<(*j)))) continue;
        Dipole_Info di(m_flavs, *i, *j, *k, subtraction_type, alphamin, alphamax, m_dipole_case);

        /* for pp>WWjj, pseudo-dipoles are not always used */
        if(m_dipole_case == DipoleCase::IDin && !ApplyIDdipole(di) )
	      di = Dipole_Info(m_flavs, *i, *j, *k, subtraction_type, alphamin, alphamax,
               DipoleCase::CS);

	    switch(di.m_split_type)
	      {
	      case SplittingType::FF:
		ret.push_back(new FF_Dipole(di)); break;
	      case SplittingType::FI:
		ret.push_back(new FI_Dipole(di)); break;
	      case SplittingType::IF:
		ret.push_back(new IF_Dipole(di)); break;
	      case SplittingType::II:
		ret.push_back(new II_Dipole(di)); break;
	      default:
		THROW(fatal_error, "Internal error");
	    }
	  }
    return ret;
  }


  ATOOLS::NLO_subevtlist RS_Process::ConstructSubevents(size_t n_dipoles) const
  {
    /* Create one NLO_subevt per CS dipole. Properties will be set
       later by Dipole_Wrapper_Processes */
    ATOOLS::NLO_subevtlist ret;
    for(size_t i(0); i<n_dipoles; i++)
      ret.push_back(new ATOOLS::NLO_subevt());

    /* Last item in list corresponds to real emission configuration
       per convention, with i=j=k=0 */
    ATOOLS::NLO_subevt* evt = new ATOOLS::NLO_subevt();
    evt->p_fl     = &Flavours()[0];

    /* Can only set this pointer to it's proper address once the
       integrator of this proc is initialized. Do this in the Init
       method */
    evt->p_mom    =  NULL;
    
    evt->m_n      =  Flavours().size();
    evt->m_i      =  0;
    evt->m_j      =  0;
    evt->m_k      =  0;
    evt->m_me     =  0.0;
    evt->m_result =  0.0;
    evt->m_trig   =  false;
    evt->p_proc   =  (void*)this;
    ret.push_back(evt);

    /* For all subevents, set  a pointer to the real emission subevent */
    for(ATOOLS::NLO_subevtlist::const_iterator it=ret.begin(); it!=ret.end(); ++it)
      (*it) -> p_real = ret.back();

    return ret;
  }


  void RS_Process::SmearSubEvents(const Dipole_Vector& dipoles,
				  ATOOLS::NLO_subevtlist& subs,
				  const double& alpha_0,
				  const double& power)
  {
    assert(subs.size() == dipoles.size() +1);
    assert(alpha_0>0.0);

    ATOOLS::NLO_subevt* realevt = subs.back();
      
    for (size_t i=0;i<dipoles.size();i++) {

      double alpha = dipoles[i]->LastKinematics()->Alpha();
      double x     = pow(alpha/alpha_0, power);

      if (alpha>alpha_0) continue;

      realevt->m_me     += (1.0-x)*subs[i]->m_me;
      realevt->m_mewgt  += (1.0-x)*subs[i]->m_mewgt;
      realevt->m_result += (1.0-x)*subs[i]->m_result;

      subs[i]->m_me     *= x;
      subs[i]->m_mewgt  *= x;
      subs[i]->m_result *= x;
    }
  }
  

  void RS_Process::ConstructRunningCouplings(MODEL::Coupling_Map& cpls,
					     const ATOOLS::NLO_subevtlist& evts,
					     const Dipole_Vector& dipoles) const
  {
    assert(dipoles.size() == evts.size()-1);

    MODEL::Coupling_Data* rqcd(cpls.Get("Alpha_QCD"));
    MODEL::Coupling_Data* rqed(cpls.Get("Alpha_QED"));

    if(!rqcd) THROW(fatal_error, "Invalid pointer");
    if(!rqed) THROW(fatal_error, "Invalid pointer");

    /* Create copies of the Coupling_Data instances, associate
       current NLO_subevt to them, and store in multimap */
    for(size_t i(0); i<dipoles.size(); i++)
      {
	MODEL::Coupling_Data* nrqcd = new MODEL::Coupling_Data(*rqcd,evts[i]);
	MODEL::Coupling_Data* nrqed = new MODEL::Coupling_Data(*rqed,evts[i]);
	cpls.insert(std::make_pair("Alpha_QCD",nrqcd));
	cpls.insert(std::make_pair("Alpha_QED",nrqed));

	/* Now tell the dipole where to get its (pre-calculated)
	   running couplings */
	dipoles[i]->CorrelatedME()->SetCouplings(nrqcd, nrqed);
      }
  }


  void RS_Process::SetScale(const PHASIC::Scale_Setter_Arguments &args)
  {
    PHASIC::Single_Process::SetScale(args);
    for(Dipole_Wrappers::const_iterator it=m_dipole_wrappers.begin();
	it!=m_dipole_wrappers.end(); ++it) (*it)->SetScaleSetter(p_scale);
  }


  void RS_Process::SetKFactor(const PHASIC::KFactor_Setter_Arguments &args)
  {
    PHASIC::Single_Process::SetKFactor(args);
    for(Dipole_Wrappers::const_iterator it=m_dipole_wrappers.begin();
	it!=m_dipole_wrappers.end(); ++it) (*it)->SetKFactor(args);
  }


  void RS_Process::SetSelector(const PHASIC::Selector_Key &key)
  {
    PHASIC::Single_Process::SetSelector(key);
    for(Dipole_Wrappers::const_iterator it=m_dipole_wrappers.begin();
	it!=m_dipole_wrappers.end(); ++it) (*it)->SetSelector(key);
  }


  void RS_Process::SetNLOMC(PDF::NLOMC_Base *const nlomc)
  {
    PHASIC::Process_Base::SetNLOMC(nlomc);
    for(Dipole_Wrappers::const_iterator it=m_dipole_wrappers.begin();
	it!=m_dipole_wrappers.end(); ++it) (*it)->SetNLOMC(nlomc);
  }


  void RS_Process::SetShower(PDF::Shower_Base *const ps)
  {
    PHASIC::Process_Base::SetShower(ps);
    for(Dipole_Wrappers::const_iterator it=m_dipole_wrappers.begin();
	it!=m_dipole_wrappers.end(); ++it) (*it)->SetShower(ps);
  }


  void RS_Process::SetGenerator(PHASIC::ME_Generator_Base *const gen)
  {
    PHASIC::Process_Base::SetGenerator(gen);
    for(Dipole_Wrappers::const_iterator it=m_dipole_wrappers.begin();
	it!=m_dipole_wrappers.end(); ++it) (*it)->SetGenerator(gen);
  }


  void RS_Process::FillProcessMap(PHASIC::NLOTypeStringProcessMap_Map *apmap)
  {
    PHASIC::Process_Base::FillProcessMap(apmap);
    for(Dipole_Wrappers::const_iterator it=m_dipole_wrappers.begin();
	it!=m_dipole_wrappers.end(); ++it) (*it)->FillProcessMap(apmap);
  }

  
  size_t RS_Process::SetMCMode(const size_t mcmode)
  {
    size_t previous_mode = PHASIC::Process_Base::SetMCMode(mcmode);
    for(Dipole_Wrappers::const_iterator it=m_dipole_wrappers.begin();
	it!=m_dipole_wrappers.end(); ++it) (*it)->SetMCMode(mcmode);
    m_mcmode = mcmode;
    return previous_mode;
  }

  
  RS_Process::~RS_Process() {
    for(Dipole_Vector::const_iterator it=m_dipoles.begin();
	it!=m_dipoles.end(); ++it)
      if(*it) delete *it;
    DeleteSubevents();
    DeleteDipoleWrappers();

//    const std::string dir = "histograms";
//    std::ofstream myfile;
//    myfile.open(dir+"/mreal_y46.yoda");
//      YODA::WriterYODA::write(myfile, m_hist_mreal_y46);
//    myfile.close();
//    myfile.open(dir+"/dipole_y46.yoda");
//      YODA::WriterYODA::write(myfile, m_hist_dipole_y46);
//    myfile.close();

//    myfile.open("alpha_min.yoda");
//      YODA::WriterYODA::write(myfile, m_hist_alpha_min);
//    myfile.close();
    m_myfile.close();
  }

  
  void RS_Process::DeleteDipoleWrappers()
  {
    for(Dipole_Wrappers::const_iterator it = m_dipole_wrappers.begin();
	it!=m_dipole_wrappers.end(); ++it)
      if((*it)) delete (*it);
    m_dipole_wrappers.clear();
  }


  void RS_Process::DeleteSubevents()
  {
    for(ATOOLS::NLO_subevtlist::const_iterator it = m_subevents.begin();
	it!=m_subevents.end(); ++it)
      if((*it)) delete (*it);
    m_subevents.clear();
  }

}
