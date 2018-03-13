#include "AddOns/DiHiggsNLO/DiHiggs_Virtual.H"
#include "AddOns/DiHiggsNLO/hhgrid.h"

#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Library_Loader.H"
#include "ATOOLS/Org/CXXFLAGS.H"

#include "PHASIC++/Process/Process_Info.H"

#include "MODEL/Main/Model_Base.H"
#include "MODEL/Main/Running_AlphaS.H"

#include <assert.h>

namespace DiHiggs {

  DiHiggs_Virtual::DiHiggs_Virtual(const PHASIC::Process_Info& pi,
				   const ATOOLS::Flavour_Vector& flavs,
				   const std::string& grid_path) :
    PHASIC::Virtual_ME2_Base(pi, flavs)
  {
    DEBUG_FUNC(this);

    m_symfac =pi.m_fi.FSSymmetryFactor();
    m_symfac*=pi.m_ii.ISSymmetryFactor();
    m_mode = 1;

    /* Use these temporary char fields to avoid compiler warnings */
    char pth    [] = "path";
    
    /* Init python and modify python search path to include
       installation directory of 'creategrid.py' */
    Py_Initialize();
    PyList_Append(PySys_GetObject(pth),
		  PyString_FromString(PYTHON_LIBS));

    if(!ATOOLS::FileExists(grid_path))
      THROW(fatal_error, "Invalid grid path " + grid_path)
	
    msg_Debugging() << "Initializing grid from file: " << grid_path << std::endl;
    p_grid = grid_initialize(grid_path.c_str());

    if(!p_grid)
      THROW(fatal_error, "Failed to initialize python grid");

    /* Get running alphaS for implementing renormalization scale dependence */
    p_as = dynamic_cast<MODEL::Running_AlphaS*>(MODEL::s_model->GetScalarFunction("alpha_S"));
    if(!p_as) THROW(fatal_error, "Failed to retrieve running alphaS");

    double NF = Flavour(kf_quark).Size()/2.0;
    if(NF!=5.0) THROW(fatal_error, "Only compatible with five flavour scheme");
    
    m_beta0 =  11. - 10.0/3.0;

  }
  

  DiHiggs_Virtual::~DiHiggs_Virtual()
  {
    Py_DECREF(p_grid);
    Py_Finalize();
  }

  
  double DiHiggs_Virtual::Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& p)
  {
    /* They pull out a factor of 4\pi^\Epsilon/\Gamma(1-\Epsilon) in
       frot of the virtual amplitude. 1/\Gamma(1-\Epsilon) is alsways
       assumed, but we need to correct for the rest */
    return 4.*M_PI;
  }


  void DiHiggs_Virtual::Calc(const ATOOLS::Vec4D_Vector& p)
  {
    THROW(fatal_error, "Invalid call");
  }

			     
  void DiHiggs_Virtual::Calc(const ATOOLS::Vec4D_Vector& p,
			     const double& born)
  {
    double s   = (p[0]+p[1]).Abs2();
    double t   = (p[0]-p[2]).Abs2();
    double mu0 = s/4.0;

    double asratio = sqr((*p_as)(m_mur2)/(*p_as)(mu0));

    /* Division by symmetry factor already included in grid, need to
       multiply in order to comply with Sherpa conventions */
    double V0 = grid_virt(p_grid,s,t)*m_symfac;

    /* Renormalization scale dependence: arXiv:1703.09252 eq. (2.6)*/
    double V = V0*asratio + 3.0*born*(sqr(log(mu0/s)) - sqr(log(m_mur2/s)));

    m_res.Finite() = V;
    m_res.IR()     = -2.0*3.0*born*log(m_mur2/s) - born*m_beta0;
    m_res.IR2()    = -2.0*3.0*born;
  }
  
}

DECLARE_VIRTUALME2_GETTER(DiHiggs::DiHiggs_Virtual,"DiHiggs_Virtual")

PHASIC::Virtual_ME2_Base *ATOOLS::Getter<PHASIC::Virtual_ME2_Base,
					 PHASIC::Process_Info,
					 DiHiggs::DiHiggs_Virtual>
::operator()(const PHASIC::Process_Info &pi) const
{
  /* Require Loop_Generator=DiHiggs in process section */
  if (pi.m_loopgenerator!="DiHiggs") return NULL;

  /* Check NLO type (allow only QCD, not EW)  */
  if (pi.m_fi.m_nloewtype ==PHASIC::nlo_type::loop) return NULL;
  if (pi.m_fi.m_nloqcdtype!=PHASIC::nlo_type::loop) return NULL;

  /* Check flavours */
  ATOOLS::Flavour_Vector flavs = pi.ExtractFlavours();
  if(!flavs[0].IsGluon() )   return NULL;
  if(!flavs[1].IsGluon() )   return NULL;
  if( flavs[2].Kfcode()!=25) return NULL;
  if( flavs[2].Kfcode()!=25) return NULL;

  /* Get the full file path to the grid */
  std::string grid_path = ATOOLS::Default_Reader().
    Get<std::string>("DIHIGGS_GRID_PATH", "Virt_full.grid");

  return new DiHiggs::DiHiggs_Virtual(pi, flavs, grid_path);
}
