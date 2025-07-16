#include "ExampleME.H"
#include "rambo.h"

#include "PHASIC++/Process/Tree_ME2_Base.H"
#include "PHASIC++/Process/Process_Info.H"
#include "ATOOLS/Org/Getter_Function.H"
#include "ATOOLS/Org/Exception.H"

#include <assert.h>

ExampleME::ExampleME(const PHASIC::Process_Info& pi) :
  Tree_ME2_Base(pi, pi.ExtractFlavours())
{
  
  m_symfac = pi.m_fi.FSSymmetryFactor();
  m_symfac*= pi.m_ii.ISSymmetryFactor();

  m_nin    = pi.m_ii.NExternal();
  m_nout   = pi.m_fi.NExternal();
  m_pinfo  = pi;

  assert(pi.m_maxcpl[0]==pi.m_mincpl[0]);
  assert(pi.m_maxcpl[1]==pi.m_mincpl[1]);

  m_oqcd = pi.m_maxcpl[0];
  m_oew  = pi.m_maxcpl[1];

  if(m_mg_proc.nprocesses!=1)
    THROW(fatal_error, "MG process contrains several subprocesses")

  // Read param_card and set parameters
  m_mg_proc.initProc("./param_card.dat");
  m_mg_proc.setInitial((int)m_flavs[0],(int)m_flavs[1]);
}

double ExampleME::RunningCouplingFactor() const
{
  if(!(p_aqed&&p_aqcd))
    THROW(fatal_error, "Running couplings not set");
  double fac(1.0);
  if (OrderQCD()!=0) fac*=pow(p_aqcd->Factor(), OrderQCD());
  if (OrderEW() !=0) fac*=pow(p_aqed->Factor(), OrderEW() );
  return fac;
}


double ExampleME::Calc(const ATOOLS::Vec4D_Vector& momenta)
{

  // Convert momenta
  std::vector< double* > p(momenta.size(), NULL);
  for(int i(0); i<momenta.size(); i++){
    p[i] = new double[4];
    p[i][0] = momenta[i][0];
    p[i][1] = momenta[i][1];
    p[i][2] = momenta[i][2];
    p[i][3] = momenta[i][3];
  }
  // Set momenta for this event
  m_mg_proc.setMomenta(p);
  // Evaluate matrix element
  m_mg_proc.sigmaKin();
  // Get Matrix element
  double result = m_mg_proc.getMatrixElements()[0];

  // Free memory that held momenta temporarily
  for(std::vector<double*>::iterator it=p.begin();
      it!=p.end(); ++it)
    delete [] (*it);

  // OL returns ME2 including 1/symfac, but Calc is supposed to return it
  // without 1/symfac, thus multiplying with symfac here
  return m_symfac*result*RunningCouplingFactor();
}

using namespace PHASIC;
DECLARE_TREEME2_GETTER(ExampleME,"ExampleME")
PHASIC::Tree_ME2_Base *ATOOLS::Getter<PHASIC::Tree_ME2_Base,
				      PHASIC::Process_Info,
				      ExampleME>::
operator()(const PHASIC::Process_Info &pi) const
{
  DEBUG_FUNC(pi);
  /* Running of coupling can only be applied here if the process has
     definite QCD and QED order */
  if(pi.m_maxcpl[0]!=pi.m_mincpl[0])
    return NULL;
  if(pi.m_maxcpl[1]!=pi.m_mincpl[1])
    return NULL;
  return new ExampleME(pi);
}
