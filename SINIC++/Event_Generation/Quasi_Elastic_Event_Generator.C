#include "SINIC++/Event_Generation/Quasi_Elastic_Event_Generator.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Random.H"

using namespace SINIC;

Quasi_Elastic_Event_Generator::
Quasi_Elastic_Event_Generator(Elastic_Event_Generator * elastic,
			      Single_Diffractive_Event_Generator * singdiff,
			      Double_Diffractive_Event_Generator * doubdiff) :
  p_elastic(elastic), p_singdiff(singdiff),p_doubdiff(doubdiff),
  m_el(p_elastic->XSec()), m_sd(p_singdiff->XSec()), m_dd(p_doubdiff->XSec()),
  m_stot(m_el+m_sd+m_dd)
{
  msg_Tracking()<<METHOD<<" with xsecs: el = "<<(m_el/1.e9)<<", "
		<<"SD = "<<(m_sd/1.e9)<<", DD = "<<(m_dd/1.e9)<<", "
		<<"and tot = "<<(m_stot/1.e9)<<" mbarns."<<std::endl;
}


bool Quasi_Elastic_Event_Generator::
QuasiElasticEvent(ATOOLS::Blob_List * blobs,const double & xsec) {
  double disc(m_stot*0.99999999*ATOOLS::ran->Get());
  disc -= m_dd;
  if (disc<=0.)
    return p_doubdiff->DoubleDiffractiveEvent(blobs,m_stot);
  disc -= m_sd;
  if (disc<=0.)
    return p_singdiff->SingleDiffractiveEvent(blobs,m_stot);
  disc -= m_el;
  if (disc>0.) 
    msg_Error()<<"Potential error in "<<METHOD<<":"<<std::endl
	       <<"   Ignore it, generate an elastic event "
	       <<"and hope for the best."<<std::endl;
  return p_elastic->ElasticEvent(blobs,m_stot);
}
