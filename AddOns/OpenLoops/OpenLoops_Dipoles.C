#include "AddOns/OpenLoops/OpenLoops_Interface.H"
#include "AddOns/OpenLoops/OpenLoops_Dipoles.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"


namespace OpenLoops {


  Dipole::Dipole(const PHASIC::Dipole_Info& di,
		 const int& id, const AmplitudeType& amptype)
    : CS_Dipole(di), m_ol_id(id), m_amptype(amptype)
  {
    switch(FlavType())
      {
      case PHASIC::FlavourType::qtoqg:
	m_prefac = 8.0*M_PI;
	break;
	
      case PHASIC::FlavourType::gtoqq:
	m_prefac = 8.0*M_PI*m_TR/m_CA;
	break;
	
      case PHASIC::FlavourType::gtogg:
	m_prefac = 16.0*M_PI;
	break;
	
      default:
	THROW(fatal_error, "Internal error");
      }

    switch (SubtractionType())
      {
      case 0: break;
      case 1: break;
      case 2: break;
      default: THROW(not_implemented, "Subtraction type "+
		     ATOOLS::ToString(SubtractionType())+
		     " not implemented")
      }
    
  }
  
  
  double Dipole::CalcCorrelator() const
  {

    double alphas = p_aqcd->Default()*p_aqcd->Factor();
    OpenLoops_Interface::SetParameter("alphas", alphas);
    
    /* <1,...,m;a,b| T_ij T_k |b,a;m,...m1> * T_ij^{-2} */
    double TijTk  = OpenLoops_Interface::EvaluateColorCorrelator(m_ol_id,
								 Momenta(),
								 BornIJ(), BornK(),
								 AmpType());
    
    double SC(0.0);
    /* <1,...,m;a,b| ptilde^\mu T_ij T_k ptilde^\nu |b,a;m,...m1> * T_ij^{-2} * ptilde^{-2} */
    if(FlavType()!=PHASIC::FlavourType::qtoqg)
      SC = OpenLoops_Interface::EvaluateSpinCorrelator(m_ol_id,
						       Momenta(),
						       CalcPtilde(),
						       BornIJ(), BornK(),
						       AmpType());
    
    return alphas * m_prefac * m_symfac * ( CalcA()*TijTk + CalcB()*SC );
  }
  

  double FF_Dipole::Calc() const
  {
    if(SubtractionType()!=0) THROW(not_implemented, "Not implemented");

    const ATOOLS::Vec4D& pi = m_kin.m_pi;
    const ATOOLS::Vec4D& pj = m_kin.m_pj;
    
    /* hep-ph/9605323v3 eq. (5.2) */
    return -1.0/(2.0*pi*pj) * CalcCorrelator();
  }


  double FI_Dipole::Calc() const
  {
    if(SubtractionType()!=0) THROW(not_implemented, "Not implemented");

    const ATOOLS::Vec4D& pi = m_kin.m_pi;
    const ATOOLS::Vec4D& pj = m_kin.m_pj;
    const double& x = m_kin.m_x;
    
    /* hep-ph/9605323v3 eq. (5.36) */
    return -1.0/(2.0*pi*pj*x) * CalcCorrelator();
  }


  double IF_Dipole::Calc() const
  {
    if(SubtractionType()!=0) THROW(not_implemented, "Not implemented");

    const ATOOLS::Vec4D& pi = m_kin.m_pi;
    const ATOOLS::Vec4D& pa = m_kin.m_pa;
    const double& x = m_kin.m_x;
    
    /* hep-ph/9605323v3 eq. (5.61) */
    return -1.0/(2.0*pi*pa*x) * CalcCorrelator();
  }

      
  double II_Dipole::Calc() const
  {
    const ATOOLS::Vec4D& pi = m_kin.m_pi;
    const ATOOLS::Vec4D& pa = m_kin.m_pa;
    const double& x = m_kin.m_x;

    /* hep-ph/9605323v3 eq. (5.136) */
    return -1.0/(2.0*pi*pa*x) * CalcCorrelator();
  }


  double FF_Dipole::CalcA() const
  {
    double zi = m_kin.m_zi;
    double zj = m_kin.m_zj;
    const double& y(m_kin.m_y);

    /* q->qg expression depends on the flavour assignment being
       i=quark, j=gluon. Need to respect that here by swapping z_i,z_j
       if neccessary. Does not affect other splittings, so can be done
       for all cases. */
    if(FlavI().IsGluon()) std::swap(zi,zj);
    
    /* Coefficients of \delta_{ss^\prime} or -g^{\mu\nu} in eq. (5.7)
       - (5.9). */
    if(FlavType()==PHASIC::FlavourType::qtoqg)
      return 2.0/(1.0-zi*(1.0-y)) - (1.+ zi);
    if(FlavType()==PHASIC::FlavourType::gtoqq)
      return 1.0;
    if(FlavType()==PHASIC::FlavourType::gtogg)
      return 1.0/(1.0-zi*(1.0-y)) + 1.0/(1-zj*(1.0-y)) - 2.0;
    
    THROW(fatal_error, "Internal error");
  }


  double FI_Dipole::CalcA() const
  {
    double zi = m_kin.m_zi;
    double zj = m_kin.m_zj;
    const double& x(m_kin.m_x);

    /* q->qg expression depends on the flavour assignment being
       i=quark, j=gluon. Need to respect that here by swapping z_i,z_j
       if neccessary. Does not affect other splittings, so can be done
       for all cases. */
    if(FlavI().IsGluon()) std::swap(zi,zj);
    
    /* Coefficients of \delta_{ss^\prime} or -g^{\mu\nu} in eq. (5.39)
       - (5.41). */
    if(FlavType()==PHASIC::FlavourType::qtoqg)
      return 2.0/(1.0-zi+(1.0-x)) - (1.+ zi);
    if(FlavType()==PHASIC::FlavourType::gtoqq)
      return 1.0;
    if(FlavType()==PHASIC::FlavourType::gtogg)
      return 1.0/(1.0-zi+(1.0-x)) + 1.0/(1-zj+(1.0-x)) - 2.0;
    
    THROW(fatal_error, "Internal error");
  }


  double IF_Dipole::CalcA() const
  {
    const double& x  = (m_kin.m_x);
    const double& ui = (m_kin.m_ui);

    /* Need this to distinguish (5.65) from (5.66) */
    const ATOOLS::Flavour& flav_a = RealFlavours()[std::min(I(),J())];
    
    /* Coefficients of \delta_{ss^\prime} or -g^{\mu\nu} in eq. (5.65)
       - (5.68). */
    if((FlavType()==PHASIC::FlavourType::qtoqg) && flav_a.IsQuark())
      return 2.0/(1.0-x+ui) - (1.+x);
    if((FlavType()==PHASIC::FlavourType::qtoqg) && flav_a.IsGluon())
      return 1.0-2.0*x*(1.0-x);
    if(FlavType()==PHASIC::FlavourType::gtoqq)
      return x;
    if(FlavType()==PHASIC::FlavourType::gtogg)
      return 1/(1.0-x+ui)-1.0+x*(1.0-x);
    
    THROW(fatal_error, "Internal error");
  }

  
  double II_Dipole::CalcA() const
  {
    const double& x  = (m_kin.m_x);
    const double& z  = (SubtractionType()==1) ? m_kin.m_x+m_kin.m_v : x;

    /* Need this to distinguish (5.145) from (5.147) */
    const ATOOLS::Flavour& flav_a = RealFlavours()[std::min(I(),J())];
    
    /* Coefficients of \delta_{ss^\prime} or -g^{\mu\nu} in eq. (5.145)
       - (5.148). */
    if((FlavType()==PHASIC::FlavourType::qtoqg) && flav_a.IsQuark())
      return 2.0/(1.0-x) - (1.+z);

    if((FlavType()==PHASIC::FlavourType::qtoqg) && flav_a.IsGluon())
      return 1.0-2.0*z*(1.0-z);

    if(FlavType()==PHASIC::FlavourType::gtoqq)
      return z;
	  
    if(FlavType()==PHASIC::FlavourType::gtogg)
      return x/(1.0-x)+z*(1.0-z);
    
    THROW(fatal_error, "Internal error");
  }

  
  ATOOLS::Vec4D FF_Dipole::CalcPtilde() const
  {
    /* \mu-\nu tensor structure in hep-ph/9605323v3 eq. (5.8), (5.9)  */
    return m_kin.m_zi*m_kin.m_pi - m_kin.m_zj*m_kin.m_pj;
  }
  
  
  ATOOLS::Vec4D FI_Dipole::CalcPtilde() const
  {
    /* \mu-\nu tensor structure in hep-ph/9605323v3 eq. (5.40), (5.41)  */
    return m_kin.m_zi*m_kin.m_pi - m_kin.m_zj*m_kin.m_pj;
  }
  
  
  ATOOLS::Vec4D IF_Dipole::CalcPtilde() const
  {
    /* \mu-\nu tensor structure in hep-ph/9605323v3 eq. (5.67), (5.68)  */
    return m_kin.m_pi/m_kin.m_ui - m_kin.m_pk/(1.0-m_kin.m_ui);
  }
  

  ATOOLS::Vec4D II_Dipole::CalcPtilde() const
  {
    /* \mu-\nu tensor structure in hep-ph/9605323v3 eq. (5.147), (5.148)  */
    return m_kin.m_pi - (m_kin.m_pi*m_kin.m_pa)/(m_kin.m_pb*m_kin.m_pa) * m_kin.m_pb;
  }


  double FF_Dipole::CalcB() const
    {
      const double& zi(m_kin.m_zi);
      const double& zj(m_kin.m_zj);

      if(FlavType()==PHASIC::FlavourType::qtoqg)
	return 0.0;
      if(FlavType()==PHASIC::FlavourType::gtoqq)
	return -4.0*zi*zj;
      if(FlavType()==PHASIC::FlavourType::gtogg)
	return +2.0*zi*zj;

      THROW(fatal_error, "Internal error");
    }


  double FI_Dipole::CalcB() const
    {
      const double& zi(m_kin.m_zi);
      const double& zj(m_kin.m_zj);

      if(FlavType()==PHASIC::FlavourType::qtoqg)
	return 1.0;
      if(FlavType()==PHASIC::FlavourType::gtoqq)
	return -4.0*zi*zj;
      if(FlavType()==PHASIC::FlavourType::gtogg)
	return +2.0*zi*zj;

      THROW(fatal_error, "Internal error");
    }


  double IF_Dipole::CalcB() const
    {
      const double& x(m_kin.m_x);

      /*signs swapped as in FI*/
      if(FlavType()==PHASIC::FlavourType::qtoqg)
	return 1.0;
      if(FlavType()==PHASIC::FlavourType::gtoqq)
	return +4.0*(1.0-x)/x;
      if(FlavType()==PHASIC::FlavourType::gtogg)
	return +2.0*(1.0-x)/x;

      THROW(fatal_error, "Internal error");
    }


  double II_Dipole::CalcB() const
    {
      const double& x(m_kin.m_x);
      const double& v(m_kin.m_v);

      /*signs swapped as in FI*/
      if(FlavType()==PHASIC::FlavourType::qtoqg)
	return 1.0;
      
      if(FlavType()==PHASIC::FlavourType::gtoqq)
	{
	  switch(SubtractionType())
	    {
	    case 0:
	      return  +4.0*(1.0-x)/x;
	    case 1:
	      return  +4.0*( (x+v)/(sqr(x+v) +v*(1-x-v)) -1);
	    case 2:
	      return  +4.0*(1.0/(x+v)-1.0);
	    default:
	      THROW(not_implemented, "Not implemented");
	    }
	}
      
      if(FlavType()==PHASIC::FlavourType::gtogg)
	{
	  switch(SubtractionType())
	    {
	    case 0:
	      return +2.0*(1.0-x)/x;
	    case 1:
	      return +2.0*( (x+v)/(sqr(x+v) +v*(1-x-v)) -1);
	    case 2:
	      return +2.0*(1.0/(x+v)-1.0);
	    default:
	      THROW(not_implemented, "Not implemented");
	    }
	}

      THROW(fatal_error, "Internal error");
    }
}

using namespace OpenLoops;


DECLARE_CSDIPOLE_GETTER(OpenLoops::Dipole, "OpenLoops::Dipole")


PHASIC::CS_Dipole* ATOOLS::Getter<PHASIC::CS_Dipole,
				  PHASIC::Dipole_Info,
				  OpenLoops::Dipole>::
operator()(const PHASIC::Dipole_Info &di) const
{

  /* TODO: proper handling of coupling setting to make sure correct
     library is loaded */
  OpenLoops_Interface::SetParameter("coupling_qcd_0", -1);
  OpenLoops_Interface::SetParameter("coupling_qcd_1", -1);
  OpenLoops_Interface::SetParameter("coupling_ew_0",  -1);
  OpenLoops_Interface::SetParameter("coupling_ew_1",  -1);

  /* Construct a Process_Info instance in order to load library */
  const ATOOLS::Flavour_Vector& flavs
    = PHASIC::CS_Dipole::ConstructBornFlavours(di.m_real_i,
					       di.m_real_j,
					       di.m_real_flavs);
  PHASIC::Subprocess_Info ii;  PHASIC::Subprocess_Info fi;
  for (size_t n=0; n<2; n++)
    ii.m_ps.push_back(PHASIC::Subprocess_Info(flavs[n]));
  for (size_t n=2; n<flavs.size(); n++)
    fi.m_ps.push_back(PHASIC::Subprocess_Info(flavs[n]));

  AmplitudeType born_types[2] = {Tree,Loop2}; int id(-1);
  for (size_t i=0; i<2; ++i)
    {
      AmplitudeType type(born_types[i]);
      id = OpenLoops_Interface::RegisterProcess(ii, fi, type);
      if (id>0)
	{
	  switch(di.m_split_type)
	    {
	    case PHASIC::SplittingType::FF:
	      return new FF_Dipole(di,id,type);
	    case PHASIC::SplittingType::FI:
	      return new FI_Dipole(di,id,type);
	    case PHASIC::SplittingType::IF:
	      return new IF_Dipole(di,id,type);
	    case PHASIC::SplittingType::II:
	      return new II_Dipole(di,id,type);
	    default:
	      THROW(fatal_error, "Internal error");
	    }
	}
    }
  
  return NULL;
}
