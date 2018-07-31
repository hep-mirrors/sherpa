#include "AddOns/Recola/Recola_Born.H"
#include "AddOns/Recola/Recola_Interface.H"
#include "PHASIC++/Process/External_ME_Args.H"

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
    double aqcd(AlphaQCD());
    // TODO: where to get this from?
    double mur=91.2;
    int defflav=Recola_Interface::GetDefaultFlav();
    set_alphas_rcl(aqcd,
		   mur,
		   defflav);

    
    double res(0.0);
    if (m_amptype==12) 
      {
	METOOLS::Divergence_Array<double> darr; double born;
	Recola_Interface::EvaluateLoop(m_recola_id, momenta, born, darr);
	res = darr.Finite();
      }
    else
      {
	THROW(not_implemented, "Not implemented");
      }
    
    // Recola returns ME2 including 1/symfac, but Calc is supposed to return it
    // without 1/symfac, thus multiplying with symfac here
    return res*m_symfac;
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

  int id = Recola_Interface::RegisterBorn(args, 12);
  if (id<=0) return NULL;
  
  return new Recola_Born(args, id, 12);
}
