#include "HT.H"

DEFINE_OBSERVABLE_GETTER(Multiplicity,Multiplicity_Getter,"Multi");
 
HT::HT(int type,double xmin,double xmax,int nbins,
			   const std::string & listname) :
  Primitive_Observable_Base(type,xmin,xmax,nbins,NULL)
{
  if (listname!="") {
    m_listname = listname;
    m_name = listname+"_HT.dat";
  }
  else
    m_name = "HT.dat";
}

void HT::Evaluate(const ATOOLS::Particle_List & pl,
		  double weight, int ncount)
{
  ATOOLS::Particle_List *jets=p_ana->GetParticleList(listname);
  double HT=0.0;
  for (ATOOLS::Particle_List::const_iterator pit=jets->begin();
       pit!=p_jets->end();++pit) {
    HT+=(*pit)->Momentum().PPerp();
  }
  p_histo->Insert(HT,weight,ncount); 
}


Primitive_Observable_Base * HT::Copy() const 
{
  return new HT(m_type,m_xmin,m_xmax,m_nbins,m_listname);
}
