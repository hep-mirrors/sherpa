#include "Gluon_Decayer.H"
#include "Hadronisation_Parameters.H"
#include "Message.H"

using namespace AHADIC;
using namespace ATOOLS;


Gluon_Decayer::Gluon_Decayer(Dipole_Splitter * splitter) :
  p_splitter(splitter)
{ 
  double norm(0.);
  for (FlavCCMap_Iterator fdit=hadpars.GetConstituents()->CCMap.begin();
       fdit!=hadpars.GetConstituents()->CCMap.end();fdit++) {
    if (hadpars.GetConstituents()->TotWeight(fdit->first)>norm)
      norm = hadpars.GetConstituents()->TotWeight(fdit->first);
  }
  DecaySpecs * decspec;
  for (FlavCCMap_Iterator fdit=hadpars.GetConstituents()->CCMap.begin();
       fdit!=hadpars.GetConstituents()->CCMap.end();fdit++) {
    if (!fdit->first.IsAnti()) {
      decspec = new DecaySpecs;
      decspec->popweight = hadpars.GetConstituents()->TotWeight(fdit->first)/norm;
      decspec->massmin   = hadpars.GetConstituents()->Mass(fdit->first);
      m_options.insert(std::make_pair(fdit->first,decspec));
      msg_Tracking()<<"Insert option : g->"<<fdit->first<<" "<<fdit->first.Bar()<<std::endl;
    }
  }
  p_splitter->SetOptions(&m_options);

  msg_Tracking()<<"------------- END OF GLUON_DECAYER --------------"<<std::endl;
  if (m_options.empty()) {
    msg_Error()<<"Error in "<<METHOD<<":"<<std::endl
	       <<"   No decay channels found for gluons, will abort."<<std::endl;
    abort();
  }
}

Gluon_Decayer::~Gluon_Decayer() {
  for (FDIter fdit=m_options.begin();fdit!=m_options.end();++fdit) {
    if (fdit->second!=NULL) { delete fdit->second; fdit->second=NULL; }
  }
  m_options.clear();
}

bool Gluon_Decayer::DecayList(Proto_Particle_List * plin)
{
  if (plin==NULL || plin->empty()) return true;

  if (!Shift(plin)) { return false; }
  FillDipoleList(plin);
  PrintDipoleList();
  DecayDipoles();
  PrintDipoleList();
  abort();

  PPL_Iterator pit;
  bool success(true);
  Proto_Particle * part, * part1, * part2;
  pit=plin->begin();
  do {
  } while (pit!=plin->end());

  return success;
}

bool Gluon_Decayer::Shift(Proto_Particle_List * pl)
{
  size_t number(pl->size());
  bool val(true);
  if (number<2) return val; 
  std::vector<Vec4D>  momenta(number);
  std::vector<double> masses(number);
  int k(0);
  Flavour flav;
  PPL_Iterator pit;
  for (pit=pl->begin();pit!=pl->end();++pit,++k) {
    flav       = (*pit)->m_flav;
    momenta[k] = (*pit)->m_mom;
    masses[k]  = hadpars.GetConstituents()->Mass(flav);
  }
  if (!hadpars.AdjustMomenta(number,&momenta.front(),&masses.front())) val=false;

  k = 0;
  for (pit=pl->begin();pit!=pl->end();++pit,++k) (*pit)->m_mom = momenta[k]; 
  return val;
}

bool Gluon_Decayer::FillDipoleList(Proto_Particle_List * plin)
{
  PPL_Iterator pit(plin->begin()), pit1(plin->begin());
  pit1++;
  Proto_Particle * begin(*pit);
  Dipole * dip;
  do {
    dip = new Dipole(*pit,*pit1);
    m_dipoles.push_back(dip);
    pit = pit1;
    pit1++;
  } while (pit1!=plin->end());
  if ((*pit)->m_flav.IsGluon()) {
    if (begin->m_flav.IsGluon()) {
      dip = new Dipole(*pit,begin);
      m_dipoles.push_back(dip);
    }
    else {
      msg_Error()<<"ERROR in "<<METHOD<<":"<<std::endl
		 <<"    Last flavour in list = "<<(*pit)->m_flav
		 <<" but first flavour = "<<begin->m_flav<<"."<<std::endl
		 <<"   Don't know what to do, abort."<<std::endl;
      abort();
    }
  }
  return true;
}

bool Gluon_Decayer::DecayDipoles() {
  DipIter dip,neighbour;
  Proto_Particle * new1, * new2;
  do { 
    dip = SelectDipole(); 
    if (dip==m_dipoles.end()) return true;
    msg_Out()<<METHOD<<" the winner is "<<sqrt((*dip)->mass2)
	     <<" ("<<((*dip)->mustdecay)<<"): "<<std::endl
	     <<"  "<<(*dip)->triplet->m_flav<<"("<<(*dip)->triplet->m_mom<<"), "
	     <<"  "<<(*dip)->antitriplet->m_flav<<"("<<(*dip)->antitriplet->m_mom<<")."
	     <<std::endl;
    if (!p_splitter->Split((*dip),1.)) {
      msg_Error()<<"ERROR in "<<METHOD<<" :"<<std::endl
		 <<"   Could not split dipole.  Redo hadronization step."<<std::endl;
      abort();
    }
    p_splitter->GetNewParticles(new1,new2);
    neighbour = dip;
    if ((*dip)->switched) {
      delete (*dip)->triplet;
      (*dip)->triplet = new2;
      neighbour--;
      (*neighbour)->antitriplet = new1;
    }
    else {
      delete (*dip)->antitriplet;
      (*dip)->antitriplet = new1;
      neighbour++;
      (*neighbour)->triplet = new2;      
    }
    (*dip)->Update();
    (*neighbour)->Update();
  } while (dip!=m_dipoles.end());
  return true;
}

DipIter Gluon_Decayer::SelectDipole() {
  double smax(-1.);
  DipIter winner=m_dipoles.end();
  for (DipIter dip=m_dipoles.begin();dip!=m_dipoles.end();dip++) {
    if ((*dip)->mustdecay && (*dip)->massbar2>smax) {
      smax   = (*dip)->massbar2;
      winner = dip;
    }
  }
  return winner;
}

void Gluon_Decayer::PrintDipoleList()
{
  for (DipIter dip=m_dipoles.begin();dip!=m_dipoles.end();dip++) {
    msg_Out()<<"Dipole("<<((*dip)->mustdecay)<<", "<<sqrt((*dip)->mass2)<<") : "<<std::endl
	     <<"  "<<(*dip)->triplet->m_flav<<"("<<(*dip)->triplet->m_mom<<"), "
	     <<" "<<hadpars.GetConstituents()->Mass((*dip)->triplet->m_flav)
	     <<" vs. "<<sqrt(Max((*dip)->triplet->m_mom.Abs2(),0.0))<<";"<<std::endl
	     <<"  "<<(*dip)->antitriplet->m_flav<<"("<<(*dip)->antitriplet->m_mom<<"),"
	     <<" "<<hadpars.GetConstituents()->Mass((*dip)->antitriplet->m_flav)
	     <<" vs. "<<sqrt(Max((*dip)->antitriplet->m_mom.Abs2(),0.0))<<"."<<std::endl;
  }
}


