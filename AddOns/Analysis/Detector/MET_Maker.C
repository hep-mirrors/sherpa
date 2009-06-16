#include "AddOns/Analysis/Detector/MET_Maker.H"
#include "AddOns/Analysis/Main/Primitive_Analysis.H"
#include "AddOns/Analysis/Detector/Detector.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"


using namespace ANALYSIS;
using namespace ATOOLS;

DECLARE_GETTER(MET_Maker_Getter,"MET_Maker",
	       Analysis_Object,Argument_Matrix);

Analysis_Object *
MET_Maker_Getter::operator()(const Argument_Matrix &parameters) const
{			
  if (parameters.size()<1) return NULL;
  //if (parameters.size()==1) abort(); // For read-in of, like 'ATLAS'


  std::string mode("ET_UP");
  MET_Maker * maker = new MET_Maker(parameters(),mode);
  double cutHad,cutEM;

  for (size_t i=0;i<parameters.size();++i) {
    const std::vector<std::string> &cur=parameters[i];
    if (cur.size()<2) continue;
    else if (cur[0]=="ECorrection") {
      if (cur[1]=="truth") {
	maker->SetECorrection(0,ATOOLS::ToType<double>(cur[2]));
      }
      else if (cur[1]=="constant") {
	maker->SetECorrection(1,ATOOLS::ToType<double>(cur[2]),0.,0.,
			      ATOOLS::ToType<double>(cur[5]),0.,0.);
      }
      else if (cur[1]=="parametrised") {
	maker->SetECorrection(2,ATOOLS::ToType<double>(cur[2]),ATOOLS::ToType<double>(cur[3]),
			      ATOOLS::ToType<double>(cur[4]),ATOOLS::ToType<double>(cur[5]),
			      ATOOLS::ToType<double>(cur[6]),ATOOLS::ToType<double>(cur[7]));
      }
    }
    else if (cur[0]=="Cuts") {
      cutHad = ATOOLS::ToType<double>(cur[1]);
      cutEM  = ATOOLS::ToType<double>(cur[2]);
      maker->SetCuts(cutHad,cutEM);
    }
  }
  return maker;
}									

void MET_Maker_Getter::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"Keyword=Acceptance   Parameters=etamin, etamax"<<std::endl; 
}

MET_Maker::MET_Maker(Primitive_Analysis * ana,const std::string mode) : 
  Object_Definition_Base(ana,"METMaker",mode),
  m_corrmode(0),
  m_cutHad(1.), m_cutEM(1.), m_deviation(0.02),
  m_Rinfty_ECal(0.78), m_ampl_ECal(0.37), m_slope_ECal(0.2),
  m_Rinfty_HCal(0.81), m_ampl_HCal(0.45), m_slope_HCal(0.2)
{  
  m_kfcode=kf_none;
  GetElements();
}

MET_Maker::~MET_Maker() {}

void MET_Maker::SetCuts(const double cutHad,const double cutEM) {
  m_cutHad = cutHad; m_cutEM = cutEM;
}

void MET_Maker::
SetECorrection(const int mode,
	       const double rinfty_ecal,const double ampl_ecal, const double slope_ecal,
	       const double rinfty_hcal,const double ampl_hcal, const double slope_hcal) {
  switch (mode) {
  case 2:
  case 1:
    m_corrmode    = 2;
    m_Rinfty_ECal = rinfty_ecal; 
    m_ampl_ECal   = ampl_ecal; 
    m_slope_ECal  = slope_ecal;
    m_Rinfty_HCal = rinfty_hcal; 
    m_ampl_HCal   = ampl_hcal; 
    m_slope_HCal  = slope_hcal;
    break;
  case 0:
  default:
    m_corrmode    = 0;
    m_deviation   = rinfty_ecal;
  }
}

void MET_Maker::ReconstructObjects(ATOOLS::Particle_List * plist,ATOOLS::Vec4D & METvector) {
  Vec4D cellmom(0.,0.,0.,0.), corrmom(0.,0.,0.,0.);
  Vec4D truemom(0.,0.,0.,0.), truecell;
  std::map<Particle *,Vec4D> parts;
  Particle * part;

  Cell * cell(NULL);
  std::list<Cell *> * cells = p_ECal->GetHitCells();

  msg_Debugging()<<METHOD<<" with ECAL:"<<std::endl;
  for (std::list<Cell *>::iterator cit(cells->begin());
       cit!=cells->end();) {
    cell = (*cit);
    if (!cell->Used() && cell->TotalDeposit()>m_cutEM) {
      cellmom += cell->TotalDeposit()*cell->Direction();
      corrmom += ReconstructedE_ECal(cell->TotalDeposit())*cell->Direction();
      msg_Debugging()<<"   Deposit = "<<cell->TotalDeposit()*cell->Direction()
	       <<" --> reconstructed = "
      	       <<ReconstructedE_ECal(cell->TotalDeposit())*cell->Direction();
      truecell = Vec4D(0.,0.,0.,0.);
      for (std::map<ATOOLS::Particle *,double>::iterator 
	     pit=cell->ParticleEntries()->begin();
	   pit!=cell->ParticleEntries()->end();pit++) {
	part = pit->first->OriginalPart();
	truecell += part->Momentum();
	if (parts.find(part)==parts.end()) {
	  truemom  += part->Momentum();
	  parts[part] = ReconstructedE_ECal(cell->TotalDeposit())*cell->Direction();
	}
	else parts[part] += ReconstructedE_ECal(cell->TotalDeposit())*cell->Direction();
      }
      msg_Debugging()<<" vs. particles = "<<truecell<<" from "<<cell->ParticleEntries()->size()
      	       <<"("<<cell->ParticleEntries()->begin()->first->Flav()<<" @ "
	       <<cell->ParticleEntries()->begin()->first->Momentum().Eta()<<")."<<std::endl;
      cell->Reset();
      cit = cells->erase(cit);
    }
    else {
      cell->Reset();
      cit++;
    }
  }

  msg_Debugging()<<METHOD<<" with HCAL:"<<std::endl;
  cells = p_HCal->GetHitCells();
  for (std::list<Cell *>::iterator cit(cells->begin());
       cit!=cells->end();) {
    cell = (*cit);
    if (!cell->Used() && cell->TotalDeposit()>m_cutHad) {
      cellmom += cell->TotalDeposit()*cell->Direction();
      corrmom += ReconstructedE_HCal(cell->TotalDeposit())*cell->Direction();
      msg_Debugging()<<"   Deposit = "<<cell->TotalDeposit()*cell->Direction()
	       <<" --> reconstructed = "
      	       <<ReconstructedE_HCal(cell->TotalDeposit())*cell->Direction();
      truecell = Vec4D(0.,0.,0.,0.);
      for (std::map<ATOOLS::Particle *,double>::iterator 
	     pit=cell->ParticleEntries()->begin();
	   pit!=cell->ParticleEntries()->end();pit++) {
	part = pit->first->OriginalPart();
	truecell += part->Momentum();
	if (parts.find(part)==parts.end()) {
	  truemom  += part->Momentum();
	  parts[part] = ReconstructedE_ECal(cell->TotalDeposit())*cell->Direction();
	}
	else parts[part] += ReconstructedE_ECal(cell->TotalDeposit())*cell->Direction();
      }
      msg_Debugging()<<" vs. particles = "<<truecell<<" from "<<cell->ParticleEntries()->size()
      	       <<"("<<cell->ParticleEntries()->begin()->first->Flav()<<")."<<std::endl;
      cell->Reset();
      cit = cells->erase(cit);
    }
    else {
      cell->Reset();
      cit++;
    }
  }
  msg_Debugging()<<"   Check particles:"<<std::endl;
  for (std::map<Particle *,Vec4D>::iterator pvit=parts.begin();pvit!=parts.end();pvit++) {
    msg_Debugging()<<"      "<<pvit->first->Flav()<<": original mom = "<<pvit->first->Momentum()
  	     <<" vs. reconstructed mom = "<<pvit->second<<std::endl;
  }
  msg_Debugging()<<"   Cellmom after HCAL: corrected = "<<corrmom<<" (from deposit = "<<cellmom<<")"
	   <<" vs. true = "<<truemom<<std::endl;
  
  std::list<Track *> tracks;
  p_chambers->GetTracks(tracks);
  Track * track(NULL);
  for (std::list<Track *>::iterator trit=tracks.begin();trit!=tracks.end();) {
    track = (*trit);
    if (!track->used && track->mom[0]>m_cutEM) {
      msg_Debugging()<<"   Add track: "<<track->flav<<" with "<<track->mom<<std::endl;
      cellmom += track->mom;
      truemom += track->mom;
      trit = tracks.erase(trit);
    }
    else {
      trit++;
    }
  }

  if (m_corrmode==2) METvector -= cellmom;
  else if (m_corrmode==0) {
    double rana,dummy;
    do { ran.Gaussian(rana,dummy); } while (dabs(rana)>2.*M_PI);
    METvector -= (1.+m_deviation*rana/M_PI)*truemom;
  }
  msg_Debugging()<<METHOD<<" reconstructed MET : "<<METvector<<std::endl;

  METvector[3] = 0.;
  METvector[0] = METvector.PPerp();

  plist->push_back(new Particle(0,Flavour(m_kfcode),METvector,'r'));
}


double MET_Maker::ReconstructedE_ECal(const double dep) {
  return dep/(m_Rinfty_ECal-m_ampl_ECal*exp(-m_slope_ECal*dep));
}

double MET_Maker::ReconstructedE_HCal(const double dep) {
  return dep/(m_Rinfty_HCal-m_ampl_HCal*exp(-m_slope_HCal*dep));
}
