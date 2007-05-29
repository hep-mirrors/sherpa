#include "MET_Maker.H"
#include "Primitive_Analysis.H"
#include "Detector.H"
#include "Message.H"
#include "MyStrStream.H"
#include "Exception.H"


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
  m_cutHad(1.), m_cutEM(1.)
{  
  m_kfcode=kf::none;
  GetElements();
}

MET_Maker::~MET_Maker() {}

void MET_Maker::SetCuts(const double cutHad,const double cutEM) {
  m_cutHad = cutHad; m_cutEM = cutEM;
}

void MET_Maker::ReconstructObjects(ATOOLS::Particle_List * plist,ATOOLS::Vec4D & METvector) {
  //std::cout<<METHOD<<" for "<<p_ECal->GetHitCells()->size()<<" hit cells in ECal and "
  //   <<p_HCal->GetHitCells()->size()<<" hit cells in HCal."<<std::endl;

  std::list<Cell *> * cells = p_ECal->GetHitCells();
  Cell * cell(NULL);

  Vec4D checkmom(0.,0.,0.,0.);
  for (std::list<Cell *>::iterator cit(cells->begin());
       cit!=cells->end();) {
    cell = (*cit);
    if (!cell->Used() && cell->TotalDeposit()>m_cutEM) {
      METvector -= cell->TotalDeposit()*cell->Direction();
      cell->Reset();
      cit = cells->erase(cit);
    }
    else {
      checkmom += cell->TotalDeposit()*cell->Direction();
      cell->Reset();
      cit++;
    }
  }
  cells = p_HCal->GetHitCells();
  for (std::list<Cell *>::iterator cit(cells->begin());
       cit!=cells->end();) {
    cell = (*cit);
    if (!cell->Used() && cell->TotalDeposit()>m_cutHad) {
      METvector -= cell->TotalDeposit()*cell->Direction();
      cell->Reset();
      cit = cells->erase(cit);
    }
    else {
      checkmom += cell->TotalDeposit()*cell->Direction();
      cell->Reset();
      cit++;
    }
  }

  std::list<Track *> tracks;
  p_chambers->GetTracks(tracks);
  Track * track(NULL);
  for (std::list<Track *>::iterator trit=tracks.begin();trit!=tracks.end();) {
    track = (*trit);
    if (!track->used && track->mom[0]>m_cutEM) {
      METvector -= track->mom;
      trit = tracks.erase(trit);
    }
    else {
      checkmom += track->mom;
      trit++;
    }
  }


  METvector[3] = 0.;
  METvector[0] = METvector.PPerp();
  Particle * part = new Particle(0,Flavour(m_kfcode),METvector,'r');
  plist->push_back(part);
}
