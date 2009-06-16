#include "AddOns/Analysis/Detector/ElMag_Calorimeter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"


using namespace ANALYSIS;
using namespace ATOOLS;

DECLARE_GETTER(ECal_Getter,"ECal",
	       Analysis_Object,Argument_Matrix);

Analysis_Object *
ECal_Getter::operator()(const Argument_Matrix &parameters) const
{			
  if (parameters.size()<1) return NULL;
  //if (parameters.size()==1) abort(); // For read-in of, like 'ATLAS'

  ElMag_Calorimeter * ecal = new ElMag_Calorimeter(parameters());
  Detector_Segment * segment(NULL);
  double   etamin(0.),etamax(0.);
  long int neta(0),nphi(0);
  std::string name("");

  for (size_t i=0;i<parameters.size();++i) {
    const std::vector<std::string> &cur=parameters[i];
    if (cur.size()<2) continue;
    else if (cur[0]=="Segment") {
      etamin  = ATOOLS::ToType<double>(cur[1]);
      etamax  = ATOOLS::ToType<double>(cur[2]);
      neta    = ATOOLS::ToType<long int>(cur[3]);
      nphi    = ATOOLS::ToType<long int>(cur[4]);
      name    = cur[5];
      segment = new Detector_Segment(etamin,etamax,neta,nphi); 
      ecal->AddDetectorSegment(segment);
    }
  }
  return ecal;
}									

void ECal_Getter::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"Keyword=Segment   Parameters=etamin, etamax, neta, nphi, name"<<std::endl; 
}

ElMag_Calorimeter::ElMag_Calorimeter(Primitive_Analysis * ana) :
  Detector_Element(ana,"ECal")
{
  m_isdet = true;
}

ElMag_Calorimeter::~ElMag_Calorimeter() { }

Analysis_Object * ElMag_Calorimeter::GetCopy() const {
  ElMag_Calorimeter * ecal = new ElMag_Calorimeter(p_ana);
  for (std::set<Detector_Segment *,DS_Order>::iterator ds=m_segments.begin();
       ds!=m_segments.end();ds++) {
    ecal->AddDetectorSegment(dynamic_cast<Detector_Segment *>((*ds)->GetCopy()));
  }
  return ecal;
}


void ElMag_Calorimeter::Reset() {
  for (std::set<Detector_Segment *,DS_Order>::iterator ds=m_segments.begin();
       ds!=m_segments.end();ds++) (*ds)->Reset();
  if (!m_hitcells.empty()) {
    for (std::list<Cell *>::iterator cit=m_hitcells.begin();cit!=m_hitcells.end();++cit)
      (*cit)->Reset();
    m_hitcells.clear();
  }
}

bool ElMag_Calorimeter::Fill(const double E,const double eta,
			     const double phi,ATOOLS::Particle * part) {
  double etamin,etamax;
  for (std::set<Detector_Segment *,DS_Order>::iterator ds=m_segments.begin();
       ds!=m_segments.end();ds++) {
    (*ds)->Dimensions(etamin,etamax);
    if ((*ds)->InEtaRange(eta)) {
      Cell * cell((*ds)->AddParticle(eta,phi,part,E));
      m_hitcells.push_back(cell);
      return true;
    }
  }
  return false;
}

bool ElMag_Calorimeter::BuildCluster(Cell * cell,std::vector<Cell *> & cluster,
				     const int dim,double & E,double & eta,double & phi) {
  cluster.clear();
  double etacell,phicell,Etest,etatest,phitest;
  int wini(-dim+1),winj(-dim+1);
  E = 0.;
  cell->Centroid(etacell,phicell);
  switch (dim) {
  case 3:
  case 2:
    for (int i=-dim+1;i<1;i++) {
      for (int j=-dim+1;j<1;j++) {
	Etest = Plaquette(cell,etatest,phitest,i,i+dim-1,j,j+dim-1);
	if (Etest>E) {
	  E    = Etest;
	  eta  = etatest; 
	  phi  = phitest;
	  wini = i;
	  winj = j;
	}  
	else if (E>0. && dabs(Etest/E-1.)<1.e-3 && 
		 sqr(etatest-etacell)+sqr(phitest-phicell)<sqr(eta-etacell)+sqr(phi-phicell)) {
	  E    = Etest;
	  eta  = etatest; 
	  phi  = phitest;
	  wini = i;
	  winj = j;
	}
      }
    }
    FillPlaquette(cell,cluster,wini,winj,dim);
    break;
  case 1:
  default:
    cell->Centroid(eta,phi);
    E = cell->TotalDeposit();
    cluster.push_back(cell);
  }
  return (!cluster.empty());
}


double ElMag_Calorimeter::Plaquette(Cell * cell,double & eta,double & phi,
				    const int strip1,const int strip2,
				    const int cell1,const int cell2) {
  Etastrip * strip_1(cell->GetEtastrip()), * strip_2(cell->GetEtastrip());
  for (int i=0;i>strip1;i--) { 
    strip_1 = strip_1->GetMinus();
    if (strip_1==NULL) return 0.; 
  }
  for (int i=0;i<strip2;i++) {
    strip_2 = strip_2->GetPlus();
    if (strip_2==NULL) return 0.; 
  }
  Cell * cell_1(cell),  * cell_2(cell), * cur(NULL);
  for (int i=0;i>cell1;i--)  cell_1  = cell_1->GetDown();
  for (int i=0;i<cell2;i++)  cell_2  = cell_2->GetUp();
  double * dim = new double[4],E(0.),entry,norm(0.),eta_entry,phi_entry;
  cell->Dimensions(dim);
  cell_1->Dimensions(dim); 
  double phimin(dim[2]);
  cell_2->Dimensions(dim); 
  double phimax(dim[3]);
  if (phimax>M_PI) phimax=M_PI*0.9999999;

  eta = 0.;
  phi = 0.;
  int ncells(cell2-cell1+1);
  for (int i=0;i<strip2-strip1+1;i++) {
    cur  = strip_1->GetZero();
    for (;;) {
      cur->Dimensions(dim);
      if (cell==cur || dim[2]>=phimin) break; 
      cur = cur->GetUp();
    }
    for (int j=0;j<ncells;j++) {
      entry = cur->TotalDeposit();
      if (entry>0.) {
	cur->Centroid(eta_entry,phi_entry);
	E    += entry;
	eta  += entry*eta_entry;
	phi  += entry*phi_entry;
	norm += entry;
      }
      cur->Dimensions(dim); 
      cur = cur->GetUp();
    }
    strip_1 = strip_1->GetPlus();
    if (strip_1==NULL) break;
  } 
  delete [] dim;

  eta /= norm; 
  phi /= norm;
  return E;
}


void ElMag_Calorimeter::FillPlaquette(Cell * cell,std::vector<Cell *> & cells,
				      const int etastrip,const int cellno,const int ncells) {
  Etastrip * strip_1(cell->GetEtastrip());
  for (int i=0;i>etastrip;i--) strip_1 = strip_1->GetMinus();
  double * dims = new double[4];
  Cell * cur(cell);
  for (int i=0;i>cellno;i--) cur = cur->GetDown();
  cur->Dimensions(dims); 
  double phimin(dims[2]);
  cur = cell;
  for (int i=0;i<cellno+ncells;i++) cur = cur->GetUp();
  cur->Dimensions(dims); 
  double phimax(dims[3]);

  for (int i=0;i<ncells;i++) {
    cur  = strip_1->GetZero();
    do { 
      cur->Dimensions(dims); 
      cur = cur->GetUp();
    } while (dims[2]<phimin);
    for (int j=0;j<ncells;j++) {
      if (cur->TotalDeposit()>0. && !cur->Used()) {
	cur->SetUsed(true);
	cells.push_back(cur);
      }
      cur->Dimensions(dims); 
      if (dims[3]>phimax) break;
      cur = cur->GetUp();
    }
    strip_1 = strip_1->GetPlus();
    if (strip_1==NULL) break;
  } 
  delete [] dims;
}
