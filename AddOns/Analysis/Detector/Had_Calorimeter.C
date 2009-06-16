#include "AddOns/Analysis/Detector/Had_Calorimeter.H"
#include "PHASIC++/Channels/Rambo.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"


using namespace ANALYSIS;
using namespace ATOOLS;

DECLARE_GETTER(HCal_Getter,"HCal",
	       Analysis_Object,Argument_Matrix);

Analysis_Object *
HCal_Getter::operator()(const Argument_Matrix &parameters) const
{			
  if (parameters.size()<1) return NULL;
  //if (parameters.size()==1) abort(); // For read-in of, like 'ATLAS'

  Had_Calorimeter * ecal = new Had_Calorimeter(parameters());
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

void HCal_Getter::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"Keyword=Segment   Parameters=etamin, etamax, neta, nphi, name"<<std::endl; 
}

Had_Calorimeter::Had_Calorimeter(Primitive_Analysis * ana) :
  Detector_Element(ana,"HCal")
{
  m_isdet = true;
}

Had_Calorimeter::~Had_Calorimeter() { }

Analysis_Object * Had_Calorimeter::GetCopy() const {
  Had_Calorimeter * ecal = new Had_Calorimeter(p_ana);
  for (std::set<Detector_Segment *,DS_Order>::iterator ds=m_segments.begin();
       ds!=m_segments.end();ds++) 
    ecal->AddDetectorSegment(dynamic_cast<Detector_Segment *>((*ds)->GetCopy()));
  return ecal;
}


void Had_Calorimeter::Reset() {
  for (std::set<Detector_Segment *,DS_Order>::iterator ds=m_segments.begin();
       ds!=m_segments.end();ds++) (*ds)->Reset();
  if (!m_hitcells.empty()) m_hitcells.clear();
}

bool Had_Calorimeter::Fill(const double E,const double eta,
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


