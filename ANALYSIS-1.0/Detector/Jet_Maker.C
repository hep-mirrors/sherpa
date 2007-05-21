#include "Jet_Maker.H"
#include "Primitive_Analysis.H"
#include "Detector.H"
#include "Message.H"
#include "MyStrStream.H"
#include "Exception.H"

using namespace ANALYSIS;
using namespace ATOOLS;

DECLARE_GETTER(Jet_Maker_Getter,"Jet_Maker",
	       Analysis_Object,Argument_Matrix);

Analysis_Object *
Jet_Maker_Getter::operator()(const Argument_Matrix &parameters) const
{			
  if (parameters.size()<1) return NULL;
  //if (parameters.size()==1) abort(); // For read-in of, like 'ATLAS'


  std::string mode("ET_UP");
  Jet_Maker * maker = new Jet_Maker(parameters(),mode);
  double Ecut,Estart,R;

  for (size_t i=0;i<parameters.size();++i) {
    const std::vector<std::string> &cur=parameters[i];
    if (cur.size()<2) continue;
    else if (cur[0]=="Simple_Cone") {
      Ecut   = ATOOLS::ToType<double>(cur[1]);
      Estart = ATOOLS::ToType<double>(cur[2]);
      R      = ATOOLS::ToType<double>(cur[3]);
      maker->SetSimpleCone(Ecut,Estart,R);
    }
    else if (cur[0]=="ECorrection") {
      maker->SetECorrection(ATOOLS::ToType<double>(cur[1]));
    }
  }
  return maker;
}									

void Jet_Maker_Getter::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"Keyword=Simple_Cone Parameters=ETcut ETstart R"<<std::endl; 
}


Jet_Maker::Jet_Maker(Primitive_Analysis * ana,const std::string mode) :
  Object_Definition_Base(ana,"JetMaker",mode), 
  p_simplecone(NULL), m_jetmode(0)
{
  GetElements();
}

Jet_Maker::~Jet_Maker() {
  if (p_simplecone) { delete p_simplecone; p_simplecone = NULL; }
}

void Jet_Maker::SetSimpleCone(const double Ecut,const double Estart,const double R) {
  p_simplecone = new Simple_Cone(Ecut,Estart,R);
  p_simplecone->SetCalorimeters(p_HCal,p_ECal);
  p_simplecone->SetMuonChambers(p_chambers);
  m_jetmode = 1;
}


void Jet_Maker::ReconstructObjects(ATOOLS::Particle_List * plist,ATOOLS::Vec4D & METvector) {
  switch (m_jetmode) {
  case 0: return;
  case 1: 
    if (!p_simplecone->ConstructJets(&m_objects)) return;
    break;
  }

  DropUsedCells();

  Particle * part;
  int iob(0);
  while (!m_objects.empty()) {
    part = m_objects.front()->CreateParticle();
    delete m_objects.front();
    m_objects.pop_front();
    plist->push_back(part);
    METvector -= part->Momentum(); 
    iob++;
  }
}


void Jet_Maker::ReconstructSimpleCones() {}

void Jet_Maker::CorrectEnergies() {}
