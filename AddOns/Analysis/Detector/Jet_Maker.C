#include "AddOns/Analysis/Detector/Jet_Maker.H"
#include "AddOns/Analysis/Main/Primitive_Analysis.H"
#include "AddOns/Analysis/Detector/Detector.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Exception.H"

#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#ifdef USING__ROOT
#include "ATOOLS/Math/Scaling.H"
#include "TH1D.h"
#include "TH2D.h"
#include "ATOOLS/Org/My_Root.H"
#endif 

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
      if (cur[1]=="truth") {
	maker->SetECorrection(ATOOLS::ToType<double>(cur[2]),0);
      }
      else if (cur[1]=="constant") {
	maker->SetECorrection(ATOOLS::ToType<double>(cur[2]),1);
      }
    }
    else if (cur[0]=="BTagging") {
      maker->SetBTagging(ATOOLS::ToType<double>(cur[1]),
			 ATOOLS::ToType<double>(cur[2]),ATOOLS::ToType<double>(cur[3]));
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
#ifdef USING__ROOT
  std::string name = std::string("E_Correction");
  (*MYROOT::myroot)(new TH2D(name.c_str(),name.c_str(),200,0.,1000.,20,0.5,1.5),name);
  name = std::string("ET_Correction");
  (*MYROOT::myroot)(new TH2D(name.c_str(),name.c_str(),200,0.,1000.,20,0.5,1.5),name);
#endif
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

void Jet_Maker::SetECorrection(const double val,const int m_correction) {
  switch (m_correction) {
  case 1:  
    m_inv    = val; break;
  case 0: 
  default: 
    m_spread = val; break;
  }
}

void Jet_Maker::SetBTagging(const double tagprob,const double rlight, const double rcharm) {
  p_btagging = new BTagging(tagprob,1./rlight,1./rcharm);
}

void Jet_Maker::ReconstructObjects(ATOOLS::Particle_List * plist,ATOOLS::Vec4D & METvector) {
  if (!p_btagging) p_btagging = new BTagging(0.5,1./100.,1./10.);
  switch (m_jetmode) {
  case 0: return;
  case 1: 
    if (!p_simplecone->ConstructJets(&m_objects)) return;
    break;
  }

  msg_Debugging()<<METHOD<<" :"<<std::endl;
  Particle * part;
  while (!m_objects.empty()) {
    m_objects.front()->SetIncludeTracks(true);
    m_objects.front()->Update();
    switch (m_correction) {
    case 1:
      m_objects.front()->CorrectE(m_inv); break;
    case 0:
    default:
      m_objects.front()->CorrectTruth(m_spread); break;
    }
    part = m_objects.front()->CreateParticle();
    msg_Debugging()<<"   Found jet : "<<m_objects.front()->Mom()
    	     <<"/"<<part->Momentum()<<" with "<<m_objects.front()->GetCells().size()
    	     <<"/"<<m_objects.front()->GetTracks().size()<<std::endl;
    delete m_objects.front();
    m_objects.pop_front();
    plist->push_back(part);
    METvector -= part->Momentum(); 
  }

  DropUsedCells();
}


void Jet_Maker::CorrectEnergies() {
  for (ObjectListIterator olit=m_objects.begin();olit!=m_objects.end();olit++) {
    switch (m_correction) {
    case 1:
      (*olit)->CorrectE(m_inv);
      break;
    case 0:
    default:
      (*olit)->CorrectTruth(m_spread);
#ifdef USING__ROOT
      std::string name = std::string("E_Correction");
      ((TH2D*)(*MYROOT::myroot)[name])->Fill((*olit)->TrueMom()[0],(*olit)->E_Correction(),1.);
      name = std::string("ET_Correction");
      // Vec4D has no member ET.
      //((TH2D*)(*MYROOT::myroot)[name])->Fill((*olit)->TrueMom().ET(),(*olit)->ET_Correction(),1.);
#endif
      break;
    }
  }  
}

void Jet_Maker::IdentifyBs() {
  for (ObjectListIterator olit=m_objects.begin();olit!=m_objects.end();olit ++) {
    if (p_btagging->Tag((*olit))) (*olit)->SetFlavour(Flavour(kf_bjet));
  }
}

