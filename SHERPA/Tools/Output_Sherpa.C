#include "SHERPA/Tools/Output_Sherpa.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace SHERPA;
using namespace ATOOLS;
using namespace std;

void Output_Sherpa::Header()
{
  m_outstream<<"# created by SHERPA "<<SHERPA_VERSION<<"."<<SHERPA_SUBVERSION
             <<endl;
}

void Output_Sherpa::Output(Blob_List* blobs, const double weight) 
{
  ATOOLS::Particle_Int_Map           P2I;
  ATOOLS::Particle_Int_Map::iterator P2Iiter;
  
  for (Blob_List::const_iterator blit=blobs->begin();blit!=blobs->end();++blit){
    for (int i=0;i<(*blit)->NInP();i++) {
      if (P2I.find((*blit)->InParticle(i))==P2I.end()) 
	P2I.insert(make_pair((*blit)->InParticle(i),P2I.size()+1));
    }
    for (int i=0;i<(*blit)->NOutP();i++) {
      if (P2I.find((*blit)->OutParticle(i))==P2I.end()) 
	P2I.insert(make_pair((*blit)->OutParticle(i),P2I.size()+1));
    }
  }
  m_outstream<<rpa.gen.NumberOfDicedEvents()-1<<" "<<P2I.size()<<" "
             <<blobs->size()<<" "<<weight<<endl;
  Particle * part;
  int kfc;
  for (P2Iiter=P2I.begin();P2Iiter!=P2I.end();P2Iiter++) {
    part = P2Iiter->first;
    kfc  = part->Flav().Kfcode(); if (part->Flav().IsAnti()) kfc=-kfc;
    m_outstream<<P2Iiter->second<<" "<<part->Status()<<" "<<part->Info()<<" "
               <<kfc<<" "
               <<" "<<part->Momentum()[0]<<" "<<part->Momentum()[1]
	       <<" "<<part->Momentum()[2]<<" "<<part->Momentum()[3]<<" \n";
  }
  for (Blob_List::const_iterator blit=blobs->begin();blit!=blobs->end();++blit){
    m_outstream<<(*blit)->Id()<<" "<<(*blit)->Status()<<" "
               <<(int)(*blit)->Type()<<" "<<(*blit)->TypeSpec()
	       <<" "<<(*blit)->NInP()<<" "<<(*blit)->NOutP()<<" \n"
	       <<" "<<(*blit)->Position()[0]<<" "<<(*blit)->Position()[1]
	       <<" "<<(*blit)->Position()[2]<<" "<<(*blit)->Position()[3]
               <<" \n";
    for (int i=0;i<(*blit)->NInP();i++)
      m_outstream<<P2I.find((*blit)->InParticle(i))->second<<" ";
    for (int i=0;i<(*blit)->NOutP();i++)
      m_outstream<<P2I.find((*blit)->OutParticle(i))->second<<" ";
    m_outstream<<" \n";
  }
}

void Output_Sherpa::Footer(std::string number)
{
  string newfile=m_basename+"."+number+m_ext;
  string footer = newfile.substr(newfile.rfind("/",newfile.size())+1);
  m_outstream<<footer<<std::endl;
}
