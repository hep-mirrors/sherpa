#include "AHADIC++/Decays/Cluster_Part.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Random.H"


using namespace AHADIC;
using namespace ATOOLS;
using namespace std;

Cluster_Part::Cluster_Part(bool ana) :
  m_ana(ana), m_fails(0), m_att(0)
{ 
  p_csplitter = new Cluster_Splitter();
}

Cluster_Part::~Cluster_Part()
{
  if (p_csplitter) delete p_csplitter;
}

bool Cluster_Part::TestDecay(Cluster * const cluster)
{
  m_att++;
  if (p_csplitter && !(*p_csplitter)(cluster)) {
    m_fails++;
    msg_Out()<<METHOD<<" fails for\n"<<(*cluster)<<"\n";
    return false;
  }
  return true;
}
