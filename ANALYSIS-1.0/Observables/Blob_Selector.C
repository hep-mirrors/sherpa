#include "Primitive_Observable_Base.H"

#include "Primitive_Analysis.H"
#include "MyStrStream.H"

using namespace ANALYSIS;

class Blob_Selector: public Primitive_Observable_Base {
private:
  ATOOLS::btp::code m_type;
  std::string m_outlist;
  int m_mode;
public:
  Blob_Selector(const int type,const int mode,const std::string &outlist):
    m_type((ATOOLS::btp::code)type), m_outlist(outlist), m_mode(mode) 
  {
    m_splitt_flag=false;
  }
  
  Primitive_Observable_Base *Copy() const 
  {
    return new Blob_Selector(m_type,m_mode,m_outlist);
  }

  void Evaluate(const ATOOLS::Blob_List &blobs,double value,int ncount)
  {
    ATOOLS::Particle_List *outlist = new ATOOLS::Particle_List();
    p_ana->AddParticleList(m_outlist,outlist);
    for (ATOOLS::Blob_List::const_iterator bit=blobs.begin();
	 bit!=blobs.end();++bit) {
      if ((*bit)->Type()==m_type) {
	if (m_mode>1) {
	  for (int i=0;i<(*bit)->NInP();++i) {
	    const ATOOLS::Particle *cur=(*bit)->ConstInParticle(i);
	    bool present=false;
	    for (ATOOLS::Particle_List::const_iterator pit=outlist->begin();
		 pit!=outlist->end();++pit) 
	      if (cur==*pit) {
		present=true;
		break;
	      }
	    if (!present) outlist->push_back(new ATOOLS::Particle(*cur));
	  }
	}
	for (int i=0;i<(*bit)->NOutP();++i) {
	  const ATOOLS::Particle *cur=(*bit)->ConstOutParticle(i);
	  if (cur->DecayBlob()!=NULL) if (m_mode<1) continue;
	  bool present=false;
	  for (ATOOLS::Particle_List::const_iterator pit=outlist->begin();
	       pit!=outlist->end();++pit) 
	    if (cur==*pit) {
	      present=true;
	      break;
	    }
	  if (!present) outlist->push_back(new ATOOLS::Particle(*cur));
	}
      }
    }
  }

};

DECLARE_GETTER(Blob_Selector_Getter,"BlobSel",
	       Primitive_Observable_Base,String_Matrix);

Primitive_Observable_Base *
Blob_Selector_Getter::operator()(const String_Matrix &parameters) const
{ 
  if (parameters.size()<1) return NULL;
  if (parameters[0].size()<2) return NULL;
  std::string outlist=parameters[0].size()>2?parameters[0][2]:"Analysed";
  return new Blob_Selector(ATOOLS::ToType<int>(parameters[0][0]),
			   ATOOLS::ToType<int>(parameters[0][1]),outlist); 
}

void Blob_Selector_Getter::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"type mode [outlist]"; 
}

