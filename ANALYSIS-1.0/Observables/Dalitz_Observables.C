#include "Primitive_Observable_Base.H"

#include "Primitive_Analysis.H"
#include "MyStrStream.H"

using namespace ANALYSIS;
using namespace ATOOLS;

class Dalitz_Observable_Base: public Primitive_Observable_Base {
private:

  Flavour m_inflav, m_outflavs[3];
  size_t m_bins;
  double m_min, m_max;
  int m_type;

  TH2D *p_histogram;

public:

  Dalitz_Observable_Base(const Flavour &in,const Flavour &out1,const Flavour &out2,const Flavour &out3,
			 const size_t &bins,const double &min,const double &max,const int &type):
    m_inflav(in), m_bins(bins), m_min(min), m_max(max), m_type(type) 
  {
    m_outflavs[0]=out1;
    m_outflavs[1]=out2;
    m_outflavs[2]=out3;
    m_splitt_flag=false;
    std::string id(ATOOLS::ToString(this));
    p_histogram=new TH2D(id.c_str(),(m_inflav.IDName()+std::string("_")+m_outflav[0].IDName()+
				     m_outflav[1].IDName()+m_outflav[2].IDName()+
				     std::string("_Dalitz")).c_str(),
			 m_bins,m_min,m_max,m_bins,m_min,m_max);
    (*MYROOT::myroot)(p_histogram,id);
  }
  
  void Evaluate(const ATOOLS::Blob_List &blobs,double value,int ncount);

  virtual void Evaluate(const Vec4D &pin,const Vec4D* pout,const double &weight,const size_t &ncount) = 0;

};// end of class

void Dalitz_Observable_Base::Evaluate(const ATOOLS::Blob_List &blobs,double value,int ncount)
{
  for (ATOOLS::Blob_List::const_iterator bit=blobs.begin();
       bit!=blobs.end();++bit) {
    if ((*bit)-NInP()==1 && (*bit)->NOutP()==3) {
      if ((*bit)->ConstInParticle(0)->Flav()==m_inflav) {
	bool cont(false);
	Vec4D pout[3];
	for (int i=0;i<3;++i) {
	  const ATOOLS::Particle *cur=(*bit)->ConstOutParticle(i);
	  if (cur->Flav()==m_outflav[0] && pout[0]==Vec4D()) pout[0]=cur->Momentum();
	  else if (cur->Flav()==m_outflav[1] && pout[1]==Vec4D()) pout[1]=cur->Momentum();
	  else if (cur->Flav()==m_outflav[2] && pout[2]==Vec4D()) pout[2]=cur->Momentum();
	  else {
	    cont=true;
	    break;
	  }
	  if (!cont) {
	    Evaluate((*bit)->ConstInParticle(0)->Momentum(),pout,value,ncount);
	  }
	}
      }
    }
  }
}

class Dalitz: public Dalitz_Observable_Base {
public:
  Dalitz(const Flavour &in,const Flavour &out1,const Flavour &out2,const Flavour &out3,
	 const size_t &bins,const double &min,const double &max,const int &type):
    Dalitz_Observable_Base(in,out1,out2,out3,bins,min,max,type) {}
  
  Primitive_Observable_Base *Copy() const 
  {
    return new Dalitz_Observable_Base(m_inflav,m_outflavs[0],m_outflavs[1],m_outflavs[2],m_bins,m_min,m_max,m_type);
  }

  void Evaluate(const Vec4D* pout);

};// end of class

void Dalitz::Evaluate(const Vec4D &pin,const Vec4D* pout,const double &weight,const size_t &ncount)
{
  double s1((pin-pout[0]).Abs2()), s2((pin-pout[1]).Abs2()), s3((pin-pout[2]).Abs2()); 
  double s0((s1+s2+s3)/3.0);
  p_histogram->Fill((s1-s2)/s0,(s3-s0)/s0,weight);
  for (size_t i(1);i<ncount;++i) p_histogram->Fill(m_min,m_min,0.0);
}

DECLARE_GETTER(Dalitz_Observable_Base_Getter,"Dalitz",
	       Primitive_Observable_Base,String_Matrix);

Primitive_Observable_Base *
Dalitz_Observable_Base_Getter::operator()(const String_Matrix &parameters) const
{ 
  if (parameters.size()<1) return NULL;
  if (parameters[0].size()<8) return NULL;
  Flavour in((kf::code)ToType<int>(parameters[0][0]));
  Flavour out1((kf::code)ToType<int>(parameters[0][1]));
  Flavour out2((kf::code)ToType<int>(parameters[0][2]));
  Flavour out3((kf::code)ToType<int>(parameters[0][3]));
  return new Dalitz(in,out1,out2,out3,ToType<int>(parameters[4]),
		    ToType<double>(parameters[0][5]),
		    ToType<double>(parameters[0][6]),
		    ToType<int>(parameters[0][7])); 
}

void Dalitz_Observable_Base_Getter::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"inflav outflav1 outflav2 outflav3 bins min max type"; 
}

