#include "Four_Particle_Observables.H"

using namespace ANALYSIS;

#include "MyStrStream.H"

template <class Class>
Primitive_Observable_Base *const GetObservable(const String_Matrix &parameters)
{									
  if (parameters.size()<1) return NULL;
  if (parameters.size()==1) {
    if (parameters[0].size()<8) return NULL;
    std::vector<ATOOLS::Flavour> f(4);
    for (short unsigned int i=0;i<4;++i) {
      int kf=ATOOLS::ToType<int>(parameters[0][i]);
      f[i]=ATOOLS::Flavour((ATOOLS::kf::code)abs(kf));
      if (kf<0) f[i]=f[i].Bar();
    }
    std::string list=parameters[0].size()>8?parameters[0][8]:"Analysed";
    return new Class(f,10*(int)(parameters[0][7]=="Log"),
		     ATOOLS::ToType<double>(parameters[0][4]),
		     ATOOLS::ToType<double>(parameters[0][5]),
		     ATOOLS::ToType<int>(parameters[0][6]),list);
  }
  else if (parameters.size()<8) return NULL;
  double min=0.0, max=1.0;
  size_t bins=100;
  std::vector<ATOOLS::Flavour> f(4);
  std::string list="Analysed", scale="Lin";
  for (size_t i=0;i<parameters.size();++i) {
    if (parameters[i].size()<2) continue;
    for (short unsigned int j=0;j<4;++j) {
      if (parameters[i][0]==std::string("FLAV")+ATOOLS::ToString(j+1)) {
	int kf=ATOOLS::ToType<int>(parameters[i][1]);
	f[j]=ATOOLS::Flavour((ATOOLS::kf::code)abs(kf));
	if (kf<0) f[j]=f[j].Bar();
      }
    }
    if (parameters[i][0]=="MIN") min=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="MAX") max=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="BINS") bins=ATOOLS::ToType<int>(parameters[i][1]);
    else if (parameters[i][0]=="SCALE") scale=parameters[i][1];
    else if (parameters[i][0]=="LIST") list=parameters[i][1];
  }
  return new Class(f,(scale=="Log")*10,min,max,bins,list);
}									

#define DEFINE_GETTER_METHOD(CLASS,NAME)				\
  Primitive_Observable_Base *					\
  NAME::operator()(const String_Matrix &parameters) const		\
  { return GetObservable<CLASS>(parameters); }

#define DEFINE_PRINT_METHOD(NAME)					\
  void NAME::PrintInfo(std::ostream &str,const size_t width) const	\
  { str<<"kf1 kf2 kf3 kf4 min max bins Lin|Log [list]"; }

#define DEFINE_OBSERVABLE_GETTER(CLASS,NAME,TAG)			\
  DECLARE_GETTER(NAME,TAG,Primitive_Observable_Base,String_Matrix);	\
  DEFINE_GETTER_METHOD(CLASS,NAME);					\
  DEFINE_PRINT_METHOD(NAME)

using namespace ATOOLS;
using namespace std;

Four_Particle_Observable_Base::Four_Particle_Observable_Base(const std::vector<Flavour> & flavs,
							     int type,double xmin,double xmax,
							     int nbins,const std::string & name) :
  Primitive_Observable_Base(type,xmin,xmax,nbins,NULL)
{
  if (flavs.size()<4) {
    msg.Error()<<"Error in Four_Particle_Observable_Base:"<<std::endl
	       <<"   No four flavours specified, try to copy flavours."<<std::endl;
  }
  std::string help = std::string("");
  Flavour fl;
  for (size_t i=0;i<4;i++) {
    if (i<flavs.size()) fl=flavs[i];
    m_flavs.push_back(fl);
    help += fl.Name();
    if (i==1) help+=std::string("--");
  }
  m_name     = name + std::string(".dat");
  m_blobtype = std::string("");
  m_blobdisc = false;
}

void Four_Particle_Observable_Base::Evaluate(double value,double weight, int ncount) {
  p_histo->Insert(value,weight,ncount); 
}

 
void Four_Particle_Observable_Base::Evaluate(int nout,const Vec4D * moms,const Flavour * flavs,
					     double weight, int ncount) 
{
  for (int i=0;i<nout;i++) { 
    if (flavs[i]==m_flavs[0]) {
      for (int j=0;j<nout;j++) { 
	if (flavs[j]==m_flavs[1] && i!=j) {
	  for (int k=0;k<nout;k++) { 
	    if (flavs[k]==m_flavs[2] && k!=j && k!=i) {
	      for (int l=0;l<nout;l++) { 
		if (flavs[l]==m_flavs[3] && l!=k && l!=j && l!=i) {
		  Evaluate(moms[i],moms[j],moms[k],moms[l],weight,ncount);
		}
	      }
	    }
	  }
	} 
      }
    }
  }
}


void Four_Particle_Observable_Base::Evaluate(const Particle_List & plist,double weight,int ncount)
{
  for (Particle_Const_Iterator plit1=plist.begin();plit1!=plist.end();++plit1) {
    if ((*plit1)->Flav()==m_flavs[0]) {
      for (Particle_Const_Iterator plit2=plist.begin();plit2!=plist.end();++plit2) {
	if ((*plit2)->Flav()==m_flavs[1] && plit1!=plit2) {
	  for (Particle_Const_Iterator plit3=plist.begin();plit3!=plist.end();++plit3) {
	    if ((*plit3)->Flav()==m_flavs[2] && plit3!=plit2 && plit3!=plit1) {
	      for (Particle_Const_Iterator plit4=plist.begin();plit4!=plist.end();++plit4) {
		if ((*plit4)->Flav()==m_flavs[3] && plit4!=plit3 && plit4!=plit2 && plit4!=plit1) {
		  Evaluate((*plit1)->Momentum(),(*plit2)->Momentum(),
			   (*plit3)->Momentum(),(*plit4)->Momentum(),weight,ncount);
		}
	      }
	    }
	  }
	}
      }
    }
  }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_OBSERVABLE_GETTER(Four_Particle_PlaneAngle,
			 Four_Particle_PlaneAngle_Getter,"PlaneAngle");

void Four_Particle_PlaneAngle::Evaluate(const Vec4D & mom1,const Vec4D & mom2,
					const Vec4D & mom3,const Vec4D & mom4,
					double weight, int ncount)
{
  Vec3D normal1 = cross(Vec3D(mom1),Vec3D(mom2));
  Vec3D normal2 = cross(Vec3D(mom3),Vec3D(mom4));
  double costh  = (normal1*normal2)/(normal1.Abs()*normal2.Abs()); 
  p_histo->Insert(costh,weight,ncount); 
}
 
Four_Particle_PlaneAngle::Four_Particle_PlaneAngle(const std::vector<Flavour> & flavs,
						   int type,double xmin,double xmax,int nbins,
						   const std::string & name) :
  Four_Particle_Observable_Base(flavs,type,xmin,xmax,nbins,name) { }

Primitive_Observable_Base * Four_Particle_PlaneAngle::Copy() const
{
  return new Four_Particle_PlaneAngle(m_flavs,m_type,m_xmin,m_xmax,m_nbins,m_name);
}

