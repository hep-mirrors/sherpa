#include "AddOns/Analysis/Observables/Four_Particle_Observables.H"
#include "AddOns/Analysis/Main/Primitive_Analysis.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"

using namespace ANALYSIS;
using namespace ATOOLS;

template <class Class>
Primitive_Observable_Base *GetObservable(const Argument_Matrix &parameters)
{									
  if (parameters.size()<1) return NULL;
  if (parameters.size()==1) {
    if (parameters[0].size()<8) return NULL;
    std::vector<ATOOLS::Flavour> f(4);
    for (short unsigned int i=0;i<4;++i) {
      int kf=ATOOLS::ToType<int>(parameters[0][i]);
      f[i]=ATOOLS::Flavour((kf_code)abs(kf));
      if (kf<0) f[i]=f[i].Bar();
    }
    std::string list=parameters[0].size()>8?parameters[0][8]:finalstate_list;
    return new Class(f,HistogramType(parameters[0][7]),
		     ATOOLS::ToType<double>(parameters[0][4]),
		     ATOOLS::ToType<double>(parameters[0][5]),
		     ATOOLS::ToType<int>(parameters[0][6]),list);
  }
  else if (parameters.size()<8) return NULL;
  double min=0.0, max=1.0;
  size_t bins=100;
  std::vector<ATOOLS::Flavour> f(4);
  std::string list=finalstate_list, scale="Lin";
  for (size_t i=0;i<parameters.size();++i) {
    if (parameters[i].size()<2) continue;
    for (short unsigned int j=0;j<4;++j) {
      if (parameters[i][0]==std::string("FLAV")+ATOOLS::ToString(j+1)) {
        int kf=ATOOLS::ToType<int>(parameters[i][1]);
        f[j]=ATOOLS::Flavour((kf_code)abs(kf));
        if (kf<0) f[j]=f[j].Bar();
      }
    }
    if (parameters[i][0]=="MIN") min=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="MAX") max=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="BINS") bins=ATOOLS::ToType<int>(parameters[i][1]);
    else if (parameters[i][0]=="SCALE") scale=parameters[i][1];
    else if (parameters[i][0]=="LIST") list=parameters[i][1];
  }
  return new Class(f,HistogramType(scale),min,max,bins,list);
}									

#define DEFINE_GETTER_METHOD(CLASS,NAME)				\
  Primitive_Observable_Base *					\
  NAME::operator()(const Argument_Matrix &parameters) const		\
  { return GetObservable<CLASS>(parameters); }

#define DEFINE_PRINT_METHOD(NAME)					\
  void NAME::PrintInfo(std::ostream &str,const size_t width) const	\
  { str<<"kf1 kf2 kf3 kf4 min max bins Lin|LinErr|Log|LogErr [list]"; }

#define DEFINE_OBSERVABLE_GETTER(CLASS,NAME,TAG)			\
  DECLARE_GETTER(NAME,TAG,Primitive_Observable_Base,Argument_Matrix);	\
  DEFINE_GETTER_METHOD(CLASS,NAME)					\
  DEFINE_PRINT_METHOD(NAME)

using namespace ATOOLS;
using namespace std;

Four_Particle_Observable_Base::Four_Particle_Observable_Base
(const std::vector<Flavour>& flavs, int type, double xmin, double xmax,
 int nbins, const std::string& listname, const std::string& name)
  : Primitive_Observable_Base(type,xmin,xmax,nbins), f_special(false) {

  if(flavs.size()<4) {
    msg_Error()<<"Error in Four_Particle_Observable_Base:"<<std::endl
	       <<"   No four flavours specified, try to copy flavours."
	       <<std::endl;
  }
  MyStrStream str;
  str<<name<<flavs[0]<<flavs[1]<<flavs[2]<<flavs[3]<<".dat";
  str>>m_name;
  Flavour fl;
  for(size_t i=0; i<4; i++) {
    if(i<flavs.size()) fl=flavs[i];
    m_flavs.push_back(fl);
  }
  m_listname = listname;
  m_blobtype = std::string("");
  m_blobdisc = false;
  if(xmin>=0.0) f_special=true;

}

void Four_Particle_Observable_Base::Evaluate(double value, double weight,
					     double ncount) {
  p_histo->Insert(value,weight,ncount); 
}

 
void Four_Particle_Observable_Base::Evaluate(int nout, const Vec4D* moms,
					     const Flavour* flavs,
					     double weight, double ncount) 
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


void Four_Particle_Observable_Base::Evaluate(const Particle_List& plist,
					     double weight, double ncount) {
  for(Particle_List::const_iterator plit1=plist.begin();
      plit1!=plist.end(); ++plit1) {
    if((*plit1)->Flav()==m_flavs[0]) {
      for(Particle_List::const_iterator plit2=plist.begin();
	  plit2!=plist.end(); ++plit2) {
	if((*plit2)->Flav()==m_flavs[1] && plit1!=plit2) {
	  for(Particle_List::const_iterator plit3=plist.begin();
	      plit3!=plist.end(); ++plit3) {
	    if((*plit3)->Flav()==m_flavs[2] && plit3!=plit2 && plit3!=plit1) {
	      for(Particle_List::const_iterator plit4=plist.begin();
		  plit4!=plist.end(); ++plit4) {
		if((*plit4)->Flav()==m_flavs[3] &&
		   plit4!=plit3 && plit4!=plit2 && plit4!=plit1) {
		  Evaluate((*plit1)->Momentum(),(*plit2)->Momentum(),
			   (*plit3)->Momentum(),(*plit4)->Momentum(),
			   weight, ncount);
		  return;
		}
	      }
	    }
	  }
	}
      }
    }
  }
  p_histo->Insert(0.0,0.0,ncount);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_OBSERVABLE_GETTER(Four_Particle_PlaneAngle,
			 Four_Particle_PlaneAngle_Getter,"PlaneAngle")

void Four_Particle_PlaneAngle::Evaluate(const Vec4D & mom1,const Vec4D & mom2,
					const Vec4D & mom3,const Vec4D & mom4,
					double weight, double ncount)
{
  Vec3D normal1 = cross(Vec3D(mom1),Vec3D(mom2));
  Vec3D normal2 = cross(Vec3D(mom3),Vec3D(mom4));
  double costh  = (normal1*normal2)/(normal1.Abs()*normal2.Abs()); 
  p_histo->Insert(costh,weight,ncount); 
}
 
Four_Particle_PlaneAngle::Four_Particle_PlaneAngle(const std::vector<Flavour> & flavs,
						   int type,double xmin,double xmax,int nbins,
						   const std::string& listname) :
  Four_Particle_Observable_Base(flavs,type,xmin,xmax,nbins,listname,"NRAngle") { }

Primitive_Observable_Base * Four_Particle_PlaneAngle::Copy() const
{
  return new Four_Particle_PlaneAngle(m_flavs,m_type,m_xmin,m_xmax,m_nbins,
				      m_listname);
}

//=============================================================================

DEFINE_OBSERVABLE_GETTER(Four_Particle_PT,
			 Four_Particle_PT_Getter,"PT4")

void Four_Particle_PT::Evaluate(const Vec4D& mom1,const Vec4D& mom2,
				const Vec4D& mom3,const Vec4D& mom4,
				double weight, double ncount)
{
  double pt = sqrt(sqr(mom1[1]+mom2[1]+mom3[1]+mom4[1]) +
		   sqr(mom1[2]+mom2[2]+mom3[2]+mom4[2]));
  p_histo->Insert(pt,weight,ncount);
}
 
Four_Particle_PT::Four_Particle_PT(const std::vector<Flavour>& flavs,
				   int type,double xmin,double xmax,int nbins,
				   const std::string & listname)
  : Four_Particle_Observable_Base(flavs,type,xmin,xmax,nbins,listname,"PT") {}

Primitive_Observable_Base* Four_Particle_PT::Copy() const
{
  return new Four_Particle_PT(m_flavs,m_type,m_xmin,m_xmax,m_nbins,m_listname);
}

//=============================================================================

DEFINE_OBSERVABLE_GETTER(Two_Partonpair_PTdiff,
			 Two_Partonpair_PTdiff_Getter,"PTdiff4")

void Two_Partonpair_PTdiff::Evaluate(const Vec4D& mom1,const Vec4D& mom2,
				     const Vec4D& mom3,const Vec4D& mom4,
				     double weight, double ncount) {
  Vec4D vecA(mom1); vecA+=mom2;
  Vec4D vecB(mom3); vecB+=mom4;
  double ptdiff=sqrt(sqr(vecA[1])+sqr(vecA[2]));
  ptdiff-=sqrt(sqr(vecB[1])+sqr(vecB[2]));
  if(f_special) ptdiff=dabs(ptdiff);
  p_histo->Insert(ptdiff, weight, ncount);

  //std::cout<<f_special<<" ptdiff ---> "<<ptdiff<<"\n";
}
 
Two_Partonpair_PTdiff::Two_Partonpair_PTdiff(const std::vector<Flavour>& flavs,
					     int type,
					     double xmin, double xmax,
					     int nbins,
					     const std::string & listname)
  : Four_Particle_Observable_Base(flavs,type,xmin,xmax,
				  nbins,listname,"PTdiff") {}

Primitive_Observable_Base* Two_Partonpair_PTdiff::Copy() const {
  return new Two_Partonpair_PTdiff(m_flavs, m_type, m_xmin, m_xmax,
				   m_nbins, m_listname);
}

//=============================================================================

DEFINE_OBSERVABLE_GETTER(Two_Partonpair_Theta,
			 Two_Partonpair_Theta_Getter,"Theta4")

void Two_Partonpair_Theta::Evaluate(const Vec4D& mom1,const Vec4D& mom2,
				    const Vec4D& mom3,const Vec4D& mom4,
				    double weight, double ncount) {
  Vec4D vecA(mom1); vecA+=mom2;
  Vec4D vecB(mom3); vecB+=mom4;
  //
  if(!f_special) {
    Vec4D plab(vecA); plab+=vecB;
    //removing the z boost effect of the considered 4 particle system
    plab[1]=plab[2]=0.0;
    if(plab.Abs2()<=0.0) {
      p_histo->Insert(-M_PI/100.0, weight, ncount);
      msg_Error()<<__PRETTY_FUNCTION__<<":\n   Warning:"
		 <<" Not able to boost the system. Insert theta=-pi/100.\n"
		 <<std::endl;
      return;
    }
    Poincare fly(plab);
    fly.Boost(vecA);
    fly.Boost(vecB);
  }
  Vec3D vec1(vecA);
  Vec3D vec2(vecB);
  double theta=acos((vec1*vec2)/(vec1.Abs()*vec2.Abs()));
  p_histo->Insert(theta, weight, ncount);

  //std::cout<<m_flavs[0]<<" : "<<mom1<<"\n";
  //std::cout<<m_flavs[1]<<" : "<<mom2<<"  :  "<<vecA<<"\n";
  //std::cout<<"----------\n";
  //std::cout<<m_flavs[2]<<" : "<<mom3<<"\n";
  //std::cout<<m_flavs[3]<<" : "<<mom4<<"  :  "<<vecB<<"\n\n";

  //if(f_special) std::cout<<" theta -------------------> "<<theta<<"\n";
  //else std::cout<<" theta boost corrected ---> "<<theta<<"\n";

}

Two_Partonpair_Theta::Two_Partonpair_Theta(const std::vector<Flavour>& flavs,
					   int type, double xmin, double xmax,
					   int nbins,
					   const std::string & listname)
  : Four_Particle_Observable_Base(flavs,type,xmin,xmax,
				  nbins,listname,"Theta") {}

Primitive_Observable_Base* Two_Partonpair_Theta::Copy() const {
  return new Two_Partonpair_Theta(m_flavs, m_type, m_xmin, m_xmax,
				  m_nbins, m_listname);
}

// ============================================================================

template <class Class>
Primitive_Observable_Base *GetObservable2(const Argument_Matrix &parameters)
{									
  if (parameters.size()<1) return NULL;
  if (parameters.size()==1) {
    if (parameters[0].size()<4) return NULL;
    std::string list=parameters[0].size()>4?parameters[0][4]:finalstate_list;
    return new Class(HistogramType(parameters[0][3]),
		     ATOOLS::ToType<double>(parameters[0][0]),
		     ATOOLS::ToType<double>(parameters[0][1]),
		     ATOOLS::ToType<int>(parameters[0][2]),list);
  }
  else if (parameters.size()<4) return NULL;
  double min=0.0, max=1.0;
  size_t bins=100;
  std::string list=finalstate_list, scale="Lin";
  for (size_t i=0;i<parameters.size();++i) {
    if (parameters[i].size()<2) continue;
    if (parameters[i][0]=="MIN") min=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="MAX") max=ATOOLS::ToType<double>(parameters[i][1]);
    else if (parameters[i][0]=="BINS") bins=ATOOLS::ToType<int>(parameters[i][1]);
    else if (parameters[i][0]=="SCALE") scale=parameters[i][1];
    else if (parameters[i][0]=="LIST") list=parameters[i][1];
  }
  return new Class(HistogramType(scale),min,max,bins,list);
}									

#define DEFINE_GETTER_METHOD2(CLASS,NAME)				\
  Primitive_Observable_Base *					\
  NAME::operator()(const Argument_Matrix &parameters) const		\
  { return GetObservable2<CLASS>(parameters); }

#define DEFINE_PRINT_METHOD2(NAME)					\
  void NAME::PrintInfo(std::ostream &str,const size_t width) const	\
  { str<<"min max bins Lin|LinErr|Log|LogErr [list]"; }

#define DEFINE_OBSERVABLE_GETTER2(CLASS,NAME,TAG)			\
  DECLARE_GETTER(NAME,TAG,Primitive_Observable_Base,Argument_Matrix);	\
  DEFINE_GETTER_METHOD2(CLASS,NAME)					\
  DEFINE_PRINT_METHOD2(NAME)


DEFINE_OBSERVABLE_GETTER2(Di_Mass,
			 Di_Mass_Getter,"DiMass")

Di_Mass::Di_Mass(unsigned int type,double xmin,double xmax,int nbins,
	       const std::string & lname) :
  Primitive_Observable_Base(type,xmin,xmax,nbins)
{
  m_listname=lname;
  m_name  = std::string("4jet_");
  if (lname!=finalstate_list) m_name=lname+std::string("_")+m_name;
  m_name += "DiMass.dat";
}

void Di_Mass::Evaluate(const ATOOLS::Blob_List & blobs,double weight, double ncount)
{
  Particle_List * pl=p_ana->GetParticleList(m_listname);

  if (pl->size()!=4) {
    p_histo->Insert(0.,0.,2*ncount);
    return;
  }

  std::vector<Vec4D> moms;
  for (Particle_List::const_iterator pit=pl->begin();pit!=pl->end();++pit) {
    moms.push_back((*pit)->Momentum());
  }

  double m1a = (moms[0]+moms[1]).Abs2();
  double m1b = (moms[2]+moms[3]).Abs2();
  double d1 = dabs(m1a-m1b);

  double m2a = (moms[0]+moms[2]).Abs2();
  double m2b = (moms[1]+moms[3]).Abs2();
  double d2 = dabs(m2a-m2b);

  double m3a = (moms[0]+moms[3]).Abs2();
  double m3b = (moms[1]+moms[2]).Abs2();
  double d3 = dabs(m3a-m3b);

  if (d1<d2 && d1<d3) {
    p_histo->Insert(sqrt(m1a),weight,ncount);
    p_histo->Insert(sqrt(m1b),weight,ncount);
  }
  else if (d2<d3) {
    p_histo->Insert(sqrt(m2a),weight,ncount);
    p_histo->Insert(sqrt(m2b),weight,ncount);
  }
  else {
    p_histo->Insert(sqrt(m3a),weight,ncount);
    p_histo->Insert(sqrt(m3b),weight,ncount);
  }
}

Primitive_Observable_Base * Di_Mass::Copy() const 
{
  return new Di_Mass(m_type,m_xmin,m_xmax,m_nbins,m_listname);
}

// ==================================================================

DEFINE_OBSERVABLE_GETTER(Four_Particle_PlaneAngleCMS,
                         Four_Particle_PlaneAngleCMS_Getter,"PlaneAngleCMS")

void Four_Particle_PlaneAngleCMS::Evaluate(const Vec4D & mom1,const Vec4D & mom2,
                                        const Vec4D & mom3,const Vec4D & mom4,
                                        double weight, double ncount)
{
  Vec4D sum = mom1+mom2+mom3+mom4;
  Poincare boost(sum);
  Vec4D p1 = boost*mom1;
  Vec4D p2 = boost*mom2;
  Vec4D p3 = boost*mom3;
  Vec4D p4 = boost*mom4;
  Vec3D normal1 = cross(Vec3D(p1),Vec3D(p3+p4));
  Vec3D normal2 = cross(Vec3D(p3),Vec3D(p3+p4));
  double costh  = (normal1*normal2)/(normal1.Abs()*normal2.Abs());
  p_histo->Insert(acos(costh),weight,ncount);
}

Four_Particle_PlaneAngleCMS::Four_Particle_PlaneAngleCMS(const std::vector<Flavour> & flavs,
                                                   int type,double xmin,double xmax,int nbins,
                                                   const std::string& listname) :
  Four_Particle_Observable_Base(flavs,type,xmin,xmax,nbins,listname,"PlaneAngleCMS") { }

Primitive_Observable_Base * Four_Particle_PlaneAngleCMS::Copy() const
{
  return new Four_Particle_PlaneAngleCMS(m_flavs,m_type,m_xmin,m_xmax,m_nbins,
                                      m_listname);
}

// ==================================================================

DEFINE_OBSERVABLE_GETTER(Four_Particle_EnergyCMS,
                         Four_Particle_EnergyCMS_Getter,"4EnergyCMS")

void Four_Particle_EnergyCMS::Evaluate(const Vec4D & mom1,const Vec4D & mom2,
                                       const Vec4D & mom3,const Vec4D & mom4,
                                       double weight, double ncount)
{
  Vec4D sum = mom1+mom2+mom3+mom4;
  Poincare boost(sum);
  Vec4D p1 = boost*mom1;
  Vec4D p2 = boost*mom2;
  Vec4D p3 = boost*mom3;
  Vec4D p4 = boost*mom4;
  double E = p1[0];
  p_histo->Insert(2.0*E/rpa->gen.Ecms(),weight,ncount);
}

Four_Particle_EnergyCMS::Four_Particle_EnergyCMS(const std::vector<Flavour> & flavs,
                                                   int type,double xmin,double xmax,int nbins,
                                                   const std::string& listname) :
  Four_Particle_Observable_Base(flavs,type,xmin,xmax,nbins,listname,"4EnergyCMS") { }

Primitive_Observable_Base * Four_Particle_EnergyCMS::Copy() const
{
  return new Four_Particle_EnergyCMS(m_flavs,m_type,m_xmin,m_xmax,m_nbins,
                                      m_listname);
}
