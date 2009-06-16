#include "AddOns/Analysis/Observables/Three_Particle_Observables.H"
#include "AddOns/Analysis/Main/Primitive_Analysis.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace ANALYSIS;

#include "ATOOLS/Org/MyStrStream.H"

template <class Class>
Primitive_Observable_Base *const GetObservable(const Argument_Matrix &parameters)
{									
  if (parameters.size()<1) return NULL;
  if (parameters.size()==1) {
    if (parameters[0].size()<7) return NULL;
    ATOOLS::Flavour f[3];
    for (short unsigned int i=0;i<3;++i) {
      int kf=ATOOLS::ToType<int>(parameters[0][i]);
      f[i]=ATOOLS::Flavour((kf_code)abs(kf));
      if (kf<0) f[i]=f[i].Bar();
    }
    std::string list=parameters[0].size()>7?parameters[0][7]:finalstate_list;
    return new Class(f[0],f[1],f[2],HistogramType(parameters[0][6]),
		     ATOOLS::ToType<double>(parameters[0][3]),
		     ATOOLS::ToType<double>(parameters[0][4]),
		     ATOOLS::ToType<int>(parameters[0][5]),list);
  }
  else if (parameters.size()<7) return NULL;
  double min=0.0, max=1.0;
  size_t bins=100;
  ATOOLS::Flavour f[3];
  std::string list=finalstate_list, scale="Lin";
  for (size_t i=0;i<parameters.size();++i) {
    if (parameters[i].size()<3) continue;
    for (short unsigned int j=0;j<3;++j) {
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
  return new Class(f[0],f[1],f[2],HistogramType(scale),min,max,bins,list);
}									

#define DEFINE_GETTER_METHOD(CLASS,NAME)				\
  Primitive_Observable_Base *					\
  NAME::operator()(const Argument_Matrix &parameters) const		\
  { return GetObservable<CLASS>(parameters); }

#define DEFINE_PRINT_METHOD(NAME)					\
  void NAME::PrintInfo(std::ostream &str,const size_t width) const	\
  { str<<"kf1 kf2 kf3 min max bins Lin|LinErr|Log|LogErr [list]"; }

#define DEFINE_OBSERVABLE_GETTER(CLASS,NAME,TAG)			\
  DECLARE_GETTER(NAME,TAG,Primitive_Observable_Base,Argument_Matrix);	\
  DEFINE_GETTER_METHOD(CLASS,NAME)					\
  DEFINE_PRINT_METHOD(NAME)

using namespace ATOOLS;
using namespace std;

Three_Particle_Observable_Base::Three_Particle_Observable_Base(const Flavour & flav1,const Flavour & flav2,
							       const Flavour & flav3,int type,double xmin,
							       double xmax,int nbins,const std::string & listname,
							       const std::string & name) :
  Primitive_Observable_Base(type,xmin,xmax,nbins), 
  m_flav1(flav1), m_flav2(flav2), m_flav3(flav3)
{
  m_listname=listname;
  MyStrStream str;
  str<<name<<m_flav1<<m_flav2<<m_flav3<<".dat";
  str>>m_name;
  m_blobtype = std::string("");
  m_blobdisc = false;
}

void Three_Particle_Observable_Base::Evaluate(const Particle_List & plist,double weight, double ncount)
{
  for (Particle_List::const_iterator plit1=plist.begin();plit1!=plist.end();++plit1) {
    if ((*plit1)->Flav()==m_flav1) {
      for (Particle_List::const_iterator plit2=plist.begin();plit2!=plist.end();++plit2) {
	if ((*plit2)->Flav()==m_flav2 && plit1!=plit2) {
	  for (Particle_List::const_iterator plit3=plist.begin();plit3!=plist.end();++plit3) {
	    if ((*plit3)->Flav()==m_flav3 && plit1!=plit3 && plit2!=plit3) {
	      Evaluate((*plit1)->Momentum(),(*plit2)->Momentum(),(*plit3)->Momentum(),weight,ncount);
	      return;
	    }
	  }
	}
      }
    }
  }
  p_histo->Insert(0.0,0.0,ncount);
}

//=============================================================================

DEFINE_OBSERVABLE_GETTER(Three_Particle_PT,
			 Three_Particle_PT_Getter,"PT3")

Three_Particle_PT::Three_Particle_PT(const Flavour& flav1, const Flavour& flav2,
				     const Flavour& flav3, int type,
				     double xmin, double xmax, int nbins,
				     const std::string & listname) :
  Three_Particle_Observable_Base(flav1,flav2,flav3,type,xmin,xmax,nbins,
				 listname,"PT") {}

void Three_Particle_PT::Evaluate(const Vec4D& mom1, const Vec4D& mom2,
				 const Vec4D& mom3, double weight, double ncount) {
  double pt = sqrt(sqr(mom1[1]+mom2[1]+mom3[1]) + sqr(mom1[2]+mom2[2]+mom3[2]));
  p_histo->Insert(pt,weight,ncount);
}

Primitive_Observable_Base* Three_Particle_PT::Copy() const {
  return new Three_Particle_PT(m_flav1,m_flav2,m_flav3,m_type,
			       m_xmin,m_xmax,m_nbins,m_listname);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_OBSERVABLE_GETTER(Three_Particle_Y,Three_Particle_Y_Getter,"Y3")

Three_Particle_Y::Three_Particle_Y(const Flavour& flav1, const Flavour& flav2,
                                   const Flavour& flav3, int type, double xmin,
                                   double xmax, int nbins,
				   const std::string& listname) :
  Three_Particle_Observable_Base(flav1,flav2,flav3,type,xmin,xmax,nbins,
				 listname,"Y") {}

void Three_Particle_Y::Evaluate(const Vec4D& mom1, const Vec4D& mom2,
				const Vec4D & mom3, double weight, double ncount) {
  double y = (mom1+mom2+mom3).Y();
  p_histo->Insert(y,weight,ncount);
}

Primitive_Observable_Base * Three_Particle_Y::Copy() const {
  return new Three_Particle_Y(m_flav1,m_flav2,m_flav3,m_type,
			      m_xmin,m_xmax,m_nbins,m_listname);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_OBSERVABLE_GETTER(Three_Particle_DEta,Three_Particle_DEta_Getter,"DEta3")

Three_Particle_DEta::Three_Particle_DEta(const Flavour & flav1,const Flavour & flav2,
					 const Flavour & flav3,int type,double xmin,
					 double xmax,int nbins,const std::string & listname) :
  Three_Particle_Observable_Base(flav1,flav2,flav3,type,xmin,xmax,nbins,listname,"deta3") 
{ 
}


void Three_Particle_DEta::Evaluate(const Vec4D & mom1,const Vec4D & mom2,const Vec4D & mom3,double weight, double ncount) 
{    
  Vec4D mother = mom1+mom2;
  double deta = abs((mother.Eta()-mom3.Eta()));
  p_histo->Insert(deta,weight,ncount); 
} 

Primitive_Observable_Base * Three_Particle_DEta::Copy() const 
{
    return new Three_Particle_DEta(m_flav1,m_flav2,m_flav3,m_type,m_xmin,m_xmax,m_nbins,m_listname);
}
    
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_OBSERVABLE_GETTER(Three_Particle_DPhi,Three_Particle_DPhi_Getter,"DPhi3")

Three_Particle_DPhi::Three_Particle_DPhi(const Flavour & flav1,const Flavour & flav2,
					 const Flavour & flav3,int type,double xmin,
					 double xmax,int nbins,const std::string & listname) :
  Three_Particle_Observable_Base(flav1,flav2,flav3,type,xmin,xmax,nbins,listname,"dphi3") 
{ 
}


void Three_Particle_DPhi::Evaluate(const Vec4D & mom1,const Vec4D & mom2,const Vec4D & mom3,double weight, double ncount) 
{ 
  
  Vec4D mother = mom1+mom2;

  double pt1=sqrt(mother[1]*mother[1]+mother[2]*mother[2]);
  double pt2=sqrt(mom3[1]*mom3[1]+mom3[2]*mom3[2]);
  double dphi=acos((mother[1]*mom3[1]+mother[2]*mom3[2])/(pt1*pt2));
  p_histo->Insert(dphi,weight,ncount); 
} 

Primitive_Observable_Base * Three_Particle_DPhi::Copy() const 
{
    return new Three_Particle_DPhi(m_flav1,m_flav2,m_flav3,m_type,m_xmin,m_xmax,m_nbins,m_listname);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_OBSERVABLE_GETTER(Three_Particle_DR,Three_Particle_DR_Getter,"DR3")

Three_Particle_DR::Three_Particle_DR(const Flavour & flav1,const Flavour & flav2,
				     const Flavour & flav3,int type, double xmin, 
				     double xmax, int nbins,const std::string & listname) :
  Three_Particle_Observable_Base(flav1,flav2,flav3,type,xmin,xmax,nbins,listname,"dr3") 
{
 }


void Three_Particle_DR::Evaluate(const Vec4D & mom1,const Vec4D & mom2,const Vec4D & mom3,double weight, double ncount) 
{ 
  
  Vec4D mother = mom1+mom2;
  
  double pt1=sqrt(mother[1]*mother[1]+mother[2]*mother[2]);
  double pt2=sqrt(mom3[1]*mom3[1]+mom3[2]*mom3[2]);
  double dphi=acos((mother[1]*mom3[1]+mother[2]*mom3[2])/(pt1*pt2));
  double c1=mother[3]/Vec3D(mother).Abs();
  double c2=mom3[3]/Vec3D(mom3).Abs();
  double deta=0.5 *log( (1 + c1)*(1 - c2)/((1-c1)*(1+c2)));
  double dr= sqrt(sqr(deta) + sqr(dphi)); 
  p_histo->Insert(dr,weight,ncount); 
} 

Primitive_Observable_Base * Three_Particle_DR::Copy() const 
{
  return new Three_Particle_DR(m_flav1,m_flav2,m_flav3,m_type,m_xmin,m_xmax,m_nbins,m_listname);
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_OBSERVABLE_GETTER(Three_Particle_3Mass2,Three_Particle_3Mass2_Getter,"3Mass2")

Three_Particle_3Mass2::Three_Particle_3Mass2(const Flavour & flav1,const Flavour & flav2,
					 const Flavour & flav3,int type,double xmin,
					 double xmax,int nbins,const std::string & listname) :
  Three_Particle_Observable_Base(flav1,flav2,flav3,type,xmin,xmax,nbins,listname,"3Mass2") 
{ 
}


void Three_Particle_3Mass2::Evaluate(const Vec4D & mom1,const Vec4D & mom2,const Vec4D & mom3,double weight, double ncount) 
{ 
  
  double mass = (mom1+mom2+mom3).Abs2();
  p_histo->Insert(mass,weight,ncount); 
} 

Primitive_Observable_Base * Three_Particle_3Mass2::Copy() const 
{
    return new Three_Particle_3Mass2(m_flav1,m_flav2,m_flav3,m_type,m_xmin,m_xmax,m_nbins,m_listname);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_OBSERVABLE_GETTER(Three_Particle_3Mass,Three_Particle_3Mass_Getter,"3Mass")

Three_Particle_3Mass::Three_Particle_3Mass(const Flavour & flav1,const Flavour & flav2,
					 const Flavour & flav3,int type,double xmin,
					 double xmax,int nbins,const std::string & listname) :
  Three_Particle_Observable_Base(flav1,flav2,flav3,type,xmin,xmax,nbins,listname,"3Mass") 
{ 
}


void Three_Particle_3Mass::Evaluate(const Vec4D & mom1,const Vec4D & mom2,const Vec4D & mom3,double weight, double ncount) 
{ 
  
  double mass = sqrt( (mom1+mom2+mom3).Abs2() );
  p_histo->Insert(mass,weight,ncount); 
} 

Primitive_Observable_Base * Three_Particle_3Mass::Copy() const 
{
    return new Three_Particle_3Mass(m_flav1,m_flav2,m_flav3,m_type,m_xmin,m_xmax,m_nbins,m_listname);
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_OBSERVABLE_GETTER(Three_Particle_3EnergyCMS,Three_Particle_3EnergyCMS_Getter,"3EnergyCMS")

Three_Particle_3EnergyCMS::Three_Particle_3EnergyCMS(const Flavour & flav1,const Flavour & flav2,
                                         const Flavour & flav3,int type,double xmin,
                                         double xmax,int nbins,const std::string & listname) :
  Three_Particle_Observable_Base(flav1,flav2,flav3,type,xmin,xmax,nbins,listname,"3EnergyCMS")
{ 
}


void Three_Particle_3EnergyCMS::Evaluate(const Vec4D & mom1,const Vec4D & mom2,const Vec4D & mom3,double weight, double ncount)
{ 
  Vec4D sum = mom1+mom2+mom3;
  Poincare boost(sum);
  Vec4D p1 = boost*mom1;
  Vec4D p2 = boost*mom2;
  Vec4D p3 = boost*mom3;
  double E = p1[0];
  p_histo->Insert(2.0*E/rpa.gen.Ecms(),weight,ncount);
} 

Primitive_Observable_Base * Three_Particle_3EnergyCMS::Copy() const
{
    return new Three_Particle_3EnergyCMS(m_flav1,m_flav2,m_flav3,m_type,m_xmin,m_xmax,m_nbins,m_listname);
}

