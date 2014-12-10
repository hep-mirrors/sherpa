#include "METOOLS/Explicit/Current.H"
#include "METOOLS/Currents/C_Spinor.H"

namespace METOOLS {

  template <typename SType>
  class CF: public Current,
	    public Current_Contractor<SType> {
  public:

    typedef std::complex<SType>   SComplex;
    typedef std::vector<SComplex> SComplex_Vector;

    typedef CSpinor<SType> CSpinorType;
    typedef std::vector<CSpinorType*> CSpinorType_Vector;

  protected:

    SComplex m_cmass2, m_cmass;

    std::string CLabel() const;

  public:

    CF(const Current_Key &key);

    void ConstructJ(const ATOOLS::Vec4D &p,const int ch,
		    const int cr,const int ca,const int mode);
    void SetGauge(const ATOOLS::Vec4D &k);

    void AddPropagator();

    void SContract
    (const Current &c,const Int_Vector &pols,
     SComplex_Vector &ress,const size_t &offset) const;

    std::string Format(const CObject *c) const;

    char Type() const;    

  };// end of class CF

}// end of namespace METOOLS

#include "METOOLS/Explicit/Vertex.H"
#include "METOOLS/Explicit/Color_Calculator.H"
#include "METOOLS/Explicit/Dipole_Kinematics.H"
#include "METOOLS/Explicit/Dipole_Color.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/STL_Tools.H"
#include "ATOOLS/Org/MyStrStream.H"

#define M_I SComplex(0.0,1.0)

using namespace METOOLS;
using namespace ATOOLS;

template <typename SType>
CF<SType>::CF(const Current_Key &key): 
  Current(key),
  m_cmass2(Complex(sqr(this->m_mass),-this->m_mass*this->m_width)), 
  m_cmass(sqrt(m_cmass2))
{
  if (key.m_n==1 && key.p_model->ScalarNumber("WidthScheme")!=1) 
    m_cmass=sqrt(m_cmass2=Complex(sqr(this->m_mass),0.0));
}

template <typename SType>
void CF<SType>::ConstructJ(const ATOOLS::Vec4D &p,const int ch,
			   const int cr,const int ca,const int mode)
{
  this->m_p=p;
  this->ResetJ();
  bool anti(this->m_fl.IsAnti());
  if (this->m_fl.Majorana()) anti=(mode&1)?this->m_dir<0:this->m_dir>0;
  if (ch>=0) {
    CSpinorType j(anti^(this->m_dir>0)?
		  CSpinorType(this->m_fl.Majorana()?-2:-1,-this->m_dir,
			      this->m_fl.Majorana()?(mode?1:-1):1,
			      p,cr,ca,0,0,sqr(this->m_mass),
			      this->m_fl.MassSign()):
		  CSpinorType(this->m_fl.Majorana()?2:1,this->m_dir,
			      this->m_fl.Majorana()?(mode?-1:1):1,
			      p,cr,ca,0,0,sqr(this->m_mass),
			      this->m_fl.MassSign()));
    j.SetH(anti^(this->m_dir>0)?1:0);
#ifdef DEBUG__BG
    msg_Debugging()<<METHOD<<"(): "<<(this->m_dir>0?'I':'O')<<"+ "<<this->m_id
		   <<" "<<j<<" "<<(this->m_dir>0?this->m_fl.Bar():this->m_fl)
		   <<", m = "<<m_cmass<<" ("<<p.Mass()<<")\n";
#endif
    CSpinorType *c(CSpinorType::New(j));
    AddJ(c);
    if (p_sub) static_cast<Dipole_Color*>
      (p_sub->In().front()->Color().front())->AddJJK(c);
  }
  if (ch<=0) {
    CSpinorType j(anti^(this->m_dir>0)?
		  CSpinorType(this->m_fl.Majorana()?-2:-1,-this->m_dir,
			      this->m_fl.Majorana()?(mode?-1:1):-1,
			      p,cr,ca,0,0,sqr(this->m_mass),
			      this->m_fl.MassSign()):
		  CSpinorType(this->m_fl.Majorana()?2:1,this->m_dir,
			      this->m_fl.Majorana()?(mode?1:-1):-1,
			      p,cr,ca,0,0,sqr(this->m_mass),
		              this->m_fl.MassSign()));
    j.SetH(anti^(this->m_dir>0)?0:1);
#ifdef DEBUG__BG
    msg_Debugging()<<METHOD<<"(): "<<(this->m_dir>0?'I':'O')<<"- "<<this->m_id
		   <<" "<<j<<" "<<(this->m_dir>0?this->m_fl.Bar():this->m_fl)
		   <<", m = "<<m_cmass<<" ("<<p.Mass()<<")\n";
#endif
    CSpinorType *c(CSpinorType::New(j));
    AddJ(c);
    if (p_sub) static_cast<Dipole_Color*>
      (p_sub->In().front()->Color().front())->AddJJK(c);
  }
#ifdef DEBUG__BG
  if (p_sub) Print();
#endif
}

template <typename SType>
void CF<SType>::SetGauge(const ATOOLS::Vec4D &k)
{
}

template <typename SType>
void CF<SType>::AddPropagator()
{
  const CSpinorType hs;
  // add propagator for off-shell leg
  SComplex prop(M_I/(SType(this->m_p.Abs2())-m_cmass2));
  if (this->m_osd) prop=SComplex(M_I);
  SComplex pp(Spinor<SType>::PPlus(this->m_p));
  SComplex pm(Spinor<SType>::PMinus(this->m_p));
  SComplex pt(Spinor<SType>::PT(this->m_p));
  SComplex ptc(Spinor<SType>::PTC(this->m_p));
#ifdef DEBUG__BG
  msg_Debugging()<<"propagator: "<<prop
		 <<" <- p^2 = "<<this->m_p.Abs2()<<", m = "<<m_cmass<<"\n";
  msg_Debugging()<<"pp = "<<pp<<", pm = "<<pm<<", pt = "<<pt<<"\n";
#endif
  for (size_t i(0);i<m_j.size();++i) {
  CSpinorType_Vector *j(m_j[i].template Get<CSpinorType>());
  for (typename CSpinorType_Vector::iterator 
	 jit(j->begin());jit!=j->end();++jit) {
    CSpinorType j((*jit)->R(),(*jit)->B(),(**jit)(0),(**jit)(1),
		  (*jit)->H(),(*jit)->S(),
		  ((*jit)->On()&1)<<1|((*jit)->On()&2)>>1);
    if ((*jit)->B()>0) {// S(-p)
      j[0]=-pm*(**jit)[2]+ptc*(**jit)[3];
      j[1]=pt*(**jit)[2]-pp*(**jit)[3];
      j[2]=-pp*(**jit)[0]-ptc*(**jit)[1];
      j[3]=-pt*(**jit)[0]-pm*(**jit)[1];
    }
    else {// S(p)
      j[0]=(**jit)[2]*pp+(**jit)[3]*pt;
      j[1]=(**jit)[2]*ptc+(**jit)[3]*pm;
      j[2]=(**jit)[0]*pm-(**jit)[1]*pt;
      j[3]=-(**jit)[0]*ptc+(**jit)[1]*pp;
    }
    if (m_fl.MassSign()>=0)
      **jit=(this->m_msv?j+**jit*m_cmass:j)*prop;
    else **jit=(this->m_msv?j-**jit*m_cmass:j)*prop;
  }
  }
}

template <typename SType> void CF<SType>::SContract
(const Current &c,const Int_Vector &pols,
 SComplex_Vector &ress,const size_t &offset) const
{
#ifdef DEBUG__BG
  msg_Debugging()<<METHOD<<"(): {\n";
  msg_Indent();
#endif
  double phase(0.0);
  const std::vector<size_t> *pm(NULL);
  if (p_sub) {
    Vertex *v(p_sub->Sub()->In().front());
    if (v->Info()->Mode()==1) {
      phase=v->Kin()->Phase(offset==1?0:1);
      pm=&v->Kin()->PM();
    }
  }
  if (c.Type()!='F') THROW(fatal_error,"Invalid current type.");
  size_t i(0);
  for (typename CObject_Matrix::const_iterator 
	 ajit1(m_j.begin());ajit1!=m_j.end();++ajit1) {	
    const CSpinorType_Vector *j(ajit1->Get<CSpinorType>());
    for (typename CObject_Matrix::const_iterator 
	   ajit2(c.J().begin());ajit2!=c.J().end();++ajit2,++i) {
      // if (!pols[i]) continue;
      const CSpinorType_Vector *cj(ajit2->Get<CSpinorType>());
      for (typename CSpinorType_Vector::const_iterator 
	     jit2(cj->begin());jit2!=cj->end();++jit2) 
	for (typename CSpinorType_Vector::const_iterator 
	       jit1(j->begin());jit1!=j->end();++jit1)
	  if ((**jit1)(0)==(**jit2)(1) && (**jit1)(1)==(**jit2)(0) &&
	      (*jit1)->S()==offset && (*jit2)->S()==offset) {
#ifdef DEBUG__BG
	    msg_Debugging()<<"Add ("<<m_hm[i]<<")"
		       <<**jit1***jit2<<" ["<<offset<<"]\n";
#endif
	    ress[m_hm[i]]+=**jit1***jit2;
	    if (offset && pm) {
#ifdef DEBUG__BG
	      msg_Debugging()<<"Add ("<<(*pm)[m_hm[i]]<<")"
			     <<**jit1***jit2*SType(phase)<<" ["
			     <<offset<<"] ( phase = "<<phase<<" )\n";
#endif
	      ress[(*pm)[m_hm[i]]]+=**jit1***jit2*SType(phase);
	    }
	  }
    }
  }
#ifdef DEBUG__BG
  msg_Debugging()<<"}\n";
#endif
}

template <typename SType>
std::string CF<SType>::Format(const CObject *c) const
{
  return ToString(*(CSpinorType*)c,6);
}

template <typename SType>
char CF<SType>::Type() const
{
  return 'F';
}

template <typename SType>
std::string CF<SType>::CLabel() const
{
  return "fermion,label.side=right,label=$"+
    (this->m_out.empty()?this->m_fl.Bar():this->m_fl).TexName()+"$";
}

DECLARE_GETTER(CF<double>,"DF",Current,Current_Key);

Current *ATOOLS::Getter<Current,Current_Key,CF<double> >::
operator()(const Current_Key &key) const
{
  if (key.m_fl.IsFermion()) return new CF<double>(key);
  return NULL;
}

void ATOOLS::Getter<Current,Current_Key,CF<double> >::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"fermion current (double)";
}
