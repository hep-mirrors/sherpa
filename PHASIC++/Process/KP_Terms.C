#include "PHASIC++/Process/KP_Terms.H"

#include "PHASIC++/Main/Process_Integrator.H"
#include "PDF/Main/ISR_Handler.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace PHASIC;
using namespace ATOOLS;

KP_Terms::KP_Terms(Process_Base *const proc,const int mode):
  m_imode(6),
  p_proc(proc), p_flkern(new Flavour_Kernels()), p_masskern(NULL)
{
  SetNC(3.0);
  m_flavs=p_proc->Flavours();
  m_massive=mode&1;
  m_cemode=mode&2;
  size_t fsgluons(0);
  for (size_t i=p_proc->NIn();i<m_flavs.size();i++) {
    if (m_flavs[i].Strong()&&m_flavs[i].IsMassive()) m_massive=1;
    if (m_flavs[i].IsGluon()) fsgluons++; 
  }
  if (m_massive||fsgluons) {
    p_masskern = new Massive_Kernels();
    if (p_masskern->Nmf()>0) m_massive=1;
    if (!m_massive) {
      delete p_masskern;
      p_masskern=NULL;
    }
    else {
      m_xpa.resize(p_masskern->Nmf()*fsgluons);
      m_xpb.resize(p_masskern->Nmf()*fsgluons);
    }
  }
  for (size_t i(0);i<p_proc->Flavours().size();++i)
    if (p_proc->Flavours()[i].Strong()) m_plist.push_back(i);
  for (int i=0;i<8;i++) m_kpca[i]=0.;
  for (int i=0;i<8;i++) m_kpcb[i]=0.;
}

KP_Terms::~KP_Terms()
{
  if (p_masskern) delete p_masskern;
  delete p_flkern;
}

void KP_Terms::SetNC(const double &nc)
{
  p_flkern->SetNC(nc);
  if (p_masskern) p_masskern->SetNC(nc);
}

void KP_Terms::SetCoupling(const MODEL::Coupling_Map *cpls)
{
  if (cpls->find("Alpha_QCD")!=cpls->end()) p_cpl=cpls->find("Alpha_QCD")->second;
  else THROW(fatal_error,"Coupling not found");
  msg_Tracking()<<"DipoleSplitting_Base:: alpha = "<<*p_cpl<<std::endl;
  m_cpldef = p_cpl->Default()/(2.*M_PI);
}

void KP_Terms::SetAlpha(const double &aff,const double &afi,
			const double &aif,const double &aii)
{
  if (!p_masskern) p_flkern->SetAlpha(aff);
  else             p_flkern->SetAlpha(1.);
  if (p_masskern)  p_masskern->SetAlpha(aff,afi,aif,aii);
}

void KP_Terms::Calculate
(const Vec4D_Vector &mom,const double &x0,const double &x1,
 const double &eta0,const double &eta1,const double &weight)
{
  bool sa=m_flavs[0].Strong();
  bool sb=m_flavs[1].Strong();
  if (!sa && !sb) return;
  if (x0<eta0 || x1<eta1) return; 
  size_t pls=1;
  if (sa&&sb) pls++;
  double muf = p_proc->ScaleSetter()->Scale(stp::fac,1);
  for (int i=0;i<8;i++) m_kpca[i]=0.;
  for (int i=0;i<8;i++) m_kpcb[i]=0.;

  if (sa) {
    double w=1.-eta0;
    int type=m_flavs[0].IntSpin();
    if (m_imode&2) {
    m_kpca[0]=-w*p_flkern->Kb1(type,x0)+p_flkern->Kb2(type)-p_flkern->Kb4(type,eta0);
    m_kpca[1]=w*(p_flkern->Kb1(type,x0)+p_flkern->Kb3(type,x0));
    m_kpca[2]=-w*p_flkern->Kb1(type+2,x0)+p_flkern->Kb2(type+2)-p_flkern->Kb4(type+2,eta0);
    m_kpca[3]=w*(p_flkern->Kb1(type+2,x0)+p_flkern->Kb3(type+2,x0));
    for (int i=0;i<4;i++) m_kpca[i]*=m_dsij[0][0];

    double t=0.;
    if (!m_massive) {
      for (size_t i=pls;i<m_plist.size();i++) {
	int itype=m_flavs[m_plist[i]].IntSpin();
	t+= m_dsij[0][i]*p_flkern->ft(itype);
      }
      m_kpca[type*2-2]+=t*(-w*p_flkern->t1(x0)+p_flkern->t2()-p_flkern->t4(eta0));
      m_kpca[type*2-1]+=t*w*p_flkern->t1(x0);
    }
    else {
      size_t xpcnt(0);
      for (size_t i=pls;i<m_plist.size();i++) {
	int spin=m_flavs[m_plist[i]].IntSpin();
	double saj=dabs(2.*mom[m_plist[0]]*mom[m_plist[i]]);
	double muq2=saj;
	double sajx = saj/x0;
	double muq2x=sajx;
	if (spin!=2) muq2=sqr(m_flavs[m_plist[i]].Mass())/saj;
	if (spin!=2) muq2x=sqr(m_flavs[m_plist[i]].Mass())/sajx;
	m_kpca[0]+=m_dsij[0][i]*
	  (-w*p_masskern->t1(type,spin,muq2,x0)-w*p_masskern->t1p(type,spin,muq2,x0)
	  +p_masskern->t2(type,spin,muq2)-p_masskern->t4(type,spin,muq2,eta0));
	m_kpca[1]+=m_dsij[0][i]*
	  (w*(p_masskern->t1(type,spin,muq2x,x0)+p_masskern->t1p(type,spin,muq2,x0)+p_masskern->t3(type,spin,muq2x,x0)));
	m_kpca[2]+=m_dsij[0][i]*
	  (-w*p_masskern->t1(type+2,spin,muq2,x0)-w*p_masskern->t1p(type+2,spin,muq2,x0)
	  +p_masskern->t2(type+2,spin,muq2)-p_masskern->t4(type+2,spin,muq2,eta0));
	m_kpca[3]+=m_dsij[0][i]*
	  (w*(p_masskern->t1(type+2,spin,muq2x,x0)+p_masskern->t1p(type+2,spin,muq2,x0)+p_masskern->t3(type+2,spin,muq2x,x0)));	
	if (spin==2) {
          for (size_t j=0;j<p_masskern->Nmf();j++) {
	    m_xpa[xpcnt].xp=1.-4.*sqr(p_masskern->FMass(j))/saj;
	    if (m_xpa[xpcnt].xp>eta0) {
	      m_kpca[0]+=m_dsij[0][i]*p_masskern->t6(type,m_xpa[xpcnt].xp);
	      m_kpca[1]+=m_dsij[0][i]*w*p_masskern->t5(type,x0,m_xpa[xpcnt].xp);
	      m_kpca[2]+=m_dsij[0][i]*p_masskern->t6(type+2,m_xpa[xpcnt].xp);
	      m_kpca[3]+=m_dsij[0][i]*w*p_masskern->t5(type+2,x0,m_xpa[xpcnt].xp);
	      
	      m_xpa[xpcnt].kpc=m_dsij[0][i]*
		(-w*p_masskern->t5(type,x0,m_xpa[xpcnt].xp)-p_masskern->t6(type,m_xpa[xpcnt].xp)-p_masskern->t7(type,eta0,m_xpa[xpcnt].xp));
	    }
	    xpcnt++;
	  }
	}
      }
    }

    if (sb) {
      if(!m_massive){
        m_kpca[0]-=m_dsij[0][1]*(-w*p_flkern->Kt1(type,x0)+p_flkern->Kt2(type)-p_flkern->Kt4(type,eta0));
        m_kpca[1]-=m_dsij[0][1]*w*(p_flkern->Kt1(type,x0)+p_flkern->Kt3(type,x0));
        m_kpca[2]-=m_dsij[0][1]*(-w*p_flkern->Kt1(type+2,x0)+p_flkern->Kt2(type+2)-p_flkern->Kt4(type+2,eta0));
        m_kpca[3]-=m_dsij[0][1]*w*(p_flkern->Kt1(type+2,x0)+p_flkern->Kt3(type+2,x0));
      }
      else{
        m_kpca[0]-=m_dsij[0][1]*(-w*p_masskern->Kt1(type,x0)+p_masskern->Kt2(type)-p_masskern->Kt4(type,eta0));
        m_kpca[1]-=m_dsij[0][1]*w*(p_masskern->Kt1(type,x0)+p_masskern->Kt3(type,x0));
        m_kpca[2]-=m_dsij[0][1]*(-w*p_masskern->Kt1(type+2,x0)+p_masskern->Kt2(type+2)-p_masskern->Kt4(type+2,eta0));
        m_kpca[3]-=m_dsij[0][1]*w*(p_masskern->Kt1(type+2,x0)+p_masskern->Kt3(type+2,x0));
      }
    }
    }
    
    if (m_imode&4) {
    double asum=0.,fsum=0.;
    for (size_t i=1;i<m_plist.size();i++) {
      fsum+=m_dsij[0][i];
      asum+=m_dsij[0][i]*log(muf/dabs(2.*mom[m_plist[0]]*mom[m_plist[i]]));
    }
    if (fsum!=0.0) asum/=fsum;
    m_kpca[4]=fsum*(-w*p_flkern->P1(type,x0)+p_flkern->P2(type)-p_flkern->P4(type,eta0));
    m_kpca[5]=fsum*w*(p_flkern->P1(type,x0)+p_flkern->P3(type,x0));
    m_kpca[6]=fsum*(-w*p_flkern->P1(type+2,x0)+p_flkern->P2(type+2)-p_flkern->P4(type+2,eta0));
    m_kpca[7]=fsum*w*(p_flkern->P1(type+2,x0)+p_flkern->P3(type+2,x0));
    m_kpca[0]+=asum*m_kpca[4];
    m_kpca[1]+=asum*m_kpca[5];
    m_kpca[2]+=asum*m_kpca[6];
    m_kpca[3]+=asum*m_kpca[7];
    }
  }
  
  if (sb) {
    double w=1.-eta1;
    int type=m_flavs[1].IntSpin();
    if (m_imode&2) {
    m_kpcb[0]=-w*p_flkern->Kb1(type,x1)+p_flkern->Kb2(type)-p_flkern->Kb4(type,eta1);
    m_kpcb[1]=w*(p_flkern->Kb1(type,x1)+p_flkern->Kb3(type,x1));
    m_kpcb[2]=-w*p_flkern->Kb1(type+2,x1)+p_flkern->Kb2(type+2)-p_flkern->Kb4(type+2,eta1);
    m_kpcb[3]=w*(p_flkern->Kb1(type+2,x1)+p_flkern->Kb3(type+2,x1));
    for (int i=0;i<4;i++) m_kpcb[i]*=m_dsij[0][0];

    double t=0.;
    if (!m_massive) {
      for (size_t i=pls;i<m_plist.size();i++) {
	int itype=m_flavs[m_plist[i]].IntSpin();
	t+= m_dsij[pls-1][i]*p_flkern->ft(itype);
      }
      m_kpcb[type*2-2]+=t*(-w*p_flkern->t1(x1)+p_flkern->t2()-p_flkern->t4(eta1));
      m_kpcb[type*2-1]+=t*w*p_flkern->t1(x1);
    }
    else {
      size_t xpcnt(0);
      for (size_t i=pls;i<m_plist.size();i++) {
	int spin=m_flavs[m_plist[i]].IntSpin();
	double saj=dabs(2.*mom[m_plist[pls-1]]*mom[m_plist[i]]);
	double muq2=saj;
	double sajx = saj/x1;
	double muq2x=sajx;
	if (spin!=2) muq2=sqr(m_flavs[m_plist[i]].Mass())/saj;// mu-tilde
	if (spin!=2) muq2x=sqr(m_flavs[m_plist[i]].Mass())/sajx;// mu
	m_kpcb[0]+=m_dsij[pls-1][i]*
	  (-w*p_masskern->t1(type,spin,muq2,x1) - w*p_masskern->t1p(type,spin,muq2,x1) 
	  +p_masskern->t2(type,spin,muq2)-p_masskern->t4(type,spin,muq2,eta1));
	m_kpcb[1]+=m_dsij[pls-1][i]*
	  (w*(p_masskern->t1(type,spin,muq2x,x1) + p_masskern->t1p(type,spin,muq2,x1) +p_masskern->t3(type,spin,muq2x,x1)));
	m_kpcb[2]+=m_dsij[pls-1][i]*
	  (-w*p_masskern->t1(type+2,spin,muq2,x1)-w*p_masskern->t1p(type+2,spin,muq2,x1)
	  +p_masskern->t2(type+2,spin,muq2)-p_masskern->t4(type+2,spin,muq2,eta1));
	m_kpcb[3]+=m_dsij[pls-1][i]*
	  (w*(p_masskern->t1(type+2,spin,muq2x,x1)+p_masskern->t1p(type+2,spin,muq2,x1)+p_masskern->t3(type+2,spin,muq2x,x1)));
	if (spin==2) {
          for (size_t j=0;j<p_masskern->Nmf();j++) {
	    m_xpb[xpcnt].xp=1.-4.*sqr(p_masskern->FMass(j))/saj;
	    if (m_xpb[xpcnt].xp>eta1) {
	      m_kpcb[0]+=m_dsij[pls-1][i]*p_masskern->t6(type,m_xpb[xpcnt].xp);
	      m_kpcb[1]+=m_dsij[pls-1][i]*w*p_masskern->t5(type,x1,m_xpb[xpcnt].xp);
	      m_kpcb[2]+=m_dsij[pls-1][i]*p_masskern->t6(type+2,m_xpb[xpcnt].xp);
	      m_kpcb[3]+=m_dsij[pls-1][i]*w*p_masskern->t5(type+2,x1,m_xpb[xpcnt].xp);

	      m_xpb[xpcnt].kpc=m_dsij[pls-1][i]*
		(-w*p_masskern->t5(type,x1,m_xpb[xpcnt].xp)-p_masskern->t6(type,m_xpb[xpcnt].xp)-p_masskern->t7(type,eta1,m_xpb[xpcnt].xp));
	    }
	    xpcnt++;
	  }
	}
      }
    }

    if (sa) {
      if (!m_massive){
        m_kpcb[0]-=m_dsij[0][1]*(-w*p_flkern->Kt1(type,x1)+p_flkern->Kt2(type)-p_flkern->Kt4(type,eta1));
        m_kpcb[1]-=m_dsij[0][1]*w*(p_flkern->Kt1(type,x1)+p_flkern->Kt3(type,x1));
        m_kpcb[2]-=m_dsij[0][1]*(-w*p_flkern->Kt1(type+2,x1)+p_flkern->Kt2(type+2)-p_flkern->Kt4(type+2,eta1));
        m_kpcb[3]-=m_dsij[0][1]*w*(p_flkern->Kt1(type+2,x1)+p_flkern->Kt3(type+2,x1));
      }
      else {
        m_kpcb[0]-=m_dsij[0][1]*(-w*p_masskern->Kt1(type,x1)+p_masskern->Kt2(type)-p_masskern->Kt4(type,eta1));
        m_kpcb[1]-=m_dsij[0][1]*w*(p_masskern->Kt1(type,x1)+p_masskern->Kt3(type,x1));
        m_kpcb[2]-=m_dsij[0][1]*(-w*p_masskern->Kt1(type+2,x1)+p_masskern->Kt2(type+2)-p_masskern->Kt4(type+2,eta1));
        m_kpcb[3]-=m_dsij[0][1]*w*(p_masskern->Kt1(type+2,x1)+p_masskern->Kt3(type+2,x1));
      }
    }
    }

    if (m_imode&4) {
    double asum=0.,fsum=0.;
    for (size_t i=0;i<m_plist.size();i++) if (i!=pls-1) {
      fsum+=m_dsij[pls-1][i];
      asum+=m_dsij[pls-1][i]*log(muf/dabs(2.*mom[m_plist[pls-1]]*mom[m_plist[i]]));
    }
    if (fsum!=0.0) asum/=fsum;
    m_kpcb[4]=fsum*(-w*p_flkern->P1(type,x1)+p_flkern->P2(type)-p_flkern->P4(type,eta1));
    m_kpcb[5]=fsum*w*(p_flkern->P1(type,x1)+p_flkern->P3(type,x1));
    m_kpcb[6]=fsum*(-w*p_flkern->P1(type+2,x1)+p_flkern->P2(type+2)-p_flkern->P4(type+2,eta1));
    m_kpcb[7]=fsum*w*(p_flkern->P1(type+2,x1)+p_flkern->P3(type+2,x1));
    m_kpcb[0]+=asum*m_kpcb[4];
    m_kpcb[1]+=asum*m_kpcb[5];
    m_kpcb[2]+=asum*m_kpcb[6];
    m_kpcb[3]+=asum*m_kpcb[7];
    }
  }

  double gfac=weight;
  if (sa) gfac/=(1.-eta0);
  if (sb) gfac/=(1.-eta1);

  if (sa) {
    for (int i=0;i<8;i++) m_kpca[i]*=gfac;
    for (size_t i=0;i<m_xpa.size();i++) m_xpa[i].kpc*=gfac;
  }
  if (sb) {
    for (int i=0;i<8;i++) m_kpcb[i]*=gfac;
    for (size_t i=0;i<m_xpb.size();i++) m_xpb[i].kpc*=gfac;
  }
}

double KP_Terms::Get(const double &x0,const double &x1,
		     const double &eta0,const double &eta1,
		     const ATOOLS::Flavour_Vector &flav,
		     const int mode)
{
  bool sa=flav[0].Strong();
  bool sb=flav[1].Strong();
  PDF::PDF_Base *pdfa(p_proc->Integrator()->ISR()->PDF(mode));
  PDF::PDF_Base *pdfb(p_proc->Integrator()->ISR()->PDF(1-mode));
  if (sa && (pdfa==NULL || !pdfa->Contains(flav[0]))) return 0.0;
  if (sb && (pdfb==NULL || !pdfb->Contains(flav[1]))) return 0.0;
  if (!sa && !sb) return 0.;
  if (x0<eta0 || x1<eta1) return 0.; 
  size_t pls=1;
  if (sa&&sb) pls++;
  Flavour gluon(kf_gluon);
  Flavour quark(kf_quark);
  double fa=0.,faq=0.,fag=0.,faqx=0.,fagx=0.;
  double fb=0.,fbq=0.,fbg=0.,fbqx=0.,fbgx=0.;
  double g2massq(0.);
  double muf = p_proc->ScaleSetter()->Scale(stp::fac);

  if (sa) {
    if (m_cemode && eta0*rpa->gen.PBeam(0)[0]<flav[0].Mass(true)) {
      msg_Tracking()<<METHOD<<"(): E < m ! ( "<<eta0*rpa->gen.PBeam(0)[0]
		    <<" vs. "<<flav[0].Mass(true)<<" )"<<std::endl;
      return 0.0;
    }
    pdfa->Calculate(eta0/x0,muf);
    fagx = pdfa->GetXPDF(gluon)/eta0;
    if (flav[0].IsQuark()) faqx = pdfa->GetXPDF(flav[0])/eta0;
    else {
      for (size_t i=0;i<quark.Size();i++) faqx+= pdfa->GetXPDF(quark[i]);
      faqx/=eta0;
    }
    pdfa->Calculate(eta0,muf);
    fa  = pdfa->GetXPDF(flav[0])/eta0;
    if ((m_cemode && IsZero(fa,1.0e-16)) || !(fa>0.)) return 0.;
    fag = pdfa->GetXPDF(gluon)/eta0;
    if (flav[0].IsQuark()) faq = fa;
    else {
      for (size_t i=0;i<quark.Size();i++) faq+= pdfa->GetXPDF(quark[i]);
      faq/=eta0;
    }

    if (m_massive) {
      for (size_t i=0;i<m_xpa.size();i++) if (m_xpa[i].xp>eta0) {
	pdfa->Calculate(eta0/m_xpa[i].xp,muf);
	g2massq+=m_xpa[i].kpc*pdfa->GetXPDF(flav[0])/eta0/fa;
      }
    }    
  }
  if (sb) {
    if (m_cemode && eta1*rpa->gen.PBeam(1)[0]<flav[1].Mass(true)) {
      msg_Tracking()<<METHOD<<"(): E < m ! ( "<<eta1*rpa->gen.PBeam(1)[0]
		    <<" vs. "<<flav[1].Mass(true)<<" )"<<std::endl;
      return 0.0;
    }
    pdfb->Calculate(eta1/x1,muf);
    fbgx = pdfb->GetXPDF(gluon)/eta1;
    if (flav[1].IsQuark()) fbqx = pdfb->GetXPDF(flav[1])/eta1;
    else {
      for (size_t i=0;i<quark.Size();i++) fbqx+= pdfb->GetXPDF(quark[i]);
      fbqx/=eta1;
    }
    pdfb->Calculate(eta1,muf);
    fb = pdfb->GetXPDF(flav[1])/eta1;
    if ((m_cemode && IsZero(fb,1.0e-16)) || !(fb>0.)) return 0.;
    fbg = pdfb->GetXPDF(gluon)/eta1;
    if (flav[1].IsQuark()) fbq = fb;
    else {
      for (size_t i=0;i<quark.Size();i++) fbq+= pdfb->GetXPDF(quark[i]);
      fbq/=eta1;
    }

    if (m_massive) {
      for (size_t i=0;i<m_xpb.size();i++) if (m_xpb[i].xp>eta1) {
	pdfb->Calculate(eta1/m_xpb[i].xp,muf);
	g2massq+=m_xpb[i].kpc*pdfb->GetXPDF(flav[1])/eta1/fb;
      }
    }    
  }

  double res=g2massq;
  if (sa) {
    res+= (m_kpca[0]*faq+m_kpca[1]*faqx+m_kpca[2]*fag+m_kpca[3]*fagx)/fa;
  }
  
  if (sb) {
    res+= (m_kpcb[0]*fbq+m_kpcb[1]*fbqx+m_kpcb[2]*fbg+m_kpcb[3]*fbgx)/fb;
  }
  return res;
}

void KP_Terms::FillMEwgts(ATOOLS::ME_Weight_Info &wgtinfo)
{
  if (wgtinfo.m_type&mewgttype::KP) {
    for (int i=0;i<4;i++) {
      wgtinfo.m_wfac[i]=m_kpca[i];
      wgtinfo.m_wfac[i+4]=m_kpcb[i];
      wgtinfo.m_wfac[i+8]=m_kpca[i+4];
      wgtinfo.m_wfac[i+12]=m_kpcb[i+4];
    }
  }
}

void KP_Terms::SetDSij(const std::vector<std::vector<double> > &ds)
{
  m_dsij.resize(ds.size());
  for (size_t i(0);i<ds.size();++i) {
    m_dsij[i].resize(ds[i].size());
    for (size_t j(0);j<ds[i].size();++j) m_dsij[i][j]=ds[i][j];
  }
}
