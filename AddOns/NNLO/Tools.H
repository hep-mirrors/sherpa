#ifndef NNLO__src__Tools_H
#define NNLO__src__Tools_H

#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "ATOOLS/Phys/NLO_Subevt.H"
#include "PDF/Main/PDF_Base.H"
#include "PHASIC++/Process/Process_Base.H"
#include "AMEGIC++/Amplitude/Zfunctions/Basic_Sfuncs.H"

namespace SHNNLO {

  extern PDF::PDF_Base *s_pdf;
  extern MODEL::Running_AlphaS *s_as;
  extern AMEGIC::Basic_Sfuncs *s_bs;

  extern double pmap[6],s_pdfmin[2];
  Complex spa(int i,int j);
  Complex spb(int i,int j);
  double sij(int i,int j);

  double Li2(double z);
  double Li3(double z);
  double Li4(double z);
  double S2p2(double z);
  double H0p0m(double z);
  double H00mp(double z);
  double H0mpp(double z);
  double H0mpm(double z);
  double H0mmp(double z);
  double Hmpmp(double z);

  double GetX(const ATOOLS::Vec4D &p,const int i);
  double PDF(int i,double x,double muf);
  double GetPDF(const ATOOLS::Flavour &fl,double x,double muf2);
  double GetXPDF(ATOOLS::Cluster_Leg *l,double muf2);
  double Beta0(const double &nf);
  double Hab(const ATOOLS::Flavour &a,const ATOOLS::Flavour &b);
  double FPab(const ATOOLS::Flavour &a,const ATOOLS::Flavour &b,const double &z);
  double SPab(const ATOOLS::Flavour &a,const ATOOLS::Flavour &b,const double &z);
  double IPab(const ATOOLS::Flavour &a,const ATOOLS::Flavour &b,const double &x);

  ATOOLS::Cluster_Leg *GetSplitter
  (const ATOOLS::Cluster_Amplitude &ampl);

  double CollinearCounterTerms
  (ATOOLS::Cluster_Leg *const l,
   const double &t1,const double &t2,
   const double &lmur2,const double &lmuf2,
   const double &ran);

  double Weight(double &ct,ATOOLS::Cluster_Amplitude *const ampl,
		const double &mur2,const double &muf2,
		const double *k0sq,const int mode);

  double NLODiffWeight(PHASIC::Process_Base *const proc,double &wgt,
		       const double &mur2,const double &muf2,
		       const double *k0sq,const int fomode,
		       const int umode,const std::string &varid);

  double NNLODeltaWeight(PHASIC::Process_Base *const proc,
			 const double &wgt,const int fomode);
  double NNLODiffWeight(PHASIC::Process_Base *const proc,double &wgt,
			const double &mur2,const double &muf2,
			const double *k0sq,const int mode,const int fomode,
			const int umode,const std::string &varid);

  extern long int s_ntrials;
  extern double   s_disc, s_p1, s_p2, s_z[4];

  double SetWeight(ATOOLS::Cluster_Amplitude *const ampl,
		   const int mode,double wgt,
		   const double &w,const double &w1,
		   const int umode,const std::string &varid);

}

#endif
