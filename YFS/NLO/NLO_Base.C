#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Message.H"
#include "YFS/NLO/NLO_Base.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include "YFS/NLO/Virtual.H"
#include "YFS/NLO/VirtualVirtual.H"

using namespace YFS;
using namespace MODEL;
using namespace ATOOLS;
using namespace std;

double massmin=2220;
double rcount=1;
double sumw =0;

std::ofstream out_recola;
std::ofstream out_sub, out_real, out_finite;

static double Lambda(double x,double y,double z)
{
  return abs(x*x+y*y+z*z-2.*x*y-2.*x*z-2.*y*z);
}

NLO_Base::NLO_Base() {
  p_yfsFormFact = new YFS::YFS_Form_Factor();
  p_nlodipoles = new YFS::Define_Dipoles();
  p_real = NULL;
  p_virt = NULL;
  p_realvirt = NULL;
  p_vv = NULL;
  m_evts = 0;
  m_recola_evts = 0;
  m_realtool = 0;
  m_realvirt = 0;
  m_looptool = 0;
  m_rrtool = 0;
  m_vvtool = 0;
  if(m_isr_debug || m_fsr_debug){
  	m_histograms2d["IFI_EIKONAL"] = new Histogram_2D(0, -1., 1., 20, 0, 5., 20 );
  	m_histograms2d["REAL_SUB"] = new Histogram_2D(0,  0, sqrt(m_s), 200, 0, sqrt(m_s)/2., 20);
  	m_histograms2d["REAL_COLL_RATIO"] = new Histogram_2D(0,  0, 2.*M_PI, 20, 0, sqrt(m_s)/2., 125);
  	m_histograms2d["REAL_COLL_RATIO"] = new Histogram_2D(0, 0, sqrt(m_s)/2., 125, 0, 10, 20);
  	m_histograms2d["REAL_RATIO"] = new Histogram_2D(0,  0, 15, 16, 0, sqrt(m_s)/2., 200);
  	m_histograms1d["Real_diff"] = new Histogram(0, -1, 1, 100);
  	m_histograms2d["Real_Flux"] = new Histogram_2D(0, 0, 1.1, 50, 80, 100, 100);
  	m_histograms1d["k_E"] = new Histogram(0, 0, sqrt(m_s)/2, sqrt(m_s)/2);
  	m_histograms1d["k_pt"] = new Histogram(0, 0, sqrt(m_s)/2, sqrt(m_s)/2);
  	m_histograms1d["dip_mass"] = new Histogram(0, 0, sqrt(m_s), sqrt(m_s));
  	m_histograms1d["k_E_pass"] = new Histogram(0, 0, sqrt(m_s)/2, sqrt(m_s)/2);
  	m_histograms1d["k_pt_pass"] = new Histogram(0, 0, sqrt(m_s)/2, sqrt(m_s)/2);
  	m_histograms1d["dip_mass_pass"] = new Histogram(0, 0, sqrt(m_s), sqrt(m_s));
  	if (!ATOOLS::DirectoryExists(m_debugDIR_NLO)) ATOOLS::MakeDir(m_debugDIR_NLO);
  }
  if(m_check_poles==1){
  	if (!ATOOLS::DirectoryExists(m_debugDIR_NLO)) ATOOLS::MakeDir(m_debugDIR_NLO);
  	m_histograms1d["SinglePoleCD"] = new Histogram(0, 0, 25, 25);
  	m_histograms1d["SinglePoleVV"] = new Histogram(0, 0, 25, 25);
  	m_histograms1d["DoublePoleVV"] = new Histogram(0, 0, 25, 25);
  	m_histograms1d["OneLoopEpsLP"] = new Histogram(0, -1.5, -0.5, 50);
  	m_histograms1d["OneLoopEpsYFS"] = new Histogram(0, -1.5, -0.5, 50);
  	m_histograms1d["RealLoopEpsLP"] = new Histogram(0, -5,  0.0, 50);
  	m_histograms1d["RealLoopEpsYFS"] = new Histogram(0, -5, 0.0, 50);
  	m_histograms1d["relativediff"] = new Histogram(0, -20., -5.0, 50);
  	m_histograms1d["RVSinglePoleCD"] = new Histogram(0, 0, 25, 25);
  	m_histograms2d["REAL_SUB"] = new Histogram_2D(0,  0, sqrt(m_s)/2., 200, 0, 2*M_PI, 20);
  	m_histograms2d["REAL"] = new Histogram_2D(0,  0, sqrt(m_s)/2., 200, 0, 2*M_PI, 20);
  }
}

NLO_Base::~NLO_Base() {
  if(m_isr_debug || m_fsr_debug || m_check_real_sub ||m_check_poles ){
		Histogram_2D * histo2d;
		string name;
		for (map<string, Histogram_2D *>::iterator hit = m_histograms2d.begin();
		        hit != m_histograms2d.end(); hit++) {
			histo2d = hit->second;
			name  = string(m_debugDIR_NLO) + "/" + hit->first + string(".dat");
			// histo2d->MPISync();
			histo2d->Finalize();
			histo2d->Output(name);
			delete histo2d;
		}
		Histogram * histo1d;
		for (map<string, Histogram *>::iterator hit = m_histograms1d.begin();
		        hit != m_histograms1d.end(); hit++) {
			histo1d = hit->second;
			name  = string(m_debugDIR_NLO) +  "/" + hit->first + string(".dat");
			histo1d->MPISync();
			histo1d->Finalize();
			histo1d->Output(name);
			delete histo1d;
		}
	}
	// PRINT_VAR(massmin)
	if(p_yfsFormFact) delete p_yfsFormFact;
	if(p_nlodipoles) delete p_nlodipoles;
	if(p_realvirt) delete p_realvirt;
	if(p_real) delete p_real;
	if(p_virt) delete p_virt;
	if(p_vv) delete p_vv;
}


void NLO_Base::InitializeVirtual(const PHASIC::Process_Info& pi) {
	if(!m_eex_virt && !m_looptool){
		p_virt = new YFS::Virtual(pi);
		m_looptool = true;
	}
}

void NLO_Base::InitializeReal(const PHASIC::Process_Info& pi) {
	if(m_coll_real || m_realtool) return;
	p_real = new YFS::Real(pi);
	m_realtool = true;
}

void NLO_Base::InitializeRealVirtual(const PHASIC::Process_Info& pi) {
	// if(m_realvirt)
	p_realvirt = new YFS::RealVirtual(pi);
	m_realvirt = true;
	// m_looptool = true;
}

void NLO_Base::InitializeRealReal(const PHASIC::Process_Info& pi) {
	p_realreal = new YFS::RealReal(pi);
	m_rrtool = true;
}

void NLO_Base::InitializeVV(const PHASIC::Process_Info& pi) {
	p_vv = new YFS::VirtualVirtual(pi);
	m_vvtool = true;
}


void NLO_Base::Init(Flavour_Vector &flavs, Vec4D_Vector &plab, Vec4D_Vector &born) {
	m_flavs = flavs;
	m_plab = plab;
	m_bornMomenta = born;
}


double NLO_Base::CalculateNLO() {
	m_failcut = false;
	double result{0.0}, virt, real, rv, rr;
	// PRINT_VAR(m_looptool);
	if(!m_real_only) result += CalculateVirtual();
	virt = result;
	if(!m_virtual_only) result += CalculateReal();
	real = result-virt;
	result += CalculateRealReal();
	rr = result - real - virt - rv;
	result += CalculateRealVirtual();
	result += CalculateVV();
	rv = result - real - virt;
	if(m_failcut) return 0;
	return result;
}


double NLO_Base::CalculateVirtual() {
	if(m_eex_virt){
		// subtract born to avoid double counting
		// already present in eex!!
		return p_dipoles->CalculateEEXVirtual()*m_born-m_born;
	}
	if (!m_looptool) return 0;
	double virt;
	double sub;
	p_dipoles->p_yfsFormFact->p_virt = p_virt->p_loop_me;
	CheckMassReg();
	if(!HasISR()) virt = p_virt->Calc(m_bornMomenta, m_born);
	else virt = p_virt->Calc(m_plab, m_born);
	if(m_check_virt_born) {
			if (!IsEqual(m_born, p_virt->p_loop_me->ME_Born(), 1e-6)) {
			msg_Error() << METHOD << "\n Warning! Loop provider's born is different! YFS Subtraction likely fails\n"
									<< "Loop Provider " << ":  "<<p_virt->p_loop_me->ME_Born()
									<< "\nSherpa" << ":  "<<m_born<<std::endl
									<<"PhaseSpace Point = ";
			for(auto _p: m_plab) msg_Error()<<_p<<std::endl;
		}
	}	
	if(p_virt->FailCut()) return 0;
	if(m_virt_sub && p_virt->p_loop_me->Mode()!=1) sub = p_dipoles->CalculateVirtualSub();
	else sub = 0;
	m_oneloop = (virt- sub * m_born/m_rescale_alpha );
	if(p_virt->p_loop_me->Mode()==1) {
		m_oneloop /= m_rescale_alpha; 
		// PRINT_VAR(m_rescale_alpha);
	}
	if(IsBad(m_oneloop) || IsBad(sub)){
		msg_Error()<<"YFS Virtual is NaN"<<std::endl
							 <<"Virtual:  "<<virt<<std::endl
							 <<"Subtraction: "<<sub*m_born<<std::endl
							 <<"PhaseSpace Point: "<<std::endl<<m_plab<<std::endl;
	}
	if(m_check_poles==1){
		if(m_virt_sub==0) sub = p_dipoles->CalculateVirtualSub();
		double p1 = p_virt->p_loop_me->ME_E1()*p_virt->m_factor;
		double yfspole = p_dipoles->Get_E1();
		int ncorrect = ::countMatchingDigits(p1, -yfspole);
		double reldiff = (p1+yfspole)/p1;
		if(!IsEqual(p1,-yfspole,1e-4)){
			msg_Error()<<"Poles do not cancel in YFS Virtuals"<<std::endl
					 <<"Correct digits =  "<<ncorrect<<std::endl
					 <<"Relative diff =  "<<reldiff<<std::endl
					 // <<"Process =  "<<p_virt->p_loop_me->Name()<<std::endl
					 <<"One-Loop Provider V eps^{-1}  = "<<p1<<std::endl
					 <<"Sherpa V eps^{-1} = "<<yfspole<<std::endl
					 <<"Sherpa/One-Loop = "<<yfspole/p1<<std::endl;
					 return 0;
		}
		else{
			int i = 0;
			msg_Debugging()<<std::setprecision(32);
			msg_Debugging()<<"Poles cancel in YFS Virtuals to "<<ncorrect<<" digits"<<std::endl
					 		<<"Relative diff =  "<<reldiff<<std::endl;
			m_histograms1d["SinglePoleCD"]->Insert(ncorrect);
			m_histograms1d["OneLoopEpsYFS"]->Insert(log10(fabs(yfspole)));
			m_histograms1d["OneLoopEpsLP"]->Insert(log10(fabs(p1)));
			m_histograms1d["relativediff"]->Insert(log10(fabs(reldiff)));
			// msg_Out()<<"PhaseSpace point: "<<std::endl;
			// for(auto &p: m_plab) {
			// 	msg_Out()<<"p["<<i<<"] = "<<p<<std::endl;
			// 	i++;
			// }
			msg_Debugging()<<std::setprecision(32)<<"One-Loop Provider V eps^{-1}  = "<<p1<<std::endl
			 			<<"Sherpa V eps^{-1}  = "<<yfspole<<std::endl;
		}
		// Check Rescaling
		// Vec4D_Vector pscale = m_plab;
		// double virtscale;
		// for (double scale = 1; scale < 100; scale+=5)
		// {
		// 	std::vector<double> masses;
		// 	for (int i = 0; i < m_plab.size(); ++i)
		// 	{
		// 		masses.push_back(m_plab[i].Mass() *scale);
		// 	}
		// 	RescaleMasses(pscale, masses);
		// 	PRINT_VAR(masses);
		// 	PRINT_VAR(pscale);
		// 	PRINT_VAR(m_plab);
		// 	virtscale = p_virt->Calc(pscale, m_born);
		// 	msg_Out()<<"Virtual/ScaledVirtual = "<<virt/virtscale<<std::endl;
		// }
	}
	return m_oneloop;
}


double NLO_Base::CalculateReal() {
	if(m_coll_real) return p_dipoles->CalculateEEX()*m_born;
	if (!m_realtool) return 0;
	double real(0);
	double collreal = p_dipoles->CalculateEEX()*m_born;
	for (auto k : m_ISRPhotons) {
		if(m_check_real_sub && !HasFSR()) {
				if(k.E() < 0.2*sqrt(m_s)) continue;
				CheckRealSub(k,0);
		}
		if(m_isr_debug || m_fsr_debug){
			double fullreal = CalculateReal(k);
			real+=fullreal;
			double coll = p_dipoles->GetDipoleII()[0].Beta1(k);
			coll /=  p_dipoles->GetDipoleII()[0].Eikonal(k);
			if(fullreal!=0) m_histograms2d["REAL_COLL_RATIO"]->Insert(k.E(), coll*m_born/fullreal);
		}
		else {
			if(k.E() > m_hardmin) real+=CalculateReal(k,0);
		}
	}
	int fsrcount(0);
	for (auto k : m_FSRPhotons) {
		if(m_check_real_sub) {
			if(k.E() < 0.2*sqrt(m_s)) continue;
				CheckRealSub(k,1);
		}
		if(k.E() > m_hardmin){
			real+=CalculateReal(k, 1);
			fsrcount++;
		}
	}
	return real;
}


double NLO_Base::CalculateReal(Vec4D k, int fsrcount) {
	// if(fsrcount!=1) return 0;
	double norm = 2.*pow(2 * M_PI, 3);
	Vec4D_Vector p(m_plab),pi(m_bornMomenta), pf(m_bornMomenta);
  	dipoletype::code fluxtype;
	Vec4D kk = k;
	m_evts+=1;
	p_nlodipoles->MakeDipoles(m_flavs,m_plab,m_plab);
	fluxtype = p_nlodipoles->WhichResonant(k);
  // if(fluxtype==dipoletype::final || fsrcount==4){
  if(fsrcount==1 || fsrcount==4){
  	if(!HasFSR()) msg_Error()<<"Wrong dipole type in "<<METHOD<<endl;
  	for (Dipole_Vector::iterator Dip = p_nlodipoles->GetDipoleFF()->begin();
       Dip != p_nlodipoles->GetDipoleFF()->end(); ++Dip) {
  		 double scalek = p_fsr->ScalePhoton(k);
  		 Dip->SetPhotonScale(scalek);
  		 Dip->AddPhotonToDipole(k);
  		 if(!Dip->BoostNLO()) return 0;
  		 int i(0);
  		 for (auto f : Dip->m_flavs) {
      	p[p_nlodipoles->m_flav_label[f]] =  Dip->GetNewMomenta(i);
      	i++;
    	}
    	// k = Dip->m_dipolePhotons[0];
  	}
  }
 	else {
 		MapMomenta(p,k);
 	}
 	p.push_back(k);
 	// if(fluxtype==dipoletype::final || fsrcount==4) MapInitial(p);
 	if(fsrcount==1 || fsrcount==4) MapInitial(p);
 	CheckMasses(p,1);
 	Vec4D_Vector pp = p;
 	pp.pop_back();
	p_nlodipoles->MakeDipolesII(m_flavs,pp,m_plab);
	p_nlodipoles->MakeDipolesIF(m_flavs,pp,m_plab);
	p_nlodipoles->MakeDipoles(m_flavs,pp,m_plab);
	double r = p_real->Calc_R(p) / norm;
	if(p_real->FailCut()) m_failcut = true;
	double flux;
	if(m_flux_mode==1) flux = p_nlodipoles->CalculateFlux(k);
	else if(m_flux_mode==2) flux = 0.5*(p_nlodipoles->CalculateFlux(kk)+p_nlodipoles->CalculateFlux(k));
	else flux = p_dipoles->CalculateFlux(kk);
	double tot,rcoll;
	double subloc = p_nlodipoles->CalculateRealSub(k);
	double subb;
	m_real = r;
	if(fsrcount==0 || fsrcount==3) subb = p_dipoles->CalculateRealSubEEX(kk);
	else subb = p_dipoles->CalculateRealSubEEX(kk);
	// if(IsZero(subb)) return 0;
	m_eikeex = subb;
	m_subloc = subloc;
	if(!CheckMomentumConservation(p)) {
		msg_Error()<<"Momentum Conservation fails in "<<METHOD<<std::endl;
		if(m_isr_debug || m_fsr_debug){
			m_histograms1d["k_E"]->Insert(k.E());
			m_histograms1d["k_pt"]->Insert(k.PPerp());
			m_histograms1d["dip_mass"]->Insert((p[2]+p[3]).Mass());
		}
		return 0;
	}
	if((p[2]+p[3]).Mass() < massmin) massmin = (p[2]+p[3]).Mass();
	if(m_isr_debug || m_fsr_debug){
			m_histograms1d["k_E_pass"]->Insert(k.E());
			m_histograms1d["k_pt_pass"]->Insert(k.PPerp());
			m_histograms1d["dip_mass_pass"]->Insert((p[2]+p[3]).Mass());
		}
	if(IsZero(r)) return 0;
	if(IsBad(r) || IsBad(flux)) {
		msg_Error()<<"Bad point for YFS Real"<<std::endl
							 <<"Real ME is : "<<r<<std::endl
							 <<"Flux is : "<<flux<<std::endl;
		return 0;
	}
	// if(fsrcount==0 || fsrcount==3) flux=1;
	m_recola_evts+=1;
	if(m_submode==submode::local) tot =  (r*flux-subloc*m_born/m_rescale_alpha)/subloc;
	else if(m_submode==submode::global) tot =  (r*flux-subloc*m_born/m_rescale_alpha)/subb;
	else if(m_submode==submode::off) tot =  (r*flux)/subb;
	else msg_Error()<<METHOD<<" Unknown YFS Subtraction Mode "<<m_submode<<std::endl;
  if(m_isr_debug || m_fsr_debug){
		double diff = ((r/subloc - m_born)-( rcoll/subb - m_born))/((r/subloc - m_born)+( rcoll/subb - m_born));
		m_histograms1d["Real_diff"]->Insert(diff);
		if(m_isr_debug) m_histograms2d["Real_Flux"]->Insert(flux,sqrt(p_dipoles->GetDipoleII()->Sprime()));
  }
  if(m_no_subtraction) return r/subloc;
  if(IsBad(tot)){
  	msg_Error()<<"NLO real is NaN"<<std::endl
  							<<"R = "<<r<<std::endl
  							<<"Local  S = "<<subloc*m_born<<std::endl
  							<<"GLobal S = "<<subb<<std::endl;
  }
	if(m_isr_debug || m_fsr_debug) {
		m_histograms2d["IFI_EIKONAL"]->Insert(k.Y(),k.PPerp(), p_nlodipoles->CalculateRealSubIF(k));
		m_histograms2d["REAL_SUB"]->Insert((p[0]+p[1]).Mass(), k.E(), tot/m_born);
		m_histograms2d["REAL"]->Insert(k.E(), k.Theta(), r);
		m_histograms2d["REAL_SUB"]->Insert(k.E(), k.Theta(), tot);
	}
	sumw+=tot;
	rcount+=1;
	double avg=sumw/rcount;
	if(rcount==1000){
		m_ravg = avg;
	}
	if(rcount>1000){
		double diff = fabs(1.-m_ravg/avg)*100;
		if(diff>10){
			msg_Debugging()<<"Large jump in Real weight for "<<k<<std::endl;
			m_ravg = avg;
			// return 0;
		}
	}
	if(fsrcount>=3) return r*flux-subloc*m_born/m_rescale_alpha;
	return tot;
}

double NLO_Base::CalculateRealVirtual() {
	if (!m_realvirt) return 0;
	double realvirtual(0);
	for (auto k : m_ISRPhotons) {
		if(m_check_rv) {
			if(k.E() < 0.3*sqrt(m_s)) continue;
			CheckRealVirtualSub(k);
		}
		if(CheckPhotonForReal(k))	realvirtual+=CalculateRealVirtual(k,0);
		// realvirtual+=CalculateRealVirtual(k,0);
	}
	for (auto k : m_FSRPhotons) {
		// if(m_check_rv) {
		// 	if(k.E() < 0.2*sqrt(m_s)) continue;
		// 	CheckRealVirtualSub(k);
		// }
		if(CheckPhotonForReal(k)) realvirtual+=CalculateRealVirtual(k, 1);
		// realvirtual+=CalculateRealVirtual(k, 1);
	}
	// if(IsZero(realvirtual)) realvirtual = p_dipoles->CalculateRealSubEEX();
	return realvirtual;
}



double NLO_Base::CalculateRealVirtual(Vec4D k, int fsrcount) {
	if (!m_realvirt) return 0;
	Vec4D_Vector p(m_plab),pi(m_bornMomenta), pf(m_bornMomenta);
	double tot(0), sub(0);
	double norm = 2.*pow(2 * M_PI, 3);
	double flux(1);
	Vec4D kk = k;
	p_nlodipoles->MakeDipoles(m_flavs,m_plab,m_plab);
	dipoletype::code fluxtype = p_nlodipoles->WhichResonant(k);
  // if(fluxtype==dipoletype::final){
  if(fsrcount){
  	if(!HasFSR()) msg_Error()<<"Wrong dipole type in "<<METHOD<<endl;
  	for (Dipole_Vector::iterator Dip = p_nlodipoles->GetDipoleFF()->begin();
       Dip != p_nlodipoles->GetDipoleFF()->end(); ++Dip) {
  		 double scalek = p_fsr->ScalePhoton(k);
  		 Dip->SetPhotonScale(scalek);
  		 Dip->AddPhotonToDipole(k);
  		 if(!Dip->BoostNLO()) {
  		 	msg_Error()<<"NLO Boost Failed"<<std::endl;
  		 	return 0;
  		 }
  		 int i(0);
  		 for (auto f : Dip->m_flavs) {
      	p[p_nlodipoles->m_flav_label[f]] =  Dip->GetNewMomenta(i);
      	i++;
    		}
    	// k = Dip->m_dipolePhotons[0];
  	}
  }
 	else {
 		MapMomenta(p,k);
 	}
 	double yfspole;
 	p.push_back(k);
 	if(fsrcount) MapInitial(p);
 	CheckMasses(p,1);
	// CheckMasses(p, 1);
	Vec4D_Vector pp = p;
 	pp.pop_back();
 	Flavour_Vector fl = m_flavs;
 	p_nlodipoles->MakeDipolesII(fl,pp,pp);
	p_nlodipoles->MakeDipolesIF(fl,pp,pp);
	p_nlodipoles->MakeDipoles(fl,pp,pp);
	p_nlodipoles->p_yfsFormFact->p_virt = p_realvirt->p_loop_me;
	double subloc = p_nlodipoles->CalculateRealVirtualSubEps(k);
	yfspole = p_nlodipoles->Get_E1();
	double aB = p_nlodipoles->CalculateRealSub(k)*m_oneloop;//CalculateVirtual();//*p_realvirt->m_factor;
	

	p_nlodipoles->MakeDipolesII(fl,pp,m_plab);
	p_nlodipoles->MakeDipolesIF(fl,pp,m_plab);
	p_nlodipoles->MakeDipoles(fl,pp,m_plab);
	// p.push_back(k);
	// m_plab = pp;
	if(m_flux_mode==1) flux = p_nlodipoles->CalculateFlux(k);
	else if(m_flux_mode==2) flux = 0.5*(p_nlodipoles->CalculateFlux(kk)+p_nlodipoles->CalculateFlux(k));
	else flux = p_dipoles->CalculateFlux(kk);
	// PRINT_VAR(yfspole);
	double subb;

	subb = (fsrcount!=1?p_dipoles->CalculateRealSubEEX(kk):p_dipoles->CalculateRealSubEEX(k));
	
	if(p.size()!=(m_flavs.size()+1)){
		msg_Error()<<"Mismatch in "<<METHOD<<std::endl;
	}
	double r = p_realvirt->Calc(p, m_born) / norm;
	if(p_realvirt->FailCut()) m_failcut = true;;
	if (IsBad(r)) {
		msg_Error()<<"Real-Virtual is "<<r<<std::endl;
		return 0;
	}
	if(IsZero(r)) {
		// PRINT_VAR(r);
		// PRINT_VAR(k);
		return 0;
	}
	// m_plab = pp;
	// double aB = subloc*m_oneloop;//CalculateVirtual();//*p_realvirt->m_factor;
	double aB = p_nlodipoles->CalculateRealSub(k)*m_oneloop;//*p_realvirt->m_factor;
	// double aB = subloc*CalculateVirtual()*p_realvirt->m_factor;
	// yfspole*=m_oneloop*p_realvirt->m_factor;
	// double aB = subloc*(p_virt->Calc(pp, m_born) - m_born*p_nlodipoles->CalculateVirtualSub());
	// yfspole*=(p_virt->p_loop_me->ME_E1()*p_virt->m_factor-m_born*p_nlodipoles->Get_E1());
	// double tot = (r-aB) / subloc;
	if(m_submode==submode::local) tot =  (r*flux-aB/m_rescale_alpha)/subloc;
	else if(m_submode==submode::global) tot =  (r*flux-aB/m_rescale_alpha)/subb;
	else if(m_submode==submode::off) tot =  (r*flux)/subb;
	if(m_check_poles==1 && r!=0){
		double pr1 = p_realvirt->p_loop_me->ME_E1()*p_realvirt->m_factor*flux/norm;
		double pr2 = p_realvirt->p_loop_me->ME_E1()*p_realvirt->m_factor;
		// yfspole = p_nlodipoles->Get_E1();
		double yfspoleV = p_dipoles->Get_E1();
		double correctdigit = ::countMatchingDigits(pr2, -p_nlodipoles->Get_E1());
		m_histograms1d["RVSinglePoleCD"]->Insert(correctdigit);
		m_histograms1d["RealLoopEpsLP"]->Insert(log10(fabs(pr1)));
		m_histograms1d["RealLoopEpsYFS"]->Insert(log10(fabs(yfspole)));
		if(!IsEqual(pr2,-yfspole,1e-4)){
			msg_Out()<<"Poles do not cancel in YFS Real-Virtuals"<<std::endl
					 <<"Process =  "<<p_realvirt->p_loop_me->Name()<<std::endl
					 <<"Correct Digits =  "<<correctdigit<<std::endl
					 <<"One-Loop Provider RV eps^{-1}  = "<<pr2<<std::endl
					 <<"Sherpa RV eps^{-1} = "<<-yfspole<<std::endl
					 <<"Sherpa/One-Loop = "<<-yfspole/pr2<<std::endl;
			return 0;
		}
		else{
			msg_Debugging()<<std::setprecision(16)<<"Poles cancel in YFS Real-Virtuals"<<std::endl
						<<"Process =  "<<p_realvirt->p_loop_me->Name()<<std::endl
						<<"Correct Digits =  "<<correctdigit<<std::endl
					 	<<"One-Loop Provider RV eps^{-1}  = "<<pr2<<std::endl
			 			<<"Sherpa RV eps^{-1} = "<<p_nlodipoles->Get_E1()<<std::endl;
		}
	}
	return tot;
}





double NLO_Base::CalculateRealReal() {
	if (!m_rrtool) return 0;
	double rr(0);
	Vec4D_Vector photons;
	for (auto k : m_ISRPhotons) photons.push_back(k);
	for (auto k : m_FSRPhotons) photons.push_back(k);
	if(photons.size()==0) return 0;
	size_t nISR = m_ISRPhotons.size();
	size_t nFSR = m_FSRPhotons.size();
	for (int i = 0; i < photons.size(); ++i) {
		for (int j = i+1; j < photons.size(); ++j) {
			Vec4D k  = photons[i];
			Vec4D kk = photons[j];
			int isFSR_i = (i >= nISR) ? 1 : 0;
			int isFSR_j = (j >= nISR) ? 1 : 0;
			rr+=CalculateRealReal(k,kk,isFSR_i,isFSR_j);
			if(m_check_rr_sub){
				// k*=2;
				// kk*=2;
				if(k.E() < 0.3*sqrt(m_s)) continue;
				if(kk.E() < 0.1*sqrt(m_s)) continue;
				if(!m_failcut) CheckRealRealSub(k, kk, isFSR_i, isFSR_j);
			}
		}
	}
	return rr;
}


double NLO_Base::CalculateRealReal(Vec4D k1, Vec4D k2, int fsr1, int fsr2){
	double norm = 2.*pow(2 * M_PI, 6);
	Vec4D_Vector p(m_plab),pi(m_bornMomenta), pf(m_bornMomenta), plab(m_plab);
 	Vec4D_Vector pp = p;
	Vec4D kk1 = k1;
	Vec4D kk2 = k2;
	// dipoletype::code fluxtype1 = p_nlodipoles->WhichResonant(k1);
	// dipoletype::code fluxtype2 = p_nlodipoles->WhichResonant(k2);
	p_nlodipoles->MakeDipolesII(m_flavs,m_plab,m_plab);
	p_nlodipoles->MakeDipolesIF(m_flavs,m_plab,m_plab);
	p_nlodipoles->MakeDipoles(m_flavs,m_plab,m_plab);
	// fsr1 = (fluxtype1==dipoletype::final?1:0);
	// fsr2 = (fluxtype2==dipoletype::final?1:0);
  if(fsr1 && !fsr2){
  	if(!HasFSR()) msg_Error()<<"Wrong dipole type in "<<METHOD<<endl;
  	for (Dipole_Vector::iterator Dip = p_nlodipoles->GetDipoleFF()->begin();
       Dip != p_nlodipoles->GetDipoleFF()->end(); ++Dip) {
  		 Dip->ClearPhotons();
  		 double scalek = p_fsr->ScalePhoton(k1);
  		 Dip->SetPhotonScale(scalek);
  		 Dip->AddPhotonToDipole(k1);
  		 if(!Dip->BoostNLO()) {
  		 	msg_Error()<<"NLO boost failed" <<endl
  		 						<<"k1 = "<<k1<<endl
  		 						<<"k2 = "<<k2<<endl;
  		 	return 0;
  		 }
  		 int i(0);
  		 for (auto f : Dip->m_flavs) {
      	p[p_nlodipoles->m_flav_label[f]] =  Dip->GetNewMomenta(i);
      	i++;
    	}
    	// k1 = Dip->m_dipolePhotons[0];
  	}
  }
  if(!fsr1 && fsr2){
  	if(!HasFSR()) msg_Error()<<"Wrong dipole type in "<<METHOD<<endl;
  	for (Dipole_Vector::iterator Dip = p_nlodipoles->GetDipoleFF()->begin();
       Dip != p_nlodipoles->GetDipoleFF()->end(); ++Dip) {
  		 double scalek = p_fsr->ScalePhoton(k2);
  		 Dip->SetPhotonScale(scalek);
  		 Dip->AddPhotonToDipole(k2);
  		 if(!Dip->BoostNLO()) {
  		 	msg_Error()<<"NLO boost failed" <<endl
  		 						<<"k1 = "<<k1<<endl
  		 						<<"k2 = "<<k2<<endl;
  		 	return 0;
  		 }
  		 int i(0);
  		 for (auto f : Dip->m_flavs) {
      	p[p_nlodipoles->m_flav_label[f]] =  Dip->GetNewMomenta(i);
      	i++;
    	}
    	// k2 = Dip->m_dipolePhotons[0];
  	}
  }
  if(fsr1 && fsr2){
  	if(!HasFSR()) msg_Error()<<"Wrong dipole type in "<<METHOD<<endl;
  	Dipole_Vector *diplo = p_dipoles->GetDipoleFF();
  	for (Dipole_Vector::iterator Dip = p_nlodipoles->GetDipoleFF()->begin();
       Dip != p_nlodipoles->GetDipoleFF()->end(); ++Dip) {
  		 double scalek = p_fsr->ScalePhoton(k1+k2);
  		 // scalek += p_fsr->ScalePhoton(k2);
  		 Dip->ClearPhotons();
  		 Dip->SetPhotonScale(scalek);
  		 Dip->AddPhotonToDipole(k1);
  		 Dip->AddPhotonToDipole(k2);
  		 if(!Dip->BoostNLO()) {
  		 	msg_Error()<<"NLO boost failed" <<endl
  		 							<<"k1 = "<<k1<<endl
  		 						<<"k2 = "<<k2<<endl;
  		 	return 0;
  		 }
  		 int i(0);
  		 for (auto f : Dip->m_flavs) {
      	p[p_nlodipoles->m_flav_label[f]] =  Dip->GetNewMomenta(i);
      	i++;
    	}
    	// k1 = Dip->m_dipolePhotons[0];
    	// k2 = Dip->m_dipolePhotons[1];
  	}
  }

  if(!fsr1 && !fsr2) MapMomenta(p, k1, k2);
 	p.push_back(k1);
 	p.push_back(k2);
 	if(fsr1 || fsr2) MapInitial(p);
 	CheckMasses(p, 2);
  Vec4D_Vector _p = p;
 	_p.pop_back();
 	_p.pop_back();
  m_plab = pp;
	p_nlodipoles->MakeDipolesII(m_flavs,_p,m_plab);
	p_nlodipoles->MakeDipolesIF(m_flavs,_p,m_plab);
	p_nlodipoles->MakeDipoles(m_flavs,_p,m_plab);
	double subloc1 = p_nlodipoles->CalculateRealSub(k1);
	double subloc2 = p_nlodipoles->CalculateRealSub(k2);
	double flux;
	if(m_flux_mode==1) flux = p_nlodipoles->CalculateFlux(k1)*p_nlodipoles->CalculateFlux(k2);
	// if(m_flux_mode==1) flux = p_nlodipoles->CalculateFlux(k1,k2);
	else flux = p_dipoles->CalculateFlux(k1)*p_dipoles->CalculateFlux(k2);
	double tot,rcoll;
	if(!CheckMomentumConservation(p)) {
		return 0;
	}
	// return 0;
	double r = p_realreal->Calc_R(p) / norm ;
	if(p_realreal->FailCut()) {
		m_failcut = true;
		return 0;
	}
	if(IsBad(r) || IsBad(flux)) {
		msg_Error()<<"Bad point for YFS Real"<<std::endl
							 <<"Real ME is : "<<r<<std::endl
							 <<"Flux is : "<<flux<<std::endl;
		return 0;
	}
	m_recola_evts+=1;
	double real1 = CalculateReal(kk1,3+fsr1);
	double real2 = CalculateReal(kk2,3+fsr2);
	double sub1 = (fsr1!=1?p_dipoles->CalculateRealSubEEX(kk1):p_dipoles->CalculateRealSubEEX(kk1));
	double sub2 = (fsr2!=1?p_dipoles->CalculateRealSubEEX(kk2):p_dipoles->CalculateRealSubEEX(kk2));
	double fullsub = (-subloc2*real1 -subloc1*real2-subloc1*subloc2*m_born);
	tot = (r*flux + fullsub)/sub1/sub2;
  if(IsBad(tot)){
  	msg_Error()<<"NNLO RR is NaN"<<std::endl;
  }
  m_plab = plab;
	return tot;
}


double NLO_Base::CalculateVV() {
	if (!m_vvtool) return 0;
	if(m_eex_virt){
		return p_dipoles->CalculateEEXVirtual()*m_born-m_born;
	}
	double virt;
	double sub;
	// CheckMassReg();
	if(!HasISR()) virt = p_vv->Calc(m_bornMomenta, m_born);
	else virt = p_vv->Calc(m_plab, m_born);
	if(m_check_virt_born) {
			if (!IsEqual(m_born, p_virt->p_loop_me->ME_Born(), 1e-6)) {
			msg_Error() << METHOD << "\n Warning! Loop provider's born is different! YFS Subtraction likely fails\n"
									<< "Loop Provider " << ":  "<<p_virt->p_loop_me->ME_Born()
									<< "\nSherpa" << ":  "<<m_born<<std::endl
									<<"PhaseSpace Point = ";
			for(auto _p: m_plab) msg_Error()<<_p<<std::endl;
		}
	}	
	if(p_vv->FailCut()) return 0;
	if(m_virt_sub && p_virt->p_loop_me->Mode()!=1) sub = p_dipoles->CalculateVirtualSub();
	else sub = 0;
	double sub2 = p_dipoles->CalculateVVSubEps();
	// m_oneloop = (virt- sub * m_born/m_rescale_alpha );
	m_oneloop = (virt- sub * CalculateVirtual()/m_rescale_alpha -0.5*sub*sub*m_born/m_rescale_alpha);
	if(p_virt->p_loop_me->Mode()==1) {
		m_oneloop /= m_rescale_alpha; 
		// PRINT_VAR(m_rescale_alpha);
	}
	if(IsBad(m_oneloop) || IsBad(sub)){
		msg_Error()<<"YFS Virtual is NaN"<<std::endl
							 <<"Virtual:  "<<virt<<std::endl
							 <<"Subtraction: "<<sub*m_born<<std::endl
							 <<"PhaseSpace Point: "<<std::endl<<m_plab<<std::endl;
	}
	double loope1 = p_vv->p_loop_me->ME_E1()*p_vv->m_factor;//*p_vv->m_factor;//+p_virt->p_loop_me->ME_E1()*p_virt->m_factor;;
	double loope2 = 2.*p_vv->p_loop_me->ME_E2()*p_vv->m_factor*p_vv->m_factor;
	double yfse1 = p_dipoles->Get_E1();
	double yfse2 = p_dipoles->GetVV_E2();
	// PRINT_VAR(::countMatchingDigits(loope1, yfse1));
	// PRINT_VAR(::countMatchingDigits(loope2, -yfse2));
	// PRINT_VAR(loope1);
	// PRINT_VAR(loope1);
	// PRINT_VAR(loope1/yfse1);
	// PRINT_VAR(loope2);
	// PRINT_VAR(yfse2);
	// PRINT_VAR(loope2/yfse2);
	// PRINT_VAR(m_born);
	// PRINT_VAR(1./m_born);
	// PRINT_VAR(p_dipoles->Get_E1()+0.5*pow(p_dipoles->Get_E1(),2));
	if(m_check_poles==1){
		if(m_virt_sub==0) sub = p_dipoles->CalculateVirtualSub();
		const double p1 =  p_vv->p_loop_me->ME_E1()*p_vv->m_factor;
		const double p2 =  2.*p_vv->p_loop_me->ME_E2()*p_vv->m_factor*p_vv->m_factor;
		const double yfspole1 =  (p_dipoles->Get_E1());
		const double yfspole2 =  p_dipoles->GetVV_E2();
		PRINT_VAR(p1/yfspole1);
		int ncorrect1 = ::countMatchingDigits(p1, yfspole1, 32);
		int ncorrect2 = ::countMatchingDigits(p2, -yfspole2, 32);
		if(!IsEqual(p2,-yfspole2,1e-6) || ncorrect1 < 10){
			msg_Error()<<"Poles do not cancel in YFS Double Virtuals"<<std::endl
					 <<"Correct digits \epsion^{-1} =  "<<ncorrect1<<std::endl
					 <<"Correct digits \epsion^{-2} =  "<<ncorrect2<<std::endl;
					 return 0;
		}
		else{
			int i = 0;
			msg_Debugging()<<std::setprecision(32);
			msg_Out()<<"Poles cancel in YFS double Virtuals to "<<ncorrect2<<" digits"<<std::endl;
			m_histograms1d["SinglePoleVV"]->Insert(ncorrect1);
			m_histograms1d["DoublePoleVV"]->Insert(ncorrect2);
		}
	}
	return 0;
}


void NLO_Base::RandomRotate(Vec4D &p){
  Vec4D t1 = p;
  // rotate around x
  p[2] = cos(m_ranTheta)*t1[2] - sin(m_ranTheta)*t1[3];
  p[3] = sin(m_ranTheta)*t1[2] + cos(m_ranTheta)*t1[3];
  t1 = p;
  // rotate around z
  p[1] = cos(m_ranPhi)*t1[1]-sin(m_ranPhi)*t1[2];
  p[2] = sin(m_ranPhi)*t1[1]+cos(m_ranPhi)*t1[2];
}

void NLO_Base::MapMomenta(Vec4D_Vector &p, Vec4D &k) {
	Vec4D Q;
	Vec4D QQ, PP;
  double s = (m_plab[0]+m_plab[1]).Abs2();
  double t = (m_plab[0]-m_plab[2]).Abs2();
  m_ranTheta = acos(1.+2.*t/s);
	m_ranPhi = ran->Get()*2.*M_PI;
	Poincare boostLab(m_bornMomenta[0]+m_bornMomenta[1]);
	for (int i = 2; i < p.size(); ++i)
	{
		Q += p[i];
	}
	Q += k;
	double sq = Q.Abs2();
	Poincare boostQ(Q);
  Poincare pRot(m_bornMomenta[0], Vec4D(0., 0., 0., 1.));
	for (int i = 0; i < p.size(); ++i) {
		pRot.RotateBack(p[i]);
		boostQ.Boost(p[i]);
		// RandomRotate(p[i]);
	}
	pRot.RotateBack(k);
	boostQ.Boost(k);
	// RandomRotate(k);
	double qx(0), qy(0), qz(0);
	for (int i = 2; i < p.size(); ++i)
	{
		qx += p[i][1];
		qy += p[i][2];
		qz += p[i][3];
	}
	if (!IsEqual(k[1], -qx, 1e-5) || !IsEqual(k[2], -qy, 1e-5) || !IsEqual(k[3], -qz, 1e-5) ) {
		if( k[1]> 1e-6 && k[2]> 1e-6 && k[3]> 1e-6 ){
			msg_Error() << "YFS Mapping has failed for ISR\n";
			msg_Error() << " Photons px = " << k[1] << "\n Qx = " << -qx << std::endl;
			msg_Error() << " Photons py = " << k[2] << "\n Qy = " << -qy << std::endl;
			msg_Error() << " Photons pz = " << k[3] << "\n Qz = " << -qz << std::endl;
		}
		}
	for (int i = 2; i < p.size(); ++i)
	{
		QQ += p[i];
	}
	QQ+=k;
	double sqq = QQ.Abs2();
	if (!IsEqual(sqq, sq, 1e-6))
	{
		msg_Error() << "YFS Real mapping not conserving momentum in " << METHOD << std::endl;
	}
	// if(m_is_isr) QQ = p[0]+p[1];
  // double zz = sqrt(sqq) / 2.;
	// double z = zz * sqrt((sqq - sqr(m_flavs[0].Mass() - m_flavs[1].Mass())) * (sqq - sqr(m_flavs[0].Mass() + m_flavs[1].Mass()))) / sqq;
	double sign_z = (p[0][3] < 0 ? -1 : 1);
	// p[0] = {zz, 0, 0, z};
	// p[1] = {zz, 0, 0, -z};
  double m1 = m_flavs[0].Mass();
  double m2 = m_flavs[1].Mass();
  double lamCM = 0.5*sqrt(Lambda(sqq,m1*m1,m2*m2)/sqq);
  double E1 = lamCM*sqrt(1+m1*m1/sqr(lamCM));
  double E2 = lamCM*sqrt(1+m2*m2/sqr(lamCM));
 	p[0] = {E1, 0, 0, sign_z*lamCM};
  p[1] = {E2, 0, 0, -sign_z*lamCM};
  Poincare pRot2(m_bornMomenta[0], Vec4D(0., 	0., 0, 1.));
	for (int i = 0; i < p.size(); ++i)
	{
		pRot2.Rotate(p[i]);
		boostLab.BoostBack(p[i]);
	}
	pRot2.Rotate(k);
	boostLab.BoostBack(k);
}



void NLO_Base::MapMomenta(Vec4D_Vector &p, Vec4D &k1, Vec4D &k2) {
	Vec4D Q;
	Vec4D QQ, PP;
  double s = (m_plab[0]+m_plab[1]).Abs2();
  double t = (m_plab[0]-m_plab[2]).Abs2();
  m_ranTheta = acos(1.+2.*t/s);
	m_ranPhi = ran->Get()*2.*M_PI;
	Poincare boostLab(m_bornMomenta[0]+m_bornMomenta[1]);
	for (int i = 2; i < p.size(); ++i)
	{
		Q += p[i];
	}
	Q += k1+k2;
	double sq = Q.Abs2();
	Poincare boostQ(Q);
  Poincare pRot(m_bornMomenta[0], Vec4D(0., 0., 0., 1.));
	for (int i = 0; i < p.size(); ++i) {
		pRot.Rotate(p[i]);
		boostQ.Boost(p[i]);
		// RandomRotate(p[i]);
	}
	// pRot.Rotate(k);
	boostQ.Boost(k1);
	boostQ.Boost(k2);
	// RandomRotate(k);
	double qx(0), qy(0), qz(0);
	for (int i = 2; i < p.size(); ++i)
	{
		qx += p[i][1];
		qy += p[i][2];
		qz += p[i][3];
	}
	// if (!IsEqual(k[1], -qx, 1e-5) || !IsEqual(k[2], -qy, 1e-5) || !IsEqual(k[3], -qz, 1e-5) ) {
	// 	if( k[1]> 1e-6 && k[2]> 1e-6 && k[3]> 1e-6 ){
	// 		msg_Error() << "YFS Mapping has failed for ISR\n";
	// 		msg_Error() << " Photons px = " << k[1] << "\n Qx = " << -qx << std::endl;
	// 		msg_Error() << " Photons py = " << k[2] << "\n Qy = " << -qy << std::endl;
	// 		msg_Error() << " Photons pz = " << k[3] << "\n Qz = " << -qz << std::endl;
	// 	}
	// 	}
	for (int i = 2; i < p.size(); ++i)
	{
		QQ += p[i];
	}
	QQ+=k1+k2;
	double sqq = QQ.Abs2();
	if (!IsEqual(sqq, sq, 1e-6))
	{
		msg_Error() << "YFS Real mapping not conserving momentum in " << METHOD << std::endl;
	}
	// if(m_is_isr) QQ = p[0]+p[1];
  // double zz = sqrt(sqq) / 2.;
	// double z = zz * sqrt((sqq - sqr(m_flavs[0].Mass() - m_flavs[1].Mass())) * (sqq - sqr(m_flavs[0].Mass() + m_flavs[1].Mass()))) / sqq;
	double sign_z = (p[0][3] < 0 ? -1 : 1);
	// p[0] = {zz, 0, 0, z};
	// p[1] = {zz, 0, 0, -z};
  double m1 = m_flavs[0].Mass();
  double m2 = m_flavs[1].Mass();
  double lamCM = 0.5*sqrt(Lambda(sqq,m1*m1,m2*m2)/sqq);
  double E1 = lamCM*sqrt(1+m1*m1/sqr(lamCM));
  double E2 = lamCM*sqrt(1+m2*m2/sqr(lamCM));
 	p[0] = {E1, 0, 0, sign_z*lamCM};
  p[1] = {E2, 0, 0, -sign_z*lamCM};
  Poincare pRot2(m_bornMomenta[0], Vec4D(0., 	0., 0, 1.));
	for (int i = 0; i < p.size(); ++i)
	{
		pRot2.Rotate(p[i]);
		boostLab.BoostBack(p[i]);
	}
	pRot2.Rotate(k1);
	pRot2.Rotate(k2);
	boostLab.BoostBack(k1);
	boostLab.BoostBack(k2);
}


void NLO_Base::MapInitial(Vec4D_Vector &p){
	Vec4D QQ;
	Vec4D_Vector born=p;
	for (int i = 2; i < p.size(); ++i)
	{
		QQ += p[i];
	}
	double sqq = QQ.Abs2();
	double sign_z = (p[0][3] > 0 ? -1 : 1);
  double m1 = m_flavs[0].Mass();
  double m2 = m_flavs[1].Mass();
  double lamCM = 0.5*sqrt(Lambda(sqq,m1*m1,m2*m2)/sqq);
  double E1 = lamCM*sqrt(1+m1*m1/sqr(lamCM));
  double E2 = lamCM*sqrt(1+m2*m2/sqr(lamCM));
 	p[0] = {E1, 0, 0, sign_z*lamCM};
  p[1] = {E2, 0, 0, -sign_z*lamCM};
	Poincare boostLab(QQ);
  Poincare pRot = Poincare(p[0], Vec4D(0., 0., 0., 1.));
	for (int i = 0; i < 2; ++i)
	{
		pRot.Rotate(p[i]);
		boostLab.BoostBack(p[i]);
		// pRot2.Rotate(p[i]);
	}
}


void NLO_Base::CheckMasses(Vec4D_Vector &p, int realmode){
	bool allonshell=true;
	std::vector<double> masses;
	Flavour_Vector flavs = m_flavs;
	if(realmode>=1) flavs.push_back(Flavour(kf_photon));
	if(realmode>=2) flavs.push_back(Flavour(kf_photon));
	if(p.size() != flavs.size()) msg_Error()<<"Mismatch between mass and flavour vectors in "<<METHOD<<std::endl;
	for (int i = 0; i < p.size(); ++i)
	{
		masses.push_back(flavs[i].Mass());
		if(!IsEqual(p[i].Mass(),flavs[i].Mass())&& flavs[i].Mass()!=0){
			msg_Debugging()<<"Wrong particle masses in YFS Mapping"<<std::endl
								 <<"Flavour = "<<flavs[i]<<", with mass = "<<flavs[i].Mass()<<std::endl
								 <<"Four momentum = "<<p[i]<<", with mass = "<<p[i].Mass()<<std::endl;
			allonshell = false;

		}
	}
	if(!allonshell) {
		m_stretcher.StretchMomenta(p, masses);
		for (int i = 0; i < p.size(); ++i){
			msg_Debugging()<<"Flavour = "<<flavs[i]<<", with mass = "<<flavs[i].Mass()<<std::endl
							 <<"Four momentum = "<<p[i]<<", with new mass = "<<p[i].Mass()<<std::endl;
		}
	}
}

void NLO_Base::RescaleMasses(Vec4D_Vector &p, std::vector<double> masses){
	bool allonshell=true;
	if(p.size() != masses.size()) msg_Error()<<"Mismatch between mass and vectors in "<<METHOD<<std::endl;
	m_stretcher.StretchMomenta(p, masses);
	// return true;
}


bool NLO_Base::CheckPhotonForReal(const Vec4D &k) {
	for (int i = 0; i < m_plab.size(); ++i)
	{
		if (m_flavs[i].IsChargedLepton()) {
			double sik = (k + m_plab[i]).Abs2();
			if (sik < m_hardmin ) {
				// msg_Out()<<"Rejecting photon k = "<<k<<std::endl
				// 				<<"sik = "<<sik<<std::endl;
				return false;
			}
			// if(m_plab[i].PPerp()< m_hardmin) return false;
		}
	}
	if(k.PPerp() < 1e-4) return false;
	return true;
}


bool NLO_Base::CheckPhotonForReal(const Vec4D &k, const Vec4D_Vector &p) {
	for (int i = 0; i < p.size(); ++i)
	{
		if (m_flavs[i].IsChargedLepton()) {
			double sik = (k + p[i]).Abs2();
			if (sik < m_hardmin*p[i].Abs2() ) {
				msg_Out()<<"Rejecting photon k = "<<k<<std::endl
								<<"sik = "<<sik<<std::endl;
				return false;
			}
			// if(p[i].PPerp() < m_hardmin) return false;
		}
	}
	// if(k.PPerp() < m_hardmin) return false;
	return true;
}

bool NLO_Base::CheckMomentumConservation(Vec4D_Vector p){
  Vec4D incoming = p[0]+p[1];
  Vec4D outgoing;
  for (int i = 2; i < p.size(); ++i)
  {
    if(p[i].E() < 0 || IsBad(p[i].E())) {
    	msg_Error()<<"Energy less than zero!: "<<p[i]<<std::endl;
    	return false;
    }
    outgoing+=p[i];
  }
  Vec4D diff = incoming - outgoing;
  if(!IsEqual(incoming,outgoing, 1e-8)){
    msg_Error()<<METHOD<<std::endl<<"Momentum not conserverd in YFS NLO"<<std::endl
               <<"Incoming momentum = "<<incoming<<std::endl
               <<"Outgoing momentum = "<<outgoing<<std::endl
               <<"Difference = "<<diff<<std::endl
               <<"Vetoing Event "<<std::endl;
    return false;
  }
  return true;
}

void NLO_Base::CheckMassReg(){
	double virt;
	if (m_check_mass_reg==1 && !m_realvirt) {
		out_sub.open("yfs-sub.txt", std::ios_base::app);
		out_recola.open("virtual-res.txt", std::ios_base::app); // append instead of overwrite
		out_finite.open("yfs-finite.txt", std::ios_base::app);
		if(!HasISR()) virt = p_virt->Calc(m_bornMomenta, m_born);
		else virt = p_virt->Calc(m_plab, m_born);
		if (!IsEqual(m_born, p_virt->p_loop_me->ME_Born()*m_rescale_alpha, 1e-6)) {
			msg_Error() << METHOD << "\n Warning! Loop provider's born is different! YFS Subtraction likely fails\n"
									<< "Loop Provider " << ":  "<<p_virt->p_loop_me->ME_Born()
									<< "Sherpa" << ":  "<<m_born;
		}
		double sub = p_dipoles->CalculateVirtualSub();
		std::cout << setprecision(15);
		out_sub<< setprecision(15) << m_photonMass << "," << -sub*m_born/m_rescale_alpha << std::endl;
		out_recola<< setprecision(15) << m_photonMass << "," << virt << std::endl;
		out_finite<< setprecision(15) << m_photonMass << "," << virt - sub*m_born/m_rescale_alpha << std::endl;
		out_sub.close();
		out_recola.close();
		exit(0);
	}
}


void NLO_Base::CheckRealSub(Vec4D k, int mode){
		// if(k.E() < 20) return;
		// k*=100;
		double real;
		std::string filename="Real_subtracted_";
		std::string filename1="Sub_term_";
		std::string filename2="Real_ME_";
		for(auto f: m_flavs) {
			filename+=f.IDName();
			filename+="_";
			filename1+=f.IDName();
			filename1+="_";
			filename2+=f.IDName();
			filename2+="_";
		}
		filename+=".txt";
		filename1+=".txt";
		filename2+=".txt";
		if(ATOOLS::FileExists(filename))  ATOOLS::Remove(filename);
		if(ATOOLS::FileExists(filename1))  ATOOLS::Remove(filename1);
		if(ATOOLS::FileExists(filename2))  ATOOLS::Remove(filename2);
		out_finite.open(filename, std::ios_base::app);
		out_sub.open(filename1, std::ios_base::app);
		out_real.open(filename2, std::ios_base::app);
		for (double i = 1; i < 20 ; i+=0.1)
		{
			k=k/i;
			real=CalculateReal(k,mode);
			if(k.E() <= 1e-16 || real==0) break;
			out_finite<<k.E()<<","<<fabs(real)/m_born<<std::endl;
			out_real<<k.E()<<","<<(m_real)*p_nlodipoles->CalculateFlux(k)<<std::endl;
			out_sub<<k.E()<<","<<m_subloc*m_born/m_rescale_alpha<<std::endl;
		}
		out_finite.close();
		out_real.close();
		out_sub.close();
		exit(0);
}

void NLO_Base::CheckRealVirtualSub(Vec4D k){
		// if(k.E() < 20) return;
		// k*=100;
		double real;
		std::string filename="RealVirtual_subtracted_";
		for(auto f: m_flavs) {
			filename+=f.IDName();
			filename+="_";
		}
		filename+=".txt";
		if(ATOOLS::FileExists(filename))  ATOOLS::Remove(filename);
		out_sub.open(filename, std::ios_base::app);
		// if(k.E() < 0.8*sqrt(m_s)/2.) return;
		for (double i = 1; i < 20 ; i+=0.005)
		{
			k=k/i;
			real=CalculateRealVirtual(k);
			PRINT_VAR(real);
			out_sub<<k.E()<<","<<fabs(real)/m_born<<std::endl;
			if(k.E() < m_isrcut ) break;
			// m_histograms2d["Real_me_sub"]->Insert(k.E(),fabs(real), 1);
		}
		out_sub.close();
		exit(0);
}

void NLO_Base::CheckRealRealSub(Vec4D k1, Vec4D k2, int fsr1, int fsr2){
		// if(k.E() < 20) return;
		// k*=100;
		double real;
		Vec4D _k1 = k1;
		Vec4D _k2 = k2;
		std::string filename1="RealReal_k1_subtracted_";
		std::string filename2="RealReal_k2_subtracted_";
		std::string filename3="RealReal_k1_k2_subtracted_";
		for(auto f: m_flavs) {
			filename1+=f.IDName();
			filename2+=f.IDName();
			filename3+=f.IDName();
			filename1+="_";
			filename2+="_";
			filename3+="_";
		}
		filename1+=".txt";
		filename2+=".txt";
		filename3+=".txt";
		if(ATOOLS::FileExists(filename1))  ATOOLS::Remove(filename1);
		if(ATOOLS::FileExists(filename2))  ATOOLS::Remove(filename2);
		if(ATOOLS::FileExists(filename3))  ATOOLS::Remove(filename3);
		out_sub.open(filename1, std::ios_base::app);
		// if(k.E() < 0.8*sqrt(m_s)/2.) return;
		for (double i = 1; i < 20 ; i+=0.02)
		{
			k1=k1/i;
			// if(k1.E()< m_isrcut*sqrt(m_s)) break;
			real=CalculateRealReal(k1,k2,fsr1, fsr2);
			out_sub<<k1.E()<<","<<fabs(real)/m_born<<std::endl;
			if(k1.E() <= 1e-16) break;
			// m_histograms2d["Real_me_sub"]->Insert(k.E(),fabs(real), 1);
		}
		out_sub.close();
		out_sub.open(filename2, std::ios_base::app);
	 	k2 = _k2;
	 	k1 = _k1;
		for (double i = 1; i < 20 ; i+=0.02)
		{
			k2=k2/i;
			// if(k2.E() <= 1e-16 ) break;
			real=CalculateRealReal(k1,k2, fsr1, fsr2);
			out_sub<<k2.E()<<","<<fabs(real)/m_born<<std::endl;
			if(k2.E()< m_isrcut*sqrt(m_s)) break;
			// if(IsZero(real)) break;
			// m_histograms2d["Real_me_sub"]->Insert(k.E(),fabs(real), 1);
		}
		out_sub.close();
		out_sub.open(filename3, std::ios_base::app);
		k2 = _k2;
	 	k1 = _k1;
		for (double i = 1; i < 20 ; i+=0.02)
		{
			k2=k2/i;
			k1=k1/i;
			real=CalculateRealReal(k1,k2,fsr1,fsr2);
			out_sub<<k1.E()<<","<<fabs(real)/m_born<<std::endl;
			if(k1.E() <= 1e-16 || real==0 && !m_failcut) break;
			// if(IsZero(real)) break;
			// m_histograms2d["Real_me_sub"]->Insert(k.E(),fabs(real), 1);
		}
		out_sub.close();
		exit(0);
}
