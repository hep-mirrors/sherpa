#include "ATOOLS/Math/Random.H"
#include "YFS/NLO/NLO_Base.H"
#include "MODEL/Main/Running_AlphaQED.H"


using namespace YFS;
using namespace MODEL;
using namespace ATOOLS;
using namespace std;

double massmin=2220;

std::ofstream out_recola;
std::ofstream out_sub, out_real, out_finite;

double SqLam(double x,double y,double z)
{
  return abs(x*x+y*y+z*z-2.*x*y-2.*x*z-2.*y*z);
}

NLO_Base::NLO_Base() {
  p_yfsFormFact = new YFS::YFS_Form_Factor();
  p_nlodipoles = new YFS::Define_Dipoles();
  m_evts = 0;
  m_recola_evts = 0;
  m_realtool = 0;
  m_realvirt = 0;
  m_looptool = 0;
  m_rrtool = 0;
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
  	if (!ATOOLS::DirectoryExists(m_debugDIR_NLO)) {
			ATOOLS::MakeDir(m_debugDIR_NLO);
		}
  }
}

NLO_Base::~NLO_Base() {
  if(m_isr_debug || m_fsr_debug || m_check_real_sub){
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
	// if(p_real) delete p_real;
	// if(p_virt) delete p_virt;
}


void NLO_Base::InitializeVirtual(const PHASIC::Process_Info& pi) {
	if(!m_eex_virt){
		p_virt = new YFS::Virtual(pi);
		m_looptool = true;
	}
}

void NLO_Base::InitializeReal(const PHASIC::Process_Info& pi) {
	if(m_coll_real) return;
	p_real = new YFS::Real(pi);
	m_realtool = true;
}

void NLO_Base::InitializeRealVirtual(const PHASIC::Process_Info& pi) {
	p_realvirt = new YFS::RealVirtual(pi);
	m_realvirt = true;
	// m_looptool = true;
}

void NLO_Base::InitializeRealReal(const PHASIC::Process_Info& pi) {
	p_realreal = new YFS::RealReal(pi);
	m_rrtool = true;
}


void NLO_Base::Init(Flavour_Vector &flavs, Vec4D_Vector &plab, Vec4D_Vector &born) {
	m_flavs = flavs;
	m_plab = plab;
	m_bornMomenta = born;
}


double NLO_Base::CalculateNLO() {
	double result{0.0};
	if(!m_real_only) result += CalculateVirtual();
	if(!m_virtual_only) result += CalculateReal();
	result += CalculateVirtual();
	result += CalculateReal();
	result += CalculateRealVirtual();
	result += CalculateRealReal();
	return result;
}


double NLO_Base::CalculateVirtual() {
	if(m_eex_virt){
		// PRINT_VAR(p_dipoles->CalculateEEXVirtual());
		// subtract born to avoid double counting
		// already present in eex!!
		return p_dipoles->CalculateEEXVirtual()*m_born-m_born;
	}
	if (!m_looptool && !m_realvirt) return 0;
	double virt;
	double sub;
	p_dipoles->p_yfsFormFact->p_virt = p_virt;
	CheckMassReg();
	// for(auto pp: m_plab) PRINT_VAR(pp.Mass());
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
	if(virt==0) return 0;
	if(m_virt_sub) sub = p_dipoles->CalculateVirtualSub();
	else {
		sub = 0;
		// virt=virt;
	}
	m_oneloop = (virt- sub * m_born/m_rescale_alpha );
	if(IsBad(m_oneloop)){
		msg_Error()<<"YFS Virtual is NaN"<<std::endl
							 <<"Virtual:  "<<m_oneloop<<std::endl
							 <<"Subtraction: "<<sub*m_born<<std::endl;
	}
	if(m_check_poles==1){
		if(m_virt_sub==0) sub = p_dipoles->CalculateVirtualSub();
		double p1 = p_virt->p_loop_me->ME_E1()*p_virt->m_factor;
		double yfspole = p_dipoles->Get_E1();
		if(!IsEqual(p1,-yfspole,1e-6)){
			msg_Error()<<"Poles do not cancel in YFS Virtuals"<<std::endl
								 <<"One-Loop Provider  = "<<p1<<std::endl
								 <<"Sherpa  = "<<yfspole<<std::endl;
		}
		{
			msg_Debugging()<<"Poles cancel in YFS Virtuals"<<std::endl
								 		 <<"One-Loop Provider  = "<<p1<<std::endl
								 			<<"Sherpa  = "<<yfspole<<std::endl;
		}
	}
	return m_oneloop;
}


double NLO_Base::CalculateReal() {
	if (!m_realtool) return 0;
	double real(0);
	m_ksum*=0;
	for(auto &k: m_FSRPhotons) m_ksum += k;
	for(auto &k: m_ISRPhotons) m_ksum += k;
	// msg->SetPrecision(16);
	double sp = p_dipoles->GetDipoleII()->Sprime();
	if(m_coll_real) return p_dipoles->CalculateEEX()*m_born;
	double collreal = p_dipoles->CalculateEEX()*m_born;
	for (auto k : m_ISRPhotons) {
		if(m_check_real_sub) {
				CheckRealSub(k);
		}
		if(m_isr_debug || m_fsr_debug){
			double fullreal = CalculateReal(k);
			real+=fullreal;
			double coll = p_dipoles->GetDipoleII()[0].Beta1(k);
			coll /=  p_dipoles->GetDipoleII()[0].Eikonal(k);
			if(fullreal!=0) m_histograms2d["REAL_COLL_RATIO"]->Insert(k.E(), coll*m_born/fullreal);
		}
		else real+=CalculateReal(k,0);
	}
	int fsrcount(0);
	for (auto k : m_FSRPhotons) {
		if(m_check_real_sub) {
			if(k.E() < 0.2*sqrt(m_s)) continue;
				CheckRealSub(k);
		}
		real+=CalculateReal(k, 1);
		fsrcount++;
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
  // if(fluxtype==dipoletype::final){
  if(fsrcount){
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
    	k = Dip->m_dipolePhotons[0];
  	}
  }
 	else {
 		MapMomenta(p,k);
 	}
 	p.push_back(k);
 	if(fsrcount) MapInitial(p);
 	Vec4D_Vector pp = p;
 	pp.pop_back();
	p_nlodipoles->MakeDipolesII(m_flavs,pp,m_plab);
	p_nlodipoles->MakeDipolesIF(m_flavs,pp,m_plab);
	p_nlodipoles->MakeDipoles(m_flavs,pp,m_plab);
	double flux;
	if(m_flux_mode==1) flux = p_nlodipoles->CalculateFlux(k);
	else if(m_flux_mode==2) flux = 0.5*(p_nlodipoles->CalculateFlux(kk)+p_nlodipoles->CalculateFlux(k));
	else flux = p_dipoles->CalculateFlux(kk);
	double tot,rcoll;
	double subloc = p_nlodipoles->CalculateRealSub(k);
	double subb;
	if(fsrcount) subb = p_dipoles->CalculateRealSubEEX(k);
	else subb = p_dipoles->CalculateRealSubEEX(k);
	if(IsZero(subb)) return 0;
	if(!CheckMomentumConservation(p)) {
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
	// CheckMasses(p,1);
	double r = p_real->Calc_R(p) / norm;
	if(IsZero(r)) return 0;
	if(IsBad(r) || IsBad(flux)) {
		msg_Error()<<"Bad point for YFS Real"<<std::endl
							 <<"Real ME is : "<<r<<std::endl
							 <<"Flux is : "<<flux<<std::endl;
		return 0;
	}
	m_recola_evts+=1;
	// if(!fsrcount) r*=flux;
	// PRINT_VAR(m_born);
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
	}
	return tot;
}

double NLO_Base::CalculateRealVirtual() {
	if (!m_realvirt) return 0;
	double real(0);
	for (auto k : m_ISRPhotons) {
		if(k.PPerp() > m_hardmin)	real+=CalculateRealVirtual(k,0);
		else real+=m_oneloop*CalculateReal(k);
	}
	for (auto k : m_FSRPhotons) {
		if(k.PPerp() > m_hardmin) real+=CalculateRealVirtual(k, 1);
		else real+=m_oneloop*CalculateReal(k,1);
	}
	return real;
}



double NLO_Base::CalculateRealVirtual(Vec4D k, int fsrcount) {
	if (!m_realvirt) return 0;
	Vec4D_Vector p(m_plab),pi(m_bornMomenta), pf(m_bornMomenta);
	double tot(0), sub(0);
	double norm = 2.*pow(2 * M_PI, 3);
	double flux;
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
  		 if(!Dip->BoostNLO()) return 0;
  		 int i(0);
  		 for (auto f : Dip->m_flavs) {
      	p[p_nlodipoles->m_flav_label[f]] =  Dip->GetNewMomenta(i);
      	i++;
    	}
    	k = Dip->m_dipolePhotons[0];
  	}
  }
 	else {
 		MapMomenta(p,k);
 	}
 	p.push_back(k);
 	if(fsrcount) MapInitial(p);
	Vec4D_Vector pp = p;
 	pp.pop_back();
	p_nlodipoles->MakeDipolesII(m_flavs,pp,m_plab);
	p_nlodipoles->MakeDipolesIF(m_flavs,pp,m_plab);
	p_nlodipoles->MakeDipoles(m_flavs,pp,m_plab);
	p.push_back(k);

	if(m_flux_mode==1) flux = p_nlodipoles->CalculateFlux(k);
	else if(m_flux_mode==2) flux = 0.5*(p_nlodipoles->CalculateFlux(kk)+p_nlodipoles->CalculateFlux(k));
	else flux = p_dipoles->CalculateFlux(kk);

	double subloc = p_nlodipoles->CalculateRealSub(k);
	double subb;

	if(fsrcount) subb = p_dipoles->CalculateRealSubEEX(k);
	else subb = p_dipoles->CalculateRealSubEEX(kk);

	double r = p_realvirt->Calc(p, m_born) / norm;
	if (r == 0 || IsBad(r)) return 0;
	double aB = subloc*CalculateVirtual();
	// double tot = (r-aB) / subloc;
	if(m_submode==submode::local) tot =  (r*flux-aB)/subloc;
	else if(m_submode==submode::global) tot =  (r*flux-aB)/subb;
	else if(m_submode==submode::off) tot =  (r*flux)/subb;
	// real += tot;

	return tot;
}





double NLO_Base::CalculateRealReal() {
	if (!m_rrtool) return 0;
	double real(0), sub(0);
	Vec4D_Vector photons, p(m_plab);
	for (auto k : m_ISRPhotons) photons.push_back(k);
	for (auto k : m_FSRPhotons) photons.push_back(k);
	Vec4D Q = m_bornMomenta[0]+m_bornMomenta[1];
	// Vec4D_Vector photons;
	double len = m_ISRPhotons.size();
	double flux;
	for (int i = 0; i < photons.size(); ++i) {
		for (int j = 0; j < i; ++j) {
			Vec4D k  = photons[i];
			Vec4D kk = photons[j];
			if(m_check_rr_sub){
				if(k.E() < 0.2*sqrt(m_s) && kk.E() < 0.2*sqrt(m_s)) continue;
				CheckRealRealSub(k,kk);
			}
			real+=CalculateRealReal(k,kk, i>(len-1)?1:0, j>(len-1)?1:0);
			// real+=CalculateRealReal(k,kk, 0, 1);
		}
	}
	return real;
}


double NLO_Base::CalculateRealReal(Vec4D k1, Vec4D k2, int fsr1, int fsr2){
	double norm = 2.*pow(2 * M_PI, 3);
	Vec4D_Vector p(m_plab),pi(m_bornMomenta), pf(m_bornMomenta);
	Vec4D kk1 = k1;
	Vec4D kk2 = k2;
	double real1 = CalculateReal(k1);
	double real2 = CalculateReal(k2);
	double sub1 = p_dipoles->CalculateRealSubEEX(k1);
	double sub2 = p_dipoles->CalculateRealSubEEX(k2);
	p_nlodipoles->MakeDipoles(m_flavs,m_plab,m_plab);
	dipoletype::code fluxtype1 = p_nlodipoles->WhichResonant(k1);
	dipoletype::code fluxtype2 = p_nlodipoles->WhichResonant(k2);
  // if(fluxtype==dipoletype::final){
  if(fsr1 && !fsr2){
  	if(!HasFSR()) msg_Error()<<"Wrong dipole type in "<<METHOD<<endl;
  	for (Dipole_Vector::iterator Dip = p_nlodipoles->GetDipoleFF()->begin();
       Dip != p_nlodipoles->GetDipoleFF()->end(); ++Dip) {
  		 double scalek = p_fsr->ScalePhoton(k1);
  		 Dip->SetPhotonScale(scalek);
  		 Dip->AddPhotonToDipole(k1);
  		 if(!Dip->BoostNLO()) {
  		 	msg_Error()<<"NLO boost failed" <<endl;
  		 	return 0;
  		 }
  		 int i(0);
  		 for (auto f : Dip->m_flavs) {
      	p[p_nlodipoles->m_flav_label[f]] =  Dip->GetNewMomenta(i);
      	i++;
    	}
    	k1 = Dip->m_dipolePhotons[0];
  	}
  	p.push_back(k1);
  	MapMomenta(p,k2);
  	p.pop_back();
  }
  if(!fsr1 && fsr2){
  	if(!HasFSR()) msg_Error()<<"Wrong dipole type in "<<METHOD<<endl;
  	Dipole_Vector *diplo = p_dipoles->GetDipoleFF();
  	for (Dipole_Vector::iterator Dip = p_nlodipoles->GetDipoleFF()->begin();
       Dip != p_nlodipoles->GetDipoleFF()->end(); ++Dip) {
  		 double scalek = p_fsr->ScalePhoton(k2);
  		 Dip->SetPhotonScale(scalek);
  		 Dip->AddPhotonToDipole(k2);
  		 if(!Dip->BoostNLO()) {
  		 	msg_Error()<<"NLO boost failed" <<endl;
  		 	return 0;
  		 }
  		 int i(0);
  		 for (auto f : Dip->m_flavs) {
      	p[p_nlodipoles->m_flav_label[f]] =  Dip->GetNewMomenta(i);
      	i++;
    	}
    	k2 = Dip->m_dipolePhotons[0];
  	}
  	p.push_back(k2);
  	MapMomenta(p,k1);
  	p.pop_back();
  }
  if(fsr1 && fsr2){
  	if(!HasFSR()) msg_Error()<<"Wrong dipole type in "<<METHOD<<endl;
  	Dipole_Vector *diplo = p_dipoles->GetDipoleFF();
  	for (Dipole_Vector::iterator Dip = p_nlodipoles->GetDipoleFF()->begin();
       Dip != p_nlodipoles->GetDipoleFF()->end(); ++Dip) {
  		 double scalek = p_fsr->ScalePhoton(k1);
  		 Dip->SetPhotonScale(scalek);
  		 Dip->AddPhotonToDipole(k1);
  		 Dip->AddPhotonToDipole(k2);
  		 if(!Dip->BoostNLO()) {
  		 	msg_Error()<<"NLO boost failed" <<endl;
  		 	return 0;
  		 }
  		 int i(0);
  		 for (auto f : Dip->m_flavs) {
      	p[p_nlodipoles->m_flav_label[f]] =  Dip->GetNewMomenta(i);
      	i++;
    	}
    	k1 = Dip->m_dipolePhotons[0];
    	k2 = Dip->m_dipolePhotons[0];
  	}
  }
  if(!fsr1 && !fsr2) MapMomenta(p, k1, k2);
 	p.push_back(k1);
 	p.push_back(k2);
 	if(fsr1 || fsr2) MapInitial(p);
 	Vec4D_Vector pp = p;
 	pp.pop_back();
 	pp.pop_back();
	p_nlodipoles->MakeDipolesII(m_flavs,pp,m_plab);
	p_nlodipoles->MakeDipolesIF(m_flavs,pp,m_plab);
	p_nlodipoles->MakeDipoles(m_flavs,pp,m_plab);
	double flux;
	flux = p_nlodipoles->CalculateFlux(k1,k2);//*p_nlodipoles->CalculateFlux(k2);
	// else if(m_flux_mode==2) flux = 0.5*(p_nlodipoles->CalculateFlux(kk)+p_nlodipoles->CalculateFlux(k));
	// else flux = p_dipoles->CalculateFlux(kk);
	double tot,rcoll;
	double subloc1 = p_nlodipoles->CalculateRealSub(k1);
	double subloc2 = p_nlodipoles->CalculateRealSub(k2);
	// if(fsrcount) subb = p_dipoles->CalculateRealSubEEX(k);
	// else subb = p_dipoles->CalculateRealSubEEX(kk);
	// if(IsZero(subb)) return 0;
	if(!CheckMomentumConservation(p)) {
		return 0;
	}
	double r = p_realreal->Calc_R(p) / norm /norm;
	if(IsZero(r)) return 0;
	if(IsBad(r) || IsBad(flux)) {
		msg_Error()<<"Bad point for YFS Real"<<std::endl
							 <<"Real ME is : "<<r<<std::endl
							 <<"Flux is : "<<flux<<std::endl;
		return 0;
	}
	m_recola_evts+=1;
	tot = (r -sub1*subloc2*real1 -sub2*subloc1*real2-subloc1*subloc2*m_born*m_rescale_alpha)/sub1/sub2;

  if(IsBad(tot)){
  	msg_Error()<<"NNLO RR is NaN"<<std::endl;
  }
	return tot;
}
	// 		MapMomenta(p, k, kk);
	// 		CheckMasses(p);
	// 		flux  = (Q-kk-k).Abs2()/(Q).Abs2();
	// 		// if(k.PPerp()<0.1 || kk.PPerp()<0.1) return 0;
	// 		if(!CheckPhotonForReal(k)||!CheckPhotonForReal(kk)) return 0;

	// 		Vec4D ksum = k + kk;
	// 		p_nlodipoles->MakeDipoles(m_flavs,p,m_plab);
	// 		p_nlodipoles->MakeDipolesII(m_flavs,p,m_bornMomenta);
	// 		p_nlodipoles->MakeDipolesIF(m_flavs,p,m_plab);
	// 		double subloc1 = p_nlodipoles->CalculateRealSub(k);
	// 		double subloc2 = p_nlodipoles->CalculateRealSub(kk);
	// 		double subb1   = p_dipoles->CalculateRealSub(k);
	// 		double subb2   = p_dipoles->CalculateRealSub(kk);	
	// 		double subb = subloc1*CalculateReal(kk,1);
	// 		subb += subloc2*CalculateReal(k,1);
	// 		subb += subloc1*subloc2*m_born;
	// 		p.push_back(k);
	// 		p.push_back(kk);
	// 		CheckMomentumConservation(p);
	// 		if(m_fsrmode==1){
	// 		for (int i = 2; i < m_plab.size(); ++i)
	// 		{
	// 			Q += m_plab[i];
	// 		}
	// 		double sq = (Q).Abs2();
	// 		double sx = (Q+k+kk).Abs2();
	// 		double sqq = (m_plab[2]+m_plab[3]).Abs2();
	// 		double sxx = (m_plab[2]+m_plab[3]+k+kk).Abs2();	
	// 		 	// subflux = (sq*sq-(MBar*conj(MBar)).real())/(sx*sx-(MBar*conj(MBar)).real())-1;
	// 		// flux = (sqr(sq-mz*mz)+sqr(gz)*sqr(sq)/sqr(mz))/(sqr(sx-mz*mz)+sqr(gz)*sqr(sx)/sqr(mz));
	// 		p_nlodipoles->CalculateFlux(k,kk);
	// 		if(m_isr_debug || m_fsr_debug) m_histograms2d["IFI_EIKONAL"]->Insert(k.Y(),k.PPerp(), p_nlodipoles->CalculateRealSubIF(k));
	// }
	// 		double r = p_realreal->Calc_R(p) / norm / norm *flux;
	// 		if(IsNan(r)) r=0;
	// 		if (r == 0) continue;
	// 		real += (r - subb)/subb1/subb2;
	// 	}
	// }



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
		pRot.Rotate(p[i]);
		boostQ.Boost(p[i]);
		// RandomRotate(p[i]);
	}
	// pRot.Rotate(k);
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
  double lamCM = 0.5*sqrt(SqLam(sqq,m1*m1,m2*m2)/sqq);
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
  double lamCM = 0.5*sqrt(SqLam(sqq,m1*m1,m2*m2)/sqq);
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
	double sign_z = (p[0][3] < 0 ? -1 : 1);
  double m1 = m_flavs[0].Mass();
  double m2 = m_flavs[1].Mass();
  double lamCM = 0.5*sqrt(SqLam(sqq,m1*m1,m2*m2)/sqq);
  double E1 = lamCM*sqrt(1+m1*m1/sqr(lamCM));
  double E2 = lamCM*sqrt(1+m2*m2/sqr(lamCM));
 	p[0] = {E1, 0, 0, sign_z*lamCM};
  p[1] = {E2, 0, 0, -sign_z*lamCM};
	Poincare boostLab(QQ);
  Poincare pRot;
	for (int i = 0; i < 2; ++i)
	{
		if(i==0) pRot = Poincare(p[i], Vec4D(0., 0., 0., 1.));
		// pRot.Rotate(p[i]);
		boostLab.BoostBack(p[i]);
		// pRot2.Rotate(p[i]);
	}
}


void NLO_Base::CheckMasses(Vec4D_Vector &p, int realmode){
	bool allonshell=true;
	std::vector<double> masses;
	Flavour_Vector flavs = m_flavs;
	if(realmode) flavs.push_back(Flavour(kf_photon));
	for (int i = 0; i < p.size(); ++i)
	{
		masses.push_back(flavs[i].Mass());
		if(!IsEqual(p[i].Mass(),flavs[i].Mass(),1e-6)){
			msg_Debugging()<<"Wrong particle masses in YFS Mapping"<<std::endl
								 <<"Flavour = "<<flavs[i]<<", with mass = "<<flavs[i].Mass()<<std::endl
								 <<"Four momentum = "<<p[i]<<", with mass = "<<p[i].Mass()<<std::endl;
			allonshell = false;

		}
	}
	if(!allonshell) m_stretcher.StretchMomenta(p, masses);
	// return true;
}

bool NLO_Base::CheckPhotonForReal(const Vec4D &k) {
	for (int i = 0; i < m_plab.size(); ++i)
	{
		if (m_flavs[i].IsChargedLepton()) {
			double sik = (k + m_plab[i]).Abs2();
			if (sik < m_hardmin ) {
				return false;
			}
		}
	}
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


void NLO_Base::CheckRealSub(Vec4D k){
		// if(k.E() < 20) return;
		// k*=100;
		double real;
		std::string filename="Real_subtracted_";
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
			real=CalculateReal(k);
			out_sub<<k.E()<<","<<fabs(real)<<std::endl;
			if(k.E() < 1e-10 || real==0) break;
			// m_histograms2d["Real_me_sub"]->Insert(k.E(),fabs(real), 1);
		}
		out_sub.close();
		exit(0);
}

void NLO_Base::CheckRealRealSub(Vec4D k1, Vec4D k2){
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
		for (double i = 1; i < 20 ; i+=0.005)
		{
			k1=k1/i;
			real=CalculateRealReal(k1,k2);
			out_sub<<k1.E()<<","<<fabs(real)<<std::endl;
			if(k1.E() < 1e-10 || real==0) break;
			// m_histograms2d["Real_me_sub"]->Insert(k.E(),fabs(real), 1);
		}
		out_sub.close();
		out_sub.open(filename2, std::ios_base::app);
	 	k2 = _k2;
	 	k1 = _k1;
		for (double i = 1; i < 20 ; i+=0.005)
		{
			k2=k2/i;
			real=CalculateRealReal(k1,k2);
			out_sub<<k2.E()<<","<<fabs(real)<<std::endl;
			if(k1.E() < 1e-10 || real==0) break;
			// m_histograms2d["Real_me_sub"]->Insert(k.E(),fabs(real), 1);
		}
		out_sub.close();
		out_sub.open(filename3, std::ios_base::app);
		k2 = _k2;
	 	k1 = _k1;
		for (double i = 1; i < 20 ; i+=0.005)
		{
			k2=k2/i;
			k1=k1/i;
			real=CalculateRealReal(k1,k2);
			out_sub<<k2.E()<<","<<fabs(real)<<std::endl;
			if(k1.E() < 1e-10 || real==0) break;
			// m_histograms2d["Real_me_sub"]->Insert(k.E(),fabs(real), 1);
		}
		out_sub.close();
		exit(0);
}
