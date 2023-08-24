#include "PHASIC++/Channels/Channel_Elements.H"
#include "ATOOLS/Math/Random.H"
#include "YFS/NLO/NLO_Base.H"


using namespace YFS;

std::ofstream out_recola;
std::ofstream out_sub, out_real, out_finite;

NLO_Base::NLO_Base() {
  p_yfsFormFact = new YFS::YFS_Form_Factor();
  p_nlodipoles = new YFS::Define_Dipoles();
  p_global_dipoles = new YFS::Define_Dipoles();
  m_evts = 0;
  m_recola_evts = 0;
  m_realrealtool = 0;
  m_realtool = 0;
  m_looptool = 0;
  m_realvirt = 0;
  if(m_isr_debug || m_fsr_debug){
  	m_histograms2d["Real_me"] = new Histogram_2D(0, -1, 1, 20, 0, sqrt(m_s)/2., 20 );
  	m_histograms1d["Real_diff"] = new Histogram(0, -1, 1, 40);
  }
}


NLO_Base::~NLO_Base() {
	// if (p_real) delete p_real;
	// if (p_virt) delete p_virt;
	// if (p_realvirt) delete p_realvirt;
	// if (p_realreal) delete p_realreal;
	// if (p_griffin)   delete p_griffin;
  if(m_isr_debug || m_fsr_debug){
		Histogram_2D * histo2d;
		string name;
		for (map<string, Histogram_2D *>::iterator hit = m_histograms2d.begin();
		        hit != m_histograms2d.end(); hit++) {
			histo2d = hit->second;
			name  = string("./YFS_NLO_Hist/") + hit->first + string(".dat");
			// histo2d->MPISync();
			histo2d->Finalize();
			histo2d->Output(name);
			delete histo2d;
		}
		Histogram * histo1d;
		for (map<string, Histogram *>::iterator hit = m_histograms1d.begin();
		        hit != m_histograms1d.end(); hit++) {
			histo1d = hit->second;
			name  = string("./YFS_NLO_Hist/") + hit->first + string(".dat");
			histo1d->MPISync();
			histo1d->Finalize();
			histo1d->Output(name);
			delete histo1d;
		}
	}
	PRINT_VAR(m_recola_evts);
	PRINT_VAR(m_evts);
	msg_Out()<<"Percentage of Recola events = "<<m_recola_evts/m_evts*100.<<"% "<<std::endl;
}


void NLO_Base::InitializeVirtual(const PHASIC::Process_Info& pi) {
	if (m_griff != 0) {
		p_griffin = new Griffin::Griffin_Interface();
		p_griffin->Initialize(pi);
		return;
	}
	p_virt = new YFS::Virtual(pi);
	m_looptool = true;
}

void NLO_Base::InitializeReal(const PHASIC::Process_Info& pi) {
	p_real = new YFS::Real(pi);
	m_realtool = true;
}

void NLO_Base::InitializeRealVirtual(const PHASIC::Process_Info& pi) {
	p_realvirt = new YFS::RealVirtual(pi);
	m_realvirt = true;
}

void NLO_Base::InitializeRealReal(const PHASIC::Process_Info& pi) {
	p_realreal = new YFS::RealReal(pi);
	m_realrealtool = true;
}


void NLO_Base::Init(Flavour_Vector &flavs, Vec4D_Vector &plab, Vec4D_Vector &born) {
	m_flavs = flavs;
	m_plab = plab;
	m_bornMomenta = born;
}


double NLO_Base::CalculateNLO() {
	double result{0.0};
	result += CalculateVirtual();
	result += CalculateReal();
	result += CalculateRealVirtual();
	result += CalculateRealReal();
	return result;
}


double NLO_Base::CalculateVirtual() {
	if (m_griff != 0) {
		double v =   p_griffin->EvaluateLoop(m_plab);
		// PRINT_VAR(v);
		return v;
	}
	if (!m_looptool) return 0;
	double virt;
	double sub;
	if (m_check_mass_reg) {
		out_sub.open("yfs-sub.txt", std::ios_base::app);
		out_recola.open("recola-res.txt", std::ios_base::app); // append instead of overwrite
		out_finite.open("yfs-finite.txt", std::ios_base::app);
		virt = p_virt->Calc(m_plab, m_born);
		if (!IsEqual(m_born, p_virt->p_loop_me->ME_Born(), 1e-4)) {
			msg_Error() << METHOD << "\n Warning! Loop provider's born is different! YFS Subtraction likely fails\n"
									<< "Loop Provider " << ":  "<<p_virt->p_loop_me->ME_Born()
									<< "Sherpa" << ":  "<<m_born;
		}
		// born = p_virt->p_loop_me->ME_Born();
		sub = p_dipoles->CalculateVirtualSub();
		std::cout << setprecision(10);
		out_sub<< setprecision(10) << m_photonMass << "," << -sub*m_born << std::endl;
		out_recola<< setprecision(10) << m_photonMass << "," << virt << std::endl;
		out_finite<< setprecision(10) << m_photonMass << "," << virt  -sub*m_born << std::endl;
		out_sub.close();
		out_recola.close();
		exit(1);
	}
	if(m_fsrmode==2) virt = p_virt->Calc(m_bornMomenta, m_born);
	else virt = p_virt->Calc(m_plab, m_born);
	if (!IsEqual(m_born, p_virt->p_loop_me->ME_Born(), 1e-4)) {
		msg_Error() << METHOD << "\n Warning! Loop provider's born is different! YFS Subtraction likely fails\n"
								<< "Loop Provider " << ":  "<<p_virt->p_loop_me->ME_Born()
								<< "\nSherpa" << ":  "<<m_born<<std::endl
								<<"PhaseSpace Point = ";
		for(auto p: m_plab) msg_Error()<<p<<std::endl;
	}
	sub = p_dipoles->CalculateVirtualSub();
	m_oneloop = (virt - sub * m_born);
	return m_oneloop;
}


double NLO_Base::CalculateReal() {
	if (!m_realtool) return 0;
	double real(0), sub(0);
	m_is_isr=1;
	for (auto k : m_ISRPhotons) real+= CalculateReal(k);
	m_is_isr = 0;
	for (auto k : m_FSRPhotons) real+= CalculateReal(k);
	// if(m_FSRPhotons.size()!=0) real+=CalculateReal(m_FSRPhotons[0]);
 	return real;
}


double NLO_Base::CalculateReal(Vec4D k) {
	double norm = 2.*pow(2 * M_PI, 3);
	m_sp = (m_plab[0]+m_plab[1]).Abs2();
	Vec4D_Vector p(m_plab),pi(m_bornMomenta), pf(m_bornMomenta);
	Vec4D kk = k;
	double B = m_born;
	MapMomenta(p, k);
	CheckMasses(p);
	double sx = (m_plab[0]+m_plab[1]).Abs2();
	double sq = (m_eikmom[3]+m_eikmom[2]+kk).Abs2();
	p_nlodipoles->MakeDipoles(m_flavs,p,m_plab);
	p_nlodipoles->MakeDipolesII(m_flavs,p,p);
	p_nlodipoles->MakeDipolesIF(m_flavs,p,m_plab);
	p_global_dipoles->MakeDipolesII(m_flavs,p,m_bornMomenta);
	p_global_dipoles->MakeDipolesIF(m_flavs,p,m_bornMomenta);
	p_global_dipoles->MakeDipoles(m_flavs,m_plab,m_bornMomenta);
	double flux;
	Vec4D Q = m_bornMomenta[0]+m_bornMomenta[1];
	flux = (Q-kk).Abs2()/(m_bornMomenta[0]+m_bornMomenta[1]).Abs2();
	double tot,colltot,rcoll;
	double subloc = p_nlodipoles->CalculateRealSub(k);
	double subb   = p_dipoles->CalculateRealSubEEX(kk);
	m_evts+=1;
	rcoll = p_dipoles->CalculateEEXReal(kk)*B;
	if (!CheckPhotonForReal(kk)) { 
		if(m_no_subtraction) return rcoll/subb;
		return ( rcoll/subb - B);
	}
	double eex = rcoll/subb - B;
	if(m_fsrmode==2) flux = sqr(sx/sq)*(sqr(sq-91.1876*91.1876)+sqr(2.4952*sq)/sqr(91.1876))/(sqr(sx-91.1876*91.1876)+sqr(2.4952*sx)/sqr(91.1876));
	p.push_back(k);
	double r = p_real->Calc_R(p) / norm * flux; 
	if(IsBad(r) || IsBad(flux)) {
		msg_Error()<<"Bad point for YFS Real"<<std::endl
							 <<"Real ME is : "<<r<<std::endl
							 <<"Flux is : "<<flux<<std::endl;
		return 0;
	}
	m_recola_evts+=1;
	tot =  (r-subloc*m_born)/subloc;
  if(m_isr_debug || m_fsr_debug){
		m_histograms2d["Real_me"]->Insert(k.CosTheta(m_bornMomenta[0]), k.E(), r/rcoll);
		double diff = (tot-colltot)/(tot+colltot);
		m_histograms1d["Real_diff"]->Insert(diff);
  }
  if(m_no_subtraction) return r/subloc;
	return tot; 
}

double NLO_Base::CalculateRealVirtual() {
	if (!m_realvirt) return 0;
	double real(0), sub(0);
	double norm = 2.*pow(2 * M_PI, 3);
	Vec4D_Vector photons;
	for (auto k : m_ISRPhotons) photons.push_back(k);
	for (auto k : m_FSRPhotons) photons.push_back(k);
	for (auto k : photons) {
		Vec4D_Vector p(m_plab);
		MapMomenta(p, k);
		if(!CheckPhotonForReal(k)) return 0;
		if(k.E()/sqrt(m_s) < 1e-3 || k.PPerp() < 0.1 ) return 0;

		p_nlodipoles->MakeDipoles(m_flavs,p,m_bornMomenta);
		p_nlodipoles->MakeDipolesII(m_flavs,p,m_bornMomenta);
		p_nlodipoles->MakeDipolesIF(m_flavs,p,m_bornMomenta);
		double subb   = p_dipoles->CalculateRealSub(k);
		double subloc = p_nlodipoles->CalculateRealSub(k);
		p.push_back(k);
		double r = p_realvirt->Calc(p, m_born);
		if (r == 0 || IsBad(r)) continue;
		double aB = p_dipoles->CalculateRealSub(k)*m_oneloop;
		// PRINT_VAR(aB);
		// PRINT_VAR(r);
		double tot = (r-aB) / p_dipoles->CalculateRealSub(k);
		// tot -= m_born * p_yfsFormFact->BVV_full(p[0], p[1], m_photonMass, sqrt(m_s) / 2., 3);
		real += tot;
	}
	return real/norm;
}


double NLO_Base::CalculateRealReal() {
	if (!m_realrealtool) return 0;
	// if (m_ISRPhotons.size() <= 1 && m_FSRPhotons.size() <= 1) return 0;
	double real(0), sub(0);
	Vec4D_Vector photons, p(m_plab);
	for (auto k : m_ISRPhotons) photons.push_back(k);
	// for (auto k : m_FSRPhotons) photons.push_back(k);
	// if (NHardPhotons(photons) == 0) return 0;
	double norm = 2.*pow(2 * M_PI, 3);
	for (int i = 0; i < photons.size(); ++i) {
		for (int j = 0; j < i; ++j) {
			Vec4D k  = photons[i];
			Vec4D kk = photons[j];
			// if(k.PPerp()<0.1 || kk.PPerp()<0.1) return 0;
			Vec4D ksum = k + kk;
			p = m_plab;
			MapMomenta(p, ksum);
			p_nlodipoles->MakeDipoles(m_flavs,p,m_bornMomenta);
			p_nlodipoles->MakeDipolesII(m_flavs,p,p);
			p_nlodipoles->MakeDipolesIF(m_flavs,p,m_bornMomenta);
			double subloc1 = p_nlodipoles->CalculateRealSub(k);
			double subloc2 = p_nlodipoles->CalculateRealSub(kk);
			double subb1   = p_dipoles->CalculateRealSub(k);
			double subb2   = p_dipoles->CalculateRealSub(kk);	
			if(k.E()/sqrt(m_s) < 1e-3 || kk.E()/sqrt(m_s) < 1e-3 ) return 0;
			if(k.PPerp() < 1e-1 || kk.PPerp() < 1e-1 ) return 0;
			// if(!CheckPhotonForReal(k)||!CheckPhotonForReal(kk)) return 0;
			double subb = subloc1*CalculateReal(kk)*subb2;
			subb += subloc2*CalculateReal(k)*subb1;
			subb += subloc1*subloc2*m_born;
			p.push_back(k);
			p.push_back(kk);
			double r = p_realreal->Calc_R(p) / norm / norm;
			if(IsNan(r)) r=0;
			if (r == 0) continue;
			real += (r - subb)/subb1/subb2;
		}
	}
	return real;
}


void NLO_Base::MapInitial(Vec4D_Vector &p, Vec4D &k){
	Vec4D Q;
  Q = p[0] + p[1] - k;
  // if(Q.Abs2() > )
  double sp = Q * Q;
  double zz = sqrt(sp) / 2.;
  double z = zz * sqrt((sp - sqr(p[0].Mass() - p[1].Mass())) * (sp - sqr(p[0].Mass() + p[1].Mass()))) / sp;
  p[0] = {zz, 0, 0, z};
  p[1] = {zz, 0, 0, -z};
  ATOOLS::Poincare poin(Q);
  Poincare pRot(p[0], Vec4D(0., 0., 0., 1.));
  for (int i = 0; i < 2; ++i) {
    pRot.Rotate(p[i]);
    poin.BoostBack(p[i]);
  }
	PHASIC::CE.Isotropic2Momenta(p[0]+p[1]-k, sqr(p[2].Mass()), sqr(p[3].Mass()), p[2], p[3], ran->Get(), ran->Get(), -1, 1);

  PRINT_VAR(p[0]+p[1]-k);
  PRINT_VAR(p[3]+p[2]);
}

void NLO_Base::MapMomenta(Vec4D_Vector &p, Vec4D &k, int mode) {
	Vec4D Q;
	Vec4D QQ, PP;
	Poincare boostLab(m_bornMomenta[0] + m_bornMomenta[1]);
	// Poincare boostLab(p[0] + p[1]);
	if(mode==0){
		for (int i = 2; i < p.size(); ++i)
		{
			Q += p[i];
		}
		Q += k;
	}
	else if (mode==1){
		for (int i = 0; i < 2; ++i)
		{
			Q += p[i];
		}
		Q -= k;	
	}
	// else Q-=k;

	double sq = Q.Abs2();
	// else Q-=k;
	Poincare boostQ(Q);
  Poincare pRot(m_bornMomenta[0], Vec4D(0., 0., 0., 1.));
  if(mode==0){
		for (int i = 2; i < p.size(); ++i) {
			pRot.Rotate(p[i]);
			boostQ.Boost(p[i]);
		}
	}
	else{
		for (int i = 0; i < 2; ++i) {
			pRot.Rotate(p[i]);
			boostQ.Boost(p[i]);
		}
	}
	pRot.Rotate(k);
	boostQ.Boost(k);
	double qx(0), qy(0), qz(0);
	if(mode==0){
		for (int i = 2; i < p.size(); ++i)
		{
			qx += p[i][1];
			qy += p[i][2];
			qz += p[i][3];
		}
		// if (!IsEqual(k[1], -qx, 1e-6) || !IsEqual(k[2], -qy, 1e-6) || !IsEqual(k[3], -qz, 1e-6) ) {
		// 	if( k[1]> 1e-8 && k[2]> 1e-8 && k[3]> 1e-8 ){
		// 		msg_Error() << "YFS Mapping has failed for ISR\n";
		// 		msg_Error() << " Photons px = " << k[1] << "\n Qx = " << -qx << std::endl;
		// 		msg_Error() << " Photons py = " << k[2] << "\n Qy = " << -qy << std::endl;
		// 		msg_Error() << " Photons pz = " << k[3] << "\n Qz = " << -qz << std::endl;
		// 	}
		// }
		for (int i = 2; i < p.size(); ++i)
		{
			QQ += p[i];
		}
		QQ+=k;
		double sqq = QQ.Abs2();
		if (!IsEqual(sqq, sq, 1e-8))
		{
			msg_Error() << "YFS Real mapping not conserving momentum in " << METHOD << std::endl;
		}
		// if(m_is_isr) QQ = p[0]+p[1];
    double zz = sqrt(sqq) / 2.;
		double z = zz * sqrt((sqq - sqr(m_flavs[0].Mass() - m_flavs[1].Mass())) * (sqq - sqr(m_flavs[0].Mass() + m_flavs[1].Mass()))) / sqq;
		p[0] = {zz, 0, 0, z};
		p[1] = {zz, 0, 0, -z};
		// PHASIC::CE.Isotropic2Momenta(QQ, sqr(m_flavs[0].Mass()), sqr(m_flavs[1].Mass()), p[0], p[1], ran->Get(), ran->Get(), -1, 1);
	}
	else {
		for (int i = 0; i < 2; ++i)
		{
			QQ += p[i];
		}
		QQ-=k;
		double sqq = QQ.Abs2();
		if (!IsEqual(sqq, sq, 1e-8))
		{
			msg_Error() << "YFS Real mapping not conserving momentum in " << METHOD << std::endl;
		}
		// if(m_is_isr) QQ = p[0]+p[1];
		// double zz = sqrt(sqq) / 2.;
		// double z = zz * sqrt((sqq - sqr(m_flavs[0].Mass() - m_flavs[1].Mass())) * (sqq - sqr(m_flavs[0].Mass() + m_flavs[1].Mass()))) / sqq;
		// p[0] = {zz, 0, 0, z};
		// p[1] = {zz, 0, 0, -z};
		PHASIC::CE.Isotropic2Momenta(QQ, sqr(m_flavs[2].Mass()), sqr(m_flavs[3].Mass()), p[2], p[3], ran->Get(), ran->Get(), -1, 1);
	}
  Poincare pRot2(p[0], Vec4D(0., 	0., 0., 1.));
	for (int i = 0; i < p.size(); ++i)
	{
		// pRot2.Rotate(p[i]);
		boostLab.Boost(p[i]);
	}
	// pRot2.Rotate(k);
	boostLab.Boost(k);
}

void NLO_Base::CheckMasses(Vec4D_Vector &p){
	bool allonshell=true;
	std::vector<double> masses;
	for (int i = 0; i < p.size(); ++i)
	{
		masses.push_back(m_flavs[i].Mass());
		if(!IsEqual(p[i].Mass(),m_flavs[i].Mass(),1e-6)){
			msg_Debugging()<<"Wrong particle masses in YFS Mapping"<<std::endl
								 <<"Flavour = "<<m_flavs[i]<<", with mass = "<<m_flavs[i].Mass()<<std::endl
								 <<"Four momentum = "<<p[i]<<", with mass = "<<p[i].Mass()<<std::endl;
			allonshell = false;

		}
	}
	if(!allonshell) m_stretcher.StretchMomenta(p, masses);
	// return true;
}

bool NLO_Base::CheckPhotonForReal(const Vec4D &k) {
	// if (k.E() < m_hardmin || 2.*k.E()/sqrt(m_s) <  1e-3 ) return false;
	if (k.E() < m_hardmin ) return false;
	// if(k.PPerp() < 0.1) return false;
	// for (int i = 0; i < m_bornMomenta.size(); ++i)
	// {
	// 	if (m_flavs[i].IsChargedLepton()) {
	// 		double sik = (k + m_plab[i]).Abs2();
	// 		if (sqrt(sik)/sqrt(m_s) < m_hardmin ) {
	// 			return false;
	// 		}
	// 	}
	// }
	return true;
}


void NLO_Base::MakeHardMomenta(){
	 // m_plab[0] = 
	m_eikmom = m_plab;
	 for (Dipole_Vector::iterator Dip = p_dipoles->GetDipoleFF()->begin();
         Dip != p_dipoles->GetDipoleFF()->end(); ++Dip) {
	 			for (int i = 0; i < Dip->m_eikmomentum.size(); ++i)
	 			{
	 				m_eikmom[i+2] = Dip->GetEikMomenta(i);
	 			}
	 }
}