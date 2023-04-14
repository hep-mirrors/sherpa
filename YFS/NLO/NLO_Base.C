#include "PHASIC++/Channels/Channel_Elements.H"
#include "ATOOLS/Math/Random.H"
#include "YFS/NLO/NLO_Base.H"




using namespace YFS;

std::ofstream out_recola;
std::ofstream out_sub, out_real, out_finite;

NLO_Base::NLO_Base() {
  p_yfsFormFact = new YFS::YFS_Form_Factor();

}


NLO_Base::~NLO_Base() {
	if (p_real) delete p_real;
	if (p_virt) delete p_virt;
	if (p_realvirt) delete p_realvirt;
	if (p_realreal) delete p_realreal;
	if (p_griffin)   delete p_griffin;
}


void NLO_Base::InitializeVirtual(const PHASIC::Process_Info& pi) {
	if (m_griff != 0) {
		p_griffin = new Griffin::Griffin_Interface();
		p_griffin->Initialize(pi);
		return;
	}
	p_virt = new YFS::Virtual(pi);
	m_looptool = true;
	m_real_only = true;
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
	double result;
	if (p_virt || m_griff) result += CalculateVirtual();
	if (p_real) result += CalculateReal();
	if (p_realvirt) result += CalculateRealVirtual();
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
		std::cout << setprecision(16);
		out_sub << m_photonMass << "," << -sub*m_born << std::endl;
		out_recola << m_photonMass << "," << virt << std::endl;
		out_finite << m_photonMass << "," << virt  -sub*m_born << std::endl;
		out_sub.close();
		out_recola.close();
		exit(1);
	}
	// if (m_check_poles) {
	// 	out_sub.open("yfs-sub.txt", std::ios_base::app);
	// 	virt = p_virt->Calc(m_plab, m_born);
	// 	DivArrC  subd = p_yfsFormFact->BVV_full_eps(m_plab[0], m_plab[1], m_photonMass, sqrt(m_s) / 2., 3);
	// 	PRINT_VAR(subd.GetEpsilon());
	// 	out_sub << rpa->gen.Ecms() << "," << -sub*m_born << std::endl;
	// 	out_sub.close();
	// 	exit(1);
	// }
	virt = p_virt->Calc(m_plab, m_born);
	if (!IsEqual(m_born, p_virt->p_loop_me->ME_Born(), 1e-4)) {
		msg_Error() << METHOD << "\n Warning! Loop provider's born is different! YFS Subtraction likely fails\n"
								<< "Loop Provider " << ":  "<<p_virt->p_loop_me->ME_Born()
								<< "\nSherpa" << ":  "<<m_born<<std::endl;
	}
	sub = p_dipoles->CalculateVirtualSub();
	m_oneloop = (virt - sub * m_born);
	return m_oneloop;
}


double NLO_Base::CalculateReal() {
	if (!m_realtool) return 0;
	double real(0), sub(0);
	Vec4D_Vector photons;
	for (auto k : m_ISRPhotons) photons.push_back(k);
	for (auto kk : m_FSRPhotons) photons.push_back(kk);
	int i(-1);
	int NISR = m_ISRPhotons.size();
	int nn(0);
	for (auto k : photons) {
		i++;
		if (!CheckPhotonForReal(k)) {
			// use collinear approx
			// if (nn < NISR) {
			// 	p_realff->SetMode(0);
			// 	p_realff->SetIncoming(p_dipoles->GetDipoleII());
			// 	p_realff->SetBorn(m_born);
			// 	real += p_realff->Beta10(k);
			// }
			// else if (m_fsrmode != 0 && nn > NISR) {
			// 	p_realff->SetMode(m_fsrmode);
			// 	for (Dipole_Vector::iterator Dip = p_dipoles->GetDipoleFF()->begin();
			// 	        Dip != p_dipoles->GetDipoleFF()->end(); ++Dip) {
			// 		p_realff->SetIncoming(Dip, m_bornMomenta, m_fsrphotonsforME, p_fsr->m_yini, p_fsr->m_zini);
			// 		p_realff->SetBorn(m_born);
			// 		real += p_realff->Beta10(k);
			// 	}
			// }
			continue;
		}
		double tot = CalculateReal(k);
		real += tot;
	}
	return real;
}


double NLO_Base::CalculateReal(Vec4D &k) {
	double norm = 2.*pow(2 * M_PI, 3);
	Vec4D_Vector p(m_plab);
	Vec4D kk = k;
	MapMomenta(p, k);
	p_dipoles->MakeDipolesIF(m_flavs, m_plab, p);
	double subb   = p_dipoles->CalculateRealSub(k);
	double subloc = p_dipoles->CalculateRealSubLocal(p, k);
	p.push_back(k);
	double r = p_real->Calc_R(p) / norm;
	if (r == 0) return 0;
	double tot =  ( r - subb * m_born) / subb;
	return tot;
}


double NLO_Base::CalculateRealVirtual() {
	if (!m_realvirt || m_ISRPhotons.size() == 0) return 0;
	double real(0), sub(0);
	double norm = 2.*pow(2 * M_PI, 3);
	Vec4D_Vector photons;
	for (auto k : m_ISRPhotons) photons.push_back(k);
	for (auto k : m_FSRPhotons) photons.push_back(k);
	for (auto k : photons) {
		Vec4D_Vector p(m_plab);
		MapMomenta(p, k);
		double subb  = p_dipoles->CalculateRealSub(k);
		double subloc = 0;
		p.push_back(k);
		double r = p_realvirt->Calc(p, m_born) / norm;
		if (r == 0 || IsBad(r)) continue;
		double aB = p_dipoles->CalculateRealVirtualSub(k)*CalculateVirtual();
		// PRINT_VAR(aB);
		// PRINT_VAR(r);
		double tot = (r-aB) / p_dipoles->CalculateRealSub(k);
		// double tot =  ( r + p_dipoles->CalculateRealSubLocal(p,k) * (CalculateVirtual())) / p_dipoles->CalculateRealSub(k);
		// tot -= m_born * p_yfsFormFact->BVV_full(p[0], p[1], m_photonMass, sqrt(m_s) / 2., 3);
		real += tot;
	}
	return real;
}


double NLO_Base::CalculateRealReal() {
	if (!m_realrealtool) return 0;
	// if (m_ISRPhotons.size() <= 1 && m_FSRPhotons.size() <= 1) return 0;
	double real(0), sub(0);
	Vec4D_Vector photons;
	for (auto k : m_ISRPhotons) photons.push_back(k);
	for (auto k : m_FSRPhotons) photons.push_back(k);
	// if (NHardPhotons(photons) == 0) return 0;
	double norm = 2.*pow(2 * M_PI, 3);
	// double real =
	// p_realff->CalculateVirt();
	// int NISR = m_ISRPhotons.size();
	int nisr(0);
	int nfsr(m_ISRPhotons.size());
	// for (auto &k : photons) {
	// for (auto &kk: photons){
	for (int i = 0; i < photons.size(); ++i) {
		for (int j = 0; j < i; ++j) {
			Vec4D k  = photons[i];
			Vec4D kk = photons[j];

			// k *= scale;
			if (!CheckPhotonForReal(k) || !CheckPhotonForReal(kk)) {
				// use collinear approx
				// if (nn < NISR) {
				//   p_realff->SetMode(0);
				//   p_realff->SetIncoming(p_dipoles->GetDipoleII());
				//   p_realff->SetBorn(m_born);
				//   real += p_realff->Beta10(k);
				// }
				// else if (m_fsrmode != 0 && nn > NISR) {
				//   p_realff->SetMode(m_fsrmode);
				//   for (Dipole_Vector::iterator Dip = p_dipoles->GetDipoleFF()->begin();
				//        Dip != p_dipoles->GetDipoleFF()->end(); ++Dip) {
				//     p_realff->SetIncoming(Dip, m_bornMomenta, m_fsrphotonsforME, p_fsr->m_yini, p_fsr->m_zini);
				//     p_realff->SetBorn(m_born);
				//     // p_realff->CalculateVirt();
				//     // PRINT_VAR(p_realff->Beta10(k));
				//     real += p_realff->Beta10(k);
				//   }
				// }
				// return real;
				continue;
			}
			// if(k.PPerp() < m_hardmin) continue;
			Vec4D_Vector p(m_reallab);
			Vec4D Q;
			Poincare boostLab(m_bornMomenta[0] + m_bornMomenta[1]);
			for (int i = 2; i < p.size(); ++i)
			{
				Q += p[i];
			}
			double sq = Q.Abs2();
			Q += k;
			Q += kk;
			Poincare boostQ(Q);
			for (int i = 0; i < p.size(); ++i) {
				boostQ.Boost(p[i]);
			}
			boostQ.Boost(k);
			boostQ.Boost(kk);
			double qx(0), qy(0), qz(0);
			for (int i = 2; i < p.size(); ++i)
			{
				qx += p[i][1];
				qy += p[i][2];
				qz += p[i][3];
			}
			Vec4D tk = k + kk;
			if (!IsEqual(tk[1], -qx, 1e-3) || !IsEqual(tk[2], -qy, 1e-3) || !IsEqual(tk[3], -qz, 1e-3) ) {
				// if( k[1]!=-qx || k[2]!=-qy || k[3]!=-qz ){
				msg_Error() << "YFS Mapping has failed for ISR\n";
				msg_Error() << " Photons px = " << tk[1] << "\n Qx = " << -qx << std::endl;
				msg_Error() << " Photons py = " << tk[2] << "\n Qy = " << -qy << std::endl;
				msg_Error() << " Photons pz = " << tk[3] << "\n Qz = " << -qz << std::endl;
				continue;
			}
			Vec4D QQ, PP;
			for (int i = 2; i < p.size(); ++i)
			{
				QQ += p[i];
			}
			double sqq = QQ.Abs2();
			if (!IsEqual(sqq, sq, 1e-4))
			{
				msg_Error() << "YFS Real mapping not conserving momentum in " << METHOD << std::endl;
			}
			QQ += k;
			QQ += kk;
			PHASIC::CE.Isotropic2Momenta(QQ, sqr(p[0].Mass()), sqr(p[1].Mass()), p[0], p[1], ran->Get(), ran->Get(), -1, 1);
			for (int ii = 0; ii < p.size(); ++ii)
			{
				boostLab.Boost(p[ii]);
			}
			for (int ij = 0; ij < 2; ++ij)
			{
				PP += p[ij];
			}
			double s1 = p_dipoles->CalculateRealSubLocal(p, k);
			double s2 = p_dipoles->CalculateRealSubLocal(p, kk);
			double subb = p_dipoles->CalculateRealSub(k) * p_dipoles->CalculateRealSub(kk);
			// subloc += Eikonal(k, p[0], p[1]) * m_realcorr[j];
			// subloc += Eikonal(kk, p[0], p[1]) * m_realcorr[i];
			// subloc += Eikonal(k, p[0], p[1]) * Eikonal(kk, p[0], p[1]) * m_born;
			// if (m_fsrmode != 0) {
			//   subloc += Eikonal(k, p[2], p[3]);
			//   subloc += Eikonal(k, p[0], p[2]);
			//   subloc -= Eikonal(k, p[0], p[3]);
			//   subloc -= Eikonal(k, p[1], p[2]);
			//   subloc += Eikonal(k, p[1], p[3]);
			// }
			// for (int i = 0; i < p.size; ++i)
			// {
			//   subloc +=
			// }
			p.push_back(k);
			p.push_back(kk);
			double r = p_realreal->Calc_R(p) / norm / norm;
			// double rc = p_realff->Beta10(k);
			if (r == 0) continue;
			// double tot =  ( r - s1 * r2 - s2 * r1 - s1 * s2 * m_born) / subb;
			// real += (r/Eikonal(k,p[0], p[1]) -  m_born)*m_sp/m_s;
			// real += tot;
		}
	}
	PRINT_VAR(real);
	return real;
}




void NLO_Base::MapMomenta(Vec4D_Vector &p, Vec4D &k) {
	Vec4D Q;
	Poincare boostLab(m_bornMomenta[0] + m_bornMomenta[1]);
	for (int i = 2; i < p.size(); ++i)
	{
		Q += p[i];
	}
	double sq = Q.Abs2();
	Q += k;
	Poincare boostQ(Q);
	for (int i = 2; i < p.size(); ++i) {
		boostQ.Boost(p[i]);
	}
	boostQ.Boost(k);
	double qx(0), qy(0), qz(0);
	for (int i = 2; i < p.size(); ++i)
	{
		qx += p[i][1];
		qy += p[i][2];
		qz += p[i][3];
	}
	if (!IsEqual(k[1], -qx, 1e-3) || !IsEqual(k[2], -qy, 1e-3) || !IsEqual(k[3], -qz, 1e-3) ) {
		// if( k[1]!=-qx || k[2]!=-qy || k[3]!=-qz ){
		msg_Error() << "YFS Mapping has failed for ISR\n";
		msg_Error() << " Photons px = " << k[1] << "\n Qx = " << -qx << std::endl;
		msg_Error() << " Photons py = " << k[2] << "\n Qy = " << -qy << std::endl;
		msg_Error() << " Photons pz = " << k[3] << "\n Qz = " << -qz << std::endl;
	}
	Vec4D QQ, PP;
	for (int i = 2; i < p.size(); ++i)
	{
		QQ += p[i];
	}
	double sqq = QQ.Abs2();
	if (!IsEqual(sqq, sq, 1e-4))
	{
		msg_Error() << "YFS Real mapping not conserving momentum in " << METHOD << std::endl;
	}
	QQ += k;
	PHASIC::CE.Isotropic2Momenta(QQ, sqr(p[0].Mass()), sqr(p[1].Mass()), p[0], p[1], ran->Get(), ran->Get(), -1, 1);
	for (int i = 0; i < p.size(); ++i)
	{
		boostLab.Boost(p[i]);
	}
	for (int i = 0; i < 2; ++i)
	{
		PP += p[i];
	}
}

bool NLO_Base::CheckPhotonForReal(const Vec4D &k) {
	for (int i = 0; i < m_bornMomenta.size(); ++i)
	{
		if (m_flavs[i].IsChargedLepton()) {
			double sik = (k + m_bornMomenta[i]).Abs2();
			if (sqrt(sik) < m_hardmin) {
				return false;
			}
			if (k.PPerp()<0.1) return false;
		}
	}
	return true;
}