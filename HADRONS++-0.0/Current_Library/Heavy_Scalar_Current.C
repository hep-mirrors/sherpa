#include "Heavy_Scalar_Current.H"

using namespace HADRONS;
using namespace ATOOLS;

void Heavy_Scalar_Current::SetModelParameters( struct GeneralModel _md )
{
  m_Vxx = _md("Vxx",0.04);
  switch( int(_md("HS_FORM_FACTOR", 4)+0.5) ) {
  case 1:
    p_ff = new ISGW( _md, p_masses );
    msg.Tracking()<<"    Using ISGW form factor model for "<<m_name<<std::endl;
    break;
  case 2:
    p_ff = new ISGW2( _md, p_masses );
    msg.Tracking()<<"    Using ISGW2 form factor model for "<<m_name<<std::endl;
    break;
  case 3:
    p_ff = new HQET( _md, p_masses );
    msg.Tracking()<<"    Using HQET form factor model for "<<m_name<<std::endl;
    break;
  case 4:
    p_ff = new HQET2( _md, p_masses );
    msg.Tracking()<<"    Using HQET2 form factor model for "<<m_name<<std::endl;
    break;
  default:
    msg.Error()<<METHOD<<": You chose a form factor model which does not "
      <<"exist for current "<<m_name<<". Aborting."<<std::endl;
    abort();
  }
}

void Heavy_Scalar_Current::Calc()
{
  Vec4D pB = p_moms[0];
  Vec4D pS = p_moms[1];
  Vec4D q = pB - pS;
  double q2 = q*q;
  double mB = p_masses[0];
  double mS = p_masses[1];

  p_ff->CalcFFs(pB,pS);
  double Fplus = p_ff->Fplus();
  double F0 = p_ff->F0();

  p_results[0] = m_Vxx * (Fplus * (pB+pS - (mB-mS)/q2*(pB-pS)) + F0 * (mB-mS)/q2*(pB-pS));
}


Heavy_Scalar_Current::HQET::HQET( GeneralModel _md, double* _masses )
  : FF_Base( _md, _masses ) 
{
  m_rho2  = _md("HQET_rho2",0.7);
  m_c     = _md("HQET_R1",0.0);
}

void Heavy_Scalar_Current::HQET::CalcFFs( Vec4D pB, Vec4D pS )
{
  Vec4D v1 = pB/m_mB;
  Vec4D v2 = pS/m_mS;
  double w = v1*v2;
  m_Fplus  = 1-m_rho2*(w-1)+m_c*(w-1)*(w-1);
  m_F0     = 0.0;
  m_calced = true;
}




Heavy_Scalar_Current::HQET2::HQET2( GeneralModel _md, double* _masses )
  : FF_Base( _md, _masses ) 
{
  m_rho2  = _md("HQET2_rho2",0.7);
  m_v1_1  = _md("HQET2_v1_1",0.98);
}

void Heavy_Scalar_Current::HQET2::CalcFFs( Vec4D pB, Vec4D pS )
{
  double w = (pB/m_mB)*(pS/m_mS);
  const double z = (sqrt(w+1)-sqrt(2.0))/(sqrt(w+1)+sqrt(2.0));
  double v1 = m_v1_1*(1.0- 8.0*m_rho2*z + (51.0*m_rho2-10.0)*z*z - (252.0*m_rho2-84.0)*z*z*z);
  m_Fplus  = v1;
  m_F0     = 0.0;
  m_calced = true;
}





Heavy_Scalar_Current::ISGW::ISGW( GeneralModel _md, double* _masses )
  : FF_Base( _md, _masses ) 
{
  // Defaults are set for B -> D
  m_msb    = _md("ISGW_msb",5.2);
  m_msd    = _md("ISGW_msd",0.33);
  m_bb2    = _md("ISGW_bb2",0.41*0.41);

  m_msq    = _md("ISGW_msq",1.82);
  m_bx2    = _md("ISGW_bx2",0.31*0.31);
  m_kap2   = _md("ISGW_kap2",0.7*0.7);

}

void Heavy_Scalar_Current::ISGW::CalcFFs( Vec4D pB, Vec4D pS )
{
  Vec4D q = pB - pS;
  double q2 = q*q;
 
  double mtb = m_msb +m_msd;
  double mtx = m_msq + m_msd;
  double mup=1.0/(1.0/m_msq+1.0/m_msb);
  double mum=1.0/(1.0/m_msq-1.0/m_msb);
  double bbx2=0.5*(m_bb2+m_bx2);
  double tm=(m_mB-m_mS)*(m_mB-m_mS);
  if ( q2>tm ) { q2=0.99*tm; }
  double F3 = sqrt(mtx/mtb)*pow(sqrt(m_bx2*m_bb2)/bbx2,1.5)*
              exp(-1.0*((m_msd*m_msd*(tm-q2)/(4.0*mtb*mtx*m_kap2*bbx2))));
  m_Fplus = F3*(1+(m_msb/(2.0*mum))-(m_msb*m_msq*m_msd*m_bb2/(4.0*mup*mum*mtx*bbx2)));
  double Fminus = F3*(1.0-(mtb+mtx)*(0.5/m_msq-(m_msd*m_bb2/(4.0*mup*mtx*bbx2))));
  m_F0 = Fminus/((m_mB*m_mB-m_mS*m_mS)/q2)+m_Fplus;
  m_calced = true;
}



Heavy_Scalar_Current::ISGW2::ISGW2( GeneralModel _md, double* _masses )
  : FF_Base( _md, _masses ) 
{
  // Defaults are set for B -> D
  m_msb = _md("ISGW2_msb",5.2);
  m_msd = _md("ISGW2_msd",0.33);
  m_bb2 = _md("ISGW2_bb2",0.185761);
  m_mbb = _md("ISGW2_mbb",5.31);
  m_nf  = _md("ISGW2_nf",4.0);

  m_msq = _md("ISGW2_msq",1.82);
  m_bx2 = _md("ISGW2_bx2",0.2025);
  m_mbx = _md("ISGW2_mbx",1.975);
  m_nfp = _md("ISGW2_nfp",3.0);
}


void Heavy_Scalar_Current::ISGW2::CalcFFs( Vec4D pB, Vec4D pS )
{
  Vec4D q = pB - pS;
  double q2 = q*q;

  double mtb = m_msb + m_msd;
  double mtx = m_msq + m_msd;  
  double mup = 1.0/(1.0/m_msq+1.0/m_msb);
  double bbx2= 0.5*(m_bb2+m_bx2);
  double tm  = (m_mB-m_mS)*(m_mB-m_mS);
  if ( q2>tm ) { q2=0.99*tm; }
  
  double mqm = 0.1;
  double r2  = 3.0/(4.0*m_msb*m_msq)+3*m_msd*m_msd/(2*m_mbb*m_mbx*bbx2) + 
    (16.0/(m_mbb*m_mbx*(33.0-2.0*m_nfp)))*
    log(Getas(mqm,mqm)/Getas(m_msq,m_msq));
  
  double f3 = sqrt(mtx/mtb)*pow(sqrt(m_bx2*m_bb2)/bbx2,1.5) /
    ((1.0+r2*(tm-q2)/12.0)*(1.0+r2*(tm-q2)/12.0));
  
  double w = 1.0 + (( tm - q2 ) / ( 2.0* m_mB * m_mS ));
  double ai = -1.0* ( 6.0/( 33.0 - 2.0*m_nf));  
  double cji = pow(( Getas( m_msb,m_msb ) / Getas( m_msq,m_msq ) ),ai);
  
  double zji = m_msq / m_msb;
  
  double gammaji = GetGammaji( zji );
  double chiji = -1.0 - ( gammaji / ( 1- zji ));
  double betaji_fppfm = gammaji - (2.0/3.0)*chiji;
  double betaji_fpmfm = gammaji + (2.0/3.0)*chiji;
  double rfppfm = cji *(1.0 + betaji_fppfm*Getas( m_msq,sqrt(m_msb*m_msq) )/M_PI);
  double rfpmfm = cji *(1.0 + betaji_fpmfm*Getas( m_msq,sqrt(m_msb*m_msq) )/M_PI);
  double f3fppfm = f3*pow(( m_mbb / mtb ),-0.5)*pow((m_mbx/mtx),0.5);
  double f3fpmfm = f3*pow(( m_mbb / mtb ),0.5)*pow((m_mbx/mtx),-0.5);
  double fppfm = f3fppfm* rfppfm * ( 2.0 - ( ( mtx/m_msq)*(1- ( (m_msd*m_msq*m_bb2)
						       /(2.0*mup*mtx*bbx2)))));
  double fpmfm = f3fpmfm* rfpmfm * ( mtb/m_msq) * ( 1 - ( ( m_msd*m_msq*m_bb2)/
						 ( 2.0*mup*mtx*bbx2)));
  
  m_Fplus = (fppfm + fpmfm)/2.0;
  double Fminus = (fppfm - fpmfm)/2.0;
  m_F0    = (Fminus/((m_mB*m_mB-m_mS*m_mS)/q2))+(m_Fplus);  
  m_calced = true;
}


double Heavy_Scalar_Current::ISGW2::Getas( double massq, double massx )
{
  double pi = std::acos(-1.0);
  double lqcd2 = 0.04;
  double nflav = 4;
  double temp = 0.6;
  
  if ( massx > 0.6 ) {
    if ( massq < 1.85 ) {
      nflav = 3.0;}
    
    temp = 12.0*pi / ( 33.0 - 2.0*nflav) /
      log( massx*massx/lqcd2);
  }
  return temp;
}


double Heavy_Scalar_Current::ISGW2::GetGammaji ( double z )
{
   double temp;
   temp = 2+((2.0*z)/(1-z))*log(z);
   temp = -1.0*temp;
   return temp;
} 



DECLARE_GETTER(Heavy_Scalar_Current_Getter, "Heavy_Scalar_Current",
               Current_Base,Flavour_Info);

Current_Base* Heavy_Scalar_Current_Getter::operator()(const Flavour_Info &parameters) const
{
  return new Heavy_Scalar_Current(parameters.flavs, parameters.nout, parameters.indices, "Heavy_Scalar_Current");
}

void Heavy_Scalar_Current_Getter::
    PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"implement me";
}

