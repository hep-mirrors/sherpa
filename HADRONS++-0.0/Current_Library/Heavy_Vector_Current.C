#include "Heavy_Vector_Current.H"
#include "Tools.H"

using namespace HADRONS;
using namespace ATOOLS;

void Heavy_Vector_Current::SetModelParameters( struct GeneralModel _md )
{
  m_Vxx = _md("Vxx",0.04);
  switch( int(_md("HV_FORM_FACTOR", 4)+0.5) ) {
  case 1:
    p_ff = new HQET( _md, p_masses );
    msg.Tracking()<<"    Using HQET form factor model for "<<m_name<<std::endl;
    break;
  case 2:
    p_ff = new ISGW( _md, p_masses );
    msg.Tracking()<<"    Using ISGW form factor model for "<<m_name<<std::endl;
    break;
  case 3:
    p_ff = new HQET2( _md, p_masses );
    msg.Tracking()<<"    Using HQET2 form factor model for "<<m_name<<std::endl;
    break;
  case 4:
    p_ff = new ISGW2( _md, p_masses );
    msg.Tracking()<<"    Using ISGW2 for "<<m_name<<std::endl;
    break;
  default:
    msg.Error()<<METHOD<<": You chose a form factor model which does not "
      <<"exist for current "<<m_name<<". Aborting."<<std::endl;
    abort();
  }
}

void Heavy_Vector_Current::Calc()
{
  Vec4D pB = p_moms[0];
  Vec4D pV = p_moms[1];
  Vec4D q = pB - pV;
  double q2 = q*q;
  double mB = p_masses[0];
  double mV = p_masses[1];

  // Get formfactors
  p_ff->CalcFFs(q2);
  double A0 = p_ff->A0();
  double A1 = p_ff->A1();
  double A2 = p_ff->A2();
  double A3 = p_ff->A3();
  double V  = p_ff->V();

  for( int h_had=0; h_had<3; h_had++) {
    ComplexVec4D eps = Tools::ComplexBosonPolarizationVectorC(p_moms[1], h_had);
    Complex i = Complex(0.0,1.0);
    p_results[h_had] = m_Vxx*(  V*2.0*i/(mB+mV) * cross(eps.Conjugate(),pV,pB)
                                - A1*(mB+mV) * eps.Conjugate()
                                + A2*(eps.Conjugate()*q)/(mB+mV) * (pB+pV)
                                + (A3-A0)*2.0*mV*(eps.Conjugate()*q)/q2*q );
  }
}


Heavy_Vector_Current::HQET::HQET( GeneralModel _md, double* _masses )
  : FF_Base( _md, _masses ) 
{
  m_rho2  = _md("rho2",0.769);
  m_R1    = _md("R1",1.328);
  m_R2    = _md("R2",0.92);
  m_RStar = ( 2.0*sqrt(m_mB*m_mV))/(m_mB+m_mV);
}

void Heavy_Vector_Current::HQET::CalcFFs( double q2 )
{
  double w  = (m_mB*m_mB + m_mV*m_mV - q2) / (2.0 * m_mB * m_mV);
  double xi = 1.0 - m_rho2*(w-1.0);

  m_A0 = 0.0;
  m_A1 = (1.0 - q2/(sqr(m_mB+m_mV)))*xi/m_RStar;
  m_A2 = m_R2/m_RStar*xi;
  m_A3 = (m_mB+m_mV)/(2.0*m_mV)*m_A1 - (m_mB-m_mV)/(2.0*m_mV)*m_A2;
  m_V  = m_R1/m_RStar*xi;
  m_calced = true;
}


Heavy_Vector_Current::HQET2::HQET2( GeneralModel _md, double* _masses )
  : FF_Base( _md, _masses ) 
{
  m_rho2  = _md("rho2",1.34);
  m_ha1_1 = _md("ha1_1",0.91);
  m_r1_1  = _md("r1_1",1.18);
  m_r2_1  = _md("r2_1",0.71);
}


void Heavy_Vector_Current::HQET2::CalcFFs( double q2 )
{
  double w  = (m_mB*m_mB + m_mV*m_mV - q2) / (2.0 * m_mB * m_mV);
  const double z = (sqrt(w+1)-sqrt(2.))/(sqrt(w+1)+sqrt(2.)); 
  double ha1 =m_ha1_1*(1.- 8.*m_rho2*z + (53.*m_rho2-15.)*z*z - (231.*m_rho2-91.)*z*z*z);
  double r1 = m_r1_1-0.12*(w-1)+0.05*(w-1)*(w-1);
  double r2 = m_r2_1+0.11*(w-1)-0.06*(w-1)*(w-1);
  double rstar = ( 2.0*sqrt(m_mB*m_mV))/(m_mB+m_mV);

  m_A1 = (1.0 - (q2/((m_mB+m_mV)*(m_mB+m_mV))))*ha1;
  m_A1 = m_A1/rstar;
  m_A2 = (r2/rstar)*ha1;
  m_V  = (r1/rstar)*ha1;
  m_A0 = 0.0;
  m_A3 = (m_mB+m_mV)/(2.0*m_mV)*m_A1 - (m_mB-m_mV)/(2.0*m_mV)*m_A2; 
  m_calced = true;;
}




Heavy_Vector_Current::ISGW::ISGW( GeneralModel _md, double* _masses )
  : FF_Base( _md, _masses ) 
{
  m_msb       = _md("msb",5.2);
  m_msd       = _md("msd",0.33);
  m_bb2       = _md("bb2",0.1681);
  m_msq       = _md("msq",1.82);
  m_bx2       = _md("bx2",0.1521);
  m_kap2      = _md("kap2",0.49);
}

void Heavy_Vector_Current::ISGW::CalcFFs( double q2 )
{
  double mtb  = m_msb+m_msd;
  double mtx  = m_msq+m_msd;  
  double mum  = 1.0/(1.0/m_msq-1.0/m_msb);
  double bbx2 = 0.5*(m_bb2+m_bx2);

  double tm   = (m_mB-m_mV)*(m_mB-m_mV);
  if ( q2 > tm ) q2 = 0.99*tm;
  double F3   = sqrt(mtx/mtb)*pow(sqrt(m_bx2*m_bb2)/bbx2,1.5)*
                exp(-1.0*((m_msd*m_msd*(tm-q2)/(4.0*mtb*mtx*m_kap2*bbx2))));

  double F  = 2.0*mtb*F3;
  double G  = 0.5*F3*((1/m_msq)-(m_msd*m_bb2/(2.0*mum*mtx*bbx2)));
  double AP = (-1.0*F3/(2.0*mtx))*(1.0+(m_msd*(m_bb2-m_bx2)/(m_msb
         *(m_bb2+m_bx2)))-(m_msd*m_msd*m_bx2*m_bx2/(4.0*mum*mtb*bbx2*bbx2)));
  double AM = 0.0;
  m_A1 = F/(m_mB+m_mV);
  m_A2 = -1.0*AP*(m_mB+m_mV);
  m_A3 = (m_mB+m_mV)/(2.0*m_mV)*m_A1 - (m_mB-m_mV)/(2.0*m_mV)*m_A2;
  m_A0 = m_A3 - q2*AM/(2.0*m_mV);
  m_V  = G*(m_mB+m_mV);
  m_calced = true;
}


Heavy_Vector_Current::ISGW2::ISGW2( GeneralModel _md, double* _masses )
  : FF_Base( _md, _masses ) 
{
  m_msb       = _md("msb",5.2);
  m_msd       = _md("msd",0.33);
  m_bb2       = _md("bb2",0.18576);
  m_mbb       = _md("mbb",5.31);
  m_nf        = _md("nf",4.0);
  m_cf        = _md("cf",0.989);
  m_msq       = _md("msq",1.82);
  m_bx2       = _md("bx2",0.14440);
  m_mbx       = _md("mbx",1.9750);
  m_nfp       = _md("nfp",3.0);
}

void Heavy_Vector_Current::ISGW2::CalcFFs( double q2 )
{
  double mtb  = m_msb+m_msd;
  double mtx  = m_msq+m_msd;  
  double mup  = 1.0/(1.0/m_msq+1.0/m_msb);
  double mum  = 1.0/(1.0/m_msq-1.0/m_msb);
  double bbx2 = 0.5*(m_bb2+m_bx2);

  double tm   = (m_mB-m_mV)*(m_mB-m_mV);
  if ( q2 > tm ) q2 = 0.99*tm;
  double wt = 1.0+(tm-q2)/(2.0*m_mbb*m_mbx);
  double mqm = 0.1;
  
  double r2 = 3.0/(4.0*m_msb*m_msq)+3*m_msd*m_msd/(2*m_mbb*m_mbx*bbx2) + 
    (16.0/(m_mbb*m_mbx*(33.0-2.0*m_nfp)))*
    log(Getas(mqm,mqm)/Getas(m_msq,m_msq));

  double w = 1.0 + ((tm - q2 ) / ( 2.0* m_mB * m_mV ));
 
  double ai = -1.0* ( 6.0/( 33.0 - 2.0*m_nf));  
  
  double cji = pow(( Getas( m_msb,m_msb ) / Getas( m_msq,m_msq ) ),ai);
  double zji = m_msq / m_msb;

  double gammaji = GetGammaji( zji );

  double chiji = -1.0 - ( gammaji / ( 1- zji ));
  
  double betaji_g = (2.0/3.0)+gammaji;
  double betaji_f = (-2.0/3.0)+gammaji;
  double betaji_appam = -1.0-chiji+(4.0/(3.0*(1.0-zji)))+
                 (2.0*(1+zji)*gammaji/(3.0*(1.0-zji)*(1.0-zji)));
  
  double betaji_apmam = (1.0/3.0)-chiji-(4.0/(3.0*(1.0-zji)))-
                 (2.0*(1+zji)*gammaji/(3.0*(1.0-zji)*(1.0-zji)))+
                 gammaji;

  double r_g = cji*(1+(betaji_g*Getas( m_msq,sqrt(m_mB*m_msq) )/M_PI));
  double r_f = cji*(1+(betaji_f*Getas( m_msq,sqrt(m_mB*m_msq) )/M_PI));
  double r_appam = cji*(1+(betaji_appam*Getas( m_msq,sqrt(m_mB*m_msq) )/M_PI));
  double r_apmam = cji*(1+(betaji_apmam*Getas( m_msq,sqrt(m_mB*m_msq) )/M_PI));

  
  double f3=sqrt(mtx/mtb)*pow(sqrt(m_bx2*m_bb2)/bbx2,1.5)/
    ((1.0+r2*(tm-q2)/12.0)*(1.0+r2*(tm-q2)/12.0));
  
  double f3f=sqrt(m_mbx*m_mbb/(mtx*mtb))*f3;
  double f3g=sqrt(mtx*mtb/(m_mbx*m_mbb))*f3;
  double f3appam=sqrt(mtb*mtb*mtb*m_mbx/(m_mbb*m_mbb*m_mbb*mtx))*f3;
  double f3apmam=sqrt(mtx*mtb/(m_mbx*m_mbb))*f3;
  double f=m_cf*mtb*(1+wt+m_msd*(wt-1)/(2*mup))*f3f*r_f;
  double g=0.5*(1/m_msq-m_msd*m_bb2/(2*mum*mtx*bbx2))*f3g*r_g;
  
  double appam=cji*(m_msd*m_bx2*(1-m_msd*m_bx2/(2*mtb*bbx2))/ 
	     ((1+wt)*m_msq*m_msb*bbx2)-
	     betaji_appam*Getas( m_msq,sqrt(m_msq*m_mB) )/
	     (mtb*M_PI))*f3appam;
  
  double apmam=-1.0*(mtb/m_msb-m_msd*m_bx2/(2*mup*bbx2)+wt*m_msd*mtb*m_bx2*
	      (1-m_msd*m_bx2/(2*mtb*bbx2))/((wt+1)*m_msq*m_msb*bbx2))*
            f3apmam*r_apmam/mtx;
  
  double ap=0.5*(appam+apmam);
  double am=0.5*(appam-apmam);  

  m_V  = g*(m_mB+m_mV);
  m_A1 = f/(m_mB+m_mV);
  m_A2 = -1.0*ap*(m_mB+m_mV);
  m_A3 = (m_mB+m_mV)/(2.0*m_mV)*m_A1 - (m_mB-m_mV)/(2.0*m_mV)*m_A2;
  m_A0 = m_A3 + ( (q2*am)/(2.0*m_mV));
  m_calced = true;
}

double Heavy_Vector_Current::ISGW2::Getas( double massq, double massx )
{
  double lqcd2 = 0.04;
  double nflav = 4;
  double temp = 0.6;
  
  if ( massx > 0.6 ) {
    if ( massq < 1.85 ) {
      nflav = 3.0;}
    
    temp = 12.0*M_PI / ( 33.0 - 2.0*nflav) /
      log( massx*massx/lqcd2);
  }
  return temp;
}
 

double Heavy_Vector_Current::ISGW2::GetGammaji ( double z )
{
   double temp;
   temp = 2+((2.0*z)/(1-z))*log(z);
   temp = -1.0*temp;
   return temp;
} 


DECLARE_GETTER(Heavy_Vector_Current_Getter, "Heavy_Vector_Current",
               Current_Base,Flavour_Info);

Current_Base* Heavy_Vector_Current_Getter::operator()(const Flavour_Info &parameters) const
{
  return new Heavy_Vector_Current(parameters.flavs, parameters.nout, parameters.indices, "Heavy_Vector_Current");
}

void Heavy_Vector_Current_Getter::
    PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"implement me";
}

