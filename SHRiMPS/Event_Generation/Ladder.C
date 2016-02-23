#include "SHRiMPS/Event_Generation/Ladder.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Random.H"

using namespace SHRIMPS;
using namespace ATOOLS;


std::ostream & SHRIMPS::
operator<<(std::ostream & s, const colour_type::code & colour) {
  if (colour==colour_type::singlet)      s<<" singlet ";
  else if (colour==colour_type::triplet) s<<" triplet ";
  else if (colour==colour_type::octet)   s<<"  octet  ";
  else                                   s<<"   none  ";
  return s;
}

std::ostream & SHRIMPS::
operator<<(std::ostream & s, const Ladder_Particle & part) {
  s<<"   "<<part.m_flav<<"  "<<part.m_mom<<" "
   <<"(y="<<part.m_mom.Y()<<", pt^2="<<part.m_mom.PPerp2()<<") "
   <<"{"<<part.GetFlow(1)<<" "<<part.GetFlow(2)<<"}"
   <<" at "<<part.m_pos<<".\n";
  return s;
}

Ladder_Particle::Ladder_Particle(const Particle * part) :
  p_part(NULL), m_flav(part->Flav()), m_mom(part->Momentum()), 
  m_flow(Flow(NULL)), m_marked(false), 
  m_IS(part->Info()=='I' || part->Info()=='G')
{
  SetFlow(1,part->GetFlow(1));
  SetFlow(2,part->GetFlow(2));
  m_beam = dabs(m_mom.Y())>5.;
  
  m_pos = part->ProductionBlob()?
    part->ProductionBlob()->Position()/(rpa->hBar()*rpa->c()):
    Vec4D(0.,0.,0.,0.);
}

Particle * Ladder_Particle::GetParticle() {
  msg_Tracking()<<METHOD<<"("<<this<<"  --> "<<p_part<<")\n";
  if (!p_part) { 
    p_part = new Particle(-1,m_flav,m_mom,m_IS?'I':'F');
    p_part->SetNumber();
    p_part->SetFlow(1,GetFlow(1,false));
    p_part->SetFlow(2,GetFlow(2,false));
  }
  return p_part;
}

std::ostream & SHRIMPS::operator<<(std::ostream & s, const T_Prop & tprop) {
  s<<" | ["<<tprop.m_col<<"] for "
   <<"q = "<<tprop.m_q<<" (qt = "<<sqrt(tprop.m_qt2)<<", q = "
   <<sqrt(dabs(tprop.m_q.Abs2()))<<")"
   <<" and Q0 = "<<sqrt(tprop.m_q02)<<" | \n";
  return s;
}

std::ostream & SHRIMPS::operator<<(std::ostream & s, const TPropList & props) {
  s<<"T propagator list ("<<props.size()<<", "<<(&props)<<"): \n";
  if (props.size()>0) {
    for (TPropList::const_iterator piter=props.begin();
	 piter!=props.end();piter++) s<<(*piter);
  }
  s<<"\n";
  return s;
}

std::ostream & SHRIMPS::
operator<<(std::ostream & s,const Ladder & ladder) {
  s<<"   ---------------------------------------------------------\n"
   <<"Ladder ("<<ladder.m_tprops.size()<<" props, "<<&ladder.m_tprops<<") "
   <<"at position "<<ladder.m_position<<" (b= "
   <<(sqrt(sqr(ladder.m_position[1])+sqr(ladder.m_position[2])))<<"), "
   <<"kt^2 = "<<ladder.MaxKT2()<<", wt = "<<ladder.Weight()<<":\n"
   <<" * "<<(*ladder.p_inpart1)<<" * "<<(*ladder.p_inpart2)<<"\n";
  int i(0);
  TPropList::const_iterator citer=ladder.m_tprops.begin();
  for (LadderMap::const_iterator yiter=ladder.m_emissions.begin();
       yiter!=ladder.m_emissions.end();yiter++) {
    s<<"  y_{"<<i<<"} = "<<yiter->first<<", k_{"<<i<<"} = "
     <<yiter->second<<"\n";
    if (citer!=ladder.m_tprops.end()) {
      s<<(*citer)<<"\n";
      citer++;
    }
    i++;
  }
  s<<"   ---------------------------------------------------------\n";
  return s;
}

Ladder::Ladder(const Vec4D & position) :
  m_position(position), m_weight(1.), 
  m_Y(0.), m_Ymin(0.), m_Ymax(0.), m_deltaY(0.), 
  m_diffractive(false), m_rescatter(false), m_harddiffractive(false),  
  m_maxkt2(0.), m_minkt2(0.), m_shat(0.), m_that(0.), m_uhat(0.), m_mu2(0.), 
  p_inpart1(NULL), p_inpart2(NULL), m_enforceup(false)
{ }

Ladder::~Ladder() {
  msg_Debugging()<<METHOD<<" delete Ladder_Particle("
		 <<p_inpart1<<" "<<p_inpart2<<").\n";
  if (p_inpart1) { delete p_inpart1; p_inpart1 = 0; }
  if (p_inpart2) { delete p_inpart2; p_inpart2 = 0; }
}


void Ladder::UpdatePropagators() {
  Vec4D q(p_inpart1->m_mom);
  double q2, qt2;
  TPropList::iterator piter=GetPropsBegin();
  for (LadderMap::iterator liter=GetEmissionsBegin();
       liter!=GetEmissionsEnd();liter++) {
    q  -= liter->second.m_mom;
    q2  = dabs(q.Abs2());
    qt2 = q.PPerp2();
    piter->m_q   = q;
    piter->m_q2  = q2;
    piter->m_qt2 = qt2; 
    piter++;
    if (piter==GetPropsEnd()) break;
  }
}

double Ladder::MRKweight() const {
  //std::cout<<METHOD<<":\n";
  double mrkweight(1.), qt2,q2;
  if (m_tprops.size()>1) {
    for (TPropList::const_iterator pit=m_tprops.begin();
	 pit!=m_tprops.end();pit++) {
      qt2 = pit->m_q.PPerp2(); 
      q2  = pit->m_q.Abs2();
      mrkweight *= qt2/ATOOLS::Max(qt2,q2);
    }
  }
  return ATOOLS::dabs(mrkweight);
}

bool Ladder::ExtractHardest() {
  UpdatePropagators();
  m_that = m_mu2 = -1.;
  m_harddiffractive = false;
  LadderMap::iterator liter=GetEmissionsBegin(), outiter1;
  TPropList::iterator initer;
  for (TPropList::iterator piter=GetPropsBegin();
       piter!=GetPropsEnd();piter++) {
     if (dabs(piter->m_q2)>m_that) {
       outiter1 = liter;
       initer   = piter;
       m_that   = dabs(initer->m_q2);
       if (piter->m_col==colour_type::singlet) m_harddiffractive = true;
       else m_harddiffractive = false;
       //msg_Tracking()<<METHOD<<" finds hardest: that = "<<m_that<<" "
       //	<<"("<<piter->m_col<<" -> "<<m_harddiffractive<<").\n";
     }
     liter++;
  }
  if (m_that==-1.) return false;
  p_hardestprop = &(*initer);
  m_propcol   = initer->m_col;
  m_that      = dabs(initer->m_q2);
  m_mu2       = initer->m_q02;
  liter       = outiter1; liter++;
  m_shat      = (outiter1->second.m_mom+liter->second.m_mom).Abs2();
  m_Yhat      = (outiter1->second.m_mom+liter->second.m_mom).Y();
  m_DeltaYhat = outiter1->second.m_mom.Y()-liter->second.m_mom.Y();
  double q12(0.), q22(0.);
  if (initer!=GetPropsBegin()) {
    initer--;
    q12 = initer->m_q2;
    initer++;
  }
  initer++;
  if (initer!=GetPropsEnd()) {
    q22 = initer->m_q2;
  }
  initer--;
  m_uhat    = m_shat-m_that+q12+q22;
  SetMEFSParticles(&outiter1->second,&liter->second);
  return true;
}
 
bool Ladder::
ReconstructMEFlavours(Flavour & i1,Flavour & i2,
		      Flavour & o1,Flavour & o2) {
  o1 = p_part1->m_flav;  o2 = p_part2->m_flav;
  //msg_Out()<<METHOD<<":"<<o1<<" + "<<o2<<" with "<<m_propcol<<".\n";
  if (o1.IsGluon() && o2.IsGluon()) {
    if (m_propcol==colour_type::singlet ||
	m_propcol==colour_type::octet) {
      i1 = i2 = Flavour(kf_gluon);
      return true;
    }
    else if (m_propcol==colour_type::triplet) {
      i1 = Flavour(kf_u);
      i2 = Flavour(kf_u).Bar();
      return true;
    }
    return false;
  }
  else if ((o1.IsGluon() && o2.IsQuark()) || (o1.IsQuark() && o2.IsGluon())) {
    if (m_propcol==colour_type::singlet ||
	m_propcol==colour_type::octet) {
      i1 = o1;
      i2 = o2;
      return true;
    }
    else if (m_propcol==colour_type::triplet) {
      i1 = o2;
      i2 = o1;
      return true;
    }
    return false;
  }
  else if (o1.IsQuark() && o2.IsQuark()) {
    if ((!o1.IsAnti() && o2.IsAnti()) || (o1.IsAnti() && !o2.IsAnti())) {
      if (m_propcol==colour_type::singlet ||
	  m_propcol==colour_type::octet) {
	i1 = o1;
	i2 = o2;
	return true;
      }
      else if (m_propcol==colour_type::triplet && o1==o2.Bar()) {
	i1 = i2 = Flavour(kf_gluon);
	return true;
      }
      return false;      
    }
    else if ((!o1.IsAnti() && !o2.IsAnti()) || (o1.IsAnti() && o2.IsAnti())) {
      if (m_propcol==colour_type::singlet ||
	  m_propcol==colour_type::octet) {
	i1 = o1;
	i2 = o2;
	return true;
      }
      return false;
    }
  }
  return false;
}


bool Ladder::CheckFourMomentum() {
  Vec4D check(p_inpart1->m_mom+p_inpart2->m_mom);
  double shat(check.Abs2());
  TPropList::iterator prop = m_tprops.begin();
  for (LadderMap::iterator liter=m_emissions.begin();
       liter!=m_emissions.end();liter++) {
    check -= liter->second.m_mom;
    if (prop!=m_tprops.end()) {
      if ((check.Perp()-prop->m_q.Perp()).Abs2()>1.e-6) {
	msg_Error()<<"-------------------------------------------\n"
		   <<METHOD<<" failed: check = "<<check<<" vs "<<prop->m_q<<"\n"
		   <<(*this)<<"\n"<<p_inpart1->m_mom<<" / "
		   <<p_inpart2->m_mom<<".\n";
      }
      prop++;
    }
  }
  if (dabs(check.Abs2())/shat>1.e-6) {
    msg_Error()<<"-------------------------------------------\n"
	       <<METHOD<<" failed: check = "<<check<<", "<<check.Abs2()<<"\n"
	       <<(*this)<<"\n"<<p_inpart1->m_mom<<" / "
	       <<p_inpart2->m_mom<<".\n";
    return false;
  }
  return true;
}



bool Ladder::SwapColourIndices() {
  if (p_inpart1->m_flav.IsQuark() || p_inpart2->m_flav.IsQuark()) return false;
  for (LadderMap::iterator liter=GetEmissionsBegin();
       liter!=GetEmissionsEnd();liter++) { 
    liter->second.m_flow.SwapColourIndices();
    if (liter->second.m_flav.IsQuark()) {
      liter->second.m_flav = liter->second.m_flav.Bar();
    }
  }
  p_inpart1->m_flow.SwapColourIndices();
  p_inpart2->m_flow.SwapColourIndices();
  return true;
}

bool Ladder::CanReplaceColour(const size_t & oldc,const size_t & newc,
			      const size_t & pos,const bool & inclIS) {
  Ladder_Particle * part;
  for (LadderMap::iterator liter=GetEmissionsBegin();
       liter!=GetEmissionsEnd();liter++) {
    part = &liter->second;
    if (part->GetFlow(pos)==oldc) {
      if (part->GetFlow(3-pos)!=newc) {
	part->SetFlow(pos,newc);
	return true;
      }
    }
  }
  if (p_inpart1->GetFlow(3-pos)==oldc && p_inpart1->GetFlow(pos)!=newc) {
    p_inpart1->SetFlow(3-pos,newc);
    return true;
  }
  if (p_inpart2->GetFlow(3-pos)==oldc && p_inpart2->GetFlow(pos)!=newc) {
    p_inpart2->SetFlow(3-pos,newc);
    return true;
  }
  return false;
}

bool Ladder::GenerateColourIndices(size_t & fix) {
  msg_Out()<<"#############################################\n"
		<<METHOD<<"(fix = "<<fix<<"):\n";
  LadderMap::iterator lbeg=GetEmissionsBegin(), * liter=&lbeg;
  LadderMap::iterator lend=GetEmissionsEnd();
  lend--;
  TPropList::iterator citer=GetPropsBegin();
  int col1(0), col2(0);
  if (!FixFirstColours(*liter,col1,col2,fix,citer)) return false;
  while (lend->first-lbeg->first>0.0001) {
    msg_Out()<<"   before intermediate ("<<(*liter)->first<<") with "
    	     <<"["<<col1<<", "<<col2<<"] for "<<citer->m_col<<".\n";
    if (!FixIntermediateColours(*liter,col1,col2,fix,citer)) return false;
  }
  msg_Out()<<"   before last ("<<(*liter)->first<<") with "
  	   <<"["<<col1<<", "<<col2<<"] for "<<citer->m_col<<".\n";
  if (!FixLastColours(*liter,col1,col2,fix,citer)) return false;
  msg_Out()<<METHOD<<"(fix = "<<fix<<"):\n"<<(*this)
		<<"#############################################\n";
  return true;
}



bool Ladder::FixFirstColours(LadderMap::iterator & liter,int & col1,int & col2,
			     const size_t & fix,TPropList::iterator & citer) {
  colour_type::code colour(citer->m_col);
  Ladder_Particle * inpart(p_inpart1);
  Ladder_Particle * outpart(&liter->second);

  msg_Out()<<METHOD<<"("<<liter->first<<"; "
  	   <<inpart->m_flav<<" & "<<outpart->m_flav<<") "
  	   <<"for "<<colour<<" fix = "<<fix<<".\n";

  if (colour==colour_type::singlet) {
    outpart->SetFlow(1,inpart->GetFlow(1));
    outpart->SetFlow(2,inpart->GetFlow(2));
    if (MoreSinglets(citer)) col1 = col2 = -1;
    else {
      if (fix==1) {
	col1 = -1;
	col2 = p_inpart2->GetFlow(2);
      }
      else if (fix==2) {
	col1 = p_inpart2->GetFlow(1);
	col2 = -1;
      }
      else {
	msg_Error()<<"Error in "<<METHOD<<":\n"
		   <<"   No fix = "<<fix<<" declared.  Don't know what to do.\n"
		   <<"   Will hope for the best.\n";
	col1 = col2 = -1;
      }
    }
  }
  else if (colour==colour_type::octet) {
    if (inpart->m_flav != outpart->m_flav || fix==0) {
      msg_Error()<<"Error in "<<METHOD<<":\n"
		 <<"   Octet propagator with flavour mismatch: "
		 <<inpart->m_flav<<" --> "<<outpart->m_flav<<".\n";
      return false;
    }
    if (fix==1) {
      outpart->SetFlow(1,-1);
      if (inpart->m_flav.IsQuark()) outpart->SetFlow(2,0);
      if (inpart->m_flav.IsGluon()) outpart->SetFlow(2,inpart->GetFlow(2));
      col1 = inpart->GetFlow(1);
      col2 = outpart->GetFlow(1);
    }
    else if (fix==2) {
      if (inpart->m_flav.IsQuark()) outpart->SetFlow(1,0);
      if (inpart->m_flav.IsGluon()) outpart->SetFlow(1,inpart->GetFlow(1));
      outpart->SetFlow(2,-1);
      col1 = outpart->GetFlow(2);
      col2 = inpart->GetFlow(2);
    }
  }
  liter++;
  citer++;
  if (citer==GetPropsEnd()) citer--;
  msg_Out()<<"["<<inpart->GetFlow(1)<<", "<<inpart->GetFlow(2)<<"] --> "
		<<"colour = "<<"["<<col1<<", "<<col2<<"] + "
		<<"["<<outpart->GetFlow(1)<<", "<<outpart->GetFlow(2)<<"] for "
		<<liter->first<<"\n";
  return true;
}

bool Ladder::
FixIntermediateColours(LadderMap::iterator & liter,int & col1,int & col2,
		       size_t & fix,TPropList::iterator & citer) {
  msg_Out()<<METHOD<<"(fix = "<<fix<<") "
		<<"with ["<<col1<<", "<<col2<<"]\n";

  colour_type::code colour(citer->m_col);
  Ladder_Particle * outpart(&liter->second);
  Ladder_Particle * refpart(p_inpart2);

  bool fromsing(false);
  if (colour==colour_type::singlet) {
    outpart->SetFlow(1,col1);
    outpart->SetFlow(2,col2);
    if (MoreSinglets(citer)) col1 = col2 = -1;
    else {
      if (fix==1) {
	col1 = -1;
	col2 = p_inpart2->GetFlow(2);
      }
      else if (fix==2) {
	col1 = p_inpart2->GetFlow(1);
	col2 = -1;
      }
      else {
	msg_Error()<<"Error in "<<METHOD<<":\n"
		   <<"   No fix = "<<fix<<" decalred.  Don't know what to do.\n"
		   <<"   Will set it to 0 & hope for the best.\n";
	fix = 0;
	col1 = col2 = -1;
      }
    }
  }
  else if (colour==colour_type::octet) {
    if (col1==-1 && col2==-1) {
      msg_Out()<<" ... from singlet.\n";
      fromsing=true;
    }
    if (outpart->m_flav.IsGluon()) {
      if (fix==1) {
	outpart->SetFlow(1,-1);
	outpart->SetFlow(2,col2);
	col1 = fromsing?outpart->GetFlow(2):col1;
	col2 = outpart->GetFlow(1);
      }
      else if (fix==2) {
	outpart->SetFlow(1,col1);
	outpart->SetFlow(2,-1);
	col1 = outpart->GetFlow(2);
	col2 = fromsing?outpart->GetFlow(1):col2;
      }
      else if (fix==0) {
	if (refpart->m_flav.IsQuark()) {
	  if (!refpart->m_flav.IsAnti()) {
	    outpart->SetFlow(1,refpart->GetFlow(1));
	    outpart->SetFlow(2,-1);
	    col1 = outpart->GetFlow(2);
	    col2 = outpart->GetFlow(1);
	    fix = 2;
	  }
	  else if (refpart->m_flav.IsAnti()) {
	    outpart->SetFlow(1,-1);
	    outpart->SetFlow(2,refpart->GetFlow(2));
	    col1 = outpart->GetFlow(2);
	    col2 = outpart->GetFlow(1);
	    fix = 1;
	  }
	}
	else if (refpart->m_flav.IsGluon()) {
	  fix = ran->Get()<0.5?1:2;
	  if (fix==2) {
	    outpart->SetFlow(1,refpart->GetFlow(1));
	    outpart->SetFlow(2,-1);
	    col1 = outpart->GetFlow(2);
	    col2 = outpart->GetFlow(1);
	  }
	  else if (fix==1) {
	    outpart->SetFlow(1,-1);
	    outpart->SetFlow(2,refpart->GetFlow(2));
	    col1 = outpart->GetFlow(2);
	    col2 = outpart->GetFlow(1);
	  }
	}
      }
    }
  }
  liter++; 
  citer++;
  if (citer==GetPropsEnd()) citer--;
  msg_Out()<<" --> "<<outpart->m_flav<<" "
  	   <<"["<<outpart->GetFlow(1)<<", "<<outpart->GetFlow(2)<<"] -> "
  	   <<colour<<" ["<<col1<<", "<<col2<<"]\n";
  return true;
}

bool Ladder::
FixLastColours(LadderMap::iterator & liter,const int & col1,const int & col2,
	       const size_t & fix,TPropList::iterator & citer) {
  Ladder_Particle * inpart(p_inpart2);
  Ladder_Particle * outpart(&liter->second);
  colour_type::code colour(citer->m_col);

  msg_Out()<<METHOD<<"(fix = "<<fix<<", "<<liter->first<<"; "
  	   <<inpart->m_flav<<" & "<<outpart->m_flav<<").\n";

  if (colour==colour_type::singlet) {
    outpart->SetFlow(1,inpart->GetFlow(1));
    outpart->SetFlow(2,inpart->GetFlow(2));
  }
  else if (colour==colour_type::triplet) {
    if (inpart->m_flav.IsGluon() && outpart->m_flav.IsQuark()) {
      if (outpart->m_flav.IsAnti()) {
	outpart->SetFlow(1,inpart->GetFlow(1));
	outpart->SetFlow(2,0);
      }
      else if (!outpart->m_flav.IsAnti()) {
	outpart->SetFlow(1,inpart->GetFlow(1));
	outpart->SetFlow(2,0);
      }
      else {
	msg_Error()<<"Error: "<<outpart->m_flav<<" for triplet propagator.\n";
	return false;
      }
    }
    else if (inpart->m_flav.IsQuark() && outpart->m_flav.IsGluon()) {
      if (inpart->m_flav.IsAnti()) {
	outpart->SetFlow(1,col2);
	outpart->SetFlow(2,inpart->GetFlow(2));
      }
      else if (!inpart->m_flav.IsAnti()) {
	outpart->SetFlow(1,inpart->GetFlow(1));
	outpart->SetFlow(2,col1);
      }
      else {
	msg_Error()<<"Error: "<<outpart->m_flav<<" for triplet propagator.\n";
	return false;
      }
    }
  }
  else if (colour==colour_type::octet) {
    if (inpart->m_flav != outpart->m_flav) {
      msg_Error()
	<<"Error in "<<METHOD<<":\n"
	<<"   Octet propagator with flavour mismatch: "
	<<inpart->m_flav<<" --> "<<outpart->m_flav<<".\n";
      return false;
    }
    else if (inpart->m_flav.IsQuark()) {
      if (inpart->m_flav.IsAnti()) {
	outpart->SetFlow(1,0);
	outpart->SetFlow(2,col2);
      }
      else if (!inpart->m_flav.IsAnti()) {
	outpart->SetFlow(1,col1);
	outpart->SetFlow(2,0);
      }
    }
    else if (inpart->m_flav.IsGluon()) {
      if (fix==1) {
	outpart->SetFlow(1,inpart->GetFlow(1));
	outpart->SetFlow(2,col2);
      }
      else if (fix==2) {
	outpart->SetFlow(1,col1);
	outpart->SetFlow(2,inpart->GetFlow(2));
      }
    }
  }
  msg_Out()<<"Out of "<<METHOD<<" with:\n"<<(*this)<<"\n";
  return true;
}
 
bool Ladder::MoreSinglets(const TPropList::iterator & citer) {
  TPropList::iterator cit(citer); cit++;
  while (cit!=GetPropsEnd()) {
    if (cit->m_col==colour_type::singlet) return true;
    cit++;
  }
  return false;
}
 
