#include "SHRiMPS/Beam_Remnants/Colour_Generator.H"
#include "SHRiMPS/Beam_Remnants/Hadron_Dissociation.H"
#include "ATOOLS/Phys/Flow.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace SHRIMPS;
using namespace ATOOLS;
using namespace std;


Colour_Generator::
Colour_Generator(vector<Hadron_Dissociation *> * hadrons) {
  for (size_t i=0;i<hadrons->size();i++) 
    m_hadrons.push_back((*hadrons)[i]);
}

Colour_Generator::~Colour_Generator() {}

bool Colour_Generator::
operator()(Ladder * ladder,Particle ** parts,const size_t & counter) {
  p_ladder = ladder;
  if (p_ladder->IsRescatter()) return Rescatter();
  return Primary(parts,counter);
}

bool Colour_Generator::Rescatter() {
  p_compensator = NULL;
  if (p_ladder->IsDiffractive()) return FixUncorrelatedIndices();
  return FixCorrelatedIndices();
}

bool Colour_Generator::FixUncorrelatedIndices() {
  if (!p_ladder->IsDiffractive()) {
    p_ladder->MakeDiffractive();
//     msg_Out()<<METHOD<<" made ladder diffractive:\n"<<(*p_ladder)<<"\n";
  }
  size_t fix(FixAColourLine());
  p_ladder->GenerateColourIndices(fix);
/*  msg_Out()<<"-----------------------------------------------------\n"
  	   <<METHOD<<":\n"<<(*p_ladder)<<"\n"
  	   <<"-----------------------------------------------------\n";*/
  return true;
}

size_t Colour_Generator::FixAColourLine() {
  size_t fix(0);
  if (p_inpart1->Flav().IsQuark()) { 
    if (p_inpart1->Flav().IsAnti()) fix = 2; else fix = 1;
  }
  else if (p_inpart2->Flav().IsQuark()) {
    if (p_inpart2->Flav().IsAnti()) fix = 1; else fix = 2;
  }
  else fix = ran->Get()>0.5?1:2;
/*  msg_Out()<<METHOD
  	   <<"("<<p_inpart1->Flav()<<" "<<p_inpart2->Flav()<<"), "
  	   <<"fix= "<<fix<<".\n";*/
  return fix;
}

bool Colour_Generator::FixCorrelatedIndices() {
  Particle * part1(p_ladder->GetIn1()->GetParticle());
  Particle * part2(p_ladder->GetIn2()->GetParticle());
  size_t fix(ColourConnected(part1,part2));
  int mod;
  //   msg_Out()<<METHOD<<" after first check fix = "<<fix<<".\n";
  switch (fix) {
  case 0: {
    /*    msg_Out()<<"********** fix = "<<fix
	  <<"["<<part1->GetFlow(1)<<" "<<part1->GetFlow(2)<<"] "
	  <<"["<<part2->GetFlow(1)<<" "<<part2->GetFlow(2)<<"].\n";*/
    fix = SelectColourReplacement(part1,part2);
    if (fix==0) {
      mod = ModifyOriginators(part1,part2);
      if (mod<0) return false;
      else fix = size_t(mod); 
    }
    if (fix==0) return FixUncorrelatedIndices();
    break;
  }
  case 3:
    fix = ran->Get()>=0.5?2:1;
    break;
  default:
    break;
  }
//   msg_Out()<<METHOD<<"(fix = "<<fix<<"):\n";
  if (part1==p_ladder->GetIn2()->GetParticle()) {
    fix = 3-fix; 
//     msg_Out()<<"   twisted initial states:modify fix = "<<fix<<".\n";
  }
  p_ladder->GenerateColourIndices(fix);
//   msg_Out()<<(*p_ladder);
  return true;
}

int Colour_Generator::ModifyOriginators(Particle * part1,Particle * part2) {
  bool beam(ran->Get()>0.5);
  Particle * part(beam?part1:part2), * spec(beam?part2:part1);
  Blob * partblob((beam?p_inpart1:p_inpart2)->ProductionBlob());
  Blob * specblob((beam?p_inpart2:p_inpart1)->ProductionBlob());
  size_t index(1+int(ran->Get()>0.5)),col, col1, adj, adj1;
  
  /*    msg_Out()<<"-----------------------------------------------------\n"
	<<METHOD<<"("<<beam<<") for "
	<<"["<<part1->GetFlow(1)<<", "<<part1->GetFlow(2)<<"] "
	<<"["<<part2->GetFlow(1)<<", "<<part2->GetFlow(2)<<"], "
	<<"blobs = "<<partblob->Id()<<" & "<<specblob->Id()<<".\n";
	msg_Out()<<(*p_ladder)<<"\n";
	msg_Out()<<" Deal with {"<<part->Number()<<", "<<spec->Number()<<"}:\n"
	<<(*part)<<"\n"<<(*spec)<<"\n";*/
  
  for (int i=0;i<2;i++) {
    index = 3-index;
    col  = part->GetFlow(index);
    adj  = part->GetFlow(3-index);
    col1 = spec->GetFlow(3-index);
    adj1 = spec->GetFlow(index);
    //     msg_Out()<<i<<"th attempt: 
    //Check to replace "<<col<<" <-> "<<col1<<".\n";
    if (col==0 || col1==0) continue;
    if (ReplaceColour(partblob,index,col,col1,adj)) {
      part->SetFlow(index,col1);
      if ((part->Flav().IsGluon() && col1==0) ||
	  (part->Flav().IsQuark() && part->Flav().IsAnti() && 
	   index==2 && col1==0) ||
	  (part->Flav().IsQuark() && !part->Flav().IsAnti() && 
	   index==1 && col1==0))
	msg_Error()<<"Error with colours here (1).\n";
      (beam?p_ladder->GetIn1():p_ladder->GetIn2())->m_flow.SetCode(index,col1);
      
      /*	msg_Out()<<"index(1) = "
		<<(beam?index:3-index)<<" for beam = "<<beam<<": "
		<<col<<" ---> "<<col1<<", "
		<<"set "<<part->Number()
		<<"(index = "<<index<<") = "<<col1<<".\n";*/
      
      if (!partblob->CheckColour()) {
	msg_Error()<<"Problem in "<<METHOD<<":\n"
		   <<"   Extra blob ("<<partblob->Id()<<") seems fishy: "
		   <<"Bad colour configuration.\n"<<(*partblob)<<"\n";
	return -1;
      }
      return beam?index:3-index;
    } 
    if (ReplaceColour(specblob,3-index,col1,col,adj1)) {
      (beam?p_ladder->GetIn2():p_ladder->GetIn1())->m_flow.SetCode(3-index,col);
      if ((spec->Flav().IsGluon() && col1==0) ||
	  (spec->Flav().IsQuark() && spec->Flav().IsAnti() && 
	   index==2 && col1==0) ||
	  (spec->Flav().IsQuark() && !spec->Flav().IsAnti() && 
	   index==1 && col1==0))
	msg_Error()<<"Error with colours here (2).\n";
      spec->SetFlow(3-index,col);
      
      /*	msg_Out()<<"index(2) = "
		<<(beam?3-index:index)<<" for beam = "<<beam<<": "
		<<col1<<" ---> "<<col<<", "
		<<"set "<<spec->Number()
		<<"(index = "<<(3-index)<<") = "<<col<<"\n"
		<<(*spec)<<"\n"<<(*part)<<"\n";*/
      
      if (!specblob->CheckColour()) {
	msg_Error()<<"Problem in "<<METHOD<<":\n"
		   <<"   Extra blob ("<<specblob->Id()<<") seems fishy: "
		   <<"Bad colour configuration.\n"<<(*specblob)<<"\n";
	return -1;
      }
      return beam?index:3-index;
    } 
  }
  return 0;
}

bool Colour_Generator::
ReplaceColour(Blob * blob,const size_t & index,
	      const size_t & col,const size_t & col1,const size_t & adj) 
{
  Particle * parti,* partj;
  for (int i=0;i<blob->NInP();i++) {
    if (blob->InParticle(i)->GetFlow(index)==col) return false;
  }
  for (int i=0;i<blob->NOutP();i++) {
    parti = blob->OutParticle(i);
    if (parti->GetFlow(index)==col &&
	parti->GetFlow(3-index)==adj) {
      if (parti->GetFlow(3-index)==col1) return false;
      for (int j=0;j<blob->NOutP();j++) {
	if (i==j) continue;
	partj = blob->OutParticle(j);
	if (partj->DecayBlob()) continue;
	if (partj->GetFlow(3-index)==col) {
	  if (partj->GetFlow(index)==col1) return false;
	  
	  /*	    msg_Out()<<METHOD
		    <<"("<<parti->Number()<<" & "<<partj->Number()<<") "
		    <<"successful.  Decayblobs = "
		    <<(parti->DecayBlob()?parti->DecayBlob()->Id():0)
		    << " & "
		    <<(partj->DecayBlob()?partj->DecayBlob()->Id():0)
		    << ".\n"<<(*blob)<<"\n";*/
	  
	  if ((parti->Flav().IsGluon() && col1==0) ||
	      (parti->Flav().IsQuark() && parti->Flav().IsAnti() && 
	       index==2 && col1==0) ||
	      (parti->Flav().IsQuark() && !parti->Flav().IsAnti() && 
	       index==1 && col1==0))
	    msg_Error()<<"Error with colours here (3i).\n";
	  if ((partj->Flav().IsGluon() && col1==0) ||
	      (partj->Flav().IsQuark() && partj->Flav().IsAnti() && 
	       index==2 && col1==0) ||
	      (partj->Flav().IsQuark() && !partj->Flav().IsAnti() && 
	       index==1 && col1==0))
	    msg_Error()<<"Error with colours here (3j).\n";
	  parti->SetFlow(index,col1);
	  partj->SetFlow(3-index,col1);
	  // 	  msg_Out()<<(*parti)<<"\n"<<(*partj)<<"\n";
	  if (!blob->CheckColour()) {
	    msg_Error()<<"Problem in "<<METHOD<<":\n"
		       <<"   Extra blob ("<<blob->Id()<<") seems fishy: "
		       <<"Bad colour configuration.\n"<<(*blob)<<"\n";
	    return false;
	  }
	  return true;
	}
      }
    }
  }
  /*  msg_Out()<<METHOD<<"("<<col<<", "<<col1<<") for blob ["<<blob->Id()<<"]: "
      <<"Colour not found in index = "<<index<<".\n";*/
  return false;
}

size_t Colour_Generator::
SelectColourReplacement(Particle * part1,Particle * part2) {
  bool dir(ran->Get()>=0.5),found(false);
  size_t beam,index,col,col1,idcol;
  for (size_t i=0;i<2;i++) {
    beam = dir?i:1-i;
    for (index = 1;index<3;index++) {
      col = (beam==1?part2:part1)->GetFlow(index);
      for (set<int>::iterator cit=m_col[beam][2-index].begin();
	   cit!=m_col[beam][2-index].end();cit++) {
	if (col==(*cit)) {
	  col1 = (beam==1?part1:part2)->GetFlow(3-index);
	  if (col1!=0) {
	    found = false;
            for (int i=0;i<m_hadrons[1-beam]->GetParticles().size();i++){ 
	      if (m_hadrons[1-beam]->GetParticles()[i]->GetFlow(index)==col && 
		  m_hadrons[1-beam]->GetParticles()[i]->GetFlow(3-index)!=col1){
		found = true;
		idcol = m_hadrons[1-beam]->GetParticles()[i]->GetFlow(3-index);
		break;
	      }
	    }	      
	    break;
	  }
	}
      }
      if (found) break;
    }
    if (found) break;
  }
  if (!found) return 0;

  p_compensator = new Blob(p_ladder->Position()*rpa->hBar()*rpa->c());
  p_compensator->SetId();
  p_compensator->SetType(btp::Soft_Collision);
  p_compensator->SetTypeSpec("ColourCompensation");
  p_compensator->SetStatus(blob_status::inactive);
  p_compensator->AddToInParticles((beam==1?part2:part1));
  (beam==1?part2:part1)->SetStatus(part_status::decayed);

  Particle * gluon = new Particle(-1,Flavour(kf_gluon),Vec4D(0.,0.,0.,0.),'R');
  gluon->SetNumber();
  gluon->SetFlow(index,col1);
  gluon->SetFlow(3-index,col);
  if (col1==0 || col1==0)
    msg_Error()<<"Error with colours here (4).\n";

  gluon->SetStatus(part_status::decayed);
  p_compensator->AddToInParticles(gluon);
  p_blob->AddToOutParticles(gluon);

  Particle * part= new Particle(*((beam==1?part2:part1)));
  part->SetNumber();
  part->SetFlow(index,col1);
  if ((part->Flav().IsGluon() && col1==0) ||
      (part->Flav().IsQuark() && part->Flav().IsAnti() && 
       index==2 && col1==0) ||
      (part->Flav().IsQuark() && !part->Flav().IsAnti() && 
       index==1 && col1==0))
    msg_Error()<<"Error with colours here (5).\n";
  (beam==1?p_ladder->GetIn2():p_ladder->GetIn1())->m_flow.SetCode(index,col1);
  (beam==1?p_ladder->GetIn2():p_ladder->GetIn1())->p_part = part;

  for (int i=0;i<p_blob->NInP();i++) {
    part = p_blob->InParticle(i);
    if (part->GetFlow(index)==col && part->GetFlow(3-index)==idcol && 
	part->Beam()==(1-beam)){
      part->SetFlow(index,col1);
      break;
    }
    if ((part->Flav().IsGluon() && col1==0) ||
	(part->Flav().IsQuark() && part->Flav().IsAnti() && 
	 index==2 && col1==0) ||
	(part->Flav().IsQuark() && !part->Flav().IsAnti() && 
	 index==1 && col1==0))
      msg_Error()<<"Error with colours here (6).\n";
  }

  m_col[beam][2-index].erase(col);
  m_col[beam][2-index].insert(col1);

  return (beam==1?3-index:index);
}

bool Colour_Generator::Primary(Particle ** parts,const size_t & counter) {
  size_t fix;
  int col[2][2], colour;
  col[0][0] = col[0][1] = col[1][0] = col[1][1] = -1;
  if (p_ladder->IsDiffractive()) fix = UncorrelatedIndices(col);
  else                           fix = CorrelatedIndices(col);
  for (size_t beam=0;beam<2;beam++) {
    for (size_t index=0;index<2;index++) {
      colour = col[beam][index]; 
      if (colour==-1) {
	msg_Error()<<"Error in "<<METHOD<<":\n"
		   <<"   Col["<<beam<<"]["<<index<<"] = -1. "
		   <<"   Will return falseand hope for the best.\n";
	return false;
      }
      p_ladder->GetIn(beam)->SetFlow(index+1,colour);
      parts[1-beam]->SetFlow(index+1,colour);
      if ((parts[1-beam]->Flav().IsGluon() && colour==0) ||
	  (parts[1-beam]->Flav().IsQuark() && 
	   parts[1-beam]->Flav().IsAnti() && 
	   index==1 && colour==0) ||
	  (parts[1-beam]->Flav().IsQuark() && 
	   !parts[1-beam]->Flav().IsAnti() && 
	   index==0 && colour==0))
	msg_Error()<<"Error with colours here (7).\n";
      m_hadrons[1-beam]->GetParticle(counter)->SetFlow(index+1,colour);
    }
  }
  if (!p_ladder->GenerateColourIndices(fix)) return false;

  return true;
}

int Colour_Generator::CorrelatedIndices(int col[2][2]) {
  int fix(-1), colour(-1);
  Flavour flav[2];
  for (size_t beam=0;beam<2;beam++) {
    flav[beam] = p_ladder->GetIn(beam)->m_flav;
    if (flav[beam].IsQuark()) {
      if (flav[beam].IsAnti())       col[beam][0] = 0;
      else if (!flav[beam].IsAnti()) col[beam][1] = 0;
    }
  }
  if (flav[0].IsGluon() && flav[1].IsGluon()) {
    size_t index    = PickIndexAndColour(colour);
    col[0][index]   = col[1][1-index] = colour;
    col[0][1-index] = PickOneColour(0,1-index,colour);
    col[1][index]   = PickOneColour(1,index,colour);
    fix = index==0?1:2;
  }
  else {
    for (size_t beam=0;beam<2;beam++) {
      if (flav[beam].IsQuark() && !flav[beam].IsAnti() && 
	  flav[1-beam].IsGluon()) {
	colour         = PickColourPair(beam,0);
	col[beam][0]   = col[1-beam][1] = colour;
	col[1-beam][0] = PickOneColour(1-beam,0,colour);
	fix = beam==0?1:2;
	break;
      }
    }
  }
  //   msg_Out()<<METHOD<<"("<<flav[0]<<" "<<flav[1]<<"): fix = "<<fix<<".\n";
  return fix;
} 

int Colour_Generator::UncorrelatedIndices(int col[2][2]) {
  int fix(-1), avoid(-1);
  Flavour flav[2];
  for (size_t beam=0;beam<2;beam++) {
    flav[beam] = p_ladder->GetIn(beam)->m_flav;
    if (flav[beam].IsQuark()) {
      if (flav[beam].IsAnti())       col[beam][0] = 0;
      else if (!flav[beam].IsAnti()) col[beam][1] = 0;
    }
  }
  if (flav[0].IsGluon() && flav[1].IsGluon()) {
    fix = ran->Get()>0?1:2;
    for (size_t beam=0;beam<2;beam++) {
      avoid = -1;
      for (size_t index=0;index<2;index++) {
	avoid = col[beam][index] = PickOneColour(beam,index,avoid);
      }
    }
  }
  else {
    for (size_t beam=0;beam<2;beam++) {
      if (flav[beam].IsQuark() && !flav[beam].IsAnti() && 
	  flav[1-beam].IsGluon()) {
	avoid = -1;
	col[beam][0]   = PickOneColour(beam,0,avoid);
	avoid          = col[1-beam][0] = PickOneColour(1-beam,0,avoid);
	col[1-beam][1] = PickOneColour(1-beam,1,avoid);
	fix = beam==0?1:2;
	break;
      }
    }
  }
  return fix;
} 

size_t Colour_Generator::PickIndexAndColour(int & col) {
  msg_Tracking()<<METHOD<<": " 
		<<m_col[0][0].size()<<" "<<m_col[1][1].size()<<" & "
		<<m_col[0][1].size()<<" "<<m_col[1][0].size();
  int max0(Max(m_col[0][0].size(),m_col[1][1].size()));
  int max1(Max(m_col[0][1].size(),m_col[1][0].size()));
  size_t index(-1);
  if (max0>max1)      index = 0;
  else if (max0<max1) index = 1;
  else                index = ran->Get()>0.5?0:1;

  col = PickColourPair(0,index);
  msg_Tracking()<<"    col["<<0<<"]["<<index<<"] ---> "<<col<<".\n";
  return index;
}

int Colour_Generator::
PickOneColour(const size_t & beam,const size_t & index,const int & avoid) {
  msg_Tracking()<<METHOD<<"(beam = "<<beam<<", index = "<<index<<", "
		<<"avoid = "<<avoid<<"): sizes = "
		<<m_col[beam][index].size()<<" "<<m_col[1-beam][1-index].size();
  int col(-1);
  for (set<int>::iterator cit=m_col[beam][index].begin();
       cit!=m_col[beam][index].end();cit++) {
    if ((*cit)!=avoid) {
      col = (*cit);
      m_col[beam][index].erase(col);
      break;
    }
  }
  if (col==-1) {
    Flow Flow;
    col = Flow.Counter();
    m_col[beam][1-index].insert(col);
  }
  msg_Tracking()<<" ---> "<<col<<".\n";
  return col;
}

int Colour_Generator::
PickColourPair(const size_t & beam,const size_t & index) {
  msg_Tracking()<<METHOD<<"(beam = "<<beam<<", index = "<<index<<"): "
		<<m_col[beam][index].size()<<" "<<m_col[1-beam][1-index].size();
  int col(-1);
  for (set<int>::iterator cit1=m_col[beam][index].begin();
       cit1!=m_col[beam][index].end();cit1++) {
    for (set<int>::iterator cit2=m_col[1-beam][1-index].begin();
	 cit2!=m_col[1-beam][1-index].end();cit2++) {
      if ((*cit1)==(*cit2)) {
	col = (*cit1);
	m_col[beam][index].erase(col);
	m_col[1-beam][1-index].erase(col);
	break;
      }
    }
    if (col!=-1) break;
  }
  if (col==-1) {
    Flow Flow;
    col = Flow.Counter();
    m_col[beam][1-index].insert(col);
    m_col[1-beam][index].insert(col);
  }
  msg_Tracking()<<" ---> "<<col<<".\n";
  return col;
}

void Colour_Generator::PickTwoColours(const size_t & beam,int * cols) {
  Flow flow;
  size_t index;
  cols[0] = cols[1] = -1;
  for (index=0;index<2;index++) {
    if (m_col[beam][index].size()==0) cols[index] = flow.Counter();
    else {
      cols[index] = (*m_col[beam][index].begin());
    }
  }
  if (cols[0]==cols[1]) {
    if (m_col[beam][0].size()==1 && m_col[beam][1].size()==1) {
      cols[ran->Get()>0.5?0:1] = flow.Counter();
    }
    else {
      index = m_col[beam][0].size()>m_col[beam][1].size()?0:1;
      set<int>::iterator cit=m_col[beam][index].begin();
      cit++;
      cols[index] = (*cit);
      m_col[beam][index].erase(cols[index]);
    }
  }
  for (index=0;index<2;index++) {
    if (cols[index]==(*m_col[beam][1-index].begin())) {
      m_col[beam][1-index].erase(cols[index]);
    }
  }
  msg_Tracking()<<METHOD<<" yields "<<cols[0]<<" "<<cols[1]<<".\n";
}

void Colour_Generator::FinalColours() {
  //   msg_Out()<<METHOD<<":\n";
  size_t flow(0),col(0);
  Particle * part;
  Flavour flav;
  for (size_t beam=0;beam<2;beam++) {
    size_t length(m_hadrons[1-beam]->Size());
    
    /*      msg_Out()<<"   Trips["<<beam<<"]:";
	    for (set<int>::iterator cit=m_col[beam][0].begin();
	    cit!=m_col[beam][0].end();cit++) msg_Out()<<" "<<(*cit);
	    msg_Out()<<"\n";
	    msg_Out()<<"   Antis["<<beam<<"]:";
	    for (set<int>::iterator cit=m_col[beam][1].begin();
	    cit!=m_col[beam][1].end();cit++) msg_Out()<<" "<<(*cit);
	    msg_Out()<<"\n";*/
    
    for (size_t i=length-2;i<length;i++) {
      part = m_hadrons[1-beam]->GetParticle(i);
      flav = part->Flav();
      if (flav.IsQuark() || flav.IsDiQuark()) {
	if ((flav.IsQuark()   && !flav.IsAnti()) ||
	    (flav.IsDiQuark() && flav.IsAnti()))    flow=0;
	if ((flav.IsQuark()   && flav.IsAnti()) ||
	    (flav.IsDiQuark() && !flav.IsAnti()))   flow=1;
	if (m_col[beam][flow].size()>0) {
	  col = (*m_col[beam][flow].begin());
	  part->SetFlow(flow+1,col);
	  m_col[beam][flow].erase(m_col[beam][flow].begin());
	}
	else {
	  part->SetFlow(flow+1,-1);
	  m_col[beam][1-flow].insert(part->GetFlow(flow+1));	  
	}
      }
      else if (flav.IsGluon()) {
	for (flow=0;flow<2;flow++) {
	  if (m_col[beam][flow].size()>0) {
	    col = (*m_col[beam][flow].begin());
	    part->SetFlow(flow+1,col);
	    m_col[beam][flow].erase(m_col[beam][flow].begin());
	  }
	  else {
	    part->SetFlow(flow+1,-1);
	    m_col[beam][1-flow].insert(part->GetFlow(flow+1));	  
	  }
	}
      }
      //       msg_Out()<<(*part)<<"\n";
    }
  }
}
