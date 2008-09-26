#include "Width_Calculator.H"

using namespace MODEL;
using namespace ATOOLS;

double Width_Calculator::Width(Decay_Channel * dec) {
  m_M  = dec->GetDecaying().PSMass();
  for (short int i=0;i<dec->NumberOfDecayProducts();i++)
    m_m[i] = dec->GetDecayProduct(i).PSMass(); 

  if (dec->GetDecaying().IsScalar()) {
    switch (dec->NumberOfDecayProducts()) {
    case 2:
      if (dec->GetDecayProduct(0).IsFermion() &&
	  dec->GetDecayProduct(1).IsFermion()) return SFF(dec);
      if (dec->GetDecayProduct(0).IsScalar() &&
	  dec->GetDecayProduct(1).IsScalar()) return SSS(dec);
      if ((dec->GetDecayProduct(0).IsVector() &&
	   dec->GetDecayProduct(1).IsScalar()) ||
	  (dec->GetDecayProduct(1).IsVector() &&
	   dec->GetDecayProduct(0).IsScalar())) return SVS(dec);
      if (dec->GetDecayProduct(0).IsVector() &&
	  dec->GetDecayProduct(1).IsVector())   return SVV(dec);
    case 3:
    default:
      break;
    }
  }
  else if (dec->GetDecaying().IsVector()) {
    switch (dec->NumberOfDecayProducts()) {
    case 2:
      if (dec->GetDecayProduct(0).IsFermion() &&
	  dec->GetDecayProduct(1).IsFermion())  return VFF(dec);
      if (dec->GetDecayProduct(0).IsScalar() &&
	  dec->GetDecayProduct(1).IsScalar())   return VSS(dec);
      if ((dec->GetDecayProduct(0).IsVector() &&
	   dec->GetDecayProduct(1).IsScalar()) ||
	  (dec->GetDecayProduct(1).IsVector() &&
	   dec->GetDecayProduct(0).IsScalar())) return VVS(dec);
      if (dec->GetDecayProduct(0).IsVector() &&
	  dec->GetDecayProduct(1).IsVector())   return VVV(dec);
    case 3:
    default:
      break;
    }
  }
  else if (dec->GetDecaying().IsFermion()) {
    switch (dec->NumberOfDecayProducts()) {
    case 2:
      if ((dec->GetDecayProduct(0).IsFermion() &&
	   dec->GetDecayProduct(1).IsVector())||
           (dec->GetDecayProduct(1).IsFermion() &&
	    dec->GetDecayProduct(0).IsVector())) return FFV(dec);
      if ((dec->GetDecayProduct(0).IsFermion() &&
	   dec->GetDecayProduct(1).IsScalar())||
          (dec->GetDecayProduct(1).IsFermion() &&
	   dec->GetDecayProduct(0).IsScalar())) return FFS(dec);
    case 3:
    default:
      break;
    }
  }

  msg_Error()<<"Error in "<<METHOD<<":"<<std::endl
	     <<"   Did not find a width calculator for the decay channel below,"<<std::endl
	     <<"   return 0 width."<<std::endl;
  dec->Output();

  return 0;
}

Single_Vertex * Width_Calculator::FindVertex(const Flavour & in,
					     const Flavour & out1,const Flavour & out2) {
  Vertex_List & vertexlist = (*p_vertextable)[in];
  for (size_t i=0;i<vertexlist.size();i++) {
    if ((vertexlist[i]->in[1]==out1 && vertexlist[i]->in[2]==out2) ||
	(vertexlist[i]->in[1]==out2 && vertexlist[i]->in[2]==out1)) {
      return vertexlist[i];
    };
  }
  return NULL;
}

double Width_Calculator::TwoBodyPref(const double M,const double m1,const double m2) {
  return sqrt(Lambda2(M,m1,m2))/(8.*M_PI*sqr(M));
}

double Width_Calculator::Lambda2(const double M,const double m1,const double m2) {
  return ((sqr(M)-sqr(m1+m2))*(sqr(M)-sqr(m1-m2)))/(4.*sqr(M));
}

double Width_Calculator::ColorFactor(Single_Vertex * vertex) 
{
  double colfac(1.);
  if (vertex->Color.size()==1){
    switch  (int(vertex->Color.back().m_type)){
    case (int(cf::D)) : colfac *= 3.; break;
    case (int(cf::T)) : colfac *= 4.; break;
    case (int(cf::F)) : colfac *=24.; break;
    }
  }
  return colfac;
}

double Width_Calculator::Norm(Flavour in)
{
  double norm(1.);
  if (in.IsFermion())                 norm/=2.;
  else if (in.IsVector())             norm/=3.;
  if (abs(in.StrongCharge())==3)      norm/=3.;
  else if (abs(in.StrongCharge())==9) norm/=8.;
  else if (abs(in.StrongCharge())==8) norm/=8.;
  return norm;
}


///////////////////////////////////////////////////////////////////////////////////////////////

double Width_Calculator::SFF(Decay_Channel * dec) {
  Single_Vertex * vertex(FindVertex(dec->GetDecaying(),
				    dec->GetDecayProduct(0),dec->GetDecayProduct(1)));
  if (vertex==NULL) return 0.;
  double M2(sqr(m_M)),m02(sqr(m_m[0])),m12(sqr(m_m[1]));
  double ME2   = (M2-m02-m12)*(std::norm(vertex->cpl[0].Value())+
			       std::norm(vertex->cpl[1].Value()))-
                 2.*m_m[0]*m_m[1]*(2.*std::abs(vertex->cpl[0].Value()*
					       vertex->cpl[1].Value()));
  double width = TwoBodyPref(m_M,m_m[0],m_m[1])*ME2*ColorFactor(vertex)*Norm(dec->GetDecaying());

  return width;
}

double Width_Calculator::SSS(Decay_Channel * dec) {
  Single_Vertex * vertex(FindVertex(dec->GetDecaying(),
				    dec->GetDecayProduct(0),dec->GetDecayProduct(1)));
  if (vertex==NULL) return 0.;
  double ME2   = std::norm(vertex->cpl[0].Value());
  if ((dec->GetDecayProduct(0))==(dec->GetDecayProduct(1)))
    ME2*=1./2; 
  double width = TwoBodyPref(m_M,m_m[0],m_m[1])*ME2*ColorFactor(vertex)*Norm(dec->GetDecaying());

  return width;
}

double Width_Calculator::SVS(Decay_Channel * dec) {
  Single_Vertex * vertex(FindVertex(dec->GetDecaying(),
				    dec->GetDecayProduct(0),dec->GetDecayProduct(1)));
  if (vertex==NULL) return 0.;
  double m02(sqr(m_M)),m12(sqr(m_m[1])),m22(sqr(m_m[0]));

  double ME2   = std::norm(vertex->cpl[0].Value());
  ME2*=(sqr(m02)-2.*(m12+m22)*m02+sqr(m12-m22))/m22;

  double width = TwoBodyPref(m_M,m_m[0],m_m[1])*ME2*ColorFactor(vertex)*Norm(dec->GetDecaying());

  return width;
}

double Width_Calculator::SVV(Decay_Channel * dec) {
  Single_Vertex * vertex(FindVertex(dec->GetDecaying(),
				    dec->GetDecayProduct(0),dec->GetDecayProduct(1)));
  if (vertex==NULL) return 0.;

  double M2(sqr(m_M)),m02(sqr(m_m[0])),m12(sqr(m_m[1]));

  double ME2 = std::norm(vertex->cpl[0].Value())*(2.*sqr(M2)-m02*m12);

  if ((dec->GetDecayProduct(0))==(dec->GetDecayProduct(1)) &&
      dec->GetDecayProduct(0).Strong()==false)
    ME2*=1./2; 
  double width = TwoBodyPref(m_M,m_m[0],m_m[1])*ME2*ColorFactor(vertex)*Norm(dec->GetDecaying());

  return width;
}


double Width_Calculator::VFF(Decay_Channel * dec) {
  Single_Vertex * vertex(FindVertex(dec->GetDecaying(),
				    dec->GetDecayProduct(0),dec->GetDecayProduct(1)));
  if (vertex==NULL) return 0.;

  double ME2   = 
    ((3.*(sqr(m_M)-sqr(m_m[0])-sqr(m_m[1]))-4.*Lambda2(m_M,m_m[0],m_m[1]))*
     (std::norm(vertex->cpl[0].Value())+std::norm(vertex->cpl[1].Value()))) +
    6.*m_m[0]*m_m[1]*(2.*std::abs(vertex->cpl[0].Value()*vertex->cpl[1].Value()));
  double width = TwoBodyPref(m_M,m_m[0],m_m[1])*ME2*ColorFactor(vertex)*Norm(dec->GetDecaying());

  std::cout<<"In "<<METHOD<<" for : "
	   <<dec->GetDecaying()<<" -> "<<dec->GetDecayProduct(0)<<"+"<<dec->GetDecayProduct(1)<<std::endl
	   <<"  = ME^2 = "<<ME2<<" * (colour = "<<ColorFactor(vertex)<<") * norm = "
	   <<Norm(dec->GetDecaying())<<std::endl;

  return width;
}

double Width_Calculator::VSS(Decay_Channel * dec) {
  Single_Vertex * vertex(FindVertex(dec->GetDecaying(),
				    dec->GetDecayProduct(0),dec->GetDecayProduct(1)));
  if (vertex==NULL) return 0.;
  double m02(sqr(m_M)),m12(sqr(m_m[0])),m22(sqr(m_m[1]));
  double ME2   =  std::norm(vertex->cpl[0].Value());
  ME2*=(sqr(m02)-2.*(m12+m22)*m02+sqr(m12-m22))/m02;

  if ((dec->GetDecayProduct(0))==(dec->GetDecayProduct(1))) ME2*=1./2;
  
  double width = TwoBodyPref(m_M,m_m[0],m_m[1])*ME2*ColorFactor(vertex)*Norm(dec->GetDecaying());

  return width;
}

double Width_Calculator::VVS(Decay_Channel * dec) {
  Single_Vertex * vertex(FindVertex(dec->GetDecaying(),
				    dec->GetDecayProduct(0),dec->GetDecayProduct(1)));
  if (vertex==NULL) return 0.;
  double M2(sqr(m_M)),m02(sqr(m_m[0])),m12(sqr(m_m[1]));
  double ME2   = (std::norm(vertex->cpl[0].Value()))*
                 (sqr(M2+m02-m12)+8.*m02*M2)/(4.*m02*M2);
  double width = TwoBodyPref(m_M,m_m[0],m_m[1])*ME2*ColorFactor(vertex)*Norm(dec->GetDecaying());

  return width;
}

double Width_Calculator::VVV(Decay_Channel * dec) {
  Single_Vertex * vertex(FindVertex(dec->GetDecaying(),
				    dec->GetDecayProduct(0),dec->GetDecayProduct(1)));
  if (vertex==NULL) return 0.;
  double m02(sqr(m_M)),m12(sqr(m_m[0])),m22(sqr(m_m[1]));
  double ME2  = (std::norm(vertex->cpl[0].Value()))*
     ((sqr(m02) - 2.*m02*m12 + sqr(m12) - 2.*m02*m22 - 2.*m12*m22 + sqr(m22))*
      (sqr(m02) + 10.*m02*m12 + sqr(m12) + 10.*m02*m22 + 10.*m12*m22 + sqr(m22)))/
      (4.*m02*m12*m22);

  if ((dec->GetDecayProduct(0))==(dec->GetDecayProduct(1))) ME2*=1./2;
 
  double width = TwoBodyPref(m_M,m_m[0],m_m[1])*ME2*ColorFactor(vertex)*Norm(dec->GetDecaying());

  return width;
}

double Width_Calculator::FFS(Decay_Channel * dec) {
  Single_Vertex * vertex(FindVertex(dec->GetDecaying(),
				    dec->GetDecayProduct(0),dec->GetDecayProduct(1)));
  if (vertex==NULL) return 0.;

  double ME2   = (2.*m_M*m_m[1]*(2.*std::abs(vertex->cpl[0].Value()*vertex->cpl[1].Value())) +
                 (sqr(m_M)+sqr(m_m[1])-sqr(m_m[0]))*
                 (std::norm(vertex->cpl[0].Value())+std::norm(vertex->cpl[1].Value())));
  double width = TwoBodyPref(m_M,m_m[1],m_m[0])*ME2*ColorFactor(vertex)*Norm(dec->GetDecaying());

  return width;
}

double Width_Calculator::FFV(Decay_Channel * dec) {
  Single_Vertex * vertex(FindVertex(dec->GetDecaying(),
				    dec->GetDecayProduct(0),dec->GetDecayProduct(1)));
  if (vertex==NULL) return 0.;
  double M2(sqr(m_M)),m02(sqr(m_m[1])),m12(sqr(m_m[0]));

  double ME2   = (M2+m02-m12 +(M2 -m02-m12)*(M2-m02+m12)/m12)*
                 (std::norm(vertex->cpl[0].Value())+std::norm(vertex->cpl[1].Value()))-
                 6.*m_M*m_m[0]*(2.*std::abs(vertex->cpl[0].Value()*vertex->cpl[1].Value())) ;
  double width = TwoBodyPref(m_M,m_m[1],m_m[0])*ME2*ColorFactor(vertex)*Norm(dec->GetDecaying());

  return width;
}
