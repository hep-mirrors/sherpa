#include "Width_Calculator.H"

using namespace MODEL;
using namespace ATOOLS;

double Width_Calculator::Width(Decay_Channel * dec) {
  if (dec->GetDecaying().IsVector()) {
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
    case 3:
    default:
      msg_Error()<<"Error in "<<METHOD<<":"<<std::endl
		 <<"   Did not find a width calculator for the decay channel below,"<<std::endl
		 <<"   return 0 width."<<std::endl;
      dec->Output();
    }
  }
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

double Width_Calculator::VFF(Decay_Channel * dec) {
  Single_Vertex * vertex(FindVertex(dec->GetDecaying(),dec->GetDecayProduct(0),dec->GetDecayProduct(1)));
  if (vertex==NULL) return 0.;

  double M(dec->GetDecaying().PSMass());
  double m1(dec->GetDecayProduct(0).PSMass()), m2(dec->GetDecayProduct(1).PSMass());

  double ME2   = 
    ((sqr(M)-sqr(m1)-sqr(m2)+4.*Lambda2(M,m1,m2))*
     (std::norm(vertex->cpl[0].Value())+std::norm(vertex->cpl[1].Value()))) -
    6.*m1*m2*(2.*std::abs(vertex->cpl[0].Value()*vertex->cpl[1].Value()));
  double width = TwoBodyPref(M,m1,m2)/3.*ME2;

  if (vertex->ncf==1 && vertex->Color->m_type==cf::D) { width*=3.; }
  return width;
}

double Width_Calculator::VSS(Decay_Channel * dec) {
  return 0;
}

double Width_Calculator::VVS(Decay_Channel * dec) {
  return 0;
}

