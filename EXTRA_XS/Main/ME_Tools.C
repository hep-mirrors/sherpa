#include "EXTRA_XS/Main/ME_Tools.H"

#include "ATOOLS/Org/Exception.H"

using namespace METOOLS;

std::vector<Vertex*> EXTRAXS::ConstructVertices(Current* cur1,
						Current* cur2,
						Current* prop)
{
  std::vector<Vertex*> ret;
  Current_Vector curs(2);
  curs[0]=cur1;
  curs[1]=cur2;

  // try first rotation
  Vertex_Key *vkey(Vertex_Key::New(curs,prop,MODEL::s_model));
  MODEL::VMIterator_Pair keyrange(MODEL::s_model->GetVertex(vkey->ID()));
  for (MODEL::Vertex_Map::const_iterator it=keyrange.first; it!=keyrange.second; ++it) {
    vkey->p_mv=it->second;
    vkey->m_p=std::string(1,'D');
    ret.push_back(new Vertex(*vkey));
    ret.back()->AddJ(vkey->m_j);
    ret.back()->SetJC(prop);
  }
  vkey->Delete();

  //try second rotation
  std::swap<Current*>(curs[0],curs[1]);
  vkey=Vertex_Key::New(curs,prop,MODEL::s_model);
  keyrange=MODEL::s_model->GetVertex(vkey->ID());
  for (MODEL::Vertex_Map::const_iterator it=keyrange.first; it!=keyrange.second; ++it) {
    vkey->p_mv=it->second;//fixme!!
    vkey->m_p=std::string(1,'D');
    ret.push_back(new Vertex(*vkey));
    ret.back()->AddJ(vkey->m_j);
    ret.back()->SetJC(prop);
  }
  vkey->Delete();

  if (ret.size()==0) THROW(fatal_error, "vertex not found: "+vkey->ID());

  return ret;
}

