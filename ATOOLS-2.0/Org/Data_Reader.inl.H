#ifndef Data_Reader_inl_H
#define Data_Reader_inl_H

namespace ATOOLS {

  inline void Data_Reader::SetTags(const std::map<std::string,std::string> &tags) 
  { 
    m_tags=tags; 
  }
  
  inline void Data_Reader::AddTags(const std::map<std::string,std::string> &tags) 
  { 
    for (std::map<std::string,std::string>::const_iterator tit=tags.begin();
	 tit!=tags.end();++tit) m_tags[tit->first]=tit->second; 
  }
  
  inline void Data_Reader::SetTag(const std::string &tag,const std::string &value) 
  { 
    m_tags.clear(); 
    m_tags[tag]=value; 
  }
  
  inline void Data_Reader::AddTag(const std::string &tag,const std::string &value) 
  { 
    m_tags[tag]=value; 
  }

  inline const std::string &Data_Reader::Tag(const std::string &tag) const 
  { 
    return m_tags.find(tag)->second;
  }
  
  inline const std::map<std::string,std::string> &Data_Reader::Tags() const
  { 
    return m_tags;
  }

}// end of namespace ATOOLS

#endif
