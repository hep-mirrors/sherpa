#ifndef ResonanceFlavour_h
#define ResonanceFlavour_h

namespace HADRONS {

  class ResonanceFlavour {

	std::string m_name;
	double m_mass, m_width;

	public:
	ResonanceFlavour() {};
	ResonanceFlavour( std::string _name, double _mass, double width );
	ResonanceFlavour( ATOOLS::kf::code _kfc, double _mass, double width );
	~ResonanceFlavour() {}

	void Set( ATOOLS::kf::code _kfc, double _mass, double width );
	 
	double Mass() { return m_mass; }
	double Width() { return m_width; }
	std::string Name() { return m_name; }

  };

}

#endif
  
