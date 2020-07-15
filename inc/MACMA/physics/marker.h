// this file is part of MACMA
// author: cecile.grigne@univ-brest.fr

#ifndef MARKER_H_
#define MARKER_H_

#include <string>
using std::string;

class Earth;
class Continent;

/* ==================================
 Markers for following special positions.
 -- Abstract class	--
 ================================== */
class Marker
{
public:  
	virtual ~Marker();
	
	inline string getId() { return _id; }
	
	inline double getPosition() { return _position; }
	inline void setPosition(double pos) { _position = pos; }
	
	inline double getInitAge() { return _initAge; }
	inline void setInitAge(double age) { _initAge = age; }
	
	inline void setContinent(Continent* continent) { _continent = continent; }
	inline Continent* getContinent() { return _continent; }
	
	virtual	double getAgeMa();
	virtual bool isA(string);
	
	virtual void updatePosition(double) = 0;	
	
	virtual string getShortId();
	
protected:
	Marker(Earth* earth, double position);
	
	string _id;
	Earth* _earth;
	double _position;
	double _initAge;
	Continent* _continent;
};

#endif
