// this file is part of MACMA
// author: cecile.grigne@univ-brest.fr

#ifndef COLLISION_H_
#define COLLISION_H_

#include <string>
#include "MACMA/physics/marker.h"

using std::string;

class Continent;

/* 
 class Collision : extends marker
 */

class Collision : public Marker
{
public:
	Collision(Earth* earth,
			  double position);
	
	~Collision();
	
	bool isA(string);
	void updatePosition(double);
	
	inline string getClass() { return "Collision"; }	
};


#endif

