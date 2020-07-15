// this file is part of MACMA
// author: cecile.grigne@univ-brest.fr

#include "MACMA/physics/earth.h"
#include "MACMA/physics/continent.h"
#include "MACMA/physics/collision.h"

using namespace std;
// =======================================================
Collision::Collision(Earth* earth, double position) : Marker(earth, position)
{
	_id = Earth::getId("Collision");
}
// =======================================================
Collision::~Collision()
{
}
// =======================================================
bool
Collision::isA(string str)
{
	return (str.compare("Marker") == 0 || 
			str.compare("Collision") == 0);
}
// =======================================================
void
Collision::updatePosition(double dt)
{
	if (!_earth->fixedConfiguration()) 
	{
		if (_continent) 
		{
			double U = _continent->getU();  // in cm/yr
			::changePosition(_position, _earth->cm_to_deg(U) * dt); 
		}
		else 
		{
			cerr << "ERROR: " << getShortId() 
			<< " does not have a continent" << endl;
		}
	}
}
// =======================================================
