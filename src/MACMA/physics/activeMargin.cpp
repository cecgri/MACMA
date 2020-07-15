// this file is part of MACMA
// author: cecile.grigne@univ-brest.fr

#include "MACMA/physics/earth.h"
#include "MACMA/physics/activeMargin.h"

using namespace std;
// =======================================================
ActiveMargin::ActiveMargin(Earth* earth, double position) : Marker(earth, position)
{
    _id = Earth::getId("ActiveMargin");
	_subduction = NULL;
    _extremity = NULL;
}
// =======================================================
ActiveMargin::~ActiveMargin()
{
}
// =======================================================
bool
ActiveMargin::isA(string str)
{
    return (str.compare("Marker") == 0 ||
            str.compare("ActiveMargin") == 0);
}
// =======================================================
void
ActiveMargin::updatePosition(double dt)
{
    if(!_earth->fixedConfiguration())
    {
        if (_subduction)
		{
			// given exactly the position of the subduction while it's active
			_position = _subduction->getPosition();
		}
		else if (_continent)
		{
			// following the accretion position within the continent when margin not active anymore
			double U = _continent->getU();  // in cm/yr
			::changePosition(_position, _earth->cm_to_deg(U) * dt); 
		}
        else
            cerr << "ERROR: " << getShortId()
            << " does not have a subduction" << endl;
    }
}
// =======================================================
