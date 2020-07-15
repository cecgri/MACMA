// this file is part of MACMA
// author: cecile.grigne@univ-brest.fr

#include "MACMA/physics/earth.h"
#include "MACMA/physics/marker.h"

using namespace std;
// =======================================================
Marker::Marker(Earth* earth, 
			   double position) : 
_earth(earth),
_position(position)
{
	_initAge = _earth->getAgeMa();
}
// =======================================================
Marker::~Marker()
{
}
// =======================================================
bool
Marker::isA(string str)
{
	return (str.compare("Marker") == 0);
}
// =======================================================
string
Marker::getShortId()
{
	string shortId = _id;
	if (this->isA("Collision")) 
	{
		string str = "CO";
		shortId.replace(0, 9, str);
	}
    else if (this->isA("ActiveMargin"))
    {
        string str = "AM";
        shortId.replace(0, 12, str);
    }
	
	return shortId;
}
// =======================================================
double
Marker::getAgeMa()
{
	return _initAge - _earth->getAgeMa();
}
// =======================================================
