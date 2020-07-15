// this file is part of MACMA
// author: cecile.grigne@univ-brest.fr

#include "MACMA/physics/earth.h"
#include "MACMA/physics/age.h"

using namespace std;

// --------------------------------------------------
Age::Age(Earth* earth, double value, double position) : 
_earth(earth),
_value(value),
_effValue(value),
_position(position)
{
	_ocean = true;
	_extra = false; 
	_originContCounter = 0;
}
// --------------------------------------------------
Age::~Age()
{
}
// --------------------------------------------------
void
Age::reset(double pos)
{
	_position = pos;
	_value = 0.0;
	_effValue = 0.0;
	_ocean = true;
	_extra = false;
}
// --------------------------------------------------
void
Age::move(double displacement)
{
	::changePosition(_position, displacement);
}
// --------------------------------------------------
void
Age::addValue(double add)
{
	_value += add;
}
// --------------------------------------------------
