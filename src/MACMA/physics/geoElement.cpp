// Copyright (C) 2011 ENIB
//   Ecole Nationale d'Ingenieurs de Brest (ENIB)
//   CS 73862 - 29238 BREST Cedex 3 - France
//   Tel: +33(0)298 05 89 89, Fax: +33(0)298 05 89 79, e-mail: combes@enib.fr
//
// This file is part of MACMA.
//
//   MACMA is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
//
//   MACMA is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with MACMA.  If not, see <http://www.gnu.org/licenses/>.
//
//   author: leyaouanq@cervval.com

#include "MACMA/physics/enum_struct.h"
#include "MACMA/physics/geoElement.h"
#include "MACMA/physics/earth.h"

using namespace std;
// ==================================================
GeoElement::GeoElement(Earth* earth, double position) :
_earth(earth),
_position(position),
_rightAge(earth, 0.0, position),
_leftAge(earth, 0.0, position)
{
	_U = 0.0;
	_selected = false;
	
	_leftSection = NULL;
	_rightSection = NULL;
}
// -------------------------------------------------------
GeoElement::~GeoElement()
{
	_leftAges.clear();
	_rightAges.clear();
}
// -------------------------------------------------------
bool
GeoElement::isA(string className)
{
	return className == "GeoElement";
}
// -------------------------------------------------------
string
GeoElement::getClass()
{
	return "GeoElement";
}
// -------------------------------------------------------
void
GeoElement::updatePosition(double dt)
{
	if(!_earth->fixedConfiguration())
		_position = getNextPosition(dt);
}
// -------------------------------------------------------
double 
GeoElement::getNextPosition(double dt)
{
	return ::addPosition(_position, _earth->cm_to_deg(_U) * dt);
}
// -------------------------------------------------------
double
GeoElement::computePlateThickness(Direction direction)
{ /* computes plate thickness on the right 
   or left of the element */
	double age;
	
	switch(direction)
    {
		case LEFT :
		{
			age = _leftAge.getEffValue();
			break;
		}
		case RIGHT :
		{
			age = _rightAge.getEffValue();
			break;
		}
		default : 
		{
			cerr << "ERROR: computePlateThickness(dir) used with no direction dir given"
			<< endl;
			cerr << _id
			<< "\t" << _position << endl;
			cerr << direction << endl;
			exit(EXIT_FAILURE);
		}
    }
	
	double H = Earth::coeff_plate_thickness * 
    sqrt(_earth->get_kappa() * ::Myr_to_sec(age));
	
	return H;
}
// --------------------------------------------------
bool
GeoElement::isStrictlyWithin(double rightPos, double leftPos)
{ /* position of element is within rightPos and leftPos
   ( > to rightPos and < to leftPos )
   */
	return ::isStrictlyWithin( _position,
							  rightPos,
							  leftPos );
}
// --------------------------------------------------
bool
GeoElement::isWithinOrEqual(double rightpos, double leftpos, Direction dir)
{
	return ::isWithinOrEqual( _position,
							 rightpos,
							 leftpos,
							 dir );
}
// --------------------------------------------------
string
GeoElement::getShortName()
{
	string shortName;
	if(this->isA("RightSubduction"))
		shortName = "RS";
	else if(this->isA("LeftSubduction"))
		shortName = "LS";
	else if(this->isA("Ridge"))
		shortName = "R";
	else if(this->isA("Staple"))
		shortName = "S";
	else if(this->isA("RightContinentExtremity"))
		shortName = "RCE";
	else if(this->isA("LeftContinentExtremity"))
		shortName = "LCE";
	
	return shortName;
}
// --------------------------------------------------
string
GeoElement::getShortId()
{
	string shortName = this->getShortName();
	unsigned int length = _id.find(".");
	string shortId = _id;
	shortId.replace(0, length, shortName);
	
	return shortId;
}
// --------------------------------------------------
bool
GeoElement::sectionAgesIsEmpty(Direction direction)
{
	switch(direction)
    {
		case RIGHT : 
			return _rightAges.empty();
			break;
			
		case LEFT :
			return _leftAges.empty();
			break;
			
		default :
			cerr << "Error: GeoElement::sectionAgesIsEmpty(dir) is used without correct dir" << endl;
			exit(EXIT_FAILURE);
    }
}
// --------------------------------------------------
void
GeoElement::setRightAge(Earth* earth, double value, double pos)
{	// used for init in macma.cpp
	_rightAge = Age(earth, value, pos);
}
// --------------------------------------------------
void
GeoElement::setLeftAge(Earth* earth, double value, double pos)
{
	_leftAge = Age(earth, value, pos);
}
// --------------------------------------------------
void
GeoElement::clearRightAges()
{
	_rightAges.clear();
    _rightAge.setPosition(_position);
    _rightAge.setValue(0.0);
    _rightAge.setEffValue(0.0);
    _rightAge.setOceanic(true); // should be by default
}
// --------------------------------------------------
void
GeoElement::clearLeftAges()
{
	_leftAges.clear();
    _leftAge.setPosition(_position);
    _leftAge.setValue(0.0);
    _leftAge.setEffValue(0.0);
    _leftAge.setOceanic(true); // should be by default
}
// --------------------------------------------------

