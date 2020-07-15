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

#ifndef GEOELEMENT_H_
#define GEOELEMENT_H_

#include <vector>
#include <string>
using std::vector;
using std::string;
#include "MACMA/physics/enum_struct.h"
#include "MACMA/physics/age.h"

class Earth;
class PlateSection;

/*
 ==================================================
 Class GeoElement : Abstract class
 -------------------------------------------------
 - Unify interfaces, staples, continent extremity
 ==================================================
 */

class GeoElement
{
public:
	virtual ~GeoElement();
	
	// --- ID
    inline string getId() const { return _id; }
    inline void setId(string id) { _id = id; }
    virtual string getClass();
    virtual string getShortName();
    virtual string getShortId(); //for shorter outputs
	
    // --- Position
    inline double getPosition() const { return _position; }
    inline void setPosition(double position) { _position = position; }
    virtual void updatePosition(double dt);
    virtual double getNextPosition(double);
	
    // --- Velocity
    inline double getU() const { return _U; }
    inline void setU(double u) { _U = u; }
	
    // --- Ages
    inline void setRightAge(Age age) { _rightAge = age; }
    inline void setLeftAge(Age age) { _leftAge = age; }
    virtual void setRightAge(Earth*, double, double);
    virtual void setLeftAge(Earth*, double, double);
	
    inline Age& getRightAge() { return _rightAge; }
    inline Age& getLeftAge() { return _leftAge; }
    inline vector<Age>& getRightAges() { return _rightAges; }
    inline vector<Age>& getLeftAges() { return _leftAges; }
	
	inline void setRightAges(vector<Age> ages) { _rightAges = ages; }
	inline void setLeftAges(vector<Age> ages) { _leftAges = ages; }	
    
	virtual void clearRightAges();
	virtual void clearLeftAges();
	
	virtual bool sectionAgesIsEmpty(Direction);
	
	// thickness of plate: related to age
	virtual double computePlateThickness(Direction); 
	
	
	// --- position
	virtual bool isStrictlyWithin(double rightPosition, 
								  double leftPosition);
	virtual bool isWithinOrEqual(double rightPosition, 
								 double leftPosition,
								 Direction dir);
	
	// --- Sections
	inline void setLeftSection(PlateSection* section) { _leftSection = section; }
	inline PlateSection* getLeftSection() { return _leftSection; }
	inline void setRightSection(PlateSection* section) { _rightSection = section; }
	inline PlateSection* getRightSection() { return _rightSection; }
	
	// --- To define
	virtual bool isA(string className);
	virtual void updateVelocity() = 0;
	
	// --- Graphics
	inline void setSelected(bool selected) { _selected = selected; }
	inline bool isSelected() { return _selected; }
	
protected:
	GeoElement(Earth* earth, double position);
	
	// --- ID
	string _id;
	
	// --- Earth
	Earth* _earth;
	
	// --- Position
	double _position;
	
	// --- Velocity
	double _U; // cm/yr
	
	// --- Ages
	vector<Age> _leftAges;
	vector<Age> _rightAges;
	Age _leftAge;
	Age _rightAge;
	
	// --- Velocities 
	double _leftU;
	double _rightU;
	
	// --- Sections
	PlateSection* _leftSection;
	PlateSection* _rightSection;
	
	// --- Graphics
	bool _selected;
	
};

#endif

