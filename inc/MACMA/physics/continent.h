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

#ifndef CONTINENT_H_
#define CONTINENT_H_

#include <vector>
#include <string>
#include <fstream>
using std::vector;
using std::string;
using std::ofstream;
#include "MACMA/physics/enum_struct.h"

// forward declarations
class Earth;
class Plate;
class PlateSection;
class LeftContinentExtremity;
class RightContinentExtremity;
class WarmingZone;
class Age;

// ==================================================
class Continent
{
public:
	Continent(Earth* earth, double position, double length, string origin);
	~Continent();
	
	inline string getId() { return _id; }
	
	virtual double getLeftPosition();
	virtual double getRightPosition();
	virtual double getMiddlePosition();
	
	void updatePosition(double dt);
	
	// --- Getters / Setters
	virtual void setPosition(double position);
	inline double getPosition() { return _position; }
	virtual void setLength(double length);
	inline double getLength() { return _length; }
	
	virtual double getLm(); // length in m
	virtual double getU(); // in cm/yr
	
	virtual string getOrigin() { return _origin; }
	
	inline void setPlateSection(PlateSection* section) { _plateSection = section; }
	inline PlateSection* getPlateSection() { return _plateSection; }
	virtual Plate* getPlate();
	
	// extremities 
	inline void setLeftExtremity(LeftContinentExtremity* extremity) { _leftExtremity = extremity; }
	inline LeftContinentExtremity* getLeftExtremity() { return _leftExtremity; }
	inline void setRightExtremity(RightContinentExtremity* extremity) { _rightExtremity = extremity; }
	inline RightContinentExtremity* getRightExtremity() { return _rightExtremity; }
	
	// warming zone and break up
	inline WarmingZone* getWarmingZone() { return _warmingZone; }
	inline void setWarmingZone(WarmingZone* wz) { _warmingZone = wz; }
	void updateBreakable();
	bool checkOpening();
	inline double getBirthDate() { return _birthDate; }
	inline bool canBreak() { return _breakable; }
	
	// ages
	inline vector<AgeGroup> getAgeGroups() { return _ageGroups; }
	void shareNeighborAges();
	void gatherAges(vector<Age>, vector<Age>);
	void setAgeGroups();
	
	// others
	unsigned int getCounter();
	
	// files
	void updateFile();
	
	// --- Graphics
	inline void setSelected(bool selected) { _selected = selected; }
	inline bool isSelected() { return _selected; }
	
protected:
	string _id;
	
	Earth* _earth;
	double _position;
	double _length;
	PlateSection* _plateSection;
	
	LeftContinentExtremity* _leftExtremity;
	RightContinentExtremity* _rightExtremity;
	
	WarmingZone* _warmingZone;
	bool _breakable;
	double _birthDate;
	string _origin;
	
	ofstream _file;
	
	vector<AgeGroup> _ageGroups;
	
	// --- Graphics
	bool _selected;
};

#endif

