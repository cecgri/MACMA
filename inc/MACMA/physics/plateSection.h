// Copyright (C) 2011 ENIB
//   Ecole Nationale d'Ingenieurs de Brest (ENIB)
//   CS 73862 - 29238 BREST Cedex 3 - France
//   Tel: +33(0)298 05 89 89, Fax: +33(0)298 05 89 79, e-mail: combes@enib.fr
//
// This file is part of MACMA.
//
//   MACMA is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as publixsshed by
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

#ifndef PLATE_SECTION_H
#define PLATE_SECTION_H

#include <string>
#include <vector>
#include "MACMA/physics/enum_struct.h"
#include "MACMA/physics/earth.h"
#include "MACMA/physics/continent.h"
#include "MACMA/physics/geoElement.h"
using std::string;
using std::vector;

class Plate;

/*
 ====================================================
 Class PlateSection
 ---------------------------------------------------
 - Represents the plate section between two elements
 - Is on a plate
 ====================================================
 */

class PlateSection
{
public:
    PlateSection(Earth* earth, GeoElement* leftElement, GeoElement* rightElement);
    ~PlateSection();
	
    // --- ID
    inline string getId() { return _id; }
	
    // --- Plate
    inline Plate* getPlate() { return _plate; }
    inline void setPlate(Plate* plate) { _plate = plate; }
	
    // --- Borders
    inline void setRightElement(GeoElement* element) { _rightElement = element; }
    inline GeoElement* getRightElement() { return _rightElement; }
    inline void setLeftElement(GeoElement* element) { _leftElement = element; }
    inline GeoElement* getLeftElement() { return _leftElement; }
	
    // --- continent
    virtual bool isContinental();
    inline void setContinent(Continent* continent) { _continent = continent; }
    inline Continent* getContinent() { return _continent; }
	
    // --- Computing
    virtual double computeLength();
    virtual double computeLm();
	
    void updateQtot();
    void updateBathymetry();
    inline double getQtot() {return _Qtot; }
    inline double getVbelowRidge_HSCM() { return _VbelowRidge_HSCM; }
    inline double getVbelowRidge_PM125() { return _VbelowRidge_PM125; }
    inline double getVbelowRidge_PM95() { return _VbelowRidge_PM95; }
    double computeMeanThickness();
	
    // --- Ages
    virtual void initAges();
    virtual void removeExtraAges();
    virtual void moveAges();
    virtual void createNewOceanAges();
    virtual void createNewContinentAges();
    virtual void removeSubductedAges();
    virtual void correctAges();
    virtual void update_U_for_Ages();
    virtual void effectiveAges();
    virtual vector<Age> splitAges(double, double, Direction);
    inline vector<AgeGroup> getAgeGroups() { return _ageGroups; }
    virtual void setAgeGroups();
	
    // --- Graphics
    inline void setSelected(bool selected) { _selected = selected; }
    inline bool isSelected() { return _selected; }
	
    // --- positions
    virtual bool contains(double position);
	
    inline double getLength() { return _length; }
    inline void setLength(double l) { _length = l; }
	
	
protected:
    // ID
    string _id;
	
    // --- Earth
    Earth* _earth;
	
    // --- Plate
    Plate* _plate;
	
    // --- heat flux
    double _Qtot;
	
    // --- volume below ridge
    double _VbelowRidge_HSCM;
    double _VbelowRidge_PM95;
    double _VbelowRidge_PM125;
	
    // --- continent
    Continent* _continent;
	
    // --- Ages
    vector<AgeGroup> _ageGroups;
	
    // --- Borders
    GeoElement* _leftElement;
    GeoElement* _rightElement;
	
    // --- length (in deg)
    double _length;
	
    // --- Graphics
    bool _selected;
};

#endif

