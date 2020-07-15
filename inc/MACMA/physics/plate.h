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

#ifndef PLATE_H_
#define PLATE_H_

#include <vector>
#include <fstream>
using std::vector;

class Earth;
class GeoElement;
class PlateSection;
class Cell;
class Continent;
class Force;

/*
==================================================
 Class Plate
 -------------------------------------------------
 - Represents the lithosphere between two interfaces
 - Is composed of plate sections
==================================================
*/
class Plate
{
 public:
  Plate(Earth* earth, vector<PlateSection*> sections);
  ~Plate();
  
  // --- ID
  inline string getId() { return _id; } // useless for plates

  // --- Computing
  virtual double computeLength(); // in degrees
  virtual double computeLm(); // in meters
  virtual double computeHorizontalViscosity(); 
  virtual double computeContinentalLength(); // in degrees

  // --- Sections
  virtual PlateSection* findSection(const char* name);
  virtual int findSection(PlateSection* section);
  inline vector<PlateSection*>& accessSections(){ return _sections; }
  void addSection(PlateSection*);

  // --- Cells
  inline vector<Cell*>& accessCells() {return _cells; }
  virtual void setCells();

  // --- Continents
  vector<Continent*> getContinents();

  // --- Interfaces
  virtual GeoElement* getLeftElement();
  virtual GeoElement* getRightElement();
  
  // --- Plates
  virtual Plate* getLeftPlate();
  virtual Plate* getRightPlate();

  // --- Velocity
  inline void setU(double v) {_U = v; }
  inline double getU() { return _U; } // in cm/yr

  // --- Position
  double getMiddlePosition();

  // --- Continent
  bool isContinental();
  bool isOceanic();
  
  // --- Forces
  void updateForces(); // compute all that's needed to get plate velocities
  void completeForces(); /* add Mantle Drag, Bending and Viscous Shear computation
			    (for outputs) */
  inline Force* & accessForces() { return _forces; }
  inline void setForces(Force* force) { _forces = force ;}

  bool isTensive(); /* check if difference of forces left-right 
		       make tensive or compressive plate */


  // -- properties of the plate
  inline bool isSubducting() { return _subducting; }

  inline void setType(string str) { _type = str; }
  inline string getType() { return _type; }
  inline bool hasType(string str) { return _type == str; }
  void updateType(); 

  // --- Graphics
  inline void setSelected(bool selected) { _selected = selected; }
  inline bool isSelected() { return _selected; }
  
 protected:
  // --- ID
  string _id; // not a real id for plates
  // --- Earth
  Earth* _earth;
  
  // --- Sections
  vector<PlateSection*> _sections;
  
  // --- Cells
  vector<Cell*> _cells;

  // --- Values
  double _U;		// cm/yr
  
  // --- Graphics
  bool _selected;
  
  // --- Forces
  Force* _forces;

  // --- others
  string _type; // SU, RU, RS, UU, RR or SS
  bool _subducting; // defined in updateType()
};

#endif
