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

#ifndef CELL_H_
#define CELL_H_

#include <string>
using std::string;

// forward declarations
class Earth;
class Plate;
class GeoElement;

/*
==================================================
 Class Cell
 -------------------------------------------------
 - Is under a plate  -> just for drawing
 ==================================================
*/

class Cell
{
 public:
  Cell(Earth* earth);
  ~Cell();
  
  // --- ID
  inline string getId() { return _id; }
  
  // --- Orientation
  void updateOrientation();
  inline double getOrientation() { return _orientation; }
  inline void setOrientation(double orient) {_orientation = orient; }


  // --- Elements
  inline void setLeftElement(GeoElement* element) { _leftElement = element; }
  inline void setRightElement(GeoElement* element) { _rightElement = element; }
  inline GeoElement* getLeftElement() { return _leftElement; }
  inline GeoElement* getRightElement() { return _rightElement; }

  // --- Position
  inline void setLeftPosition(double left) {_leftPosition = left; }
  inline void setRightPosition(double right) {_rightPosition = right; }
  inline double getLeftPosition() {return _leftPosition; }
  inline double getRightPosition() {return _rightPosition; }

  // --- Plate
  inline void setPlate(Plate* plate) {_plate = plate; }
  inline Plate* getPlate() {return _plate; }
  
  // --- Graphics
  inline void setSelected(bool selected) { _selected = selected; }
  inline bool isSelected() { return _selected; }
  virtual double getStartAngle();
  virtual double getStopAngle();
  double getSpeedFactor();

 protected:
  string _id;
  
  // --- Earth
  Earth* _earth;
  
  double _orientation;
  
  // --- Elements
  GeoElement* _leftElement;
  GeoElement* _rightElement;
  
  // --- Positions
  double _leftPosition;
  double _rightPosition;

  // --- Plate
  Plate* _plate;
  
  // --- Graphics
  bool _selected;
};

#endif
