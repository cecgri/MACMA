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

#ifndef INTERFACE_H_
#define INTERFACE_H_

#include <string>
#include <fstream>
#include "MACMA/physics/enum_struct.h"
#include "MACMA/physics/geoElement.h"

using std::string;
using std::ofstream;

class Earth;
class Plate;

/*
  ==================================================
  Class Interface : extends GeoElement
  -------------------------------------------------
  - Unify ridges and subductions
  - Has left and right plates
  ==================================================
*/

class Interface : public GeoElement
{
 public:
  virtual ~Interface();
  
  // --- GeoElement
  virtual bool isA(string className);
  virtual string getClass();
  virtual void updateVelocity() = 0;
  
  // --- Interface
  inline Plate* getLeftPlate() { return _leftPlate; }
  inline void setLeftPlate(Plate* plate) { _leftPlate = plate; }
  inline Plate* getRightPlate() { return _rightPlate; }
  inline void setRightPlate(Plate* plate) { _rightPlate = plate; }

  Direction getDirectionDrivingForces(Direction);
  Interface* findNextInterface(Direction);
  bool nearInterface(double);

 protected:
  Interface(Earth* earth, double position);
  
  Plate* _leftPlate;
  Plate* _rightPlate;

  ofstream _file;
};

#endif

