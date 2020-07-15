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

#ifndef STAPLE_H_
#define STAPLE_H_

#include <string>
#include "MACMA/physics/geoElement.h"
using std::string;

class Earth;
class Plate;
class PlateSection;
class RightContinentExtremity;
class LeftContinentExtremity;

/*
==================================================
 Class Staple : extends GeoElement
 -------------------------------------------------
 - Has left and right plate sections
 - Velocity law : U = U_plate
==================================================
*/

class Staple : public GeoElement
{
 public:
  Staple(Earth* earth, double position);
  ~Staple();
		
  // --- GeoElement
  virtual void updateVelocity();
  virtual void updatePosition(double dt);
  virtual double getNextPosition(double dt);
  virtual bool isA(string className);
  virtual string getClass();

  // --- Continents
  inline void setLeftNeighborExtremity(RightContinentExtremity* extremity) { _leftNeighborExtremity = extremity; }
  inline RightContinentExtremity* getLeftNeighborExtremity() { return _leftNeighborExtremity; }
  inline void setRightNeighborExtremity(LeftContinentExtremity* extremity) { _rightNeighborExtremity = extremity; }
  inline LeftContinentExtremity* getRightNeighborExtremity() { return _rightNeighborExtremity; }

 protected:
  RightContinentExtremity* _leftNeighborExtremity;
  LeftContinentExtremity* _rightNeighborExtremity;
};

#endif
