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

#ifndef LEFT_CONTINENT_EXTREMITY_H_
#define LEFT_CONTINENT_EXTREMITY_H_

#include "MACMA/physics/continentExtremity.h"

class Earth;
class Continent;
/*
==================================================
 Class LeftContinentExtremity : extends ContinentExtremity
 -------------------------------------------------
==================================================
*/

class LeftContinentExtremity : public ContinentExtremity
{
 public:
  LeftContinentExtremity(Earth* earth, Continent* continent);
  ~LeftContinentExtremity();
  
  // --- GeoElement
  //		virtual void updateAges(double dt);
  inline bool isA(string className) { return className == "LeftContinentExtremity" || 
      className == "ContinentExtremity" || 
      className == "GeoElement"; }

  inline string getClass() { return "LeftContinentExtremity"; }

  void modifyVelocity();
  
 protected:
  
};

#endif
