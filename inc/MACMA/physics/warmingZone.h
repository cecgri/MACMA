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

#ifndef WARMING_ZONE_H_
#define WARMING_ZONE_H_

#include <string>
#include <fstream>
#include <math.h>
using std::string;
using std::ofstream;

class Earth;
class Continent;

/*
==================================================
 Class WarmingZone
 -------------------------------------------------
==================================================
*/

class WarmingZone
{
 public:
  WarmingZone(Earth* earth, Continent* continent);
  ~WarmingZone();

  // --- ID
  inline string getId() { return _id; }

  // --- Computing
  double computeThermalPressure();
  double computeTmax();
  double computeWarmingTimeMyr();
  void updateWarming();
  double computeU();
  void updateF();
  

  inline double getLifetime() { return _lifetime; }
  virtual void updateLifetime(double dt);
		
//  inline void setWarming(double warming) { _warming = warming; }
  inline double getWarming() { return _warming; }
	
  inline double getF() { return _F; }
	
  virtual void computeLifetime0(double warming);

 protected:
  string _id;

  // --- Earth
  Earth* _earth;

  // --- Continent
  Continent* _continent;

  double _lifetime; // Myr

  // --- warming
  double _thermalPressure;
  double _warming;
  double _F;

};

#endif
