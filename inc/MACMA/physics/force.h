// This file is part of MACMA
//
// author: cecile.grigne@univ-brest.fr
//

#ifndef FORCE_H_
#define FORCE_H_

#include "MACMA/physics/enum_struct.h"

class Earth;
class GeoElement;

class Force
{
 public:
  Force(Earth*);
  ~Force();

  //  gets ////////////////////////////
  inline double getRidgePush() { return _ridgePush; }
  inline double getSlabPull() { return _slabPull; }
  inline double getSlabSuction() { return _slabSuction; }
  inline double getEtaH() { return _etaH; }
  inline double getEtaVRight() { return _etaVRight; }
  inline double getEtaVLeft() { return _etaVLeft; }
  inline double getEtaBRight() { return _etaBRight; }
  inline double getEtaBLeft() { return _etaBLeft; }
  inline double getMantleDrag() { return _mantleDrag ; }
  inline double getViscousShear() { return _viscousShear ;}
  inline double getBending() { return _bending ;}
  /*     sum of forces or viscosities*/
  inline double getDrivingForces() { return _slabPull + _slabSuction + _ridgePush; }
  inline double getEtaLeft() { return _etaVLeft + _etaBLeft; }
  inline double getEtaRight() { return _etaVRight + _etaBRight; }

  // sets /////////////////////////////
  inline void setRidgePush(double r) { _ridgePush = r; }
  inline void setSlabPull(double s) { _slabPull = s; }
  inline void setSlabSuction(double s) { _slabSuction = s; }
  inline void setEtaH(double e) { _etaH = e; }
  inline void setEtaVRight(double e) { _etaVRight = e; }
  inline void setEtaVLeft(double e) { _etaVLeft = e; }
  inline void setEtaBRight(double e) { _etaBRight = e; }
  inline void setEtaBLeft(double e) { _etaBLeft = e; }
  inline void setMantleDrag(double m) { _mantleDrag = m ; }
  inline void setViscousShear(double v) { _viscousShear = v ;}
  inline void setBending(double b) { _bending = b ;}

  // re-init
  void zeroAll();

  // computing
  void addDrivingForces(GeoElement*, Direction);
  void computeViscosities(GeoElement*, Direction);

 protected:
  Earth* _earth;

  double _ridgePush;
  double _slabPull;
  double _slabSuction;
  double _etaH;
  double _etaVRight;
  double _etaVLeft;
  double _etaBRight;
  double _etaBLeft;
  double _viscousShear;
  double _mantleDrag;
  double _bending;
  // for computation of velocity of 'virtual' plates:
  double _URight;
  double _ULeft;
};

#endif
