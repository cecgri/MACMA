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
//   author: cecile.grigne@univ-brest.fr

#ifndef TIMER_H_
#define TIMER_H_

#include <QObject>
#include "MACMA/physics/earth.h"


class EarthTimer : public QObject
{

  Q_OBJECT

 public:
  EarthTimer(Earth*);
  ~EarthTimer();

  public slots:
    void isTimeReached(double);
    inline double getNextTimeWrite() { return _nextTimeWrite; }
    inline string getId() { return _id; }

 signals:
    void timeReached();

 private:
    Earth* _earth;
    double _nextTimeWrite;
    string _id;
    static unsigned int _timerCounter;
};


#endif

