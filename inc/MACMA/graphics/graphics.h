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

#ifndef GRAPHICS_H_
#define GRAPHICS_H_

// --- QT
#include <QWidget>
#include <QPainter>
#include <QPicture>

// --- MACMA
#include "MACMA/physics/earth.h"

class Graphics : public QWidget
{
 public:
  Graphics(QWidget* parent, Earth* earth);
  ~Graphics();
  
  QPoint center();
  
  void draw(QPainter* painter);
  void drawEarth(QPainter* painter);
  void drawCell(QPainter* painter, double start, double stop, double orientation, double speedfactor);
  void drawRidge(QPainter* painter, double position, bool selected);
  void drawSubduction(QPainter* painter, string className, double position, double depth, bool selected);
  void drawStaple(QPainter* painter, double position, bool selected);
  void drawPlate(QPainter* painter, double start_angle, double stop_angle, bool selected);
  void drawSection(QPainter* painter, double start_angle, double stop_angle, bool selected);
  void drawAgeOceans(QPainter* painter, double ageMax, double start_angle, double stop_angle);
  void drawContinent(QPainter* painter, bool breakable, double start_angle, double stop_angle, bool selected);
  void drawAgeContinents(QPainter* painter, double ageMax, double start_angle, double stop_angle);
  void drawTime(QPainter* painter, double tMyr, double ageMa);
  void drawTemperature(QPainter* painter, double T);  
  void drawTemperatureModif(QPainter* painter, double T, double ageMyr);  
  void drawSurfaceRatio(QPainter* painter, double S);
  void drawFlux(QPainter* painter, double Q);  
  void drawColorScale(QPainter* painter);
  void initColorScale();

  void saveToPNG(const char* filename);
  
 protected:
  void paintEvent(QPaintEvent* event);
  
  Earth*  _earth;
  QPixmap _imgEarth;
  QRect   _rectEarth;
  
  double _heightFactor;
  
  QPicture _picColorScale;
  map<double, QColor> _ageContColor;
  map<double, QColor> _ageOceanColor;
};

#endif 
