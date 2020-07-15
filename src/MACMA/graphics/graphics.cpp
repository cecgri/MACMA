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

#include "MACMA/graphics/graphics.h"

#define ST_HEIGHT 0.05

// --- Couleurs Cellule
#define CELL_COLD_R 0.0
#define CELL_COLD_G 0.0
#define CELL_COLD_B 255.0

#define CELL_MANTLE_R 255.0
#define CELL_MANTLE_G 189
#define CELL_MANTLE_B 0.0

#define CELL_HOT_R 255.0
#define CELL_HOT_G 0.0
#define CELL_HOT_B 0.0

// --- Couleur WarmingZone
#define WARMING_R 255.0
#define WARMING_G 0.0
#define WARMING_B 0.0

#define CONTINENT_PEN_WIDTH 3.0

// --- 
#define ST_CELL_RMIN 0.54
#define ST_CELL_RMAX 0.965

#define PLATE_PEN_WIDTH 8.0
#define PEN_WIDTH 6.0

using namespace std;

Graphics::Graphics(QWidget* parent, Earth* earth) : QWidget(parent), _earth(earth)
{
  QString filename = ":/images/earth.png";
  _imgEarth.load(filename);
  _heightFactor = _imgEarth.height()*0.344;

  _rectEarth = QRect(-_imgEarth.width()*0.35, -_imgEarth.height()*0.35, 
		     _imgEarth.width()*0.7, _imgEarth.height()*0.7);


  // colors for ocean ages
  _ageOceanColor[25.0] = QColor(250, 0, 0);
  _ageOceanColor[50.0] = QColor(250, 160, 40);
  _ageOceanColor[75.0] = QColor(250, 250, 50);
  _ageOceanColor[100.0] = QColor(140, 250, 50);
  _ageOceanColor[125.0] = QColor(40, 250, 90);
  _ageOceanColor[150.0] = QColor(45, 250, 240);
  _ageOceanColor[175.0] = QColor(20, 120, 250);
  _ageOceanColor[200.0] = QColor(10, 40, 250);
  _ageOceanColor[225.0] = QColor(150, 40, 250);
  _ageOceanColor[250.0] = QColor(250, 40, 250);

  // FORMER CONTINENTS : GREEN
  //_ageContColor[4500.0] = QColor(0, 250, 0);

  // colors for continental ages
  _ageContColor[4500.0] = QColor(0, 0, 80);
  _ageContColor[4250.0] = QColor(0, 0, 120);
  _ageContColor[4000.0] = QColor(0, 0, 160);
  _ageContColor[3750.0] = QColor(0, 0, 200);
  _ageContColor[3500.0] = QColor(0, 0, 220);
  _ageContColor[3250.0] = QColor(0, 0, 250);
  _ageContColor[3000.0] = QColor(40, 40, 250);
  _ageContColor[2750.0] = QColor(90, 90, 250);
  _ageContColor[2500.0] = QColor(120, 120, 250);
  _ageContColor[2250.0] = QColor(160, 160, 250);
  _ageContColor[2000.0] = QColor(250, 250, 250);
  _ageContColor[1750.0] = QColor(250, 200, 200);
  _ageContColor[1500.0] = QColor(250, 160, 160);
  _ageContColor[1250.0] = QColor(250, 120, 120);
  _ageContColor[1000.0] = QColor(250, 80, 80);
  _ageContColor[750.0]  = QColor(250, 0, 0);
  _ageContColor[500.0]  = QColor(175, 0, 0);
  _ageContColor[250.0]  = QColor(100, 0, 0);

  initColorScale();
}

Graphics::~Graphics()
{

}
// --------------------------------------------------
void
Graphics::paintEvent(QPaintEvent* event)
{
  QPainter painter(this);
  painter.setRenderHint(QPainter::Antialiasing);
  draw(&painter);
}
// --------------------------------------------------
QPoint
Graphics::center()
{
  return QPoint(width()/2.0, height()/2.0);
}
// --------------------------------------------------
void
Graphics::draw(QPainter* painter)
{
  EarthState state = _earth->getState();

  // --- Earth
  drawEarth(painter);

  // --- Cells
  for(unsigned int i=0; i<state.cells.size(); i++)
    drawCell(painter, state.cells[i].first.first, state.cells[i].first.second, 
	     state.cells[i].second.first, state.cells[i].second.second);  
  // first.first: startAngle, first.second: stopAngle, 
  // second.first: orientation, second.second: velocity

  // --- Plates
  for(unsigned int i=0; i<state.plates.size(); i++)
    drawPlate(painter, state.plates[i].first, state.plates[i].second, state.selectedPlate == (int)i);

  // --- Sections
  for(unsigned int i=0; i<state.sections.size(); i++)
    drawSection(painter, state.sections[i].first, state.sections[i].second, state.selectedSection == (int)i);

  // --- ages of oceans
  if(state.drawAgeOceans)
    {
      for(unsigned int i=0; i<state.ageOceans.size(); i++)
	{
	  drawAgeOceans(painter, state.ageOceans[i].first,
			state.ageOceans[i].second.first,
			state.ageOceans[i].second.second);
	}
    }

  // --- Continents
   for(unsigned int i=0; i<state.continents.size(); i++)
     {
       drawContinent(painter, state.continents[i].first, 
		     state.continents[i].second.first, state.continents[i].second.second, 
		     state.selectedContinent == (int)i);
       // first: breakable, second.first: startAngle, second.second: stopAngle
     }
   for(unsigned int i=0; i<state.ageContinents.size(); i++)
     {
       drawAgeContinents(painter, state.ageContinents[i].first,
			 state.ageContinents[i].second.first,
			 state.ageContinents[i].second.second);
       // first: ageMax; second.first: rightPosition; second.second: leftPosition
     }

  // --- Interfaces
  for(unsigned int i=0; i<state.interfaces.size(); i++)
    {
      if(state.interfaces[i].first == "Ridge")
	drawRidge(painter, state.interfaces[i].second.first, false);
      if(state.interfaces[i].first == "LeftSubduction" || state.interfaces[i].first == "RightSubduction")
	drawSubduction(painter, state.interfaces[i].first, state.interfaces[i].second.first, state.interfaces[i].second.second, false);
      if(state.interfaces[i].first == "Staple")
	drawStaple(painter, state.interfaces[i].second.first, false);
    }
  if(state.selectedInterface != -1)
    {
      int i = state.selectedInterface;
      if(state.interfaces[i].first == "Ridge")
	drawRidge(painter, state.interfaces[i].second.first, true);
      if(state.interfaces[i].first == "LeftSubduction" || state.interfaces[i].first == "RightSubduction")
	drawSubduction(painter, state.interfaces[i].first, state.interfaces[i].second.first, state.interfaces[i].second.second, true);
      if(state.interfaces[i].first == "Staple")
	drawStaple(painter, state.interfaces[i].second.first, true);
    }

  // --- Time
  drawTime(painter, state.tMyr, state.ageMa);

  // --- Temperature
  drawTemperature(painter, state.T);

  // --- color scale
  if(_earth->runStarted() && _earth->growingContinents())
    {
      drawColorScale(painter);
      drawSurfaceRatio(painter, state.ScontRatio);
    }

  // --- Heat flux
  if(_earth->runStarted())
    drawFlux(painter, state.Q);
}
// --------------------------------------------------
void
Graphics::drawEarth(QPainter* painter)
{
  // --- Background
  QRadialGradient grad(center(), height()/1.5);
  grad.setColorAt(0.0, Qt::white);
  grad.setColorAt(1.0, Qt::black);
  painter->fillRect(0, 0, width(), height(), grad);
    
  // --- Earth
  painter->save();
  painter->translate(width()/2.0, height()/2.0);
  painter->drawPixmap(_rectEarth, _imgEarth); 
  painter->restore();
}
// --------------------------------------------------
void 
Graphics::drawCell(QPainter* painter, double start_angle, double stop_angle, double orientation, double speedfactor)
{
  if(stop_angle <= start_angle)
    stop_angle += 360.0;
  
  start_angle += 0.03 * 180.0 / M_PI;
  stop_angle -= 0.03 * 180.0 / M_PI;

  double deg_start = start_angle;
  double deg_stop = stop_angle;

  double start = - start_angle * M_PI / 180.0;
  double stop = - stop_angle * M_PI / 180.0;

  double trans = 0.2;
  double deg_trans = 180.0 * trans / M_PI;

  bool withMiddle = (deg_trans >= (stop_angle-start_angle)/2.0);
  double deg_end = (deg_stop-deg_start-2.*deg_trans) / 360.;

  QConicalGradient red_orange_blue;
  QConicalGradient blue_orange_red;

  if(!withMiddle)
    {
      red_orange_blue = QConicalGradient(0, 0, deg_start+deg_trans-2.0) ;
      red_orange_blue.setColorAt(0.0, QColor(CELL_HOT_R, CELL_HOT_G, CELL_HOT_B));
      red_orange_blue.setColorAt(0.2*deg_end, QColor(CELL_MANTLE_R, CELL_MANTLE_G, CELL_MANTLE_B));
      red_orange_blue.setColorAt(0.8*deg_end, QColor(CELL_MANTLE_R, CELL_MANTLE_G, CELL_MANTLE_B));
      red_orange_blue.setColorAt(deg_end, QColor(CELL_COLD_R, CELL_COLD_G, CELL_COLD_B));
      
      blue_orange_red = QConicalGradient(0, 0, deg_start+deg_trans-2.0) ;
      blue_orange_red.setColorAt(0.0, QColor(CELL_COLD_R, CELL_COLD_G, CELL_COLD_B));
      blue_orange_red.setColorAt(0.2*deg_end, QColor(CELL_MANTLE_R, CELL_MANTLE_G, CELL_MANTLE_B));
      blue_orange_red.setColorAt(0.8*deg_end, QColor(CELL_MANTLE_R, CELL_MANTLE_G, CELL_MANTLE_B));
      blue_orange_red.setColorAt(deg_end, QColor(CELL_HOT_R, CELL_HOT_G, CELL_HOT_B));
    }
 
  painter->save();
  painter->translate(width()/2.0, height()/2.0);

  /* - pathStart: inside vertical right border of cell 
     - pathEnd: inside vertical left border of cell
     (red if ascending, blue if descending)

     - pathMiddleBottom and pathMiddleTop: inside
     horizontal borders (gradient, if there's room = !withMiddle)

     - pathFullMiddle: dashes inside
     
  */


  // --- Start
  QPainterPath pathStart;

  if(!withMiddle)
    {
      pathStart.moveTo(_heightFactor*ST_CELL_RMIN*cos(start-trans),
		       _heightFactor*ST_CELL_RMIN*sin(start-trans));
      pathStart.cubicTo( _heightFactor*ST_CELL_RMIN*cos(start),
			 _heightFactor*ST_CELL_RMIN*sin(start),
			 _heightFactor*ST_CELL_RMIN*cos(start),
			 _heightFactor*ST_CELL_RMIN*sin(start),
			 _heightFactor*((ST_CELL_RMAX-ST_CELL_RMIN)/2.0+ST_CELL_RMIN)*cos(start),
			 _heightFactor*((ST_CELL_RMAX-ST_CELL_RMIN)/2.0+ST_CELL_RMIN)*sin(start));
      pathStart.cubicTo( _heightFactor*ST_CELL_RMAX*cos(start),
			 _heightFactor*ST_CELL_RMAX*sin(start),
			 _heightFactor*ST_CELL_RMAX*cos(start),
			 _heightFactor*ST_CELL_RMAX*sin(start),
			 _heightFactor*ST_CELL_RMAX*cos(start-trans),
			 _heightFactor*ST_CELL_RMAX*sin(start-trans));
    }
  else
    {
      pathStart.moveTo(_heightFactor*ST_CELL_RMIN*cos(start+(stop-start)/2.0),
		       _heightFactor*ST_CELL_RMIN*sin(start+(stop-start)/2.0));
      pathStart.cubicTo( _heightFactor*ST_CELL_RMIN*cos(start),
			 _heightFactor*ST_CELL_RMIN*sin(start),
			 _heightFactor*ST_CELL_RMIN*cos(start),
			 _heightFactor*ST_CELL_RMIN*sin(start),
			 _heightFactor*((ST_CELL_RMAX-ST_CELL_RMIN)/2.0+ST_CELL_RMIN)*cos(start),
			 _heightFactor*((ST_CELL_RMAX-ST_CELL_RMIN)/2.0+ST_CELL_RMIN)*sin(start));
      pathStart.cubicTo( _heightFactor*ST_CELL_RMAX*cos(start),
			 _heightFactor*ST_CELL_RMAX*sin(start),
			 _heightFactor*ST_CELL_RMAX*cos(start),
			 _heightFactor*ST_CELL_RMAX*sin(start),
			 _heightFactor*ST_CELL_RMAX*cos(start+(stop-start)/2.0),
			 _heightFactor*ST_CELL_RMAX*sin(start+(stop-start)/2.0));
    }

  // --- Middle
  // Bottom
  QPainterPath pathMiddleBottom;
  if(!withMiddle)
    {
      pathMiddleBottom.moveTo(_heightFactor*ST_CELL_RMIN*cos(start-trans),
			      _heightFactor*ST_CELL_RMIN*sin(start-trans));
      pathMiddleBottom.arcTo(0.98*_rectEarth.topLeft().x()*ST_CELL_RMIN,
			     0.98*_rectEarth.topLeft().y()*ST_CELL_RMIN,
			     0.98*_rectEarth.width()*ST_CELL_RMIN,
			     0.98*_rectEarth.height()*ST_CELL_RMIN, 
			     deg_start+deg_trans, deg_stop-deg_start-2.0*deg_trans);
    }

  // Top
  QPainterPath pathMiddleTop;
  if(!withMiddle)
    {
      pathMiddleTop.moveTo(_heightFactor*ST_CELL_RMAX*cos(start-trans),
			   _heightFactor*ST_CELL_RMAX*sin(start-trans));
      pathMiddleTop.arcTo(0.98*_rectEarth.topLeft().x()*ST_CELL_RMAX,
			  0.98*_rectEarth.topLeft().y()*ST_CELL_RMAX,
			  0.98*_rectEarth.width()*ST_CELL_RMAX,
			  0.98*_rectEarth.height()*ST_CELL_RMAX, 
			  deg_start+deg_trans, deg_stop-deg_start-2.0*deg_trans);
    }
    
  // --- End
  QPainterPath pathEnd;

  if(!withMiddle)
    {
      pathEnd.moveTo(_heightFactor*ST_CELL_RMAX*cos(stop+trans),
		     _heightFactor*ST_CELL_RMAX*sin(stop+trans));
      pathEnd.cubicTo( _heightFactor*ST_CELL_RMAX*cos(stop),
		       _heightFactor*ST_CELL_RMAX*sin(stop),
		       _heightFactor*ST_CELL_RMAX*cos(stop),
		       _heightFactor*ST_CELL_RMAX*sin(stop),
		       _heightFactor*((ST_CELL_RMAX-ST_CELL_RMIN)/2.0+ST_CELL_RMIN)*cos(stop),
		       _heightFactor*((ST_CELL_RMAX-ST_CELL_RMIN)/2.0+ST_CELL_RMIN)*sin(stop));
      pathEnd.cubicTo( _heightFactor*ST_CELL_RMIN*cos(stop),
		       _heightFactor*ST_CELL_RMIN*sin(stop),
		       _heightFactor*ST_CELL_RMIN*cos(stop),
		       _heightFactor*ST_CELL_RMIN*sin(stop),
		       _heightFactor*ST_CELL_RMIN*cos(stop+trans),
		       _heightFactor*ST_CELL_RMIN*sin(stop+trans));
    }
  else
    {
      pathEnd.moveTo(_heightFactor*ST_CELL_RMAX*cos(start+(stop-start)/2.0),
		     _heightFactor*ST_CELL_RMAX*sin(start+(stop-start)/2.0));
      pathEnd.cubicTo( _heightFactor*ST_CELL_RMAX*cos(stop),
		       _heightFactor*ST_CELL_RMAX*sin(stop),
		       _heightFactor*ST_CELL_RMAX*cos(stop),
		       _heightFactor*ST_CELL_RMAX*sin(stop),
		       _heightFactor*((ST_CELL_RMAX-ST_CELL_RMIN)/2.0+ST_CELL_RMIN)*cos(stop),
		       _heightFactor*((ST_CELL_RMAX-ST_CELL_RMIN)/2.0+ST_CELL_RMIN)*sin(stop));
      pathEnd.cubicTo( _heightFactor*ST_CELL_RMIN*cos(stop),
		       _heightFactor*ST_CELL_RMIN*sin(stop),
		       _heightFactor*ST_CELL_RMIN*cos(stop),
		       _heightFactor*ST_CELL_RMIN*sin(stop),
		       _heightFactor*ST_CELL_RMIN*cos(start+(stop-start)/2.0),
		       _heightFactor*ST_CELL_RMIN*sin(start+(stop-start)/2.0));
    }

  // --- Full Middle
  QPainterPath pathFullMiddle;

  // Start
  if(!withMiddle)
    {
      pathFullMiddle.moveTo(_heightFactor*ST_CELL_RMIN*cos(start-trans),
			    _heightFactor*ST_CELL_RMIN*sin(start-trans));
      pathFullMiddle.cubicTo( _heightFactor*ST_CELL_RMIN*cos(start),
			      _heightFactor*ST_CELL_RMIN*sin(start),
			      _heightFactor*ST_CELL_RMIN*cos(start),
			      _heightFactor*ST_CELL_RMIN*sin(start),
			      _heightFactor*((ST_CELL_RMAX-ST_CELL_RMIN)/2.0+ST_CELL_RMIN)*cos(start),
			      _heightFactor*((ST_CELL_RMAX-ST_CELL_RMIN)/2.0+ST_CELL_RMIN)*sin(start));
      pathFullMiddle.cubicTo( _heightFactor*ST_CELL_RMAX*cos(start),
			      _heightFactor*ST_CELL_RMAX*sin(start),
			      _heightFactor*ST_CELL_RMAX*cos(start),
			      _heightFactor*ST_CELL_RMAX*sin(start),
			      _heightFactor*ST_CELL_RMAX*cos(start-trans),
			      _heightFactor*ST_CELL_RMAX*sin(start-trans));
    }
  else
    {
      pathFullMiddle.moveTo(_heightFactor*ST_CELL_RMIN*cos(start+(stop-start)/2.0),
			    _heightFactor*ST_CELL_RMIN*sin(start+(stop-start)/2.0));
      pathFullMiddle.cubicTo( _heightFactor*ST_CELL_RMIN*cos(start),
			      _heightFactor*ST_CELL_RMIN*sin(start),
			      _heightFactor*ST_CELL_RMIN*cos(start),
			      _heightFactor*ST_CELL_RMIN*sin(start),
			      _heightFactor*((ST_CELL_RMAX-ST_CELL_RMIN)/2.0+ST_CELL_RMIN)*cos(start),
			      _heightFactor*((ST_CELL_RMAX-ST_CELL_RMIN)/2.0+ST_CELL_RMIN)*sin(start));
      pathFullMiddle.cubicTo( _heightFactor*ST_CELL_RMAX*cos(start),
			      _heightFactor*ST_CELL_RMAX*sin(start),
			      _heightFactor*ST_CELL_RMAX*cos(start),
			      _heightFactor*ST_CELL_RMAX*sin(start),
			      _heightFactor*ST_CELL_RMAX*cos(start+(stop-start)/2.0),
			      _heightFactor*ST_CELL_RMAX*sin(start+(stop-start)/2.0));
    }

  // Top
  if(!withMiddle)
    {
      pathFullMiddle.arcTo(0.98*_rectEarth.topLeft().x()*ST_CELL_RMAX,
			   0.98*_rectEarth.topLeft().y()*ST_CELL_RMAX,
			   0.98*_rectEarth.width()*ST_CELL_RMAX,
			   0.98*_rectEarth.height()*ST_CELL_RMAX, 
			   deg_start+deg_trans, deg_stop-deg_start-2.0*deg_trans);
    }

  // End
  if(!withMiddle)
    {
      pathFullMiddle.cubicTo( _heightFactor*ST_CELL_RMAX*cos(stop),
			      _heightFactor*ST_CELL_RMAX*sin(stop),
			      _heightFactor*ST_CELL_RMAX*cos(stop),
			      _heightFactor*ST_CELL_RMAX*sin(stop),
			      _heightFactor*((ST_CELL_RMAX-ST_CELL_RMIN)/2.0+ST_CELL_RMIN)*cos(stop),
			      _heightFactor*((ST_CELL_RMAX-ST_CELL_RMIN)/2.0+ST_CELL_RMIN)*sin(stop));
      pathFullMiddle.cubicTo( _heightFactor*ST_CELL_RMIN*cos(stop),
			      _heightFactor*ST_CELL_RMIN*sin(stop),
			      _heightFactor*ST_CELL_RMIN*cos(stop),
			      _heightFactor*ST_CELL_RMIN*sin(stop),
			      _heightFactor*ST_CELL_RMIN*cos(stop+trans),
			      _heightFactor*ST_CELL_RMIN*sin(stop+trans));
    }
  else
    {
      pathFullMiddle.moveTo(_heightFactor*ST_CELL_RMAX*cos(start+(stop-start)/2.0),
			    _heightFactor*ST_CELL_RMAX*sin(start+(stop-start)/2.0));
      pathFullMiddle.cubicTo( _heightFactor*ST_CELL_RMAX*cos(stop),
			      _heightFactor*ST_CELL_RMAX*sin(stop),
			      _heightFactor*ST_CELL_RMAX*cos(stop),
			      _heightFactor*ST_CELL_RMAX*sin(stop),
			      _heightFactor*((ST_CELL_RMAX-ST_CELL_RMIN)/2.0+ST_CELL_RMIN)*cos(stop),
			      _heightFactor*((ST_CELL_RMAX-ST_CELL_RMIN)/2.0+ST_CELL_RMIN)*sin(stop));
      pathFullMiddle.cubicTo( _heightFactor*ST_CELL_RMIN*cos(stop),
			      _heightFactor*ST_CELL_RMIN*sin(stop),
			      _heightFactor*ST_CELL_RMIN*cos(stop),
			      _heightFactor*ST_CELL_RMIN*sin(stop),
			      _heightFactor*ST_CELL_RMIN*cos(start+(stop-start)/2.0),
			      _heightFactor*ST_CELL_RMIN*sin(start+(stop-start)/2.0));
    }

  // Bottom
  if(!withMiddle)
    {
      pathFullMiddle.arcTo(0.98*_rectEarth.topLeft().x()*ST_CELL_RMIN,
			   0.98*_rectEarth.topLeft().y()*ST_CELL_RMIN,
			   0.98*_rectEarth.width()*ST_CELL_RMIN,
			   0.98*_rectEarth.height()*ST_CELL_RMIN, 
			   deg_stop-deg_trans, -(deg_stop-deg_start-2.0*deg_trans));
    }

  // --- Full Ext
  QPainterPath pathFullExt;
  double offset = -0.02; 

  // Start
  if(!withMiddle)
    {
      pathFullExt.moveTo(0.98*_heightFactor*ST_CELL_RMIN*cos(start-trans),
                         0.98*_heightFactor*ST_CELL_RMIN*sin(start-trans));
      pathFullExt.cubicTo( 0.98*_heightFactor*ST_CELL_RMIN*cos(start-offset),
                           0.98*_heightFactor*ST_CELL_RMIN*sin(start-offset),
                           0.98*_heightFactor*ST_CELL_RMIN*cos(start-offset),
                           0.98*_heightFactor*ST_CELL_RMIN*sin(start-offset),
                           _heightFactor*((ST_CELL_RMAX-ST_CELL_RMIN)/2.0+ST_CELL_RMIN)*cos(start-offset),
                           _heightFactor*((ST_CELL_RMAX-ST_CELL_RMIN)/2.0+ST_CELL_RMIN)*sin(start-offset));
      pathFullExt.cubicTo( 1.02*_heightFactor*ST_CELL_RMAX*cos(start-offset),
                           1.02*_heightFactor*ST_CELL_RMAX*sin(start-offset),
                           _heightFactor*ST_CELL_RMAX*cos(start-offset),
                           _heightFactor*ST_CELL_RMAX*sin(start-offset),
                           1.02*_heightFactor*ST_CELL_RMAX*cos(start-trans),
                           1.02*_heightFactor*ST_CELL_RMAX*sin(start-trans));
    }
  else
    {
      pathFullExt.moveTo(0.98*_heightFactor*ST_CELL_RMIN*cos(start+(stop-start)/2.0),
			 0.98*_heightFactor*ST_CELL_RMIN*sin(start+(stop-start)/2.0));
      pathFullExt.cubicTo( 0.98*_heightFactor*ST_CELL_RMIN*cos(start-offset),
                           0.98*_heightFactor*ST_CELL_RMIN*sin(start-offset),
                           0.98*_heightFactor*ST_CELL_RMIN*cos(start-offset),
                           0.98*_heightFactor*ST_CELL_RMIN*sin(start-offset),
                           _heightFactor*((ST_CELL_RMAX-ST_CELL_RMIN)/2.0+ST_CELL_RMIN)*cos(start-offset),
                           _heightFactor*((ST_CELL_RMAX-ST_CELL_RMIN)/2.0+ST_CELL_RMIN)*sin(start-offset));
      pathFullExt.cubicTo( 1.02*_heightFactor*ST_CELL_RMAX*cos(start-offset),
                           1.02*_heightFactor*ST_CELL_RMAX*sin(start-offset),
                           _heightFactor*ST_CELL_RMAX*cos(start-offset),
                           _heightFactor*ST_CELL_RMAX*sin(start-offset),
                           1.02*_heightFactor*ST_CELL_RMAX*cos(start+(stop-start)/2.0),
                           1.02*_heightFactor*ST_CELL_RMAX*sin(start+(stop-start)/2.0));
    }

  // Top
  if(!withMiddle)
    {
      pathFullExt.arcTo(_rectEarth.topLeft().x()*ST_CELL_RMAX,
			_rectEarth.topLeft().y()*ST_CELL_RMAX,
			_rectEarth.width()*ST_CELL_RMAX,
			_rectEarth.height()*ST_CELL_RMAX, 
			deg_start+deg_trans, deg_stop-deg_start-2.0*deg_trans);
    }

  // End
  if(!withMiddle)
    {
      pathFullExt.cubicTo( 1.01*_heightFactor*ST_CELL_RMAX*cos(stop+offset),
                           1.01*_heightFactor*ST_CELL_RMAX*sin(stop+offset),
                           _heightFactor*ST_CELL_RMAX*cos(stop+offset),
                           _heightFactor*ST_CELL_RMAX*sin(stop+offset),
                           _heightFactor*((ST_CELL_RMAX-ST_CELL_RMIN)/2.0+ST_CELL_RMIN)*cos(stop+offset),
                           _heightFactor*((ST_CELL_RMAX-ST_CELL_RMIN)/2.0+ST_CELL_RMIN)*sin(stop+offset));
      pathFullExt.cubicTo( 0.98*_heightFactor*ST_CELL_RMIN*cos(stop+offset),
                           0.98*_heightFactor*ST_CELL_RMIN*sin(stop+offset),
                           0.98*_heightFactor*ST_CELL_RMIN*cos(stop+offset),
                           0.98*_heightFactor*ST_CELL_RMIN*sin(stop+offset),
                           0.98*_heightFactor*ST_CELL_RMIN*cos(stop+trans),
                           0.98*_heightFactor*ST_CELL_RMIN*sin(stop+trans));
    }
  else
    {
      pathFullExt.cubicTo( 1.01*_heightFactor*ST_CELL_RMAX*cos(stop+offset),
                           1.01*_heightFactor*ST_CELL_RMAX*sin(stop+offset),
                           _heightFactor*ST_CELL_RMAX*cos(stop+offset),
                           _heightFactor*ST_CELL_RMAX*sin(stop+offset),
                           _heightFactor*((ST_CELL_RMAX-ST_CELL_RMIN)/2.0+ST_CELL_RMIN)*cos(stop+offset),
                           _heightFactor*((ST_CELL_RMAX-ST_CELL_RMIN)/2.0+ST_CELL_RMIN)*sin(stop+offset));
      pathFullExt.cubicTo( 0.98*_heightFactor*ST_CELL_RMIN*cos(stop+offset),
                           0.98*_heightFactor*ST_CELL_RMIN*sin(stop+offset),
                           0.98*_heightFactor*ST_CELL_RMIN*cos(stop+offset),
                           0.98*_heightFactor*ST_CELL_RMIN*sin(stop+offset),
                           0.98*_heightFactor*ST_CELL_RMIN*cos(start+(stop-start)/2.0),
                           0.98*_heightFactor*ST_CELL_RMIN*sin(start+(stop-start)/2.0));
    }

  // Bottom
  if(!withMiddle)
    {
      pathFullExt.arcTo(0.96*_rectEarth.topLeft().x()*ST_CELL_RMIN,
			0.96*_rectEarth.topLeft().y()*ST_CELL_RMIN,
			0.96*_rectEarth.width()*ST_CELL_RMIN,
			0.96*_rectEarth.height()*ST_CELL_RMIN, 
			deg_stop-deg_trans, -(deg_stop-deg_start-2.0*deg_trans));
    }

  // --- Full Int
  QPainterPath pathFullInt;

  // Start
  offset = 0.02; 
  if(!withMiddle)
    {
      pathFullInt.moveTo(1.02*_heightFactor*ST_CELL_RMIN*cos(start-trans),
                         1.02*_heightFactor*ST_CELL_RMIN*sin(start-trans));
      pathFullInt.cubicTo( 1.02*_heightFactor*ST_CELL_RMIN*cos(start-offset),
                           1.02*_heightFactor*ST_CELL_RMIN*sin(start-offset),
                           1.02*_heightFactor*ST_CELL_RMIN*cos(start-offset),
                           1.02*_heightFactor*ST_CELL_RMIN*sin(start-offset),
                           _heightFactor*((ST_CELL_RMAX-ST_CELL_RMIN)/2.0+ST_CELL_RMIN)*cos(start-offset),
                           _heightFactor*((ST_CELL_RMAX-ST_CELL_RMIN)/2.0+ST_CELL_RMIN)*sin(start-offset));
      pathFullInt.cubicTo( 0.98*_heightFactor*ST_CELL_RMAX*cos(start-offset),
                           0.98*_heightFactor*ST_CELL_RMAX*sin(start-offset),
                           0.985*_heightFactor*ST_CELL_RMAX*cos(start-offset),
                           0.985*_heightFactor*ST_CELL_RMAX*sin(start-offset),
                           0.985*_heightFactor*ST_CELL_RMAX*cos(start-trans),
                           0.985*_heightFactor*ST_CELL_RMAX*sin(start-trans));
    }
  else
    {
      pathFullInt.moveTo(1.02*_heightFactor*ST_CELL_RMIN*cos(start+(stop-start)/2.0),
                         1.02*_heightFactor*ST_CELL_RMIN*sin(start+(stop-start)/2.0));
      pathFullInt.cubicTo( 1.02*_heightFactor*ST_CELL_RMIN*cos(start-offset),
                           1.02*_heightFactor*ST_CELL_RMIN*sin(start-offset),
                           1.02*_heightFactor*ST_CELL_RMIN*cos(start-offset),
                           1.02*_heightFactor*ST_CELL_RMIN*sin(start-offset),
                           _heightFactor*((ST_CELL_RMAX-ST_CELL_RMIN)/2.0+ST_CELL_RMIN)*cos(start-offset),
                           _heightFactor*((ST_CELL_RMAX-ST_CELL_RMIN)/2.0+ST_CELL_RMIN)*sin(start-offset));
      pathFullInt.cubicTo( 0.98*_heightFactor*ST_CELL_RMAX*cos(start-offset),
                           0.98*_heightFactor*ST_CELL_RMAX*sin(start-offset),
                           0.985*_heightFactor*ST_CELL_RMAX*cos(start-offset),
                           0.985*_heightFactor*ST_CELL_RMAX*sin(start-offset),
                           0.985*_heightFactor*ST_CELL_RMAX*cos(start+(stop-start)/2.0),
                           0.985*_heightFactor*ST_CELL_RMAX*sin(start+(stop-start)/2.0));
    }
        
  // Top
  if(!withMiddle)
    {
      pathFullInt.arcTo(0.965*_rectEarth.topLeft().x()*ST_CELL_RMAX,
			0.965*_rectEarth.topLeft().y()*ST_CELL_RMAX,
			0.965*_rectEarth.width()*ST_CELL_RMAX,
			0.965*_rectEarth.height()*ST_CELL_RMAX, 
			deg_start+deg_trans, deg_stop-deg_start-2.0*deg_trans);
    }

  // End
  if(!withMiddle)
    {
      pathFullInt.cubicTo(0.98*_heightFactor*ST_CELL_RMAX*cos(stop+offset),
			  0.98*_heightFactor*ST_CELL_RMAX*sin(stop+offset),
			  _heightFactor*ST_CELL_RMAX*cos(stop+offset),
			  _heightFactor*ST_CELL_RMAX*sin(stop+offset),
			  _heightFactor*((ST_CELL_RMAX-ST_CELL_RMIN)/2.0+ST_CELL_RMIN)*cos(stop+offset),
			  _heightFactor*((ST_CELL_RMAX-ST_CELL_RMIN)/2.0+ST_CELL_RMIN)*sin(stop+offset));
      pathFullInt.cubicTo( 1.02*_heightFactor*ST_CELL_RMIN*cos(stop+offset),
                           1.02*_heightFactor*ST_CELL_RMIN*sin(stop+offset),
                           1.02*_heightFactor*ST_CELL_RMIN*cos(stop+offset),
                           1.02*_heightFactor*ST_CELL_RMIN*sin(stop+offset),
                           1.02*_heightFactor*ST_CELL_RMIN*cos(stop+trans),
                           1.02*_heightFactor*ST_CELL_RMIN*sin(stop+trans));
    }
  else
    {
      pathFullInt.cubicTo(0.98*_heightFactor*ST_CELL_RMAX*cos(stop+offset),
			  0.98*_heightFactor*ST_CELL_RMAX*sin(stop+offset),
			  _heightFactor*ST_CELL_RMAX*cos(stop+offset),
			  _heightFactor*ST_CELL_RMAX*sin(stop+offset),
			  _heightFactor*((ST_CELL_RMAX-ST_CELL_RMIN)/2.0+ST_CELL_RMIN)*cos(stop+offset),
			  _heightFactor*((ST_CELL_RMAX-ST_CELL_RMIN)/2.0+ST_CELL_RMIN)*sin(stop+offset));
      pathFullInt.cubicTo( 1.02*_heightFactor*ST_CELL_RMIN*cos(stop+offset),
                           1.02*_heightFactor*ST_CELL_RMIN*sin(stop+offset),
                           1.02*_heightFactor*ST_CELL_RMIN*cos(stop+offset),
                           1.02*_heightFactor*ST_CELL_RMIN*sin(stop+offset),
                           1.02*_heightFactor*ST_CELL_RMIN*cos(start+(stop-start)/2.0),
                           1.02*_heightFactor*ST_CELL_RMIN*sin(start+(stop-start)/2.0));
    }

  // Bottom
  if(!withMiddle)
    {
      pathFullInt.arcTo(1.01*_rectEarth.topLeft().x()*ST_CELL_RMIN,
			1.01*_rectEarth.topLeft().y()*ST_CELL_RMIN,
			1.01*_rectEarth.width()*ST_CELL_RMIN,
			1.01*_rectEarth.height()*ST_CELL_RMIN, 
			deg_stop-deg_trans, -(deg_stop-deg_start-2.0*deg_trans));
    }
    
  // Draw
  {
    // Start
    QPen pen;
    pen.setWidth(1.5*PEN_WIDTH);
    if(orientation > 0)
      pen.setColor(QColor(CELL_HOT_R, CELL_HOT_G, CELL_HOT_B));
    else
      pen.setColor(QColor(CELL_COLD_R, CELL_COLD_G, CELL_COLD_B));
    painter->setPen(pen);
    painter->drawPath(pathStart);

    // End
    if(orientation > 0)
      pen.setColor(QColor(CELL_COLD_R, CELL_COLD_G, CELL_COLD_B));
    else
      pen.setColor(QColor(CELL_HOT_R, CELL_HOT_G, CELL_HOT_B));
    painter->setPen(pen);
    painter->drawPath(pathEnd);
  }

  if(!withMiddle)
    {
      // Bottom + Top
      QPen mPen;
      if(orientation > 0)
	mPen.setBrush(QBrush(red_orange_blue));
      else
	mPen.setBrush(QBrush(blue_orange_red));
      mPen.setWidth(1.5*PEN_WIDTH);
      painter->setPen(mPen);
      painter->drawPath(pathMiddleBottom);
      painter->drawPath(pathMiddleTop);
    }

  {
    // Full Middle
    QPen pen; 
    pen.setStyle(Qt::DotLine);
    pen.setWidth(1.5*PEN_WIDTH);
    pen.setBrush(Qt::black);
    QVector<qreal> dashes;
    dashes << 0.02 << 10.0;
    pen.setDashPattern(dashes);
    double dashOffset = -orientation * fmod(speedfactor * ((double) (clock())/((double)CLOCKS_PER_SEC)), 10.0);
    pen.setDashOffset(dashOffset);
    pen.setJoinStyle(Qt::RoundJoin);
    painter->setPen(pen);
    painter->setBrush(Qt::NoBrush);
    painter->drawPath(pathFullMiddle);
  }

  {
    // Full Ext
    QPen pen; 
    pen.setWidth(0.5*PEN_WIDTH);
    pen.setBrush(Qt::black);
    pen.setJoinStyle(Qt::RoundJoin);
    painter->setPen(pen);
    painter->setBrush(Qt::NoBrush);
    painter->drawPath(pathFullExt);
  }

  {
    // Full Int
    QPen pen; 
    pen.setWidth(0.5*PEN_WIDTH);
    pen.setBrush(Qt::black);
    pen.setJoinStyle(Qt::RoundJoin);
    painter->setPen(pen);
    painter->setBrush(QBrush(QColor(CELL_MANTLE_R, CELL_MANTLE_G, CELL_MANTLE_B)));
    pathFullInt.setFillRule(Qt::WindingFill);

    painter->drawPath(pathFullInt);
  }

  {
    // Full Int
    QPen pen; 
    pen.setWidth(0.5*PEN_WIDTH);
    pen.setBrush(Qt::black);
    pen.setJoinStyle(Qt::RoundJoin);
    painter->setPen(pen);
    painter->setBrush(QPixmap(":/images/magma.png"));
    pathFullInt.setFillRule(Qt::WindingFill);

    painter->drawPath(pathFullInt);
  }

  painter->restore();
}
// --------------------------------------------------
void 
Graphics::drawRidge(QPainter* painter, double position, bool selected)
{
  position = -position * M_PI / 180.0;

  if(selected)
    {
      QPen pen(Qt::green);
      pen.setWidth(1.);
      painter->setPen(pen);
      painter->setBrush(QBrush(Qt::green));
    }
  else
    {
      QPen pen(Qt::black);
      pen.setWidth(1.);
      painter->setPen(pen);
      painter->setBrush(QBrush(Qt::black));
    }

  {
    painter->save();
    painter->translate(width()/2.0, height()/2.0);
    QPolygon polygon;
    polygon << QPoint(_heightFactor*0.98*cos(position-0.008), _heightFactor*0.98*sin(position-0.008))
	    << QPoint(cos(position-0.04)*_heightFactor*1.02, sin(position-0.04)*_heightFactor*1.02) 
	    << QPoint(cos(position-0.008)*1.06*_heightFactor, sin(position-0.008)*1.06*_heightFactor);
    painter->drawPolygon(polygon);   
    painter->restore();
  }

  {
    painter->save();
    painter->translate(width()/2.0, height()/2.0);
    QPolygon polygon;
    polygon << QPoint(_heightFactor*0.98*cos(position+0.008), _heightFactor*0.98*sin(position+0.008))
	    << QPoint(cos(position+0.04)*_heightFactor*1.02, sin(position+0.04)*_heightFactor*1.02) 
	    << QPoint(cos(position+0.008)*1.06*_heightFactor, sin(position+0.008)*1.06*_heightFactor);
    painter->drawPolygon(polygon);   
    painter->restore();
  }
}
// --------------------------------------------------
void 
Graphics::drawSubduction(QPainter* painter, string className, double position, double depth, bool selected)
{
  position = -position * M_PI / 180.0;
  
  if(selected)
    {
      QPen pen(Qt::green);
      pen.setWidth(PLATE_PEN_WIDTH);
      painter->setPen(pen);
      painter->setBrush(Qt::NoBrush);
    }
  else
    {
      QPen pen(Qt::black);
      pen.setWidth(PLATE_PEN_WIDTH);
      painter->setPen(pen);
      painter->setBrush(Qt::NoBrush);
    }
  
  double depthPerCent = min(depth / 800e3, 1.0);
  depthPerCent *= 0.1;
  
  painter->save();
  painter->translate(width()/2.0, height()/2.0);
  if(className == "RightSubduction")
    {
      QPainterPath path;
      path.moveTo(1.01*_heightFactor*cos(position), 1.01*_heightFactor*sin(position));
      path.cubicTo(1.01*_heightFactor*cos(position), 1.01*_heightFactor*sin(position),
		   _heightFactor*0.99*cos(position+0.07), _heightFactor*0.99*sin(position+0.07),
		   _heightFactor*(0.96-depthPerCent)*cos(position+0.10), 
		   _heightFactor*(0.96-depthPerCent)*sin(position+0.10));
      painter->drawPath(path);
    }
  else // LEFT
    {
      QPainterPath path;
      path.moveTo(1.01*_heightFactor*cos(position), 1.01*_heightFactor*sin(position));
      path.cubicTo(1.01*_heightFactor*cos(position), 1.01*_heightFactor*sin(position),
		   _heightFactor*0.99*cos(position-0.07), _heightFactor*0.99*sin(position-0.07),
		   _heightFactor*(0.96-depthPerCent)*cos(position-0.10), 
		   _heightFactor*(0.96-depthPerCent)*sin(position-0.10));
      painter->drawPath(path);
    }
  painter->restore();
}
// --------------------------------------------------
void 
Graphics::drawStaple(QPainter* painter, double position, bool selected)
{
  position = -position * M_PI / 180.0;

  if(selected)
    {
      QPen pen(Qt::green);
      pen.setWidth(PEN_WIDTH);
      painter->setPen(pen);
      painter->setBrush(QBrush(Qt::green));
    }
  else
    {
      QPen pen(Qt::black);
      pen.setWidth(PEN_WIDTH);
      painter->setPen(pen);
      painter->setBrush(QBrush(Qt::black));
    }

  painter->save();
  painter->translate(width()/2.0, height()/2.0);
  QPointF c(1.02*_heightFactor*cos(position), 1.02*_heightFactor*sin(position));
  painter->drawEllipse(c, PEN_WIDTH*0.8, PEN_WIDTH*0.8);
  painter->restore();
}
// --------------------------------------------------
void 
Graphics::drawPlate(QPainter* painter, double start_angle, double stop_angle, bool selected)
{
  if(selected)
    {
      QPen pen(Qt::green);
      pen.setWidth(PLATE_PEN_WIDTH);
      painter->setPen(pen);
    }
  else
    {
      QPen pen(Qt::black);
      pen.setWidth(PLATE_PEN_WIDTH);
      painter->setPen(pen);
    }
  
  if(stop_angle <= start_angle)
    stop_angle += 360.0;
  
  start_angle = start_angle * 16.0; // 
  stop_angle = stop_angle * 16.0 - start_angle;
  
  painter->save();
  painter->translate(width()/2.0, height()/2.0);
  painter->drawArc(_rectEarth, start_angle, stop_angle);
  painter->restore();
}
// --------------------------------------------------
void 
Graphics::drawSection(QPainter* painter, double start_angle, double stop_angle, bool selected)
{
  if(selected)
    {
      QPen pen(Qt::red);
      pen.setWidth(PLATE_PEN_WIDTH);
      painter->setPen(pen);
      
      if(stop_angle <= start_angle)
	stop_angle += 360.0;
      
      start_angle = start_angle * 16.0; // drawArc uses 1/16th of degrees
      stop_angle = stop_angle * 16.0 - start_angle;
      
      painter->save();
      painter->translate(width()/2.0, height()/2.0);
      painter->drawArc(_rectEarth, start_angle, stop_angle);
      painter->restore();
    }
}
// --------------------------------------------------
void 
Graphics::drawContinent(QPainter* painter, bool breakable, double start_angle, double stop_angle, bool selected)
{
  if(stop_angle <= start_angle)
    stop_angle += 360.0; // in deg
  
  double start = start_angle * M_PI / 180.0;  // in rad
  double stop = stop_angle * M_PI / 180.0; 
  
  if(selected)
    {
      QPen pen(Qt::red);
      pen.setWidth(PEN_WIDTH);
      painter->setPen(pen);
    }
  else
    {
      QPen pen(Qt::black);
      pen.setWidth(CONTINENT_PEN_WIDTH);
      painter->setPen(pen);
    }
  //  painter->setBrush(QBrush(QColor(CONTINENT_R, CONTINENT_G, CONTINENT_B)));
  //painter->setBrush(QBrush(QColor(AGE4500_R, AGE4500_G, AGE4500_B)));  
  painter->setBrush(QBrush(_ageContColor.find(4500.0)->second));

  painter->save();
  painter->translate(width()/2.0, height()/2.0);
  
  QPainterPath path;
  path.moveTo(1.01*_heightFactor*cos(start), -1.01*_heightFactor*sin(start));
  path.arcTo(_rectEarth.topLeft().x()*0.99,
	     _rectEarth.topLeft().y()*0.99,
	     _rectEarth.width()*0.99,
	     _rectEarth.height()*0.99, 
	     start_angle, stop_angle-start_angle);
  path.lineTo(1.05*_heightFactor*cos(stop), -1.05*_heightFactor*sin(stop));
  path.arcTo(_rectEarth.topLeft().x()*1.05,
	     _rectEarth.topLeft().y()*1.05,
	     _rectEarth.width()*1.05,
	     _rectEarth.height()*1.05, 
	     stop_angle, start_angle-stop_angle);
  path.lineTo(1.01*_heightFactor*cos(start), -1.01*_heightFactor*sin(start));
  painter->drawPath(path);
  painter->restore();
  
  // Warming Zone
  if(breakable)
    {
      QPen p(QColor(255, 0, 0, 50));
      p.setCapStyle(Qt::FlatCap);
      p.setWidth(2.0);
      painter->setPen(p);
      painter->setBrush(QColor(255, 0, 0, 100));

      painter->save();
      painter->translate(width()/2.0, height()/2.0);
      QPainterPath path2;
      path2.moveTo(0.9182*_heightFactor*cos(start), -0.9182*_heightFactor*sin(start));
      path2.arcTo(_rectEarth.topLeft().x()*0.9,
      		  _rectEarth.topLeft().y()*0.9,
      		  _rectEarth.width()*.9,
      		  _rectEarth.height()*0.9, 
      		  start_angle, stop_angle-start_angle);
      path2.lineTo(0.97*_heightFactor*cos(stop), -0.97*_heightFactor*sin(stop));
      
      path2.arcTo(_rectEarth.topLeft().x()*0.97,
      		  _rectEarth.topLeft().y()*0.97,
      		  _rectEarth.width()*0.97,
      		  _rectEarth.height()*0.97, 
      		  stop_angle, start_angle-stop_angle);
      path2.closeSubpath();
      painter->drawPath(path2);
      painter->restore();
    }
}
// --------------------------------------------------
void
Graphics::drawAgeOceans(QPainter* painter, double ageMax, double start_angle, double stop_angle)
{
  if(stop_angle == start_angle)
    return;
  if(stop_angle < start_angle)
    stop_angle += 360.0; // in deg

  double start = start_angle * 16.0;  // in 1/16th of degrees
  double span = stop_angle * 16.0 - start;
 
  double age = 25.0;
  double maxAge = 250.0;
  double step = 25.0;

  QPen p;
  p.setWidth(PLATE_PEN_WIDTH);

  bool found = false;
  while(true)
    {
      if(age >= ageMax) 
	{
	  p.setColor(_ageOceanColor.find(age)->second);
	  found = true;
	  break;
	}

      age += step;
      if(age > maxAge)
	break;
    }

  if(!found)
    p.setColor(_ageOceanColor.find(maxAge)->second);
  
  painter->save();

  painter->setPen(p);

  painter->translate(width()/2.0, height()/2.0);
  painter->drawArc(_rectEarth, start, span);

  painter->restore();
}
// --------------------------------------------------
void
Graphics::drawAgeContinents(QPainter* painter, double ageMax, double start_angle, double stop_angle)
{
  if(stop_angle == start_angle)
    return;
  if(stop_angle < start_angle)
    stop_angle += 360.0; // in deg
  
  double start = start_angle * M_PI / 180.0;  // in rad
  double stop = stop_angle * M_PI / 180.0; 
  
  QPen p(Qt::NoPen);
  painter->setPen(p);
    
  double age = 4500.0;
  double minAge = 500.0;
  double step = 250.0;
  
  bool found = false;
  while(true)
    {
      if(ageMax > age)
	{
	  painter->setBrush(_ageContColor.find(age)->second);
	  found = true;
	  break;
	}

      age -= step;
      if(age < minAge)
	break;
    }
  if(!found)
    painter->setBrush(_ageContColor.find(250.0)->second);
    
  painter->save();
  painter->translate(width()/2.0, height()/2.0);
  
  QPainterPath path;
  path.moveTo(1.01*_heightFactor*cos(start), -1.01*_heightFactor*sin(start));
  path.arcTo(_rectEarth.topLeft().x()*0.99,
	     _rectEarth.topLeft().y()*0.99,
	     _rectEarth.width()*0.99,
	     _rectEarth.height()*0.99, 
	     start_angle, stop_angle-start_angle);
  path.lineTo(1.05*_heightFactor*cos(stop), -1.05*_heightFactor*sin(stop));
  path.arcTo(_rectEarth.topLeft().x()*1.05,
	     _rectEarth.topLeft().y()*1.05,
	     _rectEarth.width()*1.05,
	     _rectEarth.height()*1.05, 
	     stop_angle, start_angle-stop_angle);
  path.lineTo(1.01*_heightFactor*cos(start), -1.01*_heightFactor*sin(start));
  painter->drawPath(path);
  painter->restore();
}
// --------------------------------------------------
void
Graphics::initColorScale()
{
  QPainter paintColorScale;
  paintColorScale.begin(&_picColorScale);
  QPen pen(Qt::NoPen);
  paintColorScale.setPen(pen);
  
  QRectF rectColor;
  QRectF rectText;

  double x = 0.0;
  double y = 0.0;
  double widthRect = _rectEarth.width() * 0.020;
  double widthText = _rectEarth.width() * 0.050;
  double tickYPos = y + widthRect;
  double tickLength = _rectEarth.width() * 0.0025;
  double textPos = tickYPos + 1.5 * tickLength;
  double age = 4500.0;
  while(true)
    {
      // colored square
      if(age > 0.0)
	{
	  paintColorScale.setBrush(QBrush(_ageContColor.find(age)->second));
	  rectColor = QRect(x, y, 
			    widthRect, widthRect);
	  paintColorScale.drawRect(rectColor);
	}

      // tick and text for age
      rectText = QRect(x - 0.5*widthText, tickYPos-2.0*tickLength,
		       widthText, widthText);

      paintColorScale.save();
      if( age == 4000.0 || age == 3000.0 
	  || age == 2000.0 || age == 1000.0
	  || age == 0.0 )
	{
	  QPen p(Qt::white);
	  // tick
	  p.setWidth(2);
	  paintColorScale.setPen(p);
	  QLineF line(x, tickYPos, x, tickYPos + tickLength);
	  paintColorScale.drawLine(line);
	  // text
	  paintColorScale.setFont(QFont("Arial", 12));
	  QString textAge(::doubleToString(age/1000.0, 1).c_str());
	  paintColorScale.drawText(rectText, Qt::AlignCenter, textAge);
	  paintColorScale.restore();
	}
      paintColorScale.restore();

      age -= 250.0;
      if(age < 0)
	break;

      x += widthRect;
    }
  
  paintColorScale.end();

}
// --------------------------------------------------
void
Graphics::drawColorScale(QPainter* painter)
{
  painter->save();
  painter->translate(-width()*0.08, height()*0.900);
  painter->setPen(Qt::white);
  painter->setFont(QFont("Arial", 16));
  QString textAge("continent age (Ga)");
  painter->drawText(rect(), Qt::AlignRight, textAge);
  painter->restore();

  painter->drawPicture(0.630*width(), 0.935*height(), _picColorScale);
}
// --------------------------------------------------
void 
Graphics::drawTime(QPainter* painter, double tMyr, double ageMa)
{  
  unsigned int precision = 2;

  painter->save();
  painter->translate(width()*0.05, height()*0.025);
  painter->setPen(Qt::white);
  painter->setFont(QFont("Arial", 25));
  string prefix = "age = ";
  string suffix = " Ma";
  QString textAge(::doubleToString(ageMa, precision, prefix, suffix).c_str());
  painter->drawText(rect(), Qt::AlignLeft, textAge);
  painter->restore();

  painter->save();
  painter->translate(width()*0.05, height()*0.075);
  painter->setPen(Qt::gray);
  painter->setFont(QFont("Arial", 18));
  prefix = "run = ";
  suffix = " Myr";
  QString textt(::doubleToString(tMyr, precision, prefix, suffix).c_str());
  painter->drawText(rect(), Qt::AlignLeft, textt);
  painter->restore();
}
// --------------------------------------------------
void 
Graphics::drawTemperature(QPainter* painter, double T)
{  
  painter->save();
  painter->translate(-width()*0.05, height()*0.025);
  painter->setPen(Qt::yellow);
  painter->setFont(QFont("Arial", 25));
  string prefix = "T = ";
  string suffix = " K";
  QString textAge(::doubleToString(T, 1, prefix, suffix).c_str());
  painter->drawText(rect(), Qt::AlignRight, textAge);
  painter->restore();
}

// --------------------------------------------------
void 
Graphics::drawTemperatureModif(QPainter* painter, double T, double ageMyr)
{  
  painter->save();
  painter->translate(-width()*0.05, height()*0.025);
  painter->setPen(Qt::yellow);
  painter->setFont(QFont("Arial", 25));
  string prefix = "T = ";
  string suffix = " C";
  double Tmod = T - 273.0;
  QString textAge(::doubleToString(Tmod, 1, prefix, suffix).c_str());
  painter->drawText(rect(), Qt::AlignRight, textAge);
  painter->restore();
}


// --------------------------------------------------
void
Graphics::drawFlux(QPainter* painter, double Q)
{  
  painter->save();
  painter->translate(-width()*0.05, height()*0.075);
  painter->setPen(Qt::cyan);
  painter->setFont(QFont("Arial", 22));
  string prefix = "Q = ";
  string suffix = " TW";
  QString textAge(::doubleToString(Q, 1, prefix, suffix).c_str());
  painter->drawText(rect(), Qt::AlignRight, textAge);
  painter->restore();
}
// --------------------------------------------------
void
Graphics::drawSurfaceRatio(QPainter* painter, double S)
{  
  painter->save();
  painter->translate(+width()*0.025, height()*0.940);
  painter->setPen(Qt::white);
  painter->setFont(QFont("Arial", 18));
  string prefix = "continent ratio = ";
  string suffix = " %";
  S *= 100.0;
  QString textS(::doubleToString(S, 1, prefix, suffix).c_str());
  painter->drawText(rect(), Qt::AlignLeft, textS);
  painter->restore();
}
// --------------------------------------------------
void
Graphics::saveToPNG(const char* filename)
{
  QImage img(width(), height(), QImage::Format_RGB16);
  QPainter painter(&img);
  painter.setRenderHint(QPainter::Antialiasing);
  render(&painter);
  painter.end();
  img.save(QString(filename));

  //    cerr << "Saving " << filename << endl;
}
// --------------------------------------------------
#undef ST_CELL_RMIN
#undef ST_CELL_RMAX

#undef CELL_COLD_R
#undef CELL_COLD_G
#undef CELL_COLD_B

#undef CELL_MANTLE_R
#undef CELL_MANTLE_G
#undef CELL_MANTLE_B

#undef CELL_HOT_R
#undef CELL_HOT_G
#undef CELL_HOT_B

#undef ST_HEIGHT
