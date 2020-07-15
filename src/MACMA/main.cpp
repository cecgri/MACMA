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
//   author: leyaouanq@cervval.com;
//   modified by cecile.grigne@univ-brest.fr (all modifications since 2013)

#include <QApplication>

#include "MACMA/graphics/macma.h"
#include "MACMA/graphics/macmaGui.h"
#include "MACMA/graphics/macmaNonGui.h"

using namespace std;

void 
messageOutput(QtMsgType type, const char *msg)
{
}

int 
main(int argc, char *argv[])
{

  if(argc > 1)
    {
      string input = argv[1];

      MacmaNonGui* macma = new MacmaNonGui();
      macma->run(input);

      return 0;
    }
  else
    {
      qInstallMsgHandler(messageOutput);
  
      QApplication app(argc, argv);
      
      MacmaGui* window = new MacmaGui();
     
      // window->showMaximized();
      window->showNormal();

      return app.exec();
    }

}

