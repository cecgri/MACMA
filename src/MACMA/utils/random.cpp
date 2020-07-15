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

#include "MACMA/utils/random.h"

using namespace std;

Random::Random(unsigned int iseed)
{
  // --- Initialization
  //  srand((unsigned)(time(0)));
  srand(iseed);
}

Random::~Random()
{

}

// --- int
int
Random::nextInt(void)
{
  return rand();
}

int
Random::nextInt(int max)
{
  return (int)((double)rand() / ((double)RAND_MAX + 1) * max);
}

int
Random::nextInt(int min, int max)
{
  if (min>max)
    {
      cerr << "Error Random::nextInt : min > max" << endl;
      return -1;
    }
  else
    { // CG correction
      return min + (int)((double)rand() / ((double)RAND_MAX + 1) * (max-min));
    }
}

int
Random::nextEvenInt(int min, int max)
{
  int n = 1;
  while(((n = nextInt(min, max)) % 2) != 0);
  return n;
}

int
Random::nextOddInt(int min, int max)
{
  int n = 0;
  while(((n = nextInt(min, max)) % 2) == 0);
  return n;
}

// --- float
float
Random::nextFloat(void)
{
  return rand()/((float)(RAND_MAX)+1.0);
}

float
Random::nextFloat(float max)
{
  return nextFloat()*max;
}

float
Random::nextFloat(float min, float max)
{
  if (min>max)
    {
      return nextFloat()*(min-max)+max;
    }
  else
    {
      return nextFloat()*(max-min)+min;
    }
}

// --- double
double
Random::nextDouble(void)
{
  return rand()/((double)(RAND_MAX)+1);
}

double
Random::nextDouble(double max)
{
  return nextDouble()*max;
}

double
Random::nextDouble(double min, double max)
{
  if (min>max)
    {
      return nextDouble()*(min-max)+max;
    }
  else
    {
      return nextDouble()*(max-min)+min;
    }
}
