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

#ifndef RANDOM_H_
#define RANDOM_H_

#include <iostream>
#include <stdlib.h>
#include <time.h>

class Random
{
	public:
		Random(unsigned int);
		~Random();

		//generates a pseudo-random integer between 0 and 32767
		int nextInt(void);

		//generates a pseudo-random integer between 0 and max
		int nextInt(int max);

		//generates a pseudo-random integer between min and max
		int nextInt(int min, int max);

		//generates a pseudo-random even integer between min and max
		int nextEvenInt(int min, int max);

		//generates a pseudo-random odd integer between min and max
		int nextOddInt(int min, int max);

		//generates a pseudo-random float between 0.0 and 0.999...
		float nextFloat(void);

		//generates a pseudo-random float between 0.0 and max
		float nextFloat(float max);

		//generates a pseudo-random float between min and max
		float nextFloat(float min, float max);

		//generates a pseudo-random double between 0.0 and 0.999...
		double nextDouble(void);

		//generates a pseudo-random double between 0.0 and max
		double nextDouble(double max);

		//generates a pseudo-random double between min and max
		double nextDouble(double min, double max);
};


#endif /* RANDOM_H_ */
