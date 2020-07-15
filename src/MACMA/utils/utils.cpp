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

#include "MACMA/utils/utils.h"

// avec la fonction "force_rmfile()" pour supprimer un fichier
void 
force_rmfile(char *filename)
{
	FILE *f = fopen(filename, "r") ;
	if (! f)
		return ;
	fclose(f) ;
	chmod(filename, S_IWRITE) ;
	remove(filename);
}    

// Fonction "force_rmdir()" pour supprimer un repertoire
void 
force_rmdir(char *dirname)
{
	DIR *dir;
	struct dirent *ent;
	dir = opendir(dirname) ;
	if (! dir)
	{
		force_rmfile(dirname) ;
		return ;
	} 
	while((ent = readdir(dir)) != NULL)
	{
		if (ent->d_name[0] == '.')
			continue;
		char *filename = (char*)malloc(strlen(dirname) + 1 + strlen(ent->d_name) + 1 + 1);
		strcpy(filename, dirname);
		strcat(filename, "/");
		strcat(filename, ent->d_name); 
		force_rmdir(filename);
		free(filename);
	}
	rmdir(dirname);
}
