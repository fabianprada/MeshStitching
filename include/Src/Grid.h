/*
Copyright (c) 2017, Michael Kazhdan and Fabian Prada
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of
conditions and the following disclaimer. Redistributions in binary form must reproduce
the above copyright notice, this list of conditions and the following disclaimer
in the documentation and/or other materials provided with the distribution.

Neither the name of the Johns Hopkins University nor the names of its contributors
may be used to endorse or promote products derived from this software without specific
prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO THE IMPLIED WARRANTIES
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE  GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
*/


#ifndef SIMPLE_GRID_INCLUDED
#define SIMPLE_GRID_INCLUDED
#include <functional>
#include <algorithm>
#include <omp.h>

template<class T>
class Grid3D{
public:
	Grid3D(){}
	Grid3D(int p_x_dim, int p_y_dim, int p_z_dim){
		x_dim = p_x_dim;
		y_dim = p_y_dim;
		z_dim = p_z_dim;
		xy_dim = x_dim*y_dim;
		xyz_dim = xy_dim*z_dim;
		data = new T[xyz_dim];
	}
	Grid3D(const Grid3D<T> & p_grid){
		x_dim = p_grid.x_dim;
		y_dim = p_grid.y_dim;
		z_dim = p_grid.z_dim;
		xy_dim = x_dim*y_dim;
		xyz_dim = xy_dim*z_dim;
		data = new T[xyz_dim];
		for (int i = 0; i < xyz_dim; i++) data[i] = p_grid.data[i];
	}
	~Grid3D(){
		delete data;
	}
	void copy(const Grid3D<T> & p_grid){
		x_dim = p_grid.x_dim;
		y_dim = p_grid.y_dim;
		z_dim = p_grid.z_dim;
		xy_dim = x_dim*y_dim;
		xyz_dim = xy_dim*z_dim;
		if (data) delete data;
		data = new T[xyz_dim];
		for (int i = 0; i < xyz_dim; i++) data[i] = p_grid.data[i];
	}
	void resize(int p_x_dim, int p_y_dim, int p_z_dim){
		x_dim = p_x_dim;
		y_dim = p_y_dim;
		z_dim = p_z_dim;
		xy_dim = x_dim*y_dim;
		xyz_dim = xy_dim*z_dim;
		if (data) delete data;
		data = new T[xyz_dim];
	}

	int Read(const char * filename);
	void Save(const char * filename) const;

	T& operator () (int x, int y, int z){ return data[xy_dim*z + x_dim*y + x]; }
	const T& operator() (int x, int y, int z) const { return data[xy_dim*z + x_dim*y + x]; }
	T * data = 0;
	int x_dim, y_dim, z_dim, xy_dim, xyz_dim;
};

template<class T>
int Grid3D<T>::Read(const char * filename){
	FILE * file;
	file = fopen(filename, "rb");
	if (!file) { printf("Unable to read %s \n", filename); return 0; }
	int p_x_dim, p_y_dim, p_z_dim;
	fread(&p_x_dim, sizeof(int), 1, file);
	fread(&p_y_dim, sizeof(int), 1, file);
	fread(&p_z_dim, sizeof(int), 1, file);
	resize(p_x_dim, p_y_dim, p_z_dim);
	fread(&data[0], sizeof(T), xyz_dim, file);
	fclose(file);
	return 1;
}

template<class T>
void Grid3D<T>::Save(const char * filename) const{
	FILE * file;
	file = fopen(filename, "wb");
	fwrite(&x_dim, sizeof(int), 1, file);
	fwrite(&y_dim, sizeof(int), 1, file);
	fwrite(&z_dim, sizeof(int), 1, file);
	fwrite(&data[0], sizeof(T), xyz_dim, file);
	fclose(file);
}
#endif// SIMPLE_GRID_INCLUDED