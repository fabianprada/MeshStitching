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


#ifndef VECTOR_IO_INCLUDED
#define VECTOR_IO_INCLUDED
#include<cstdlib>
#include<cstdio>
#include<vector>
#include<unordered_map>

template<typename T>
int ReadVector(std::vector<T> & vec, const char * fileName){

	FILE * file;
	file = fopen(fileName, "rb");
	if (!file) { printf("Unable to read %s \n", fileName); return 0; }
	int vecSize;
	fread(&vecSize, sizeof(int), 1, file);
	vec.resize(vecSize);
	fread(&vec[0], sizeof(T), vecSize, file);
	fclose(file);
	return 1;
}

template<typename T>
void WriteVector(const std::vector<T> & vec, const char * fileName){

	FILE * file;
	file = fopen(fileName, "wb");
	int vecSize = (int)vec.size();
	fwrite(&vecSize, sizeof(int), 1, file);
	fwrite(&vec[0], sizeof(T), vecSize, file);
	fclose(file);
}

template<typename T, typename U>
void WriteMap(const std::unordered_map<T,U> & map, const char * prefix){

	std::vector<T> mapKeys(map.size());
	std::vector<U> mapValues(map.size());
	int counter = 0;
	for (auto iter = map.begin(); iter != map.end(); iter++){
		mapKeys[counter] = (*iter).first;
		mapValues[counter] = (*iter).second;
		counter++;
	}

	char keyFileName[256];
	sprintf(keyFileName, "%s.mapkey", prefix);
	char valueFileName[256];
	sprintf(valueFileName, "%s.mapvalue", prefix);

	WriteVector(mapKeys, keyFileName);
	WriteVector(mapValues, valueFileName);
}

template<typename T, typename U>
int ReadMap(std::unordered_map<T, U> & map, const char * prefix){

	std::vector<T> mapKeys;
	std::vector<U> mapValues;

	char keyFileName[256];
	sprintf(keyFileName, "%s.mapkey", prefix);
	char valueFileName[256];
	sprintf(valueFileName, "%s.mapvalue", prefix);

	if (!ReadVector(mapKeys, keyFileName) || !ReadVector(mapValues, valueFileName)){
		return -1;
	}
	
	if (mapKeys.size() != mapValues.size()){
		printf("Keys and values size does not match! \n");
		return -1;
	}

	for (int i = 0; i < mapKeys.size(); i++) map[mapKeys[i]] = mapValues[i];
}

#endif //VECTOR_IO_INCLUDED