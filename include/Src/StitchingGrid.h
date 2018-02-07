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


//#include "SimpleMesh.h"
#include "BijectiveMap.h"
#include "VectorIO.h"
#include "Grid.h"
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <set>

//Voxel elements indexig
//Triangle elements indexng

class ChainEdge{
public:
	ChainEdge(){
		sourceVertex = targetVertex = -1;
		previousEdge = nextEdge = -1;
		sourceVertexType = targetVertexType = -1;
		targetVertexOrthogonalEdge = sourceVertexOrthogonalEdge = -1;
	}
	unsigned long long sourceVertex;
	unsigned long long targetVertex;
	unsigned int previousEdge;
	unsigned int nextEdge;
	unsigned int targetVertexOrthogonalEdge;
	unsigned int sourceVertexOrthogonalEdge;
	unsigned char sourceVertexType;
	unsigned char targetVertexType;
	unsigned long faceKey;
	unsigned int chainIndex;
};


class TriangleMap{
public:
	std::vector<std::pair<unsigned long long, unsigned long long>> interiorEdges;
	std::vector<std::pair<double, unsigned long long>> intersectionTime;
};

class FaceMap{
public:
	std::unordered_map<unsigned long long, unsigned long long> forwardVertex;
	std::unordered_map<unsigned long long, unsigned long long> backwardVertex;
	std::unordered_map<unsigned long long, unsigned char> vertexType;
};


class EdgeIntersections{
public:
	std::unordered_set<unsigned long long> vertices;
};


// old_v = new_v*scale +  translation
void CenterMesh(SimpleMesh & mesh, double & scale, Point3D<double> & translation){
	Point3D< double > min_v = mesh.vertices[0];
	Point3D< double > max_v = mesh.vertices[0];
	for (int v = 0; v < mesh.vertices.size(); v++) for (int c = 0; c < 3; c++){
		min_v[c] = std::min<double>(min_v[c], mesh.vertices[v][c]);
		max_v[c] = std::max<double>(max_v[c], mesh.vertices[v][c]);
	}

	Point3D< double > axis_lengths = max_v - min_v;
	double largestAxis = std::max<double>(std::max<double>(axis_lengths[0], axis_lengths[1]), axis_lengths[2])*1.1;


	Point3D< double > center = (max_v + min_v) / 2;
	for (int v = 0; v < mesh.vertices.size(); v++){
		mesh.vertices[v] = ((mesh.vertices[v] - center) / largestAxis) + Point3D< double >(0.5, 0.5, 0.5);
	}

	scale = largestAxis;
	translation = center - Point3D< double >(0.5, 0.5, 0.5)*largestAxis;
}

void CenterMeshes(SimpleMesh & sourceMesh, SimpleMesh & targetMesh, double & scale, Point3D<double> & translation){
	Point3D< double > min_v = sourceMesh.vertices[0];
	Point3D< double > max_v = sourceMesh.vertices[0];
	for (int v = 0; v < sourceMesh.vertices.size(); v++) for (int c = 0; c < 3; c++){
		min_v[c] = std::min<double>(min_v[c], sourceMesh.vertices[v][c]);
		max_v[c] = std::max<double>(max_v[c], sourceMesh.vertices[v][c]);
	}

	for (int v = 0; v < targetMesh.vertices.size(); v++) for (int c = 0; c < 3; c++){
		min_v[c] = std::min<double>(min_v[c], targetMesh.vertices[v][c]);
		max_v[c] = std::max<double>(max_v[c], targetMesh.vertices[v][c]);
	}

	Point3D< double > axis_lengths = max_v - min_v;
	double largestAxis = std::max<double>(std::max<double>(axis_lengths[0], axis_lengths[1]), axis_lengths[2])*1.1;


	Point3D< double > center = (max_v + min_v) / 2;
	for (int v = 0; v < sourceMesh.vertices.size(); v++){
		sourceMesh.vertices[v] = ((sourceMesh.vertices[v] - center) / largestAxis) + Point3D< double >(0.5, 0.5, 0.5);
	}
	for (int v = 0; v < targetMesh.vertices.size(); v++){
		targetMesh.vertices[v] = ((targetMesh.vertices[v] - center) / largestAxis) + Point3D< double >(0.5, 0.5, 0.5);
	}

	scale = largestAxis;
	translation = center - Point3D< double >(0.5, 0.5, 0.5)*largestAxis;
}


unsigned long long SetIntersectionKey(const unsigned long tElement, const unsigned long vElement, const unsigned long intersectionType){
	return ((((unsigned long long)intersectionType << 63) & 0x8000000000000000) | (((unsigned long long)tElement << 32) & 0x7FFFFFFF00000000) | ((unsigned long long)vElement & 0x00000000FFFFFFFF));
}

void GetIntersectionKey(unsigned long long key, unsigned long & tElement, unsigned long & vElement, unsigned long & intersectionType){
	vElement = static_cast<unsigned long>(key & 0x00000000FFFFFFFF);
	tElement = static_cast<unsigned long>((key >> 32) & 0x000000007FFFFFFF);
	intersectionType = static_cast<unsigned long>((key >> 63) & 0x0000000000000001);
}

unsigned long long SetMeshEdgeKey(const unsigned long i0, const unsigned long i1){
	return ((((unsigned long long)i0 << 32) & 0xFFFFFFFF00000000) | ((unsigned long long)i1 & 0x00000000FFFFFFFF));
}

void GetMeshEdgeIndices(unsigned long long key, unsigned long & i0, unsigned long & i1){
	i1 = static_cast<unsigned long>(key & 0x00000000FFFFFFFF);
	i0 = static_cast<unsigned long>((key >> 32) & 0x00000000FFFFFFFF);
}

unsigned long long InvertMeshEdgeKey(unsigned long long key){
	unsigned long i0, i1;
	GetMeshEdgeIndices(key, i0, i1);
	return SetMeshEdgeKey(i1, i0);
}

void SetMeshEdgeKeyMap(std::unordered_map<unsigned long long, unsigned long > & inverseEdgeKeyMap, std::vector<unsigned long long> & edgeKeyMap, const SimpleMesh & mesh){
	unsigned long edgeCounter = 0;

	for (int i = 0; i < mesh.triangles.size(); i++){
		for (int j = 0; j < 3; j++){
			int v[2] = { mesh.triangles[i][j], mesh.triangles[i][(j + 1) % 3] };
			if (v[0]> v[1]) std::swap(v[0], v[1]);
			unsigned long long edgeKey = SetMeshEdgeKey(v[0], v[1]);
			if (inverseEdgeKeyMap.find(edgeKey) == inverseEdgeKeyMap.end()){
				edgeKeyMap.push_back(edgeKey);
				inverseEdgeKeyMap[edgeKey] = edgeCounter;
				edgeCounter++;
			}
		}
	}
}

int SetMeshEdgeKeyTriangleMap(std::unordered_map<unsigned long long, int > & edgeKeyTriangleMap, const SimpleMesh & mesh){
	
	for (int i = 0; i < mesh.triangles.size(); i++){
		for (int j = 0; j < 3; j++){
			unsigned long long edgeKey = SetMeshEdgeKey(mesh.triangles[i][j], mesh.triangles[i][(j + 1) % 3]);
			if (edgeKeyTriangleMap.find(edgeKey) == edgeKeyTriangleMap.end()){
				edgeKeyTriangleMap[edgeKey] = i;
			}
			else{
				printf("Edge adjacent to multiple triangles! \n");
				return 0;
			}
		}
	}
	return 1;
}


unsigned long SetVoxelElementKey(const unsigned long axialIndices[3], const unsigned long cIndex){
	return	(((cIndex << 30) & 0xC0000000) | ((axialIndices[0] << 20) & 0x3FF00000) | ((axialIndices[1] << 10) & 0x000FFC00) | (axialIndices[2] & 0x000003FF));
}

unsigned long SetPointVoxelKey(const Point3D<double> p){
	unsigned long axialIndices[3] = { static_cast<unsigned long>(floor(p[0])), static_cast<unsigned long>(floor(p[1])), static_cast<unsigned long>(floor(p[2])) };
	return SetVoxelElementKey(axialIndices, 0);
}

void GetVoxelElementKeyIndices(unsigned long  key, unsigned long axialIndices[3], unsigned long & cIndex){
	cIndex = static_cast<unsigned long>((key >> 30) & 0x00000003);
	axialIndices[0] = static_cast<unsigned long>((key >> 20) & 0x000003FF);
	axialIndices[1] = static_cast<unsigned long>((key >> 10) & 0x000003FF);
	axialIndices[2] = static_cast<unsigned long>((key & 0x000003FF));
}

void PrintIntersectionInfo(unsigned long long key, const std::vector<unsigned long long> & edgeKeyMap, const std::vector<TriangleIndex> & triangles){
	unsigned long vElement = static_cast<unsigned long>(key & 0x00000000FFFFFFFF);
	unsigned long tElement = static_cast<unsigned long>((key >> 32) & 0x000000007FFFFFFF);
	unsigned long intersectionType = static_cast<unsigned long>((key >> 63) & 0x0000000000000001);

	if (intersectionType == 0){
		//Triangle Edge - Voxel Face
		unsigned long long meshEdgeKey = edgeKeyMap[tElement];
		unsigned long v0, v1;
		GetMeshEdgeIndices(meshEdgeKey, v0, v1);
		unsigned long axialIndices[3];
		unsigned long cIndex;
		GetVoxelElementKeyIndices(vElement, axialIndices, cIndex);
		printf("Mesh edge (%lu,%lu) with voxel face (%lu,%lu,%lu,%lu) \n", v0, v1, axialIndices[0], axialIndices[1], axialIndices[2], cIndex);
	}
	else if (intersectionType == 1){
		//Voxel Edge - Triangle Face
		unsigned long axialIndices[3];
		unsigned long cIndex;
		GetVoxelElementKeyIndices(vElement, axialIndices, cIndex);
		printf("Mesh triangle %lu = (%lu,%lu,%lu) with voxel edge (%lu,%lu,%lu,%lu) \n", tElement, triangles[tElement][0], triangles[tElement][1], triangles[tElement][2], axialIndices[0], axialIndices[1], axialIndices[2], cIndex);
	}
	else{
		printf("Unidentified intersection! \n");
	}
}


int TriangleVoxelFaceIntersection(const SquareMatrix<double, 2> & inverse_metric, const Point3D<double> & normal, const double level, const SimpleMesh & mesh, const unsigned long tIndex, std::unordered_map<unsigned long long, unsigned long > & meshInverseEdgeKeyMap, const std::vector<unsigned long long> & edgeKeyMap,
	const unsigned long axialIndices[3], const unsigned long cIndex, const int res,
	std::unordered_map<unsigned long long, unsigned long long> & adjacentVertexMap, std::unordered_map<unsigned long long, SamplePoint> & vertexSample, std::unordered_map<unsigned long, EdgeIntersections> & edgeIntersections){

	//Fast elimination
	double voxelWidth = 1.0 / double(res);
	Point3D<double> faceCorners[4];
	faceCorners[0][cIndex] = double(axialIndices[cIndex]) * voxelWidth;
	faceCorners[0][(cIndex + 1) % 3] = double(axialIndices[(cIndex + 1) % 3]) * voxelWidth;
	faceCorners[0][(cIndex + 2) % 3] = double(axialIndices[(cIndex + 2) % 3]) * voxelWidth;

	faceCorners[1][cIndex] = double(axialIndices[cIndex]) * voxelWidth;
	faceCorners[1][(cIndex + 1) % 3] = double(axialIndices[(cIndex + 1) % 3] + 1) * voxelWidth;
	faceCorners[1][(cIndex + 2) % 3] = double(axialIndices[(cIndex + 2) % 3]) * voxelWidth;

	faceCorners[2][cIndex] = double(axialIndices[cIndex]) * voxelWidth;
	faceCorners[2][(cIndex + 1) % 3] = double(axialIndices[(cIndex + 1) % 3] + 1) * voxelWidth;
	faceCorners[2][(cIndex + 2) % 3] = double(axialIndices[(cIndex + 2) % 3] + 1) * voxelWidth;

	faceCorners[3][cIndex] = double(axialIndices[cIndex]) * voxelWidth;
	faceCorners[3][(cIndex + 1) % 3] = double(axialIndices[(cIndex + 1) % 3]) * voxelWidth;
	faceCorners[3][(cIndex + 2) % 3] = double(axialIndices[(cIndex + 2) % 3] + 1) * voxelWidth;

	double cornerOffset[4];
	for (int i = 0; i < 4; i++) cornerOffset[i] = Point3D<double>::Dot(faceCorners[i], normal) - level;

	if ((cornerOffset[0] > 0 && cornerOffset[1] > 0 && cornerOffset[2] > 0 && cornerOffset[3] > 0) || (cornerOffset[0] < 0 && cornerOffset[1] < 0 && cornerOffset[2] < 0 && cornerOffset[3] < 0)){
		return 0;
	}

	//Compute edge intersections

	//Triangle edges with voxel faces
	unsigned long long intersectionsIndices[2] = { static_cast<unsigned long long>(-1), static_cast<unsigned long long>(-1) };
	Point3D<double> intersectionsPoints[2];
	int intersectionsCounter = 0;
	for (int j = 0; j < 3; j++){
		int v[2] = { mesh.triangles[tIndex][j], mesh.triangles[tIndex][(j + 1) % 3] };
		bool reverseOrder = false;
		if (v[0]> v[1]) {
			std::swap(v[0], v[1]);
			reverseOrder = true;
		}
		double axialStart = mesh.vertices[v[0]][cIndex];
		double axialEnd = mesh.vertices[v[1]][cIndex];
		double alpha = (double(axialIndices[cIndex]) * voxelWidth - axialStart) / (axialEnd - axialStart);
		if (alpha > 0.0  && alpha < 1.0){
			Point3D<double> isect = mesh.vertices[v[0]] * (1.0 - alpha) + mesh.vertices[v[1]] * alpha;
			if (isect[(cIndex + 1) % 3] > double(axialIndices[(cIndex + 1) % 3]) * voxelWidth && isect[(cIndex + 1) % 3] <  double(axialIndices[(cIndex + 1) % 3] + 1) * voxelWidth &&
				isect[(cIndex + 2) % 3] >  double(axialIndices[(cIndex + 2) % 3]) * voxelWidth && isect[(cIndex + 2) % 3] < double(axialIndices[(cIndex + 2) % 3] + 1) * voxelWidth){
				if (intersectionsCounter < 2){
					unsigned long long meshEdgeKey = SetMeshEdgeKey(v[0], v[1]);
					unsigned long voxelFaceKey = SetVoxelElementKey(axialIndices, cIndex);
					unsigned long tElement = meshInverseEdgeKeyMap[meshEdgeKey];
					unsigned long long intersectionKey = SetIntersectionKey(tElement, voxelFaceKey, static_cast<unsigned long>(0));
					intersectionsIndices[intersectionsCounter] = intersectionKey;
					intersectionsPoints[intersectionsCounter] = isect;
					intersectionsCounter++;

					//Mod Nov 11/2016
					Point3D<double> baricentric3D;
					baricentric3D[j] = (1.0 - alpha);
					baricentric3D[(j + 1) % 3] = alpha;
					if (reverseOrder) std::swap(baricentric3D[j], baricentric3D[(j + 1) % 3]);
					SamplePoint sample;
					sample.tIndex = tIndex;
					sample.baricentricCoord = Point2D<double>(baricentric3D[1], baricentric3D[2]);
					vertexSample[intersectionKey] = sample;
				}
				else{
					printf("More than two intersection between triangle and voxel grid face \n");
				}
			}
		}
	}

	//Voxel edges with triangle

	Point3D<double> d[2] = { mesh.vertices[mesh.triangles[tIndex][1]] - mesh.vertices[mesh.triangles[tIndex][0]], mesh.vertices[mesh.triangles[tIndex][2]] - mesh.vertices[mesh.triangles[tIndex][0]] };

	for (int j = 0; j < 4; j++){
		double alpha = cornerOffset[j] / (cornerOffset[j] - cornerOffset[(j + 1) % 4]);
		if (alpha > 0.0  && alpha < 1.0){
			Point3D<double> isect = faceCorners[j] * (1.0 - alpha) + faceCorners[(j + 1) % 4] * alpha;
			Point3D<double> projection = isect - mesh.vertices[mesh.triangles[tIndex][0]];
			Point2D<double> projection_coeffs = Point2D<double>(Point3D<double>::Dot(projection, d[0]), Point3D<double>::Dot(projection, d[1]));
			Point2D<double> baricentricCoord = inverse_metric*projection_coeffs;
			if (baricentricCoord[0] > 0.0 && baricentricCoord[1] > 0.0 && (baricentricCoord[0] + baricentricCoord[1]) < 1.0){
				if (intersectionsCounter < 2){

					unsigned long edgeAxialIndices[3] = { axialIndices[0], axialIndices[1], axialIndices[2] };
					unsigned long voxelEdgeKey;
					if (j == 0){
						voxelEdgeKey = SetVoxelElementKey(edgeAxialIndices, (cIndex + 1) % 3);
					}
					if (j == 1){
						edgeAxialIndices[(cIndex + 1) % 3]++;
						voxelEdgeKey = SetVoxelElementKey(edgeAxialIndices, (cIndex + 2) % 3);
					}
					if (j == 2){
						edgeAxialIndices[(cIndex + 2) % 3]++;
						voxelEdgeKey = SetVoxelElementKey(edgeAxialIndices, (cIndex + 1) % 3);
					}
					if (j == 3){
						voxelEdgeKey = SetVoxelElementKey(edgeAxialIndices, (cIndex + 2) % 3);
					}
					unsigned long long intersectionKey = SetIntersectionKey(tIndex, voxelEdgeKey, static_cast<unsigned long>(1));
					intersectionsIndices[intersectionsCounter] = intersectionKey;
					edgeIntersections[voxelEdgeKey].vertices.insert(intersectionKey);
					intersectionsPoints[intersectionsCounter] = isect;

					SamplePoint sample;
					sample.tIndex = tIndex;
					sample.baricentricCoord = baricentricCoord;
					vertexSample[intersectionKey] = sample;

					intersectionsCounter++;
				}
				else{
					printf("More than two intersection between triangle and voxel grid face%d! \n");
				}
			}
		}
	}
	if (intersectionsCounter == 2){
		//Check oriention
		Point3D<double> intersectionRay = intersectionsPoints[1] - intersectionsPoints[0];
		Point3D<double> axialRay(0, 0, 0);
		axialRay[cIndex] = 1.0;
		Point3D<double> rotIntersectionRay = Point3D<double>::CrossProduct(axialRay, intersectionRay);
		double orientation = Point3D<double>::Dot(rotIntersectionRay, normal);
		if (orientation > 0){
			std::swap(intersectionsIndices[0], intersectionsIndices[1]);
		}


		if (adjacentVertexMap.find(intersectionsIndices[0]) != adjacentVertexMap.end()){
			printf("ERROR : Multiple mappings from a single vertex %llu!", intersectionsIndices[0]);
			return -1;
		}
		else{
			adjacentVertexMap[intersectionsIndices[0]] = intersectionsIndices[1];
		}
		return 2;
	}
	else if (intersectionsCounter == 1){
		printf("ERROR : Unexpected single intersection \n");
		for (int i = 0; i < 4; i++){
			printf("%g %g %g \n", faceCorners[i][0], faceCorners[i][1], faceCorners[i][2]);
		}
		for (int i = 0; i < 3; i++){
			printf("%g %g %g \n", mesh.vertices[mesh.triangles[tIndex][i]][0], mesh.vertices[mesh.triangles[tIndex][i]][1], mesh.vertices[mesh.triangles[tIndex][i]][2]);
		}
		return -1;
	}
	return 0;
}

class TriangleEdgeVoxelIntersection{
public:
	TriangleEdgeVoxelIntersection(){
		intersectionTime[0] = intersectionTime[1] = -1.0;
	}
	unsigned long long intersectionKey[2];
	double intersectionTime[2];
};

int TriangleVoxelFaceIntersection(const SquareMatrix<double, 2> & inverse_metric, const Point3D<double> & normal, const double level, const SimpleMesh & mesh, const unsigned long tIndex, std::unordered_map<unsigned long long, unsigned long > & meshInverseEdgeKeyMap, const std::vector<unsigned long long> & edgeKeyMap,
	const unsigned long axialIndices[3], const unsigned long cIndex, const int res,
	std::unordered_map<unsigned long long, unsigned long long> & adjacentVertexMap, TriangleEdgeVoxelIntersection triangleEdgeIntersections[3], std::unordered_map<unsigned long long, Point3D<double>> & vertexPos, bool positiveOrientation){

	//Fast elimination
	double voxelWidth = 1.0 / double(res);
	Point3D<double> faceCorners[4];
	faceCorners[0][cIndex] = double(axialIndices[cIndex]) * voxelWidth;
	faceCorners[0][(cIndex + 1) % 3] = double(axialIndices[(cIndex + 1) % 3]) * voxelWidth;
	faceCorners[0][(cIndex + 2) % 3] = double(axialIndices[(cIndex + 2) % 3]) * voxelWidth;

	faceCorners[1][cIndex] = double(axialIndices[cIndex]) * voxelWidth;
	faceCorners[1][(cIndex + 1) % 3] = double(axialIndices[(cIndex + 1) % 3] + 1) * voxelWidth;
	faceCorners[1][(cIndex + 2) % 3] = double(axialIndices[(cIndex + 2) % 3]) * voxelWidth;

	faceCorners[2][cIndex] = double(axialIndices[cIndex]) * voxelWidth;
	faceCorners[2][(cIndex + 1) % 3] = double(axialIndices[(cIndex + 1) % 3] + 1) * voxelWidth;
	faceCorners[2][(cIndex + 2) % 3] = double(axialIndices[(cIndex + 2) % 3] + 1) * voxelWidth;

	faceCorners[3][cIndex] = double(axialIndices[cIndex]) * voxelWidth;
	faceCorners[3][(cIndex + 1) % 3] = double(axialIndices[(cIndex + 1) % 3]) * voxelWidth;
	faceCorners[3][(cIndex + 2) % 3] = double(axialIndices[(cIndex + 2) % 3] + 1) * voxelWidth;

	double cornerOffset[4];
	for (int i = 0; i < 4; i++) cornerOffset[i] = Point3D<double>::Dot(faceCorners[i], normal) - level;

	if ((cornerOffset[0] > 0 && cornerOffset[1] > 0 && cornerOffset[2] > 0 && cornerOffset[3] > 0) || (cornerOffset[0] < 0 && cornerOffset[1] < 0 && cornerOffset[2] < 0 && cornerOffset[3] < 0)){
		return 0;
	}

	//Compute edge intersections

	//Triangle edges with voxel faces
	unsigned long long intersectionsIndices[2] = { static_cast<unsigned long long>(-1), static_cast<unsigned long long>(-1) };
	Point3D<double> intersectionsPoints[2];
	int intersectionsCounter = 0;
	for (int j = 0; j < 3; j++){
		int v[2] = { mesh.triangles[tIndex][j], mesh.triangles[tIndex][(j + 1) % 3] };
		bool preservedOrientation = true;
		if (v[0]> v[1]){
			std::swap(v[0], v[1]);
			preservedOrientation = false;
		}

		double axialStart = mesh.vertices[v[0]][cIndex];
		double axialEnd = mesh.vertices[v[1]][cIndex];
		double alpha = (double(axialIndices[cIndex]) * voxelWidth - axialStart) / (axialEnd - axialStart);
		if (alpha > 0.0  && alpha < 1.0){
			Point3D<double> isect = mesh.vertices[v[0]] * (1.0 - alpha) + mesh.vertices[v[1]] * alpha;
			if (isect[(cIndex + 1) % 3] > double(axialIndices[(cIndex + 1) % 3]) * voxelWidth && isect[(cIndex + 1) % 3] <  double(axialIndices[(cIndex + 1) % 3] + 1) * voxelWidth &&
				isect[(cIndex + 2) % 3] >  double(axialIndices[(cIndex + 2) % 3]) * voxelWidth && isect[(cIndex + 2) % 3] < double(axialIndices[(cIndex + 2) % 3] + 1) * voxelWidth){
				if (intersectionsCounter < 2){
					unsigned long long meshEdgeKey = SetMeshEdgeKey(v[0], v[1]);
					unsigned long voxelFaceKey = SetVoxelElementKey(axialIndices, cIndex);
					unsigned long tElement = meshInverseEdgeKeyMap[meshEdgeKey];
					unsigned long long intersectionKey = SetIntersectionKey(tElement, voxelFaceKey, static_cast<unsigned long>(0));
					intersectionsIndices[intersectionsCounter] = intersectionKey;
					intersectionsPoints[intersectionsCounter] = isect;
					intersectionsCounter++;

					vertexPos[intersectionKey] = isect;

					TriangleEdgeVoxelIntersection & triangleEdgeIntersection = triangleEdgeIntersections[j];
					double isectTime = preservedOrientation ? alpha : 1.0 - alpha;
					if (triangleEdgeIntersection.intersectionTime[0] == -1){
						triangleEdgeIntersection.intersectionTime[0] = isectTime;
						triangleEdgeIntersection.intersectionKey[0] = intersectionKey;
					}
					else if (triangleEdgeIntersection.intersectionTime[1] == -1){
						triangleEdgeIntersection.intersectionTime[1] = isectTime;
						triangleEdgeIntersection.intersectionKey[1] = intersectionKey;
					}
					else{
						printf("Unexpected condition: more than two intersection on an edge! \n");
						return 0;
					}
				}
				else{
					printf("More than two intersection between triangle and voxel grid face \n");
				}
			}
		}
	}

	//Voxel edges with triangle

	Point3D<double> d[2] = { mesh.vertices[mesh.triangles[tIndex][1]] - mesh.vertices[mesh.triangles[tIndex][0]], mesh.vertices[mesh.triangles[tIndex][2]] - mesh.vertices[mesh.triangles[tIndex][0]] };

	for (int j = 0; j < 4; j++){
		double alpha = cornerOffset[j] / (cornerOffset[j] - cornerOffset[(j + 1) % 4]);
		if (alpha > 0.0  && alpha < 1.0){
			Point3D<double> isect = faceCorners[j] * (1.0 - alpha) + faceCorners[(j + 1) % 4] * alpha;
			Point3D<double> projection = isect - mesh.vertices[mesh.triangles[tIndex][0]];
			Point2D<double> projection_coeffs = Point2D<double>(Point3D<double>::Dot(projection, d[0]), Point3D<double>::Dot(projection, d[1]));
			Point2D<double> baricentricCoord = inverse_metric*projection_coeffs;
			if (baricentricCoord[0] > 0.0 && baricentricCoord[1] > 0.0 && (baricentricCoord[0] + baricentricCoord[1]) < 1.0){
				if (intersectionsCounter < 2){

					unsigned long edgeAxialIndices[3] = { axialIndices[0], axialIndices[1], axialIndices[2] };
					unsigned long voxelEdgeKey;
					if (j == 0){
						voxelEdgeKey = SetVoxelElementKey(edgeAxialIndices, (cIndex + 1) % 3);
					}
					if (j == 1){
						edgeAxialIndices[(cIndex + 1) % 3]++;
						voxelEdgeKey = SetVoxelElementKey(edgeAxialIndices, (cIndex + 2) % 3);
					}
					if (j == 2){
						edgeAxialIndices[(cIndex + 2) % 3]++;
						voxelEdgeKey = SetVoxelElementKey(edgeAxialIndices, (cIndex + 1) % 3);
					}
					if (j == 3){
						voxelEdgeKey = SetVoxelElementKey(edgeAxialIndices, (cIndex + 2) % 3);
					}
					unsigned long long intersectionKey = SetIntersectionKey(tIndex, voxelEdgeKey, static_cast<unsigned long>(1));
					intersectionsIndices[intersectionsCounter] = intersectionKey;
					intersectionsPoints[intersectionsCounter] = isect;

					vertexPos[intersectionKey] = isect;

					intersectionsCounter++;
				}
				else{
					printf("More than two intersection between triangle and voxel grid face%d! \n");
				}
			}
		}
	}
	if (intersectionsCounter == 2){
		//Check oriention
		Point3D<double> intersectionRay = intersectionsPoints[1] - intersectionsPoints[0];
		Point3D<double> axialRay(0, 0, 0);
		axialRay[cIndex] = positiveOrientation ? 1.0 : -1.0;
		Point3D<double> rotIntersectionRay = Point3D<double>::CrossProduct(axialRay, intersectionRay);
		double orientation = Point3D<double>::Dot(rotIntersectionRay, normal);
		if (orientation > 0){
			std::swap(intersectionsIndices[0], intersectionsIndices[1]);
			//std::swap(intersectionsPoints[0], intersectionsPoints[1]);
		}

		if (adjacentVertexMap.find(intersectionsIndices[0]) != adjacentVertexMap.end()){
			printf("ERROR : Multiple mappings from a single vertex %llu!", intersectionsIndices[0]);
			return -1;
		}
		else{
			adjacentVertexMap[intersectionsIndices[0]] = intersectionsIndices[1];
		}
		return 2;
	}
	else if (intersectionsCounter == 1){
		printf("ERROR : Unexpected single intersection \n");
		for (int i = 0; i < 4; i++){
			printf("%g %g %g \n", faceCorners[i][0], faceCorners[i][1], faceCorners[i][2]);
		}
		for (int i = 0; i < 3; i++){
			printf("%g %g %g \n", mesh.vertices[mesh.triangles[tIndex][i]][0], mesh.vertices[mesh.triangles[tIndex][i]][1], mesh.vertices[mesh.triangles[tIndex][i]][2]);
		}
		return -1;
	}
	return 0;
}

class OrderedIntersection{
public:
	unsigned long long vertex;
	double time;
};

class OrderedIntersectionComparison{
public:
	bool operator()(const OrderedIntersection & i0, const OrderedIntersection & i1){
		if (i0.vertex == i1.vertex){
			return false;
		}
		else return i0.time < i1.time;
	}
};

int TriangleVoxelFaceIntersection(const SquareMatrix<double, 2> & inverse_metric, const Point3D<double> & normal, const double level, const SimpleMesh & mesh, const unsigned long tIndex, std::unordered_map<unsigned long long, unsigned long > & meshInverseEdgeKeyMap, const std::vector<unsigned long long> & edgeKeyMap,
	const unsigned long axialIndices[3], const unsigned long cIndex, const int res,
	std::unordered_map<unsigned long long, unsigned long long> & adjacentVertexMap, TriangleEdgeVoxelIntersection triangleEdgeIntersections[3], std::unordered_map<unsigned long long, Point3D<double>> & vertexPos, bool positiveOrientation, std::unordered_set<unsigned long long> & splitEdgeKeys, std::unordered_map< unsigned long long, std::set<OrderedIntersection, OrderedIntersectionComparison>>& splitEdgeIntersections){

	//Fast elimination
	double voxelWidth = 1.0 / double(res);
	Point3D<double> faceCorners[4];
	faceCorners[0][cIndex] = double(axialIndices[cIndex]) * voxelWidth;
	faceCorners[0][(cIndex + 1) % 3] = double(axialIndices[(cIndex + 1) % 3]) * voxelWidth;
	faceCorners[0][(cIndex + 2) % 3] = double(axialIndices[(cIndex + 2) % 3]) * voxelWidth;

	faceCorners[1][cIndex] = double(axialIndices[cIndex]) * voxelWidth;
	faceCorners[1][(cIndex + 1) % 3] = double(axialIndices[(cIndex + 1) % 3] + 1) * voxelWidth;
	faceCorners[1][(cIndex + 2) % 3] = double(axialIndices[(cIndex + 2) % 3]) * voxelWidth;

	faceCorners[2][cIndex] = double(axialIndices[cIndex]) * voxelWidth;
	faceCorners[2][(cIndex + 1) % 3] = double(axialIndices[(cIndex + 1) % 3] + 1) * voxelWidth;
	faceCorners[2][(cIndex + 2) % 3] = double(axialIndices[(cIndex + 2) % 3] + 1) * voxelWidth;

	faceCorners[3][cIndex] = double(axialIndices[cIndex]) * voxelWidth;
	faceCorners[3][(cIndex + 1) % 3] = double(axialIndices[(cIndex + 1) % 3]) * voxelWidth;
	faceCorners[3][(cIndex + 2) % 3] = double(axialIndices[(cIndex + 2) % 3] + 1) * voxelWidth;

	double cornerOffset[4];
	for (int i = 0; i < 4; i++) cornerOffset[i] = Point3D<double>::Dot(faceCorners[i], normal) - level;

	if ((cornerOffset[0] > 0 && cornerOffset[1] > 0 && cornerOffset[2] > 0 && cornerOffset[3] > 0) || (cornerOffset[0] < 0 && cornerOffset[1] < 0 && cornerOffset[2] < 0 && cornerOffset[3] < 0)){
		return 0;
	}

	//Compute edge intersections

	//Triangle edges with voxel faces
	unsigned long long intersectionsIndices[2] = { static_cast<unsigned long long>(-1), static_cast<unsigned long long>(-1) };
	Point3D<double> intersectionsPoints[2];
	int intersectionsCounter = 0;
	for (int j = 0; j < 3; j++){
		int v[2] = { mesh.triangles[tIndex][j], mesh.triangles[tIndex][(j + 1) % 3] };
		bool preservedOrientation = true;
		if (v[0]> v[1]){
			std::swap(v[0], v[1]);
			preservedOrientation = false;
		}

		double axialStart = mesh.vertices[v[0]][cIndex];
		double axialEnd = mesh.vertices[v[1]][cIndex];
		double alpha = (double(axialIndices[cIndex]) * voxelWidth - axialStart) / (axialEnd - axialStart);
		if (alpha > 0.0  && alpha < 1.0){
			Point3D<double> isect = mesh.vertices[v[0]] * (1.0 - alpha) + mesh.vertices[v[1]] * alpha;
			if (isect[(cIndex + 1) % 3] > double(axialIndices[(cIndex + 1) % 3]) * voxelWidth && isect[(cIndex + 1) % 3] <  double(axialIndices[(cIndex + 1) % 3] + 1) * voxelWidth &&
				isect[(cIndex + 2) % 3] >  double(axialIndices[(cIndex + 2) % 3]) * voxelWidth && isect[(cIndex + 2) % 3] < double(axialIndices[(cIndex + 2) % 3] + 1) * voxelWidth){
				if (intersectionsCounter < 2){
					unsigned long long meshEdgeKey = SetMeshEdgeKey(v[0], v[1]);
					unsigned long voxelFaceKey = SetVoxelElementKey(axialIndices, cIndex);
					unsigned long tElement = meshInverseEdgeKeyMap[meshEdgeKey];
					unsigned long long intersectionKey = SetIntersectionKey(tElement, voxelFaceKey, static_cast<unsigned long>(0));
					intersectionsIndices[intersectionsCounter] = intersectionKey;
					intersectionsPoints[intersectionsCounter] = isect;
					intersectionsCounter++;

					vertexPos[intersectionKey] = isect;

					TriangleEdgeVoxelIntersection & triangleEdgeIntersection = triangleEdgeIntersections[j];
					double isectTime = preservedOrientation ? alpha : 1.0 - alpha;
					if (triangleEdgeIntersection.intersectionTime[0] == -1){
						triangleEdgeIntersection.intersectionTime[0] = isectTime;
						triangleEdgeIntersection.intersectionKey[0] = intersectionKey;
					}
					else if (triangleEdgeIntersection.intersectionTime[1] == -1){
						triangleEdgeIntersection.intersectionTime[1] = isectTime;
						triangleEdgeIntersection.intersectionKey[1] = intersectionKey;
					}
					else{
						printf("Unexpected condition: more than two intersection on an edge! \n");
						return 0;
					}

					if (splitEdgeKeys.find(meshEdgeKey) != splitEdgeKeys.end()){
						OrderedIntersection orderedItersection;
						orderedItersection.time = 1.0 - isectTime;
						orderedItersection.vertex = intersectionKey;
						splitEdgeIntersections[meshEdgeKey].insert(orderedItersection);
					}
				}
				else{
					printf("More than two intersection between triangle and voxel grid face \n");
				}
			}
		}
	}

	//Voxel edges with triangle

	Point3D<double> d[2] = { mesh.vertices[mesh.triangles[tIndex][1]] - mesh.vertices[mesh.triangles[tIndex][0]], mesh.vertices[mesh.triangles[tIndex][2]] - mesh.vertices[mesh.triangles[tIndex][0]] };

	for (int j = 0; j < 4; j++){
		double alpha = cornerOffset[j] / (cornerOffset[j] - cornerOffset[(j + 1) % 4]);
		if (alpha > 0.0  && alpha < 1.0){
			Point3D<double> isect = faceCorners[j] * (1.0 - alpha) + faceCorners[(j + 1) % 4] * alpha;
			Point3D<double> projection = isect - mesh.vertices[mesh.triangles[tIndex][0]];
			Point2D<double> projection_coeffs = Point2D<double>(Point3D<double>::Dot(projection, d[0]), Point3D<double>::Dot(projection, d[1]));
			Point2D<double> baricentricCoord = inverse_metric*projection_coeffs;
			if (baricentricCoord[0] > 0.0 && baricentricCoord[1] > 0.0 && (baricentricCoord[0] + baricentricCoord[1]) < 1.0){
				if (intersectionsCounter < 2){

					unsigned long edgeAxialIndices[3] = { axialIndices[0], axialIndices[1], axialIndices[2] };
					unsigned long voxelEdgeKey;
					if (j == 0){
						voxelEdgeKey = SetVoxelElementKey(edgeAxialIndices, (cIndex + 1) % 3);
					}
					if (j == 1){
						edgeAxialIndices[(cIndex + 1) % 3]++;
						voxelEdgeKey = SetVoxelElementKey(edgeAxialIndices, (cIndex + 2) % 3);
					}
					if (j == 2){
						edgeAxialIndices[(cIndex + 2) % 3]++;
						voxelEdgeKey = SetVoxelElementKey(edgeAxialIndices, (cIndex + 1) % 3);
					}
					if (j == 3){
						voxelEdgeKey = SetVoxelElementKey(edgeAxialIndices, (cIndex + 2) % 3);
					}
					unsigned long long intersectionKey = SetIntersectionKey(tIndex, voxelEdgeKey, static_cast<unsigned long>(1));
					intersectionsIndices[intersectionsCounter] = intersectionKey;
					intersectionsPoints[intersectionsCounter] = isect;

					vertexPos[intersectionKey] = isect;

					intersectionsCounter++;
				}
				else{
					printf("More than two intersection between triangle and voxel grid face%d! \n");
				}
			}
		}
	}
	if (intersectionsCounter == 2){
		//Check oriention
		Point3D<double> intersectionRay = intersectionsPoints[1] - intersectionsPoints[0];
		Point3D<double> axialRay(0, 0, 0);
		axialRay[cIndex] = positiveOrientation ? 1.0 : -1.0;
		Point3D<double> rotIntersectionRay = Point3D<double>::CrossProduct(axialRay, intersectionRay);
		double orientation = Point3D<double>::Dot(rotIntersectionRay, normal);
		if (orientation > 0){
			std::swap(intersectionsIndices[0], intersectionsIndices[1]);
			//std::swap(intersectionsPoints[0], intersectionsPoints[1]);
		}

		if (adjacentVertexMap.find(intersectionsIndices[0]) != adjacentVertexMap.end()){
			printf("ERROR : Multiple mappings from a single vertex %llu!", intersectionsIndices[0]);
			return -1;
		}
		else{
			adjacentVertexMap[intersectionsIndices[0]] = intersectionsIndices[1];
		}
		return 2;
	}
	else if (intersectionsCounter == 1){
		printf("ERROR : Unexpected single intersection \n");
		for (int i = 0; i < 4; i++){
			printf("%g %g %g \n", faceCorners[i][0], faceCorners[i][1], faceCorners[i][2]);
		}
		for (int i = 0; i < 3; i++){
			printf("%g %g %g \n", mesh.vertices[mesh.triangles[tIndex][i]][0], mesh.vertices[mesh.triangles[tIndex][i]][1], mesh.vertices[mesh.triangles[tIndex][i]][2]);
		}
		return -1;
	}
	return 0;
}

static int functionCall = 0;

int TriangleVoxelIntersection(const SquareMatrix<double, 2> & inverse_metric, const Point3D<double> & normal, const double level, const SimpleMesh & mesh, const unsigned long tIndex, std::unordered_map<unsigned long long, unsigned long > & meshInverseEdgeKeyMap, const std::vector<unsigned long long> & edgeKeyMap,
	const unsigned long axialIndices[3], const int res,
	std::unordered_map<unsigned long long, unsigned long long> & adjacentVertexMap, std::unordered_map<unsigned long long, Point3D<double>> & vertexPos, std::unordered_set<unsigned long long> & splitEdgeKeys, std::unordered_map< unsigned long long, std::set<OrderedIntersection, OrderedIntersectionComparison>> & splitEdgeIntersections){

	functionCall++;

	TriangleEdgeVoxelIntersection triangleEdgeIntersections[3];

	//Check Interior Vertices
	double voxelWidth = 1.0 / double(res);
	Point3D<double> voxelMinCorner = Point3D<double>(axialIndices[0], axialIndices[1], axialIndices[2])* voxelWidth;
	Point3D<double> voxelMaxCorner = Point3D<double>(axialIndices[0] + 1, axialIndices[1] + 1, axialIndices[2] + 1)* voxelWidth;

	for (int i = 0; i < 3; i++){
		Point3D<double> vertex = mesh.vertices[mesh.triangles[tIndex][i]];
		bool isInterior = true;
		for (int c = 0; c < 3; c++) if (!(voxelMinCorner[c] < vertex[c] && vertex[c] < voxelMaxCorner[c])) isInterior = false;
		if (isInterior) {
			if (triangleEdgeIntersections[i].intersectionTime[0] == -1){
				triangleEdgeIntersections[i].intersectionTime[0] = 0;
				triangleEdgeIntersections[i].intersectionKey[0] = mesh.triangles[tIndex][i];
			}
			else{
				triangleEdgeIntersections[i].intersectionTime[1] = 0;
				triangleEdgeIntersections[i].intersectionKey[1] = mesh.triangles[tIndex][i];
			}

			if (triangleEdgeIntersections[(i + 2) % 3].intersectionTime[0] == -1){
				triangleEdgeIntersections[(i + 2) % 3].intersectionTime[0] = 1.0;
				triangleEdgeIntersections[(i + 2) % 3].intersectionKey[0] = mesh.triangles[tIndex][i];
			}
			else{
				triangleEdgeIntersections[(i + 2) % 3].intersectionTime[1] = 1.0;
				triangleEdgeIntersections[(i + 2) % 3].intersectionKey[1] = mesh.triangles[tIndex][i];
			}
		}
	}

	for (unsigned long dir = 0; dir < 3; dir++){
		for (int offset = 0; offset < 2; offset++){
			unsigned long faceAxialIndices[3] = { axialIndices[0], axialIndices[1], axialIndices[2] };
			if (offset == 1) faceAxialIndices[dir]++;
			TriangleVoxelFaceIntersection(inverse_metric, normal, level, mesh, tIndex, meshInverseEdgeKeyMap, edgeKeyMap, faceAxialIndices, dir, res, adjacentVertexMap, triangleEdgeIntersections, vertexPos, offset == 0, splitEdgeKeys,splitEdgeIntersections);
		}
	}


	for (int i = 0; i < 3; i++){
		if (triangleEdgeIntersections[i].intersectionTime[0] != -1){
			if (triangleEdgeIntersections[i].intersectionTime[1] == -1){
				printf("Triangle edge only have a single intersection! \n");
				return 0; 
			}
			else{
				if (triangleEdgeIntersections[i].intersectionTime[1] < triangleEdgeIntersections[i].intersectionTime[0]){
					std::swap(triangleEdgeIntersections[i].intersectionTime[0], triangleEdgeIntersections[i].intersectionTime[1]);
					std::swap(triangleEdgeIntersections[i].intersectionKey[0], triangleEdgeIntersections[i].intersectionKey[1]);
				}
				adjacentVertexMap[triangleEdgeIntersections[i].intersectionKey[0]] = triangleEdgeIntersections[i].intersectionKey[1];
			}
		}
	}
}


int TriangleVoxelIntersection(const SquareMatrix<double, 2> & inverse_metric, const Point3D<double> & normal, const double level, const SimpleMesh & mesh, const unsigned long tIndex, std::unordered_map<unsigned long long, unsigned long > & meshInverseEdgeKeyMap, const std::vector<unsigned long long> & edgeKeyMap,
	const unsigned long axialIndices[3], const int res,
	std::unordered_map<unsigned long long, unsigned long long> & adjacentVertexMap, std::unordered_map<unsigned long long, Point3D<double>> & vertexPos){

	functionCall++;

	TriangleEdgeVoxelIntersection triangleEdgeIntersections[3];

	//Check Interior Vertices
	double voxelWidth = 1.0 / double(res);
	Point3D<double> voxelMinCorner = Point3D<double>(axialIndices[0], axialIndices[1], axialIndices[2])* voxelWidth;
	Point3D<double> voxelMaxCorner = Point3D<double>(axialIndices[0] + 1, axialIndices[1] + 1, axialIndices[2] + 1)* voxelWidth;

	for (int i = 0; i < 3; i++){
		Point3D<double> vertex = mesh.vertices[mesh.triangles[tIndex][i]];
		bool isInterior = true;
		for (int c = 0; c < 3; c++) if (!(voxelMinCorner[c] < vertex[c] && vertex[c] < voxelMaxCorner[c])) isInterior = false;
		if (isInterior) {
			if (triangleEdgeIntersections[i].intersectionTime[0] == -1){
				triangleEdgeIntersections[i].intersectionTime[0] = 0;
				triangleEdgeIntersections[i].intersectionKey[0] = mesh.triangles[tIndex][i];
			}
			else{
				triangleEdgeIntersections[i].intersectionTime[1] = 0;
				triangleEdgeIntersections[i].intersectionKey[1] = mesh.triangles[tIndex][i];
			}

			if (triangleEdgeIntersections[(i + 2) % 3].intersectionTime[0] == -1){
				triangleEdgeIntersections[(i + 2) % 3].intersectionTime[0] = 1.0;
				triangleEdgeIntersections[(i + 2) % 3].intersectionKey[0] = mesh.triangles[tIndex][i];
			}
			else{
				triangleEdgeIntersections[(i + 2) % 3].intersectionTime[1] = 1.0;
				triangleEdgeIntersections[(i + 2) % 3].intersectionKey[1] = mesh.triangles[tIndex][i];
			}
		}
	}

	for (unsigned long dir = 0; dir < 3; dir++){
		for (int offset = 0; offset < 2; offset++){
			unsigned long faceAxialIndices[3] = { axialIndices[0], axialIndices[1], axialIndices[2] };
			if (offset == 1) faceAxialIndices[dir]++;
			TriangleVoxelFaceIntersection(inverse_metric, normal, level, mesh, tIndex, meshInverseEdgeKeyMap, edgeKeyMap, faceAxialIndices, dir, res, adjacentVertexMap, triangleEdgeIntersections, vertexPos, offset == 0);
		}
	}


	for (int i = 0; i < 3; i++){
		if (triangleEdgeIntersections[i].intersectionTime[0] != -1){
			if (triangleEdgeIntersections[i].intersectionTime[1] == -1){
				printf("Triangle edge only have a single intersection! \n");
				return 0;
			}
			else{
				if (triangleEdgeIntersections[i].intersectionTime[1] < triangleEdgeIntersections[i].intersectionTime[0]){
					std::swap(triangleEdgeIntersections[i].intersectionTime[0], triangleEdgeIntersections[i].intersectionTime[1]);
					std::swap(triangleEdgeIntersections[i].intersectionKey[0], triangleEdgeIntersections[i].intersectionKey[1]);
				}
				adjacentVertexMap[triangleEdgeIntersections[i].intersectionKey[0]] = triangleEdgeIntersections[i].intersectionKey[1];
			}
		}
	}
}

// Vertex type indexing
// - - 2 - - 
// |		|
// 3   4	1
// |		|
// - - 0 - -


int SetEdgeProperties(ChainEdge & edge, unsigned long cIndex){

	unsigned long vElementKeys[2];

	unsigned long s_tElement, s_intersectionType;
	GetIntersectionKey(edge.sourceVertex, s_tElement, vElementKeys[0], s_intersectionType);
	unsigned long s_axialIndices[3];
	unsigned long s_cIndex;
	GetVoxelElementKeyIndices(vElementKeys[0], s_axialIndices, s_cIndex);

	unsigned long t_tElement, t_intersectionType;
	GetIntersectionKey(edge.targetVertex, t_tElement, vElementKeys[1], t_intersectionType);
	unsigned long t_axialIndices[3];
	unsigned long t_cIndex;
	GetVoxelElementKeyIndices(vElementKeys[1], t_axialIndices, t_cIndex);

	unsigned long faceAxialIndices[3] = { std::min<unsigned long>(s_axialIndices[0], t_axialIndices[0]), std::min<unsigned long>(s_axialIndices[1], t_axialIndices[1]), std::min<unsigned long>(s_axialIndices[2], t_axialIndices[2]) };

	unsigned long keys[5];
	{//j == 4
		keys[4] = SetVoxelElementKey(faceAxialIndices, cIndex);
	}
	{// j == 0
		keys[0] = SetVoxelElementKey(faceAxialIndices, (cIndex + 1) % 3);
	}
	{// j == 3
		keys[3] = SetVoxelElementKey(faceAxialIndices, (cIndex + 2) % 3);
	}
	{// j == 1
		faceAxialIndices[(cIndex + 1) % 3]++;
		keys[1] = SetVoxelElementKey(faceAxialIndices, (cIndex + 2) % 3);
		faceAxialIndices[(cIndex + 1) % 3]--;
	}
	{// j == 2
		faceAxialIndices[(cIndex + 2) % 3]++;
		keys[2] = SetVoxelElementKey(faceAxialIndices, (cIndex + 1) % 3);
		faceAxialIndices[(cIndex + 2) % 3]--;
	}


	unsigned char vertexType[2];
	for (int i = 0; i < 2; i++){
		bool foundKey = false;
		for (unsigned char j = 0; j < 5; j++){
			if (vElementKeys[i] == keys[j]){
				vertexType[i] = j;
				foundKey = true;
			}
		}
		if (!foundKey){
			printf("Undefined vertex type! \n");
			return 0;
		}
	}

	edge.faceKey = keys[4];
	edge.sourceVertexType = vertexType[0];
	edge.targetVertexType = vertexType[1];
	return 1;
}

int AddFace(std::unordered_map<unsigned long long, unsigned long long> & forwardMap, const FaceMap & faceMap, bool positiveOrientedFace){
	const std::unordered_map<unsigned long long, unsigned long long> & forwardVertexMap = faceMap.forwardVertex;
	for (auto forwardVertexIter = forwardVertexMap.begin(); forwardVertexIter != forwardVertexMap.end(); forwardVertexIter++){
		unsigned long long sourceVertex = (*forwardVertexIter).first;
		unsigned long long targetVertex = (*forwardVertexIter).second;
		if (positiveOrientedFace){
			if (forwardMap.find(sourceVertex) != forwardMap.end()){
				printf("Vertex already added! \n");
				return 0;
			}
			else{
				forwardMap[sourceVertex] = targetVertex;
			}

		}
		else{
			if (forwardMap.find(targetVertex) != forwardMap.end()){
				printf("Vertex already added! \n");
				return 0;
			}
			else{
				forwardMap[targetVertex] = sourceVertex;
			}
		}
	}
	return 1;
}

int RemoveFace(std::unordered_map<unsigned long long, unsigned long long> & forwardMap, const FaceMap & faceMap, bool positiveOrientedFace){
	const std::unordered_map<unsigned long long, unsigned long long> & forwardVertexMap = faceMap.forwardVertex;
	for (auto forwardVertexIter = forwardVertexMap.begin(); forwardVertexIter != forwardVertexMap.end(); forwardVertexIter++){
		unsigned long long sourceVertex = (*forwardVertexIter).first;
		unsigned long long targetVertex = (*forwardVertexIter).second;
		if (!positiveOrientedFace){
			if (forwardMap.find(sourceVertex) == forwardMap.end()){
				printf("Vertex to delete not found! \n");
				return 0;
			}
			else{
				forwardMap.erase(sourceVertex);
			}
		}
		else {
			if (forwardMap.find(targetVertex) == forwardMap.end()){
				printf("Vertex to delete not found! \n");
				return 0;
			}
			else{
				forwardMap.erase(targetVertex);
			}
		}
	}
	return 1;
}


class ComponentMap{
public:
	std::unordered_map<unsigned long long, unsigned long long> forwardMap;
	std::unordered_map<unsigned long long, unsigned long long> backwardMap;
	std::unordered_map<unsigned long, unsigned char> boundaryEdgesDegree;
	std::unordered_map<unsigned long, bool> boundaryFaceOrientation;
	void clear(){
		forwardMap.clear();
		backwardMap.clear();
		boundaryEdgesDegree.clear();
		boundaryFaceOrientation.clear();
	}
};

int AddFace(ComponentMap & refMap, const FaceMap & faceMap, bool positiveOrientedFace){
	const std::unordered_map<unsigned long long, unsigned long long> & forwardVertexMap = faceMap.forwardVertex;
	for (auto forwardVertexIter = forwardVertexMap.begin(); forwardVertexIter != forwardVertexMap.end(); forwardVertexIter++){
		unsigned long long sourceVertex = (*forwardVertexIter).first;
		unsigned long long targetVertex = (*forwardVertexIter).second;
		if (positiveOrientedFace){
			if (refMap.forwardMap.find(sourceVertex) != refMap.forwardMap.end()){
				printf("source vertex already added in forward map! %llu \n", sourceVertex);
				return 0;
			}
			else{
				refMap.forwardMap[sourceVertex] = targetVertex;
			}

		}
		else{
			if (refMap.forwardMap.find(targetVertex) != refMap.forwardMap.end()){
				printf("target vertex already added in forward map! %llu \n", targetVertex);
				return 0;
			}
			else{
				refMap.forwardMap[targetVertex] = sourceVertex;
			}
		}
	}

	const std::unordered_map<unsigned long long, unsigned long long> & backwardVertexMap = faceMap.backwardVertex;
	for (auto backwardVertexIter = backwardVertexMap.begin(); backwardVertexIter != backwardVertexMap.end(); backwardVertexIter++){
		unsigned long long sourceVertex = (*backwardVertexIter).first;
		unsigned long long targetVertex = (*backwardVertexIter).second;
		if (positiveOrientedFace){
			if (refMap.backwardMap.find(sourceVertex) != refMap.backwardMap.end()){
				printf("source vertex already added in backward map! %llu \n", sourceVertex);
				return 0;
			}
			else{
				refMap.backwardMap[sourceVertex] = targetVertex;
			}

		}
		else{
			if (refMap.backwardMap.find(targetVertex) != refMap.backwardMap.end()){
				printf("target vertex already added in backward map! %llu \n", targetVertex);
				return 0;
			}
			else{
				refMap.backwardMap[targetVertex] = sourceVertex;
			}
		}
	}
	return 1;
}


int RemoveFace(ComponentMap & refMap, const FaceMap & faceMap, bool positiveOrientedFace){
	const std::unordered_map<unsigned long long, unsigned long long> & forwardVertexMap = faceMap.forwardVertex;
	for (auto forwardVertexIter = forwardVertexMap.begin(); forwardVertexIter != forwardVertexMap.end(); forwardVertexIter++){
		unsigned long long sourceVertex = (*forwardVertexIter).first;
		unsigned long long targetVertex = (*forwardVertexIter).second;
		if (positiveOrientedFace){
			if (refMap.forwardMap.find(sourceVertex) == refMap.forwardMap.end()){
				printf("Vertex to delete not found! \n");
				return 0;
			}
			else{
				refMap.forwardMap.erase(sourceVertex);
			}
		}
		else{
			if (refMap.forwardMap.find(targetVertex) == refMap.forwardMap.end()){
				printf("Vertex to delete not found! \n");
				return 0;
			}
			else{
				refMap.forwardMap.erase(targetVertex);
			}
		}
	}

	const std::unordered_map<unsigned long long, unsigned long long> & backwardVertexMap = faceMap.backwardVertex;
	for (auto backwardVertexIter = backwardVertexMap.begin(); backwardVertexIter != backwardVertexMap.end(); backwardVertexIter++){
		unsigned long long sourceVertex = (*backwardVertexIter).first;
		unsigned long long targetVertex = (*backwardVertexIter).second;
		if (positiveOrientedFace){
			if (refMap.backwardMap.find(sourceVertex) == refMap.backwardMap.end()){
				printf("Vertex to delete not found! \n");
				return 0;
			}
			else{
				refMap.backwardMap.erase(sourceVertex);
			}
		}
		else {
			if (refMap.backwardMap.find(targetVertex) == refMap.backwardMap.end()){
				printf("Vertex to delete not found! \n");
				return 0;
			}
			else{
				refMap.backwardMap.erase(targetVertex);
			}
		}
	}
	return 1;
}

void GetNeigbourEdgeKeys(unsigned long faceKey, unsigned long edgeKeys[4]){
	unsigned long axialIndices[3];
	unsigned long cIndex;
	GetVoxelElementKeyIndices(faceKey, axialIndices, cIndex);
	for (unsigned long offset = 0; offset < 2; offset++){
		for (unsigned long dir = 0; dir < 2; dir++){
			unsigned long _axialIndices[3] = { axialIndices[0], axialIndices[1], axialIndices[2] };
			if (offset == 1) _axialIndices[(cIndex + 1 + dir) % 3]++;
			edgeKeys[offset * 2 + dir] = SetVoxelElementKey(_axialIndices, (cIndex + 2 - dir) % 3);
		}
	}
}

void GetNeigbourFaceKeys(unsigned long edgeKey, unsigned long faceKeys[4]){
	unsigned long axialIndices[3];
	unsigned long cIndex;
	GetVoxelElementKeyIndices(edgeKey, axialIndices, cIndex);
	for (unsigned long offset = 0; offset < 2; offset++){
		for (unsigned long dir = 0; dir < 2; dir++){
			unsigned long _axialIndices[3] = { axialIndices[0], axialIndices[1], axialIndices[2] };
			if (offset == 1) _axialIndices[(cIndex + 1 + dir) % 3]--;
			faceKeys[offset * 2 + dir] = SetVoxelElementKey(_axialIndices, (cIndex + 2 - dir) % 3);
		}
	}
}

int ConstructVoxelMap(ComponentMap & voxelMap, std::unordered_map<unsigned long, FaceMap> & faceMaps, unsigned long axialIndices[3]){
	voxelMap.clear();
	for (unsigned long c = 0; c < 3; c++){
		for (unsigned long offset = 0; offset < 2; offset++){
			unsigned long faceAxialIndices[3] = { axialIndices[0], axialIndices[1], axialIndices[2] };
			if (offset) faceAxialIndices[c]++;
			unsigned long faceKey = SetVoxelElementKey(faceAxialIndices, c);
			if(!AddFace(voxelMap, faceMaps[faceKey], offset)) return 0;
			voxelMap.boundaryFaceOrientation[faceKey] = offset;
			unsigned long edgeKeys[4];
			GetNeigbourEdgeKeys(faceKey, edgeKeys);
			for (int e = 0; e < 4; e++) voxelMap.boundaryEdgesDegree[edgeKeys[e]] = 1;
		}
	}
	return 1;
}

int ReplicateEdgeIntersections(std::unordered_map<unsigned long, FaceMap> & faceMaps, EdgeIntersections & edge, ComponentMap & refMap, ComponentMap & addMap, unsigned long faceKeys[4], int changeFaceIndex, std::unordered_map<unsigned long long, SamplePoint> & vertexSample, bool verbose){
	for (auto vertexIter = edge.vertices.begin(); vertexIter != edge.vertices.end(); vertexIter++){
		unsigned long long vertexKey = (*vertexIter);
		unsigned long vElement = static_cast<unsigned long>(vertexKey & 0x00000000FFFFFFFF);
		unsigned long tElement = static_cast<unsigned long>((vertexKey >> 32) & 0x000000007FFFFFFF);
		unsigned long intersectionType = static_cast<unsigned long>((vertexKey >> 63) & 0x0000000000000001);
		unsigned long triangleOffset = 2 << 29;
		unsigned long long copyVertexKey = SetIntersectionKey(triangleOffset + tElement, vElement, intersectionType);

		bool faceOrientation[2] = { refMap.boundaryFaceOrientation[faceKeys[changeFaceIndex]], addMap.boundaryFaceOrientation[faceKeys[(changeFaceIndex + 1) % 4]] };
		if(verbose)printf("Face Orientations %d %d \n", faceOrientation[0], faceOrientation[1]);
		if(verbose)printf("Original Vertex %llu \n", vertexKey);
		if(verbose)printf("Vertex Copy %llu \n", copyVertexKey);
		bool isInFowardMap[2];
		for (int i = 0; i < 2; i++){
			FaceMap & faceMap = faceMaps[faceKeys[(changeFaceIndex + i) % 4]];
			
			unsigned char vertexType = faceMap.vertexType[vertexKey];
			faceMap.vertexType.erase(vertexKey);
			faceMap.vertexType[copyVertexKey] = vertexType;

			if (faceMap.forwardVertex.find(vertexKey) != faceMap.forwardVertex.end()){
				unsigned long long forward = faceMap.forwardVertex[vertexKey];
				faceMap.forwardVertex[copyVertexKey] = forward;
				faceMap.backwardVertex[forward] = copyVertexKey;
				faceMap.forwardVertex.erase(vertexKey);
				if (i == 0){
					if (faceOrientation[i]){
						refMap.forwardMap[copyVertexKey] = forward;
						refMap.forwardMap.erase(vertexKey);
						if(verbose) printf("Deleting vertex from forward map %llu \n", vertexKey);
						refMap.backwardMap[forward] = copyVertexKey;
					}
					else{
						refMap.backwardMap[copyVertexKey] = forward;
						refMap.backwardMap.erase(vertexKey);
						if(verbose) printf("Deleting vertex from backward map %llu \n", vertexKey);
						refMap.forwardMap[forward] = copyVertexKey;
					}
				}
				else{
					addMap.forwardMap[copyVertexKey] = forward;
					addMap.forwardMap.erase(vertexKey);
					addMap.backwardMap[forward] = copyVertexKey;
					//addMap.backwardMap[copyVertexKey] = vertexKey;
				}
				isInFowardMap[i] = faceOrientation[i];
			}
			else if (faceMap.backwardVertex.find(vertexKey) != faceMap.backwardVertex.end()){
				unsigned long long backward = faceMap.backwardVertex[vertexKey];
				faceMap.backwardVertex[copyVertexKey] = backward;
				faceMap.forwardVertex[backward] = copyVertexKey;
				faceMap.backwardVertex.erase(vertexKey);
				if (i == 0){
					if (faceOrientation[i]){
						refMap.backwardMap[copyVertexKey] = backward;
						refMap.backwardMap.erase(vertexKey);
						if(verbose) printf("Deleting vertex from backward map %llu \n", vertexKey);
						refMap.forwardMap[backward] = copyVertexKey;
						//refMap.forwardMap[copyVertexKey] = vertexKey;
					}
					else{
						refMap.forwardMap[copyVertexKey] = backward;
						refMap.forwardMap.erase(vertexKey);
						if(verbose) printf("Deleting vertex from forward map %llu \n", vertexKey);
						refMap.backwardMap[backward] = copyVertexKey;
					}
				}
				else{
					addMap.backwardMap[copyVertexKey] = backward;
					addMap.backwardMap.erase(vertexKey);
					addMap.forwardMap[backward] = copyVertexKey;
					//addMap.forwardMap[copyVertexKey] = vertexKey;
				}
				isInFowardMap[i] = !faceOrientation[i];
			}
			else{
				printf("Not found vertex map!\n");
				return 0;
			}
		}
		if (isInFowardMap[0] == isInFowardMap[1]){
			printf("Adjacent faces map is invalid !\n");
			return 0;
		}
		vertexSample[copyVertexKey] = vertexSample[vertexKey];
	}
}

int SingularEdgeIntersection(ComponentMap & refMap, ComponentMap & addMap, std::unordered_map<unsigned long, FaceMap> & faceMaps, std::unordered_map<unsigned long, EdgeIntersections> & sourceEdgeIntersections, std::unordered_map<unsigned long, EdgeIntersections> & targetEdgeIntersections){
	for (auto addBoundaryEdgesIter = addMap.boundaryEdgesDegree.begin(); addBoundaryEdgesIter != addMap.boundaryEdgesDegree.end(); addBoundaryEdgesIter++){
		unsigned long edgeKey = (*addBoundaryEdgesIter).first;
		unsigned long addDegree = (*addBoundaryEdgesIter).second;
		if (refMap.boundaryEdgesDegree.find(edgeKey) != refMap.boundaryEdgesDegree.end()){
			unsigned long refDegree = refMap.boundaryEdgesDegree[edgeKey];
			if (addDegree == 1 && refDegree == 1){
				unsigned long faceKeys[4];
				bool isInAddBoundary[4];
				bool isInRefBoundary[4];
				int isInAddBoundaryCount = 0;
				int isInRefBoundaryCount = 0;
				GetNeigbourFaceKeys(edgeKey, faceKeys);
				bool criticalEdge = true;
				for (int j = 0; j < 4; j++){
					unsigned long currentFace = faceKeys[j];
					if (addMap.boundaryFaceOrientation.find(currentFace) == addMap.boundaryFaceOrientation.end()){
						isInAddBoundary[j] = false;
					}
					else{
						isInAddBoundary[j] = true;
						isInAddBoundaryCount++;
					}
					if (refMap.boundaryFaceOrientation.find(currentFace) == refMap.boundaryFaceOrientation.end()){
						isInRefBoundary[j] = false;
					}
					else{
						isInRefBoundary[j] = true;
						isInRefBoundaryCount++;
					}
					if (isInAddBoundary[j] == isInRefBoundary[j]){
						//printf("common face!");
						criticalEdge = false;
						//break;
					}
				}

				if (criticalEdge && (sourceEdgeIntersections[edgeKey].vertices.size() > 0 || targetEdgeIntersections[edgeKey].vertices.size() > 0)){
					if (isInRefBoundaryCount != 2 || isInAddBoundaryCount != 2){
						printf("Unexpected face count! \n");
						return 0;
					}
					if (isInRefBoundary[0] == isInRefBoundary[2] || isInAddBoundary[0] == isInAddBoundary[2]){
						printf("Unexpected face adjacency! \n");
						return 0;
					}
					int changeFaceIndex = -1;
					for (int i = 0; i < 4; i++){
						if (isInRefBoundary[i] && isInAddBoundary[(i + 1) % 4]){
							changeFaceIndex = i;
						}
					}

					if (changeFaceIndex == -1){
						printf("Not found delete face index! \n");
						return 0;
					}
					return 2;
				}
			}
		}
	}
	return 1;
}

int MergeComponents(ComponentMap & refMap, ComponentMap & addMap, std::unordered_map<unsigned long, FaceMap> & faceMaps, std::unordered_map<unsigned long, EdgeIntersections> & edgeIntersections, std::unordered_map<unsigned long long, SamplePoint> & vertexSample, bool verbose){
	
	// Process Edges
	
	std::unordered_set<unsigned long> edgesToAdd;
	std::unordered_set<unsigned long> edgesToEdit;
	for (auto addBoundaryEdgesIter = addMap.boundaryEdgesDegree.begin(); addBoundaryEdgesIter != addMap.boundaryEdgesDegree.end(); addBoundaryEdgesIter++){
		unsigned long edgeKey = (*addBoundaryEdgesIter).first;
		unsigned long addDegree = (*addBoundaryEdgesIter).second;
		if (refMap.boundaryEdgesDegree.find(edgeKey) != refMap.boundaryEdgesDegree.end()){
			unsigned long refDegree = refMap.boundaryEdgesDegree[edgeKey];
			if (addDegree == 1 && refDegree == 1){
				unsigned long faceKeys[4];
				bool isInAddBoundary[4];
				bool isInRefBoundary[4];
				int isInAddBoundaryCount = 0;
				int isInRefBoundaryCount = 0;
				GetNeigbourFaceKeys(edgeKey, faceKeys);
				bool criticalEdge = true;
				for (int j = 0; j < 4; j++){
					unsigned long currentFace = faceKeys[j];
					if (addMap.boundaryFaceOrientation.find(currentFace) == addMap.boundaryFaceOrientation.end()){
						isInAddBoundary[j] = false;
					}
					else{
						isInAddBoundary[j] = true;
						isInAddBoundaryCount++;
					}
					if (refMap.boundaryFaceOrientation.find(currentFace) == refMap.boundaryFaceOrientation.end()){
						isInRefBoundary[j] = false;
					}
					else{
						isInRefBoundary[j] = true;
						isInRefBoundaryCount++;
					}
					if (isInAddBoundary[j] == isInRefBoundary[j]){
						//printf("common face!");
						criticalEdge = false;
						//break;
					}
				}

				if (criticalEdge && edgeIntersections[edgeKey].vertices.size()){
					if (isInRefBoundaryCount != 2 || isInAddBoundaryCount != 2){
						printf("Unexpected face count! \n");
						return 0;
					}
					if (isInRefBoundary[0] == isInRefBoundary[2] || isInAddBoundary[0] == isInAddBoundary[2]){
						printf("Unexpected face adjacency! \n");
						return 0;
					}
					int changeFaceIndex = -1;
					for (int i = 0; i < 4; i++){
						if (isInRefBoundary[i] && isInAddBoundary[(i + 1) % 4]){
							changeFaceIndex = i;
						}
					}

					if (changeFaceIndex == -1){
						printf("Not found delete face index! \n");
						return 0;
					}
					
					if (verbose){
						printf("Replicating edge intersections! \n");
						printf("Edge %llu \n", edgeKey);
						printf("Add %d %d %d %d\n", isInAddBoundary[0], isInAddBoundary[1], isInAddBoundary[2], isInAddBoundary[3]);
						printf("Ref %d %d %d %d \n", isInRefBoundary[0], isInRefBoundary[1], isInRefBoundary[2], isInRefBoundary[3]);
						printf("Face 0 %llu \n", faceKeys[0]);
						printf("Face 1 %llu \n", faceKeys[1]);
						printf("Face 2 %llu \n", faceKeys[2]);
						printf("Face 3 %llu \n", faceKeys[3]);
						printf("changeFaceIndex %d \n", changeFaceIndex);
					}
					if (!ReplicateEdgeIntersections(faceMaps, edgeIntersections[edgeKey], refMap, addMap, faceKeys, changeFaceIndex, vertexSample, verbose)) return 0;
				}
			}
			edgesToEdit.insert(edgeKey);
		}
		else{
			edgesToAdd.insert(edgeKey);
		}
	}
	for (auto edgeIter = edgesToEdit.begin(); edgeIter != edgesToEdit.end(); edgeIter++){
		unsigned long edgeKey = (*edgeIter);
		unsigned long refDegree = refMap.boundaryEdgesDegree[edgeKey];
		unsigned long addDegree = addMap.boundaryEdgesDegree[edgeKey];
		if ((addDegree + refDegree) == 4){
			refMap.boundaryEdgesDegree.erase(edgeKey);
		}
		else{
			refMap.boundaryEdgesDegree[edgeKey] = addDegree + refDegree;
		}
		
	}
	for (auto edgeIter = edgesToAdd.begin(); edgeIter != edgesToAdd.end(); edgeIter++){
		refMap.boundaryEdgesDegree[*edgeIter] = addMap.boundaryEdgesDegree[*edgeIter];
	}

	//Process Faces
	std::unordered_set<unsigned long> facesToAdd;
	for (auto addFaceIter = addMap.boundaryFaceOrientation.begin(); addFaceIter != addMap.boundaryFaceOrientation.end(); addFaceIter++){
		unsigned long faceKey = (*addFaceIter).first;
		if (refMap.boundaryFaceOrientation.find(faceKey) != refMap.boundaryFaceOrientation.end()){
			if(!RemoveFace(refMap, faceMaps[faceKey], refMap.boundaryFaceOrientation[faceKey]))return 0;
			refMap.boundaryFaceOrientation.erase(faceKey);
		}
		else{
			facesToAdd.insert(faceKey);
		}
	}

	for (auto faceIter = facesToAdd.begin(); faceIter != facesToAdd.end(); faceIter++){
		if (verbose) {
			printf("Adding face %llu  with orientation %d \n", *faceIter, addMap.boundaryFaceOrientation[*faceIter]);
		}
		if(!AddFace(refMap, faceMaps[*faceIter], addMap.boundaryFaceOrientation[*faceIter]))return 0;
		refMap.boundaryFaceOrientation[*faceIter] = addMap.boundaryFaceOrientation[*faceIter];
	}

	return 1;
}

int AddVoxel(ComponentMap & refMap, std::unordered_map<unsigned long, FaceMap> & faceMaps, std::unordered_map<unsigned long, EdgeIntersections> & edgeIntersections, std::unordered_map<unsigned long long, SamplePoint> & vertexSample, unsigned long axialIndices[3], bool verbose, bool allowReplication){
	ComponentMap voxelMap;
	ConstructVoxelMap(voxelMap, faceMaps, axialIndices);
	if (allowReplication){
		if (!MergeComponents(refMap, voxelMap, faceMaps, edgeIntersections, vertexSample, verbose)) return 0;
		return 1;
	}
	else{
		int testReplication = SingularEdgeIntersection(refMap, voxelMap, faceMaps, edgeIntersections, edgeIntersections);
		if (testReplication == 1){
			if (!MergeComponents(refMap, voxelMap, faceMaps, edgeIntersections, vertexSample, verbose)) return 0;
			return 1;
		}
		else return testReplication;
	}
}

int AddVoxel(ComponentMap & sourceRefMap, std::unordered_map<unsigned long, FaceMap> & sourceFaceMaps, std::unordered_map<unsigned long, EdgeIntersections> & sourceEdgeIntersections, std::unordered_map<unsigned long long, SamplePoint> & sourceVertexSample,
	ComponentMap & targetRefMap, std::unordered_map<unsigned long, FaceMap> & targetFaceMaps, std::unordered_map<unsigned long, EdgeIntersections> & targetEdgeIntersections, std::unordered_map<unsigned long long, SamplePoint> & targetVertexSample,
			 unsigned long axialIndices[3], bool verbose, bool allowReplication){
	
	ComponentMap sourceVoxelMap;
	ConstructVoxelMap(sourceVoxelMap, sourceFaceMaps, axialIndices);
	ComponentMap targetVoxelMap;
	ConstructVoxelMap(targetVoxelMap, targetFaceMaps, axialIndices);
	if (allowReplication){
		if (!MergeComponents(sourceRefMap, sourceVoxelMap, sourceFaceMaps, sourceEdgeIntersections, sourceVertexSample, verbose)) return 0;
		if (!MergeComponents(targetRefMap, targetVoxelMap, targetFaceMaps, targetEdgeIntersections, targetVertexSample, verbose)) return 0;
		return 1;
	}
	else{
		int testSourceReplication = SingularEdgeIntersection(sourceRefMap, sourceVoxelMap, sourceFaceMaps, sourceEdgeIntersections, targetEdgeIntersections);
		if (testSourceReplication == 1){
			if (!MergeComponents(sourceRefMap, sourceVoxelMap, sourceFaceMaps, sourceEdgeIntersections, sourceVertexSample, verbose)) return 0;
			if (!MergeComponents(targetRefMap, targetVoxelMap, targetFaceMaps, targetEdgeIntersections, targetVertexSample, verbose)) return 0;
			return 1;
		}
		else return testSourceReplication;
	}
}


int LoopRepresentatives(std::unordered_map<unsigned long long, unsigned long long> & forwardMap, std::vector<unsigned long long> & loopRepresentatives){
	loopRepresentatives.clear();
	std::unordered_set<unsigned long long> alreadyVisitedVertex;
	unsigned int loopCounter = 0;
	unsigned int edgeCounter = 0;
	for (auto iter = forwardMap.begin(); iter != forwardMap.end(); iter++){
		unsigned long long sourceVertex = (*iter).first;
		if (alreadyVisitedVertex.find(sourceVertex) == alreadyVisitedVertex.end()){
			loopRepresentatives.push_back(sourceVertex);
			unsigned long long currentVertex = sourceVertex;
			bool terminate = false;
			unsigned int startEdgeCounter = edgeCounter;
			do{
				if (alreadyVisitedVertex.find(currentVertex) == alreadyVisitedVertex.end()){
					alreadyVisitedVertex.insert(currentVertex);
					auto mappedVertex = forwardMap.find(currentVertex);
					if (mappedVertex != forwardMap.end()){
						unsigned long long nextVertex = (*mappedVertex).second;
						edgeCounter++;
						currentVertex = nextVertex;
					}
					else{
						printf("Vertex to dead end! \n");
						return -1;
					}
				}
				else{
					if (currentVertex != sourceVertex){
						printf("Non simple loop! \n");
						return -1;
					}
					terminate = true;
				}
			} while (!terminate);
			loopCounter++;
		}
	}
}


int LoopVertices(std::unordered_map<unsigned long long, unsigned long long> & forwardMap, std::vector<std::vector<unsigned long long>> & loopVertices){
	loopVertices.clear();
	std::unordered_set<unsigned long long> alreadyVisitedVertex;
	unsigned int loopCounter = 0;
	unsigned int edgeCounter = 0;
	for (auto iter = forwardMap.begin(); iter != forwardMap.end(); iter++){
		unsigned long long sourceVertex = (*iter).first;
		if (alreadyVisitedVertex.find(sourceVertex) == alreadyVisitedVertex.end()){
			std::vector< unsigned long long> currentLoop;
			unsigned long long currentVertex = sourceVertex;
			bool terminate = false;
			unsigned int startEdgeCounter = edgeCounter;
			do{
				if (alreadyVisitedVertex.find(currentVertex) == alreadyVisitedVertex.end()){
					alreadyVisitedVertex.insert(currentVertex);
					currentLoop.push_back(currentVertex);
					auto mappedVertex = forwardMap.find(currentVertex);
					if (mappedVertex != forwardMap.end()){
						unsigned long long nextVertex = (*mappedVertex).second;
						edgeCounter++;
						currentVertex = nextVertex;
					}
					else{
						printf("Vertex to dead end! \n");
						return -1;
					}
				}
				else{
					if (currentVertex != sourceVertex){
						printf("Non simple loop! \n");
						return -1;
					}
					terminate = true;
				}
			} while (!terminate);
			loopCounter++;
			loopVertices.push_back(currentLoop);
		}
	}
}



bool SimilarFaces(FaceMap & sourceFace, FaceMap & targetFace){
	if (sourceFace.forwardVertex.size() == 0 && targetFace.forwardVertex.size() == 0){
		//printf("Similar Face!");
		return true;
	}
	else if (sourceFace.forwardVertex.size() == 1 && targetFace.forwardVertex.size() == 1){

		unsigned long long sourceVertex = (*sourceFace.forwardVertex.begin()).first;
		unsigned long long targetVertex = (*targetFace.forwardVertex.begin()).first;
		unsigned long s_tElement, s_vElement, s_intersectionType;
		unsigned long t_tElement, t_vElement, t_intersectionType;
		GetIntersectionKey(sourceVertex, s_tElement, s_vElement, s_intersectionType);
		GetIntersectionKey(targetVertex, t_tElement, t_vElement, t_intersectionType);
		if (s_vElement == t_vElement){
			sourceVertex = (*sourceFace.forwardVertex.begin()).second;
			targetVertex = (*targetFace.forwardVertex.begin()).second;
			GetIntersectionKey(sourceVertex, s_tElement, s_vElement, s_intersectionType);
			GetIntersectionKey(targetVertex, t_tElement, t_vElement, t_intersectionType);
			if (s_vElement == t_vElement){
				return true;
			}
			else return false;
		}
		else return false;
	}
	else{
		return false;
	}
}

void FindSimilarFaces(std::unordered_map<unsigned long, FaceMap> & sourceFaceMaps, std::unordered_map<unsigned long, FaceMap> & targetFaceMaps, std::unordered_set<unsigned long> & similarFaces){
	for (auto faceIter = sourceFaceMaps.begin(); faceIter != sourceFaceMaps.end(); faceIter++){
		unsigned long faceKey = (*faceIter).first;
		if (targetFaceMaps.find(faceKey) != targetFaceMaps.end()){
			FaceMap & sourceFace = (*faceIter).second;
			FaceMap & targetFace = targetFaceMaps[faceKey];
			if (SimilarFaces(sourceFace, targetFace)){
				similarFaces.insert(faceKey);
			}
		}
	}
}

void FindWeakFaces(std::unordered_map<unsigned long, FaceMap> & sourceFaceMaps, std::unordered_map<unsigned long long, double > & sourceVertexConfidence,
	std::unordered_map<unsigned long, FaceMap> & targetFaceMaps, std::unordered_map<unsigned long long, double > & targetVertexConfidence,
	std::unordered_set<unsigned long> & weakFaces, const double threshold){

	for (auto faceIter = sourceFaceMaps.begin(); faceIter != sourceFaceMaps.end(); faceIter++){
		unsigned long faceKey = (*faceIter).first;
		FaceMap & face = (*faceIter).second;
		for (auto vertexIter = face.vertexType.begin(); vertexIter != face.vertexType.end(); vertexIter++){
			unsigned long long vertexKey = (*vertexIter).first;
			if (sourceVertexConfidence[vertexKey] < threshold) weakFaces.insert(faceKey);
		}
	}

	for (auto faceIter = targetFaceMaps.begin(); faceIter != targetFaceMaps.end(); faceIter++){
		unsigned long faceKey = (*faceIter).first;
		FaceMap & face = (*faceIter).second;
		for (auto vertexIter = face.vertexType.begin(); vertexIter != face.vertexType.end(); vertexIter++){
			unsigned long long vertexKey = (*vertexIter).first;
			if (targetVertexConfidence[vertexKey] < threshold) weakFaces.insert(faceKey);
		}
	}
}




int GrowVoxelRegion(const int gridRes, std::unordered_map<unsigned long, FaceMap> & faceMaps, std::unordered_map<unsigned long, EdgeIntersections> & edgeIntersections, std::unordered_map<unsigned long long, SamplePoint> & vertexSample, std::unordered_set<unsigned long> & voxelKeys, const std::vector<unsigned long> & voxelSeeds, ComponentMap & refMap, std::vector<unsigned long> & addedVoxels, unsigned int growingRadius, bool allowReplication, bool stopAtNonWeakFaces, std::unordered_set<unsigned long> & weakFaces){
	//Grow voxel region
	std::unordered_map<unsigned long, unsigned int> addedVoxelsDepth;
	std::unordered_set<unsigned long> addedFaces;
	std::queue<unsigned long> voxelsQueue;
	std::unordered_map<unsigned long, unsigned int> skippingCount;
	int iter = 0;

	for (int i = 0; i < voxelSeeds.size(); i++){
		voxelsQueue.push(voxelSeeds[i]);
		addedVoxelsDepth[voxelSeeds[i]] = 0;
	}
	do{
		iter++;
		unsigned long currentVoxel = voxelsQueue.front();
		voxelsQueue.pop();
		unsigned long axialIndices[3];
		unsigned long _c;
		GetVoxelElementKeyIndices(currentVoxel, axialIndices, _c);
		if(0) printf("Adding voxel(%d,%d,%d) \n", axialIndices[0], axialIndices[1], axialIndices[2]);
		int addingOutCome = AddVoxel(refMap, faceMaps, edgeIntersections, vertexSample, axialIndices, false, allowReplication);
		if (addingOutCome == 0) return 0;
		else if (addingOutCome == 1){
			addedVoxels.push_back(currentVoxel);
			unsigned int currentDepth = addedVoxelsDepth[currentVoxel];
			if (currentDepth < growingRadius){
				for (unsigned long c = 0; c < 3; c++){
					for (unsigned long offset = 0; offset < 2; offset++){
						unsigned long neigbourVoxelIndices[3] = { axialIndices[0], axialIndices[1], axialIndices[2] };
						unsigned long neigbourFaceIndices[3] = { axialIndices[0], axialIndices[1], axialIndices[2] };
						if (neigbourVoxelIndices[c] > 0 && offset == 0){
							neigbourVoxelIndices[c]--;
						}
						else if (neigbourVoxelIndices[c] < gridRes - 1 && offset == 1){
							neigbourVoxelIndices[c]++;
							neigbourFaceIndices[c]++;
						}
						unsigned long neighbourVoxelKey = SetVoxelElementKey(neigbourVoxelIndices, 0);
						unsigned long neighbourFaceKey = SetVoxelElementKey(neigbourFaceIndices, c);
						if (voxelKeys.find(neighbourVoxelKey) != voxelKeys.end()){
							if (addedVoxelsDepth.find(neighbourVoxelKey) == addedVoxelsDepth.end()){
								if (!stopAtNonWeakFaces || weakFaces.find(neighbourFaceKey) != weakFaces.end()){
									voxelsQueue.push(neighbourVoxelKey);
									addedVoxelsDepth[neighbourVoxelKey] = currentDepth + 1;
								}
							}
						}
					}
				}
			}
		}
		else if(addingOutCome == 2){
			printf("Skipped voxel %llu! \n", currentVoxel);
			if (skippingCount.find(currentVoxel) == skippingCount.end())skippingCount[currentVoxel] = 0;
			else skippingCount[currentVoxel]++;
			if (skippingCount[currentVoxel] == 5) printf("Voxel %llu being skipped multiple times! \n", currentVoxel);
			else voxelsQueue.push(currentVoxel);
		}
		else{
			printf("Unknown adding outcome! \n", currentVoxel);
		}
	} while (!voxelsQueue.empty());
	if (0)printf("Voxels added %d \n", addedVoxelsDepth.size());
	return 1;

}


int GrowVoxelRegion(const int gridRes, 
	std::unordered_map<unsigned long, FaceMap> & sourceFaceMaps, std::unordered_map<unsigned long, EdgeIntersections> & sourceEdgeIntersections, std::unordered_map<unsigned long long, SamplePoint> & sourceVertexSample, ComponentMap & sourceRefMap,
	std::unordered_map<unsigned long, FaceMap> & targetFaceMaps, std::unordered_map<unsigned long, EdgeIntersections> & targetEdgeIntersections, std::unordered_map<unsigned long long, SamplePoint> & targetVertexSample, ComponentMap & targetRefMap,
	std::unordered_set<unsigned long> & voxelKeys, const std::vector<unsigned long> & voxelSeeds, std::unordered_set<unsigned long> & addedVoxels, int maxGrowingRadius, bool allowReplication, bool stopAtNonWeakFaces, std::unordered_set<unsigned long> & weakFaces){
	
	//Grow voxel region
	std::unordered_map<unsigned long, unsigned int> addedVoxelsDepth;
	std::unordered_set<unsigned long> addedFaces;
	std::queue<unsigned long> voxelsQueue;
	int iter = 0;

	std::unordered_map<unsigned long, unsigned int> skippingCount;

	unsigned int currentDepth = 0;

	for (int i = 0; i < voxelSeeds.size(); i++){
		voxelsQueue.push(voxelSeeds[i]);
		addedVoxelsDepth[voxelSeeds[i]] = 0;
	}
	do{
		iter++;
		unsigned long currentVoxel = voxelsQueue.front();
		voxelsQueue.pop();
		unsigned long axialIndices[3];
		unsigned long _c;
		GetVoxelElementKeyIndices(currentVoxel, axialIndices, _c);
		if (0) printf("Adding voxel(%d,%d,%d) \n", axialIndices[0], axialIndices[1], axialIndices[2]);
		int addingOutCome = AddVoxel(sourceRefMap, sourceFaceMaps, sourceEdgeIntersections, sourceVertexSample, targetRefMap, targetFaceMaps, targetEdgeIntersections, targetVertexSample, axialIndices, false, allowReplication);
		if (addingOutCome == 0) return 0;
		else if (addingOutCome == 1){
			addedVoxels.insert(currentVoxel);
			currentDepth = addedVoxelsDepth[currentVoxel];
			if (currentDepth < maxGrowingRadius){
				for (unsigned long c = 0; c < 3; c++){
					for (unsigned long offset = 0; offset < 2; offset++){
						unsigned long neigbourVoxelIndices[3] = { axialIndices[0], axialIndices[1], axialIndices[2] };
						unsigned long neigbourFaceIndices[3] = { axialIndices[0], axialIndices[1], axialIndices[2] };
						if (neigbourVoxelIndices[c] > 0 && offset == 0){
							neigbourVoxelIndices[c]--;
						}
						else if (neigbourVoxelIndices[c] < gridRes - 1 && offset == 1){
							neigbourFaceIndices[c]++;
							neigbourVoxelIndices[c]++;
						}
						unsigned long neighbourVoxelKey = SetVoxelElementKey(neigbourVoxelIndices, 0);
						unsigned long neighbourFaceKey = SetVoxelElementKey(neigbourFaceIndices, c);
						if (voxelKeys.find(neighbourVoxelKey) != voxelKeys.end()){
							if (addedVoxelsDepth.find(neighbourVoxelKey) == addedVoxelsDepth.end()){
								if ( (!stopAtNonWeakFaces || weakFaces.find(neighbourFaceKey) != weakFaces.end()) && !SimilarFaces(sourceFaceMaps[neighbourFaceKey], targetFaceMaps[neighbourFaceKey])){
									voxelsQueue.push(neighbourVoxelKey);
									addedVoxelsDepth[neighbourVoxelKey] = currentDepth + 1;
								}
							}
						}
					}
				}
			}
		}
		else if (addingOutCome == 2){
			printf("Skipped voxel %llu! \n", currentVoxel);
			if (skippingCount.find(currentVoxel) == skippingCount.end())skippingCount[currentVoxel] = 0;
			else skippingCount[currentVoxel]++;
			if (skippingCount[currentVoxel] == 5) printf("Voxel %llu being skipped multiple times! \n", currentVoxel);
			else voxelsQueue.push(currentVoxel);
		}
		else{
			printf("Unknown adding outcome! \n", currentVoxel);
		}
	} while (!voxelsQueue.empty());
	if (1)printf("Voxels added %d \n", addedVoxelsDepth.size());
	if (1)printf("Current depth  %d. Max depth %d \n", currentDepth, maxGrowingRadius);
	return 1;

}


int AddComponent(std::unordered_map<unsigned long, int> & vertexComponent, unsigned long vIndex, int currentComponent, std::unordered_map<unsigned long, std::unordered_set<unsigned long>> & neighbours, std::unordered_set<unsigned long> & componentElements){
	vertexComponent[vIndex] = currentComponent;
	std::queue<unsigned long> visitingQueue;
	visitingQueue.push(vIndex);
	componentElements.insert(vIndex);

	while (!visitingQueue.empty()){
		unsigned long currentVertex = visitingQueue.front();
		visitingQueue.pop();
		const std::unordered_set<unsigned long> & vertexNeighbours = neighbours[currentVertex];
		for (auto neighbourIter = vertexNeighbours.begin(); neighbourIter != vertexNeighbours.end(); neighbourIter++){
			unsigned long neighbourVertex = *neighbourIter;
			if (vertexComponent.find(neighbourVertex) == vertexComponent.end()){
				vertexComponent[neighbourVertex] = currentComponent;
				visitingQueue.push(neighbourVertex);
				componentElements.insert(neighbourVertex);

			}
			else if (vertexComponent[neighbourVertex] == currentComponent){}
			else{
				printf("Unexpected Condition On A Connected Component. Expected %d. Obtained %d\n", currentComponent, vertexComponent[neighbourVertex]);
				return 0;
			}
		}
	}
	return 1;
}

int ComputeComponentBoundary(unsigned long vIndex, std::unordered_map<unsigned long, std::unordered_set<unsigned long>> & neighbours, std::unordered_set<unsigned long> & boundaryFaces){
	std::queue<unsigned long > visitingQueue;
	visitingQueue.push(vIndex);
	std::unordered_set<unsigned long > componentElements;
	componentElements.insert(vIndex);

	while (!visitingQueue.empty()){
		unsigned long currentVertex = visitingQueue.front();
		visitingQueue.pop();
		//Add faces
		unsigned long axialIndices[3];
		unsigned long _c;
		GetVoxelElementKeyIndices(currentVertex, axialIndices, _c);
		for (unsigned long c = 0; c < 3; c++){
			for (unsigned long offset = 0; offset < 2; offset++){
				unsigned long faceAxialIndices[3] = { axialIndices[0], axialIndices[1], axialIndices[2] };
				if (offset) faceAxialIndices[c]++;
				unsigned long faceKey = SetVoxelElementKey(faceAxialIndices, c);
				if (boundaryFaces.find(faceKey) == boundaryFaces.end()) boundaryFaces.insert(faceKey);
				else boundaryFaces.erase(faceKey);
			}
		}
		//Add neighbours to queue
		const std::unordered_set<unsigned long> & vertexNeighbours = neighbours[currentVertex];
		for (auto neighbourIter = vertexNeighbours.begin(); neighbourIter != vertexNeighbours.end(); neighbourIter++){
			unsigned long neighbourVertex = *neighbourIter;
			if (componentElements.find(neighbourVertex) == componentElements.end()){
				componentElements.insert(neighbourVertex);
				visitingQueue.push(neighbourVertex);
			}
		}
	}
	return 1;
}

int FillRegionHoles(const int gridRes,
	std::unordered_map<unsigned long, FaceMap> & sourceFaceMaps, std::unordered_map<unsigned long, EdgeIntersections> & sourceEdgeIntersections, std::unordered_map<unsigned long long, SamplePoint> & sourceVertexSample, ComponentMap & sourceRefMap,
	std::unordered_map<unsigned long, FaceMap> & targetFaceMaps, std::unordered_map<unsigned long, EdgeIntersections> & targetEdgeIntersections, std::unordered_map<unsigned long long, SamplePoint> & targetVertexSample, ComponentMap & targetRefMap,
	std::unordered_set<unsigned long> & voxelKeys, std::unordered_set<unsigned long> & addedVoxelKeys, const int maxComponentSize){

	std::unordered_map<unsigned long, bool> addBoundaryFaces = sourceRefMap.boundaryFaceOrientation;

	std::unordered_map<unsigned long, std::unordered_set<unsigned long>> voxelAdjacency;

	for (auto voxelIter = voxelKeys.begin(); voxelIter != voxelKeys.end(); voxelIter++){
		unsigned long currentVoxelKey = *voxelIter;
		if (addedVoxelKeys.find(currentVoxelKey) == addedVoxelKeys.end()){
			unsigned long axialIndices[3];
			unsigned long _c;
			GetVoxelElementKeyIndices(currentVoxelKey, axialIndices, _c);
			unsigned long neigbourVoxelIndices[3] = { axialIndices[0], axialIndices[1], axialIndices[2] };
			for (unsigned long c = 0; c < 3; c++){
				for (unsigned long offset = 0; offset < 2; offset++){
					unsigned long neigbourVoxelIndices[3] = { axialIndices[0], axialIndices[1], axialIndices[2] };
					if (neigbourVoxelIndices[c] > 0 && offset == 0) neigbourVoxelIndices[c]--;
					else if (neigbourVoxelIndices[c] < gridRes - 1 && offset == 1) neigbourVoxelIndices[c]++;
					unsigned long neighbourVoxelKey = SetVoxelElementKey(neigbourVoxelIndices, 0);
					if (voxelKeys.find(neighbourVoxelKey) != voxelKeys.end()){
						if (addedVoxelKeys.find(neighbourVoxelKey) == addedVoxelKeys.end()){
							voxelAdjacency[currentVoxelKey].insert(neighbourVoxelKey);
						}
					}
				}
			}
		}
	} 

	std::unordered_map<unsigned long, int> vertexComponent;
	int currentComponent = -1;

	std::vector<unsigned long> newVoxelsToAdd;

	for (auto voxelIter = voxelKeys.begin(); voxelIter != voxelKeys.end(); voxelIter++){
		unsigned long currentVoxelKey = *voxelIter;
		if (vertexComponent.find(currentVoxelKey) == vertexComponent.end()){
			currentComponent++;
			std::unordered_set<unsigned long> componentElements;
			AddComponent(vertexComponent, currentVoxelKey, currentComponent, voxelAdjacency, componentElements);
			//printf("Component size %d \n", componentElements.size());
			if (componentElements.size() < maxComponentSize){
				//Find boundary faces
				std::unordered_set<unsigned long> boundaryFaces;
				for (auto elementIter = componentElements.begin(); elementIter != componentElements.end(); elementIter++){
					unsigned long currentElement = *elementIter;
					unsigned long axialIndices[3];
					unsigned long _c;
					GetVoxelElementKeyIndices(currentElement, axialIndices, _c);
					for (unsigned long c = 0; c < 3; c++){
						for (unsigned long offset = 0; offset < 2; offset++){
							unsigned long faceAxialIndices[3] = { axialIndices[0], axialIndices[1], axialIndices[2] };
							if (offset) faceAxialIndices[c]++;
							unsigned long faceKey = SetVoxelElementKey(faceAxialIndices, c);
							if (boundaryFaces.find(faceKey) == boundaryFaces.end()) boundaryFaces.insert(faceKey);
							else boundaryFaces.erase(faceKey);
						}
					}
				}
				//Check boundary face similarity
				bool allSimilarFaces = true;
				bool adjacentToBoundary = false;
				for (auto faceIter = boundaryFaces.begin(); faceIter != boundaryFaces.end(); faceIter++){
					unsigned long faceKey = *faceIter;
					if (!SimilarFaces(sourceFaceMaps[faceKey], targetFaceMaps[faceKey])){
						allSimilarFaces = false;
						break;
					}
					if (!adjacentToBoundary){
						if (addBoundaryFaces.find(faceKey) != addBoundaryFaces.end()){
							adjacentToBoundary = true;
						}
					}
				}
				if (allSimilarFaces && adjacentToBoundary){
					if(0) printf("Found valid component!. Size %d \n", componentElements.size());
					for (auto elementIter = componentElements.begin(); elementIter != componentElements.end(); elementIter++){
						newVoxelsToAdd.push_back(*elementIter);
					}
				}
			}
		}
	}

	for (int i = 0; i < newVoxelsToAdd.size(); i++){
		unsigned long currentVoxel = newVoxelsToAdd[i];
		unsigned long axialIndices[3];
		unsigned long _c;
		GetVoxelElementKeyIndices(currentVoxel, axialIndices, _c);
		if (!AddVoxel(sourceRefMap, sourceFaceMaps, sourceEdgeIntersections, sourceVertexSample, targetRefMap, targetFaceMaps, targetEdgeIntersections, targetVertexSample, axialIndices, false, true)) return 0;
		addedVoxelKeys.insert(currentVoxel);
	}
	return 1;
}

int TraverseFaceEdges(FaceMap & faceMap, std::unordered_map<unsigned long long, unsigned long long> & adjacentVertexMap, std::vector<unsigned long long> & vertexSequence, bool reverseOrder){
	vertexSequence.clear();
	if (faceMap.forwardVertex.size() > 1){
		printf("Non simple face! \n");
		return 0;
	}
	else if (faceMap.forwardVertex.size() == 1){
		std::vector<unsigned long long> _vertexSequence;
		unsigned long long startVertex = (*faceMap.forwardVertex.begin()).first;
		unsigned long long endVertex = (*faceMap.forwardVertex.begin()).second;
		unsigned long long currentVertex = startVertex;
		do{
			_vertexSequence.push_back(currentVertex);
			currentVertex = adjacentVertexMap[currentVertex];
		} while (currentVertex != endVertex);
		_vertexSequence.push_back(currentVertex);

		if (reverseOrder){
			int length = _vertexSequence.size();
			vertexSequence.resize(length);
			for (int i = 0; i < length; i++)vertexSequence[i] = _vertexSequence[length - 1 - i];
		}
		else vertexSequence = _vertexSequence;
	}
	return 1;
}


void SetLoopPoints(std::unordered_map<unsigned long long, unsigned long long> & forwardMap, std::unordered_map<unsigned long long, SamplePoint> & vertexSample, const SimpleMesh & mesh, std::vector<std::vector<Point3D<double>>> & loopPoints){
	std::vector<unsigned long long> loopRepresentatives;
	LoopRepresentatives(forwardMap, loopRepresentatives);

	printf("Total loops = %d \n", loopRepresentatives.size());

	loopPoints.clear();
	loopPoints.resize(loopRepresentatives.size());
	for (int i = 0; i < loopRepresentatives.size(); i++){
		unsigned long long startVertex = loopRepresentatives[i];
		unsigned long long currentVertex = startVertex;
		do{
			SamplePoint sample = vertexSample[currentVertex];
			Point3D<double> pos = mesh.vertices[mesh.triangles[sample.tIndex][0]] * (1.0 - sample.baricentricCoord[0] - sample.baricentricCoord[1]) + mesh.vertices[mesh.triangles[sample.tIndex][1]] * sample.baricentricCoord[0] + mesh.vertices[mesh.triangles[sample.tIndex][2]] * sample.baricentricCoord[1];
			loopPoints[i].push_back(pos);
			currentVertex = forwardMap[currentVertex];
		} while (currentVertex != startVertex);
		printf("Loop %d size = %d \n", i, loopPoints[i].size());
	}
}

int InitializeVoxelKeys(std::unordered_map<unsigned long, FaceMap> & faceMaps, std::unordered_set<unsigned long> & voxelKeys){
	
	for (auto faceIter = faceMaps.begin(); faceIter != faceMaps.end(); faceIter++){
		unsigned long faceKey = (*faceIter).first;
		unsigned long axialIndices[3];
		unsigned long cIndex;
		GetVoxelElementKeyIndices(faceKey, axialIndices, cIndex);
		unsigned long voxelKey = SetVoxelElementKey(axialIndices, 0);
		voxelKeys.insert(voxelKey);
		axialIndices[cIndex]--;
		voxelKey = SetVoxelElementKey(axialIndices, 0);
		voxelKeys.insert(voxelKey);
	}
	printf("Intersected Voxels %d \n", voxelKeys.size());
	return 1;
} 

unsigned long long SetWedgeEdgeKey(const unsigned long i0, const unsigned long i1){
	return ((((unsigned long long)i0 << 32) & 0xFFFFFFFF00000000) | ((unsigned long long)i1 & 0x00000000FFFFFFFF));
}

void GetWedgeEdgeIndices(unsigned long long key, unsigned long & i0, unsigned long & i1){
	i1 = static_cast<unsigned long>(key & 0x00000000FFFFFFFF);
	i0 = static_cast<unsigned long>((key >> 32) & 0x00000000FFFFFFFF);
}

//Map wedges to triangles
int SetWedgeKeyMap(std::unordered_map<unsigned long long, unsigned long > & inverseEdgeKeyMap, std::unordered_map<unsigned long long,int> & wedgeMap, const SimpleMesh & mesh){
	unsigned long edgeCounter = 0;

	for (int i = 0; i < mesh.triangles.size(); i++){
		unsigned long reducedEdgeKey[3];
		for (int j = 0; j < 3; j++){
			int v[2] = { mesh.triangles[i][j], mesh.triangles[i][(j + 1) % 3] };
			if (v[0]> v[1]) std::swap(v[0], v[1]);
			unsigned long long  edgeKey = SetMeshEdgeKey(v[0], v[1]);
			if (inverseEdgeKeyMap.find(edgeKey) == inverseEdgeKeyMap.end()){
				printf("Not found edge key! \n");
				return 0;
			}
			else{
				reducedEdgeKey[j] = inverseEdgeKeyMap[edgeKey];
			}
		}
		for (int j = 0; j < 3; j++){
			unsigned long long  wedgeKey = SetWedgeEdgeKey(reducedEdgeKey[j], reducedEdgeKey[(j + 1) % 3]);
			wedgeMap[wedgeKey] = i;
			wedgeKey = SetWedgeEdgeKey(reducedEdgeKey[(j + 1) % 3], reducedEdgeKey[j]);
			wedgeMap[wedgeKey] = i;
		}
	}
	return 1;
}

double SegmentVoxelFaceIntersection(Point3D<double> vertices[2], unsigned long voxelFaceKey, const double voxelWidth, bool knownIntersection){
	unsigned long axialIndices[3];
	unsigned long cIndex;
	GetVoxelElementKeyIndices(voxelFaceKey, axialIndices, cIndex);

	double axialStart = vertices[0][cIndex];
	double axialEnd =vertices[1][cIndex];
	double alpha = (double(axialIndices[cIndex]) * voxelWidth - axialStart) / (axialEnd - axialStart);
	if (alpha > 0.0  && alpha < 1.0){
		Point3D<double> isect = vertices[0] * (1.0 - alpha) + vertices[1]* alpha;
		if (isect[(cIndex + 1) % 3] > double(axialIndices[(cIndex + 1) % 3]) * voxelWidth && isect[(cIndex + 1) % 3] <  double(axialIndices[(cIndex + 1) % 3] + 1) * voxelWidth &&
			isect[(cIndex + 2) % 3] >  double(axialIndices[(cIndex + 2) % 3]) * voxelWidth && isect[(cIndex + 2) % 3] < double(axialIndices[(cIndex + 2) % 3] + 1) * voxelWidth){
			return alpha;
		}
	}
	else{
		if (knownIntersection) printf("ERROR: voxel face intersection not found!\n");
		return -1.0;
	}
}

class IndexedTriangleEdge{
public:
	int triangleIndex;
	int edgeOffset;
};

int SetIndexedTriangleEdgeMap(std::unordered_map<unsigned long long, IndexedTriangleEdge > & indexedTriangleEdgeMap, const SimpleMesh & mesh){

	for (int i = 0; i < mesh.triangles.size(); i++){
		for (int j = 0; j < 3; j++){
			unsigned long long edgeKey = SetMeshEdgeKey(mesh.triangles[i][j], mesh.triangles[i][(j + 1) % 3]);
			if (indexedTriangleEdgeMap.find(edgeKey) == indexedTriangleEdgeMap.end()){
				IndexedTriangleEdge indexedTriangleEdge;
				indexedTriangleEdge.triangleIndex = i;
				indexedTriangleEdge.edgeOffset = j;
				indexedTriangleEdgeMap[edgeKey] = indexedTriangleEdge;
			}
			else{
				printf("Edge adjacent to multiple triangles! \n");
				return 0;
			}
		}
	}
	return 1;
}


int InitalizeTriangleMaps(const int depth, const SimpleMesh & mesh, std::unordered_map<int, TriangleMap> & triangleMaps, std::unordered_map<unsigned long, FaceMap> & faceMaps, std::unordered_map<unsigned long long, unsigned long long> adjacentVertexMap[3], std::unordered_map<unsigned long, bool> boundaryFaceOrientation){
	int res = 2 << depth;
	double voxelWidth = 1.0 / double(res);
	
	std::unordered_map<unsigned long long, unsigned long > inverseEdgeKeyMap;
	std::vector<unsigned long long> edgeKeyMap;
	SetMeshEdgeKeyMap(inverseEdgeKeyMap, edgeKeyMap, mesh);

	std::unordered_map<unsigned long long, int>  wedgeMap;
	SetWedgeKeyMap(inverseEdgeKeyMap, wedgeMap, mesh);

	std::unordered_map<unsigned long long, IndexedTriangleEdge > indexedTriangleEdgeMap;
	SetIndexedTriangleEdgeMap(indexedTriangleEdgeMap, mesh);

	for (auto faceIter = boundaryFaceOrientation.begin(); faceIter != boundaryFaceOrientation.end(); faceIter++){
		unsigned long faceKey = (*faceIter).first;
		bool orientation = (*faceIter).second;
		FaceMap & sourceFace = faceMaps[faceKey];
		unsigned long axialIndices[3];
		unsigned long cIndex;
		GetVoxelElementKeyIndices(faceKey, axialIndices, cIndex);

		std::vector<unsigned long long> vertexSequence;
		if (!TraverseFaceEdges(sourceFace, adjacentVertexMap[cIndex], vertexSequence, !orientation)) return 0;

		if (vertexSequence.size() > 1){
			unsigned long long perviousVertex = vertexSequence[0];
			unsigned long p_vElement, p_tElement, p_intersectionType;
			GetIntersectionKey(perviousVertex, p_tElement, p_vElement, p_intersectionType);
			int previousVertexTriangleIndex = (p_intersectionType == 1) ? p_tElement % (2 << 29) : -1;
			for (int vertexIter = 1; vertexIter < vertexSequence.size(); vertexIter++){
				unsigned long long currentVertex = vertexSequence[vertexIter];
				unsigned long c_vElement, c_tElement, c_intersectionType;
				GetIntersectionKey(currentVertex, c_tElement, c_vElement, c_intersectionType);
				int currentVertexTriangleIndex = -1;
				if (c_intersectionType == 1){
					currentVertexTriangleIndex = c_tElement % (2 << 29);
				}
				else{//Compute intersection time
					unsigned long long edgeKey = edgeKeyMap[c_tElement];
					unsigned long i0, i1;
					GetMeshEdgeIndices(edgeKey, i0, i1);
					Point3D<double> vertices[2] = { mesh.vertices[i0], mesh.vertices[i1] };
					double alpha = SegmentVoxelFaceIntersection(vertices, c_vElement, voxelWidth, true);
					if (alpha <0 || alpha > 1.0){
						printf("Not found segment voxel face intersection!\n");
						return 0;
					}
					{
						TriangleMap & triangleMap = triangleMaps[indexedTriangleEdgeMap[edgeKey].triangleIndex];
						triangleMap.intersectionTime.push_back(std::pair<double, unsigned long long>(double(indexedTriangleEdgeMap[edgeKey].edgeOffset) + alpha, currentVertex));
					}
					{
						unsigned long long oppositeEdgeKey = SetMeshEdgeKey(i1, i0);
						TriangleMap & triangleMap = triangleMaps[indexedTriangleEdgeMap[oppositeEdgeKey].triangleIndex];
						triangleMap.intersectionTime.push_back(std::pair<double, unsigned long long>(double(indexedTriangleEdgeMap[oppositeEdgeKey].edgeOffset) + (1.0 - alpha), currentVertex));
					}
				}
				if (previousVertexTriangleIndex != -1 && currentVertexTriangleIndex != -1){
					if (previousVertexTriangleIndex != currentVertexTriangleIndex){
						printf("Unexpected condition: interior points must belong to same triangle ; %d %d! \n", previousVertexTriangleIndex, currentVertexTriangleIndex);
						return 0;
					}
				}

				if (previousVertexTriangleIndex != -1){
					TriangleMap & triangleMap = triangleMaps[previousVertexTriangleIndex];
					triangleMap.interiorEdges.push_back(std::pair<unsigned long long, unsigned long long>(perviousVertex, currentVertex));
				}
				else if (currentVertexTriangleIndex != -1){
					TriangleMap & triangleMap = triangleMaps[currentVertexTriangleIndex];
					triangleMap.interiorEdges.push_back(std::pair<unsigned long long, unsigned long long>(perviousVertex, currentVertex));
				}
				else{//triangle edge to triangle edge segment
					unsigned long long wedgeKey = SetWedgeEdgeKey(p_tElement, c_tElement);
					if (wedgeMap.find(wedgeKey) == wedgeMap.end()){
						printf("Wedge key not found! \n");
						return 0;
					}
					else{
						int tIndex = wedgeMap[wedgeKey];
						TriangleMap & triangleMap = triangleMaps[tIndex];
						triangleMap.interiorEdges.push_back(std::pair<unsigned long long, unsigned long long>(perviousVertex, currentVertex));
					}
				}

				p_vElement = c_vElement;
				p_tElement = c_tElement;
				p_intersectionType = c_intersectionType;
				previousVertexTriangleIndex = currentVertexTriangleIndex;
				perviousVertex = currentVertex;
			}
		}
	}
}

#define ANSI_DECLARATORS
#define TRILIBRARY
#define NO_TIMER

#include <Triangulation/triangle.c>
void ConstrainedTriangulation(const std::vector<Point2D<double>> & vertices, const std::vector<std::pair<int, int>> & constrainedEdges, const std::vector<bool> & isBoundaryVertex, std::vector<TriangleIndex> & outputTriangles){
	
	struct triangulateio in,out;

	/* Define input points. */

	in.numberofpoints = vertices.size();
	
	in.pointlist = (REAL *)malloc(in.numberofpoints * 2 * sizeof(REAL));
	in.pointmarkerlist = (int *)malloc(in.numberofpoints * sizeof(int));
	for (int i = 0; i < in.numberofpoints; i++) {
		in.pointlist[2 * i] = vertices[i][0];
		in.pointlist[2 * i + 1] = vertices[i][1];
		in.pointmarkerlist[i] = isBoundaryVertex[i] ? 1 : 0; //Check boundary markers documentation
	}

	in.numberofsegments = constrainedEdges.size();
	in.segmentlist = (int *)malloc(in.numberofsegments * 2 * sizeof(int));
	for (int i = 0; i < in.numberofsegments; i++) {
		in.segmentlist[2 * i] = constrainedEdges[i].first + 1;
		in.segmentlist[2 * i + 1] = constrainedEdges[i].second + 1;
	}
	in.numberofholes = 0;
	in.numberofregions = 0;
	in.numberofpointattributes = 0;
	in.segmentmarkerlist = (int *)NULL;

	out.pointlist = (REAL *)NULL;
	out.trianglelist = (int *)NULL;
	out.segmentlist = (int *)NULL;
	out.pointmarkerlist = (int *)NULL;
	out.triangleattributelist = (REAL *)NULL;
	out.segmentmarkerlist = (int *)NULL;
	/* Refine the triangulation according to the attached */
	/*   triangle area constraints.                       */

	triangulate("pQ", &in, &out, (struct triangulateio *) NULL);

	outputTriangles.resize(out.numberoftriangles);
	for (int i = 0; i < out.numberoftriangles; i++) outputTriangles[i] = TriangleIndex(out.trianglelist[3 * i] - 1, out.trianglelist[3 * i + 1] - 1, out.trianglelist[3 * i + 2] - 1);

	free(in.pointlist);
	free(in.pointmarkerlist);
	free(in.segmentlist);
	free(in.segmentmarkerlist);
	
	free(out.pointlist);
	free(out.trianglelist);
	free(out.segmentlist);
	free(out.pointmarkerlist);
	free(out.triangleattributelist);
	free(out.segmentmarkerlist);
}
#undef REAL
#undef dest
#undef ANSI_DECLARATORS
#undef TRILIBRARY
#undef NO_TIMER

unsigned long long GetNonReplicatedKey(unsigned long long inputKey){
	unsigned long vElement, tElement, intersectionType;
	GetIntersectionKey(inputKey, tElement, vElement, intersectionType);
	if (intersectionType == 1){//Triangle interior intersection
		tElement = tElement % (2 << 29);
		return SetIntersectionKey(tElement, vElement, intersectionType);
	}
	else{
		return inputKey;
	}
	
}

int ProcessTriangleMap(const SimpleMesh & mesh, std::unordered_map<int, TriangleMap> & triangleMaps, std::unordered_map<unsigned long long, SamplePoint> & vertexSample, std::vector<unsigned long long > & subdivisionTriangles){
	
	printf("Triangle to subudivide %d \n", triangleMaps.size());
	
	for (auto tIter = triangleMaps.begin(); tIter != triangleMaps.end(); tIter++){
		int tIndex = (*tIter).first;
		TriangleMap & tMap = (*tIter).second;
		
		//Reduce vertices

		std::unordered_map< unsigned long long, int> reducedMap;
		std::vector<unsigned long long> invReducedMap;
		int vertexCounter = 0;

		for (int segIter = 0; segIter < tMap.interiorEdges.size(); segIter++){
			unsigned long long prevVertex = tMap.interiorEdges[segIter].first;
			unsigned long long nextVertex = tMap.interiorEdges[segIter].second;
			if (reducedMap.find(prevVertex) == reducedMap.end()){
				invReducedMap.push_back(prevVertex);
				reducedMap[prevVertex] = vertexCounter;
				vertexCounter++;
			}
			if (reducedMap.find(nextVertex) == reducedMap.end()){
				invReducedMap.push_back(nextVertex);
				reducedMap[nextVertex] = vertexCounter;
				vertexCounter++;
			}
		}

		for (int k = 0; k < 3; k++){
			unsigned long long vIndex = mesh.triangles[tIndex][k];
			invReducedMap.push_back(vIndex);
			reducedMap[vIndex] = vertexCounter;
			vertexCounter++;
		}

		//Set planar vertices
		//Map 3D triangle point to the unit 2D right triangle 
		std::vector<Point2D<double>> planarVertices(vertexCounter);
		std::vector<Point3D<double>> embeddingVertices(vertexCounter);

		Point3D<double> origin = mesh.vertices[mesh.triangles[tIndex][0]];
		Point3D<double> d[2] = { mesh.vertices[mesh.triangles[tIndex][1]] - mesh.vertices[mesh.triangles[tIndex][0]], mesh.vertices[mesh.triangles[tIndex][2]] - mesh.vertices[mesh.triangles[tIndex][0]] };
		SquareMatrix<double, 2> metric_tensor;
		for (int k = 0; k < 2; k++)for (int l = 0; l < 2; l++) metric_tensor(k, l) = Point3D<double>::Dot(d[k], d[l]);
		SquareMatrix<double, 2> inverse_metric = metric_tensor.inverse();
		for (int i = 0; i < vertexCounter; i++){
			if (vertexSample.find(invReducedMap[i]) == vertexSample.end()){
				printf("Not found vertex with key %llu \n", invReducedMap[i]);
				return 0;
			}
			SamplePoint sample = vertexSample[invReducedMap[i]];
			Point3D<double> vertex3DPos = mesh.vertices[mesh.triangles[sample.tIndex][0]] * (1.0 - sample.baricentricCoord[0] - sample.baricentricCoord[1]) + mesh.vertices[mesh.triangles[sample.tIndex][1]] * sample.baricentricCoord[0] + mesh.vertices[mesh.triangles[sample.tIndex][2]] * sample.baricentricCoord[1];
			Point2D<double> rhs;
			for (int k = 0; k < 2; k++) rhs[k] = Point3D<double>::Dot(vertex3DPos - origin, d[k]);
			planarVertices[i] = inverse_metric*rhs;
			embeddingVertices[i] = vertex3DPos;
		}
		
		//Set edges
		std::vector<std::pair<int,int>> constrainedEdges;
		std::vector<bool> isBoundaryVertex(vertexCounter, false);

		tMap.intersectionTime.push_back(std::pair<double, unsigned long long>(0.0, mesh.triangles[tIndex][0]));
		tMap.intersectionTime.push_back(std::pair<double, unsigned long long>(1.0, mesh.triangles[tIndex][1]));
		tMap.intersectionTime.push_back(std::pair<double, unsigned long long>(2.0, mesh.triangles[tIndex][2]));
		std::sort(tMap.intersectionTime.begin(), tMap.intersectionTime.end());
		
		int boundaryVertexCount = tMap.intersectionTime.size();

		for (int i = 0; i < boundaryVertexCount; i++){
			unsigned long long prevVertex = tMap.intersectionTime[i].second;
			unsigned long long nextVertex = tMap.intersectionTime[(i + 1) % boundaryVertexCount].second;
			int reducedPrevVertex = reducedMap[prevVertex];
			int reducedNextVertex = reducedMap[nextVertex];
			constrainedEdges.push_back(std::pair<int, int>(reducedPrevVertex, reducedNextVertex));
			isBoundaryVertex[reducedPrevVertex] = true;
		}

		//Interior edges
		for (int segIter = 0; segIter < tMap.interiorEdges.size(); segIter++){
			unsigned long long prevVertex = tMap.interiorEdges[segIter].first;
			unsigned long long nextVertex = tMap.interiorEdges[segIter].second;
			int reducedPrevVertex = reducedMap[prevVertex];
			int reducedNextVertex = reducedMap[nextVertex];
			constrainedEdges.push_back(std::pair<int, int>(reducedPrevVertex, reducedNextVertex));
		}

		std::vector<TriangleIndex> triangulation;
		ConstrainedTriangulation(planarVertices, constrainedEdges,isBoundaryVertex, triangulation);

		for (int i = 0; i < triangulation.size(); i++){
			for (int k = 0; k < 3; k++){
				unsigned long long originalVIndex = invReducedMap[triangulation[i][k]];
				subdivisionTriangles.push_back(originalVIndex);
			}
		}

		if (0){
			SimpleMesh triangulationMesh;
			triangulationMesh.vertices = embeddingVertices;
			triangulationMesh.triangles = triangulation;
			char outputMesh[256];
			sprintf(outputMesh, "Triangulation-T%06d.ply", tIndex);
			WriteSimpleMesh(triangulationMesh, outputMesh);
		}
	}
}

//Initial mesh is already scaled and centralized to [0,1]^3
int InitializeFaceMaps(const SimpleMesh & mesh, int depth, std::unordered_map<unsigned long, FaceMap> & faceMaps, std::unordered_map<unsigned long, EdgeIntersections> & edgeIntersections, std::unordered_map<unsigned long long, SamplePoint> & vertexSample, std::unordered_map<unsigned long long, unsigned long long> adjacentVertexMap[3]){
	
	int res = 2 << depth;
	printf("Grid Resolution %d \n", res);
	double voxelWidth = 1.0 / double(res);

	std::unordered_map<unsigned long long, unsigned long > inverseEdgeKeyMap;
	std::vector<unsigned long long> edgeKeyMap;
	SetMeshEdgeKeyMap(inverseEdgeKeyMap, edgeKeyMap, mesh);

	std::vector<ChainEdge> edges[3];
	std::vector<unsigned int> chains[3];
	std::unordered_map<unsigned long long, unsigned int> vertexToEdgeMap[3];
	
	for (unsigned long c = 0; c < 3; c++){

		printf("Channel %d : \n", c);

		std::unordered_map<unsigned long long, unsigned long long> & _adjacentVertexMap = adjacentVertexMap[c];

		for (int t = 0; t < mesh.triangles.size(); t++){
			//Compute triangle BBox
			Point3D< double > min_v = mesh.vertices[mesh.triangles[t][0]];
			Point3D< double > max_v = mesh.vertices[mesh.triangles[t][0]];
			for (int k = 1; k < 3; k++) for (int c = 0; c < 3; c++){
				min_v[c] = std::min<double>(min_v[c], mesh.vertices[mesh.triangles[t][k]][c]);
				max_v[c] = std::max<double>(max_v[c], mesh.vertices[mesh.triangles[t][k]][c]);
			}
			SquareMatrix<double, 2> metric_tensor;
			Point3D<double> d[2] = { mesh.vertices[mesh.triangles[t][1]] - mesh.vertices[mesh.triangles[t][0]], mesh.vertices[mesh.triangles[t][2]] - mesh.vertices[mesh.triangles[t][0]] };
			for (int k = 0; k < 2; k++)for (int l = 0; l < 2; l++) metric_tensor(k, l) = Point3D<double>::Dot(d[k], d[l]);
			SquareMatrix<double, 2> inverse_metric = metric_tensor.inverse();

			//Compute normal and level
			Point3D<double> normal = Point3D<double>::CrossProduct(mesh.vertices[mesh.triangles[t][1]] - mesh.vertices[mesh.triangles[t][0]], mesh.vertices[mesh.triangles[t][2]] - mesh.vertices[mesh.triangles[t][0]]);
			normal /= Point3D<double>::Length(normal);
			Point3D<double> baricenter = (mesh.vertices[mesh.triangles[t][0]] + mesh.vertices[mesh.triangles[t][1]] + mesh.vertices[mesh.triangles[t][2]]) / 3.0;
			double level = Point3D<double>::Dot(baricenter, normal);

			//Compute candidate face indices 
			unsigned long min_axial_index[3];
			unsigned long max_axial_index[3];

			for (int k = 0; k < 3; k++){
				min_axial_index[k] = static_cast<unsigned long>(floor(min_v[k] * double(res)));
				max_axial_index[k] = static_cast<unsigned long>(ceil(max_v[k] * double(res)));
			}

			for (unsigned long iter0 = min_axial_index[0]; iter0 <= max_axial_index[0]; iter0++){
				for (unsigned long iter1 = min_axial_index[1]; iter1 <= max_axial_index[1]; iter1++){
					for (unsigned long iter2 = min_axial_index[2]; iter2 <= max_axial_index[2]; iter2++){
						//Compute chain edge intersection
						unsigned long axialIndices[3] = { iter0, iter1, iter2 };
						if (TriangleVoxelFaceIntersection(inverse_metric, normal, level, mesh, t, inverseEdgeKeyMap, edgeKeyMap, axialIndices, c, res, _adjacentVertexMap, vertexSample, edgeIntersections) == -1){ return -1; }
					}
				}
			}
		}

		printf("\t Edges %d \n", _adjacentVertexMap.size());

		std::vector<ChainEdge> & _edges = edges[c];
		_edges.resize(_adjacentVertexMap.size());

		unsigned int edgeCounter = 0;
		unsigned int iterCounter = 0;
		std::unordered_set<unsigned long long> alreadyVisitedVertex;
		std::unordered_map<unsigned long long, unsigned int> & _vertexToEdgeMap = vertexToEdgeMap[c];
		std::vector<unsigned int> & _chains = chains[c];
		unsigned int chainCount = 0;

		for (auto iter = _adjacentVertexMap.begin(); iter != _adjacentVertexMap.end(); iter++){
			iterCounter++;
			unsigned long long sourceVertex = (*iter).first;
			if (alreadyVisitedVertex.find(sourceVertex) == alreadyVisitedVertex.end()){
				unsigned long long currentVertex = sourceVertex;
				bool terminate = false;
				unsigned int startEdgeCounter = edgeCounter;
				do{
					if (alreadyVisitedVertex.find(currentVertex) == alreadyVisitedVertex.end()){
						alreadyVisitedVertex.insert(currentVertex);
						auto mappedVertex = _adjacentVertexMap.find(currentVertex);
						if (mappedVertex != _adjacentVertexMap.end()){
							unsigned long long nextVertex = (*mappedVertex).second;
							_edges[edgeCounter].sourceVertex = currentVertex;
							_edges[edgeCounter].targetVertex = nextVertex;
							_edges[edgeCounter].chainIndex = chainCount;
							if(!SetEdgeProperties(_edges[edgeCounter], c)) return -1;

							_vertexToEdgeMap[currentVertex] = edgeCounter;
							if (edgeCounter > startEdgeCounter){
								_edges[edgeCounter].previousEdge = edgeCounter - 1;
								_edges[edgeCounter - 1].nextEdge = edgeCounter;
							}
							edgeCounter++;
							currentVertex = nextVertex;
						}
						else{
							printf("Vertex to dead end (%d)  %llu \n", c, currentVertex);
							return -1;
						}
					}
					else{
						if (edgeCounter > startEdgeCounter){
							if (_edges[_vertexToEdgeMap[currentVertex]].previousEdge != -1){
								printf("Multiple incoming mappings to the same vertex (%d) %llu! \n", c, currentVertex);
								return -1;
							}
							else{
								_edges[_vertexToEdgeMap[currentVertex]].previousEdge = edgeCounter - 1;
								_edges[edgeCounter -1].nextEdge = _vertexToEdgeMap[currentVertex];
							}
						}
						terminate = true;
					}
				} while (!terminate);

				//Add Chain
				_chains.push_back(startEdgeCounter);
				chainCount++;
			}
		}

		printf("\t Chains %d \n", chainCount);

		for (int chainIter = 0; chainIter <_chains.size(); chainIter++){
			//Traverse the current chain
			unsigned int startEdge = _chains[chainIter];
			unsigned int currentEdge = startEdge;

			//Find first edge intersection
			bool foundIntersection = false;

			do{
				const ChainEdge & currentEdgeData = _edges[currentEdge];
				if (currentEdgeData.sourceVertexType != 4) foundIntersection = true;
				else currentEdge = currentEdgeData.nextEdge;
			} while (!foundIntersection && currentEdge != startEdge);

			if (!foundIntersection){
				unsigned long long sourceVertex = _edges[startEdge].sourceVertex;
				unsigned long sourceFaceKey = _edges[startEdge].faceKey;
				faceMaps[sourceFaceKey].forwardVertex[sourceVertex] = sourceVertex;
				faceMaps[sourceFaceKey].backwardVertex[sourceVertex] = sourceVertex;
				faceMaps[sourceFaceKey].vertexType[sourceVertex] = _edges[startEdge].sourceVertexType;
				if (_edges[startEdge].sourceVertexType != 4){
					printf("Unexpected vertex type! \n");
					return 0;
				}
			}

			else{
				startEdge = currentEdge;
				unsigned int sourceEdge = currentEdge;
				unsigned char sourceIntersectionType = _edges[currentEdge].sourceVertexType;
				unsigned long sourceFaceKey = _edges[currentEdge].faceKey;
				do{
					const ChainEdge & currentEdgeData = _edges[currentEdge];
					if (currentEdgeData.targetVertexType != 4){
						if (currentEdgeData.faceKey != sourceFaceKey){
							printf("Source and target face key does not match! \n");
							return -1;
						}

						unsigned long long sourceVertex, targetVertex;
						sourceVertex = _edges[sourceEdge].sourceVertex;
						targetVertex = _edges[currentEdge].targetVertex;
						faceMaps[sourceFaceKey].forwardVertex[sourceVertex] = targetVertex;
						faceMaps[sourceFaceKey].backwardVertex[targetVertex] = sourceVertex;
						faceMaps[sourceFaceKey].vertexType[sourceVertex] = _edges[sourceEdge].sourceVertexType;
						faceMaps[sourceFaceKey].vertexType[targetVertex] = _edges[sourceEdge].targetVertexType;

						currentEdge = currentEdgeData.nextEdge;

						sourceEdge = currentEdge;

						const ChainEdge & nextEdgeData = _edges[currentEdge];

						if (nextEdgeData.faceKey == sourceFaceKey){
							printf("Unexpected equal faces!\n");
							return -1;
						}
						sourceFaceKey = nextEdgeData.faceKey;

						if (nextEdgeData.sourceVertexType == 4){
							printf("Unexpected interior vertex!\n");
							return -1;
						}
						sourceIntersectionType = nextEdgeData.sourceVertexType;
					}
					else currentEdge = currentEdgeData.nextEdge;
				} while (currentEdge != startEdge);
			}
		}
	}

	printf("Intersected Grid Faces %d \n", faceMaps.size());
	return 1;
}


int MeshClipping(const SimpleMesh & mesh, int depth, SimpleMesh & clippedMesh, std::unordered_map<unsigned long long, int> & vertexIndex, std::vector<unsigned long> & triangleVoxelMap){
	
	triangleVoxelMap.clear();
	vertexIndex.clear();
	
	int res = 2 << depth;
	printf("Grid Resolution %d \n", res);
	double voxelWidth = 1.0 / double(res);

	std::unordered_map<unsigned long long, unsigned long > inverseEdgeKeyMap;
	std::vector<unsigned long long> edgeKeyMap;
	SetMeshEdgeKeyMap(inverseEdgeKeyMap, edgeKeyMap, mesh);

	std::vector<TriangleIndex> clippedTriangles;
	std::vector<Point3D<double>> clippedVertices;

	for (int i = 0; i < mesh.vertices.size(); i++){
		vertexIndex[i] = i;
		clippedVertices.push_back(mesh.vertices[i]);
	}
	int lastVertexIndex = mesh.vertices.size();

	for (int t = 0; t < mesh.triangles.size(); t++){
		//Compute triangle BBox
		Point3D< double > min_v = mesh.vertices[mesh.triangles[t][0]];
		Point3D< double > max_v = mesh.vertices[mesh.triangles[t][0]];
		for (int k = 1; k < 3; k++) for (int c = 0; c < 3; c++){
			min_v[c] = std::min<double>(min_v[c], mesh.vertices[mesh.triangles[t][k]][c]);
			max_v[c] = std::max<double>(max_v[c], mesh.vertices[mesh.triangles[t][k]][c]);
		}
		SquareMatrix<double, 2> metric_tensor;
		Point3D<double> d[2] = { mesh.vertices[mesh.triangles[t][1]] - mesh.vertices[mesh.triangles[t][0]], mesh.vertices[mesh.triangles[t][2]] - mesh.vertices[mesh.triangles[t][0]] };
		for (int k = 0; k < 2; k++)for (int l = 0; l < 2; l++) metric_tensor(k, l) = Point3D<double>::Dot(d[k], d[l]);
		SquareMatrix<double, 2> inverse_metric = metric_tensor.inverse();

		//Compute normal and level
		Point3D<double> normal = Point3D<double>::CrossProduct(mesh.vertices[mesh.triangles[t][1]] - mesh.vertices[mesh.triangles[t][0]], mesh.vertices[mesh.triangles[t][2]] - mesh.vertices[mesh.triangles[t][0]]);
		normal /= Point3D<double>::Length(normal);
		Point3D<double> baricenter = (mesh.vertices[mesh.triangles[t][0]] + mesh.vertices[mesh.triangles[t][1]] + mesh.vertices[mesh.triangles[t][2]]) / 3.0;
		double level = Point3D<double>::Dot(baricenter, normal);

		//Compute candidate face indices 
		unsigned long min_axial_index[3];
		unsigned long max_axial_index[3];

		for (int k = 0; k < 3; k++){
			min_axial_index[k] = static_cast<unsigned long>(floor(min_v[k] * double(res)));
			max_axial_index[k] = static_cast<unsigned long>(ceil(max_v[k] * double(res)));
		}


		for (unsigned long iter0 = min_axial_index[0]; iter0 <= max_axial_index[0]; iter0++){
			for (unsigned long iter1 = min_axial_index[1]; iter1 <= max_axial_index[1]; iter1++){
				for (unsigned long iter2 = min_axial_index[2]; iter2 <= max_axial_index[2]; iter2++){
					//Compute chain edge intersection
					unsigned long axialIndices[3] = { iter0, iter1, iter2 };
					std::unordered_map<unsigned long long, unsigned long long> forwardMap;
					std::unordered_map<unsigned long long, Point3D<double>> vertexPos;
					unsigned long voxelKey = SetVoxelElementKey(axialIndices, 0);
					if (TriangleVoxelIntersection(inverse_metric, normal, level, mesh, t, inverseEdgeKeyMap, edgeKeyMap, axialIndices, res, forwardMap,vertexPos) == -1){ return -1; }
					if (forwardMap.size()){
						if (forwardMap.size() < 3){
							printf("Invalid triangle - voxel intersection!\n");
							return 0;
						}
						else{
							if (1){//Check loop conditions
								std::vector<unsigned long long> loopRepresentatives;
								if (!LoopRepresentatives(forwardMap, loopRepresentatives))return 0;
								if (loopRepresentatives.size() != 1){
									printf("triangle voxel intersection with more than a single loop! \n");
									return 0;
								}
							}

							unsigned long long startVertex = (*forwardMap.begin()).first;
							unsigned long long currentVertex = startVertex;
							int loopPos = 0;
							int zeroIndex;
							int pIndex;
							do{
								int vIndex = -1;
								//Add vertices
								if (vertexIndex.find(currentVertex) == vertexIndex.end()){
									vertexIndex[currentVertex] = lastVertexIndex;
									vIndex = lastVertexIndex;
									clippedVertices.push_back(vertexPos[currentVertex]);
									lastVertexIndex++;
								}
								else vIndex = vertexIndex[currentVertex];

								if (loopPos == 0) zeroIndex = vIndex;
								else if (loopPos == 1) pIndex = vIndex;
								else{
									clippedTriangles.push_back(TriangleIndex(zeroIndex, pIndex, vIndex));
									triangleVoxelMap.push_back(voxelKey);
									pIndex = vIndex;
								}
								loopPos++;

								currentVertex = forwardMap[currentVertex];
							} while (currentVertex != startVertex);
						}
					}
				}
			}
		}
	}

	clippedMesh.triangles = clippedTriangles;
	clippedMesh.vertices = clippedVertices;

	return 1;

}

int GetFaceAdjacentTriangles(std::unordered_map<unsigned long, FaceMap> & sourceFaceMaps, unsigned long faceKey, bool orientation, std::unordered_map<unsigned long long, unsigned long long> sourceAdjacentVertexMap[3], std::unordered_map<unsigned long long, int> & sourceClippedVertexIndex, std::unordered_map<unsigned long long, int> & sourceClippedEdgeTriangleMap, std::vector<int> & sourceClippedVertexSequence, std::vector<int> & sourceClippedTriangleSequence, std::vector<int> & sourceClippedOppositeVertexSequence, std::vector<double> & sourceCumLengthVector, SimpleMesh & sourceClippedMesh){
	
	FaceMap & sourceFace = sourceFaceMaps[faceKey];

	unsigned long axialIndices[3];
	unsigned long cIndex;
	GetVoxelElementKeyIndices(faceKey, axialIndices, cIndex);

	std::vector<unsigned long long> sourceFaceVertexSequence;
	if (!TraverseFaceEdges(sourceFace, sourceAdjacentVertexMap[cIndex], sourceFaceVertexSequence, !orientation)) return 0;
	int sourceSequenceLenght = sourceFaceVertexSequence.size();

	if (sourceSequenceLenght > 1){

		sourceClippedVertexSequence.resize(sourceSequenceLenght);
		for (int i = 0; i < sourceSequenceLenght; i++)sourceClippedVertexSequence[i] = sourceClippedVertexIndex[sourceFaceVertexSequence[i]];

		sourceClippedTriangleSequence.resize(sourceSequenceLenght - 1);
		for (int i = 0; i < sourceSequenceLenght - 1; i++){
			unsigned long long edgeKey = SetMeshEdgeKey(sourceClippedVertexSequence[i], sourceClippedVertexSequence[i + 1]);
			if (sourceClippedEdgeTriangleMap.find(edgeKey) == sourceClippedEdgeTriangleMap.end()){
				printf("Not found opposite edge tp %llu \n", edgeKey);
				return 0;
			}
			sourceClippedTriangleSequence[i] = sourceClippedEdgeTriangleMap[edgeKey];
		}

		sourceClippedOppositeVertexSequence.resize(sourceSequenceLenght - 1, -1);
		for (int i = 0; i < sourceSequenceLenght - 1; i++){
			int tIndex = sourceClippedTriangleSequence[i];
			int diffVerticesCount = 0;
			for (int k = 0; k < 3; k++){
				int vIndex = sourceClippedMesh.triangles[tIndex][k];
				if (vIndex != sourceClippedVertexSequence[i] && vIndex != sourceClippedVertexSequence[i + 1]){
					sourceClippedOppositeVertexSequence[i] = vIndex;
				}
				else{
					diffVerticesCount++;
				}
			}
			if (diffVerticesCount != 2){
				printf("Unexpected triangle indices!\n");
				return 0;
			}
		}

		sourceCumLengthVector.resize(sourceSequenceLenght, 0.0);
		double sourceCumLength = 0.0;
		for (int i = 1; i < sourceSequenceLenght; i++){
			double length = Point3D<double>::Length(sourceClippedMesh.vertices[sourceClippedVertexSequence[i]] - sourceClippedMesh.vertices[sourceClippedVertexSequence[i-1]]);
			sourceCumLength += length;
			sourceCumLengthVector[i] = sourceCumLength;
		}
		for (int i = 0; i < sourceSequenceLenght; i++) sourceCumLengthVector[i] /= sourceCumLength;
	}
	else{

	}
	return 1;
}



int IdentifyTrianglesToSplit(const std::vector<TriangleIndex> & triangles, std::unordered_set<int> & trianglesToClip, std::unordered_set<int> & trianglesToSplit, std::unordered_set<unsigned long long> & splitEdgeKeys){
	
	//Initialize edges of triangleToClip
	std::unordered_map<unsigned long long, int> edgeToTriangle;
	
	for (auto triangleIter = trianglesToClip.begin(); triangleIter != trianglesToClip.end(); triangleIter++){
		int t = (*triangleIter);
		for (int k = 0; k < 3; k++){
			unsigned long long edgeKey = SetMeshEdgeKey(triangles[t][k], triangles[t][(k + 1) % 3]);
			if (edgeToTriangle.find(edgeKey) == edgeToTriangle.end()) edgeToTriangle[edgeKey] = t;
			else{
				printf("Non manifold mesh!! \n");
				return 0;
			}
		}
	}


	//Find non clipped triangles adjacent to triangle to clip
	for (int i = 0; i < triangles.size(); i++){
		if (trianglesToClip.find(i) == trianglesToClip.end()){
			for (int k = 0; k < 3; k++){
				unsigned long long edgeKey = SetMeshEdgeKey(triangles[i][(k + 1) % 3], triangles[i][k]);
				if (edgeToTriangle.find(edgeKey) != edgeToTriangle.end()){
					int v[2] = { triangles[i][(k + 1) % 3], triangles[i][k] };
					if (v[0]> v[1]) std::swap(v[0], v[1]);
					unsigned long long orderedEdgeKey = SetMeshEdgeKey(v[0],v[1]);
					splitEdgeKeys.insert(orderedEdgeKey);
					trianglesToSplit.insert(i);
				}
			}
		}
	}
	return 1;
}

int MeshSubsetClippling(const SimpleMesh & mesh, std::unordered_map<int, TriangleMap> & triangleMaps, std::unordered_map<unsigned long long, SamplePoint> & vertexSample, const std::vector<unsigned long long > & subdivisionTriangles, SimpleMesh & clippedMesh, std::unordered_map<unsigned long long, int> & vertexIndex, std::vector<int> & originalToClippedTriangleMap){
	originalToClippedTriangleMap.resize(mesh.triangles.size(), -1);
	vertexIndex.clear();

	int lastVertexIndex = 0;
	int lastTriangleIndex = 0;

	//Add preserved triangles
	for (int t = 0; t < mesh.triangles.size(); t++){
		if (triangleMaps.find(t) == triangleMaps.end()){
			int reducedVIndex[3];
			for (int k = 0; k < 3; k++){
				int unreducedVIndex = mesh.triangles[t][k];
				if (vertexIndex.find(unreducedVIndex) == vertexIndex.end()){
					SamplePoint sample = vertexSample[unreducedVIndex];
					Point3D<double> vertex3DPos = mesh.vertices[mesh.triangles[sample.tIndex][0]] * (1.0 - sample.baricentricCoord[0] - sample.baricentricCoord[1]) + mesh.vertices[mesh.triangles[sample.tIndex][1]] * sample.baricentricCoord[0] + mesh.vertices[mesh.triangles[sample.tIndex][2]] * sample.baricentricCoord[1];
					clippedMesh.vertices.push_back(vertex3DPos);
					vertexIndex[unreducedVIndex] = lastVertexIndex;
					lastVertexIndex++;
				}
				reducedVIndex[k] = vertexIndex[unreducedVIndex];
			}
			clippedMesh.triangles.push_back(TriangleIndex(reducedVIndex[0], reducedVIndex[1], reducedVIndex[2]));
			originalToClippedTriangleMap[t] = lastTriangleIndex;
			lastTriangleIndex++;
		}
	}

	//Add subdivided triangles
	for (int t = 0; t < subdivisionTriangles.size() / 3; t++){
		int reducedVIndex[3];
		for (int k = 0; k < 3; k++){
			unsigned long long unreducedVIndex = subdivisionTriangles[3*t + k];
			if (vertexIndex.find(unreducedVIndex) == vertexIndex.end()){
				SamplePoint sample = vertexSample[unreducedVIndex];
				Point3D<double> vertex3DPos = mesh.vertices[mesh.triangles[sample.tIndex][0]] * (1.0 - sample.baricentricCoord[0] - sample.baricentricCoord[1]) + mesh.vertices[mesh.triangles[sample.tIndex][1]] * sample.baricentricCoord[0] + mesh.vertices[mesh.triangles[sample.tIndex][2]] * sample.baricentricCoord[1];
				clippedMesh.vertices.push_back(vertex3DPos);
				vertexIndex[unreducedVIndex] = lastVertexIndex;
				lastVertexIndex++;
			}
			reducedVIndex[k] = vertexIndex[unreducedVIndex];
		}
		clippedMesh.triangles.push_back(TriangleIndex(reducedVIndex[0], reducedVIndex[1], reducedVIndex[2]));
	}

}

int MeshSubsetClippling(const SimpleMesh & mesh, int depth, SimpleMesh & clippedMesh, std::unordered_map<unsigned long long, int> & vertexIndex, std::unordered_set<int> & trianglesToClip, std::vector<int> & originalToClippedTriangleMap){


	originalToClippedTriangleMap.resize(mesh.triangles.size(), -1);

	//Identify triangles to split
	std::unordered_set<unsigned long long> splitEdgeKeys;
	std::unordered_set<int> trianglesToSplit;
	IdentifyTrianglesToSplit(mesh.triangles, trianglesToClip, trianglesToSplit, splitEdgeKeys);
	std::unordered_map< unsigned long long, std::set<OrderedIntersection,OrderedIntersectionComparison>> splitEdgeIntersections;

	//triangleVoxelMap.clear();
	vertexIndex.clear();

	int res = 2 << depth;
	printf("Grid Resolution %d \n", res);
	double voxelWidth = 1.0 / double(res);

	std::unordered_map<unsigned long long, unsigned long > inverseEdgeKeyMap;
	std::vector<unsigned long long> edgeKeyMap;
	SetMeshEdgeKeyMap(inverseEdgeKeyMap, edgeKeyMap, mesh);

	std::vector<TriangleIndex> clippedTriangles;
	std::vector<Point3D<double>> clippedVertices;

	for (int i = 0; i < mesh.vertices.size(); i++){
		vertexIndex[i] = i;
		clippedVertices.push_back(mesh.vertices[i]);
	}
	int initialNumVertices = mesh.vertices.size();
	int lastVertexIndex = mesh.vertices.size();


	for (auto triangleIter = trianglesToClip.begin(); triangleIter != trianglesToClip.end(); triangleIter++){

		int t = (*triangleIter);

		//Compute triangle BBox
		Point3D< double > min_v = mesh.vertices[mesh.triangles[t][0]];
		Point3D< double > max_v = mesh.vertices[mesh.triangles[t][0]];
		for (int k = 1; k < 3; k++) for (int c = 0; c < 3; c++){
			min_v[c] = std::min<double>(min_v[c], mesh.vertices[mesh.triangles[t][k]][c]);
			max_v[c] = std::max<double>(max_v[c], mesh.vertices[mesh.triangles[t][k]][c]);
		}
		SquareMatrix<double, 2> metric_tensor;
		Point3D<double> d[2] = { mesh.vertices[mesh.triangles[t][1]] - mesh.vertices[mesh.triangles[t][0]], mesh.vertices[mesh.triangles[t][2]] - mesh.vertices[mesh.triangles[t][0]] };
		for (int k = 0; k < 2; k++)for (int l = 0; l < 2; l++) metric_tensor(k, l) = Point3D<double>::Dot(d[k], d[l]);
		SquareMatrix<double, 2> inverse_metric = metric_tensor.inverse();

		//Compute normal and level
		Point3D<double> normal = Point3D<double>::CrossProduct(mesh.vertices[mesh.triangles[t][1]] - mesh.vertices[mesh.triangles[t][0]], mesh.vertices[mesh.triangles[t][2]] - mesh.vertices[mesh.triangles[t][0]]);
		normal /= Point3D<double>::Length(normal);
		Point3D<double> baricenter = (mesh.vertices[mesh.triangles[t][0]] + mesh.vertices[mesh.triangles[t][1]] + mesh.vertices[mesh.triangles[t][2]]) / 3.0;
		double level = Point3D<double>::Dot(baricenter, normal);

		//Compute candidate face indices 
		unsigned long min_axial_index[3];
		unsigned long max_axial_index[3];

		for (int k = 0; k < 3; k++){
			min_axial_index[k] = static_cast<unsigned long>(floor(min_v[k] * double(res)));
			max_axial_index[k] = static_cast<unsigned long>(ceil(max_v[k] * double(res)));
		}


		for (unsigned long iter0 = min_axial_index[0]; iter0 <= max_axial_index[0]; iter0++){
			for (unsigned long iter1 = min_axial_index[1]; iter1 <= max_axial_index[1]; iter1++){
				for (unsigned long iter2 = min_axial_index[2]; iter2 <= max_axial_index[2]; iter2++){
					//Compute chain edge intersection
					unsigned long axialIndices[3] = { iter0, iter1, iter2 };
					std::unordered_map<unsigned long long, unsigned long long> forwardMap;
					std::unordered_map<unsigned long long, Point3D<double>> vertexPos;
					unsigned long voxelKey = SetVoxelElementKey(axialIndices, 0);
					if (TriangleVoxelIntersection(inverse_metric, normal, level, mesh, t, inverseEdgeKeyMap, edgeKeyMap, axialIndices, res, forwardMap, vertexPos, splitEdgeKeys, splitEdgeIntersections) == -1){ return -1; }
					if (forwardMap.size()){
						if (forwardMap.size() < 3){
							printf("Invalid triangle - voxel intersection!\n");
							return 0;
						}
						else{
							if (1){//Check loop conditions
								std::vector<unsigned long long> loopRepresentatives;
								if (!LoopRepresentatives(forwardMap, loopRepresentatives))return 0;
								if (loopRepresentatives.size() != 1){
									printf("triangle voxel intersection with more than a single loop! \n");
									return 0;
								}
							}

							unsigned long long startVertex = (*forwardMap.begin()).first;
							unsigned long long currentVertex = startVertex;
							int loopPos = 0;
							int zeroIndex;
							int pIndex;
							do{
								int vIndex = -1;
								//Add vertices
								if (vertexIndex.find(currentVertex) == vertexIndex.end()){
									vertexIndex[currentVertex] = lastVertexIndex;
									vIndex = lastVertexIndex;
									clippedVertices.push_back(vertexPos[currentVertex]);
									lastVertexIndex++;
								}
								else vIndex = vertexIndex[currentVertex];

								if (loopPos == 0) zeroIndex = vIndex;
								else if (loopPos == 1) pIndex = vIndex;
								else{
									clippedTriangles.push_back(TriangleIndex(zeroIndex, pIndex, vIndex));
									//triangleVoxelMap.push_back(voxelKey);
									pIndex = vIndex;
								}
								loopPos++;

								currentVertex = forwardMap[currentVertex];
							} while (currentVertex != startVertex);
						}
					}
				}
			}
		}
	}

	for (auto splitTriangleIter = trianglesToSplit.begin(); splitTriangleIter != trianglesToSplit.end(); splitTriangleIter++){
		int tIndex = (*splitTriangleIter);
		std::vector<unsigned long long> triangleBoundaryVertices;
		for (int k = 0; k < 3; k++){
			triangleBoundaryVertices.push_back(mesh.triangles[tIndex][k]);
			int v[2] = { mesh.triangles[tIndex][k], mesh.triangles[tIndex][(k + 1) % 3] };
			if (v[0]> v[1]) std::swap(v[0], v[1]);
			unsigned long long orderedEdgeKey = SetMeshEdgeKey(v[0], v[1]);
			if (splitEdgeIntersections.find(orderedEdgeKey) != splitEdgeIntersections.end()){
				std::set<OrderedIntersection, OrderedIntersectionComparison> & edgeIntersections = splitEdgeIntersections[orderedEdgeKey];
				//for (int j = 0; j < edgeIntersections.size(); j++){
				//	triangleBoundaryVertices.push_back(edgeIntersections[j].second);
				//}
				for (auto edgeIntersectionsIter = edgeIntersections.begin(); edgeIntersectionsIter != edgeIntersections.end(); edgeIntersectionsIter++){
					OrderedIntersection isect = *edgeIntersectionsIter;
					triangleBoundaryVertices.push_back(isect.vertex);
				}
			}
		}
		if (triangleBoundaryVertices.size() == 3){
			clippedTriangles.push_back(TriangleIndex(vertexIndex[triangleBoundaryVertices[0]], vertexIndex[triangleBoundaryVertices[1]], vertexIndex[triangleBoundaryVertices[2]]));
			originalToClippedTriangleMap[tIndex] = clippedTriangles.size() - 1;
			//triangleVoxelMap.push_back(-1);
		}
		else if (triangleBoundaryVertices.size() > 3){
			Point3D<double> baricenter = (mesh.vertices[mesh.triangles[tIndex][0]] + mesh.vertices[mesh.triangles[tIndex][1]] + mesh.vertices[mesh.triangles[tIndex][2]]) / 3.0;
			int baricenterReducedIndex = lastVertexIndex;
			int baricenterUnreducedIndex = initialNumVertices + tIndex;
			vertexIndex[baricenterUnreducedIndex] = baricenterReducedIndex;
			clippedVertices.push_back(baricenter);
			lastVertexIndex++;

			for (int j = 0; j < triangleBoundaryVertices.size(); j++){
				clippedTriangles.push_back(TriangleIndex(baricenterReducedIndex, vertexIndex[triangleBoundaryVertices[j]], vertexIndex[triangleBoundaryVertices[(j + 1) % triangleBoundaryVertices.size()]]));
				//triangleVoxelMap.push_back(-1);
			}
 		}
 	}

	for (int t = 0; t < mesh.triangles.size(); t++){
		if (trianglesToClip.find(t) == trianglesToClip.end() && trianglesToSplit.find(t) == trianglesToSplit.end()){
			clippedTriangles.push_back(TriangleIndex(mesh.triangles[t][0], mesh.triangles[t][1], mesh.triangles[t][2]));
			originalToClippedTriangleMap[t] = clippedTriangles.size() - 1;
		}
	}

	clippedMesh.triangles = clippedTriangles;
	clippedMesh.vertices = clippedVertices;

	if (1){
		ColoredMesh coloredMesh;
		coloredMesh.vertices.resize(mesh.triangles.size() * 3);
		coloredMesh.colors.resize(mesh.triangles.size() * 3);
		coloredMesh.triangles.resize(mesh.triangles.size());
		for (int i = 0; i < mesh.triangles.size(); i++){
			coloredMesh.triangles[i] = TriangleIndex(3 * i, 3 * i + 1, 3 * i + 2);
			for (int k = 0; k < 3; k++) coloredMesh.vertices[3 * i + k] = mesh.vertices[mesh.triangles[i][k]];
			if (trianglesToClip.find(i) != trianglesToClip.end() && trianglesToSplit.find(i) != trianglesToSplit.end()){
				printf("Simultaneous split and clip triangle! \n");
				return 0;
			}
			else if(trianglesToClip.find(i) != trianglesToClip.end()){
				for (int k = 0; k < 3; k++) coloredMesh.colors[3 * i + k] = Point3D<double>(255.0,0.0,0.0);
			}
			else if (trianglesToSplit.find(i) != trianglesToSplit.end()){
				for (int k = 0; k < 3; k++) coloredMesh.colors[3 * i + k] = Point3D<double>(0.0, 0.0, 255.0);
			}
			else{
				for (int k = 0; k < 3; k++) coloredMesh.colors[3 * i + k] = Point3D<double>(0.0, 255.0, 0.0);
			}
		}
		WriteColoredMesh(coloredMesh, "Classification.ply");
	}
	return 1;
}

int IdentifyBoundaryAdjacentTriangles(const SimpleMesh & mesh, std::unordered_map<unsigned long, bool> & boundaryFaceOrientation, std::unordered_map<unsigned long, FaceMap> & faceMaps, std::unordered_map<unsigned long long, unsigned long long> adjacentVertexMap[3], std::unordered_set<int> & boundaryTriangles){

	std::unordered_map<unsigned long long, unsigned long > inverseEdgeKeyMap;
	std::vector<unsigned long long> edgeKeyMap;
	SetMeshEdgeKeyMap(inverseEdgeKeyMap, edgeKeyMap, mesh);
	std::unordered_map<unsigned long long, int > edgeKeyTriangleMap;
	SetMeshEdgeKeyTriangleMap(edgeKeyTriangleMap, mesh);

	for (auto faceIter = boundaryFaceOrientation.begin(); faceIter != boundaryFaceOrientation.end(); faceIter++){
		unsigned long faceKey = (*faceIter).first;
		bool orientation = (*faceIter).second;

		FaceMap & faceMap = faceMaps[faceKey];

		unsigned long axialIndices[3];
		unsigned long cIndex;
		GetVoxelElementKeyIndices(faceKey, axialIndices, cIndex);

		std::vector<unsigned long long> vertexSequence;
		if (!TraverseFaceEdges(faceMap, adjacentVertexMap[cIndex], vertexSequence, false)) return 0;
		for (int i = 0; i < vertexSequence.size(); i++){
			unsigned long long vertexKey = vertexSequence[i];
			unsigned long vElement, tElement, intersectionType;
			GetIntersectionKey(vertexKey, tElement, vElement, intersectionType);
			if (intersectionType == 1){//Triangle interior intersection
				int tIndex = tElement % (2 << 29);
				boundaryTriangles.insert(tIndex);
			}
			else if (intersectionType == 0){//Triangle edge intersection
				unsigned long long edgeKey = edgeKeyMap[tElement];
				unsigned long long oppEdgeKey = InvertMeshEdgeKey(edgeKey);
				if (edgeKeyTriangleMap.find(edgeKey) != edgeKeyTriangleMap.end()){
					boundaryTriangles.insert(edgeKeyTriangleMap[edgeKey]);
				}
				if (edgeKeyTriangleMap.find(oppEdgeKey) != edgeKeyTriangleMap.end()){
					boundaryTriangles.insert(edgeKeyTriangleMap[oppEdgeKey]);
				}
			}
		}
	}
	return 1;
}

int IdentifyBoundaryEdges(std::unordered_map<unsigned long, bool> & boundaryFaceOrientation, std::unordered_map<unsigned long, FaceMap> & faceMaps, std::unordered_map<unsigned long long, unsigned long long> adjacentVertexMap[3], std::unordered_map<unsigned long long, int> & vertexIndex, std::unordered_set<unsigned long long> & boundaryEdges, bool boundaryOrientation){
	printf("Num vertices clipped mesh %d \n", vertexIndex.size());
	for (auto faceIter = boundaryFaceOrientation.begin(); faceIter != boundaryFaceOrientation.end(); faceIter++){
		unsigned long faceKey = (*faceIter).first;
		bool faceOrientation = (*faceIter).second;

		FaceMap & faceMap = faceMaps[faceKey];

		unsigned long axialIndices[3];
		unsigned long cIndex;
		GetVoxelElementKeyIndices(faceKey, axialIndices, cIndex);

		std::vector<unsigned long long> vertexSequence;
		if (!TraverseFaceEdges(faceMap, adjacentVertexMap[cIndex], vertexSequence, faceOrientation)) return 0;
		if (vertexSequence.size() > 1){
			for (int i = 0; i < vertexSequence.size() - 1; i++){
				if (vertexIndex.find(vertexSequence[i]) == vertexIndex.end()){
					printf("Not found vertex! \n");
					return 0;
				}
				if (vertexIndex.find(vertexSequence[i + 1]) == vertexIndex.end()){
					printf("Not found vertex! \n");
					return 0;
				}

				unsigned long long edgeKey = boundaryOrientation ? SetMeshEdgeKey(vertexIndex[vertexSequence[i]], vertexIndex[vertexSequence[i + 1]]) : SetMeshEdgeKey(vertexIndex[vertexSequence[i + 1]], vertexIndex[vertexSequence[i]]);
				boundaryEdges.insert(edgeKey);
			}
		}
	}
	return 1;
}


int BoundaryConstrainedPropagation(const SimpleMesh & mesh, std::unordered_set<unsigned long long> & boundaryEdges, std::unordered_set<int> & markedTriangles){

	std::unordered_map<unsigned long long, int > edgeKeyTriangleMap;
	SetMeshEdgeKeyTriangleMap(edgeKeyTriangleMap, mesh);

	std::queue<int> triangleQueue;
	for (auto boundaryEdgeIter = boundaryEdges.begin(); boundaryEdgeIter != boundaryEdges.end(); boundaryEdgeIter++){
		if (edgeKeyTriangleMap.find(*boundaryEdgeIter) != edgeKeyTriangleMap.end()){
			int triangleId = edgeKeyTriangleMap[*boundaryEdgeIter];
			if (markedTriangles.find(triangleId) == markedTriangles.end()){
				markedTriangles.insert(triangleId);
				triangleQueue.push(triangleId);
			}
		}
		else{
			printf("Not found edge! \n");
			return 0;
		}

	}

	while (!triangleQueue.empty()){
		int triangleId = triangleQueue.front();
		triangleQueue.pop();
		for (int k = 0; k < 3; k++){
			unsigned long long  edgeKey = SetMeshEdgeKey(mesh.triangles[triangleId][k], mesh.triangles[triangleId][(k + 1) % 3]);
			if (boundaryEdges.find(edgeKey) == boundaryEdges.end()){
				unsigned long long oppositeEdgeKey = InvertMeshEdgeKey(edgeKey);
				if (edgeKeyTriangleMap.find(oppositeEdgeKey) != edgeKeyTriangleMap.end()){
					int neighbourId = edgeKeyTriangleMap[oppositeEdgeKey];
					if (markedTriangles.find(neighbourId) == markedTriangles.end()){
						markedTriangles.insert(neighbourId);
						triangleQueue.push(neighbourId);
					}
				}
			}
		}
	}

	return 1;
}

int VoxelBasedSegmentation(const int gridRes, const SimpleMesh & mesh, std::unordered_set<unsigned long> & addedVoxels, std::unordered_set<int> & markedTriangles, bool markExterior){
	for (int i = 0; i < mesh.triangles.size(); i++){
		Point3D<double> baricenter = Point3D<double>(mesh.vertices[mesh.triangles[i][0]] + mesh.vertices[mesh.triangles[i][1]] + mesh.vertices[mesh.triangles[i][2]])*double(gridRes) / 3.0;
		unsigned long axial_index[3] = { floor(baricenter[0]), floor(baricenter[1]), floor(baricenter[2]) };
		unsigned long  voxelKey = SetVoxelElementKey(axial_index, 0);
		if (addedVoxels.find(voxelKey) == addedVoxels.end()){
			if (markExterior){
				markedTriangles.insert(i);
			}
		}
		else{
			if (!markExterior){
				markedTriangles.insert(i);
			}
		}
	}
	return 1;
}

int SplitFaceAdjacentTriangles(std::unordered_map<unsigned long long, int> & sourceClippedEdgeTriangleMap, std::vector<int> & sourceClippedVertexSequence, std::vector<int> & sourceClippedTriangleSequence, std::vector<int> & sourceClippedOppositeVertexSequence, std::vector<double> & sourceCumLengthVector, std::vector<double> & targetCumLengthVector, SimpleMesh & sourceClippedMesh, bool orientation){
	int targetOffset = 1;
	for (int i = 0; i < sourceClippedTriangleSequence.size(); i++){
		int tIndex = sourceClippedTriangleSequence[i];
		int oppositeVertex = sourceClippedOppositeVertexSequence[i];
		int edgeStartVertex = sourceClippedVertexSequence[i];
		int edgeEndVertex = sourceClippedVertexSequence[i + 1];

		std::vector<double> splittingRatio;
		while (sourceCumLengthVector[i] < targetCumLengthVector[targetOffset] && targetCumLengthVector[targetOffset] < sourceCumLengthVector[i + 1] && targetOffset < targetCumLengthVector.size() - 1){
			double splitRatio = (targetCumLengthVector[targetOffset] - sourceCumLengthVector[i]) / (sourceCumLengthVector[i + 1] - sourceCumLengthVector[i]);
			splittingRatio.push_back(splitRatio);
			targetOffset++;
		}

		if (splittingRatio.size()){//Split triangle
			int oldVertexCount = sourceClippedMesh.vertices.size();
			//Add vertices
			for (int j = 0; j < splittingRatio.size(); j++){
				Point3D<double> splitVertexPos = sourceClippedMesh.vertices[edgeStartVertex] * (1.0 - splittingRatio[j]) + sourceClippedMesh.vertices[edgeEndVertex] * splittingRatio[j];
				sourceClippedMesh.vertices.push_back(splitVertexPos);
			}

			//Remove initial triangle
			for (int k = 0; k < 3; k++){
				unsigned long long edgeKey = SetMeshEdgeKey(sourceClippedMesh.triangles[tIndex][k], sourceClippedMesh.triangles[tIndex][(k + 1) % 3]);
				sourceClippedEdgeTriangleMap.erase(edgeKey);
			}

			//Add triangles
			sourceClippedMesh.triangles[tIndex] = TriangleIndex(oppositeVertex, edgeStartVertex, oldVertexCount);
			for (int k = 0; k < 3; k++){
				unsigned long long edgeKey = SetMeshEdgeKey(sourceClippedMesh.triangles[tIndex][k], sourceClippedMesh.triangles[tIndex][(k + 1) % 3]);
				sourceClippedEdgeTriangleMap[edgeKey] = tIndex;
			}

			for (int j = 0; j < splittingRatio.size() - 1; j++){
				sourceClippedMesh.triangles.push_back(TriangleIndex(oppositeVertex, oldVertexCount + j, oldVertexCount + j + 1));
				int lastTriangleIndex = sourceClippedMesh.triangles.size() - 1;
				for (int k = 0; k < 3; k++){
					unsigned long long edgeKey = SetMeshEdgeKey(sourceClippedMesh.triangles[lastTriangleIndex][k], sourceClippedMesh.triangles[lastTriangleIndex][(k + 1) % 3]);
					sourceClippedEdgeTriangleMap[edgeKey] = lastTriangleIndex;
				}
			}
			int lastVertexIndex = sourceClippedMesh.vertices.size() - 1;
			sourceClippedMesh.triangles.push_back(TriangleIndex(oppositeVertex, lastVertexIndex, edgeEndVertex));
			int lastTriangleIndex = sourceClippedMesh.triangles.size() - 1;
			for (int k = 0; k < 3; k++){
				unsigned long long edgeKey = SetMeshEdgeKey(sourceClippedMesh.triangles[lastTriangleIndex][k], sourceClippedMesh.triangles[lastTriangleIndex][(k + 1) % 3]);
				sourceClippedEdgeTriangleMap[edgeKey] = lastTriangleIndex;
			}
		}
	}
}

int SubdivideFaceAdjacentTriangles(std::unordered_map<unsigned long, FaceMap> & sourceFaceMaps, std::unordered_map<unsigned long long, unsigned long long> sourceAdjacentVertexMap[3], std::unordered_map<unsigned long long, int> & sourceClippedVertexIndex, std::unordered_map<unsigned long long, int> & sourceClippedEdgeTriangleMap, SimpleMesh & sourceClippedMesh,
	std::unordered_map<unsigned long, FaceMap> & targetFaceMaps, std::unordered_map<unsigned long long, unsigned long long> targetAdjacentVertexMap[3], std::unordered_map<unsigned long long, int> & targetClippedVertexIndex, std::unordered_map<unsigned long long, int> & targetClippedEdgeTriangleMap, SimpleMesh & targetClippedMesh,
	std::unordered_map<int, int> & matchedVertices, unsigned long faceKey, bool orientation){

	std::vector<int> sourceClippedVertexSequence;
	std::vector<int> sourceClippedTriangleSequence;
	std::vector<int> sourceClippedOppositeVertexSequence;
	std::vector<double> sourceCumLengthVector;
	if (!GetFaceAdjacentTriangles(sourceFaceMaps, faceKey, orientation, sourceAdjacentVertexMap, sourceClippedVertexIndex, sourceClippedEdgeTriangleMap, sourceClippedVertexSequence, sourceClippedTriangleSequence, sourceClippedOppositeVertexSequence, sourceCumLengthVector, sourceClippedMesh)) return 0;
	int sourceSequenceLenght = sourceClippedVertexSequence.size();

	std::vector<int> targetClippedVertexSequence;
	std::vector<int> targetClippedTriangleSequence;
	std::vector<int> targetClippedOppositeVertexSequence;
	std::vector<double> targetCumLengthVector;
	if (!GetFaceAdjacentTriangles(targetFaceMaps, faceKey, !orientation, targetAdjacentVertexMap, targetClippedVertexIndex, targetClippedEdgeTriangleMap, targetClippedVertexSequence, targetClippedTriangleSequence, targetClippedOppositeVertexSequence, targetCumLengthVector, targetClippedMesh)) return 0;
	int targetSequenceLenght = targetClippedVertexSequence.size();

	//sourceClippedVertexSequence and targetClippedVertexSequence are positive oriented (according to the outward surface normal!)

	if (sourceSequenceLenght < 2 && targetSequenceLenght >= 2 || sourceSequenceLenght < 2 && targetSequenceLenght >= 2){
		printf("Non equivalent faces!\n");
		return 0;
	}

	if (targetSequenceLenght >= 2 || sourceSequenceLenght >= 2){

		std::vector<double> inverseSourceCumLengthVector(sourceSequenceLenght);
		for (int i = 0; i < sourceSequenceLenght; i++){
			inverseSourceCumLengthVector[i] = 1.0 - sourceCumLengthVector[sourceSequenceLenght - 1 - i];
		}

		std::vector<double> inverseTargetCumLengthVector(targetSequenceLenght);
		for (int i = 0; i < targetSequenceLenght; i++){
			inverseTargetCumLengthVector[i] = 1.0 - targetCumLengthVector[targetSequenceLenght - 1 - i];
		}

		int sourceInitialVcount = sourceClippedMesh.vertices.size();
		int targetInitialVcount = targetClippedMesh.vertices.size();

		if (targetSequenceLenght > 2)SplitFaceAdjacentTriangles(sourceClippedEdgeTriangleMap, sourceClippedVertexSequence, sourceClippedTriangleSequence, sourceClippedOppositeVertexSequence, sourceCumLengthVector, inverseTargetCumLengthVector, sourceClippedMesh, orientation);
		if (sourceSequenceLenght > 2)SplitFaceAdjacentTriangles(targetClippedEdgeTriangleMap, targetClippedVertexSequence, targetClippedTriangleSequence, targetClippedOppositeVertexSequence, targetCumLengthVector, inverseSourceCumLengthVector, targetClippedMesh, orientation);

		int sourceFinalVcount = sourceClippedMesh.vertices.size();
		int targetFinalVcount = targetClippedMesh.vertices.size();

		if (sourceFinalVcount != (sourceInitialVcount + targetSequenceLenght - 2)){
			printf("Unexpected vertex count! \n");
			return 0;
		}
		if (targetFinalVcount != (targetInitialVcount + sourceSequenceLenght - 2)){
			printf("Unexpected vertex count! \n");
			return 0;
		}

		matchedVertices[targetClippedVertexSequence[targetSequenceLenght - 1]] = sourceClippedVertexSequence[0];
		matchedVertices[targetClippedVertexSequence[0]] = sourceClippedVertexSequence[sourceSequenceLenght - 1];

		for (int i = 1; i < sourceSequenceLenght - 1; i++){
			matchedVertices[targetInitialVcount + i - 1] = sourceClippedVertexSequence[sourceSequenceLenght - 1 - i];
		}

		for (int i = 1; i < targetSequenceLenght - 1; i++){
			matchedVertices[targetClippedVertexSequence[i]] = sourceFinalVcount - i;
		}
	}

	return 1;
}

int StitchClippedMeshes(SimpleMesh & sourceClippedMesh, const std::unordered_set<int> & sourceTrianglesToRemove, std::vector<int> & sourceClippedToStitchedTriangleMap, const SimpleMesh & targetClippedMesh, const std::unordered_set<int> & targetTrianglesToRemove, std::vector<int> & targetClippedToStitchedTriangleMap, SimpleMesh & mergeMesh, std::unordered_map<int, int> & matchedVertices){
	
	std::vector<Point3D<double>> & mergedMeshVertices = mergeMesh.vertices;
	std::vector<TriangleIndex> & mergedMeshTriangles = mergeMesh.triangles;

	mergedMeshTriangles.reserve(sourceClippedMesh.triangles.size() + targetClippedMesh.triangles.size());
	mergedMeshVertices.reserve(sourceClippedMesh.vertices.size() + targetClippedMesh.vertices.size());

	//Average matched vertices
	for (auto vertexIter = matchedVertices.begin(); vertexIter != matchedVertices.end(); vertexIter++){
		int targetVertexIndex = (*vertexIter).first;
		int sourceVertexIndex = (*vertexIter).second;
		sourceClippedMesh.vertices[sourceVertexIndex] = (sourceClippedMesh.vertices[sourceVertexIndex] + targetClippedMesh.vertices[targetVertexIndex]) / 2.0;
	}

	int lastVertexIndex = 0;
	int lastTriangleIndex = 0;

	sourceClippedToStitchedTriangleMap.resize(sourceClippedMesh.triangles.size(),-1);

	std::vector<int> mergedSourceVertexIndex(sourceClippedMesh.vertices.size(), -1);
	for (int i = 0; i < sourceClippedMesh.triangles.size(); i++){
		if (sourceTrianglesToRemove.find(i) == sourceTrianglesToRemove.end()){
			TriangleIndex mergedIndices;
			for (int k = 0; k < 3; k++){
				int vIndex = sourceClippedMesh.triangles[i][k];
				if (mergedSourceVertexIndex[vIndex] == -1){
					mergedMeshVertices.push_back(sourceClippedMesh.vertices[vIndex]);
					mergedSourceVertexIndex[vIndex] = lastVertexIndex;
					lastVertexIndex++;
				}
				mergedIndices[k] = mergedSourceVertexIndex[vIndex];
			}
			mergedMeshTriangles.push_back(mergedIndices);
			sourceClippedToStitchedTriangleMap[i] = lastTriangleIndex;
			lastTriangleIndex++;
		}
	}

	targetClippedToStitchedTriangleMap.resize(targetClippedMesh.triangles.size(), -1);

	std::vector<int> mergedTargetVertexIndex(targetClippedMesh.vertices.size(), -1);
	for (int i = 0; i < targetClippedMesh.triangles.size(); i++){
		if (targetTrianglesToRemove.find(i) == targetTrianglesToRemove.end()){
			TriangleIndex mergedIndices;
			for (int k = 0; k < 3; k++){
				int vIndex = targetClippedMesh.triangles[i][k];
				if (mergedTargetVertexIndex[vIndex] == -1){
					if (matchedVertices.find(vIndex) != matchedVertices.end()){
						if (mergedSourceVertexIndex[matchedVertices[vIndex]] == -1){
							printf("Matched source vertex not indexed! \n");
							return 0;
						}
						else mergedTargetVertexIndex[vIndex] = mergedSourceVertexIndex[matchedVertices[vIndex]];
					}
					else{
						mergedMeshVertices.push_back(targetClippedMesh.vertices[vIndex]);
						mergedTargetVertexIndex[vIndex] = lastVertexIndex;
						lastVertexIndex++;
					}
				}
				mergedIndices[k] = mergedTargetVertexIndex[vIndex];
			}
			mergedMeshTriangles.push_back(mergedIndices);
			targetClippedToStitchedTriangleMap[i] = lastTriangleIndex;
			lastTriangleIndex++;
		}
	}
	return 1;
}

int Stitching(std::unordered_map<unsigned long, FaceMap> & sourceFaceMaps, std::unordered_map<unsigned long long, unsigned long long> sourceAdjacentVertexMap[3], std::unordered_map<unsigned long long, int> & sourceClippedVertexIndex, std::unordered_map<unsigned long long, int> & sourceClippedEdgeTriangleMap, std::unordered_set<int> & sourceTrianglesToRemove, std::vector<int> & sourceClippedToStitchedTriangleMap, SimpleMesh & sourceClippedMesh,
	std::unordered_map<unsigned long, FaceMap> & targetFaceMaps, std::unordered_map<unsigned long long, unsigned long long> targetAdjacentVertexMap[3], std::unordered_map<unsigned long long, int> & targetClippedVertexIndex, std::unordered_map<unsigned long long, int> & targetClippedEdgeTriangleMap, std::unordered_set<int> & targetTrianglesToRemove, std::vector<int> & targetClippedToStitchedTriangleMap,  SimpleMesh & targetClippedMesh,
			  std::unordered_map<unsigned long, bool> & boundaryFaceOrientation,
			  SimpleMesh & mergedMesh){

	std::unordered_map<int, int> matchedVertices;
	//Subdivide edges in boundary faces
	for (auto faceIter = boundaryFaceOrientation.begin(); faceIter != boundaryFaceOrientation.end(); faceIter++){
		unsigned long faceKey = (*faceIter).first;
		bool orientation = (*faceIter).second;
		if (!SubdivideFaceAdjacentTriangles(sourceFaceMaps, sourceAdjacentVertexMap, sourceClippedVertexIndex, sourceClippedEdgeTriangleMap,sourceClippedMesh,
			targetFaceMaps, targetAdjacentVertexMap, targetClippedVertexIndex, targetClippedEdgeTriangleMap,targetClippedMesh,
			matchedVertices, faceKey, orientation)) return 0;
	}

	//Merge
	if (!StitchClippedMeshes(sourceClippedMesh, sourceTrianglesToRemove, sourceClippedToStitchedTriangleMap, targetClippedMesh, targetTrianglesToRemove, targetClippedToStitchedTriangleMap,  mergedMesh, matchedVertices)) return 0;
	return 1;
}
