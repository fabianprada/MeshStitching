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


#include "SimpleMesh.h"
#include "RayTracer.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <queue>

int AddComponent(std::vector<int> & vertexComponent, int vIndex, int currentComponent, const std::vector<std::vector<int>> & neighbours, const std::vector<double> & vertexWeights, double & componentWeight, int & componentSize){
	componentWeight = 0;
	componentSize = 0;
	vertexComponent[vIndex] = currentComponent;
	std::queue<int> visitingQueue;
	visitingQueue.push(vIndex);

	componentWeight += vertexWeights[vIndex];
	componentSize++;

	while (!visitingQueue.empty()){
		int currentVertex = visitingQueue.front();
		visitingQueue.pop();
		const std::vector<int> & vertexNeighbours = neighbours[currentVertex];
		for (int i = 0; i < vertexNeighbours.size(); i++){
			if (vertexComponent[vertexNeighbours[i]] == -1) {
				vertexComponent[vertexNeighbours[i]] = currentComponent;
				visitingQueue.push(vertexNeighbours[i]);
				componentWeight += vertexWeights[vertexNeighbours[i]];
				componentSize++;
			}
			else if (vertexComponent[vertexNeighbours[i]] == currentComponent){}
			else{
				printf("Unexpected Condition On A Connected Component. Expected %d. Obtained %d\n", currentComponent, vertexComponent[vertexNeighbours[i]]);
				return 0;
			}
		}
	}
	return 1;
}

void ComputeNormalizedVertexWeights(const SimpleMesh & mesh, std::vector<double> & vertexWeights){
	vertexWeights.resize(mesh.vertices.size(), 0.0);
	double cumWeights = 0.0;
	for (int i = 0; i < mesh.triangles.size(); i++){
		float area = Point3D<double>::Length(Point3D<double>::CrossProduct(mesh.vertices[mesh.triangles[i][1]] - mesh.vertices[mesh.triangles[i][0]], mesh.vertices[mesh.triangles[i][2]] - mesh.vertices[mesh.triangles[i][0]])) / 2.0;
		for (int v = 0; v < 3; v++){
			vertexWeights[mesh.triangles[i][v]] += area;
			cumWeights += area;
		}
	}
	for (int i = 0; i < mesh.vertices.size(); i++) vertexWeights[i] /= cumWeights;
}

void RemoveSmallComponents(const SimpleMesh & mesh, std::vector<bool> & markedForRemeshing, const int componentSizeThreshold = 5, const double componentWeightThreshold = 0.0001)
{
	std::vector<double> vertexWeights;
	ComputeNormalizedVertexWeights(mesh, vertexWeights);

	int vCount = mesh.vertices.size();
	int tCount = mesh.triangles.size();


	std::vector<std::vector<int>> neighbours(vCount);
	for (int i = 0; i < tCount; i++){
		int vIndices[3];
		for (int k = 0; k < 3; k++) vIndices[k] = mesh.triangles[i][k];
		for (int k = 0; k < 3; k++){
			if (markedForRemeshing[vIndices[k]] && markedForRemeshing[vIndices[(k + 1) % 3]]){
				neighbours[vIndices[k]].push_back(vIndices[(k + 1) % 3]);
				neighbours[vIndices[(k + 1) % 3]].push_back(vIndices[k]);
			}
		}
	}

	std::vector<int> vertexComponent;
	vertexComponent.resize(vCount);
	for (int v = 0; v < vCount; v++) vertexComponent[v] = -1;
	int currentComponent = -1;

	std::vector<double> componentsWeight;
	std::vector<int> componentsSize;
	for (int v = 0; v < vCount; v++){
		if (vertexComponent[v] == -1){
			currentComponent++;
			double currentComponentWeight;
			int currentComponentSize;
			AddComponent(vertexComponent, v, currentComponent, neighbours, vertexWeights, currentComponentWeight, currentComponentSize);
			componentsWeight.push_back(currentComponentWeight);
			componentsSize.push_back(currentComponentSize);
		}
	}

	int removedVertices = 0;
	int preservedVertices = 0;
	for (int v = 0; v < vCount; v++){
		if (markedForRemeshing[v]){
			int vComponent = vertexComponent[v];
			if (!(componentsWeight[vComponent] > componentWeightThreshold && componentsSize[vComponent] > componentSizeThreshold)){
				markedForRemeshing[v] = false;
				removedVertices++;
			}
			else{
				preservedVertices++;
			}
		}
	}
	if (1) printf("Preserved to remesh %d. Removed to remesh %d \n", preservedVertices, removedVertices);
}

class SamplePoint{
public:
	SamplePoint(){
		tIndex = -1;
		baricentricCoord = Point2D<double>(0, 0);
	}
	int tIndex;
	Point2D<double> baricentricCoord;
};

void BiProjectionMap(const SimpleMesh & sourceMesh, const SimpleMesh & targetMesh, const int vertexIndex, float distCutOff, float angleCutOff, SamplePoint & mapping, bool verbose = false){
	MeshRayIntersection  targetRayQuery;
	targetRayQuery.Init(targetMesh.vertices, targetMesh.triangles);

	MeshRayIntersection  sourceRayQuery;
	sourceRayQuery.Init(sourceMesh.vertices, sourceMesh.triangles);

	mapping.tIndex = -1;

	IntersectionData sourceToTargetIntersection;
	if (targetRayQuery.IntersectLine(sourceMesh.vertices[vertexIndex] + sourceMesh.normals[vertexIndex]* 0, sourceMesh.normals[vertexIndex], distCutOff, sourceToTargetIntersection, targetMesh, sourceMesh.normals[vertexIndex], angleCutOff, verbose)){
		const int tId = sourceToTargetIntersection.tId;
		Point3D<double> intersectionPosition = targetMesh.vertices[targetMesh.triangles[tId][0]] * (1.f - sourceToTargetIntersection.u - sourceToTargetIntersection.v) + targetMesh.vertices[targetMesh.triangles[tId][1]] * sourceToTargetIntersection.u + targetMesh.vertices[targetMesh.triangles[tId][2]] * sourceToTargetIntersection.v;
		if (verbose){
			double intersectionDistance = Point3D<double>::Length(sourceMesh.vertices[vertexIndex] - intersectionPosition);
			printf("First Intersection Distance %f \n", intersectionDistance);
		}
		Point3D<double> intersectionNormal = targetMesh.normals[targetMesh.triangles[tId][0]] * (1.f - sourceToTargetIntersection.u - sourceToTargetIntersection.v) + targetMesh.normals[targetMesh.triangles[tId][1]] * sourceToTargetIntersection.u + targetMesh.normals[targetMesh.triangles[tId][2]] * sourceToTargetIntersection.v;
		intersectionNormal /= Point3D<double>::Length(intersectionNormal);

		IntersectionData targetToSourceIntersection;
		if (sourceRayQuery.IntersectLine(intersectionPosition + intersectionNormal* 0, intersectionNormal, distCutOff, targetToSourceIntersection, sourceMesh, intersectionNormal, angleCutOff, verbose)){
			mapping.tIndex = targetToSourceIntersection.tId;
			mapping.baricentricCoord = Point2D<double>(targetToSourceIntersection.u, targetToSourceIntersection.v);
			if (verbose){
				const int tId = targetToSourceIntersection.tId;
				Point3D<double> intersectionPosition = sourceMesh.vertices[sourceMesh.triangles[tId][0]] * (1.f - targetToSourceIntersection.u - targetToSourceIntersection.v) + sourceMesh.vertices[sourceMesh.triangles[tId][1]] * targetToSourceIntersection.u + sourceMesh.vertices[sourceMesh.triangles[tId][2]] * targetToSourceIntersection.v;
				double intersectionDistance = Point3D<double>::Length(sourceMesh.vertices[vertexIndex] - intersectionPosition);
				printf("Second Intersection Distance %f \n", intersectionDistance);
			}
		}
	}
}

void ProjectionMap(const SimpleMesh & sourceMesh, const SimpleMesh & targetMesh, const int vertexIndex, float distCutOff, float angleCutOff, SamplePoint & mapping, bool verbose = false){
	MeshRayIntersection  targetRayQuery;
	targetRayQuery.Init(targetMesh.vertices, targetMesh.triangles);

	mapping.tIndex = -1;

	IntersectionData sourceToTargetIntersection;
	if (targetRayQuery.IntersectLine(sourceMesh.vertices[vertexIndex] + sourceMesh.normals[vertexIndex] * 0, sourceMesh.normals[vertexIndex], distCutOff, sourceToTargetIntersection, targetMesh, sourceMesh.normals[vertexIndex], angleCutOff, verbose)){
		mapping.tIndex = sourceToTargetIntersection.tId;
		mapping.baricentricCoord = Point2D<double>(sourceToTargetIntersection.u, sourceToTargetIntersection.v);
	}
}

void ProjectionMap(const SimpleMesh & _sourceMesh, const SimpleMesh & _targetMesh, float distCutOff, float angleCutOff, std::vector<SamplePoint> & mapping, bool verbose = false){
	
	SimpleMesh sourceMesh;
	sourceMesh.vertices = _sourceMesh.vertices;
	sourceMesh.triangles = _sourceMesh.triangles;
	Point3D< double > center;
	double area = 0.f;
	for (int i = 0; i<sourceMesh.triangles.size(); i++)
	{
		Point3D< double > n = Point3D< double >::CrossProduct(sourceMesh.vertices[sourceMesh.triangles[i][1]] - sourceMesh.vertices[sourceMesh.triangles[i][0]], sourceMesh.vertices[sourceMesh.triangles[i][2]] - sourceMesh.vertices[sourceMesh.triangles[i][0]]);
		Point3D< double > c = (sourceMesh.vertices[sourceMesh.triangles[i][0]] + sourceMesh.vertices[sourceMesh.triangles[i][1]] + sourceMesh.vertices[sourceMesh.triangles[i][2]]) / 3.0;
		double a = (double)Length(n);
		center += c*a, area += a;
	}
	center /= area;
	double max = 0;
	for (int i = 0; i< sourceMesh.vertices.size(); i++) max = std::max< double >(max, Point3D< double >::Length(sourceMesh.vertices[i] - center));
	for (int i = 0; i < sourceMesh.vertices.size(); i++) sourceMesh.vertices[i] = (sourceMesh.vertices[i] - center) / max;
	UpdateNormals(sourceMesh);

	SimpleMesh targetMesh;
	targetMesh.vertices = _targetMesh.vertices;
	targetMesh.triangles = _targetMesh.triangles;
	for (int i = 0; i < targetMesh.vertices.size(); i++) targetMesh.vertices[i] = (targetMesh.vertices[i] - center) / max;
	UpdateNormals(targetMesh);

	MeshRayIntersection  targetRayQuery;
	targetRayQuery.Init(targetMesh.vertices, targetMesh.triangles);

	mapping.resize(sourceMesh.vertices.size());
	for (int i = 0; i < sourceMesh.vertices.size(); i++){
		IntersectionData sourceToTargetIntersection;
		if (targetRayQuery.IntersectLine(sourceMesh.vertices[i] + sourceMesh.normals[i] * 0, sourceMesh.normals[i], distCutOff/max, sourceToTargetIntersection, targetMesh, sourceMesh.normals[i], angleCutOff, verbose)){
			mapping[i].tIndex = sourceToTargetIntersection.tId;
			mapping[i].baricentricCoord = Point2D<double>(sourceToTargetIntersection.u, sourceToTargetIntersection.v);
		}
	}
}

void BiProjectionMap(const SimpleMesh & _sourceMesh, const SimpleMesh & _targetMesh, float distCutOff, float angleCutOff, std::vector<SamplePoint> & mapping){
	
	SimpleMesh sourceMesh;
	sourceMesh.vertices = _sourceMesh.vertices;
	sourceMesh.triangles = _sourceMesh.triangles;
	sourceMesh.normals = _sourceMesh.normals;
	Point3D< double > center;
	double area = 0.f;
	for (int i = 0; i<sourceMesh.triangles.size(); i++)
	{
		Point3D< double > n = Point3D< double >::CrossProduct(sourceMesh.vertices[sourceMesh.triangles[i][1]] - sourceMesh.vertices[sourceMesh.triangles[i][0]], sourceMesh.vertices[sourceMesh.triangles[i][2]] - sourceMesh.vertices[sourceMesh.triangles[i][0]]);
		Point3D< double > c = (sourceMesh.vertices[sourceMesh.triangles[i][0]] + sourceMesh.vertices[sourceMesh.triangles[i][1]] + sourceMesh.vertices[sourceMesh.triangles[i][2]]) / 3.0;
		double a = (double)Length(n);
		center += c*a, area += a;
	}
	center /= area;
	double max = 0;
	for (int i = 0; i< sourceMesh.vertices.size(); i++) max = std::max< double >(max, Point3D< double >::Length(sourceMesh.vertices[i] - center));
	for (int i = 0; i < sourceMesh.vertices.size(); i++) sourceMesh.vertices[i] = (sourceMesh.vertices[i] - center) / max;

	SimpleMesh targetMesh;
	targetMesh.vertices = _targetMesh.vertices;
	targetMesh.triangles = _targetMesh.triangles;
	targetMesh.normals = _targetMesh.normals;
	for (int i = 0; i < targetMesh.vertices.size(); i++) targetMesh.vertices[i] = (targetMesh.vertices[i] - center) / max;
	

	MeshRayIntersection  targetRayQuery;
	targetRayQuery.Init(targetMesh.vertices, targetMesh.triangles);

	MeshRayIntersection  sourceRayQuery;
	sourceRayQuery.Init(sourceMesh.vertices, sourceMesh.triangles);

	mapping.resize(sourceMesh.vertices.size());

	for (int i = 0; i < sourceMesh.vertices.size(); i++){
		IntersectionData sourceToTargetIntersection;
		if (targetRayQuery.IntersectLine(sourceMesh.vertices[i] + sourceMesh.normals[i] * 0, sourceMesh.normals[i], distCutOff/max, sourceToTargetIntersection, targetMesh, sourceMesh.normals[i], angleCutOff)){

			const int tId = sourceToTargetIntersection.tId;
			Point3D<double> intersectionPosition = targetMesh.vertices[targetMesh.triangles[tId][0]] * (1.f - sourceToTargetIntersection.u - sourceToTargetIntersection.v) + targetMesh.vertices[targetMesh.triangles[tId][1]] * sourceToTargetIntersection.u + targetMesh.vertices[targetMesh.triangles[tId][2]] * sourceToTargetIntersection.v;
			Point3D<double> intersectionNormal = targetMesh.normals[targetMesh.triangles[tId][0]] * (1.f - sourceToTargetIntersection.u - sourceToTargetIntersection.v) + targetMesh.normals[targetMesh.triangles[tId][1]] * sourceToTargetIntersection.u + targetMesh.normals[targetMesh.triangles[tId][2]] * sourceToTargetIntersection.v;
			intersectionNormal /= Point3D<double>::Length(intersectionNormal);

			IntersectionData targetToSourceIntersection;
			if (sourceRayQuery.IntersectLine(intersectionPosition + intersectionNormal* 0, intersectionNormal, distCutOff/max, targetToSourceIntersection, sourceMesh, intersectionNormal, angleCutOff)){
				mapping[i].tIndex = targetToSourceIntersection.tId;
				mapping[i].baricentricCoord = Point2D<double>(targetToSourceIntersection.u, targetToSourceIntersection.v);
			}
		}
	}
}

void BiProjectionDistance(const SimpleMesh & sourceMesh, const SimpleMesh & targetMesh, const  std::vector<Eigen::VectorXf> & coordinates, float distCutOff, float angleCutOff, std::vector<float> & distance){
	std::vector<SamplePoint> mapping;
	BiProjectionMap(sourceMesh, targetMesh, distCutOff, angleCutOff, mapping);
	distance.resize(sourceMesh.vertices.size(), FLT_MAX);
	for (int i = 0; i < sourceMesh.vertices.size(); i++){
		int tIndex = mapping[i].tIndex;		
		if (tIndex != -1){
			Point2D<double> baricentric = mapping[i].baricentricCoord;
			Eigen::VectorXf mappedCoordinate = coordinates[sourceMesh.triangles[tIndex][0]] * (1.f - baricentric[0] - baricentric[1]) + coordinates[sourceMesh.triangles[tIndex][1]] * baricentric[0] + coordinates[sourceMesh.triangles[tIndex][2]] * baricentric[1];
			float infDistance = (coordinates[i] - mappedCoordinate).lpNorm<Eigen::Infinity>();
			distance[i] = infDistance;
		}
	}
}

void BiProjectionDistance(const SimpleMesh & sourceMesh, const SimpleMesh & targetMesh, float distCutOff, float angleCutOff, std::vector<float> & distance){
	std::vector<SamplePoint> mapping;
	BiProjectionMap(sourceMesh, targetMesh, distCutOff, angleCutOff, mapping);
	distance.resize(sourceMesh.vertices.size(), FLT_MAX);
	for (int i = 0; i < sourceMesh.vertices.size(); i++){
		int tIndex = mapping[i].tIndex;
		if (tIndex != -1){
			Point2D<double> baricentric = mapping[i].baricentricCoord;
			Point3D<double> projection = sourceMesh.vertices[sourceMesh.triangles[tIndex][0]] * (1.f - baricentric[0] - baricentric[1]) + sourceMesh.vertices[sourceMesh.triangles[tIndex][1]] * baricentric[0] + sourceMesh.vertices[sourceMesh.triangles[tIndex][2]] * baricentric[1];
			distance[i] = Point3D<double>::Length(projection - sourceMesh.vertices[i]);
		}
	}
}

void ProjectionDistance(const SimpleMesh & sourceMesh, const SimpleMesh & targetMesh, const  std::vector<Eigen::VectorXf> & coordinates, float distCutOff, float angleCutOff, std::vector<float> & distance){
	std::vector<SamplePoint> bimapping;
	BiProjectionMap(sourceMesh, targetMesh, distCutOff, angleCutOff, bimapping);

	std::vector<SamplePoint> mapping;
	ProjectionMap(sourceMesh, targetMesh, distCutOff, angleCutOff, mapping);

	distance.resize(sourceMesh.vertices.size(), -1.f);
	for (int i = 0; i < sourceMesh.vertices.size(); i++){
		int tIndex = bimapping[i].tIndex;
		if (tIndex != -1){
			tIndex = mapping[i].tIndex;
			Point2D<double> baricentric = mapping[i].baricentricCoord;
			Point3D<double> projection = targetMesh.vertices[targetMesh.triangles[tIndex][0]] * (1.f - baricentric[0] - baricentric[1]) + targetMesh.vertices[targetMesh.triangles[tIndex][1]] * baricentric[0] + targetMesh.vertices[targetMesh.triangles[tIndex][2]] * baricentric[1];
			distance[i] = Point3D<double>::Length(projection - sourceMesh.vertices[i]);
		}
	}
}

void ProjectionDistance(const SimpleMesh & sourceMesh, const SimpleMesh & targetMesh, float distCutOff, float angleCutOff, std::vector<float> & distance){
	std::vector<SamplePoint> mapping;
	ProjectionMap(sourceMesh, targetMesh, distCutOff, angleCutOff, mapping);

	distance.resize(sourceMesh.vertices.size(), -1.f);
	for (int i = 0; i < sourceMesh.vertices.size(); i++){
		int tIndex = mapping[i].tIndex;
		if (tIndex != -1){
			Point2D<double> baricentric = mapping[i].baricentricCoord;
			Point3D<double> projection = targetMesh.vertices[targetMesh.triangles[tIndex][0]] * (1.f - baricentric[0] - baricentric[1]) + targetMesh.vertices[targetMesh.triangles[tIndex][1]] * baricentric[0] + targetMesh.vertices[targetMesh.triangles[tIndex][2]] * baricentric[1];
			distance[i] = Point3D<double>::Length(projection - sourceMesh.vertices[i]);
		}
		else{
			distance[i] = FLT_MAX;
		}
	}
}


void ThresholdedProjectionMap(const SimpleMesh & sourceMesh, const SimpleMesh & targetMesh, const std::vector<Eigen::VectorXf> & coordinates, float distCutOff, float angleCutOff, std::vector<SamplePoint> & mapping, float biProjectionDistanceThreshold){
	std::vector<SamplePoint> bimapping;
	BiProjectionMap(sourceMesh, targetMesh, distCutOff, angleCutOff, bimapping);
	ProjectionMap(sourceMesh, targetMesh, distCutOff, angleCutOff, mapping);
	for (int i = 0; i < sourceMesh.vertices.size(); i++) {
		if (bimapping[i].tIndex == -1) mapping[i].tIndex = -1;
		else{
			int tIndex = bimapping[i].tIndex;
			Point2D<double> baricentric = bimapping[i].baricentricCoord;
			Eigen::VectorXf mappedCoordinate = coordinates[sourceMesh.triangles[tIndex][0]] * (1.f - baricentric[0] - baricentric[1]) + coordinates[sourceMesh.triangles[tIndex][1]] * baricentric[0] + coordinates[sourceMesh.triangles[tIndex][2]] * baricentric[1];
			float infDistance = (coordinates[i] - mappedCoordinate).lpNorm<Eigen::Infinity>();
			if (infDistance > biProjectionDistanceThreshold) mapping[i].tIndex = -1;
		}
	}
}

void ThresholdedProjectionMap(const SimpleMesh & sourceMesh, const SimpleMesh & targetMesh, float distCutOff, float angleCutOff, std::vector<SamplePoint> & mapping, float biProjectionDistanceThreshold){
	std::vector<SamplePoint> bimapping;
	BiProjectionMap(sourceMesh, targetMesh, distCutOff, angleCutOff, bimapping);
	ProjectionMap(sourceMesh, targetMesh, distCutOff, angleCutOff, mapping);
	int validPoints = 0;
	for (int i = 0; i < sourceMesh.vertices.size(); i++) {
		if (bimapping[i].tIndex == -1) mapping[i].tIndex = -1;
		else{
			int tIndex = bimapping[i].tIndex;
			Point2D<double> baricentric = bimapping[i].baricentricCoord;
			Point3D<double> mappedPos = sourceMesh.vertices[sourceMesh.triangles[tIndex][0]] * (1.f - baricentric[0] - baricentric[1]) + sourceMesh.vertices[sourceMesh.triangles[tIndex][1]] * baricentric[0] + sourceMesh.vertices[sourceMesh.triangles[tIndex][2]] * baricentric[1];
			double distance = Point3D<double>::Length(sourceMesh.vertices[i] - mappedPos);
			if (distance > biProjectionDistanceThreshold) mapping[i].tIndex = -1;
			else validPoints++;
		}
	}
	if(0) printf("Valid points %d of %d \n", validPoints, sourceMesh.vertices.size());
}

void AmbientVisibilityScore(const SimpleMesh & sourceMesh, float distCutOff, std::vector<float> & visibilityScore, double coneAngle, const int numSamples){
	visibilityScore.resize(sourceMesh.vertices.size(), 0);

	std::vector<Point3D< double >>normals(sourceMesh.vertices.size(), Point3D<double>(0.0, 0.0, 0.0));
	for (int t = 0; t < sourceMesh.triangles.size(); t++){
		Point3D<double> d01 = sourceMesh.vertices[sourceMesh.triangles[t][1]] - sourceMesh.vertices[sourceMesh.triangles[t][0]];
		Point3D<double> d02 = sourceMesh.vertices[sourceMesh.triangles[t][2]] - sourceMesh.vertices[sourceMesh.triangles[t][0]];
		Point3D<double> n = Point3D<double>::CrossProduct(d01, d02);
		for (int v = 0; v < 3; v++) normals[sourceMesh.triangles[t][v]] += n;
	}

	for (int i = 0; i < normals.size(); i++){
		if (Point3D<double>::Length(normals[i])> 0) normals[i] /= Point3D<double>::Length(normals[i]);
	}


	MeshRayIntersection  sourceRayQuery;
	sourceRayQuery.Init(sourceMesh.vertices, sourceMesh.triangles);

	for (int i = 0; i < sourceMesh.vertices.size(); i++){
		//Point3D< double > n = sourceMesh.normals[i];
		Point3D< double > n = normals[i];
		Point3D< double > p = sourceMesh.vertices[i];
		Point3D<double> intersectionPosition;
		Point3D<double> intersectionNormal;
		//Check if is interior
		
		bool isExterior = true;
		if (sourceRayQuery.IntersectRay(p + n * 0.1f, n, distCutOff, intersectionPosition, intersectionNormal, sourceMesh, false, false)){
			if (Point3D<double>::Dot(intersectionNormal, n) > 0){
				isExterior = false;
				visibilityScore[i] = 0.0;
			}
		}

		if (isExterior){//////////Compute ambient oclusion
			//Generate a random ortogonal vector
			Point3D<double> randVector(double(rand()) / double(RAND_MAX), double(rand()) / double(RAND_MAX), double(rand()) / double(RAND_MAX));
			randVector /= Point3D<double>::Length(randVector);
			
			Point3D<double> orth_vector1 = Point3D<double>::CrossProduct(n, randVector);
			orth_vector1 /= Point3D<double>::Length(orth_vector1);

			Point3D<double> orth_vector2 = Point3D<double>::CrossProduct(n, orth_vector1);
			orth_vector2 /= Point3D<double>::Length(orth_vector2);
			int intersectionCount = 0;

			for (int j = 0; j < numSamples; j++){
				double angle = 2.0*M_PI*double(j) / numSamples;
				Point3D<double> planeDirection = orth_vector1 * cos(angle) + orth_vector2 * sin(angle);
				Point3D<double> rayDirection = n * cos(coneAngle) + planeDirection * sin(coneAngle);
				if (sourceRayQuery.IntersectRay(p + n * 0.1f, rayDirection, distCutOff, intersectionPosition, intersectionNormal, sourceMesh, false, false)){
					intersectionCount++;
					//printf("Angle %f! \n", angle); 
				}
			}
			visibilityScore[i] = double(numSamples - intersectionCount) / double(numSamples);
			//printf("%f \n", visibilityScore[i]);
		}
	}
}

void VisibilityMap(const SimpleMesh & sourceMesh, float distCutOff, std::vector<float> & visibilityScore, bool markBackVertices, int & markedBackVertices){

	visibilityScore.resize(sourceMesh.vertices.size(), distCutOff);

	MeshRayIntersection  sourceRayQuery;
	sourceRayQuery.Init(sourceMesh.vertices, sourceMesh.triangles);
	
	for (int i = 0; i < sourceMesh.triangles.size(); i++){
		Point3D< double > n = Point3D< double >::CrossProduct(sourceMesh.vertices[sourceMesh.triangles[i][1]] - sourceMesh.vertices[sourceMesh.triangles[i][0]], sourceMesh.vertices[sourceMesh.triangles[i][2]] - sourceMesh.vertices[sourceMesh.triangles[i][0]]);
		n /= Point3D<double>::Length(n);
		Point3D< double > c = (sourceMesh.vertices[sourceMesh.triangles[i][0]] + sourceMesh.vertices[sourceMesh.triangles[i][1]] + sourceMesh.vertices[sourceMesh.triangles[i][2]]) / 3.0;
		Point3D<double> intersectionPosition;
		Point3D<double> intersectionNormal;
		double score = FLT_MAX;
		if (sourceRayQuery.IntersectRay(c + n * 0.01f, n, distCutOff, intersectionPosition, intersectionNormal, sourceMesh, false,false)){
			if (Point3D<double>::Dot(intersectionNormal, n) < 0){
				score = Point3D<double>::Length(c - intersectionPosition);
			}
			else{
				score = 0.0;
			}
		}
		for (int k = 0; k < 3; k++){
			visibilityScore[sourceMesh.triangles[i][k]] = std::min<float>(score, visibilityScore[sourceMesh.triangles[i][k]]);
		}
	}

	for (int i = 0; i < sourceMesh.vertices.size(); i++) visibilityScore[i] /= distCutOff;

	markedBackVertices = 0;

	if (markBackVertices){
	
		std::vector<bool> allAdjacentFacesAreBack(sourceMesh.vertices.size(), true);
		for (int i = 0; i < sourceMesh.triangles.size(); i++){
			Point3D< double > n = Point3D< double >::CrossProduct(sourceMesh.vertices[sourceMesh.triangles[i][1]] - sourceMesh.vertices[sourceMesh.triangles[i][0]], sourceMesh.vertices[sourceMesh.triangles[i][2]] - sourceMesh.vertices[sourceMesh.triangles[i][0]]);
			n /= Point3D<double>::Length(n);
			Point3D< double > c = (sourceMesh.vertices[sourceMesh.triangles[i][0]] + sourceMesh.vertices[sourceMesh.triangles[i][1]] + sourceMesh.vertices[sourceMesh.triangles[i][2]]) / 3.0;
			Point3D<double> intersectionPosition;
			Point3D<double> intersectionNormal;
			double score = FLT_MAX;
			if (sourceRayQuery.IntersectRay(c - n * 0.01f, -n, distCutOff, intersectionPosition, intersectionNormal, sourceMesh, false, false)){
				for (int k = 0; k < 3; k++){
					allAdjacentFacesAreBack[sourceMesh.triangles[i][k]] = false;
				}
			}
		}

		ColoredMesh coloredBack;
		coloredBack.vertices = sourceMesh.vertices;
		coloredBack.triangles = sourceMesh.triangles;
		coloredBack.colors.resize(sourceMesh.vertices.size());

		RemoveSmallComponents(sourceMesh, allAdjacentFacesAreBack);

		int backVerticesCount = 0;
		for (int i = 0; i < allAdjacentFacesAreBack.size(); i++){
			if (allAdjacentFacesAreBack[i]){
				backVerticesCount++;
				visibilityScore[i] = FLT_MAX;
				coloredBack.colors[i] = Point3D<double>(0,0,255);
			}
			else{
				coloredBack.colors[i] = Point3D<double>(255, 0, 0);
			}
		}
		if (1) WriteColoredMesh(coloredBack, "Mesh_BackVertices.ply");
		printf("Back Vertices %d \n", backVerticesCount);
		markedBackVertices = backVerticesCount;
	}

}

int ComputeBackVertices(const SimpleMesh & sourceMesh, float distCutOff, std::vector<float> & backVertexScore){

	backVertexScore.resize(sourceMesh.vertices.size(), 0);

	MeshRayIntersection  sourceRayQuery;
	sourceRayQuery.Init(sourceMesh.vertices, sourceMesh.triangles);

	std::vector<bool> allAdjacentFacesAreBack(sourceMesh.vertices.size(), true);
	for (int i = 0; i < sourceMesh.triangles.size(); i++){
		Point3D< double > n = Point3D< double >::CrossProduct(sourceMesh.vertices[sourceMesh.triangles[i][1]] - sourceMesh.vertices[sourceMesh.triangles[i][0]], sourceMesh.vertices[sourceMesh.triangles[i][2]] - sourceMesh.vertices[sourceMesh.triangles[i][0]]);
		n /= Point3D<double>::Length(n);
		Point3D< double > c = (sourceMesh.vertices[sourceMesh.triangles[i][0]] + sourceMesh.vertices[sourceMesh.triangles[i][1]] + sourceMesh.vertices[sourceMesh.triangles[i][2]]) / 3.0;
		Point3D<double> intersectionPosition;
		Point3D<double> intersectionNormal;
		double score = FLT_MAX;
		if (sourceRayQuery.IntersectRay(c - n * 0.01f, -n, distCutOff, intersectionPosition, intersectionNormal, sourceMesh, false, false)){
			for (int k = 0; k < 3; k++){
				allAdjacentFacesAreBack[sourceMesh.triangles[i][k]] = false;
			}
		}
	}

	ColoredMesh coloredBack;
	coloredBack.vertices = sourceMesh.vertices;
	coloredBack.triangles = sourceMesh.triangles;
	coloredBack.colors.resize(sourceMesh.vertices.size());

	RemoveSmallComponents(sourceMesh, allAdjacentFacesAreBack);

	int backVerticesCount = 0;
	for (int i = 0; i < allAdjacentFacesAreBack.size(); i++){
		if (allAdjacentFacesAreBack[i]){
			backVerticesCount++;
			backVertexScore[i] = FLT_MAX;
			coloredBack.colors[i] = Point3D<double>(0, 0, 255);
		}
		else{
			coloredBack.colors[i] = Point3D<double>(255, 0, 0);
		}
	}
	if (1) WriteColoredMesh(coloredBack, "Mesh_BackVertices.ply");
	printf("Back Vertices %d \n", backVerticesCount);
	return backVerticesCount;
}