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


#ifndef SURFACE_ARAP_INCLUDED
#define SURFACE_ARAP_INCLUDED

#include "SimpleMesh.h"
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/StdVector>
#include <Eigen/Sparse>
#include <Eigen/SVD>


template<typename InClass, typename OutClass>
void CopyPointClasses(const std::vector<InClass> & in, std::vector<OutClass> & out) {
	out.clear();
	out.resize(in.size());
	for (int i = 0; i < in.size(); i++)out[i] = OutClass(in[i][0], in[i][1], in[i][2]);
}

namespace SurfaceARAP{

	Eigen::Matrix2d TriMetricTensor(const Eigen::Vector3d v0, const  Eigen::Vector3d v1, const  Eigen::Vector3d v2){
		//Eigen::Vector3d d[2] = { v1 - v0, v2 - v0};
		Eigen::Vector3d d[2];
		d[0] = v1 - v0;
		d[1] = v2 - v0;
		Eigen::Matrix2d g;
		for (int i = 0; i < 2; i++)for (int j = 0; j < 2; j++){
			g(i, j) = d[i].dot(d[j]);
		}
		return g;
	}

	Eigen::Matrix3d TriStiffnesMatrix(const Eigen::Matrix2d & g){//This the negative stiffness matrix
		Eigen::Vector2d d[3];
		d[0] = Eigen::Vector2d(-1.0, -1.0);
		d[1] = Eigen::Vector2d(1.0, 0.0);
		d[2] = Eigen::Vector2d(0.0, 1.0);
		Eigen::Matrix3d s;
		Eigen::Matrix2d  g_inv = g.inverse();
		double triArea = sqrt(g.determinant());
		for (int i = 0; i < 3; i++)for (int j = 0; j < 3; j++) s(i, j) = -(d[i].dot(g_inv*d[j]))*triArea;
		return s;
	}

	void TriangleCotangentWeights(const std::vector<Eigen::Vector3d> & vertices, const std::vector<TriangleIndex> & triangles, std::vector<double> & cotangentWeights){
		cotangentWeights.resize(triangles.size() * 3);
		for (int t = 0; t < triangles.size(); t++){
			Eigen::Matrix2d g = TriMetricTensor(vertices[triangles[t][0]], vertices[triangles[t][1]], vertices[triangles[t][2]]);
			if (g.determinant() <= 0.0){
				printf("Non Positive Area Triangle!. %f\n", g.determinant());
				printf("Mass matrix assigned to 0.000001*Id\n");
				g(0, 0) = g(1, 1) = 0.000001;
				g(0, 1) = g(1, 0) = 0.0;
			}
			Eigen::Matrix3d s = TriStiffnesMatrix(g);
			for (int i = 0; i < 3; i++) cotangentWeights[3 * t + i] = s(i, (i + 1) % 3);
		}
	}

	int Intrinsic_ARAP_Setup(const int n, const std::vector<TriangleIndex> & triangles, const std::vector<double> & cotangentWeights, const std::vector<int> & fixedIndices, const std::vector<int> & softIndices, const std::vector<float> & softWeights, std::vector<int> & varIndex, Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> & ARAPCholesky, Eigen::SparseMatrix<double> & stiffnessMatrix){

		printf("Total Vertices %d \n", n);
		printf("Fixed Vertices %d \n", fixedIndices.size());
		printf("Soft Vertices %d \n", softIndices.size());

		std::vector<bool> isFixed(n, false);
		for (int i = 0; i < fixedIndices.size(); i++)isFixed[fixedIndices[i]] = true;
		for (int i = 0; i < softIndices.size(); i++) if (isFixed[softIndices[i]]){
			printf("Soft variable can not be fixed! \n");
			return 0;
		}

		varIndex.clear();
		varIndex.resize(n);
		for (int i = 0; i < n; i++) varIndex[i] = -1;

		int varCounter = 0;
		for (int i = 0; i < n; i++){
			if (!isFixed[i]){
				varIndex[i] = varCounter;
				varCounter++;
			}
		}

		if (varCounter != (n - fixedIndices.size())){
			printf("Variable counters unexpected! \n");
			return 0;
		}

		std::vector<Eigen::Triplet<double>> stiffnessTriplets;
		stiffnessTriplets.reserve(6 * triangles.size());

		for (int t = 0; t < triangles.size(); t++){
			for (int j = 0; j < 3; j++){
				int vP = varIndex[triangles[t][j]];
				int vN = varIndex[triangles[t][(j + 1) % 3]];
				if (vP != -1){
					stiffnessTriplets.push_back(Eigen::Triplet<double>(vP, vP, cotangentWeights[3 * t + j]));
					if (vN != -1) stiffnessTriplets.push_back(Eigen::Triplet<double>(vP, vN, -cotangentWeights[3 * t + j]));
				}
				if (vN != -1){
					stiffnessTriplets.push_back(Eigen::Triplet<double>(vN, vN, cotangentWeights[3 * t + j]));
					if (vP != -1) stiffnessTriplets.push_back(Eigen::Triplet<double>(vN, vP, -cotangentWeights[3 * t + j]));
				}
			}
		}

		stiffnessMatrix.resize(varCounter, varCounter);
		stiffnessMatrix.setFromTriplets(stiffnessTriplets.begin(), stiffnessTriplets.end());


		std::vector<Eigen::Triplet<double>> softWeighTriplets;
		softWeighTriplets.reserve(softIndices.size());
		for (int i = 0; i < softIndices.size(); i++){
			int vi = varIndex[softIndices[i]];
			softWeighTriplets.push_back(Eigen::Triplet<double>(vi, vi, softWeights[i]));
		}

		Eigen::SparseMatrix<double> softWeightsMatrix;
		softWeightsMatrix.resize(varCounter, varCounter);
		softWeightsMatrix.setFromTriplets(softWeighTriplets.begin(), softWeighTriplets.end());

		ARAPCholesky.analyzePattern(stiffnessMatrix + softWeightsMatrix);
		ARAPCholesky.factorize(stiffnessMatrix + softWeightsMatrix);

		return 1;
	}

	int Intrinsic_ARAP_Solve(const std::vector<TriangleIndex> & triangles, const std::vector<Eigen::Vector3d> & restVertices, const std::vector<double> & cotangentWeights, const std::vector<int> & fixedIndices, const std::vector<int> & softIndices, const std::vector<float> & softWeights, const std::vector<Eigen::Vector3d> & softConstraints, const std::vector<int> & varIndex, const Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> & ARAPCholesky, std::vector<Eigen::Vector3d> & newVertices){
		int tCount = (int)triangles.size();
		std::vector<Eigen::Matrix3d> rotationMatrices(tCount);
		int threads = omp_get_num_procs();
		for (int t = 0; t < tCount; t++) for (int k = 0; k < 3; k++)for (int l = 0; l < 3; l++) rotationMatrices[t](k, l) = 0.0;

#pragma omp parallel for num_threads( threads )
		for (int t = 0; t < tCount; t++){
			Eigen::Matrix3d scatterMat;
			for (int k = 0; k < 3; k++)for (int l = 0; l < 3; l++) scatterMat(k, l) = 0.0;
			double averageEdge = 0.0;
			double averageCotangent = 0.0;
			for (int i = 0; i < 3; i++){
				Eigen::Vector3d restEdge = restVertices[triangles[t][(i + 1) % 3]] - restVertices[triangles[t][i]];
				Eigen::Vector3d newEdge = newVertices[triangles[t][(i + 1) % 3]] - newVertices[triangles[t][i]];
				for (int k = 0; k < 3; k++)for (int l = 0; l < 3; l++) scatterMat(k, l) += restEdge[k] * newEdge[l] * cotangentWeights[3 * t + i];
				averageCotangent += cotangentWeights[3 * t + i];
				averageEdge += restEdge.norm();
			}
			averageEdge /= 3.0;
			averageCotangent /= 3.0;
			Eigen::Vector3d restEdge = (restVertices[triangles[t][1]] - restVertices[triangles[t][0]]).cross(restVertices[triangles[t][2]] - restVertices[triangles[t][0]]);
			restEdge *= (averageEdge / restEdge.norm());
			Eigen::Vector3d newEdge = (newVertices[triangles[t][1]] - newVertices[triangles[t][0]]).cross(newVertices[triangles[t][2]] - newVertices[triangles[t][0]]);
			newEdge *= (averageEdge / newEdge.norm());
			for (int k = 0; k < 3; k++)for (int l = 0; l < 3; l++) scatterMat(k, l) += restEdge[k] * newEdge[l] * averageCotangent;

			Eigen::JacobiSVD<Eigen::Matrix3d> mSVD(scatterMat, Eigen::ComputeFullU | Eigen::ComputeFullV);
			Eigen::Matrix3d U = mSVD.matrixU();
			Eigen::Matrix3d V = mSVD.matrixV();
			if (U.determinant()*V.determinant() < 0.0) for (int k = 0; k < 3; k++) U(k, 2) *= -1.0;
			rotationMatrices[t] = V*U.transpose();
		}

		Eigen::MatrixXd rhs;
		int numFree = (int)restVertices.size() - (int)fixedIndices.size();
		rhs.resize(numFree, 3);
		for (int k = 0; k < numFree; k++) for (int l = 0; l < 3; l++) rhs(k, l) = 0.0;
		for (int t = 0; t < triangles.size(); t++){
			for (int j = 0; j < 3; j++){
				int vP = varIndex[triangles[t][j]];
				int vN = varIndex[triangles[t][(j + 1) % 3]];
				Eigen::Vector3d restEdge = restVertices[triangles[t][j]] - restVertices[triangles[t][(j + 1) % 3]];
				Eigen::Vector3d rotatedEdge = (rotationMatrices[t] * restEdge) * cotangentWeights[3 * t + j];
				if (vP != -1){
					rhs.row(vP) += rotatedEdge;
					if (vN == -1) rhs.row(vP) += newVertices[triangles[t][(j + 1) % 3]] * cotangentWeights[3 * t + j];
				}
				if (vN != -1){
					rhs.row(vN) -= rotatedEdge;
					if (vP == -1) rhs.row(vN) += newVertices[triangles[t][j]] * cotangentWeights[3 * t + j];
				}
			}
		}

#pragma omp parallel for num_threads( threads )
		for (int i = 0; i < softIndices.size(); i++){
			int vi = varIndex[softIndices[i]];
			rhs.row(vi) += softConstraints[i] * softWeights[i];
		}

		Eigen::MatrixXd solution = ARAPCholesky.solve(rhs);
		for (int i = 0; i < restVertices.size(); i++){
			int vi = varIndex[i];
			if (vi != -1) for (int k = 0; k < 3; k++) newVertices[i][k] = solution(vi, k);
		}
	}

	void Intrinsic_ARAP_Energy(const std::vector<TriangleIndex> & triangles, const std::vector<Eigen::Vector3d> & restVertices, const std::vector<Eigen::Vector3d> & newVertices, const std::vector<double> & cotangentWeights, std::vector<double> & energy){
		int tCount = (int)triangles.size();

		std::vector<Eigen::Matrix3d> rotationMatrices(tCount);
		int threads = omp_get_num_procs();
		for (int t = 0; t < tCount; t++) for (int k = 0; k < 3; k++)for (int l = 0; l < 3; l++) rotationMatrices[t](k, l) = 0.0;

#pragma omp parallel for num_threads( threads )
		for (int t = 0; t < tCount; t++){
			Eigen::Matrix3d scatterMat;
			for (int k = 0; k < 3; k++)for (int l = 0; l < 3; l++) scatterMat(k, l) = 0.0;
			for (int i = 0; i < 3; i++){
				Eigen::Vector3d restEdge = restVertices[triangles[t][(i + 1) % 3]] - restVertices[triangles[t][i]];
				Eigen::Vector3d newEdge = newVertices[triangles[t][(i + 1) % 3]] - newVertices[triangles[t][i]];
				for (int k = 0; k < 3; k++)for (int l = 0; l < 3; l++) scatterMat(k, l) += restEdge[k] * newEdge[l] * cotangentWeights[3 * t + i];
			}
			Eigen::JacobiSVD<Eigen::Matrix3d> mSVD(scatterMat, Eigen::ComputeFullU | Eigen::ComputeFullV);
			Eigen::Matrix3d U = mSVD.matrixU();
			Eigen::Matrix3d V = mSVD.matrixV();
			if (U.determinant()*V.determinant() < 0.0) for (int k = 0; k < 3; k++) U(k, 2) *= -1.0;
			rotationMatrices[t] = V*U.transpose();
		}

		energy.clear();
		energy.resize(tCount, 0.0);

		for (int t = 0; t < triangles.size(); t++){
			for (int j = 0; j < 3; j++){
				Eigen::Vector3d newEdge = newVertices[triangles[t][j]] - newVertices[triangles[t][(j + 1) % 3]];
				Eigen::Vector3d restEdge = restVertices[triangles[t][j]] - restVertices[triangles[t][(j + 1) % 3]];
				Eigen::Vector3d rotatedEdge = (rotationMatrices[t] * restEdge);
				energy[t] += (newEdge - rotatedEdge).squaredNorm() * cotangentWeights[3 * t + j];
			}
		}

		if (1){ //Debug
			for (int t = 0; t < tCount; t++){
				if (energy[t] < 0.0){
					printf("Negative Triangle Energy!!\n");
				}
			}
		}
	}


	int Tri_ARAP_Setup(const int n, const std::vector<TriangleIndex> & triangles, const std::vector<double> & cotangentWeights, const std::vector<int> & fixedIndices, const std::vector<int> & softIndices, const std::vector<float> & softWeights, std::vector<int> & varIndex, Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> & ARAPCholesky, Eigen::SparseMatrix<double> & stiffnessMatrix){

		printf("Total Vertices %d \n", n);
		printf("Fixed Vertices %d \n", fixedIndices.size());
		printf("Soft Vertices %d \n", softIndices.size());

		std::vector<bool> isFixed(n, false);
		for (int i = 0; i < fixedIndices.size(); i++)isFixed[fixedIndices[i]] = true;
		for (int i = 0; i < softIndices.size(); i++) if (isFixed[softIndices[i]]){
			printf("Soft variable can not be fixed! \n");
			return 0;
		}

		varIndex.clear();
		varIndex.resize(n);
		for (int i = 0; i < n; i++) varIndex[i] = -1;

		int varCounter = 0;
		for (int i = 0; i < n; i++){
			if (!isFixed[i]){
				varIndex[i] = varCounter;
				varCounter++;
			}
		}

		if (varCounter != (n - fixedIndices.size())){
			printf("Variable counters unexpected! \n");
			return 0;
		}

		std::vector<Eigen::Triplet<double>> stiffnessTriplets;
		stiffnessTriplets.reserve(18 * triangles.size());

		for (int t = 0; t < triangles.size(); t++){
			for (int i = 0; i < 3; i++){
				for (int j = 0; j < 3; j++){
					int vP = varIndex[triangles[t][j]];
					int vN = varIndex[triangles[t][(j + 1) % 3]];
					if (vP != -1){
						stiffnessTriplets.push_back(Eigen::Triplet<double>(vP, vP, cotangentWeights[3 * t + j]));
						if (vN != -1) stiffnessTriplets.push_back(Eigen::Triplet<double>(vP, vN, -cotangentWeights[3 * t + j]));
					}
					if (vN != -1){
						stiffnessTriplets.push_back(Eigen::Triplet<double>(vN, vN, cotangentWeights[3 * t + j]));
						if (vP != -1) stiffnessTriplets.push_back(Eigen::Triplet<double>(vN, vP, -cotangentWeights[3 * t + j]));
					}
				}
			}
		}


		stiffnessMatrix.resize(varCounter, varCounter);
		stiffnessMatrix.setFromTriplets(stiffnessTriplets.begin(), stiffnessTriplets.end());


		std::vector<Eigen::Triplet<double>> softWeighTriplets;
		softWeighTriplets.reserve(softIndices.size());
		for (int i = 0; i < softIndices.size(); i++){
			int vi = varIndex[softIndices[i]];
			softWeighTriplets.push_back(Eigen::Triplet<double>(vi, vi, softWeights[i]));
		}

		Eigen::SparseMatrix<double> softWeightsMatrix;
		softWeightsMatrix.resize(varCounter, varCounter);
		softWeightsMatrix.setFromTriplets(softWeighTriplets.begin(), softWeighTriplets.end());

		ARAPCholesky.analyzePattern(stiffnessMatrix + softWeightsMatrix);
		ARAPCholesky.factorize(stiffnessMatrix + softWeightsMatrix);

		return 1;
	}

	int Tri_ARAP_Setup(const int n, const std::vector<TriangleIndex> & triangles, const std::vector<double> & cotangentWeightsSource, const std::vector<double> & cotangentWeightsTarget, const float sourceWeight, const std::vector<int> & fixedIndices, const std::vector<int> & softIndices, const std::vector<float> & softWeights, std::vector<int> & varIndex, Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> & ARAPCholesky, Eigen::SparseMatrix<double> & stiffnessMatrix){
		if (0) printf("Total Vertices %d \n", n);
		if (0) printf("Fixed Vertices %d \n", fixedIndices.size());
		if (0) printf("Soft Vertices %d \n", softIndices.size());

		std::vector<bool> isFixed(n, false);
		for (int i = 0; i < fixedIndices.size(); i++)isFixed[fixedIndices[i]] = true;
		for (int i = 0; i < softIndices.size(); i++) if (isFixed[softIndices[i]]){
			printf("Soft variable can not be fixed! \n");
			return 0;
		}

		varIndex.clear();
		varIndex.resize(n);
		for (int i = 0; i < n; i++) varIndex[i] = -1;

		int varCounter = 0;
		for (int i = 0; i < n; i++){
			if (!isFixed[i]){
				varIndex[i] = varCounter;
				varCounter++;
			}
		}

		if (varCounter != (n - fixedIndices.size())){
			printf("Variable counters unexpected! \n");
			return 0;
		}

		std::vector<Eigen::Triplet<double>> stiffnessTriplets;
		stiffnessTriplets.reserve(18 * triangles.size());

		for (int t = 0; t < triangles.size(); t++){
			for (int i = 0; i < 3; i++){
				for (int j = 0; j < 3; j++){
					int vP = varIndex[triangles[t][j]];
					int vN = varIndex[triangles[t][(j + 1) % 3]];
					double edgeWeight = cotangentWeightsSource[3 * t + j] * sourceWeight + cotangentWeightsTarget[3 * t + j] * (1.f - sourceWeight);
					if (vP != -1){
						stiffnessTriplets.push_back(Eigen::Triplet<double>(vP, vP, edgeWeight));
						if (vN != -1) stiffnessTriplets.push_back(Eigen::Triplet<double>(vP, vN, -edgeWeight));
					}
					if (vN != -1){
						stiffnessTriplets.push_back(Eigen::Triplet<double>(vN, vN, edgeWeight));
						if (vP != -1) stiffnessTriplets.push_back(Eigen::Triplet<double>(vN, vP, -edgeWeight));
					}
				}
			}
		}


		stiffnessMatrix.resize(varCounter, varCounter);
		stiffnessMatrix.setFromTriplets(stiffnessTriplets.begin(), stiffnessTriplets.end());


		std::vector<Eigen::Triplet<double>> softWeighTriplets;
		softWeighTriplets.reserve(softIndices.size());
		for (int i = 0; i < softIndices.size(); i++){
			int vi = varIndex[softIndices[i]];
			softWeighTriplets.push_back(Eigen::Triplet<double>(vi, vi, softWeights[i]));
		}

		Eigen::SparseMatrix<double> softWeightsMatrix;
		softWeightsMatrix.resize(varCounter, varCounter);
		softWeightsMatrix.setFromTriplets(softWeighTriplets.begin(), softWeighTriplets.end());

		ARAPCholesky.analyzePattern(stiffnessMatrix + softWeightsMatrix);
		ARAPCholesky.factorize(stiffnessMatrix + softWeightsMatrix);

		return 1;
	}

	void Tri_ARAP_Energy(const std::vector<TriangleIndex> & triangles, const std::vector<Eigen::Vector3d> & restVertices, const std::vector<Eigen::Vector3d> & newVertices, const std::vector<double> & cotangentWeights, std::vector<double> & energy){
		int tCount = (int)triangles.size();
		int vCount = (int)restVertices.size();
		
		std::vector<Eigen::Matrix3d> rotationMatrices(vCount);
		int threads = omp_get_num_procs();
		for (int v = 0; v < vCount; v++) for (int k = 0; k < 3; k++)for (int l = 0; l < 3; l++) rotationMatrices[v](k, l) = 0.0;

		for (int t = 0; t < tCount; t++){
			Eigen::Matrix3d scatterMat;
			for (int k = 0; k < 3; k++)for (int l = 0; l < 3; l++) scatterMat(k, l) = 0.0;
			for (int i = 0; i < 3; i++){
				Eigen::Vector3d restEdge = restVertices[triangles[t][(i + 1) % 3]] - restVertices[triangles[t][i]];
				Eigen::Vector3d newEdge = newVertices[triangles[t][(i + 1) % 3]] - newVertices[triangles[t][i]];
				for (int k = 0; k < 3; k++)for (int l = 0; l < 3; l++) scatterMat(k, l) += restEdge[k] * newEdge[l] * cotangentWeights[3 * t + i];
			}
			for (int i = 0; i < 3; i++){
				rotationMatrices[triangles[t][i]] += scatterMat;
			}
		}

#pragma omp parallel for num_threads( threads )
		for (int v = 0; v < vCount; v++){
			Eigen::JacobiSVD<Eigen::Matrix3d> mSVD(rotationMatrices[v], Eigen::ComputeFullU | Eigen::ComputeFullV);
			Eigen::Matrix3d U = mSVD.matrixU();
			Eigen::Matrix3d V = mSVD.matrixV();
			if (U.determinant()*V.determinant() < 0.0) for (int k = 0; k < 3; k++) U(k, 2) *= -1.0;
			rotationMatrices[v] = V*U.transpose();
		}
		
		energy.clear();
		energy.resize(vCount,0.0);

		for (int t = 0; t < triangles.size(); t++){
			for (int i = 0; i < 3; i++){
				for (int j = 0; j < 3; j++){
					Eigen::Vector3d newEdge = newVertices[triangles[t][j]] - newVertices[triangles[t][(j + 1) % 3]];
					Eigen::Vector3d restEdge = restVertices[triangles[t][j]] - restVertices[triangles[t][(j + 1) % 3]];
					Eigen::Vector3d rotatedEdge = (rotationMatrices[triangles[t][i]] * restEdge);
					energy[triangles[t][i]] += (newEdge - rotatedEdge).squaredNorm() * cotangentWeights[3 * t + j];
				}
			}
		}

		if (0){ //Debug
			for (int v = 0; v < vCount; v++){
				if (energy[v] < 0.0){
					printf("Negative Vertex Energy!!\n");
				}
			}
		}
	}

	void Tri_ASAP_Energy(const std::vector<TriangleIndex> & triangles, const std::vector<Eigen::Vector3d> & restVertices, const std::vector<Eigen::Vector3d> & newVertices, const std::vector<double> & cotangentWeights, std::vector<double> & energy, double similarityBound = 2.f){
		int tCount = (int)triangles.size();
		int vCount = (int)restVertices.size();
		
		double maxSimilarityScale = cbrt(similarityBound);
		
		std::vector<Eigen::Matrix3d> rotationMatrices(vCount);
		std::vector<double> referenceScale(vCount ,0.0);
		int threads = omp_get_num_procs();
		for (int v = 0; v < vCount; v++) for (int k = 0; k < 3; k++)for (int l = 0; l < 3; l++) rotationMatrices[v](k, l) = 0.0;

		for (int t = 0; t < tCount; t++){
			Eigen::Matrix3d scatterMat;
			double triangleScale = 0.0;
			for (int k = 0; k < 3; k++)for (int l = 0; l < 3; l++) scatterMat(k, l) = 0.0;
			for (int i = 0; i < 3; i++){
				Eigen::Vector3d restEdge = restVertices[triangles[t][(i + 1) % 3]] - restVertices[triangles[t][i]];
				Eigen::Vector3d newEdge = newVertices[triangles[t][(i + 1) % 3]] - newVertices[triangles[t][i]];
				for (int k = 0; k < 3; k++)for (int l = 0; l < 3; l++) scatterMat(k, l) += restEdge[k] * newEdge[l] * cotangentWeights[3 * t + i];
				triangleScale += restEdge.squaredNorm() * cotangentWeights[3 * t + i];
			}
			for (int i = 0; i < 3; i++){
				rotationMatrices[triangles[t][i]] += scatterMat;
				referenceScale[triangles[t][i]] += triangleScale;
			}
		}

#pragma omp parallel for num_threads( threads )
		for (int v = 0; v < vCount; v++){
			Eigen::Matrix3d initialRotation = rotationMatrices[v];
			Eigen::JacobiSVD<Eigen::Matrix3d> mSVD(rotationMatrices[v], Eigen::ComputeFullU | Eigen::ComputeFullV);
			Eigen::Matrix3d U = mSVD.matrixU();
			Eigen::Matrix3d V = mSVD.matrixV();
			if (U.determinant()*V.determinant() < 0.0) for (int k = 0; k < 3; k++) U(k, 2) *= -1.0;
			rotationMatrices[v] = V*U.transpose();
			double trace = (rotationMatrices[v] * initialRotation).trace();
			rotationMatrices[v] *= (trace / referenceScale[v]);
			if (rotationMatrices[v].determinant() > similarityBound){
				rotationMatrices[v] = V*U.transpose()*maxSimilarityScale;
				//printf("Max Similarity Clamped! \n");
			}
			//printf("Similarity %g \n", rotationMatrices[v].determinant());
		}

		energy.clear();
		energy.resize(vCount, 0.0);

		for (int t = 0; t < triangles.size(); t++){
			for (int i = 0; i < 3; i++){
				for (int j = 0; j < 3; j++){
					Eigen::Vector3d newEdge = newVertices[triangles[t][j]] - newVertices[triangles[t][(j + 1) % 3]];
					Eigen::Vector3d restEdge = restVertices[triangles[t][j]] - restVertices[triangles[t][(j + 1) % 3]];
					Eigen::Vector3d rotatedEdge = (rotationMatrices[triangles[t][i]] * restEdge);
					energy[triangles[t][i]] += (newEdge - rotatedEdge).squaredNorm() * cotangentWeights[3 * t + j];
				}
			}
		}

		if (0){ //Debug
			for (int v = 0; v < vCount; v++){
				if (energy[v] < 0.0){
					printf("Negative Vertex Energy!!\n");
				}
			}
		}
	}

	int Tri_ASAP_Solve(const std::vector<TriangleIndex> & triangles, const std::vector<Eigen::Vector3d> & restVertices, const std::vector<double> & cotangentWeights, const std::vector<int> & fixedIndices, const std::vector<int> & softIndices, const std::vector<float> & softWeights, const std::vector<Eigen::Vector3d> & softConstraints, const std::vector<int> & varIndex, const Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> & ARAPCholesky, std::vector<Eigen::Vector3d> & newVertices, double similarityBound = 2.f){
		double maxSimilarityScale = cbrt(similarityBound);
		
		int tCount = (int)triangles.size();
		int vCount = (int)restVertices.size();
		std::vector<Eigen::Matrix3d> rotationMatrices(vCount);
		int threads = omp_get_num_procs();
		for (int v = 0; v < vCount; v++) for (int k = 0; k < 3; k++)for (int l = 0; l < 3; l++) rotationMatrices[v](k, l) = 0.0;

		std::vector<double> referenceScale(vCount, 0.0);
		for (int t = 0; t < tCount; t++){
			Eigen::Matrix3d scatterMat;
			double triangleScale = 0.0;
			for (int k = 0; k < 3; k++)for (int l = 0; l < 3; l++) scatterMat(k, l) = 0.0;
			for (int i = 0; i < 3; i++){
				Eigen::Vector3d restEdge = restVertices[triangles[t][(i + 1) % 3]] - restVertices[triangles[t][i]];
				Eigen::Vector3d newEdge = newVertices[triangles[t][(i + 1) % 3]] - newVertices[triangles[t][i]];
				for (int k = 0; k < 3; k++)for (int l = 0; l < 3; l++) scatterMat(k, l) += restEdge[k] * newEdge[l] * cotangentWeights[3 * t + i];
				triangleScale += restEdge.squaredNorm() * cotangentWeights[3 * t + i];
			}
			for (int i = 0; i < 3; i++){
				rotationMatrices[triangles[t][i]] += scatterMat;
				referenceScale[triangles[t][i]] += triangleScale;
			}
		}

#pragma omp parallel for num_threads( threads )
		for (int v = 0; v < vCount; v++){
			Eigen::Matrix3d initialRotation = rotationMatrices[v];
			Eigen::JacobiSVD<Eigen::Matrix3d> mSVD(rotationMatrices[v], Eigen::ComputeFullU | Eigen::ComputeFullV);
			Eigen::Matrix3d U = mSVD.matrixU();
			Eigen::Matrix3d V = mSVD.matrixV();
			if (U.determinant()*V.determinant() < 0.0) for (int k = 0; k < 3; k++) U(k, 2) *= -1.0;
			rotationMatrices[v] = V*U.transpose();
			double trace = (rotationMatrices[v] * initialRotation).trace();
			rotationMatrices[v] *= (trace / referenceScale[v]);
			if (rotationMatrices[v].determinant() > similarityBound){
				rotationMatrices[v] = V*U.transpose()*maxSimilarityScale;
				//printf("Max Similarity Clamped! \n");
			}
			//printf("Similarity %g \n", rotationMatrices[v].determinant());
		}

		Eigen::MatrixXd rhs;
		int numFree = (int)restVertices.size() - (int)fixedIndices.size();
		rhs.resize(numFree, 3);
		for (int k = 0; k < numFree; k++) for (int l = 0; l < 3; l++) rhs(k, l) = 0.0;
		for (int t = 0; t < triangles.size(); t++){
			for (int i = 0; i < 3; i++){
				for (int j = 0; j < 3; j++){
					int vP = varIndex[triangles[t][j]];
					int vN = varIndex[triangles[t][(j + 1) % 3]];
					Eigen::Vector3d restEdge = restVertices[triangles[t][j]] - restVertices[triangles[t][(j + 1) % 3]];
					Eigen::Vector3d rotatedEdge = (rotationMatrices[triangles[t][i]] * restEdge) * cotangentWeights[3 * t + j];
					if (vP != -1){
						rhs.row(vP) += rotatedEdge;
						if (vN == -1) rhs.row(vP) += newVertices[triangles[t][(j + 1) % 3]] * cotangentWeights[3 * t + j];
					}
					if (vN != -1){
						rhs.row(vN) -= rotatedEdge;
						if (vP == -1) rhs.row(vN) += newVertices[triangles[t][j]] * cotangentWeights[3 * t + j];
					}
				}
			}
		}

#pragma omp parallel for num_threads( threads )
		for (int i = 0; i < softIndices.size(); i++){
			int vi = varIndex[softIndices[i]];
			rhs.row(vi) += softConstraints[i] * softWeights[i];
		}

		Eigen::MatrixXd solution = ARAPCholesky.solve(rhs);
		for (int i = 0; i < restVertices.size(); i++){
			int vi = varIndex[i];
			if (vi != -1) for (int k = 0; k < 3; k++) newVertices[i][k] = solution(vi, k);
		}
	}


	int Tri_ARAP_Solve(const std::vector<TriangleIndex> & triangles, const std::vector<Eigen::Vector3d> & restVertices, const std::vector<double> & cotangentWeights, const std::vector<int> & fixedIndices, const std::vector<int> & softIndices, const std::vector<float> & softWeights, const std::vector<Eigen::Vector3d> & softConstraints, const std::vector<int> & varIndex, const Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> & ARAPCholesky, std::vector<Eigen::Vector3d> & newVertices){
		int tCount = (int)triangles.size();
		int vCount = (int)restVertices.size();
		std::vector<Eigen::Matrix3d> rotationMatrices(vCount);
		int threads = omp_get_num_procs();
		for (int v = 0; v < vCount; v++) for (int k = 0; k < 3; k++)for (int l = 0; l < 3; l++) rotationMatrices[v](k, l) = 0.0;

		for (int t = 0; t < tCount; t++){
			Eigen::Matrix3d scatterMat;
			for (int k = 0; k < 3; k++)for (int l = 0; l < 3; l++) scatterMat(k, l) = 0.0;
			for (int i = 0; i < 3; i++){
				Eigen::Vector3d restEdge = restVertices[triangles[t][(i + 1) % 3]] - restVertices[triangles[t][i]];
				Eigen::Vector3d newEdge = newVertices[triangles[t][(i + 1) % 3]] - newVertices[triangles[t][i]];
				for (int k = 0; k < 3; k++)for (int l = 0; l < 3; l++) scatterMat(k, l) += restEdge[k] * newEdge[l] * cotangentWeights[3 * t + i];
			}
			for (int i = 0; i < 3; i++){
				rotationMatrices[triangles[t][i]] += scatterMat;
			}
		}

#pragma omp parallel for num_threads( threads )
		for (int v = 0; v < vCount; v++){
			Eigen::JacobiSVD<Eigen::Matrix3d> mSVD(rotationMatrices[v], Eigen::ComputeFullU | Eigen::ComputeFullV);
			Eigen::Matrix3d U = mSVD.matrixU();
			Eigen::Matrix3d V = mSVD.matrixV();
			if (U.determinant()*V.determinant() < 0.0) for (int k = 0; k < 3; k++) U(k, 2) *= -1.0;
			rotationMatrices[v] = V*U.transpose();
		}

		Eigen::MatrixXd rhs;
		int numFree = (int)restVertices.size() - (int)fixedIndices.size();
		rhs.resize(numFree, 3);
		for (int k = 0; k < numFree; k++) for (int l = 0; l < 3; l++) rhs(k, l) = 0.0;
		for (int t = 0; t < triangles.size(); t++){
			for (int i = 0; i < 3; i++){
				for (int j = 0; j < 3; j++){
					int vP = varIndex[triangles[t][j]];
					int vN = varIndex[triangles[t][(j + 1) % 3]];
					Eigen::Vector3d restEdge = restVertices[triangles[t][j]] - restVertices[triangles[t][(j + 1) % 3]];
					Eigen::Vector3d rotatedEdge = (rotationMatrices[triangles[t][i]] * restEdge) * cotangentWeights[3 * t + j];
					if (vP != -1){
						rhs.row(vP) += rotatedEdge;
						if (vN == -1) rhs.row(vP) += newVertices[triangles[t][(j + 1) % 3]]*cotangentWeights[3 * t + j];
					}
					if (vN != -1){
						rhs.row(vN) -= rotatedEdge;
						if (vP == -1) rhs.row(vN) += newVertices[triangles[t][j]]*cotangentWeights[3 * t + j];
					}
				}
			}
		}

#pragma omp parallel for num_threads( threads )
		for (int i = 0; i < softIndices.size(); i++){
			int vi = varIndex[softIndices[i]];
			rhs.row(vi) += softConstraints[i]*softWeights[i];
		}

		Eigen::MatrixXd solution = ARAPCholesky.solve(rhs);
		for (int i = 0; i < restVertices.size(); i++){
			int vi = varIndex[i];
			if (vi != -1) for (int k = 0; k < 3; k++) newVertices[i][k] = solution(vi, k);
		}
	}

	int Tri_ARAP_Solve(const std::vector<TriangleIndex> & triangles, const std::vector<Eigen::Vector3d> & restSourceVertices, const std::vector<Eigen::Vector3d> & restTargetVertices, const std::vector<double> & cotangentWeightsSource, const std::vector<double> & cotangentWeightsTarget, const float sourceWeight, const std::vector<int> & fixedIndices, const std::vector<int> & softIndices, const std::vector<float> & softWeights, const std::vector<Eigen::Vector3d> & softConstraints, const std::vector<int> & varIndex, const Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> & ARAPCholesky, std::vector<Eigen::Vector3d> & newVertices){
		int tCount = (int)triangles.size();
		int vCount = (int)restSourceVertices.size();
		int threads = omp_get_num_procs();

		//Upate Rotations
		std::vector<Eigen::Matrix3d> rotationMatricesSource(vCount);
		std::vector<Eigen::Matrix3d> rotationMatricesTarget(vCount);

		{
			//Source Rotations
			for (int v = 0; v < vCount; v++) for (int k = 0; k < 3; k++)for (int l = 0; l < 3; l++) rotationMatricesSource[v](k, l) = 0.0;

			for (int t = 0; t < tCount; t++){
				Eigen::Matrix3d scatterMat;
				for (int k = 0; k < 3; k++)for (int l = 0; l < 3; l++) scatterMat(k, l) = 0.0;
				for (int i = 0; i < 3; i++){
					Eigen::Vector3d restEdge = restSourceVertices[triangles[t][(i + 1) % 3]] - restSourceVertices[triangles[t][i]];
					Eigen::Vector3d newEdge = newVertices[triangles[t][(i + 1) % 3]] - newVertices[triangles[t][i]];
					for (int k = 0; k < 3; k++)for (int l = 0; l < 3; l++) scatterMat(k, l) += restEdge[k] * newEdge[l] * cotangentWeightsSource[3 * t + i];
				}
				for (int i = 0; i < 3; i++){
					rotationMatricesSource[triangles[t][i]] += scatterMat;
				}
			}

#pragma omp parallel for num_threads( threads )
			for (int v = 0; v < vCount; v++){
				Eigen::JacobiSVD<Eigen::Matrix3d> mSVD(rotationMatricesSource[v], Eigen::ComputeFullU | Eigen::ComputeFullV);
				Eigen::Matrix3d U = mSVD.matrixU();
				Eigen::Matrix3d V = mSVD.matrixV();
				if (U.determinant()*V.determinant() < 0.0) for (int k = 0; k < 3; k++) U(k, 2) *= -1.0;
				rotationMatricesSource[v] = V*U.transpose();
			}

			//Target Rotations
			for (int v = 0; v < vCount; v++) for (int k = 0; k < 3; k++)for (int l = 0; l < 3; l++) rotationMatricesTarget[v](k, l) = 0.0;

			for (int t = 0; t < tCount; t++){
				Eigen::Matrix3d scatterMat;
				for (int k = 0; k < 3; k++)for (int l = 0; l < 3; l++) scatterMat(k, l) = 0.0;
				for (int i = 0; i < 3; i++){
					Eigen::Vector3d restEdge = restTargetVertices[triangles[t][(i + 1) % 3]] - restTargetVertices[triangles[t][i]];
					Eigen::Vector3d newEdge = newVertices[triangles[t][(i + 1) % 3]] - newVertices[triangles[t][i]];
					for (int k = 0; k < 3; k++)for (int l = 0; l < 3; l++) scatterMat(k, l) += restEdge[k] * newEdge[l] * cotangentWeightsTarget[3 * t + i];
				}
				for (int i = 0; i < 3; i++){
					rotationMatricesTarget[triangles[t][i]] += scatterMat;
				}
			}

#pragma omp parallel for num_threads( threads )
			for (int v = 0; v < vCount; v++){
				Eigen::JacobiSVD<Eigen::Matrix3d> mSVD(rotationMatricesTarget[v], Eigen::ComputeFullU | Eigen::ComputeFullV);
				Eigen::Matrix3d U = mSVD.matrixU();
				Eigen::Matrix3d V = mSVD.matrixV();
				if (U.determinant()*V.determinant() < 0.0) for (int k = 0; k < 3; k++) U(k, 2) *= -1.0;
				rotationMatricesTarget[v] = V*U.transpose();
			}
		}

		//Update Vertices
		Eigen::MatrixXd rhs;
		int numFree = (int)restSourceVertices.size() - (int)fixedIndices.size();
		rhs.resize(numFree, 3);
		for (int k = 0; k < numFree; k++) for (int l = 0; l < 3; l++) rhs(k, l) = 0.0;
		for (int t = 0; t < triangles.size(); t++){
			for (int i = 0; i < 3; i++){
				for (int j = 0; j < 3; j++){
					int vP = varIndex[triangles[t][j]];
					int vN = varIndex[triangles[t][(j + 1) % 3]];
					Eigen::Vector3d restEdgeSource = restSourceVertices[triangles[t][j]] - restSourceVertices[triangles[t][(j + 1) % 3]];
					Eigen::Vector3d rotatedEdgeSource = (rotationMatricesSource[triangles[t][i]] * restEdgeSource) * cotangentWeightsSource[3 * t + j];
					
					Eigen::Vector3d restEdgeTarget = restTargetVertices[triangles[t][j]] - restTargetVertices[triangles[t][(j + 1) % 3]];
					Eigen::Vector3d rotatedEdgeTarget = (rotationMatricesTarget[triangles[t][i]] * restEdgeTarget) * cotangentWeightsTarget[3 * t + j];

					Eigen::Vector3d rotatedEdge = rotatedEdgeSource*sourceWeight + rotatedEdgeTarget*(1.f - sourceWeight);
					double edgeWeight = cotangentWeightsSource[3 * t + j] * sourceWeight + cotangentWeightsTarget[3 * t + j] * (1.f - sourceWeight);

					if (vP != -1){
						rhs.row(vP) += rotatedEdge;
						if (vN == -1) rhs.row(vP) += newVertices[triangles[t][(j + 1) % 3]]*edgeWeight;
					}
					if (vN != -1){
						rhs.row(vN) -= rotatedEdge;
						if (vP == -1) rhs.row(vN) += newVertices[triangles[t][j]]*edgeWeight;
					}
				}
			}
		}

#pragma omp parallel for num_threads( threads )
		for (int i = 0; i < softIndices.size(); i++){
			int vi = varIndex[softIndices[i]];
			rhs.row(vi) += softConstraints[i]*softWeights[i];
		}

		Eigen::MatrixXd solution = ARAPCholesky.solve(rhs);
		for (int i = 0; i < restSourceVertices.size(); i++){
			int vi = varIndex[i];
			if (vi != -1) for (int k = 0; k < 3; k++) newVertices[i][k] = solution(vi, k);
		}
	}


	void SetupARAPEnergy(const SimpleMesh & referenceMesh, double & softScale, Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> & ARAPCholesky, Eigen::SparseMatrix<double> & stiffnessMatrix, std::vector<double> & cotangentWeights,
		std::vector<int> & softIndices, std::vector<float> & softWeights, std::vector<float> & vertexWeights, std::vector<int> & freeVariableIndex, std::vector<int> & fixedIndices) {

		const std::vector<Point3D<double>> & _referenceVertices = referenceMesh.vertices;
		std::vector<Eigen::Vector3d> referenceVertices;
		CopyPointClasses(_referenceVertices, referenceVertices);

		const std::vector<TriangleIndex> & referenceTriangles = referenceMesh.triangles;

		int surfaceVertices = referenceVertices.size();
		std::vector<float> triangleWeights;

		{//Correct softScale (softScale = 10 is the default for a mesh of radius 1)
			Eigen::Vector3d center(0.f, 0.f, 0.f);
			double area = 0.f;
			for (int i = 0; i < referenceTriangles.size(); i++) {
				Eigen::Vector3d n = (referenceVertices[referenceTriangles[i][1]] - referenceVertices[referenceTriangles[i][0]]).cross(referenceVertices[referenceTriangles[i][2]] - referenceVertices[referenceTriangles[i][0]]);
				Eigen::Vector3d c = Eigen::Vector3d(referenceVertices[referenceTriangles[i][0]] + referenceVertices[referenceTriangles[i][1]] + referenceVertices[referenceTriangles[i][2]]) / 3.0;
				double a = (n).norm();
				center += c*a, area += a;
			}
			center /= area;
			double scale = 0;
			for (int i = 0; i < referenceVertices.size(); i++) scale = std::max<double>(scale, (referenceVertices[i] - center).norm());
			softScale *= scale;
		}
		{ //Set Geometry
			triangleWeights.resize(referenceTriangles.size());
			vertexWeights.resize(referenceVertices.size(), 0.f);
			for (int i = 0; i < referenceTriangles.size(); i++) {
				float area = ((referenceVertices[referenceTriangles[i][1]] - referenceVertices[referenceTriangles[i][0]]).cross(referenceVertices[referenceTriangles[i][2]] - referenceVertices[referenceTriangles[i][0]])).norm() / 2.f;
				triangleWeights[i] = area;
				for (int v = 0; v < 3; v++) vertexWeights[referenceTriangles[i][v]] += area;
			}
			float maxWeight = 0.f;
			for (int i = 0; i < referenceVertices.size(); i++) maxWeight = std::max<float>(maxWeight, vertexWeights[i]);
			for (int i = 0; i < referenceVertices.size(); i++) vertexWeights[i] /= maxWeight;

			TriangleCotangentWeights(referenceVertices, referenceTriangles, cotangentWeights);
		}
		{//Set Soft Vertices
			softIndices.resize(surfaceVertices);
			for (int i = 0; i < surfaceVertices; i++) softIndices[i] = i;

			softWeights.resize(surfaceVertices);
			for (int i = 0; i < surfaceVertices; i++) softWeights[i] = vertexWeights[i] * softScale;
		}
		Tri_ARAP_Setup((int)referenceVertices.size(), referenceTriangles, cotangentWeights, fixedIndices, softIndices, softWeights, freeVariableIndex, ARAPCholesky, stiffnessMatrix);
	}

}

#endif //SURFACE_ARAP_INCLUDED