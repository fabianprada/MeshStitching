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


#ifndef SIMPLE_FEM_INCLUDED
#define SIMPLE_FEM_INCLUDED 
#include <Eigen/Sparse>
class FEMData{
public:
	FEMData(const std::vector<TriangleIndex> & triangles, const std::vector<Point3D<double>> & vertices){
		InitializeFEM(triangles, vertices,mass,stiffness,g);
	}
	std::vector<SquareMatrix<double, 2>> g;
	Eigen::SparseMatrix<double> mass;
	Eigen::SparseMatrix<double> stiffness;
	void InitializeFEM(const std::vector<TriangleIndex> & triangles, const std::vector<Point3D<double>> & vertices,Eigen::SparseMatrix<double> & massMatrix,Eigen::SparseMatrix<double> & stiffnessMatrix, std::vector<SquareMatrix<double, 2>> & g, bool normalizeMetric = true, bool lumpMass = false);
};


SquareMatrix<double,2> GetMetricFromEmbedding(const Point3D<double> & corners_0, const Point3D<double> & corners_1, const Point3D<double> & corners_2){
	SquareMatrix<double,2>  g;
	g(0, 0) = Point3D<double>::Dot(corners_1 - corners_0,corners_1 - corners_0);
	g(1, 0) = g(0, 1) = Point3D<double>::Dot(corners_1 - corners_0,corners_2 - corners_0);
	g(1, 1) = Point3D<double>::Dot(corners_2 - corners_0,corners_2 - corners_0);
	return g;
}

SquareMatrix<double,3> GetMassMatrix(const SquareMatrix<double,2>  & g, bool lump = false) {
	double area = sqrt(g.determinant()) / 2.0;
	SquareMatrix<double,3> mass;
	if (lump) {
		mass(0, 0) = mass(1, 1) = mass(2, 2) = area*(1.0 / 3.0);
		mass(0, 1) = mass(0, 2) = mass(1, 0) = mass(1, 2) = mass(2, 0) = mass(2, 1) = 0.0;
	}
	else {
		mass(0, 0) = mass(1, 1) = mass(2, 2) = area*(1.0 / 6.0);
		mass(0, 1) = mass(0, 2) = mass(1, 0) = mass(1, 2) = mass(2, 0) = mass(2, 1) = area*(1.0 / 12.0);
	}
	return mass;
}

void GetVertexArea(const SimpleMesh & mesh, std::vector<double> & vertexArea) {
	vertexArea.resize(mesh.vertices.size(), 0.0);
	for (int i = 0; i < mesh.triangles.size(); i++) {
		double area = Point3D<double>::Length(Point3D<double>::CrossProduct((mesh.vertices[mesh.triangles[i][1]] - mesh.vertices[mesh.triangles[i][0]]),(mesh.vertices[mesh.triangles[i][2]] - mesh.vertices[mesh.triangles[i][0]]))) / 2.0;
		for (int j = 0; j < 3; j++) vertexArea[mesh.triangles[i][j]] += area / 3.0;
	}
}

double GetMeshArea(const SimpleMesh & mesh) {
	double area = 0.0;
	for (int i = 0; i < mesh.triangles.size(); i++) {
		area += Point3D<double>::Length(Point3D<double>::CrossProduct((mesh.vertices[mesh.triangles[i][1]] - mesh.vertices[mesh.triangles[i][0]]),(mesh.vertices[mesh.triangles[i][2]] - mesh.vertices[mesh.triangles[i][0]]))) / 2.0;
	}
	return area;
}


SquareMatrix<double,3> GetStiffnessMatrix(const SquareMatrix<double,2>  & g) {
	double area = sqrt(g.determinant()) / 2.0;
	SquareMatrix<double,3> stiffness;
	Point2D<double> d[3];
	d[0] = Point2D<double>(-1.0, -1.0);
	d[1] = Point2D<double>(1.0, 0.0);
	d[2] = Point2D<double>(0.0, 1.0);
	SquareMatrix<double,2>  g_inv = g.inverse();
	for (int i = 0; i < 3; i++)for (int j = 0; j < 3; j++) stiffness(i, j) = Point2D<double>::Dot(d[i],g_inv*d[j])*area;
	return stiffness;
}

void FEMData::InitializeFEM(const std::vector<TriangleIndex> & triangles, const std::vector<Point3D<double>> & vertices,Eigen::SparseMatrix<double> & massMatrix,Eigen::SparseMatrix<double> & stiffnessMatrix, std::vector<SquareMatrix<double, 2>> & g, bool normalizeMetric, bool lumpMass){
	

	std::vector<Eigen::Triplet<float>> massTriplets;
	std::vector<Eigen::Triplet<float>> stiffnessTriplets;
	double cumArea = 0;


	for (int t = 0; t <triangles.size(); t++) {
		SquareMatrix<double, 2>  metric = GetMetricFromEmbedding(vertices[triangles[t][0]],vertices[triangles[t][1]],vertices[triangles[t][2]]);
		cumArea += sqrt(metric.determinant()) / 2.0;

		SquareMatrix<double, 3> mass = GetMassMatrix(metric, lumpMass);
		SquareMatrix<double, 3> stiffness = GetStiffnessMatrix(metric);

		int vIndex[3] = {triangles[t][0],triangles[t][1],triangles[t][2] };
		for (int k = 0; k < 3; k++)for (int l = 0; l < 3; l++) {
			massTriplets.push_back(Eigen::Triplet<float>(vIndex[k], vIndex[l], mass(k, l)));
			stiffnessTriplets.push_back(Eigen::Triplet<float>(vIndex[k], vIndex[l], stiffness(k, l)));
		}
	}

	massMatrix.resize(vertices.size(),vertices.size());
	massMatrix.setFromTriplets(massTriplets.begin(), massTriplets.end());

	stiffnessMatrix.resize(vertices.size(),vertices.size());
	stiffnessMatrix.setFromTriplets(stiffnessTriplets.begin(), stiffnessTriplets.end());

	if (normalizeMetric){
		massMatrix /= cumArea;
	}
}

void SmoothSignal(FEMData * fem, Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> * solver, std::vector<Point3D<double>> & signal, bool normalize = false) {
	int vCount = signal.size();
	Eigen::MatrixXd x = Eigen::MatrixXd::Zero(vCount, 3);
	Eigen::MatrixXd b = Eigen::MatrixXd::Zero(vCount, 3);

#pragma omp parallel for num_threads( omp_get_num_procs() )
	for (int c = 0; c < 3; c++) for (int i = 0; i < vCount; i++) x(i, c) = signal[i][c];

	b = fem->mass * x;
	x = solver->solve(b);

#pragma omp parallel for num_threads( omp_get_num_procs() )
	for (int c = 0; c < 3; c++) for (int i = 0; i < vCount; i++) signal[i][c] = x(i, c);

	if (normalize) for (int i = 0; i < signal.size(); i++) signal[i] /= Point3D<double>::Length(signal[i]);
}

void SmoothSignal(FEMData * fem, Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> * solver, std::vector<double> & signal, bool compensateDC = false) {
	int vCount = signal.size();
	Eigen::VectorXd x = Eigen::VectorXd::Zero(vCount);
	Eigen::VectorXd b = Eigen::VectorXd::Zero(vCount);
	for (int i = 0; i < vCount; i++) x[i] = signal[i];
	b = fem->mass * x;
	double old_dc = 0;
	for (int i = 0; i < vCount; i++) old_dc += b[i];
	x = solver->solve(b);
	double new_dc = 0;
	for (int i = 0; i < vCount; i++) new_dc += x[i];
	double dc_compensation = (old_dc - new_dc) / vCount;
	for (int i = 0; i < vCount; i++) signal[i] = x[i];
	if (compensateDC) for (int i = 0; i < signal.size(); i++) signal[i] += dc_compensation;
}


#endif //SIMPLE_FEM_INCLUDED