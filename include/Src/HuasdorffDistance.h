#pragma once

#include <Misha/Geometry.h>
#include <igl/point_mesh_squared_distance.h>
#include <igl/AABB.h>
#include <igl/barycentric_coordinates.h>

void HausdorffDistance(const std::vector<Point3D<double>> & verticesSource, const std::vector<Point3D<double>> & verticesTarget, const std::vector<TriangleIndex> trianglesTarget, double & meanMin, double & hausdorff, std::vector<double> & distance) {

	igl::AABB<Eigen::MatrixXf, 3> kdTree;
	Eigen::MatrixXf vertexMatrix;
	Eigen::MatrixXi triangleMatrix;

	vertexMatrix.resize(verticesTarget.size(), 3);
	for (int i = 0; i < verticesTarget.size(); i++)for (int c = 0; c < 3; c++) vertexMatrix(i, c) = verticesTarget[i][c];
	triangleMatrix.resize(trianglesTarget.size(), 3);
	for (int i = 0; i < trianglesTarget.size(); i++)for (int c = 0; c < 3; c++) triangleMatrix(i, c) = trianglesTarget[i][c];
	kdTree.init(vertexMatrix, triangleMatrix);


	Eigen::VectorXf sqrD;
	Eigen::VectorXi I;
	Eigen::MatrixXf C;
	Eigen::MatrixXf P;
	P.resize(verticesSource.size(), 3);
	for (int i = 0; i < verticesSource.size(); i++)P.row(i) = Eigen::Vector3f(verticesSource[i][0], verticesSource[i][1], verticesSource[i][2]);
	kdTree.squared_distance(vertexMatrix, triangleMatrix, P, sqrD, I, C);

	distance.resize(verticesSource.size());
	hausdorff = 0;
	meanMin = 0;
	for (int i = 0; i < verticesSource.size(); i++) {
		double currentDistance = sqrt(sqrD[i]);
		distance[i] = currentDistance;
		hausdorff = std::max<double>(hausdorff, currentDistance);
		meanMin += currentDistance;
	}
	meanMin /= double(verticesSource.size());
}
