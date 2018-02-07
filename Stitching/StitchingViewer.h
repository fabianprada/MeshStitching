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


#include <Src/StitchingGrid.h>
#include "SurfaceVisualization.inl"
#include <Src/SimpleFEM.h>
#include <Src/SurfaceARAP.h>
#include <Misha/Timer.h>
#include <Src/HuasdorffDistance.h>

#define USE_CONSTRAINED_DELANAUY_CLIPPING 1


enum {
	SELECTION_STAGE,
	CLIPPING_STAGE,
	STAGE_COUNT
};

enum {
	AUTOMATIC_SELECTION,
	MANUAL_SELECTION, 
	SELECTION_MODE_COUNT
};

enum{
	REMESHING_CANDIDATES_VISUALIZATION,
	GROWING_REGION_VISUALIZATION,
	VISUALIZATION_COUNT
};

enum{
	SOURCE_MESH,
	TARGET_MESH,
	STITCHED_MESH
};

enum{
	STOP_AT_SIMILAR_FACES,
	STOP_AT_NON_WEAK_FACES,
	STOP_COUNT
};

enum{
	SINGLE_INTERSECTION,
	COMMON_INTESECTION,
	INTERSECTION_COUNT
};

enum {
	STITCHED_SIMPLE_MESH,
	STITCHED_SEGMENTED_MESH
};

class FaceCoords{
public:
	FaceCoords(double p_coords[9], int p_index){
		index = p_index;
		for (int i = 0; i < 9; i++) coords[i] = p_coords[i];
	}
	double coords[9];
	int index;
};

class FaceCoordsComparison{
public:
	bool operator() (const FaceCoords & face1, const FaceCoords & face2) const{
		for (int i = 0; i < 9; i++){
			if (face1.coords[i] < face2.coords[i]){
				return true;
			}
			else if (face2.coords[i] < face1.coords[i]){
				return false;
			}
		}
		return false;
	}
};

struct StitchingViewer
{
	static int currentSystemStage;
	static SurfaceVisualization sv;
	static SimpleMesh inputMesh[2];
	static SimpleMesh stitchedMesh;
	static ColoredMesh coloredStitchedMesh;
	static int meshMode;
	static int meshCount;

	//Seed selection parameters 
	static int selectionMode;
	static double maxHaussdorffDistance;
	static double maxVisibleHaussdorffDistance;
	static double minVisibleHaussdorffDistance;
	static double maxMeanCurvature;


	static double gridScaling;
	static Point3D<double> gridTranslation;
	static int gridRes;
	static int gridDepth;
	static Point3D<float> visualizationCenter;
	static float visualizationRadius;
	static unsigned long referenceVoxel[3];
	static int voxelIntersectionRadius;

	static float hausdorffThreshold;
	static float visibilityThreshold;

	static std::unordered_map<unsigned long, EdgeIntersections> edgeIntersections[2];
	static std::unordered_map<unsigned long, FaceMap> faceMaps[2];
	static std::unordered_set<unsigned long> voxelKeys[2];
	static std::unordered_map<unsigned long long, SamplePoint> vertexSample[2];
	static std::unordered_map<unsigned long long, unsigned long long> adjacentVertexMap[2][3];
	static std::unordered_set<int> remeshingVertices[2];
	static std::vector<unsigned long> voxelSeeds;

	static SimpleMesh clippedMesh[2];
	static std::vector<unsigned long> clippedTriangleVoxelMap[2];
	static std::unordered_map<unsigned long long, int>  clippedMeshVertexIndex[2];
	static std::unordered_map<unsigned long long, int>  clippedEdgeTriangleMap[2];

	static std::vector<std::vector<Point3D<double>>> loopPoints[2];
	static std::unordered_map<unsigned long, bool> roinBoundaryFaceOrientation;

	static std::unordered_set<unsigned long> addedVoxels;
	static ComponentMap refMap[2];
	static std::unordered_map<unsigned long long, SamplePoint>  _vertexSample[2];

	static FEMData * FEM[2];
	static Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> * normalSolver[2];
	static Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> * offsetSolver;
	static double normalSmoothWeight;
	static double offsetSmoothWeight;

	static std::unordered_set<unsigned long> commonVoxelKeys;
	static std::unordered_set<unsigned long> commonWeakFaces;
	static std::unordered_set<unsigned long> commonSimilarFaces;
	static std::vector<int> forwardMap[2];
	static std::vector<int> backwardMap[2];
	
	static void UpdateRemeshingVertices();
	static void UpdateVoxelSeeds();
	static void PerformRegionGrowing();
	static void ClippingAndStitching();
	static void SetForwardAndBackwardMaps();
	static void SmoothNonFixedVertices();
	static void ExportMesh(const char * outputName, int exportMode);

	static void NormalProjection();
	static int Init(const char* sourceName, const char* targetName, const int p_gridDepth, const int geodesicMarkers, bool disableNormalOffset, bool vertexJitter);
	static int Process(const char * outputMeshName);
	static void ToggleSelectionModeCallBack(Visualization*, const char*);
	static void ToggleMeshModeCallBack(Visualization*, const char*);
	static void PerformRegionGrowingCallBack(Visualization*, const char* prompt);
	static void StitchingCallBack(Visualization*, const char* prompt);
	static void ExportMeshCallBack(Visualization*, const char* prompt);
	static void MaxHausdorffDistanceCallBack(Visualization*, const char* prompt);
	static void MaxVisibleHausdorffDistanceCallBack(Visualization*, const char* prompt);
	static void MinVisibleHausdorffDistanceCallBack(Visualization*, const char* prompt);
	static void MaxMeanCurvatureCallBack(Visualization*, const char* prompt);
	static void Idle        ( void );
	static void KeyboardFunc( unsigned char key , int x , int y );
	static void SpecialFunc ( int key, int x, int y );
	static void Display     ( void );
	static void Reshape     ( int w , int h );
	static void MouseFunc   ( int button , int state , int x , int y );
	static void MotionFunc  ( int x , int y );
	static void SetVisualizationMesh(const SimpleMesh & mesh);
	static void SetVisualizationMesh(const ColoredMesh & mesh);
protected:
	static int _busyCursor;
	static std::vector< int > _valence;
};
int																				StitchingViewer::currentSystemStage = SELECTION_STAGE;
SurfaceVisualization															StitchingViewer::sv;

SimpleMesh																		StitchingViewer::inputMesh[2];
SimpleMesh																		StitchingViewer::stitchedMesh;
ColoredMesh																		StitchingViewer::coloredStitchedMesh;

std::unordered_map<unsigned long, FaceMap>										StitchingViewer::faceMaps[2];
std::unordered_set<unsigned long>												StitchingViewer::voxelKeys[2];
std::unordered_map<unsigned long, EdgeIntersections>							StitchingViewer::edgeIntersections[2];
std::unordered_map<unsigned long long, SamplePoint>								StitchingViewer::vertexSample[2];
std::unordered_map<unsigned long long, unsigned long long>						StitchingViewer::adjacentVertexMap[2][3];
std::unordered_set<int>															StitchingViewer::remeshingVertices[2];
std::vector<unsigned long>														StitchingViewer::voxelSeeds;
SimpleMesh																		StitchingViewer::clippedMesh[2];
std::unordered_map<unsigned long long, int>										StitchingViewer::clippedMeshVertexIndex[2];
std::unordered_map<unsigned long long, int>										StitchingViewer::clippedEdgeTriangleMap[2];
std::vector<unsigned long>														StitchingViewer::clippedTriangleVoxelMap[2];

Point3D<float>																	StitchingViewer::visualizationCenter;
float																			StitchingViewer::visualizationRadius;

int																				StitchingViewer::_busyCursor = GLUT_CURSOR_WAIT;
int																				StitchingViewer::meshMode = SOURCE_MESH;
int																				StitchingViewer::meshCount = 2;

std::unordered_set<unsigned long>												StitchingViewer::commonVoxelKeys;
std::unordered_set<unsigned long>												StitchingViewer::commonWeakFaces;
std::unordered_set<unsigned long>												StitchingViewer::commonSimilarFaces;

double																			StitchingViewer::gridScaling;
Point3D<double>																	StitchingViewer::gridTranslation;
int																				StitchingViewer::gridRes;
int																				StitchingViewer::gridDepth;
unsigned long																	StitchingViewer::referenceVoxel[3];
int																				StitchingViewer::voxelIntersectionRadius = 4;


FEMData *																		StitchingViewer::FEM[2];
Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> *							StitchingViewer::normalSolver[2];
Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> *							StitchingViewer::offsetSolver;

double																			StitchingViewer::normalSmoothWeight = (double)1e-4;
double																			StitchingViewer::offsetSmoothWeight = (double)1e-4; 

float																			StitchingViewer::hausdorffThreshold = 20.0;
float																			StitchingViewer::visibilityThreshold = 40.0;

//std::vector<unsigned long> StitchingViewer::roiVoxels;
std::unordered_map<unsigned long, bool>											StitchingViewer::roinBoundaryFaceOrientation;

std::unordered_set<unsigned long>												StitchingViewer::addedVoxels;
ComponentMap																	StitchingViewer::refMap[2];
std::unordered_map<unsigned long long, SamplePoint>								StitchingViewer::_vertexSample[2];

std::vector<int>																StitchingViewer::forwardMap[2];
std::vector<int>																StitchingViewer::backwardMap[2];
std::vector<std::vector<Point3D<double>>>										StitchingViewer::loopPoints[2];

int																				StitchingViewer::selectionMode = MANUAL_SELECTION;
double																			StitchingViewer::maxHaussdorffDistance = 100.0;
double																			StitchingViewer::maxVisibleHaussdorffDistance = 50.0;
double																			StitchingViewer::minVisibleHaussdorffDistance = 30.0;
double																			StitchingViewer::maxMeanCurvature = 0.02; 

void StitchingViewer::ToggleSelectionModeCallBack(Visualization*, const char*) {
	selectionMode = (selectionMode + 1) % SELECTION_MODE_COUNT;

	if (selectionMode == MANUAL_SELECTION) {
		sv.info[4][0] = 0;
		sv.info[3][0] = 0;
		sv.info[2][0] = 0;
		sv.info[1][0] = 0;
		sprintf(sv.info[0], "Manual vertex selection");

		remeshingVertices[0].clear();
		remeshingVertices[1].clear();
		sv.selectedVertices.clear();
	}
	else if (selectionMode == AUTOMATIC_SELECTION){
		UpdateRemeshingVertices();
	}

}
void  StitchingViewer::ToggleMeshModeCallBack(Visualization*, const char*){
	meshMode = (meshMode + 1) % meshCount;
	if(meshMode == SOURCE_MESH || meshMode == TARGET_MESH) SetVisualizationMesh(inputMesh[meshMode]);
	else if (meshMode == STITCHED_MESH) SetVisualizationMesh(coloredStitchedMesh);
	sv.updateMesh();
} 
 

void StitchingViewer::NormalProjection(){
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ARAPCholesky;
	Eigen::SparseMatrix<double> stiffnessMatrix;
	std::vector<double> cotangentWeights;
	std::vector<int> softIndices;
	std::vector<float> softWeights;
	std::vector<float> vertexWeights;
	std::vector<int> freeVariableIndex;
	std::vector<int> fixedIndices;
	double softScale = 10.f;

	SurfaceARAP::SetupARAPEnergy(inputMesh[0], softScale, ARAPCholesky, stiffnessMatrix, cotangentWeights, softIndices, softWeights, vertexWeights, freeVariableIndex, fixedIndices);

	SimpleMesh deformationMesh = inputMesh[0];
	std::vector<Eigen::Vector3d> referenceVertices;
	CopyPointClasses(inputMesh[0].vertices, referenceVertices);

	int arapIters = 10;
	int maxProjIters = 10;
	for (int projIter = 0; projIter < maxProjIters; projIter++) {
		UpdateNormals(deformationMesh);
		SmoothSignal(FEM[0], normalSolver[0], deformationMesh.normals, true);
		printf("Projection Iter %d of %d \r", projIter, maxProjIters);
		std::vector<SamplePoint> mapping;
		ThresholdedProjectionMap(deformationMesh, inputMesh[1], 200.f, 0.f, mapping, 2.f);

		std::vector<Eigen::Vector3d> softConstraints;
		softConstraints.resize(deformationMesh.vertices.size());

		std::vector<float> updateSoftWeights(deformationMesh.vertices.size());

		for (int i = 0; i < deformationMesh.vertices.size(); i++) {
			int tIndex = mapping[i].tIndex;
			if (tIndex != -1) {
				tIndex = mapping[i].tIndex;
				Point2D<double> baricentric = mapping[i].baricentricCoord;
				Point3D<double> projectedPosition = inputMesh[1].vertices[inputMesh[1].triangles[tIndex][0]] * (1.f - baricentric[0] - baricentric[1]) + inputMesh[1].vertices[inputMesh[1].triangles[tIndex][1]] * baricentric[0] + inputMesh[1].vertices[inputMesh[1].triangles[tIndex][2]] * baricentric[1];
				if (Point3D<double>::Length(projectedPosition - deformationMesh.vertices[i]) < 40.f) {
					softConstraints[i] = Eigen::Vector3d(projectedPosition[0], projectedPosition[1], projectedPosition[2]);
					updateSoftWeights[i] = vertexWeights[i] * softScale;
				}
			}
			else {
				updateSoftWeights[i] = 0.f;
			}
		}

		std::vector<Eigen::Triplet<double>> softWeighTriplets;
		softWeighTriplets.reserve(softIndices.size());
		for (int i = 0; i < softIndices.size(); i++) softWeighTriplets.push_back(Eigen::Triplet<double>(i, i, updateSoftWeights[i]));

		Eigen::SparseMatrix<double> softWeightsMatrix;
		softWeightsMatrix.resize(deformationMesh.vertices.size(), deformationMesh.vertices.size());
		softWeightsMatrix.setFromTriplets(softWeighTriplets.begin(), softWeighTriplets.end());

		ARAPCholesky.factorize(stiffnessMatrix + softWeightsMatrix);

		std::vector<Eigen::Vector3d> currentVertices;
		CopyPointClasses(deformationMesh.vertices, currentVertices);

		for (int arapIter = 0; arapIter < arapIters; arapIter++) {
			SurfaceARAP::Tri_ARAP_Solve(deformationMesh.triangles, referenceVertices, cotangentWeights, fixedIndices, softIndices, updateSoftWeights, softConstraints, freeVariableIndex, ARAPCholesky, currentVertices);
			for (int i = 0; i < currentVertices.size(); i++) softConstraints[i] = currentVertices[softIndices[i]];
		}

		CopyPointClasses(currentVertices, deformationMesh.vertices);
	}
	inputMesh[0].vertices = deformationMesh.vertices;

	UpdateNormals(inputMesh[0]);
	SmoothSignal(FEM[0], normalSolver[0], inputMesh[0].normals, true);
}

int StitchingViewer::Init(const char* sourceName, const char* targetName, const int p_gridDepth, const int geodesicMarkers, bool disableNormalOffset, bool vertexJitter)
{
	printf("Reading Input\n");
	if (!ReadSimpleMesh(inputMesh[0], sourceName)){
		printf("Unable to read source data : %s\n", sourceName);
		return 0;
	}
	if (!ReadSimpleMesh(inputMesh[1], targetName)){
		printf("Unable to read target data : %s\n", targetName);
		return 0;
	}

	if (vertexJitter){
		printf("Applying vertex jitter \n");
		for (int i = 0; i < 2; i++){
			Point3D< double > min_v = inputMesh[i].vertices[0];
			Point3D< double > max_v = inputMesh[i].vertices[0];
			for (int v = 0; v < inputMesh[i].vertices.size(); v++) for (int c = 0; c < 3; c++){
				min_v[c] = std::min<double>(min_v[c], inputMesh[i].vertices[v][c]);
				max_v[c] = std::max<double>(max_v[c], inputMesh[i].vertices[v][c]);
			}
			double bboxDiagonal = Point3D<double>::Length(max_v - min_v);
			for (int v = 0; v < inputMesh[i].vertices.size(); v++){
				Point3D<double> randVector(double(rand()) / double(RAND_MAX) - 0.5, double(rand()) / double(RAND_MAX) - 0.5, double(rand()) / double(RAND_MAX) - 0.5);
				randVector /= Point3D<double>::Length(randVector);
				inputMesh[i].vertices[v] += (randVector*bboxDiagonal*0.000001);
			}
		}
		
	}

	for (int i = 0; i < 2; i++){
		FEM[i] = new FEMData(inputMesh[i].triangles, inputMesh[i].vertices);
		normalSolver[i] = new Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>(FEM[i]->mass + FEM[i]->stiffness * normalSmoothWeight);
	}

	for (int i = 0; i < 2; i++){
		UpdateNormals(inputMesh[i]);
		SmoothSignal(FEM[i], normalSolver[i], inputMesh[i].normals, true);
	}

	if(!disableNormalOffset){
		NormalProjection();
	}

	gridDepth = p_gridDepth;
	gridRes = 2 << gridDepth;
	CenterMeshes(inputMesh[0], inputMesh[1], gridScaling, gridTranslation);

	Timer t_FaceMaps;
	
	printf("Computing face maps.. \n");
	for (int i = 0; i < 2; i++){
		if (!InitializeFaceMaps(inputMesh[i], gridDepth, faceMaps[i], edgeIntersections[i], vertexSample[i], adjacentVertexMap[i])) return 0;
		if (!InitializeVoxelKeys(faceMaps[i], voxelKeys[i])) return 0;
	}

	printf("Face maps construction time %.2f(s)\n", t_FaceMaps.elapsed());

	commonVoxelKeys = voxelKeys[0];
	for (auto voxelIter = voxelKeys[1].begin(); voxelIter != voxelKeys[1].end(); voxelIter++) commonVoxelKeys.insert(*voxelIter);

	referenceVoxel[0] = referenceVoxel[1] = referenceVoxel[2] = -1;

	FindSimilarFaces(faceMaps[0], faceMaps[1], commonSimilarFaces);
	std::vector<unsigned long> similarFacesVector(commonSimilarFaces.begin(), commonSimilarFaces.end());
	sv.faces = similarFacesVector;

	{
		Point3D< float > center;
		float area = 0.f;
		for (int i = 0; i<inputMesh[0].triangles.size(); i++)
		{
			Point3D< float > n = Point3D< float >::CrossProduct(inputMesh[0].vertices[inputMesh[0].triangles[i][1]] - inputMesh[0].vertices[inputMesh[0].triangles[i][0]], inputMesh[0].vertices[inputMesh[0].triangles[i][2]] - inputMesh[0].vertices[inputMesh[0].triangles[i][0]]);
			Point3D< float > c = (inputMesh[0].vertices[inputMesh[0].triangles[i][0]] + inputMesh[0].vertices[inputMesh[0].triangles[i][1]] + inputMesh[0].vertices[inputMesh[0].triangles[i][2]]) / 3.f;
			float a = (float)Length(n);
			center += c*a, area += a;
		}
		center /= area;
		float max = 0.f;
		for (int i = 0; i<inputMesh[0].vertices.size(); i++) max = std::max< float >(max, (float)Point3D< float >::Length(inputMesh[0].vertices[i] - center));
		visualizationCenter = center;
		visualizationRadius = max;
		printf("Visualization Center %f %f %f \n", visualizationCenter[0], visualizationCenter[1], visualizationCenter[2]);
		printf("Visualization Radius %f \n",visualizationRadius);
		SetVisualizationMesh(inputMesh[0]);
	}

	sv.integerGridScaling = 1.0 / (float(gridRes)*visualizationRadius);
	sv.integerGridTranslation = -visualizationCenter/visualizationRadius;

	sv.info.resize( 6 );
	for( int i=0 ; i<sv.info.size() ; i++ ) sv.info[i] = new char[512] , sv.info[i][0] = 0;

	sprintf(sv.info[0],"Manual vertex selection");

	sv.callBacks.push_back(Visualization::KeyboardCallBack(&sv, 'g', "grow", PerformRegionGrowingCallBack));
	sv.callBacks.push_back(Visualization::KeyboardCallBack(&sv, 's', "stitch", StitchingCallBack));
	sv.callBacks.push_back(Visualization::KeyboardCallBack(&sv, 't', "toggle mesh", ToggleMeshModeCallBack));
	sv.callBacks.push_back(Visualization::KeyboardCallBack(&sv, 'k', "toggle selection mode", ToggleSelectionModeCallBack));
	sv.callBacks.push_back(Visualization::KeyboardCallBack(&sv, 'N', "max visible distance", "Value",MaxVisibleHausdorffDistanceCallBack));
	sv.callBacks.push_back(Visualization::KeyboardCallBack(&sv, 'n', "min visible distance", "Value", MinVisibleHausdorffDistanceCallBack));
	sv.callBacks.push_back(Visualization::KeyboardCallBack(&sv, 'Y', "max hausdorff distance", "Value", MaxHausdorffDistanceCallBack));
	sv.callBacks.push_back(Visualization::KeyboardCallBack(&sv, 'C', "max mean curvature", "Value", MaxMeanCurvatureCallBack));
	sv.callBacks.push_back(Visualization::KeyboardCallBack(&sv, 'o', "export mesh", "File name", ExportMeshCallBack));

	return 1;
}


void StitchingViewer::SetVisualizationMesh(const SimpleMesh & mesh){
	sv.vertices.clear();
	sv.triangles.clear();
	sv.colors.clear();

	sv.vertices.resize(mesh.vertices.size());
	for (int i = 0; i < sv.vertices.size(); i++)sv.vertices[i] = Point3D<float>(mesh.vertices[i]);
	sv.triangles = mesh.triangles;
	sv.colors.resize(3 * sv.triangles.size());
	for (int i = 0; i<3 * sv.triangles.size(); i++) sv.colors[i] = meshMode == SOURCE_MESH  ? Point3D< float >(1.f, 0.6f, 0.6f) : Point3D< float >(0.6f, 0.6f, 1.f);

	sv.scale = 1.f;
	sv.translate = Point3D< float >(0.f, 0.f, 0.f);

	//Normalize the mesh
	for (int i = 0; i<sv.vertices.size(); i++) sv.vertices[i] = (sv.vertices[i] - visualizationCenter) / visualizationRadius;

	sv.selectedVertices.clear();
	sv.selectedVertices.insert(sv.selectedVertices.begin(), remeshingVertices[meshMode].begin(), remeshingVertices[meshMode].end());
	
	sv.loopPoints.clear();
	sv.loopPoints.resize(loopPoints[meshMode].size());
	for (int i = 0; i < loopPoints[meshMode].size(); i++) {
		sv.loopPoints[i].resize(loopPoints[meshMode][i].size());
		for (int j = 0; j < loopPoints[meshMode][i].size(); j++) {
			sv.loopPoints[i][j] = (loopPoints[meshMode][i][j] - visualizationCenter) / visualizationRadius;
		}
	}


	//sv.selectedVertex = -1;
	sv.mappedTriangle = -1;
	sv.transferedTriangle = -1;
}


void StitchingViewer::SetVisualizationMesh(const ColoredMesh & mesh){
	sv.selectedVertices.clear();
	sv.loopPoints.clear();
	sv.vertices.clear();
	sv.triangles.clear();
	sv.colors.clear();

	sv.vertices.resize(mesh.vertices.size());
	for (int i = 0; i < sv.vertices.size(); i++)sv.vertices[i] = Point3D<float>(mesh.vertices[i]);
	sv.triangles = mesh.triangles;
	sv.colors.resize(3 * sv.triangles.size());
	for (int i = 0; i<3 * sv.triangles.size(); i++) sv.colors[i] = mesh.colors[i];

	sv.scale = 1.f;
	sv.translate = Point3D< float >(0.f, 0.f, 0.f);

	//Normalize the mesh
	for (int i = 0; i<sv.vertices.size(); i++) sv.vertices[i] = (sv.vertices[i] - visualizationCenter) / visualizationRadius;
}

void StitchingViewer::ClippingAndStitching(){
	Timer t_Clipping;
	//roiVoxels = _addedVoxels;
	roinBoundaryFaceOrientation = refMap[0].boundaryFaceOrientation;

	std::unordered_set<int> trianglesToRemove[2];
	std::vector<int> originalToClippedTriangleMap[2];

	for (int i = 0; i < 2; i++) {

	#if USE_CONSTRAINED_DELANAUY_CLIPPING
	//Mod Nov 16
	//Augment vertex sample with the mesh vertices
	for (int tIndex = 0; tIndex < inputMesh[i].triangles.size(); tIndex++) {
		for (int k = 0; k < 3; k++) {
			unsigned long long vIndex = inputMesh[i].triangles[tIndex][k];
			Point3D<double> baricentric3D;
			baricentric3D[k] = 1.0;
			SamplePoint sample;
			sample.tIndex = tIndex;
			sample.baricentricCoord = Point2D<double>(baricentric3D[1], baricentric3D[2]);
			_vertexSample[i][vIndex] = sample;
		}
	}

	std::unordered_map<int, TriangleMap> triangleMaps;
	InitalizeTriangleMaps(gridDepth, inputMesh[i], triangleMaps, faceMaps[i], adjacentVertexMap[i], refMap[i].boundaryFaceOrientation);
	std::vector<unsigned long long > subdivisionTriangles;
	ProcessTriangleMap(inputMesh[i], triangleMaps, _vertexSample[i], subdivisionTriangles);
	MeshSubsetClippling(inputMesh[i], triangleMaps, _vertexSample[i], subdivisionTriangles, clippedMesh[i], clippedMeshVertexIndex[i], originalToClippedTriangleMap[i]);

	#else
	std::unordered_set<int> boundaryTriangles;
	IdentifyBoundaryAdjacentTriangles(inputMesh[i], refMap[i].boundaryFaceOrientation, faceMaps[i], adjacentVertexMap[i], boundaryTriangles);
	if (1) {
		char outputName[256];
		sprintf(outputName, "Scaled-%02d.ply", i);
		WriteSimpleMesh(inputMesh[i], outputName);
	}

	std::unordered_set<int>  trianglesToClip;
	MeshSubsetClippling(inputMesh[i], gridDepth, clippedMesh[i], clippedMeshVertexIndex[i], boundaryTriangles, originalToClippedTriangleMap[i]);
	#endif

	if (0) {
		char outputName[256];
		sprintf(outputName, "Clipped-%02d.ply", i);
		SimpleMesh scaledClipped = clippedMesh[i];
		for (int v = 0; v < scaledClipped.vertices.size(); v++) scaledClipped.vertices[v] = scaledClipped.vertices[v] * gridScaling + gridTranslation;
		WriteSimpleMesh(scaledClipped, outputName);
	}

	#if 0
	std::unordered_set<unsigned long long> boundaryEdges;
	IdentifyBoundaryEdges(refMap[i].boundaryFaceOrientation, faceMaps[i], adjacentVertexMap[i], clippedMeshVertexIndex[i], boundaryEdges, i == 0);

	BoundaryConstrainedPropagation(clippedMesh[i], boundaryEdges, trianglesToRemove[i]);
	#else
	VoxelBasedSegmentation(gridRes, clippedMesh[i], addedVoxels, trianglesToRemove[i], i == 1);
	#endif

	SetMeshEdgeKeyTriangleMap(clippedEdgeTriangleMap[i], clippedMesh[i]);


	if (0) {
		ColoredMesh coloredMesh;
		coloredMesh.vertices.resize(clippedMesh[i].triangles.size() * 3);
		coloredMesh.colors.resize(clippedMesh[i].triangles.size() * 3);
		coloredMesh.triangles.resize(clippedMesh[i].triangles.size());
		for (int j = 0; j < clippedMesh[i].triangles.size(); j++) {
			coloredMesh.triangles[j] = TriangleIndex(3 * j, 3 * j + 1, 3 * j + 2);
			for (int k = 0; k < 3; k++) coloredMesh.vertices[3 * j + k] = clippedMesh[i].vertices[clippedMesh[i].triangles[j][k]] * gridScaling + gridTranslation;
			if (trianglesToRemove[i].find(j) != trianglesToRemove[i].end()) {
				for (int k = 0; k < 3; k++) coloredMesh.colors[3 * j + k] = Point3D<double>(0.0, 0.0, 0.0);
			}
			else {
				for (int k = 0; k < 3; k++) coloredMesh.colors[3 * j + k] = i == 0 ? Point3D<double>(255.0, 0.0, 0.0) : Point3D<double>(0.0, 255.0, 0.0);
			}
		}
		char outputName[256];
		sprintf(outputName, "Segmentation-%02d.ply", i);
		WriteColoredMesh(coloredMesh, outputName);
	}
	}

	printf("Clipping time %.2f(s) \n", t_Clipping.elapsed());

	//SimpleMesh stitchedMesh;

	stitchedMesh.vertices.clear();
	stitchedMesh.triangles.clear();
	stitchedMesh.normals.clear();

	std::vector<int> clippedToStitchedTriangleMap[2];

	Timer t_Stitching;

	Stitching(faceMaps[0], adjacentVertexMap[0], clippedMeshVertexIndex[0], clippedEdgeTriangleMap[0], trianglesToRemove[0], clippedToStitchedTriangleMap[0], clippedMesh[0],
	faceMaps[1], adjacentVertexMap[1], clippedMeshVertexIndex[1], clippedEdgeTriangleMap[1], trianglesToRemove[1], clippedToStitchedTriangleMap[1], clippedMesh[1],
	roinBoundaryFaceOrientation, stitchedMesh);

	printf("Stitching time %.2f(s) \n", t_Stitching.elapsed());

	if (0) {
		SimpleMesh scaledStitched = stitchedMesh;
		for (int v = 0; v < scaledStitched.vertices.size(); v++) scaledStitched.vertices[v] = scaledStitched.vertices[v] * gridScaling + gridTranslation;
		WriteSimpleMesh(scaledStitched, "Stitched.ply");
	}
}


void StitchingViewer::SetForwardAndBackwardMaps(){
	forwardMap[0].clear();
	forwardMap[1].clear();
	backwardMap[0].clear();
	backwardMap[1].clear();

	for (int i = 0; i < 2; i++) {
		forwardMap[i].resize(inputMesh[i].triangles.size(), -1);
		int matchedFaces = 0;
#if 0
		for (int j = 0; j < originalToClippedTriangleMap[i].size(); j++) {
			int clippedIndex = originalToClippedTriangleMap[i][j];
			if (clippedIndex != -1) {
				int stitchedIndex = clippedToStitchedTriangleMap[i][clippedIndex];
				if (stitchedIndex != -1) {
					forwardMap[i][j] = stitchedIndex;
					matchedFaces++;
				}
			}
		}
#else
		std::set<FaceCoords, FaceCoordsComparison> initalFaces;
		for (int j = 0; j < inputMesh[i].triangles.size(); j++) {
			double coords[9] = { inputMesh[i].vertices[inputMesh[i].triangles[j][0]][0], inputMesh[i].vertices[inputMesh[i].triangles[j][0]][1], inputMesh[i].vertices[inputMesh[i].triangles[j][0]][2],
				inputMesh[i].vertices[inputMesh[i].triangles[j][1]][0], inputMesh[i].vertices[inputMesh[i].triangles[j][1]][1], inputMesh[i].vertices[inputMesh[i].triangles[j][1]][2],
				inputMesh[i].vertices[inputMesh[i].triangles[j][2]][0], inputMesh[i].vertices[inputMesh[i].triangles[j][2]][1], inputMesh[i].vertices[inputMesh[i].triangles[j][2]][2] };
			initalFaces.insert(FaceCoords(coords, j));
		}

		for (int j = 0; j < stitchedMesh.triangles.size(); j++) {
			double coords[9] = { stitchedMesh.vertices[stitchedMesh.triangles[j][0]][0], stitchedMesh.vertices[stitchedMesh.triangles[j][0]][1], stitchedMesh.vertices[stitchedMesh.triangles[j][0]][2],
				stitchedMesh.vertices[stitchedMesh.triangles[j][1]][0], stitchedMesh.vertices[stitchedMesh.triangles[j][1]][1], stitchedMesh.vertices[stitchedMesh.triangles[j][1]][2],
				stitchedMesh.vertices[stitchedMesh.triangles[j][2]][0], stitchedMesh.vertices[stitchedMesh.triangles[j][2]][1], stitchedMesh.vertices[stitchedMesh.triangles[j][2]][2] };
			FaceCoords simplifiedFace(coords, j);
			if (initalFaces.find(simplifiedFace) != initalFaces.end()) {
				FaceCoords matchedFace = *initalFaces.find(simplifiedFace);
				int matchedFaceIndex = matchedFace.index;
				if (forwardMap[i][matchedFaceIndex] != -1) {
					printf("ERROR: Non injective mapping!\n");
				}
				else {
					forwardMap[i][matchedFaceIndex] = j;
					matchedFaces++;
				}
			}
		}
#endif
		printf("Matched faces %d of %d \n", matchedFaces, inputMesh[i].triangles.size());

		if (0) {
			char vectorName[256];
			sprintf(vectorName, "ForwardMap-%02d.vec", i);
			WriteVector(forwardMap[i], vectorName);
		}

		backwardMap[i].resize(stitchedMesh.triangles.size(), -1);
		for (int j = 0; j < forwardMap[i].size(); j++) {
			int stitchedIndex = forwardMap[i][j];
			if (stitchedIndex != -1) {
				backwardMap[i][stitchedIndex] = j;
			}
		}

		if (0) {

			char vectorName[256];
			sprintf(vectorName, "BackwardMap-%02d.vec", i);
			WriteVector(backwardMap[i], vectorName);

			ColoredMesh coloredMesh;
			coloredMesh.vertices.resize(stitchedMesh.triangles.size() * 3);
			coloredMesh.colors.resize(stitchedMesh.triangles.size() * 3);
			coloredMesh.triangles.resize(stitchedMesh.triangles.size());
			for (int j = 0; j < stitchedMesh.triangles.size(); j++) {
				coloredMesh.triangles[j] = TriangleIndex(3 * j, 3 * j + 1, 3 * j + 2);
				for (int k = 0; k < 3; k++) coloredMesh.vertices[3 * j + k] = stitchedMesh.vertices[stitchedMesh.triangles[j][k]] * gridScaling + gridTranslation;
				if (backwardMap[i][j] != -1) {
					for (int k = 0; k < 3; k++) coloredMesh.colors[3 * j + k] = Point3D<double>(255.0, 0.0, 0.0);
				}
				else {
					for (int k = 0; k < 3; k++) coloredMesh.colors[3 * j + k] = Point3D<double>(0.0, 0.0, 0.0);
				}
			}
			char outputName[256];
			sprintf(outputName, "Preserved-%02d.ply", i);
			WriteColoredMesh(coloredMesh, outputName);
		}
	}
}


//Smooth non fixed vertices
void StitchingViewer::SmoothNonFixedVertices(){

	std::unordered_set<int> fixedVertices;
	for (int t = 0; t < stitchedMesh.triangles.size(); t++) {
		if (backwardMap[0][t] != -1 || backwardMap[1][t] != -1) {
			for (int k = 0; k < 3; k++) {
				fixedVertices.insert(stitchedMesh.triangles[t][k]);
			}
		}
	}

	std::unordered_map<int, std::unordered_set<int>> nonFixedVerticesNeighbours;
	for (int t = 0; t < stitchedMesh.triangles.size(); t++) {
		for (int k = 0; k < 3; k++) {
			if (fixedVertices.find(stitchedMesh.triangles[t][k]) == fixedVertices.end()) {
				nonFixedVerticesNeighbours[stitchedMesh.triangles[t][k]].insert(stitchedMesh.triangles[t][(k + 1) % 3]);
				nonFixedVerticesNeighbours[stitchedMesh.triangles[t][k]].insert(stitchedMesh.triangles[t][(k + 2) % 3]);
			}
		}
	}
	for (int smoothIter = 0; smoothIter < 10; smoothIter++) {
		std::unordered_map<int, Point3D<double>> newVertexPos;
		for (auto vertexIter = nonFixedVerticesNeighbours.begin(); vertexIter != nonFixedVerticesNeighbours.end(); vertexIter++) {
			int vIndex = (*vertexIter).first;
			std::unordered_set<int> & neighbours = (*vertexIter).second;
			Point3D<double> averagePos(0.0, 0.0, 0.0);
			int nCount = 0;
			for (auto neighbourIter = neighbours.begin(); neighbourIter != neighbours.end(); neighbourIter++) {
				int nIndex = *neighbourIter;
				averagePos += stitchedMesh.vertices[nIndex];
				nCount++;
			}
			newVertexPos[vIndex] = (averagePos / nCount);
		}
		for (auto vertexIter = newVertexPos.begin(); vertexIter != newVertexPos.end(); vertexIter++) {
			stitchedMesh.vertices[(*vertexIter).first] = (*vertexIter).second;
		}
	}
}

void StitchingViewer::ExportMesh(const char * outputName, const int exportMode) {
	
	if (exportMode == STITCHED_SIMPLE_MESH) {
		SimpleMesh rescaledStitched;
		rescaledStitched.vertices.resize(stitchedMesh.vertices.size());
		rescaledStitched.triangles = stitchedMesh.triangles;

		for (int v = 0; v < stitchedMesh.vertices.size(); v++) {
			rescaledStitched.vertices[v] = stitchedMesh.vertices[v] * gridScaling + gridTranslation;
		}
		WriteSimpleMesh(rescaledStitched, outputName);
	}
	else if (exportMode == STITCHED_SEGMENTED_MESH) {
		ColoredMesh coloredMesh;
		coloredMesh.vertices.resize(stitchedMesh.triangles.size() * 3);
		coloredMesh.colors.resize(stitchedMesh.triangles.size() * 3);
		coloredMesh.triangles.resize(stitchedMesh.triangles.size());
		for (int j = 0; j < stitchedMesh.triangles.size(); j++) {
			coloredMesh.triangles[j] = TriangleIndex(3 * j, 3 * j + 1, 3 * j + 2);
			for (int k = 0; k < 3; k++) coloredMesh.vertices[3 * j + k] = stitchedMesh.vertices[stitchedMesh.triangles[j][k]] * gridScaling + gridTranslation;
			if (backwardMap[0][j] != -1) {
				for (int k = 0; k < 3; k++) coloredMesh.colors[3 * j + k] = Point3D<double>(255.0, 0.0, 0.0);
			}
			else if (backwardMap[1][j] != -1) {
				for (int k = 0; k < 3; k++) coloredMesh.colors[3 * j + k] = Point3D<double>(0.0, 255.0, 0.0);
			}
			else {
				for (int k = 0; k < 3; k++) coloredMesh.colors[3 * j + k] = Point3D<double>(0.0, 0.0, 255.0);
			}
		}
		WriteColoredMesh(coloredMesh, outputName);
	}
}

int StitchingViewer::Process(const char * outputMeshName){
	
	// (1) Pick remeshing vertices using the default parameters
	UpdateRemeshingVertices();

	// (2) Mark voxel containing the remeshing vertices as seeds.
	UpdateVoxelSeeds();

	if (voxelSeeds.size()){
		
		// (3) Grown seeds
		PerformRegionGrowing();

		// (4) Clip and sticth
		ClippingAndStitching();

		// (5) Define forward and backward maps
		SetForwardAndBackwardMaps();
		
		// (6) Smooth non fixed vertices
		SmoothNonFixedVertices();

		// (7) Export mesh
		ExportMesh(outputMeshName, STITCHED_SIMPLE_MESH);

		if (1){
			char labeledMeshName[256];
			sprintf(labeledMeshName, "SegmentedMesh.ply");
			ExportMesh(labeledMeshName, STITCHED_SEGMENTED_MESH);
		}
	}
	else{
		printf("Nothing to remesh! \n");
	}
	return 1;
}
void StitchingViewer::UpdateRemeshingVertices(){

	selectionMode = AUTOMATIC_SELECTION;
	sprintf(sv.info[4], "Automatic vertex selection");
	sprintf(sv.info[3], "Max Hausdorff Distance(mm)  :  %f", maxHaussdorffDistance);
	sprintf(sv.info[2], "Min Visible Hausdorff Distance(mm)  :  %f", minVisibleHaussdorffDistance);
	sprintf(sv.info[1], "Max Visible Hausdorff Distance(mm)  :  %f", maxVisibleHaussdorffDistance);
	sprintf(sv.info[0], "Max Mean Curvature (1/mm)  :  %f", maxMeanCurvature);
	
	double _maxHaussdorffDistance = maxHaussdorffDistance/gridScaling;
	double _maxVisibleHaussdorffDistance = maxVisibleHaussdorffDistance/gridScaling;
	double _minVisibleHaussdorffDistance = minVisibleHaussdorffDistance/gridScaling;
	double _maxMeanCurvature = maxMeanCurvature*gridScaling;

	if(0)printf("Normalized max mean curvature %f \n", _maxMeanCurvature);
	if(0)printf("Grid scaling %f \n", gridScaling);
	//Update remeshing vertices
	remeshingVertices[0].clear();
	remeshingVertices[1].clear();

	Timer t_seedSelection;

	for (int i = 0; i < 2; i++){
		UpdateNormals(inputMesh[i]);
		SmoothSignal(FEM[i], normalSolver[i], inputMesh[i].normals, true);
	}


	//Energy = h(1 + p_0*v + p_1*v*m)
	double p_0 = (_maxHaussdorffDistance / _maxVisibleHaussdorffDistance) - 1.0;
	double p_1 = (_maxHaussdorffDistance / _minVisibleHaussdorffDistance) - (_maxHaussdorffDistance / _maxVisibleHaussdorffDistance);

	for (int i = 0; i < 2; i++) {

		std::vector<double>meanCurvature(inputMesh[i].vertices.size());

		{//Compute Mean Curvature
			double meshArea = 0.0;
			for (int t = 0; t < inputMesh[i].triangles.size(); t++) meshArea += Point3D<double>::Length(Point3D<double>::CrossProduct(inputMesh[i].vertices[inputMesh[i].triangles[t][1]] - inputMesh[i].vertices[inputMesh[i].triangles[t][0]], inputMesh[i].vertices[inputMesh[i].triangles[t][2]] - inputMesh[i].vertices[inputMesh[i].triangles[t][0]])) / 2.0;
			double mcFlowTimeStep = 0.001;
			Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>  * MCFlowSmoother = new Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>(FEM[i]->mass + FEM[i]->stiffness * mcFlowTimeStep);
			std::vector<Point3D<double>> vextexMCFlow = inputMesh[i].vertices;
			SmoothSignal(FEM[i], MCFlowSmoother, vextexMCFlow);

			if (0) {
				SimpleMesh mcFlowMesh;
				mcFlowMesh.vertices = vextexMCFlow;
				mcFlowMesh.triangles = inputMesh[i].triangles;
				char meanCurvatureFlowMeshName[256];
				sprintf(meanCurvatureFlowMeshName, "Mesh_MeanCurvatureFlow_%d.ply", i);
				WriteSimpleMesh(mcFlowMesh, meanCurvatureFlowMeshName);
			}

			double maxMC = 0;
			double minMC = DBL_MAX;
			for (int v = 0; v < inputMesh[i].vertices.size(); v++) {
				meanCurvature[v] = Point3D<double>::Length(vextexMCFlow[v] - inputMesh[i].vertices[v])/ mcFlowTimeStep;
				maxMC = std::max<double>(maxMC, meanCurvature[v]);
				minMC = std::min<double>(minMC, meanCurvature[v]);
			}

			for (int v = 0; v < inputMesh[i].vertices.size(); v++) meanCurvature[v] = std::min<double>(meanCurvature[v] / _maxMeanCurvature, 1.0);

			if(0)printf("Mesh %d . Max mean curvature %f. Min mean curvature %f \n", i, maxMC, minMC);
			if (0) {
				ColoredMesh coloredMeanCurvature;
				coloredMeanCurvature.vertices = inputMesh[i].vertices;
				coloredMeanCurvature.triangles = inputMesh[i].triangles;
				coloredMeanCurvature.colors.resize(inputMesh[i].vertices.size());
				for (int v = 0; v < inputMesh[i].vertices.size(); v++) {
					double factor = meanCurvature[v];
					coloredMeanCurvature.colors[v] = Point3D<double>(255, 0, 0)*factor + Point3D<double>(0, 0, 255)*(1.0 - factor);
				}
				char meanCurvatureMeshName[256];
				sprintf(meanCurvatureMeshName, "Mesh_MeanCurvature_%d.ply", i);
				WriteColoredMesh(coloredMeanCurvature, meanCurvatureMeshName);
			}
		}

		std::vector<float> ambientVisibilityScore;

		{//Compute Ambient Visibility

			double coneAngle = M_PI / 6.0;
			int numSamples = 16;
			AmbientVisibilityScore(inputMesh[i], 2000.f, ambientVisibilityScore, coneAngle, numSamples);
			if (0) {
				ColoredMesh coloredAmbientVisibility;
				coloredAmbientVisibility.vertices = inputMesh[i].vertices;
				coloredAmbientVisibility.triangles = inputMesh[i].triangles;
				coloredAmbientVisibility.colors.resize(inputMesh[i].vertices.size());
				for (int v = 0; v < inputMesh[i].vertices.size(); v++) {
					double factor = ambientVisibilityScore[v];
					coloredAmbientVisibility.colors[v] = Point3D<double>(255, 0, 0)*factor + Point3D<double>(0, 0, 255)*(1.0 - factor);
				}
				char AmbientVisibilityMeshName[256];
				sprintf(AmbientVisibilityMeshName, "Mesh_AmbientVisibility_%d.ply", i);
				WriteColoredMesh(coloredAmbientVisibility, AmbientVisibilityMeshName);
			}
		}

		std::vector<double> hausdorffDistance;
		double meanMin, hausdorff;

		HausdorffDistance(inputMesh[i].vertices, inputMesh[(i + 1) % 2].vertices, inputMesh[(i + 1) % 2].triangles, meanMin, hausdorff, hausdorffDistance);

		if (0) {
			ColoredMesh coloredHaussdorfDistance;
			coloredHaussdorfDistance.vertices = inputMesh[i].vertices;
			coloredHaussdorfDistance.triangles = inputMesh[i].triangles;
			coloredHaussdorfDistance.colors.resize(inputMesh[i].vertices.size());
			for (int v = 0; v < inputMesh[i].vertices.size(); v++) {
				double factor = 1.f - hausdorffDistance[v] / _maxVisibleHaussdorffDistance;
				coloredHaussdorfDistance.colors[v] = Point3D<double>(255, 0, 0)*factor + Point3D<double>(0, 0, 255)*(1.0 - factor);
			}
			char HaussdorfDistanceMeshName[256];
			sprintf(HaussdorfDistanceMeshName, "Mesh_HaussdorfDistance_%d.ply", i);
			WriteColoredMesh(coloredHaussdorfDistance, HaussdorfDistanceMeshName);
		}

		std::vector<bool> markedForRemeshBasedOnVisibility(inputMesh[i].vertices.size(), false);
		for (int v = 0; v < inputMesh[i].vertices.size(); v++) {
			if ((hausdorffDistance[v] * (1 + p_0*ambientVisibilityScore[v] + p_1*ambientVisibilityScore[v] * meanCurvature[v])) > _maxHaussdorffDistance) {
				markedForRemeshBasedOnVisibility[v] = true;
			}
		}
		if (0) RemoveSmallComponents(inputMesh[i], markedForRemeshBasedOnVisibility, 5);

		std::vector<bool> markedForRemeshing(inputMesh[i].vertices.size(), false);

		int remeshingCandidatesByVisibility = 0;
		for (int v = 0; v < inputMesh[i].vertices.size(); v++) {
			if (markedForRemeshBasedOnVisibility[v]) {
				remeshingCandidatesByVisibility++;
				markedForRemeshing[v] = true;
			}
		}

		if(0)printf("Marked by visibility %d \n", remeshingCandidatesByVisibility);


		for (int j = 0; j < inputMesh[i].vertices.size(); j++) if (markedForRemeshing[j]) remeshingVertices[i].insert(j);
		if (0) {
			ColoredMesh coloredMesh;
			coloredMesh.vertices = inputMesh[i].vertices;
			coloredMesh.colors.resize(inputMesh[i].vertices.size(), Point3D<double>(191.25, 191.25, 191.25));
			coloredMesh.triangles = inputMesh[i].triangles;
			for (auto j = remeshingVertices[i].begin(); j != remeshingVertices[i].end(); j++) {
				coloredMesh.colors[*j] = Point3D<double>(255, 201, 14);
			}
			char outputName[256];
			sprintf(outputName, "MarkedForRemesh-%02d.ply", i);
			WriteColoredMesh(coloredMesh, outputName);
		}
	}

	sv.selectedVertices.clear();
	sv.selectedVertices.insert(sv.selectedVertices.begin(), remeshingVertices[meshMode].begin(), remeshingVertices[meshMode].end());

	if(0)printf("Seed selection time %.2f(s)\n", t_seedSelection.elapsed());
}

void StitchingViewer::UpdateVoxelSeeds(){
	//Update voxel seeds
	voxelSeeds.clear();

	std::unordered_set<unsigned long> _voxelSeeds;
	for (int i = 0; i < 2; i++) {
		for (auto j = remeshingVertices[i].begin(); j != remeshingVertices[i].end(); j++) {
			//Point3D<double> gridPoint = ((inputMesh[i].vertices[*j] - gridTranslation) / gridScaling)*double(gridRes);
			Point3D<double> gridPoint = inputMesh[i].vertices[*j]*double(gridRes);
			unsigned long seed = SetPointVoxelKey(gridPoint);
			if (commonVoxelKeys.find(seed) == commonVoxelKeys.end()) {
				printf("Not found voxel seed! \n");
			}
			else _voxelSeeds.insert(seed);
		}
	}

	for (auto voxelIter = _voxelSeeds.begin(); voxelIter != _voxelSeeds.end(); voxelIter++) voxelSeeds.push_back(*voxelIter);

	if (0) {//Enable to use Wojtan seed selection
		_voxelSeeds.clear();
		voxelSeeds.clear();
		//Mark as seeds voxels with distinct faces
		for (auto voxelIter = commonVoxelKeys.begin(); voxelIter != commonVoxelKeys.end(); voxelIter++) {
			unsigned long currentVoxel = *voxelIter;
			unsigned long axialIndices[3];
			unsigned long _c;
			GetVoxelElementKeyIndices(currentVoxel, axialIndices, _c);
			for (unsigned long c = 0; c < 3; c++) {
				for (unsigned long offset = 0; offset < 2; offset++) {
					unsigned long neigbourFaceIndices[3] = { axialIndices[0], axialIndices[1], axialIndices[2] };
					if (offset == 1) {
						neigbourFaceIndices[c]++;
					}
					unsigned long neighbourFaceKey = SetVoxelElementKey(neigbourFaceIndices, c);
					if (!SimilarFaces(faceMaps[0][neighbourFaceKey], faceMaps[1][neighbourFaceKey])) {
						_voxelSeeds.insert(currentVoxel);
					}
				}
			}
		}
		for (auto voxelIter = _voxelSeeds.begin(); voxelIter != _voxelSeeds.end(); voxelIter++) voxelSeeds.push_back(*voxelIter);
	}
}

void StitchingViewer::PerformRegionGrowing() {
	addedVoxels.clear();
	refMap[0].clear();
	refMap[1].clear();
	if (voxelSeeds.size()) {
		std::unordered_map<unsigned long, FaceMap> _faceMaps[2];
		std::unordered_map<unsigned long, EdgeIntersections> _edgeIntersections[2];

		for (int i = 0; i < 2; i++){
			_faceMaps[i] = faceMaps[i];
			_vertexSample[i] = vertexSample[i];
			_edgeIntersections[i] = edgeIntersections[i];
		}
		std::unordered_set<unsigned long> _voxelKeys = commonVoxelKeys;
		std::unordered_set<unsigned long> _weakFaces = commonWeakFaces;
		Timer t_RegionGrowing;

		GrowVoxelRegion(gridRes,
			_faceMaps[0], _edgeIntersections[0], _vertexSample[0], refMap[0],
			_faceMaps[1], _edgeIntersections[1], _vertexSample[1], refMap[1],
			_voxelKeys, voxelSeeds, addedVoxels, INT_MAX, true, false, _weakFaces);

		printf("Region growing time %.2f(s) \n", t_RegionGrowing.elapsed());

		for (int i = 0; i < 2; i++){
			loopPoints[i].clear();
			SetLoopPoints(refMap[i].forwardMap, _vertexSample[i], inputMesh[i], loopPoints[i]);
		}
	}
}

void StitchingViewer::StitchingCallBack(Visualization*, const char* prompt) {

	// (4) Clip and sticth
	ClippingAndStitching();

	// (5) Define forward and backward maps
	SetForwardAndBackwardMaps();

	coloredStitchedMesh.vertices.resize(stitchedMesh.triangles.size() * 3);
	coloredStitchedMesh.colors.resize(stitchedMesh.triangles.size() * 3);
	coloredStitchedMesh.triangles.resize(stitchedMesh.triangles.size());
	for (int j = 0; j < stitchedMesh.triangles.size(); j++){
		coloredStitchedMesh.triangles[j] = TriangleIndex(3 * j, 3 * j + 1, 3 * j + 2);
		for (int k = 0; k < 3; k++) coloredStitchedMesh.vertices[3 * j + k] = stitchedMesh.vertices[stitchedMesh.triangles[j][k]];
		if (backwardMap[0][j] != -1){
			for (int k = 0; k < 3; k++) coloredStitchedMesh.colors[3 * j + k] = Point3D<double>(1.0, 0.6, 0.6);
		}
		else if (backwardMap[1][j] != -1){
			for (int k = 0; k < 3; k++) coloredStitchedMesh.colors[3 * j + k] = Point3D<double>(0.6, 0.6, 1.0);
		}
		else{
			for (int k = 0; k < 3; k++) coloredStitchedMesh.colors[3 * j + k] = Point3D<double>(0.6, 1.0, 0.6);
		}
	}

	sv.showSelectedVertices = false;

	SetVisualizationMesh(coloredStitchedMesh);
	meshMode = STITCHED_MESH;
	meshCount = 3;
	sv.updateMesh();

	currentSystemStage = CLIPPING_STAGE;

	glutPostRedisplay();
}

void StitchingViewer::PerformRegionGrowingCallBack(Visualization*, const char* prompt){
	UpdateVoxelSeeds();
	if (voxelSeeds.size()){
		PerformRegionGrowing();
		sv.voxels.clear();
		sv.voxels.insert(sv.voxels.begin(), addedVoxels.begin(), addedVoxels.end());

		sv.loopPoints.clear();
		sv.loopPoints.resize(loopPoints[meshMode].size());
		for (int i = 0; i < loopPoints[meshMode].size(); i++) {
			sv.loopPoints[i].resize(loopPoints[meshMode][i].size());
			for (int j = 0; j < loopPoints[meshMode][i].size(); j++) {
				sv.loopPoints[i][j] = (loopPoints[meshMode][i][j] - visualizationCenter) / visualizationRadius;
			}
		}
	}
	glutPostRedisplay();
}


void StitchingViewer::ExportMeshCallBack(Visualization*, const char* prompt) {
	glutSetCursor(_busyCursor);
	ExportMesh(prompt, STITCHED_SIMPLE_MESH);
	glutSetCursor(GLUT_CURSOR_INHERIT);
}

void StitchingViewer::MaxHausdorffDistanceCallBack(Visualization*, const char* prompt) {
	glutSetCursor(_busyCursor);
	maxHaussdorffDistance = atof(prompt);
	UpdateRemeshingVertices();
	glutSetCursor(GLUT_CURSOR_INHERIT);
}

void StitchingViewer::MaxVisibleHausdorffDistanceCallBack(Visualization*, const char* prompt) {
	glutSetCursor(_busyCursor);
	maxVisibleHaussdorffDistance = atof(prompt);
	UpdateRemeshingVertices();
	glutSetCursor(GLUT_CURSOR_INHERIT);
}

void StitchingViewer::MinVisibleHausdorffDistanceCallBack(Visualization*, const char* prompt) {
	glutSetCursor(_busyCursor);
	minVisibleHaussdorffDistance = atof(prompt);
	UpdateRemeshingVertices();
	glutSetCursor(GLUT_CURSOR_INHERIT);
}

void StitchingViewer::MaxMeanCurvatureCallBack(Visualization*, const char* prompt) {
	glutSetCursor(_busyCursor);
	maxMeanCurvature = atof(prompt);
	UpdateRemeshingVertices();
	glutSetCursor(GLUT_CURSOR_INHERIT);
}


void StitchingViewer::Idle( void ){ sv.Idle();}
void StitchingViewer::KeyboardFunc( unsigned char key , int x , int y ){
	sv.KeyboardFunc( key , x , y ); 
}
void StitchingViewer::SpecialFunc( int key , int x , int y ){
	switch (key)
	{
	case KEY_LEFTARROW: referenceVoxel[0] = std::max<unsigned long>(referenceVoxel[0] - 1, 0); break;
	case KEY_RIGHTARROW: referenceVoxel[0] = std::min<unsigned long>(referenceVoxel[0] + 1, gridRes -1); break;
	case KEY_DOWNARROW: referenceVoxel[1] = std::max<unsigned long>(referenceVoxel[1] - 1, 0); break;
	case KEY_UPARROW: referenceVoxel[1] = std::min<unsigned long>(referenceVoxel[1] + 1, gridRes - 1); break;
	case KEY_PGUP: referenceVoxel[2] = std::max<unsigned long>(referenceVoxel[2] - 1, 0); break;
	case KEY_PGDN: referenceVoxel[2] = std::min<unsigned long>(referenceVoxel[2] + 1, gridRes - 1); break;
	}
	glutPostRedisplay();
	//sv.SpecialFunc( key , x ,  y ); 
}
void StitchingViewer::Display( void ){ sv.Display(); }
void StitchingViewer::Reshape( int w , int h ){ sv.Reshape( w , h ); }
void StitchingViewer::MouseFunc( int button , int state , int x , int y )
{
	if (glutGetModifiers() & GLUT_ACTIVE_SHIFT && state == GLUT_DOWN && currentSystemStage == SELECTION_STAGE){
		int selectedVertex = sv.selectVertex(x, y);
		if (selectedVertex != -1){
			if (remeshingVertices[meshMode].find(selectedVertex) == remeshingVertices[meshMode].end()) {
				remeshingVertices[meshMode].insert(selectedVertex);
			}
			else {
				remeshingVertices[meshMode].erase(selectedVertex);
			}
			sv.selectedVertices.clear();
			sv.selectedVertices.insert(sv.selectedVertices.begin(), remeshingVertices[meshMode].begin(), remeshingVertices[meshMode].end());
		}
	}
	else sv.MouseFunc( button , state , x , y );

	glutPostRedisplay();
}
void StitchingViewer::MotionFunc( int x , int y ){ sv.MotionFunc( x , y ); }