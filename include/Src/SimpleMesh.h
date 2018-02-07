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


#ifndef SIMPLE_MESH_INCLUDED
#define SIMPLE_MESH_INCLUDED

#include <Misha/Ply.h>
#include <Misha/Image.h>

class SimpleMesh{
public:
	std::vector<Point3D<double>> vertices;
	std::vector<Point3D<double>> normals;
	std::vector<TriangleIndex> triangles;
};

class ColoredMesh : public SimpleMesh{
public:
	std::vector<Point3D<double>> colors;
};


class TexturedMesh : public SimpleMesh{
public:
	std::vector<Point2D<double>> textureCoordinates;
	Image<Point3D<double>> texture;
};

void UpdateNormals(SimpleMesh & mesh){
	mesh.normals.clear();
	mesh.normals.resize(mesh.vertices.size(), Point3D<double>(0.0, 0.0, 0.0));
	for (int t = 0; t < mesh.triangles.size(); t++){
		Point3D<double> d01 = mesh.vertices[mesh.triangles[t][1]] - mesh.vertices[mesh.triangles[t][0]];
		Point3D<double> d02 = mesh.vertices[mesh.triangles[t][2]] - mesh.vertices[mesh.triangles[t][0]];
		Point3D<double> n = Point3D<double>::CrossProduct(d01,d02);
		for (int v = 0; v < 3; v++) mesh.normals[mesh.triangles[t][v]] += n;
	}

	for (int i = 0; i < mesh.normals.size(); i++){
		if (Point3D<double>::Length(mesh.normals[i])> 0) mesh.normals[i] /= Point3D<double>::Length(mesh.normals[i]);
	}
}

int ReadSimpleMesh(SimpleMesh & mesh, const char * fileName){
	mesh.vertices.clear();
	mesh.triangles.clear();
	int file_type;
	std::vector< PlyVertex< double > > ply_vertices;
	bool readFlags[PlyVertex< double >::ReadComponents];
	if (!PlyReadTriangles(fileName, ply_vertices, mesh.triangles, PlyVertex< double >::ReadProperties, readFlags, PlyVertex< double >::ReadComponents, file_type)) return 0;
	mesh.vertices.resize(ply_vertices.size());
	for (int i = 0; i < ply_vertices.size(); i++) mesh.vertices[i] = Point3D<double>(ply_vertices[i].point[0], ply_vertices[i].point[1], ply_vertices[i].point[2]);
	UpdateNormals(mesh);
	return 1;
}

void WriteSimpleMesh(SimpleMesh & mesh, const char * fileName){
	std::vector< PlyVertex< float > > ply_vertices(mesh.vertices.size());
	for (int i = 0; i<mesh.vertices.size(); i++) ply_vertices[i].point = Point3D<float>(mesh.vertices[i][0], mesh.vertices[i][1], mesh.vertices[i][2]);
	PlyWriteTriangles(fileName, ply_vertices, mesh.triangles, PlyVertex< float >::WriteProperties, PlyVertex< float >::WriteComponents, PLY_BINARY_NATIVE);
}

void WriteColoredMesh(ColoredMesh & mesh, const char * fileName){
	std::vector< PlyColorVertex< float > > ply_vertices(mesh.vertices.size());
	for (int i = 0; i<mesh.vertices.size(); i++) ply_vertices[i].point = Point3D<float>(mesh.vertices[i][0], mesh.vertices[i][1], mesh.vertices[i][2]), ply_vertices[i].color = Point3D<float>(mesh.colors[i][0], mesh.colors[i][1], mesh.colors[i][2]);
	PlyWriteTriangles(fileName, ply_vertices, mesh.triangles, PlyColorVertex< float >::WriteProperties, PlyColorVertex< float >::WriteComponents, PLY_BINARY_NATIVE);
}

void WriteTexturedMesh(TexturedMesh & mesh, const char * fileName, const char * atlasName = NULL)
{
	std::vector< PlyTexturedFace< float > > texturedTriangles;
	texturedTriangles.resize(mesh.triangles.size());
	for (int i = 0; i < mesh.triangles.size(); i++){
		texturedTriangles[i].resize(3);
		for (int j = 0; j < 3; j++){
			texturedTriangles[i][j] = mesh.triangles[i][j];
			texturedTriangles[i].texture(j) = Point2D<float>(mesh.textureCoordinates[3 * i + j][0], mesh.textureCoordinates[3 * i + j][1]);
		}
	}

	std::vector< PlyVertex< float > > vertices(mesh.vertices.size());
	for (int i = 0; i < mesh.vertices.size(); i++)vertices[i].point = Point3D<float>(mesh.vertices[i][0], mesh.vertices[i][1], mesh.vertices[i][2]);

	if (atlasName != NULL){
		char ** comments = new char *[1];
		char atlas_comment[256];
		sprintf(atlas_comment, "TextureFile %s", atlasName);
		comments[0] = atlas_comment;
		PlyWritePolygons(fileName, vertices, texturedTriangles, PlyVertex< float >::WriteProperties, PlyVertex< float >::WriteComponents, PlyTexturedFace< float >::WriteProperties, PlyTexturedFace< float >::WriteComponents, PLY_ASCII, comments, 1);
	}
	else{
		PlyWritePolygons(fileName, vertices, texturedTriangles, PlyVertex< float >::WriteProperties, PlyVertex< float >::WriteComponents, PlyTexturedFace< float >::WriteProperties, PlyTexturedFace< float >::WriteComponents, PLY_ASCII);
	}
}

int ReadTexturedMesh(TexturedMesh & mesh, const char * meshName, const char* atlasName){
	mesh.vertices.clear();
	mesh.triangles.clear();
	mesh.textureCoordinates.clear();
	int file_type;
	std::vector< PlyVertex< double > > ply_vertices;
	std::vector< PlyTexturedFace< double > > ply_faces;
	if (!PlyReadPolygons(meshName, ply_vertices, ply_faces, PlyVertex< double >::ReadProperties, NULL, PlyVertex< double >::ReadComponents, PlyTexturedFace< double >::ReadProperties, NULL, PlyTexturedFace< double >::ReadComponents, file_type)) return 0;
	
	mesh.vertices.resize(ply_vertices.size());
	for (int i = 0; i < ply_vertices.size(); i++) mesh.vertices[i] = ply_vertices[i].point;

	mesh.triangles.resize(ply_faces.size());
	mesh.textureCoordinates.resize(3*ply_faces.size());
	for (int i = 0; i < ply_faces.size(); i++){
		for (int j = 0; j < 3; j++){
			mesh.triangles[i][j] = ply_faces[i][j];
			mesh.textureCoordinates[3 * i + j] = ply_faces[i].texture(j);
		}
	}
	UpdateNormals(mesh);
	mesh.texture.read(atlasName);

	return 1;
}

template< class Data >
int MeshSubdivide(const std::vector<TriangleIndex> & inTriangles, const std::vector<Data> & inVertexData, std::vector<TriangleIndex> & outTriangles, Data** outVertexData)
{
#define EDGE_KEY( i1 , i2 ) ( (i1)>(i2) ? ( ( (long long) (i1) )<<32 ) | ( (long long) (i2) ) : ( ( (long long) (i2) )<<32 ) | ( (long long) (i1) ) )

	std::vector< Data > vertexData;
	int vCount = (int)inVertexData.size();
	vertexData.resize(vCount);
	for (int i = 0; i<vCount; i++) vertexData[i] = inVertexData[i];

	std::unordered_map< long long, int > vMap;
	std::vector< Data > _vertexData = vertexData;
	std::vector< TriangleIndex > _triangles;
	for (int i = 0; i<inTriangles.size(); i++)
	{
		long long keys[] = { EDGE_KEY(inTriangles[i][0], inTriangles[i][1]), EDGE_KEY(inTriangles[i][1], inTriangles[i][2]), EDGE_KEY(inTriangles[i][2], inTriangles[i][0]) };
		int eIndex[3];
		for (int j = 0; j<3; j++)
			if (vMap.find(keys[j]) == vMap.end()) vMap[keys[j]] = eIndex[j] = (int)_vertexData.size(), _vertexData.push_back((vertexData[inTriangles[i][j]] + vertexData[inTriangles[i][(j + 1) % 3]]) / 2.f);
			else eIndex[j] = vMap[keys[j]];
			_triangles.push_back(TriangleIndex(eIndex[0], eIndex[1], eIndex[2]));
			_triangles.push_back(TriangleIndex(inTriangles[i][0], eIndex[0], eIndex[2]));
			_triangles.push_back(TriangleIndex(inTriangles[i][1], eIndex[1], eIndex[0]));
			_triangles.push_back(TriangleIndex(inTriangles[i][2], eIndex[2], eIndex[1]));
	}

	*outVertexData = new Data[_vertexData.size()];
	for (int i = 0; i<_vertexData.size(); i++) (*outVertexData)[i] = _vertexData[i];
	outTriangles.resize(_triangles.size());
	for (int i = 0; i<_triangles.size(); i++) outTriangles[i] = _triangles[i];
	return (int)_vertexData.size();
#undef EDGE_KEY
}

void SimpleMeshSubdivide(SimpleMesh & mesh, int numIters){
	for (int i = 0; i < numIters; i++){
		std::vector<TriangleIndex> outTriangles;
		Point3D<double>* tVertices;
		int vCount = MeshSubdivide<Point3D<double>>(mesh.triangles, mesh.vertices, outTriangles, &tVertices);
		mesh.vertices.resize(vCount);
		for (int i = 0; i<vCount; i++) mesh.vertices[i] = tVertices[i];
		delete[] tVertices;
		mesh.triangles = outTriangles;
	}
}


#endif//SIMPLE_MESH_INCLUDED
