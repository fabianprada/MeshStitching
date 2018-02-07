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


#ifndef RAY_TRACER_INCLUDED
#define RAY_TRACER_INCLUDED

#include <embree2/rtcore.h>
#include <embree2/rtcore_ray.h>
#include "SimpleMesh.h"

struct RTVertex   { float x, y, z, a; };
struct RTTriangle { int v0, v1, v2; };

class IntersectionData{
public:
	IntersectionData(){
		tId = -1;
		u = v = -1.f;
		time = FLT_MAX;
	}
	int tId;
	float u;
	float v;
	float time;
};

class MeshRayIntersection{
public:
	void Init(const std::vector<Point3D<double>> & vertices, const std::vector< TriangleIndex > & triangles);
	void Terminate();
	bool IntersectRay(const Point3D<double> & origin, const Point3D<double> & direction, float maxDistance, IntersectionData & intersection, const SimpleMesh & targetMesh, const Point3D<double> & expectedNormal, float angleCutOff, bool verbose = false) const;
	bool IntersectLine(const Point3D<double> & origin, const Point3D<double> & direction, float maxDistance, IntersectionData & intersection, const SimpleMesh & targetMesh, const Point3D<double> & expectedNormal, float angleCutOff, bool verbose = false) const;
	bool IntersectRay(const Point3D<double> & origin, const Point3D<double> & direction, float maxDistance, Point3D<double> & intersectionPosition, Point3D<double> & intersectionNormal, const SimpleMesh & targetMesh, bool verbose, bool smoothedNormal) const;
	RTCDevice device;
	RTCScene scene;
	unsigned meshID;
};

void MeshRayIntersection::Init(const std::vector<Point3D<double>> & vertices, const std::vector< TriangleIndex > & triangles){

	device = rtcNewDevice(NULL);
	scene = rtcDeviceNewScene(device, RTC_SCENE_STATIC, RTC_INTERSECT1);
	meshID = rtcNewTriangleMesh(scene, RTC_GEOMETRY_STATIC, triangles.size(), vertices.size());

	RTVertex* RTvertices = (RTVertex*)rtcMapBuffer(scene, meshID, RTC_VERTEX_BUFFER);
	for (int i = 0; i < vertices.size(); i++) RTvertices[i].x = vertices[i][0], RTvertices[i].y = vertices[i][1], RTvertices[i].z = vertices[i][2];
	rtcUnmapBuffer(scene, meshID, RTC_VERTEX_BUFFER);

	RTTriangle* RTtriangles = (RTTriangle*)rtcMapBuffer(scene, meshID, RTC_INDEX_BUFFER);
	for (int i = 0; i < triangles.size(); i++) RTtriangles[i].v0 = triangles[i][0], RTtriangles[i].v1 = triangles[i][1], RTtriangles[i].v2 = triangles[i][2];
	rtcUnmapBuffer(scene, meshID, RTC_INDEX_BUFFER);

	rtcCommit(scene);
}

void MeshRayIntersection::Terminate(){
	 rtcDeleteScene(scene);
	 rtcDeleteDevice(device);
 }

bool MeshRayIntersection::IntersectRay(const Point3D<double> & origin, const Point3D<double> & direction, float maxDistance, Point3D<double> & intersectionPosition, Point3D<double> & intersectionNormal,const SimpleMesh & targetMesh, bool verbose, bool smoothedNormal) const{
	RTCRay ray;
	ray.org[0] = origin[0];
	ray.org[1] = origin[1];
	ray.org[2] = origin[2];

	ray.dir[0] = direction[0];
	ray.dir[1] = direction[1];
	ray.dir[2] = direction[2];

	ray.tnear = 0.f;
	ray.tfar = maxDistance;
	ray.geomID = RTC_INVALID_GEOMETRY_ID;
	ray.primID = RTC_INVALID_GEOMETRY_ID;
	ray.instID = RTC_INVALID_GEOMETRY_ID;
	ray.mask = 0xFFFFFFFF;
	ray.time = 0.f;
	rtcIntersect(scene, ray);

	if (ray.geomID == meshID){
		intersectionPosition = targetMesh.vertices[targetMesh.triangles[ray.primID][0]] * (1.f - ray.u - ray.v) + targetMesh.vertices[targetMesh.triangles[ray.primID][1]] * ray.u + targetMesh.vertices[targetMesh.triangles[ray.primID][2]] * ray.v;
		if (smoothedNormal)	intersectionNormal = targetMesh.normals[targetMesh.triangles[ray.primID][0]] * (1.f - ray.u - ray.v) + targetMesh.normals[targetMesh.triangles[ray.primID][1]] * ray.u + targetMesh.normals[targetMesh.triangles[ray.primID][2]] * ray.v;
		else intersectionNormal = Point3D< double >::CrossProduct(targetMesh.vertices[targetMesh.triangles[ray.primID][1]] - targetMesh.vertices[targetMesh.triangles[ray.primID][0]], targetMesh.vertices[targetMesh.triangles[ray.primID][2]] - targetMesh.vertices[targetMesh.triangles[ray.primID][0]]);
		intersectionNormal /= Point3D<double>::Length(intersectionNormal);
		return true;
	}
	else{
		if (verbose) printf("No Ray Intersection \n");
		return false;
	}
}

bool MeshRayIntersection::IntersectRay(const Point3D<double> & origin, const Point3D<double> & direction, float maxDistance, IntersectionData & intersection, const SimpleMesh & targetMesh, const Point3D<double> & expectedNormal, float angleCutOff, bool verbose) const{
	RTCRay ray;
	ray.org[0] = origin[0];
	ray.org[1] = origin[1];
	ray.org[2] = origin[2];

	ray.dir[0] = direction[0];
	ray.dir[1] = direction[1];
	ray.dir[2] = direction[2];

	ray.tnear = 0.f;
	ray.tfar = maxDistance;
	ray.geomID = RTC_INVALID_GEOMETRY_ID;
	ray.primID = RTC_INVALID_GEOMETRY_ID;
	ray.instID = RTC_INVALID_GEOMETRY_ID;
	ray.mask = 0xFFFFFFFF;
	ray.time = 0.f;
	rtcIntersect(scene, ray);

	intersection.tId = -1;
	intersection.u = intersection.v = -1.f;
	intersection.time = FLT_MAX;
	
	if (ray.geomID == meshID){
		Point3D<double> intersectionNormal = targetMesh.normals[targetMesh.triangles[ray.primID][0]] * (1.f - ray.u - ray.v) + targetMesh.normals[targetMesh.triangles[ray.primID][1]] * ray.u + targetMesh.normals[targetMesh.triangles[ray.primID][2]] * ray.v;
		intersectionNormal /= Point3D<double>::Length(intersectionNormal);
		if (Point3D<double>::Dot(expectedNormal, intersectionNormal) > angleCutOff){
			intersection.tId = ray.primID;
			intersection.u = ray.u;
			intersection.v = ray.v;
			Point3D<double> intersectionPosition = targetMesh.vertices[targetMesh.triangles[ray.primID][0]] * (1.f - ray.u - ray.v) + targetMesh.vertices[targetMesh.triangles[ray.primID][1]] * ray.u + targetMesh.vertices[targetMesh.triangles[ray.primID][2]] * ray.v;
			intersection.time = Point3D<double>::Length(intersectionPosition - origin);
			return true;
		}
		else{
			if (verbose){
				printf("Failed normal test!\n");
				printf("Value %f Threshold %f \n", Point3D<double>::Dot(expectedNormal, intersectionNormal), angleCutOff);
				printf("Expected %f %f %f \n", expectedNormal[0], expectedNormal[1], expectedNormal[2]);
				printf("Intersection %f %f %f \n", intersectionNormal[0], intersectionNormal[1], intersectionNormal[2]);
			}
			return false;
		}
	}
	else{
		if(verbose) printf("No Ray Intersection \n");
		return false;
	}
}

bool MeshRayIntersection::IntersectLine(const Point3D<double> & origin, const Point3D<double> & direction, float maxDistance, IntersectionData & intersection, const SimpleMesh & targetMesh, const Point3D<double> & expectedNormal, float angleCutOff, bool verbose) const{

	IntersectionData _intersections[2];
	bool _intesects[2];
	_intesects[0] = IntersectRay(origin, direction, maxDistance, _intersections[0], targetMesh, expectedNormal, angleCutOff, verbose);
	_intesects[1] = IntersectRay(origin, -direction, maxDistance, _intersections[1], targetMesh, expectedNormal, angleCutOff, verbose);
	if (!_intesects[0] && !_intesects[1]){
		intersection.tId = -1;
		intersection.time = intersection.u = intersection.v = -1.f;
		if (verbose) printf("No intersection! \n");
		return false;
	}
	else{
		if (_intesects[0] && _intesects[1]){
			if (verbose) printf("Two intersections at %f and %f! \n", _intersections[0].time, _intersections[1].time);
			if (_intersections[0].time < _intersections[1].time){
				intersection = _intersections[0];
				if (verbose) printf("Positive intersection is closer! \n");
			}
			else{
				intersection = _intersections[1];
				intersection.time *= -1.f;
				if (verbose) printf("Negative intersection is closer! \n");
			}
		}
		else if(_intesects[0]){
			intersection = _intersections[0];
			if (verbose) printf("Intersection in the positive direction only! \n");
		}
		else{
			intersection = _intersections[1];
			intersection.time *= -1.f;
			if (verbose) printf("Intersection in the negative direction only! \n");
		}
		return true;
	}
}

#endif //RAY_TRACER_INCLUDED