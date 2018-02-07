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


#undef ARRAY_DEBUG
 
#include <stdlib.h>
#include <GL/glew.h>
#include <GL/glut.h>
#include "StitchingViewer.h"

cmdLineParameter<  char* > Source("source");
cmdLineParameter<  char* > Target("target");
cmdLineParameter<  char* > Output("output");
cmdLineParameter<  int   > Depth("depth", 7);
cmdLineParameter<  int   > GeodesicMarkers("markers", 5);
cmdLineReadable DisableNormalOffset("noNormalOffset");
cmdLineReadable DisableJitter("noJitter");
cmdLineReadable* params[] =
{
	&Source, &Target, &Output, &Depth, &GeodesicMarkers, &DisableNormalOffset, &DisableJitter,
	NULL
};

void ShowUsage(const char* ex)
{
	printf("Usage %s:\n", ex);
	printf("\t --%s <source>\n", Source.name);
	printf("\t --%s <target>\n", Target.name);
	printf("\t --%s <output>\n", Output.name);
	printf("\t --%s <grid depth> [%d]\n", Depth.name,Depth.value);
	printf("\t --%s <geodesic markers> [%d]\n", GeodesicMarkers.name, GeodesicMarkers.value);
	printf("\t --%s <disable normal offset>\n", DisableNormalOffset.name);
	printf("\t --%s <disable jitter>\n", DisableJitter.name);
}

int main(int argc, char* argv[])
{
	cmdLineParse(argc - 1, argv + 1, params);
	if (!Source.set || !Target.set){
		ShowUsage(argv[0]); return EXIT_FAILURE;
	}

	if (!StitchingViewer::Init(Source.value, Target.value, Depth.value, GeodesicMarkers.value, DisableNormalOffset.set,!DisableJitter.set)) return 0;

	if (!Output.set){
		glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
		glutInitWindowSize(StitchingViewer::sv.screenWidth, StitchingViewer::sv.screenHeight);
		glutInit(&argc, argv);
		char windowName[1024];
		sprintf(windowName, "Stitching Viewer");
		glutCreateWindow(windowName);

		if (glewInit() != GLEW_OK) fprintf(stderr, "[ERROR] glewInit failed\n"), exit(0);
		glutIdleFunc(StitchingViewer::Idle);
		glutDisplayFunc(StitchingViewer::Display);
		glutReshapeFunc(StitchingViewer::Reshape);
		glutMouseFunc(StitchingViewer::MouseFunc);
		glutMotionFunc(StitchingViewer::MotionFunc);
		glutKeyboardFunc(StitchingViewer::KeyboardFunc);
		glutSpecialFunc(StitchingViewer::SpecialFunc);
		glutMainLoop();
	}
	else{
		return StitchingViewer::Process(Output.value);
	}
	return 0;
}
