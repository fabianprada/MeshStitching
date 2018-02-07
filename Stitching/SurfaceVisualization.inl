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


#include <Misha/Timer.h>
#include <Misha/Visualization.h>
#include <Misha/Camera.h> 
#include <Misha/Image.h>
#ifndef M_PI
#define M_PI		3.14159265358979323846
#endif // M_PI

struct SurfaceVisualization : public Visualization
{

	GLuint offscreen_depth_texture = 0;
	GLuint offscreen_color_texture = 0;
	GLuint offscreen_framebuffer_handle = 0;
	int offscreen_frame_width, offscreen_frame_height;
	void SetupOffScreenBuffer();
	void RenderOffScreenBuffer(Image<Point3D<float>> & image);
	
	Point3D< float > materialColor;

	bool showFaces;
	std::vector<unsigned long> faces;
	std::vector<unsigned long> edges;
	bool showVoxels;
	std::vector<unsigned long> voxels;
	bool showReferenceVoxel;
	unsigned long referenceVoxel;
	float integerGridScaling;
	Point3D<float> integerGridTranslation;

	std::vector< TriangleIndex > triangles;
	std::vector< Point3D< float > > vertices,oldVertices,colors;
	std::vector<std::vector< Point3D< float > >> loopPoints;
	std::vector< int > referenceVertices;
	std::vector< Point3D< float > > exponetialPoints;
	int mappedTriangle;
	Point2D< float > mappedBaricenter;
	int transferedTriangle;
	Point2D< float > transferedBaricenter;

	Camera camera;
	float zoom;
	Point3D< float > translate;
	float scale;
	GLuint vbo, ebo;
	GLfloat lightAmbient[4], lightDiffuse[4], lightSpecular[4], shapeSpecular[4], shapeSpecularShininess;
	int oldX, oldY, newX, newY;
	float imageZoom, imageOffset[2];
	bool rotating, scaling, panning;
	bool useLight, showEdges, showColor, hasColor;
	bool showLoopPoints;
	bool showMesh;
	//int selectedVertex;

	bool showSelectedVertices;
	std::vector<int> selectedVertices;
	int selectVertex(int x, int y);

	void DrawSphere(unsigned long voxelKey);
	void DrawSphere(Point3D<float> corner);
	void DrawVoxel(unsigned long voxelKey);
	void DrawVoxel(Point3D<float> corner);
	void DrawEdge(unsigned long voxelKey);
	void DrawEdge(Point3D<float> corner, unsigned long c);

	void DrawFace(unsigned long voxelKey);
	void DrawFace(Point3D<float> corner, unsigned long c);

	SurfaceVisualization(void);
	void updateMesh();
	bool select(int x, int y, Point3D< float >& p);
	void display(void);
	void mouseFunc(int button, int state, int x, int y);
	void motionFunc(int x, int y);

	static void		ScreenshotCallBack(Visualization* v, const char* prompt);
	static void		WriteSceneConfigurationCallBack(Visualization* v, const char* prompt);
	static void		ReadSceneConfigurationCallBack(Visualization* v, const char* prompt);
	static void     ToggleLightCallBack(Visualization* v, const char*){ ((SurfaceVisualization*)v)->useLight = !((SurfaceVisualization*)v)->useLight; }
	static void     ToggleEdgesCallBack(Visualization* v, const char*){ ((SurfaceVisualization*)v)->showEdges = !((SurfaceVisualization*)v)->showEdges; }
	static void     ToggleColorCallBack(Visualization* v, const char*){ ((SurfaceVisualization*)v)->showColor = !((SurfaceVisualization*)v)->showColor; }
	static void     ToggleMeshCallBack(Visualization* v, const char*){ ((SurfaceVisualization*)v)->showMesh = !((SurfaceVisualization*)v)->showMesh; }
	static void     ToggleVoxelsCallBack(Visualization* v, const char*){ ((SurfaceVisualization*)v)->showVoxels = !((SurfaceVisualization*)v)->showVoxels; }
	static void     ToggleFacesCallBack(Visualization* v, const char*){ ((SurfaceVisualization*)v)->showFaces = !((SurfaceVisualization*)v)->showFaces; }
	bool setPosition(int x, int y, Point3D< double >& p);
	bool setPosition(int x, int y, Point3D< float >& p);
};


SurfaceVisualization::SurfaceVisualization(void)
{
	zoom = 1.05f;
	materialColor[0] = materialColor[1] = materialColor[2] = 0.75;

	lightAmbient[0] = lightAmbient[1] = lightAmbient[2] = 0.25f, lightAmbient[3] = 1.f;
	lightDiffuse[0] = lightDiffuse[1] = lightDiffuse[2] = 0.70f, lightDiffuse[3] = 1.f;
	lightSpecular[0] = lightSpecular[1] = lightSpecular[2] = 0.35f, lightSpecular[3] = 1.f;
	shapeSpecular[0] = shapeSpecular[1] = shapeSpecular[2] = 0.35f, shapeSpecular[3] = 1.f;
	shapeSpecularShininess = 128;

	oldX, oldY, newX, newY;
	imageZoom = 1.f, imageOffset[0] = imageOffset[1] = 0.f;
	rotating = scaling = panning = false;
	useLight = true;
	showEdges = false;
	showMesh = true;
	vbo = ebo = 0;
	showSelectedVertices = true;
	showColor = true;
	showLoopPoints = true;
	showVoxels = true;
	showFaces = true;

	offscreen_frame_width = 1200;
	offscreen_frame_height = 1200;

	callBacks.push_back(KeyboardCallBack(this, 'K', "screenshot", "File Name", ScreenshotCallBack));
	callBacks.push_back(KeyboardCallBack(this, 'U', "read camera", "File Name", ReadSceneConfigurationCallBack));
	callBacks.push_back(KeyboardCallBack(this, 'u', "save camera", "File Name", WriteSceneConfigurationCallBack));
	callBacks.push_back(KeyboardCallBack(this, 'l', "toggle light", ToggleLightCallBack));
	callBacks.push_back(KeyboardCallBack(this, 'c', "toggle color", ToggleColorCallBack));
	callBacks.push_back(KeyboardCallBack(this, 'e', "toggle edges", ToggleEdgesCallBack));
	callBacks.push_back(KeyboardCallBack(this, 'm', "toggle mesh", ToggleMeshCallBack));
	callBacks.push_back(KeyboardCallBack(this, 'v', "toggle voxels", ToggleVoxelsCallBack));
	callBacks.push_back(KeyboardCallBack(this, 'f', "toggle faces", ToggleFacesCallBack));
	
}

void SurfaceVisualization::WriteSceneConfigurationCallBack(Visualization* v, const char* prompt){
	const SurfaceVisualization* av = (SurfaceVisualization*)v;
	FILE * file;
	file = fopen(prompt, "wb");
	fwrite(&av->screenWidth, sizeof(int), 1, file);
	fwrite(&av->screenHeight, sizeof(int), 1, file);
	fwrite(&av->camera.position, sizeof(Point3D<double>), 1, file);
	fwrite(&av->camera.forward, sizeof(Point3D<double>), 1, file);
	fwrite(&av->camera.right, sizeof(Point3D<double>), 1, file);
	fwrite(&av->camera.up, sizeof(Point3D<double>), 1, file);
	fwrite(&av->zoom, sizeof(float), 1, file);
	fclose(file);
}

void SurfaceVisualization::ReadSceneConfigurationCallBack(Visualization* v, const char* prompt){
	SurfaceVisualization* av = (SurfaceVisualization*)v;
	FILE * file;
	file = fopen(prompt, "rb");
	if (!file) {
		printf("Camera Configuration File Not Valid \n");
	}
	else{
		fread(&av->screenWidth, sizeof(int), 1, file);
		fread(&av->screenHeight, sizeof(int), 1, file);
		fread(&av->camera.position, sizeof(Point3D<double>), 1, file);
		fread(&av->camera.forward, sizeof(Point3D<double>), 1, file);
		fread(&av->camera.right, sizeof(Point3D<double>), 1, file);
		fread(&av->camera.up, sizeof(Point3D<double>), 1, file);
		fread(&av->zoom, sizeof(float), 1, file);
		//av->offscreen_frame_height = av->screenHeight;
		//av->offscreen_frame_width = av->screenWidth;
		fclose(file);
	}
	glutPostRedisplay();
}

void SurfaceVisualization::SetupOffScreenBuffer(){
	// The depth buffer texture
	glGenTextures(1, &offscreen_depth_texture);
	glBindTexture(GL_TEXTURE_2D, offscreen_depth_texture);
	glTexStorage2D(GL_TEXTURE_2D, 1, GL_DEPTH_COMPONENT24, offscreen_frame_width, offscreen_frame_height);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);

	// The color buffer texture
	glGenTextures(1, &offscreen_color_texture);
	glBindTexture(GL_TEXTURE_2D, offscreen_color_texture);
	glTexStorage2D(GL_TEXTURE_2D, 1, GL_RGBA8, offscreen_frame_width, offscreen_frame_height);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);


	// Create and set up the FBO
	glGenFramebuffers(1, &offscreen_framebuffer_handle);
	glBindFramebuffer(GL_FRAMEBUFFER, offscreen_framebuffer_handle);
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, offscreen_depth_texture, 0);
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, offscreen_color_texture, 0);
	GLenum drawBuffers[] = { GL_COLOR_ATTACHMENT0 };
	glDrawBuffers(1, drawBuffers);

	GLenum result = glCheckFramebufferStatus(GL_FRAMEBUFFER);
	if (result == GL_FRAMEBUFFER_COMPLETE) {
		printf("Framebuffer is complete.\n");
	}
	else {
		printf("Framebuffer is not complete.\n");
	}
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

void SurfaceVisualization::RenderOffScreenBuffer(Image<Point3D<float>> & image){
	if (!offscreen_framebuffer_handle) SetupOffScreenBuffer();
	glViewport(0, 0, offscreen_frame_width, offscreen_frame_height);
	glBindFramebuffer(GL_FRAMEBUFFER, offscreen_framebuffer_handle);
	int windowScreenWidth = screenWidth;
	int windowScreenHeight = screenHeight;
	screenWidth = offscreen_frame_width;
	screenHeight = offscreen_frame_height;

	glClearColor(1, 1, 1, 1);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	display();
	screenWidth = windowScreenWidth;
	screenHeight = windowScreenHeight;
	glFlush();

	//Save color buffer to image
	Pointer(float) GLColorBuffer = AllocPointer< float >(sizeof(float)* 3 * offscreen_frame_width * offscreen_frame_height);
	glReadBuffer(GL_COLOR_ATTACHMENT0);
	glReadPixels(0, 0, offscreen_frame_width, offscreen_frame_height, GL_RGB, GL_FLOAT, GLColorBuffer);
	glFinish();

	image.resize(offscreen_frame_width, offscreen_frame_height);
	for (int i = 0; i<offscreen_frame_width; i++) for (int j = 0; j<offscreen_frame_height; j++)  for (int c = 0; c<3; c++){
		image(i, j)[c] = GLColorBuffer[c + i * 3 + (offscreen_frame_height - 1 - j) * offscreen_frame_width * 3];
	}
	FreePointer(GLColorBuffer);
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
	glViewport(0, 0, screenWidth, screenHeight);
}

void SurfaceVisualization::ScreenshotCallBack(Visualization* v, const char* prompt){
	Image<Point3D<float>> image;
	SurfaceVisualization* av = (SurfaceVisualization*)v;
	av->RenderOffScreenBuffer(image);
	image.write(prompt);
}

int SurfaceVisualization::selectVertex(int x, int y){
	Point3D<float> point;
	if (select(x, y, point)){
		int vCount = vertices.size();
		int minIndex = -1;
		float minDistance = FLT_MAX;
		for (int i = 0; i < vCount; i++) if (Point3D<float>::SquareNorm(vertices[i] - point) < minDistance){
			minDistance = Point3D<float>::SquareNorm(vertices[i] - point);
			minIndex = i;
		}
		return minIndex;

		glutPostRedisplay();
	}
	else{
		return -1;
	}
}

bool SurfaceVisualization::setPosition(int x, int y, Point3D< double >& p)
{
	double _x = (double)x / screenWidth - 0.5, _y = 1. - (double)y / screenHeight - 0.5;
	_x *= 2., _y *= 2;
	_x *= zoom, _y *= zoom;
	double r = _x*_x + _y*_y;
	if (r<1)
	{
		p = camera.forward * (-sqrt(1 - r)) + camera.right * _x + camera.up * _y;
		return true;
	}
	return false;
}
bool SurfaceVisualization::setPosition(int x, int y, Point3D< float >& p)
{
	Point3D< double > _p;
	bool ret = setPosition(x, y, _p);
	p = Point3D< float >((float)_p[0], (float)_p[1], (float)_p[2]);
	return ret;
}

void SurfaceVisualization::updateMesh()
{
	TriangleIndex *_triangles = new TriangleIndex[triangles.size()];

	for (int i = 0; i<triangles.size(); i++) for (int j = 0; j<3; j++) _triangles[i][j] = 3 * i + j;

	if (!glIsBuffer(ebo))glGenBuffers(1, &ebo);
	if (!glIsBuffer(vbo))glGenBuffers(1, &vbo);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, triangles.size() * sizeof(int) * 3, _triangles, GL_STATIC_DRAW);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

	delete[] _triangles;

	Point3D< float > *_vertices = new Point3D< float >[triangles.size() * 9];
	Point3D< float > *_normals = _vertices + triangles.size() * 3;
	Point3D< float > *_colors = _vertices + triangles.size() * 6;

	for (int i = 0; i<triangles.size(); i++)
	{
		Point3D< float > n = Point3D< float >::CrossProduct(vertices[triangles[i][1]] - vertices[triangles[i][0]], vertices[triangles[i][2]] - vertices[triangles[i][0]]);
		n /= Length(n);
		for (int j = 0; j<3; j++)
		{
			_vertices[3 * i + j] = vertices[triangles[i][j]];
			_colors[3 * i + j] = colors[3*i + j];
			_normals[3 * i + j] = n;
		}
	}
	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glBufferData(GL_ARRAY_BUFFER, 9 * triangles.size() * sizeof(Point3D< float >), _vertices, GL_DYNAMIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	delete[] _vertices;
	glutPostRedisplay();
}

bool SurfaceVisualization::select(int x, int  y, Point3D< float >& out)
{
	bool ret = false;
	Pointer(float) depthBuffer = AllocPointer< float >(sizeof(float)* screenWidth * screenHeight);
	glReadPixels(0, 0, screenWidth, screenHeight, GL_DEPTH_COMPONENT, GL_FLOAT, depthBuffer);
	float ar = (float)screenWidth / (float)screenHeight;
	float _screenWidth, _screenHeight;
	if (screenWidth>screenHeight) _screenWidth = screenWidth * ar, _screenHeight = screenHeight;
	else                           _screenWidth = screenWidth, _screenHeight = screenHeight / ar;
	{
		double _x = (double)x / screenWidth - 0.5, _y = 1. - (double)y / screenHeight - 0.5, _z;
		if (screenWidth>screenHeight) _x *= zoom*ar, _y *= zoom;
		else                           _x *= zoom, _y *= zoom / ar;
		_x *= 2., _y *= 2;
		int x1 = (int)floor(x), y1 = (int)floor(y), x2 = x1 + 1, y2 = y1 + 1;
		float dx = x - x1, dy = y - y1;
		x1 = std::max< int >(0.f, std::min< int >(x1, screenWidth - 1));
		y1 = std::max< int >(0.f, std::min< int >(y1, screenHeight - 1));
		x2 = std::max< int >(0.f, std::min< int >(x2, screenWidth - 1));
		y2 = std::max< int >(0.f, std::min< int >(y2, screenHeight - 1));
		_z =
			depthBuffer[(screenHeight - 1 - y1)*screenWidth + x1] * (1.f - dx) * (1.f - dy) +
			depthBuffer[(screenHeight - 1 - y1)*screenWidth + x2] * (dx)* (1.f - dy) +
			depthBuffer[(screenHeight - 1 - y2)*screenWidth + x1] * (1.f - dx) * (dy)+
			depthBuffer[(screenHeight - 1 - y2)*screenWidth + x2] * (dx)* (dy);
		if (_z<1) out = Point3D< float >(camera.forward * (-1.5 + 3. * _z) + camera.right * _x + camera.up * _y + camera.position), ret = true;
	}
	FreePointer(depthBuffer);
	return ret;
}

void DrawUnitSphere(){
	GLUquadric* quad = gluNewQuadric();
	glPushMatrix();
	glTranslatef(0.5f,0.5f,0.5f);
	gluSphere(quad, 0.3f, 20, 20);
	glPopMatrix();
	gluDeleteQuadric(quad);
}


void DrawPolygonVoxelFace(unsigned long c){
	if (c == 0){
		glBegin(GL_TRIANGLES);
		glVertex3f(0.f, 0.f, 0.f);
		glVertex3f(0.f, 1.f, 0.f);
		glVertex3f(0.f, 1.f, 1.f);

		glVertex3f(0.f, 0.f, 0.f);
		glVertex3f(0.f, 1.f, 1.f);
		glVertex3f(0.f, 0.f, 1.f);
		glEnd();
	}
	else if (c == 1){
		glBegin(GL_TRIANGLES);
		glVertex3f(0.f, 0.f, 0.f);
		glVertex3f(1.f, 0.f, 0.f);
		glVertex3f(1.f, 0.f, 1.f);

		glVertex3f(0.f, 0.f, 0.f);
		glVertex3f(1.f, 0.f, 1.f);
		glVertex3f(0.f, 0.f, 1.f);
		glEnd();
	}
	else if (c == 2){
		glBegin(GL_TRIANGLES);
		glVertex3f(0.f, 0.f, 0.f);
		glVertex3f(1.f, 0.f, 0.f);
		glVertex3f(1.f, 1.f, 0.f);

		glVertex3f(0.f, 0.f, 0.f);
		glVertex3f(1.f, 1.f, 0.f);
		glVertex3f(0.f, 1.f, 0.f);
		glEnd();
	}
}

void DrawWireFrameBox(){
	glBegin(GL_LINE_LOOP);
	glVertex3f(0.f, 0.f, 0.f);
	glVertex3f(1.f, 0.f, 0.f);
	glVertex3f(1.f, 1.f, 0.f);
	glVertex3f(0.f, 1.f, 0.f);
	glEnd();

	glBegin(GL_LINE_LOOP);
	glVertex3f(0.f, 0.f, 1.f);
	glVertex3f(1.f, 0.f, 1.f);
	glVertex3f(1.f, 1.f, 1.f);
	glVertex3f(0.f, 1.f, 1.f);
	glEnd();

	glBegin(GL_LINES);
	glVertex3f(0.f, 0.f, 0.f);
	glVertex3f(0.f, 0.f, 1.f);

	glVertex3f(1.f, 0.f, 0.f);
	glVertex3f(1.f, 0.f, 1.f);

	glVertex3f(1.f, 1.f, 0.f);
	glVertex3f(1.f, 1.f, 1.f);

	glVertex3f(0.f, 1.f, 0.f);
	glVertex3f(0.f, 1.f, 1.f);
	glEnd();
}

void DrawOriginBasedEdge(unsigned long c){
	Point3D<float> edge(0, 0, 0);
	edge[c] = 1.f;
	glBegin(GL_LINES);
	glVertex3f(0,0,0);
	glVertex3f(edge[0], edge[1], edge[2]);
	glEnd();
}

void SurfaceVisualization::DrawEdge(Point3D<float> corner, unsigned long c){
	glMatrixMode(GL_MODELVIEW_MATRIX);
	glPushMatrix();
	glTranslatef(integerGridTranslation[0], integerGridTranslation[1], integerGridTranslation[2]);
	glScalef(integerGridScaling, integerGridScaling, integerGridScaling);
	glTranslatef(corner[0], corner[1], corner[2]);
	DrawOriginBasedEdge(c);
	glPopMatrix();
}

void SurfaceVisualization::DrawEdge(unsigned long key){
	unsigned long axialIndices[3] = { static_cast<unsigned long>((key >> 20) & 0x000003FF), static_cast<unsigned long>((key >> 10) & 0x000003FF), static_cast<unsigned long>((key & 0x000003FF)) };
	unsigned long cIndex = static_cast<unsigned long>((key >> 30) & 0x00000003);
	DrawEdge(Point3D<float>(axialIndices[0], axialIndices[1], axialIndices[2]), cIndex);
}

void SurfaceVisualization::DrawFace(Point3D<float> corner, unsigned long c){
	glMatrixMode(GL_MODELVIEW_MATRIX);
	glPushMatrix();
	glTranslatef(integerGridTranslation[0], integerGridTranslation[1], integerGridTranslation[2]);
	glScalef(integerGridScaling, integerGridScaling, integerGridScaling);
	glTranslatef(corner[0], corner[1], corner[2]);
	DrawPolygonVoxelFace(c);
	glPopMatrix();
}


void SurfaceVisualization::DrawFace(unsigned long key){
	unsigned long axialIndices[3] = { static_cast<unsigned long>((key >> 20) & 0x000003FF), static_cast<unsigned long>((key >> 10) & 0x000003FF), static_cast<unsigned long>((key & 0x000003FF)) };
	unsigned long cIndex = static_cast<unsigned long>((key >> 30) & 0x00000003);
	DrawFace(Point3D<float>(axialIndices[0], axialIndices[1], axialIndices[2]), cIndex);
}


void SurfaceVisualization::DrawVoxel(Point3D<float> corner){
	glMatrixMode(GL_MODELVIEW_MATRIX);
	glPushMatrix();
	glTranslatef(integerGridTranslation[0], integerGridTranslation[1], integerGridTranslation[2]);
	glScalef(integerGridScaling, integerGridScaling, integerGridScaling);
	glTranslatef(corner[0], corner[1], corner[2]);
	DrawWireFrameBox();
	glPopMatrix();
}

void SurfaceVisualization::DrawVoxel(unsigned long key){
	unsigned long axialIndices[3] = { static_cast<unsigned long>((key >> 20) & 0x000003FF), static_cast<unsigned long>((key >> 10) & 0x000003FF), static_cast<unsigned long>((key & 0x000003FF)) };
	DrawVoxel(Point3D<float>(axialIndices[0], axialIndices[1], axialIndices[2]));
}

void SurfaceVisualization::DrawSphere(unsigned long key){
	unsigned long axialIndices[3] = { static_cast<unsigned long>((key >> 20) & 0x000003FF), static_cast<unsigned long>((key >> 10) & 0x000003FF), static_cast<unsigned long>((key & 0x000003FF)) };
	DrawSphere(Point3D<float>(axialIndices[0], axialIndices[1], axialIndices[2]));
}

void SurfaceVisualization::DrawSphere(Point3D<float> corner){
	glMatrixMode(GL_MODELVIEW_MATRIX);
	glPushMatrix();
	glTranslatef(integerGridTranslation[0], integerGridTranslation[1], integerGridTranslation[2]);
	glScalef(integerGridScaling, integerGridScaling, integerGridScaling);
	glTranslatef(corner[0], corner[1], corner[2]);
	DrawUnitSphere();
	glPopMatrix();
}

void SurfaceVisualization::display(void)
{
	if (!vbo && !ebo) updateMesh();
	glEnable(GL_CULL_FACE);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	float ar = (float)screenWidth / (float)screenHeight, ar_r = 1.f / ar;
	if (screenWidth>screenHeight) glOrtho(-ar*zoom, ar*zoom, -zoom, zoom, -1.5, 1.5);
	else                           glOrtho(-zoom, zoom, -ar_r*zoom, ar_r*zoom, -1.5, 1.5);
	glMatrixMode(GL_MODELVIEW);

	glLoadIdentity();
	camera.draw();

	GLfloat lPosition[4];

	{
		Point3D< float > d = camera.up + camera.right - camera.forward * 5;
		lPosition[0] = d[0], lPosition[1] = d[1], lPosition[2] = d[2];
	}
	lPosition[3] = 0.0;
	glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_FALSE);
	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
	glLightfv(GL_LIGHT0, GL_AMBIENT, lightAmbient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, lightDiffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, lightSpecular);
	glLightfv(GL_LIGHT0, GL_POSITION, lPosition);
	glEnable(GL_LIGHT0);
	if (useLight) glEnable(GL_LIGHTING);
	else           glDisable(GL_LIGHTING);

	glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
	glEnable(GL_COLOR_MATERIAL);
	glColor3f(materialColor[0], materialColor[1], materialColor[2]);
	glEnable(GL_DEPTH_TEST);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, shapeSpecular);
	glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, shapeSpecularShininess);

	glEnable(GL_DEPTH_TEST);

	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
	glEnableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(3, GL_FLOAT, 0, NULL);
	glNormalPointer(GL_FLOAT, 0, (GLubyte*)NULL + sizeof(Point3D< float >) * triangles.size() * 3);
	glColorPointer(3, GL_FLOAT, 0, (GLubyte*)NULL + sizeof(Point3D< float >) * triangles.size() * 6);
	glEnableClientState(GL_NORMAL_ARRAY);
	if (showColor) glEnableClientState(GL_COLOR_ARRAY);
	if (showMesh) glDrawElements(GL_TRIANGLES, (GLsizei)(triangles.size() * 3), GL_UNSIGNED_INT, NULL);
	glDisableClientState(GL_NORMAL_ARRAY);
	if (showColor) glDisableClientState(GL_COLOR_ARRAY);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

	if (showSelectedVertices) {
		glColor3f(0.f, 1.f, 0.f);
		GLUquadric* quad = gluNewQuadric();
		for (int i = 0; i < selectedVertices.size(); i++) {
			glPushMatrix();
			glTranslatef(vertices[selectedVertices[i]][0], vertices[selectedVertices[i]][1], vertices[selectedVertices[i]][2]);
			gluSphere(quad, 0.01f, 20, 20);
			glPopMatrix();
		}
		gluDeleteQuadric(quad);
	}


	if (showEdges)
	{
		GLint src, dst;
		glGetIntegerv(GL_BLEND_SRC, &src);
		glGetIntegerv(GL_BLEND_DST, &dst);
		Point3D< float > f = camera.forward / 256;
		glPushMatrix();
		glTranslatef(-f[0], -f[1], -f[2]);
		glColor3f(0., 0., 0.);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glEnable(GL_BLEND);
		glEnable(GL_LINE_SMOOTH);
		glLineWidth(1.0f);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glBindBuffer(GL_ARRAY_BUFFER, vbo);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
		glEnableClientState(GL_VERTEX_ARRAY);
		glVertexPointer(3, GL_FLOAT, 0, NULL);
		glDrawElements(GL_TRIANGLES, (GLsizei)(triangles.size() * 3), GL_UNSIGNED_INT, NULL);
		glBindBuffer(GL_ARRAY_BUFFER, 0);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		glDisable(GL_LINE_SMOOTH);
		glPopMatrix();
		glDisable(GL_BLEND);
		glBlendFunc(src, dst);
	}

	if (showLoopPoints){
		GLint src, dst;
		glGetIntegerv(GL_BLEND_SRC, &src);
		glGetIntegerv(GL_BLEND_DST, &dst);
		Point3D< float > f = camera.forward / 256;
		glPushMatrix();
		glTranslatef(-f[0], -f[1], -f[2]);
		glColor3f(1.f, 0.f, 0.f);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glEnable(GL_BLEND);
		glEnable(GL_LINE_SMOOTH);
		glLineWidth(2.0f);
		for (int i = 0; i < loopPoints.size(); i++){
			glBegin(GL_LINE_LOOP);
			for (int j = 0; j < loopPoints[i].size(); j++){
				glVertex3f(loopPoints[i][j][0], loopPoints[i][j][1], loopPoints[i][j][2]);
			}
			glEnd();
		}
		glDisable(GL_LINE_SMOOTH);
		glPopMatrix();
		glDisable(GL_BLEND);
		glBlendFunc(src, dst);
	}

	if (showVoxels){
		glDisable(GL_LIGHTING);
		glColor3f(0.f, 1.f, 0.f);
		glEnable(GL_BLEND);
		glEnable(GL_LINE_SMOOTH);
		glLineWidth(1.5f);
		for (int i = 0; i < voxels.size(); i++){
			DrawVoxel(voxels[i]);
		}
		glDisable(GL_LINE_SMOOTH);
		glDisable(GL_BLEND);
		if (useLight) glEnable(GL_LIGHTING);
	}

	if (showFaces){
		glDisable(GL_CULL_FACE);
		glDisable(GL_LIGHTING);
		glColor3f(1.f, 0.f, 1.f);
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		for (int i = 0; i < faces.size(); i++){
			DrawFace(faces[i]);
		}
		if (useLight) glEnable(GL_LIGHTING);
		glEnable(GL_CULL_FACE);
	}

}


void SurfaceVisualization::mouseFunc(int button, int state, int x, int y)
{
	newX = x; newY = y;

	rotating = scaling = panning = false;
	if (button == GLUT_LEFT_BUTTON)
	if (glutGetModifiers() & GLUT_ACTIVE_CTRL) panning = true;
	else                                        rotating = true;
	else if (button == GLUT_RIGHT_BUTTON) scaling = true;
}
void SurfaceVisualization::motionFunc(int x, int y)
{
	oldX = newX, oldY = newY, newX = x, newY = y;

	int imageSize = std::min< int >(screenWidth, screenHeight);
	float rel_x = (newX - oldX) / (float)imageSize * 2;
	float rel_y = (newY - oldY) / (float)imageSize * 2;
	float pRight = -rel_x * zoom, pUp = rel_y * zoom;
	float sForward = rel_y * 4;
	float rRight = rel_y, rUp = rel_x;

	if (rotating) camera.rotateUp(rUp), camera.rotateRight(rRight);
	else if (scaling) zoom *= (float)pow(0.9, (double)sForward);
	else if (panning) camera.translate(camera.right * pRight + camera.up * pUp);

	glutPostRedisplay();
}
