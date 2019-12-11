#ifndef _SCENE_H_
#define _SCENE_H_

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#include "Camera.h"
#include "Color.h"
#include "Model.h"
#include "Rotation.h"
#include "Scaling.h"
#include "Translation.h"
#include "Triangle.h"
#include "Vec3.h"
#include "Vec4.h"
#include  "Matrix4.h"
using namespace std;

class Scene
{
public:
	Color backgroundColor;
	bool cullingEnabled;
	int projectionType;
	vector< vector<Color> > image;
	vector< Camera* > cameras;
	vector< Vec3* > vertices;
	vector< Color* > colorsOfVertices;
	vector< Scaling* > scalings;
	vector< Rotation* > rotations;
	vector< Translation* > translations;
	vector< Model* > models;
	Scene(const char *xmlPath);
	void initializeImage(Camera* camera);
	void forwardRenderingPipeline(Camera* camera);
	int makeBetweenZeroAnd255(double value);
	void writeImageToPPMFile(Camera* camera);
	void convertPPMToPNG(string ppmFileName, int osType);

	// Implemented as hw
	vector < vector<Vec4*> > verticesAssembled;
	vector<Vec4*> newVertices;
	void transformation(Matrix4 transformationMatrix, Camera* camera);
	void modelingTransformation();
	Vec3* getVector3(Vec4 vector);
	Vec4* getVector4(Vec3 vector);
	void setModelw();
	void lineRasterization(int x_0, int y_0, Color* c_0, int x_1, int y_1, Color* c_1, Camera* camera);
	void triangleRasterization(Camera* camera,
		int x_0, int y_0, Color* c_0,
		int x_1, int y_1, Color* c_1,
		int x_2, int y_2, Color* c_2);
	void rasterization(Camera* camera);
	void clipColors(Vec4* v0, Vec4* v1, Vec4* v_clipped);
	bool isVisible(double d, double num, double* tEnter, double* tLeave);
	void clipping(Vec4* v0, Vec4* v1, Vec4* v0_clipped, Vec4* v1_clipped, Camera* camera);
	void clippingModels(Camera* camera);
	void backfaceCulling(Camera* camera);
	Matrix4 getCameraTransformMatrix(Camera* camera);
	Matrix4 getOrthographicProjectionMatrix(Camera* camera);
	Matrix4 getPerspectiveProjectionMatrix(Camera* camera);
	Matrix4 getViewportProjectionMatrix(Camera* camera);
	Matrix4 getTranslationMatrix(Translation* translation);
	Matrix4 getScalingMatrix(Scaling* scaling);
	Matrix4 getRotationMatrix(Rotation* rotation);
};

#endif
