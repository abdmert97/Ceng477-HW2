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
	// Implemented as hw
	vector < vector<Vec3*> > verticesAssembled;

	Scene(const char *xmlPath);
	
	void initializeImage(Camera* camera);
	void forwardRenderingPipeline(Camera* camera);
	int makeBetweenZeroAnd255(double value);
	void writeImageToPPMFile(Camera* camera);
	void convertPPMToPNG(string ppmFileName, int osType);

	// Implemented as hw
	Matrix4 getOrthographicProjectionMatrix(Camera* camera);
	Matrix4 getPerspectiveProjectionMatrix(Camera* camera);
	Matrix4 getViewportProjectionMatrix(Camera* camera);
	Matrix4 getCameraTransformMatrix(Camera* camera);
	Vec4* getVector4(Vec3 vector);
	Vec3* getVector3(Vec4 vector);
	void modelTransformation(Matrix4 worldMatrix,Camera *camera);
	void clipping(Vec3 *v0, Vec3 *v1, Camera* camera);
	bool isVisible(float d,float num,float *tEnter,float *tLeave);
	Vec3 getPointAtt(Vec3 v0, Vec3 v1, float t);
	void rasterization(Camera* camera);
	void transformation(Matrix4 transformationMatrix, Camera* camera);
	void lineRasterization(int x_0, int y_0, Color* c_0, int x_1, int y_1, Color* c_1);
	void modelingTransformation();
	Matrix4 getTranslationMatrix(Translation* translation);
	Matrix4 getScalingMatrix(Scaling* scaling);
	void triangleRasterization(Camera* camera, int x_0, int y_0, int x_1, int y_1, int x_2, int y_2);
};

#endif
