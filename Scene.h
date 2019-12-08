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
	Matrix4 OrthographicProjection(Camera* camera);
	Matrix4 PerspectiveProjection(Camera* camera);
	Matrix4 ViewportProjection(Camera* camera);
	Matrix4 getCameraTransformMatrix(Camera* camera);
	Vec4* getVector4(Vec3 vector);
	Vec3* getVector3(Vec4 vector);
	void modelTransformation(Matrix4 worldMatrix,Camera *camera);
	void clipping(Vec3 *v0, Vec3 *v1, Camera* camera);
	bool isVisible(float d,float num,float *tEnter,float *tLeave);
	Vec3 getPointAtt(Vec3 v0, Vec3 v1, float t);
	
};

#endif
