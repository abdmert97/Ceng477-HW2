#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <fstream>
#include <cmath>

#include "Scene.h"
#include "Camera.h"
#include "Color.h"
#include "Model.h"
#include "Rotation.h"
#include "Scaling.h"
#include "Translation.h"
#include "Triangle.h"
#include "Vec3.h"
#include "tinyxml2.h"
#include "Helpers.h"
#include <algorithm>

using namespace tinyxml2;
using namespace std;


/********************************************
 * GENERAL HELPER FUNCTIONS					*
 ********************************************/
void Scene::transformation(Matrix4 transformationMatrix, Camera* camera)
{
	for (int i = 0; i < models.size(); i++)
	{
		for (int j = 0; j < models[i]->numberOfTriangles; j++)
		{
			Vec4 firstVertice = *verticesAssembled[i][j * 3];
			Vec4 secondVertice = *verticesAssembled[i][j * 3 + 1];
			Vec4 thirdVertice = *verticesAssembled[i][j * 3 + 2];


			Vec4 first = multiplyMatrixWithVec4(transformationMatrix, firstVertice);
			Vec4 second = multiplyMatrixWithVec4(transformationMatrix, secondVertice);
			Vec4 third = multiplyMatrixWithVec4(transformationMatrix, thirdVertice);

			*verticesAssembled[i][j * 3] = first;
			*verticesAssembled[i][j * 3 + 1] = second;
			*verticesAssembled[i][j * 3 + 2] = third;
		}
	}

	for (int i = 0; i < newVertices.size(); i++)
	{
		*newVertices[i] = multiplyMatrixWithVec4(transformationMatrix, *newVertices[i]);
	}

}

void Scene::modelingTransformation()
{
	for (int i = 0; i < models.size(); i++)
	{
		for (int j = 0; j < models[i]->numberOfTransformations; j++)
		{
			if (models[i]->transformationTypes[j] == 't')
			{
				// TRANSLATION
				Matrix4 translationMatrix = getTranslationMatrix(translations[models[i]->transformationIds[j] - 1]);
				for (int k = 0; k < models[i]->numberOfTriangles; k++)
				{
					Vec4 firstVertice = *verticesAssembled[i][k * 3];
					Vec4 secondVertice = *verticesAssembled[i][k * 3 + 1];
					Vec4 thirdVertice = *verticesAssembled[i][k * 3 + 2];
					Vec4 first = multiplyMatrixWithVec4(translationMatrix, firstVertice);
					Vec4 second = multiplyMatrixWithVec4(translationMatrix, secondVertice);
					Vec4 third = multiplyMatrixWithVec4(translationMatrix, thirdVertice);
					*verticesAssembled[i][k * 3] = first;
					*verticesAssembled[i][k * 3 + 1] = second;
					*verticesAssembled[i][k * 3 + 2] = third;
				}
			}
			else if (models[i]->transformationTypes[j] == 's')
			{
				// SCALING
				Matrix4 scalingMatrix = getScalingMatrix(scalings[models[i]->transformationIds[j] - 1]);
				
				for (int k = 0; k < models[i]->numberOfTriangles; k++)
				{
					Vec4 firstVertice = *verticesAssembled[i][k * 3];
					Vec4 secondVertice = *verticesAssembled[i][k * 3 + 1];
					Vec4 thirdVertice = *verticesAssembled[i][k * 3 + 2];

					Vec4 first = multiplyMatrixWithVec4(scalingMatrix, firstVertice);
					Vec4 second = multiplyMatrixWithVec4(scalingMatrix, secondVertice);
					Vec4 third = multiplyMatrixWithVec4(scalingMatrix, thirdVertice);
					*verticesAssembled[i][k * 3] = first;
					*verticesAssembled[i][k * 3 + 1] = second;
					*verticesAssembled[i][k * 3 + 2] = third;
				}
			}
			else if (models[i]->transformationTypes[j] == 'r')
			{
				// ROTATING
				Matrix4 rotationMatrix = getRotationMatrix(rotations[models[i]->transformationIds[j] - 1]);
		
				for (int k = 0; k < models[i]->numberOfTriangles; k++)
				{
					Vec4 firstVertice = *verticesAssembled[i][k * 3];
					Vec4 secondVertice = *verticesAssembled[i][k * 3 + 1];
					Vec4 thirdVertice = *verticesAssembled[i][k * 3 + 2];

					Vec4 first = multiplyMatrixWithVec4(rotationMatrix, firstVertice);
					Vec4 second = multiplyMatrixWithVec4(rotationMatrix, secondVertice);
					Vec4 third = multiplyMatrixWithVec4(rotationMatrix, thirdVertice);
					*verticesAssembled[i][k * 3] = first;
					*verticesAssembled[i][k * 3 + 1] = second;
					*verticesAssembled[i][k * 3 + 2] = third;
				}
			}
		}

	}
}

Vec3* Scene::getVector3(Vec4 vector)
{
	return new Vec3(vector.x, vector.y, vector.z, vector.colorId);
}

Vec4* Scene::getVector4(Vec3 vector)
{
	Vec4* vec = new Vec4(vector.x, vector.y, vector.z, 1, vector.colorId);
	return  vec;
}

void Scene::setModelw()
{
	for (int i = 0; i < models.size(); i++)
	{
		for (int j = 0; j < models[i]->numberOfTriangles; j++)
		{
			Vec4* v0 = verticesAssembled[i][j * 3];
			Vec4* v1 = verticesAssembled[i][j * 3 + 1];
			Vec4* v2 = verticesAssembled[i][j * 3 + 2];

			v0->x /= v0->t;
			v0->y /= v0->t;
			v0->z /= v0->t;
			v0->t /= v0->t;

			v1->x /= v1->t;
			v1->y /= v1->t;
			v1->z /= v1->t;
			v1->t /= v1->t;

			v2->x /= v2->t;
			v2->y /= v2->t;
			v2->z /= v2->t;
			v2->t /= v2->t;
		}
	}

}


/********************************************
 * HELPER FUNCTIONS FOR RASTERIZATION		*
 ********************************************/
void Scene::lineRasterization(int x_0, int y_0, Color* c_0, int x_1, int y_1, Color* c_1, Camera* camera)
{
	double dx = x_1 - x_0;
	double dy = y_1 - y_0;
	// TODO: dx == 0
	double slope = dy / (dx+EPSILON);

	int xlimit = max(x_0, x_1);
	int ylimit = max(y_0, y_1);
	int xmin = min(x_0, x_1);
	int ymin = min(y_0, y_1);
	
	if (slope < -1)
	{
		int	x = max(x_0, x_1);
		
		int d = 2 * abs(dx) - abs(dy);

		for (int y = ymin; y < ylimit; y++)
		{
			if (x >= 0 && y >= 0 && x < camera->horRes && y < camera->verRes)
			{
				double colorScale_v0 = abs((y_1 - y) / dy);
				double colorScale_v1 = abs((y_0 - y) / dy);
				this->image[x][y] = Color(
					colorScale_v0 * c_0->r + colorScale_v1 * c_1->r,
					colorScale_v0 * c_0->g + colorScale_v1 * c_1->g,
					colorScale_v0 * c_0->b + colorScale_v1 * c_1->b);
			}
			if (d <= 0)
			{
				d += 2 * (abs(dx));
			}
			else
			{
				d += 2 * (abs(dx) - abs(dy));
				x--;
			}
		}
	}
	else if (slope >= -1 && slope <= 0 )
	{
		int y = max(y_0, y_1);

		int d = 2 * abs(dy) - abs(dx);

		for (int x = xmin; x < xlimit; x++)
		{
			if (x >= 0 && y >= 0 && x < camera->horRes && y < camera->verRes)
			{
				double colorScale_v0 = abs((x_1 - x) / dx);
				double colorScale_v1 = abs((x_0 - x) / dx);
				this->image[x][y] = Color(
					colorScale_v0 * c_0->r + colorScale_v1 * c_1->r,
					colorScale_v0 * c_0->g + colorScale_v1 * c_1->g,
					colorScale_v0 * c_0->b + colorScale_v1 * c_1->b);
			}
			if (d <= 0)
			{
				d += 2 * (abs(dy));
			}
			else
			{
				d += 2 * (abs(dy) - abs(dx));
				y--;
			}
		}
	}

	else if (0 < slope && slope < 1)
	{
	
		int y = min(y_0, y_1);
		int d = 2 * abs(dy) - abs(dx);
		

		for (int x = xmin; x < xlimit; x++)
		{
			if (x >= 0 && y >= 0 && x < camera->horRes && y < camera->verRes)
			{
				double colorScale_v0 = abs((x_1 - x) / dx);
				double colorScale_v1 = abs((x_0 - x) / dx);
				this->image[x][y] = Color(
					colorScale_v0 * c_0->r + colorScale_v1 * c_1->r,
					colorScale_v0 * c_0->g + colorScale_v1 * c_1->g,
					colorScale_v0 * c_0->b + colorScale_v1 * c_1->b);
			}
			if (d <= 0)
			{
				d += 2 * (abs(dy));
			}
			else
			{
				d += 2 * (abs(dy) - abs(dx));
				y++;
			}
		}
	}
	else if (slope >= 1)
	{
		int	x = min(x_0, x_1);


		int d = 2 * abs(dx) - abs(dy);

		for (int y = ymin; y < ylimit; y++)
		{
			if (x >= 0 && y >= 0 && x < camera->horRes && y < camera->verRes)
			{
				double colorScale_v0 = abs((y_1 - y) / dy);
				double colorScale_v1 = abs((y_0 - y) / dy);
				this->image[x][y] = Color(
					colorScale_v0 * c_0->r + colorScale_v1 * c_1->r,
					colorScale_v0 * c_0->g + colorScale_v1 * c_1->g,
					colorScale_v0 * c_0->b + colorScale_v1 * c_1->b);
			}
			if (d <= 0)
			{
				d += 2 * (abs(dx));
			}
			else
			{
				d += 2 * (abs(dx) - abs(dy));
				x++;
			}
		}
	}
	
	
	/*
	int y = y_0;
	float d = 2 * (y_0 - y_1) + (x_1 - x_0);
	Vec3* c_0v = new Vec3(c_0->r, c_0->g, c_0->b, -1);
	Vec3* c_1v = new Vec3(c_1->r, c_1->g, c_1->b, -1);
	Vec3* c = new Vec3(*c_0v);
	Vec3* d_c = new Vec3(255, 255, 255, -1);
	if(x_0 - x_1 != 0)
	{
		*d_c = multiplyVec3WithScalar(addVec3(*c_1v, multiplyVec3WithScalar(*c_0v, -1)), 1 / (x_1 - x_0));
	}
	for (int x = x_0; x <= x_1; x++)
	{
		if(x >= 0 && y >= 0 && x < camera->horRes && y < camera->horRes)
		this->image[x][y] = Color(round(c->x), round(c->y), round(c->z));
		//this->image[x][y] = Color(255, 255, 255);
		if (d < 0)
		{
			y += 1;
			d += 2 * ((y_0 - y_1) + (x_1 - x_0));
		}
		else
		{
			d += 2 * (y_0 - y_1);
		}
		*c = addVec3(*c, *d_c);
	}
	*/
}

void Scene::triangleRasterization(Camera* camera,
	int x_0, int y_0, Color* c_0,
	int x_1, int y_1, Color* c_1,
	int x_2, int y_2, Color* c_2)
{

	int x_min = min(min(x_0, x_1), x_2);
	int y_min = min(min(y_0, y_1), y_2);

	int x_max = max(max(x_0, x_1), x_2);
	int y_max = max(max(y_0, y_1), y_2);

	for (int y = y_min; y < y_max; y++)
	{
		for (int x = x_min; x < x_max; x++)
		{

			double f_01_xy = x * (y_0 - y_1) + y * (x_1 - x_0) + x_0 * y_1 - y_0 * x_1;
			double f_12_xy = x * (y_1 - y_2) + y * (x_2 - x_1) + x_1 * y_2 - y_1 * x_2;
			double f_20_xy = x * (y_2 - y_0) + y * (x_0 - x_2) + x_2 * y_0 - y_2 * x_0;
			double f_01_x2y2 = x_2 * (y_0 - y_1) + y_2 * (x_1 - x_0) + x_0 * y_1 - y_0 * x_1;
			double f_12_x0y0 = x_0 * (y_1 - y_2) + y_0 * (x_2 - x_1) + x_1 * y_2 - y_1 * x_2;
			double f_20_x1y1 = x_1 * (y_2 - y_0) + y_1 * (x_0 - x_2) + x_2 * y_0 - y_2 * x_0;

			double alpha = f_12_xy / f_12_x0y0;
			double beta = f_20_xy / f_20_x1y1;
			double gamma = f_01_xy / f_01_x2y2;
			if (alpha >= 0 && beta >= 0 && gamma >= 0)
			{
				double r = alpha * c_0->r + beta * c_1->r + gamma * c_2->r;
				double g = alpha * c_0->g + beta * c_1->g + gamma * c_2->g;
				double b = alpha * c_0->b + beta * c_1->b + gamma * c_2->b;
				
				if (x >= 0 && y >= 0 && x < camera->horRes && y < camera->verRes)
				{
					image[x][y].r = r;
					image[x][y].g = g;
					image[x][y].b = b;
				}
			}

		}
	}
}

void Scene::rasterization(Camera* camera)
{

	for (int i = 0; i < models.size(); i++)
	{

		for (int j = 0; j < models[i]->numberOfTriangles; j++)
		{
			Vec4 firstVertice = *verticesAssembled[i][j * 3];
			Vec4 secondVertice = *verticesAssembled[i][j * 3 + 1];
			Vec4 thirdVertice = *verticesAssembled[i][j * 3 + 2];

			// LINE RASTERIZATION

			if (models[i]->type == 0)
			{
			/*	continue;
				int x_0, y_0, x_1, y_1;
				// vertice1 -> vertice2


				if (firstVertice.t == 1 && secondVertice.t == 1)
				{
					x_0 = round(firstVertice.x);
					y_0 = round(firstVertice.y);
					x_1 = round(secondVertice.x);
					y_1 = round(secondVertice.y);
					lineRasterization(x_0, y_0, colorsOfVertices[firstVertice.colorId - 1], x_1, y_1, colorsOfVertices[secondVertice.colorId - 1], camera);
				}


				// vertice2 -> vertice3

				if (firstVertice.t == 1 && secondVertice.t == 1)
				{
					x_0 = round(secondVertice.x);
					y_0 = round(secondVertice.y);
					x_1 = round(thirdVertice.x);
					y_1 = round(thirdVertice.y);
					lineRasterization(x_0, y_0, colorsOfVertices[secondVertice.colorId - 1], x_1, y_1, colorsOfVertices[thirdVertice.colorId - 1], camera);
				}

				// vertice3 -> vertice1

				if (firstVertice.t == 1 && secondVertice.t == 1)
				{
					x_0 = round(thirdVertice.x);
					y_0 = round(thirdVertice.y);
					x_1 = round(firstVertice.x);
					y_1 = round(firstVertice.y);
					lineRasterization(x_0, y_0, colorsOfVertices[thirdVertice.colorId - 1], x_1, y_1, colorsOfVertices[firstVertice.colorId - 1], camera);
				}
				*/
			}
			else
			{
				// TRIANGLE RASTERIZATION
				if ((firstVertice.t == 1 && secondVertice.t == 1 && thirdVertice.t == 1) || !cullingEnabled)
					triangleRasterization(camera,
						ceil(firstVertice.x), ceil(firstVertice.y), colorsOfVertices[firstVertice.colorId - 1],
						ceil(secondVertice.x), ceil(secondVertice.y), colorsOfVertices[secondVertice.colorId - 1],
						ceil(thirdVertice.x), ceil(thirdVertice.y), colorsOfVertices[thirdVertice.colorId - 1]);
			}

			
		}
	}


	if (newVertices.size() == 0) return;
	for (int i = 0; i < newVertices.size() - 1; i += 2)
	{
		Vec4 firstVertice =  *newVertices[i];
		Vec4 secondVertice = *newVertices[i + 1];
		double x_0, y_0, x_1, y_1;
		if (firstVertice.t == 1 && secondVertice.t == 1)
		{
			x_0 = firstVertice.x;
			y_0 = firstVertice.y;
			x_1 = secondVertice.x;
			y_1 = secondVertice.y;
			lineRasterization(x_0, y_0, colorsOfVertices[firstVertice.colorId - 1], x_1, y_1, colorsOfVertices[secondVertice.colorId - 1], camera);
			
		}
	}
}


/********************************************
 * HELPER FUNCTIONS FOR CLIPPING			*
 ********************************************/
void Scene::clipColors(Vec4* v0, Vec4* v1, Vec4* v_clipped)
{
	Color* v0_color = colorsOfVertices[v0->colorId - 1];
	Color* v1_color = colorsOfVertices[v1->colorId - 1];
	double v0_color_scale, v1_color_scale;

	v0_color_scale = abs((v_clipped->x - v1->x) / (v1->x - v0->x));
	v1_color_scale = abs((v0->x - v_clipped->x) / (v1->x - v0->x));
	Color* color_v0_clipped = new Color(
		v0_color_scale * v0_color->r + v1_color_scale * v1_color->r,
		v0_color_scale * v0_color->g + v1_color_scale * v1_color->g,
		v0_color_scale * v0_color->b + v1_color_scale * v1_color->b);
	colorsOfVertices.push_back(color_v0_clipped);
	v_clipped->colorId = colorsOfVertices.size();
}

bool Scene::isVisible(double d, double num, double* tEnter, double* tLeave)
{
	double t = 0;
	if (d > 0)
	{
		t = num / d;
		if (t > * tLeave)
		{
			return false;
		}
		if (t > * tEnter)
		{
			*tEnter = t;
		}
	}
	else if (d < 0)
	{
		t = num / d;
		if (t < *tEnter)
		{
			return false;
		}
		if (t < *tLeave)
		{
			*tLeave = t;
		}
	}
	else if (num > 0)
		return  false;
	return  true;
}

void Scene::clipping(Vec4* v0, Vec4* v1, Vec4* v0_clipped, Vec4* v1_clipped, Camera* camera)
{
	double* tEnter = new double;
	double* tLeave = new double;
	*tLeave = 1;
	*tEnter = 0;

	const double dx = v1->x - v0->x;
	const double dy = v1->y - v0->y;
	const double dz = v1->z - v0->z;
		  
	const double xmin = -1;
	const double xmax = 1;
	const double ymin = -1;
	const double ymax = 1;
	const double zmin = -1;
	const double zmax = 1;

	if (isVisible(dx, xmin - v0->x, tEnter, tLeave))
	{
		if (isVisible(-dx, v0->x - xmax, tEnter, tLeave))
		{
			if (isVisible(dy, ymin - v0->y, tEnter, tLeave))
			{
				if (isVisible(-dy, v0->y - ymax, tEnter, tLeave))
				{
					if (isVisible(dz, zmin - v0->z, tEnter, tLeave))
					{
						if (isVisible(-dz, v0->z - zmax, tEnter, tLeave))
						{
							if (*tLeave < 1)
							{
								v1_clipped->x = v0->x + *tLeave * dx;
								v1_clipped->y = v0->y + *tLeave * dy;
								v1_clipped->z = v0->z + *tLeave * dz;
							}
							if (*tEnter > 0)
							{
								v0_clipped->x = v0->x + *tEnter * dx;
								v0_clipped->y = v0->y + *tEnter * dy;
								v0_clipped->z = v0->z + *tEnter * dz;
							}
							clipColors(v0, v1, v0_clipped);
							clipColors(v0, v1, v1_clipped);
						}

					}

				}
			}
		}
	}

}

void Scene::clippingModels(Camera* camera)
{
	for (int i = 0; i < models.size(); i++)
	{
		for (int j = 0; j < models[i]->numberOfTriangles; j++)
		{
			Vec4* v0 = verticesAssembled[i][j * 3];
			Vec4* v1 = verticesAssembled[i][j * 3 + 1];
			Vec4* v2 = verticesAssembled[i][j * 3 + 2];

			if (models[i]->type == 0)
			{
				Vec4* v0_c1 = new Vec4(*v0);
				Vec4* v1_c1 = new Vec4(*v1);
				Vec4* v2_c2 = new Vec4(*v2);
				Vec4* v0_c3 = new Vec4(*v0);
				Vec4* v1_c2 = new Vec4(*v1);
				Vec4* v2_c3 = new Vec4(*v2);
				
				clipping(v0, v1, v0_c1, v1_c1, camera);
				clipping(v1, v2, v1_c2, v2_c2, camera);
				clipping(v2, v0, v2_c3, v0_c3, camera);

				v0->t = -1;
				v1->t = -1;
				v2->t = -1;

				newVertices.push_back(v0_c1);
				newVertices.push_back(v1_c1);
				newVertices.push_back(v1_c2);
				newVertices.push_back(v2_c2);
				newVertices.push_back(v2_c3);
				newVertices.push_back(v0_c3);
			}
		}
	}
}


/********************************************
 * HELPER FUNCTIONS FOR CULLING				*
 ********************************************/
void Scene::backfaceCulling(Camera* camera)
{
	
	for (int i = 0; i < models.size(); i++)
	{
		for (int j = 0; j < models[i]->numberOfTriangles; j++)
		{
			const Vec4 v1 = *verticesAssembled[i][j * 3];
			const Vec4 v2 = *verticesAssembled[i][j * 3 + 1];
			const Vec4 v3 = *verticesAssembled[i][j * 3 + 2];

			const Vec3 v1_3 = *getVector3(v1);
			const Vec3 v2_3 = *getVector3(v2);
			const Vec3 v3_3 = *getVector3(v3);
			const Vec3 normal = crossProductVec3(subtractVec3(v2_3, v1_3), subtractVec3(v3_3, v1_3));
			Vec3 v = subtractVec3(v1_3, camera->v);

			if (dotProductVec3(normal, v) > 0)
			{
				verticesAssembled[i][j * 3]->t = 1;
				verticesAssembled[i][j * 3 + 1]->t = 1;
				verticesAssembled[i][j * 3 + 2]->t = 1;
			}
			else
			{
				verticesAssembled[i][j * 3]->t = -1;
				verticesAssembled[i][j * 3 + 1]->t = -1;
				verticesAssembled[i][j * 3 + 2]->t = -1;
			}
		}
	}
}


/********************************************
 * HELPER FUNCTIONS FOR TRANSFORM MATRICES	*
 ********************************************/
Matrix4 Scene::getCameraTransformMatrix(Camera* camera)
{
	Matrix4 cameraTranslate = getIdentityMatrix();
	
	cameraTranslate.val[0][3] = -camera->pos.x;
	cameraTranslate.val[1][3] = -camera->pos.y;
	cameraTranslate.val[2][3] = -camera->pos.z;
	
	Matrix4 cameraRotation = getIdentityMatrix();

	cameraRotation.val[0][0] = camera->u.x;
	cameraRotation.val[0][1] = camera->u.y;
	cameraRotation.val[0][2] = camera->u.z;

	cameraRotation.val[1][0] = camera->v.x;
	cameraRotation.val[1][1] = camera->v.y;
	cameraRotation.val[1][2] = camera->v.z;

	cameraRotation.val[2][0] = camera->w.x;
	cameraRotation.val[2][1] = camera->w.y;
	cameraRotation.val[2][2] = camera->w.z;

	return multiplyMatrixWithMatrix(cameraRotation, cameraTranslate);

}

Matrix4 Scene::getOrthographicProjectionMatrix(Camera* camera)
{
	Matrix4 projectionMatrix = getIdentityMatrix();
	const double xDistance = camera->right - camera->left;
	const double yDistance = camera->top - camera->bottom;
	const double zDistance = camera->far - camera->near;

	projectionMatrix.val[0][0] = 2 / xDistance;
	projectionMatrix.val[1][1] = 2 / yDistance;
	projectionMatrix.val[2][2] = -2 / zDistance;

	projectionMatrix.val[0][3] = -(camera->right + camera->left) / xDistance;
	projectionMatrix.val[1][3] = -(camera->top + camera->bottom) / yDistance;
	projectionMatrix.val[2][3] = -(camera->far + camera->near) / zDistance;

	return  projectionMatrix;
}

Matrix4 Scene::getPerspectiveProjectionMatrix(Camera* camera)
{
	Matrix4 projectionMatrix = getIdentityMatrix();
	const double xDistance = camera->right - camera->left;
	const double yDistance = camera->top - camera->bottom;
	const double zDistance = camera->far - camera->near;
	const double two_n = 2 * camera->near;

	projectionMatrix.val[0][0] = two_n / xDistance;
	projectionMatrix.val[1][1] = two_n / yDistance;
	projectionMatrix.val[2][2] = -(camera->far + camera->near) / zDistance;

	projectionMatrix.val[0][2] = (camera->right + camera->left) / xDistance;
	projectionMatrix.val[1][2] = (camera->top + camera->bottom) / yDistance;

	projectionMatrix.val[2][3] = -2 * (camera->far * camera->near) / zDistance;

	projectionMatrix.val[3][2] = -1;
	projectionMatrix.val[3][3] = 0;
	return  projectionMatrix;
}

Matrix4 Scene::getViewportProjectionMatrix(Camera* camera)
{
	Matrix4 projectionMatix = getIdentityMatrix();
	projectionMatix.val[0][0] = camera->horRes / 2.0;
	projectionMatix.val[1][1] = camera->verRes / 2.0;
	projectionMatix.val[2][2] = 0.5;

	projectionMatix.val[0][3] = (camera->horRes - 1) / 2.0;
	projectionMatix.val[1][3] = (camera->verRes - 1) / 2.0;
	projectionMatix.val[2][3] = 0.5;

	return projectionMatix;
}

Matrix4 Scene::getTranslationMatrix(Translation* translation)
{
	Matrix4 translationMatrix = getIdentityMatrix();

	translationMatrix.val[0][3] = translation->tx;
	translationMatrix.val[1][3] = translation->ty;
	translationMatrix.val[2][3] = translation->tz;

	return translationMatrix;
}

Matrix4 Scene::getScalingMatrix(Scaling* scaling)
{
	Matrix4 scalingMatrix = getIdentityMatrix();

	scalingMatrix.val[0][0] = scaling->sx;
	scalingMatrix.val[1][1] = scaling->sy;
	scalingMatrix.val[2][2] = scaling->sz;

	return scalingMatrix;
}

Matrix4 Scene::getRotationMatrix(Rotation* rotation)
{
	const double cost = cos(rotation->angle * M_PI / 180);
	const double sint = sin(rotation->angle * M_PI / 180);
	Vec3* u = new Vec3(rotation->ux, rotation->uy, rotation->uz, -1);
	*u = normalizeVec3(*u);
	Matrix4 rotationMatrix = getIdentityMatrix();

	rotationMatrix.val[0][0] = cost + u->x * u->x * (1 - cost);
	rotationMatrix.val[0][1] = u->x * u->y * (1 - cost) - u->z * sint;
	rotationMatrix.val[0][2] = u->x * u->z * (1 - cost) + u->y * sint;

	rotationMatrix.val[1][0] = u->y * u->x * (1 - cost) + u->z * sint;
	rotationMatrix.val[1][1] = cost + u->y * u->y * (1 - cost);
	rotationMatrix.val[1][2] = u->y * u->z * (1 - cost) - u->x * sint;

	rotationMatrix.val[2][0] = u->z * u->x * (1 - cost) - u->y * sint;
	rotationMatrix.val[2][1] = u->z * u->y * (1 - cost) + u->x * sint;
	rotationMatrix.val[2][2] = cost + u->z * u->z * (1 - cost);

	return rotationMatrix;
}




/*
	Transformations, clipping, culling, rasterization are done here.
	You can define helper functions inside Scene class implementation.
*/
void Scene::forwardRenderingPipeline(Camera* camera)
{
	newVertices.clear();
	verticesAssembled.clear();

	for (int i = 0; i < models.size(); i++)
	{
		vector<Vec4*> modelVertices;
		for (int j = 0; j < models[i]->numberOfTriangles; j++)
		{
			Vec3 v1 = *vertices[models[i]->triangles[j].getFirstVertexId() - 1];
			Vec3 v2 = *vertices[models[i]->triangles[j].getSecondVertexId() - 1];
			Vec3 v3 = *vertices[models[i]->triangles[j].getThirdVertexId() - 1];
			modelVertices.push_back(new Vec4(v1.x, v1.y, v1.z, 1, v1.colorId));
			modelVertices.push_back(new Vec4(v2.x, v2.y, v2.z, 1, v2.colorId));
			modelVertices.push_back(new Vec4(v3.x, v3.y, v3.z, 1, v3.colorId));
		}

		verticesAssembled.push_back(modelVertices);
	}
	const Matrix4 cameraMatrix = getCameraTransformMatrix(camera);

	Matrix4 projectionMatrix;
	if (projectionType == 0) // orthographic
	{
		projectionMatrix = getOrthographicProjectionMatrix(camera);
	}
	else  // perspective
	{
		projectionMatrix = getPerspectiveProjectionMatrix(camera);
	}
	const Matrix4 viewPortMatrix = getViewportProjectionMatrix(camera);


	modelingTransformation();
	transformation(cameraMatrix, camera);
	transformation(projectionMatrix, camera);

	setModelw();


	if (cullingEnabled)
		backfaceCulling(camera);
	clippingModels(camera);

	transformation(viewPortMatrix, camera);
	rasterization(camera);


}

/*
	Parses XML file
*/
Scene::Scene(const char* xmlPath)
{
	const char* str;
	XMLDocument xmlDoc;
	XMLElement* pElement;

	xmlDoc.LoadFile(xmlPath);

	XMLNode* pRoot = xmlDoc.FirstChild();

	// read background color
	pElement = pRoot->FirstChildElement("BackgroundColor");
	str = pElement->GetText();
	sscanf(str, "%lf %lf %lf", &backgroundColor.r, &backgroundColor.g, &backgroundColor.b);

	// read culling
	pElement = pRoot->FirstChildElement("Culling");
	if (pElement != NULL)
		pElement->QueryBoolText(&cullingEnabled);

	// read projection type
	pElement = pRoot->FirstChildElement("ProjectionType");
	if (pElement != NULL)
		pElement->QueryIntText(&projectionType);

	// read cameras
	pElement = pRoot->FirstChildElement("Cameras");
	XMLElement* pCamera = pElement->FirstChildElement("Camera");
	XMLElement* camElement;
	while (pCamera != NULL)
	{
		Camera* cam = new Camera();

		pCamera->QueryIntAttribute("id", &cam->cameraId);

		camElement = pCamera->FirstChildElement("Position");
		str = camElement->GetText();
		sscanf(str, "%lf %lf %lf", &cam->pos.x, &cam->pos.y, &cam->pos.z);

		camElement = pCamera->FirstChildElement("Gaze");
		str = camElement->GetText();
		sscanf(str, "%lf %lf %lf", &cam->gaze.x, &cam->gaze.y, &cam->gaze.z);

		camElement = pCamera->FirstChildElement("Up");
		str = camElement->GetText();
		sscanf(str, "%lf %lf %lf", &cam->v.x, &cam->v.y, &cam->v.z);

		cam->gaze = normalizeVec3(cam->gaze);
		cam->u = crossProductVec3(cam->gaze, cam->v);
		cam->u = normalizeVec3(cam->u);

		cam->w = inverseVec3(cam->gaze);
		cam->v = crossProductVec3(cam->u, cam->gaze);
		cam->v = normalizeVec3(cam->v);

		camElement = pCamera->FirstChildElement("ImagePlane");
		str = camElement->GetText();
		sscanf(str, "%lf %lf %lf %lf %lf %lf %d %d",
			&cam->left, &cam->right, &cam->bottom, &cam->top,
			&cam->near, &cam->far, &cam->horRes, &cam->verRes);

		camElement = pCamera->FirstChildElement("OutputName");
		str = camElement->GetText();
		cam->outputFileName = string(str);

		cameras.push_back(cam);

		pCamera = pCamera->NextSiblingElement("Camera");
	}

	// read vertices
	pElement = pRoot->FirstChildElement("Vertices");
	XMLElement* pVertex = pElement->FirstChildElement("Vertex");
	int vertexId = 1;

	while (pVertex != NULL)
	{
		Vec3* vertex = new Vec3();
		Color* color = new Color();

		vertex->colorId = vertexId;

		str = pVertex->Attribute("position");
		sscanf(str, "%lf %lf %lf", &vertex->x, &vertex->y, &vertex->z);

		str = pVertex->Attribute("color");
		sscanf(str, "%lf %lf %lf", &color->r, &color->g, &color->b);

		vertices.push_back(vertex);
		colorsOfVertices.push_back(color);

		pVertex = pVertex->NextSiblingElement("Vertex");

		vertexId++;
	}

	// read translations
	pElement = pRoot->FirstChildElement("Translations");
	XMLElement* pTranslation = pElement->FirstChildElement("Translation");
	while (pTranslation != NULL)
	{
		Translation* translation = new Translation();

		pTranslation->QueryIntAttribute("id", &translation->translationId);

		str = pTranslation->Attribute("value");
		sscanf(str, "%lf %lf %lf", &translation->tx, &translation->ty, &translation->tz);

		translations.push_back(translation);

		pTranslation = pTranslation->NextSiblingElement("Translation");
	}

	// read scalings
	pElement = pRoot->FirstChildElement("Scalings");
	XMLElement* pScaling = pElement->FirstChildElement("Scaling");
	while (pScaling != NULL)
	{
		Scaling* scaling = new Scaling();

		pScaling->QueryIntAttribute("id", &scaling->scalingId);
		str = pScaling->Attribute("value");
		sscanf(str, "%lf %lf %lf", &scaling->sx, &scaling->sy, &scaling->sz);

		scalings.push_back(scaling);

		pScaling = pScaling->NextSiblingElement("Scaling");
	}

	// read rotations
	pElement = pRoot->FirstChildElement("Rotations");
	XMLElement* pRotation = pElement->FirstChildElement("Rotation");
	while (pRotation != NULL)
	{
		Rotation* rotation = new Rotation();

		pRotation->QueryIntAttribute("id", &rotation->rotationId);
		str = pRotation->Attribute("value");
		sscanf(str, "%lf %lf %lf %lf", &rotation->angle, &rotation->ux, &rotation->uy, &rotation->uz);

		rotations.push_back(rotation);

		pRotation = pRotation->NextSiblingElement("Rotation");
	}

	// read models
	pElement = pRoot->FirstChildElement("Models");

	XMLElement* pModel = pElement->FirstChildElement("Model");
	XMLElement* modelElement;
	while (pModel != NULL)
	{
		Model* model = new Model();

		pModel->QueryIntAttribute("id", &model->modelId);
		pModel->QueryIntAttribute("type", &model->type);

		// read model transformations
		XMLElement* pTransformations = pModel->FirstChildElement("Transformations");
		XMLElement* pTransformation = pTransformations->FirstChildElement("Transformation");

		pTransformations->QueryIntAttribute("count", &model->numberOfTransformations);

		while (pTransformation != NULL)
		{
			char transformationType;
			int transformationId;

			str = pTransformation->GetText();
			sscanf(str, "%c %d", &transformationType, &transformationId);

			model->transformationTypes.push_back(transformationType);
			model->transformationIds.push_back(transformationId);

			pTransformation = pTransformation->NextSiblingElement("Transformation");
		}

		// read model triangles
		XMLElement* pTriangles = pModel->FirstChildElement("Triangles");
		XMLElement* pTriangle = pTriangles->FirstChildElement("Triangle");

		pTriangles->QueryIntAttribute("count", &model->numberOfTriangles);

		while (pTriangle != NULL)
		{
			int v1, v2, v3;

			str = pTriangle->GetText();
			sscanf(str, "%d %d %d", &v1, &v2, &v3);

			model->triangles.push_back(Triangle(v1, v2, v3));

			pTriangle = pTriangle->NextSiblingElement("Triangle");
		}

		models.push_back(model);

		pModel = pModel->NextSiblingElement("Model");
	}
}

/*
	Initializes image with background color
*/
void Scene::initializeImage(Camera* camera)
{
	if (this->image.empty())
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			vector<Color> rowOfColors;

			for (int j = 0; j < camera->verRes; j++)
			{
				rowOfColors.push_back(this->backgroundColor);
			}

			this->image.push_back(rowOfColors);
		}
	}
	// if image is filled before, just change color rgb values with the background color
	else
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			for (int j = 0; j < camera->verRes; j++)
			{
				this->image[i][j].r = this->backgroundColor.r;
				this->image[i][j].g = this->backgroundColor.g;
				this->image[i][j].b = this->backgroundColor.b;
			}
		}
	}
}



/*
	If given value is less than 0, converts value to 0.
	If given value is more than 255, converts value to 255.
	Otherwise returns value itself.
*/
int Scene::makeBetweenZeroAnd255(double value)
{
	if (value >= 255.0)
		return 255;
	if (value <= 0.0)
		return 0;
	return (int)(value);
}

/*
	Writes contents of image (Color**) into a PPM file.
*/
void Scene::writeImageToPPMFile(Camera* camera)
{
	ofstream fout;

	fout.open(camera->outputFileName.c_str());

	fout << "P3" << endl;
	fout << "# " << camera->outputFileName << endl;
	fout << camera->horRes << " " << camera->verRes << endl;
	fout << "255" << endl;

	for (int j = camera->verRes - 1; j >= 0; j--)
	{
		for (int i = 0; i < camera->horRes; i++)
		{
			fout << makeBetweenZeroAnd255(this->image[i][j].r) << " "
				<< makeBetweenZeroAnd255(this->image[i][j].g) << " "
				<< makeBetweenZeroAnd255(this->image[i][j].b) << " ";
		}
		fout << endl;
	}
	fout.close();
}

/*
	Converts PPM image in given path to PNG file, by calling ImageMagick's 'convert' command.
	os_type == 1 		-> Ubuntu
	os_type == 2 		-> Windows
	os_type == other	-> No conversion
*/
void Scene::convertPPMToPNG(string ppmFileName, int osType)
{
	string command;

	// call command on Ubuntu
	if (osType == 1)
	{
		command = "convert " + ppmFileName + " " + ppmFileName + ".png";
		system(command.c_str());
	}

	// call command on Windows
	else if (osType == 2)
	{
		command = "magick convert " + ppmFileName + " " + ppmFileName + ".png";
		system(command.c_str());
	}

	// default action - don't do conversion
	else
	{
	}
}