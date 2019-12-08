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

using namespace tinyxml2;
using namespace std;


void Scene::lineRasterization(int x_0, int y_0, Color* c_0, int x_1, int y_1, Color* c_1)
{
	int y = y_0;
	int distance = 2 * (y_0 - y_1) + (x_1 - x_0);
	Vec3* c_0v = new Vec3(c_0->r, c_0->g, c_0->b, 0);
	cout << x_1 << " - " << y_1<<endl;
	Vec3* c_1v = new Vec3(c_1->r, c_1->g, c_1->b, 0);
	Vec3* c = new Vec3(*c_0v);
	Vec3 d_c = addVec3(*c_1v, multiplyVec3WithScalar(*c_0v, -1));
	for (int x = x_0; x < x_1; x++)
	{
		this->image[x][y] = Color(round(c->x), round(c->x), round(c->x));

		if (distance < 0)
		{
			y += 1;
			distance += 2 * ((y_0 - y_1) + (x_1 - x_0));
		}
		else
		{
			distance += 2 * (y_0 - y_1);
		}
		*c = addVec3(*c, d_c);
	}
}

void Scene::rasterization(Camera* camera)
{
	for (int i = 0; i < models.size(); i++)
	{
		for (int j = 0; j < models[i]->numberOfTriangles; j++)
		{


			Vec3 firstVertice = *vertices[models[i]->triangles[j].getFirstVertexId() - 1];
			Vec3 secondVertice = *vertices[models[i]->triangles[j].getSecondVertexId() - 1];
			Vec3 thirdVertice = *vertices[models[i]->triangles[j].getThirdVertexId() - 1];

			// LINE RASTERIZATION

			int x_0, y_0, x_1, y_1;
			// vertice1 -> vertice2
			x_0 = round(firstVertice.x);
			y_0 = round(firstVertice.y);
			x_1 = round(secondVertice.x);
			y_1 = round(secondVertice.y);
			if (x_0 <= camera->horRes && y_0 <= camera->verRes && x_1 <= camera->horRes && y_1 <= camera->verRes)
			{
				if (x_0 >= 0 && y_0 >= 0 && x_1 >= 0 && y_1>= 0)
				{
					lineRasterization(x_0, y_0, colorsOfVertices[firstVertice.colorId - 1], x_1, y_1, colorsOfVertices[secondVertice.colorId - 1]);
				}
			}

			// vertice2 -> vertice3
			x_0 = round(secondVertice.x);
			y_0 = round(secondVertice.y);
			x_1 = round(thirdVertice.x);
			y_1 = round(thirdVertice.y);
			if (x_0 <= camera->horRes && y_0 <= camera->verRes && x_1 <= camera->horRes && y_1 <= camera->verRes)
			{
				if (x_0 >= 0 && y_0 >= 0 && x_1 >= 0 && y_1 >= 0)
				{
					lineRasterization(x_0, y_0, colorsOfVertices[secondVertice.colorId - 1], x_1, y_1, colorsOfVertices[thirdVertice.colorId - 1]);
				}
			}
			
			// vertice3 -> vertice1
			x_0 = round(thirdVertice.x);
			y_0 = round(thirdVertice.y);
			x_1 = round(firstVertice.x);
			y_1 = round(firstVertice.y);
			if (x_0 <= camera->horRes && y_0 <= camera->verRes && x_1 <= camera->horRes && y_1 <= camera->verRes)
			{
				if (x_0 >= 0 && y_0 >= 0 && x_1 >= 0 && y_1 >= 0)
				{
					lineRasterization(x_0, y_0, colorsOfVertices[thirdVertice.colorId - 1], x_1, y_1, colorsOfVertices[firstVertice.colorId - 1]);
				}
			}

			// TRIANGLE RASTERIZATION
			// TODO: Implement
		}
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
				cout << "Translating..." << endl;
				Matrix4 translationMatrix = getTranslationMatrix(translations[models[i]->transformationIds[j]]);
				for (int k = 0; k < models[i]->numberOfTriangles; k++)
				{
					Vec3 firstVertice = *vertices[models[i]->triangles[j].getFirstVertexId() - 1];
					Vec3 secondVertice = *vertices[models[i]->triangles[j].getSecondVertexId() - 1];
					Vec3 thirdVertice = *vertices[models[i]->triangles[j].getThirdVertexId() - 1];

					Vec4 first = multiplyMatrixWithVec4(translationMatrix, *getVector4(firstVertice));
					Vec4 second = multiplyMatrixWithVec4(translationMatrix, *getVector4(secondVertice));
					Vec4 third = multiplyMatrixWithVec4(translationMatrix, *getVector4(thirdVertice));
				}
			}
		}
		
	}
}

void Scene::transformation(Matrix4 transformationMatrix, Camera* camera)
{
	for (int i = 0; i < models.size(); i++)
	{
		for (int j = 0; j < models[i]->numberOfTriangles; j++)
		{
			Vec3 firstVertice = *vertices[models[i]->triangles[j].getFirstVertexId() - 1];
			Vec3 secondVertice = *vertices[models[i]->triangles[j].getSecondVertexId() - 1];
			Vec3 thirdVertice = *vertices[models[i]->triangles[j].getThirdVertexId() - 1];


			Vec4 first = multiplyMatrixWithVec4(transformationMatrix, *getVector4(firstVertice));
			Vec4 second = multiplyMatrixWithVec4(transformationMatrix, *getVector4(secondVertice));
			Vec4 third = multiplyMatrixWithVec4(transformationMatrix, *getVector4(thirdVertice));
			
			// TODO: set new vertices
		}
	}
}


void Scene::modelTransformation(Matrix4 worldMatrix,Camera *camera)
{
	for (int i = 0 ; i < models.size();i++)
	{
		for (int j = 0 ; j < models[i]->numberOfTriangles;j++)
		{

		
			Vec3 firstVertice =	*vertices[models[i]->triangles[j].getFirstVertexId() - 1];
			Vec3 secondVertice = *vertices[models[i]->triangles[j].getSecondVertexId() - 1];
			Vec3 thirdVertice = *vertices[models[i]->triangles[j].getThirdVertexId() - 1];


			Vec4 first = multiplyMatrixWithVec4(worldMatrix, *getVector4(firstVertice));
			Vec4 second = multiplyMatrixWithVec4(worldMatrix, *getVector4(secondVertice));
			Vec4 third = multiplyMatrixWithVec4(worldMatrix, *getVector4(thirdVertice));

			Vec3* first1 = getVector3(first);
			Vec3* second1 = getVector3(second);
			Vec3* third1 = getVector3(third);
			clipping(first1, second1, camera);
			clipping(second1, third1, camera);
			clipping(third1, first1, camera);
		}
	}
}
Vec3* Scene::getVector3(Vec4 vector)
{
	return new Vec3(vector.x, vector.y, vector.z, -1);
}
void Scene::clipping(Vec3 *v0, Vec3 *v1, Camera* camera)
{
	float *tEnter = new float;
	float *tLeave = new float;
	*tLeave = 1;
	*tEnter = 0;
	bool visible = false;

	float dx = v1->x - v0->x;
	float dy = v1->y - v0->y;
	float dz = v1->z - v0->z;

	float xmin = camera->left;
	float xmax = camera->right;
	float ymin = camera->bottom;
	float ymax = camera->top;
	float zmin = camera->near;
	float zmax = camera->far;

	visible = isVisible(dx, xmin -v0->x, tEnter, tLeave);
	visible = isVisible(-dx, v0->x - xmax, tEnter, tLeave);

	visible = isVisible(dy, ymin - v0->y, tEnter, tLeave);
	visible = isVisible(-dy, v0->y - ymax, tEnter, tLeave);

	visible = isVisible(dz, zmin - v0->z, tEnter, tLeave);
	visible = isVisible(-dz, v0->z - zmax, tEnter, tLeave);
	if(visible)
	{
		if(*tLeave<1)
		{
			v1->x = v0->x + dx * (*tLeave);
			v1->y = v0->y + dy * (*tLeave);
			v1->z = v0->z + dz * (*tLeave);
		}
		if (*tEnter < 1)
		{
			v0->x = v0->x + dx * (*tEnter);
			v0->y = v0->y + dy * (*tEnter);
			v0->z = v0->z + dz * (*tEnter);
		}
	}

	
}

bool Scene::isVisible(float d, float num, float *tEnter, float *tLeave)
{
	float t = 0;
	if(d>0)
	{
		t = num / d;
	}
	if(t> *tLeave)
	{
		return false;
	}
	if (t > * tEnter)
	{
		*tEnter = t;
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


Vec3 Scene::getPointAtt(Vec3 v0, Vec3 v1, float t)
{
	Vec3 p = addVec3(v0, multiplyVec3WithScalar(subtractVec3(v0, v1), t));
	return p;
}


/*
	Transformations, clipping, culling, rasterization are done here.
	You can define helper functions inside Scene class implementation.
*/
void Scene::forwardRenderingPipeline(Camera *camera)
{
	Matrix4 cameraMatrix = getCameraTransformMatrix(camera);

	Matrix4 projectionMatrix;
	if(projectionType == 0) // orthographic
	{
		projectionMatrix = getOrthographicProjectionMatrix(camera);
	}
	else  // perspective
	{
		projectionMatrix = getPerspectiveProjectionMatrix(camera);
	}
	Matrix4 viewPortMatrix = getViewportProjectionMatrix(camera);


	modelingTransformation();
	//transformation(cameraMatrix, camera);
	//transformation(projectionMatrix, camera);
	//transformation(viewPortMatrix, camera);

	rasterization(camera);
	//Matrix4 worldMatrix = multiplyMatrixWithMatrix(viewPortMatrix, multiplyMatrixWithMatrix(projectionMatrix, cameraMatrix));
	//modelTransformation(worldMatrix,camera);

	
}
Vec4 * Scene::getVector4(Vec3 vector)
{
	Vec4* vec = new Vec4(vector.x, vector.y, vector.z, 1, 0);
	return  vec;
}
Matrix4 Scene::getCameraTransformMatrix(Camera* camera)
{
	Matrix4 cameraTranslate = getIdentityMatrix();
	Vec4 *cameraPosition = new Vec4(camera->pos.x, camera->pos.y, camera->pos.z, 1, 0);
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
	Matrix4 projectionMatix = getIdentityMatrix();
	float xDistance = camera->right - camera->left;
	float yDistance = camera->top - camera->bottom;
	float zDistance = camera->far - camera->near;
	//TODO: Handle
	if(xDistance == 0 || yDistance == 0 ||zDistance ==0)
	{
		cout << "division by zero error" << endl;
	}
	projectionMatix.val[0][0] = 2/xDistance;
	projectionMatix.val[1][1] = 2/yDistance;
	projectionMatix.val[2][2] = -2/zDistance;

	projectionMatix.val[0][3] = -(camera->right+camera->left) / xDistance;
	projectionMatix.val[1][3] = -(camera->top+camera->bottom) / yDistance;
	projectionMatix.val[2][3] = -(camera->far+camera->near) / zDistance;

	return  projectionMatix;
}
Matrix4 Scene::getPerspectiveProjectionMatrix(Camera* camera)
{
	Matrix4 projectionMatix = getIdentityMatrix();
	float xDistance = camera->right - camera->left;
	float yDistance = camera->top - camera->bottom;
	float zDistance = camera->far - camera->near;
	float two_n = 2 * camera->near;
	
	projectionMatix.val[0][0] = two_n / xDistance;
	projectionMatix.val[1][1] = two_n / yDistance;
	projectionMatix.val[2][2] = -(camera->far + camera->near) / zDistance;

	projectionMatix.val[0][2] = (camera->right + camera->left) / xDistance;
	projectionMatix.val[1][2] = (camera->top + camera->bottom) / yDistance;
	
	projectionMatix.val[2][3] = -2*(camera->far * camera->near) / zDistance;
	
	projectionMatix.val[3][2] = -1;
	projectionMatix.val[3][3] = 0;
	return  projectionMatix;
}
Matrix4 Scene::getViewportProjectionMatrix(Camera* camera)
{
	Matrix4 projectionMatix = getIdentityMatrix();
	projectionMatix.val[0][0] = camera->horRes / 2;
	projectionMatix.val[1][1] = camera->verRes / 2;
	projectionMatix.val[2][2] = 0.5;

	projectionMatix.val[0][3] = (camera->horRes -1) / 2;
	projectionMatix.val[1][3] = (camera->verRes - 1) / 2;
	projectionMatix.val[2][3] = 0.5;

	
	projectionMatix.val[3][3] = 0;

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
/*
	Parses XML file
*/
Scene::Scene(const char *xmlPath)
{
	const char *str;
	XMLDocument xmlDoc;
	XMLElement *pElement;

	xmlDoc.LoadFile(xmlPath);

	XMLNode *pRoot = xmlDoc.FirstChild();

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
	XMLElement *pCamera = pElement->FirstChildElement("Camera");
	XMLElement *camElement;
	while (pCamera != NULL)
	{
		Camera *cam = new Camera();

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
	XMLElement *pVertex = pElement->FirstChildElement("Vertex");
	int vertexId = 1;

	while (pVertex != NULL)
	{
		Vec3 *vertex = new Vec3();
		Color *color = new Color();

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
	XMLElement *pTranslation = pElement->FirstChildElement("Translation");
	while (pTranslation != NULL)
	{
		Translation *translation = new Translation();

		pTranslation->QueryIntAttribute("id", &translation->translationId);

		str = pTranslation->Attribute("value");
		sscanf(str, "%lf %lf %lf", &translation->tx, &translation->ty, &translation->tz);

		translations.push_back(translation);

		pTranslation = pTranslation->NextSiblingElement("Translation");
	}

	// read scalings
	pElement = pRoot->FirstChildElement("Scalings");
	XMLElement *pScaling = pElement->FirstChildElement("Scaling");
	while (pScaling != NULL)
	{
		Scaling *scaling = new Scaling();

		pScaling->QueryIntAttribute("id", &scaling->scalingId);
		str = pScaling->Attribute("value");
		sscanf(str, "%lf %lf %lf", &scaling->sx, &scaling->sy, &scaling->sz);

		scalings.push_back(scaling);

		pScaling = pScaling->NextSiblingElement("Scaling");
	}

	// read rotations
	pElement = pRoot->FirstChildElement("Rotations");
	XMLElement *pRotation = pElement->FirstChildElement("Rotation");
	while (pRotation != NULL)
	{
		Rotation *rotation = new Rotation();

		pRotation->QueryIntAttribute("id", &rotation->rotationId);
		str = pRotation->Attribute("value");
		sscanf(str, "%lf %lf %lf %lf", &rotation->angle, &rotation->ux, &rotation->uy, &rotation->uz);

		rotations.push_back(rotation);

		pRotation = pRotation->NextSiblingElement("Rotation");
	}

	// read models
	pElement = pRoot->FirstChildElement("Models");

	XMLElement *pModel = pElement->FirstChildElement("Model");
	XMLElement *modelElement;
	while (pModel != NULL)
	{
		Model *model = new Model();

		pModel->QueryIntAttribute("id", &model->modelId);
		pModel->QueryIntAttribute("type", &model->type);

		// read model transformations
		XMLElement *pTransformations = pModel->FirstChildElement("Transformations");
		XMLElement *pTransformation = pTransformations->FirstChildElement("Transformation");

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
		XMLElement *pTriangles = pModel->FirstChildElement("Triangles");
		XMLElement *pTriangle = pTriangles->FirstChildElement("Triangle");

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
void Scene::initializeImage(Camera *camera)
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
void Scene::writeImageToPPMFile(Camera *camera)
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