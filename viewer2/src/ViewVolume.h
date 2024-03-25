#ifndef ViewVolume_h
#define ViewVolume_h 1

#include <fstream>
#include <string.h>
#include "FaceContainer.h"
#include "Hexahedron.h"
#include "Prism.h"
#include "Pyramid.h"
#include "PyramidBase5.h"
#include "Tetrahedron.h"
#include "Element.h"
#include "Point3D.h"

using namespace std;
using Clobscode::Point3D;

class ViewVolume{
	
  public:
	
	ViewVolume();
	
	virtual ~ViewVolume();
	
	virtual bool execute(char *input, const bool &edge, const unsigned int &type);
	
  protected:
	
	virtual bool yesOrNo(string &question);
	
	virtual bool askRemove(FaceContainer &tfc);
    
	virtual void removeFaces(int &axe, double start, double end, FaceContainer &tfc);
	
	virtual void clean(list<Point3D> &tpoints, list<Face> &tfaces, FaceContainer &tfc);
	
	virtual void cleanUnusedNodes(list<Point3D> &tpoints, list<Face> &tfaces,
                                  FaceContainer &tfc);
	
	virtual int askAxe();
	
	virtual double readDouble();
	
	virtual bool newTmpFile(string &stmp, ofstream &output);
	
	virtual bool readData(char *name, const unsigned int &type);
	
	virtual void headGeom(ofstream &output, int np, int nf);
	
	virtual void dataGeom(ofstream &output, list<Point3D> &tpoints, list<Face> &tfaces);
	
	virtual void printBbox();
    
    virtual bool readM3d(string name);
	
    virtual bool readMixed(string name);
    
    virtual bool readVtk(string name);
    
  protected:
	
	char op;
	vector<Point3D> points;
	FaceContainer all_faces;
	vector<Element *> elements;
	Point3D min, max;
};

#endif
