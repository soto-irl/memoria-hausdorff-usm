#ifndef ViewSurface_h
#define ViewSurface_h 1

#include <fstream>
#include <string.h>
#include <vector>
#include "Point3D.h"

using namespace std;
using Clobscode::Point3D;

class ViewSurface{
	
  public:
	
	ViewSurface();
	
	virtual ~ViewSurface();
	
	virtual bool execute(char *input);
	
  protected:
	
	virtual bool yesOrNo(string &question);
	
	
	virtual bool askRemove(list<Point3D> &tpoints, list<vector<unsigned int> > &tfaces);
	
	virtual void removeFaces(int &axe, double start, double end,
							 list<Point3D> &tpoints, list<vector<unsigned int> > &tfaces);
	
    virtual void clean(list<Point3D> &tpoints, list<vector<unsigned int> > &tfaces);
	
	virtual int askAxe();
	
	virtual double readDouble();
	
	
	
	virtual bool newTmpFile(string &stmp,ofstream &output);
    
    virtual bool readMdl(string name);
	
	virtual void headGeom(ofstream &output, int np, int nf);
	
	virtual void dataGeom(ofstream &output, list<Point3D> &tpoints,
						  list<vector<unsigned int> > &tfaces);
	
    virtual void dataGeom(ofstream &output);
    
	virtual void printBbox(); 
	
  protected:
	
	char op;
	vector<Point3D> points;
    vector<vector<unsigned int> > faces;
	Point3D min, max;
    vector<bool> removedfaces;
};

#endif
