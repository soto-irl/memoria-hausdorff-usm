#ifndef View_h
#define View_h 1

#include <fstream>
#include "../../geometry/FaceContainer.h"
#include "../../geometry/PointM3d.h"
#include "Element.h"
#include "../../read_write/ReadM3d.h"

using namespace std;

class View{
	
  public:

      View();

      virtual ~View();

      virtual bool execute(char *input);

  protected:

      virtual bool yesOrNo(string &question);

      virtual bool askRemove(list<Point> &tpoints, list<Face> &tfaces); 

      virtual void removeFaces(int &axe, double start, double end,
		     list<Point> &tpoints, list<Face> &tfaces);

      virtual void clean(list<Point> &tpoints, list<Face> &tfaces, 
		     FaceContainer *tfc);

      virtual int askAxe();
      
      virtual double readDouble();

      virtual bool newTmpFile(string &stmp,ofstream &output);

      virtual bool readData(char *name);

      virtual void headGeom(ofstream &output, int np, int nf);

      virtual void dataGeom(ofstream &output, list<Point> &tpoints, 
		     list<Face> &tfaces);

      virtual void printBbox(); 

  protected:

      char op;
      vector<Point> points;
      FaceContainer *fc;
      vector<EnhancedElement> elements;
      Point min, max;
};

#endif
