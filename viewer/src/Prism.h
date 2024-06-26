/**************************************************************************/
/*                                                                        */
/*                  Jens                                                  */
/*                                                                        */
/**************************************************************************/
/* Allows to display Jacobian Information of a mesh in .m3d format        */
/* Written by Claudio Lobos (clobos@inf.utfsm.cl) 2015					  */
/* UNIVERSIDAD TECNICA FEDERICO SANTA MARIA								  */
/**************************************************************************/

/*
 <JENS: this program computes the Element Normalized Scaled Jacobian (Jens)>
 
 Copyright (C) <2015>  <Claudio Lobos>
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/gpl.txt>
 */

#ifndef Prism_h
#define Prism_h 1

#include "Element.h"

using namespace std;

class Prism : public Element {
    
public:
    Prism(vector<unsigned int> &pts);
    
    virtual ~Prism();
    
    virtual double getAspectRatio(vector<Point3D> &points);
    
    virtual unsigned int numberOfFaces();
    
    virtual vector<unsigned int> getFacePoints(const unsigned int &i);

    
protected:
    
    virtual vector<unsigned int> getTetra(int index);
    
    virtual vector<vector<unsigned int> > getTetras();
    
};

inline unsigned int Prism::numberOfFaces() {
    return 5;
}


#endif

