#include "ViewVolume.h"

//--------------------------------------------------------------------
//--------------------------------------------------------------------

ViewVolume::ViewVolume(){
	op='\0';
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------

ViewVolume::~ViewVolume(){
    for (unsigned int i=0; i<elements.size(); i++) {
        delete elements[i];
    }
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------
bool ViewVolume::readM3d(string name){
    
    
    char word [256];
    unsigned int cant,type;
    double x,y,z;
    
    FILE *file = fopen(name.c_str(),"r");
    
    while(true){
        if(!fscanf(file,"%s",word)){
            fclose(file);
            return false;
        }
        if(!strcmp(word,"ARRAY1<STRING>]\0"))
            break;
    }
    
    fscanf(file,"%u",&cant);
    points.reserve(cant);
    
    
    if(cant<=0){
        cerr << "warning: no nodes were found\n";
        fclose(file);
        return true;
    }
    
    //read the first point (to initialize the bounds)
	fscanf(file,"%u",&type);
	fscanf(file,"%s",word);
	x=atof(word);
	fscanf(file,"%s",word);
	y=atof(word);
	fscanf(file,"%s",word);
	z=atof(word);
	Point3D p (x,y,z);

	if(type == 3 || type == 2){
		fscanf(file,"%s",word);
		fscanf(file,"%s",word);
		fscanf(file,"%s",word);
	}
    
	points.push_back(p);
	
	min[0]=max[0]=x;
	min[1]=max[1]=y;
	min[2]=max[2]=z;
    
    
    for(int i=1;i<cant;i++){
        fscanf(file,"%u",&type);
        fscanf(file,"%s",word);
        x=atof(word);
        fscanf(file,"%s",word);
        y=atof(word);
        fscanf(file,"%s",word);
        z=atof(word);
        Point3D p (x,y,z);
        
        if(type == 3 || type == 2){
            fscanf(file,"%s",word);
            fscanf(file,"%s",word);
            fscanf(file,"%s",word);
        }
        points.push_back(p);
        
        if(min[0]>x)
			min[0]=x;
		else if(max[0]<x)
			max[0]=x;
		if(min[1]>y)
			min[1]=y;
		else if(max[1]<y)
			max[1]=y;
		if(min[2]>z)
			min[2]=z;
		else if(max[2]<z)
			max[2]=z;
    }
    
    while(true){
        if(!fscanf(file,"%s",word)){
            fclose(file);
            return false;
        }
        if(!strcmp(word,"[Elements,\0"))
            break;
    }
    
    fscanf(file,"%s",word);
    fscanf(file,"%u",&cant);
    elements.reserve(cant);
    
    unsigned int idx;
    
    for(int i=0;i<cant;i++){
        
        Element *element;
		fscanf(file,"%s",word);
		if (word[0]=='H') {
			vector<unsigned int> epts;
			epts.reserve(8);
			for(int j=0;j<8;j++){
				fscanf(file,"%u",&idx);
				epts.push_back(idx);
			}
			element = new Hexahedron(epts);
		}
		else if(word[0]=='T'){
            vector<unsigned int> epts;
			epts.reserve(4);
			for(int j=0;j<4;j++){
				fscanf(file,"%u",&idx);
				epts.push_back(idx);
			}
			element = new Tetrahedron(epts);
		}
		else if(word[0]=='P'){
            vector<unsigned int> epts;
			epts.reserve(5);
			for(int j=0;j<5;j++){
				fscanf(file,"%u",&idx);
				epts.push_back(idx);
			}
			element = new Pyramid(epts);
		}
		else if(word[0]=='R'){
            vector<unsigned int> epts;
			epts.reserve(6);
			for(int j=0;j<6;j++){
				fscanf(file,"%u",&idx);
				epts.push_back(idx);
			}
			element = new Prism(epts);
		}
		else{
			continue;
		}
        
		elements.push_back(element);
        
		for(int i=0;i<3;i++)
			fscanf(file,"%s",word);
    }
    
    bool e_colored = false;
    
    //see if there are colored elements
    while(true){
        if(fscanf(file,"%s",word)==EOF) {
            fclose(file);
            return true;
        }
        if(!strcmp(word,"colored\0")){
            e_colored = true;
            break;
        }
    }
    
    if (!e_colored) {
        fclose(file);
        return true;
    }
    
    while (fscanf(file,"%s",word)!=EOF) {
        if (!strcmp(word,"range\0")) {
            unsigned int start=0, end=0;
            float red=0, green=0, blue=0;
            fscanf(file,"%u",&start);
            fscanf(file,"%u",&end);
            fscanf(file,"%f",&red);
            fscanf(file,"%f",&green);
            fscanf(file,"%f",&blue);
            
            std::ostringstream ss;
            ss << " " << red << " " << green << " " << blue << " ";
            std::string color(ss.str());
            
            //cout << " coloring from " << start << " to " << end << " width " << color << "\n";
            
            for (unsigned int i=start; i<end; i++) {
                elements[i]->setColor(color);
            }
        }
    }
    
    fclose(file);
    return true;
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------
bool ViewVolume::readVtk(string name){
    
    
    char word [256];
    unsigned int cant,type;
    double x,y,z;
    
    FILE *file = fopen(name.c_str(),"r");
    
    while(true){
        if(!fscanf(file,"%s",word)){
            fclose(file);
            return false;
        }
        if(!strcmp(word,"UNSTRUCTURED_GRID\0"))
            break;
    }
    
    //read word Points
    fscanf(file,"%s",word);
    
    fscanf(file,"%u",&cant);
    points.reserve(cant);
    
    //read word float
    fscanf(file,"%s",word);
    
    if(cant<=0){
        cerr << "warning: no nodes were found\n";
        fclose(file);
        return true;
    }
    
    //read the first point (to initialize the bounds)
    fscanf(file,"%s",word);
    x=atof(word);
    fscanf(file,"%s",word);
    y=atof(word);
    fscanf(file,"%s",word);
    z=atof(word);
    Point3D p (x,y,z);
    
    points.push_back(p);
    
    min[0]=max[0]=x;
    min[1]=max[1]=y;
    min[2]=max[2]=z;
    
    for(int i=1;i<cant;i++){
        fscanf(file,"%s",word);
        x=atof(word);
        fscanf(file,"%s",word);
        y=atof(word);
        fscanf(file,"%s",word);
        z=atof(word);
        Point3D p (x,y,z);
    
        points.push_back(p);
        
        if(min[0]>x)
            min[0]=x;
        else if(max[0]<x)
            max[0]=x;
        if(min[1]>y)
            min[1]=y;
        else if(max[1]<y)
            max[1]=y;
        if(min[2]>z)
            min[2]=z;
        else if(max[2]<z)
            max[2]=z;
    }
    
    while(true){
        if(!fscanf(file,"%s",word)){
            fclose(file);
            return false;
        }
        if(!strcmp(word,"CELLS\0"))
            break;
    }
    
    fscanf(file,"%u",&cant);
    elements.reserve(cant);
    fscanf(file,"%s",word);
    
    unsigned int idx;
    
    for(int i=0;i<cant;i++){
        
        Element *element;
        fscanf(file,"%s",word);
        if (word[0]=='8') {
            vector<unsigned int> epts;
            epts.reserve(8);
            for(int j=0;j<8;j++){
                fscanf(file,"%u",&idx);
                epts.push_back(idx);
            }
            element = new Hexahedron(epts);
        }
        else if(word[0]=='4'){
            vector<unsigned int> epts;
            epts.reserve(4);
            for(int j=0;j<4;j++){
                fscanf(file,"%u",&idx);
                epts.push_back(idx);
            }
            element = new Tetrahedron(epts);
        }
        else if(word[0]=='5'){
            vector<unsigned int> epts;
            epts.reserve(5);
            for(int j=0;j<5;j++){
                fscanf(file,"%u",&idx);
                epts.push_back(idx);
            }
            element = new Pyramid(epts);
        }
        else if(word[0]=='6'){
            vector<unsigned int> epts;
            epts.reserve(6);
            for(int j=0;j<6;j++){
                fscanf(file,"%u",&idx);
                epts.push_back(idx);
            }
            element = new Prism(epts);
        }
        else{
            continue;
        }
        
        elements.push_back(element);
    }
    
    fclose(file);
    return true;
}


//--------------------------------------------------------------------
//--------------------------------------------------------------------
bool ViewVolume::readMixed(string name){
    
    
    char word [256];
    unsigned int nodeQ=0, elementQ=0;
    double x,y,z;
    
    FILE *file = fopen(name.c_str(),"r");
    
    while(true){
        if(!fscanf(file,"%s",word)){
            fclose(file);
            return false;
        }
        if(!strcmp(word,"MIXED\0"))
            break;
    }
    
    
    //read node and element quantity:
    fscanf(file,"%u",&nodeQ);
    points.reserve(nodeQ);
    
    fscanf(file,"%u",&elementQ);
    elements.reserve(elementQ);
    
    
    if(nodeQ==0){
        cerr << "warning: no nodes were found\n";
        fclose(file);
        return true;
    }
    
    //read the first point (to initialize the bounds)
    fscanf(file,"%s",word);
    x=atof(word);
    fscanf(file,"%s",word);
    y=atof(word);
    fscanf(file,"%s",word);
    z=atof(word);
    Point3D p (x,y,z);
    
    points.push_back(p);
    
    min[0]=max[0]=x;
    min[1]=max[1]=y;
    min[2]=max[2]=z;
    
    
    for (unsigned int i=1; i<nodeQ; i++) {
        fscanf(file,"%s",word);
        x=atof(word);
        fscanf(file,"%s",word);
        y=atof(word);
        fscanf(file,"%s",word);
        z=atof(word);
        Point3D p (x,y,z);
        points.push_back(p);
        
        if(min[0]>x)
            min[0]=x;
        else if(max[0]<x)
            max[0]=x;
        if(min[1]>y)
            min[1]=y;
        else if(max[1]<y)
            max[1]=y;
        if(min[2]>z)
            min[2]=z;
        else if(max[2]<z)
            max[2]=z;
    }
    
    elements.reserve(elementQ);
    
    unsigned int idx;
    int noden=0;
    
    for(unsigned int i=0;i<elementQ;i++){
        
        Element *element;
        vector<unsigned int> epts;
        fscanf(file,"%s",word);
        
        switch (word[0]) {
            case 'H':
                epts.reserve(8);
                for(int j=0;j<8;j++){
                    fscanf(file,"%u",&idx);
                    epts.push_back(idx);
                }
                element = new Hexahedron(epts);
                break;
            case 'R':
                epts.reserve(6);
                for(int j=0;j<6;j++){
                    fscanf(file,"%u",&idx);
                    epts.push_back(idx);
                }
                element = new Prism(epts);
                break;
            case 'P':
                noden = 0;
                if (word[1]=='\0') {
                    noden = 5;
                }
                else {
                    noden = word[1]-'0';
                }
                epts.reserve(noden);
                for(int j=0;j<noden;j++){
                    fscanf(file,"%u",&idx);
                    epts.push_back(idx);
                }
                switch (noden) {
                    case 4:
                        element = new Tetrahedron(epts);
                        break;
                    case 5:
                        element = new Pyramid(epts);
                        break;
                    case 6:
                        element = new PyramidBase5(epts);
                        break;
                    default:
                        break;
                }
                break;
            case 'T':
                epts.reserve(4);
                for(int j=0;j<4;j++){
                    fscanf(file,"%u",&idx);
                    epts.push_back(idx);
                }
                element = new Tetrahedron(epts);
                break;
            default:
                continue;
                break;
        }
        
        elements.push_back(element);
    }
    
    bool e_colored = false;
    
    //see if there are colored elements
    while(true){
        if(fscanf(file,"%s",word)==EOF) {
            fclose(file);
            return true;
        }
        if(!strcmp(word,"colored\0")){
            e_colored = true;
            break;
        }
    }
    
    if (!e_colored) {
        fclose(file);
        return true;
    }

    float red=0, green=0, blue=0;

    while (fscanf(file,"%s",word)!=EOF) {

        red = atof(word);
//        fscanf(file,"%f",&red);
        fscanf(file,"%f",&green);
        fscanf(file,"%f",&blue);
        
        fscanf(file,"%s",word);

        if (!strcmp(word,"range\0")) {
            int start=0, end=0;
            fscanf(file,"%i",&start);
            fscanf(file,"%i",&end);
            
            std::ostringstream ss;
            ss << " " << red << " " << green << " " << blue << " ";
            std::string color(ss.str());
            
            //cout << " coloring from " << start << " to " << end << " width " << color << "\n";
            
            for (unsigned int i=start; i<end; i++) {
                elements[i]->setColor(color);
            }
        }
        if (!strcmp(word,"list\0")) {
            int quantity = 0,index =0;
            fscanf(file,"%i",&quantity);
            list<int> celem;
            for (unsigned int i=0; i<quantity; i++) {
                fscanf(file,"%i",&index);
                celem.push_back(index);
            }
            
            std::ostringstream ss;
            ss << " " << red << " " << green << " " << blue << " ";
            std::string color(ss.str());
            
            //cout << " coloring from width " << color << " elements:";
            
            list<int>::iterator eit;
            for (eit=celem.begin(); eit!=celem.end(); eit++) {
                elements[*eit]->setColor(color);
                //cout << " " << *eit;
            }
            //cout << "\n";
        }

    }
    
    fclose(file);
    return true;
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------

bool ViewVolume::readData(char *name, const unsigned int &type){
	
    
    switch (type) {
        case 0:
            readM3d(name);
            break;
        case 1:
            readMixed(name);
            break;
        case 2:
            readVtk(name);
            break;
        default:
            return false;
    }
    
	if ((int)points.size()<3) {
		cout << "cannot display a mesh with less than 3 nodes\n";
		return false;
	}
	
	if ((int)elements.size()<1) {
		cout << "mesh has no elements: nothing to display\n";
		return false;
	}
	
	//save faces with neigbhoring information
	for(unsigned int i=0; i<elements.size(); i++){
		//An EnhancedElement is an Element with optimizations for
		//removing and visualizing.
		elements[i]->computeBbox(points);
		
		//get the element face index points and initialize it
		for(int j=0;j<elements[i]->numberOfFaces();j++){
            //modified loop
            vector<vector<unsigned int>> fptsVec = elements[i]->getFacePoints(j);
            for(int k = 0; k < fptsVec.size(); k++) {
                vector<unsigned int> fpts = fptsVec.at(k);
                Face face(fpts);
                int fid = all_faces.addFace(face);
                all_faces.getFace(fid).addElement(i);
                //add the face to the EnhancedElement
                elements[i]->addFace(fid);
            }

            //original for loop
			// vector<unsigned int> fpts = elements[i]->getFacePoints(j);
            // Face face(fpts);

			// int fid = all_faces.addFace(face);
			// all_faces.getFace(fid).addElement(i);
			// //add the face to the EnhancedElement
			// elements[i]->addFace(fid);
		}
	}
	
	return true;
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------

bool ViewVolume::execute(char *input, const bool &edge, const unsigned int &type){
	string command, file, question;
	int sys=0;
	bool cont = true;
	list<Point3D> tmppoints;
	list<Face> tmpfaces;
	
	if (!readData(input,type)) {
		return false;
	}
	
	//detect surface elements (elements with at least one face
	//not shared by two elements).
	//this method will change the color of surface faces.
	for(unsigned int i=0; i<elements.size(); i++){
		elements[i]->applyColor(all_faces);
	}
	
	question = "continue with visualization\0";
	
	while(cont){
		
		ofstream fgeom;
		//the instruction that will be executed.
		//-b is for setting the background color to white
		command = "geomview -b 1 1 1 \0";
		
		if(!newTmpFile(file,fgeom)){
			cout << "  Error: can't generate a tmp file for visualization\n";
			return false;
		}
        
        FaceContainer temp_fc;
		
		//askRemove returns true if the user removed some
		//faces. In that case, temporal face container (tfp)
		//will be employed, instead of all face container
		//(all_faces), to creat the list of faces to display.
		if (askRemove(temp_fc)) {
			if (!edge) {
				clean(tmppoints,tmpfaces,temp_fc);
			}
			else {
				cleanUnusedNodes(tmppoints,tmpfaces,temp_fc);
			}
		}
		else {
			if (!edge) {
                //en esta funcion cambia tmpfaces y tmppoints
				clean(tmppoints,tmpfaces,all_faces);
			}
			else {
				cleanUnusedNodes(tmppoints,tmpfaces,all_faces);
			}
		}
		
        //very useful
        //falta ver como dividir un cuadrado en 2
		headGeom(fgeom,tmppoints.size(),tmpfaces.size());
        //this ones gotta change
		dataGeom(fgeom,tmppoints,tmpfaces);
        
		fgeom.close();
		//execute geomview
		command+=file;
		//i think this is useless
		sys=system(command.c_str());
		
		if(sys!=0){
			return false;
		}
		
		command.clear();
		command = "rm \0";
		command += file;
		system(command.c_str());
		cont = yesOrNo(question);
		command.clear();
		file.clear();
		
		tmpfaces.clear();
		tmppoints.clear();
	}
	
	return true;
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------

bool ViewVolume::askRemove(FaceContainer &tfc){//list<Point> &tpoints, list<Face> &tfaces){
	int axe;
	bool removed = false;
	string question;
	double start=0,end=0;
	
	question = "remove some faces\0";
	removed = yesOrNo(question);
	if(removed){
		printBbox();
		axe = askAxe();
		cout << "indicate range for elements to be viewed\n";
		cout << "(starting) ";
		start = readDouble();
		cout << "(ending) ";
		end = readDouble();
		removeFaces(axe,start,end,tfc);//tpoints,tfaces);
	}
	return removed;
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------

void ViewVolume::printBbox(){
	cout << "the geometry Bounding box is:\n";
	cout << "X in [" << min[0] << "," << max[0] << "]\n";
	cout << "Y in [" << min[1] << "," << max[1] << "]\n";
	cout << "Z in [" << min[2] << "," << max[2] << "]\n";
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------

//save in tfaces only the faces that belong to one element from the given
//FaceContainer. Then it only save in tpoints the points used by the
//selected faces.
void ViewVolume::clean(list<Point3D> &tpoints, list<Face> &tfaces,
					FaceContainer &tfc){
	
    /*
    tpoints es tmppoints (inicialmente 0)
    tfaces es tmpfaces (inicialmente 0)
    tfc es all_faces (inicialmente idk pero creo que no sirve ese numero)
    */

	vector<int> ptsidx (points.size(),-1);
	int usedpts = 0;
	for(int i=0;i<tfc.length();i++){
		if(tfc.getFace(i).numberOfElements()!=1)
			continue;
		vector<unsigned int> fpts = tfc.getFace(i).getPoints();
		for(int j=0;j<fpts.size();j++){
			if(ptsidx[fpts[j]]==-1){
				tpoints.push_back(points[fpts[j]]);
				ptsidx[fpts[j]]=usedpts;
				fpts[j]=usedpts++;
			}
			else{
				fpts[j]=ptsidx[fpts[j]];
			}
		}

        //changes in clean()

        // if (fpts.size() == 4) {
        //     vector<unsigned int> fpts1 = fpts, fpts2 = fpts;
        //     fpts1.erase(fpts1.begin() + 3);
        //     fpts2.erase(fpts2.begin() + 1);
        //     // cout << fpts.at(0) << " " << fpts.at(1) << " " << fpts.at(2) << " " << fpts.at(3) << endl;
        //     // cout << fpts1.at(0) << " " << fpts1.at(1) << " " << fpts1.at(2) << endl;
        //     // cout << fpts2.at(0) << " " << fpts2.at(1) << " " << fpts2.at(2) << endl;
        //     Face f1(fpts1);
        //     f1.setColor(tfc.getFace(i).getColor());
        //     tfaces.push_back(f1);
        //     Face f2(fpts2);
        //     f2.setColor(tfc.getFace(i).getColor());
        //     tfaces.push_back(f2);

        // }
        // else {
        //     Face f(fpts);
		//     f.setColor(tfc.getFace(i).getColor());
		//     tfaces.push_back(f);
        // }

        //original clean()
        Face f(fpts);
        f.setColor(tfc.getFace(i).getColor());
        tfaces.push_back(f);
		
	}
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------

void ViewVolume::cleanUnusedNodes(list<Point3D> &tpoints,
							   list<Face> &tfaces,
							   FaceContainer &tfc){
	
	vector<int> ptsidx (points.size(),-1);
	int usedpts = 0;
	for(int i=0;i<tfc.length();i++){
		vector<unsigned int> fpts = tfc.getFace(i).getPoints();
		for(int j=0;j<fpts.size();j++){
			if(ptsidx[fpts[j]]==-1){
				tpoints.push_back(points[fpts[j]]);
				ptsidx[fpts[j]]=usedpts;
				fpts[j]=usedpts++;
			}
			else{
				fpts[j]=ptsidx[fpts[j]];
			}
		}
		Face f(fpts);
		f.setColor(tfc.getFace(i).getColor());
		tfaces.push_back(f);
	}
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------


void ViewVolume::removeFaces(int &axe, double start, double end,
                             FaceContainer &tfc){
    //list<Point> &tpoints, list<Face> &tfaces){
	
	cout << "removing elements outside range [";
	cout << start << "," << end <<"] in axis " << axe << " ... ";
	
	//detect elements outside the range to remove them
	for(int i=0;i<elements.size();i++){
		elements[i]->remove(axe,start,end);
	}
	cout << "done\n";
	
	//FaceContainer *tfc = new FaceContainer();
	for(int i=0;i<elements.size();i++){
		if(elements[i]->wasRemoved())
			continue;
		list<unsigned int>::iterator iter;
		list<unsigned int> fele = elements[i]->getFaces();
		for(iter=fele.begin();iter!=fele.end();iter++){
			vector<unsigned int> tfp = all_faces.getFace(*iter).getPoints();
			Face nface (tfp);
			nface.addElement(i);
			nface.setColor(all_faces.getFace(*iter).getColor());
			tfc.addFace(nface);
		}
	}
	//clean(tpoints,tfaces,tfc);
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------

bool ViewVolume::yesOrNo(string &question){
	char cont;
	cout << question << "[y/n]? ";
	cin >> cont;
	while(cont!='y' && cont!='n'){
		cout << question << "[y/n]? ";
		cin >> cont;
	}
	if(cont=='y')
		return true;
	return false;
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------

int ViewVolume::askAxe(){
	char axe;
	cout << "indicate axe[x/y/z]: ";
	cin >> axe;
	while(true){
		if(axe=='x')
			return 0;
		else if(axe=='y')
			return 1;
		else if(axe=='z')
			return 2;
		cout << "invalid. indicate axe[x/y/z]: ";
		cin >> axe;
	}
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------

double ViewVolume::readDouble(){
	double coord;
	cout << "indicate coord: ";
	cin >> coord;
	return coord;
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------
//escribe en el archivo la cantidad de puntos y caras de la malla
void ViewVolume::headGeom(ofstream &output, int np, int nf){
	
	output << "OFF\n";
    output << np << " " << nf << " 0\n\n";
}

//escribe en el archivo los puntos de la malla
//en el segundo for escribe en el archivo las caras de la malla (triangulos y cuadrados)
void ViewVolume::dataGeom(ofstream &output, list<Point3D> &tpoints,
					   list<Face> &tfaces){
	
	list<Point3D>::iterator iterp;
	list<Face>::iterator iterf;
	//print each point
    //este no me importa, esta bien como esta
	for(iterp=tpoints.begin();iterp!=tpoints.end();iterp++){
		Point3D p = *iterp;
        //cout << p[0] << " " << p[1] << " " << p[2] << endl;
		output << p[0] << " ";
		output << p[1] << " ";
		output << p[2] << "\n";
	}
	
	output << endl;
	
	//print each face
	for(iterf=tfaces.begin();iterf!=tfaces.end();iterf++){
        //cout << *iterf << " " << (*iterf).getColor() << endl;
        //this ones gotta change
		output << *iterf << " ";
        //this ones fine
		output << (*iterf).getColor() << endl;
	}
}

bool ViewVolume::newTmpFile(string &stmp,ofstream &output){
    char tmp [16];
    int fd;
	
    strcpy(tmp,"geomviewXXXXXX\0");
    
    fd=mkstemp(tmp);
    
    if(fd==-1){
		cout << "  Error: Can't create a tmp file. Try saving";
		cout << " before viewing.\n";
		return false;
    }
    
    stmp+=tmp;
    output.open(tmp);
    return true;
}
