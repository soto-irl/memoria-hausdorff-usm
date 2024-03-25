#include "ViewSurface.h"

ViewSurface::ViewSurface(){
    op='\0';
}

ViewSurface::~ViewSurface(){

}

bool ViewSurface::readMdl(string name) {
    
    char word [256];
    int cant=0;
    double x=0,y=0,z=0;
    unsigned int line = 0;
    
    FILE *file = fopen(name.c_str(),"r");
    
    if (file==NULL) {
        std::cout << "File " << name << " doesn't exist\n";
        return false;
    }
    
    //read number of nodes
    while(true){
        if(std::fscanf(file,"%s",word)==EOF){
            fclose(file);
            return false;
        }
        if(!strcmp(word,"ARRAY1<POINT3D>]\0"))
            break;
    }
    
    std::fscanf(file,"%i",&cant);
    
    if(cant<=0){
        fclose(file);
        return false;
    }

    points.reserve(cant);
    
    //scan the first node to have something to compare with
    
    std::fscanf(file,"%s",word);
    x=atof(word);
    std::fscanf(file,"%s",word);
    y=atof(word);
    std::fscanf(file,"%s",word);
    z=atof(word);
    Point3D p (x,y,z);
    points.push_back(p);
    
    
    min[0] = max[0] = x;
    min[1] = max[1] = y;
    min[2] = max[2] = z;
    
    //read each node
    for( int i=1;i<cant;i++){
        std::fscanf(file,"%s",word);
        x=atof(word);
        std::fscanf(file,"%s",word);
        y=atof(word);
        std::fscanf(file,"%s",word);
        z=atof(word);
        Point3D p (x,y,z);
        points.push_back(p);
        
        //compute Bbox
        if(min[0]>x)
            min[0]=x;
        if(max[0]<x)
            max[0]=x;
        if(min[1]>y)
            min[1]=y;
        if(max[1]<y)
            max[1]=y;
        if(min[2]>z)
            min[2]=z;
        if(max[2]<z)
            max[2]=z;
        
    }
    
    //read number of faces
    cant = 0;
    while(1){
        if(std::fscanf(file,"%s",word) == EOF){
            std::cout << "didn't find faces\n";
            fclose(file);
            return false;
        }
        
        if(!strcmp(word,"ARRAY1<STRING>]\0")){
            //std::fscanf(file,"%s",word);
            std::fscanf(file,"%i",&cant);
            break;
        }
    }
    
    if(cant<=0){
        fclose(file);
        return false;
    }
    
    //this vector will hold the information of removed
    //faces. This is to avoid creating an enhancedFace object
    //just to know this info.
    removedfaces.assign(cant,false);
    
    faces.reserve(cant);
    int dust;
    
    //read each face
    for( int i=0;i<cant;i++){
        vector<unsigned int> face (3,0);
        std::fscanf(file,"%u",&face[0]);
        std::fscanf(file,"%u",&face[1]);
        std::fscanf(file,"%u",&face[2]);
        
        faces.push_back(face);
        
        //read some unnecessary integers
        for(unsigned int j=0;j<3;j++)
            std::fscanf(file,"%s",word);
    }
    fclose(file);
    
    return true;
    
}


/*bool ViewSurface::readMdlData(char *name){
    //ReadSurface *rm;
    
    //read point and element data
    //rm = new ReadSurface(name);
    //rm->execute(points,fc);
    //this vector will hold the information of removed
    //faces. This is to avoid creating an enhancedFace object
    //just to know this info.
    //removedfaces.assign(fc->length(),false);
    
    if ((int)points.size()<3) {
        cout << "cannot display a mesh with less than 3 nodes\n";
        return false;
    }
    
    if (fc->length()<1) {
        cout << "mesh has no faces: nothing to display\n";
        return false;
    }
    
    min[0] = max[0] = points[0][0];
    min[1] = max[1] = points[0][1];
    min[2] = max[2] = points[0][2];
    
    for(int i=1;i<points.size();i++){
        
        double x = points[i][0];
        double y = points[i][1];
        double z = points[i][2];
        
        //compute Bbox
        if(min[0]>x)
            min[0]=x;
        if(max[0]<x)
            max[0]=x;
        if(min[1]>y)
            min[1]=y;
        if(max[1]<y)
            max[1]=y;
        if(min[2]>z)
            min[2]=z;
        if(max[2]<z)
            max[2]=z;
    }
    
    return true;
}*/


bool ViewSurface::execute(char *input){
    string command,file,question;
    int sys=0;
    bool cont = true;
    list<Point3D> tmppoints;
    list<vector<unsigned int> > tmpfaces;
    //FaceContainer *tmpf;
    
    if (!readMdl(input)) {
        return false;
    }
    
    question = "continue with visualization\0";
    
    while(cont){
        
        ofstream fgeom;
        //the instruction that will be executed.
        command = "geomview -b 1 1 1 \0";
        
        if(!newTmpFile(file,fgeom)){
            cerr << "  Error: can't generate a tmp file for visualization\n";
            return false;
        }
        
        if (askRemove(tmppoints,tmpfaces)) {
            headGeom(fgeom,tmppoints.size(),tmpfaces.size());
            dataGeom(fgeom,tmppoints,tmpfaces);
        }
        else {
            headGeom(fgeom,points.size(),faces.size());
            dataGeom(fgeom);
        }
        
        
        fgeom.close();
        
        //execute geomview
        command+=file;
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
        removedfaces.assign(faces.size(),false);
    }
    
    return true;
}

bool ViewSurface::askRemove(list<Point3D> &tpoints, list<vector<unsigned int> > &tfaces){
    int axe;
    bool removed = false;
    string question;
    double start=0,end=0;
    
    question = "remove some faces\0";
    removed = yesOrNo(question);
    if(removed){
        printBbox();
        axe = askAxe();
        cout << "indicate range for faces to be viewed\n";
        cout << "(starting) ";
        start = readDouble();
        cout << "(ending) ";
        end = readDouble();
        removeFaces(axe,start,end,tpoints,tfaces);
    }

    return removed;
}

//save in tfaces only the faces that belong to one element from the given
//FaceContainer. Then it only save in tpoints the points used by the
//selected faces.
void ViewSurface::clean(list<Point3D> &tpoints, list<vector<unsigned int> > &tfaces){
    
    vector<unsigned int> ptsidx (points.size(),0);
    vector<bool> used (points.size(),false);
    int usedpts = 0;
    list<vector<unsigned int> >::iterator isf;
    
    for (unsigned int i=0; i<removedfaces.size(); i++) {
        if (!removedfaces[i]) {
            vector<unsigned int> fpts = faces[i];
            vector<unsigned int> nfpts (fpts.size(),0);
            for (unsigned int j=0; j<fpts.size(); j++){
                if(!used[fpts[j]]){
                    tpoints.push_back(points[fpts[j]]);
                    ptsidx[fpts[j]]=usedpts;
                    used[fpts[j]] = true;
                    nfpts[j]=usedpts++;
                }
                else{
                    nfpts[j]=ptsidx[fpts[j]];
                }
            }
            tfaces.push_back(nfpts);
        }
    }
}

void ViewSurface::removeFaces(int &axe, double start, double end,
                              list<Point3D> &tpoints, list<vector<unsigned int> > &tfaces){
    
    cout << "removing faces outside range [";
    cout << start << "," << end <<"] in axis " << axe << " ... ";
    
    bool removed = false;
    
    //detect faces outside the range to remove them
    for (unsigned int i=0; i<faces.size(); i++){
        removed = false;
        vector<unsigned int> fpts = faces[i];
        for (unsigned int j=0; j<fpts.size(); j++) {
            if (points[fpts[j]][axe]<start || points[fpts[j]][axe]>end) {
                removedfaces[i] = true;
                break;
            }
        }
    }
    
    clean(tpoints, tfaces);
    cout << "done\n";
}


void ViewSurface::printBbox(){
    cout << "the geometry Bounding box is:\n";
    cout << "X in [" << min[0] << "," << max[0] << "]\n";
    cout << "Y in [" << min[1] << "," << max[1] << "]\n";
    cout << "Z in [" << min[2] << "," << max[2] << "]\n";
}

bool ViewSurface::yesOrNo(string &question){
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

int ViewSurface::askAxe(){
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

double ViewSurface::readDouble(){
    double coord;
    cout << "indicate coord: ";
    cin >> coord;
    return coord;
}

void ViewSurface::headGeom(ofstream &output, int np, int nf){
    output << "OFF\n";
    output << np << " " << nf << " 0\n\n";
}

void ViewSurface::dataGeom(ofstream &output, list<Point3D> &tpoints,
                           list<vector<unsigned int> > &tfaces){
    
    list<Point3D>::iterator iterp;
    list<vector<unsigned int> >::iterator iterf;
    //print each point
    for(iterp=tpoints.begin();iterp!=tpoints.end();iterp++){
        Point3D p = *iterp;
        output << p[0] << " ";
        output << p[1] << " ";
        output << p[2] << "\n";
    }
    
    output << endl;
    
    //print each face
    for(iterf=tfaces.begin();iterf!=tfaces.end();iterf++){
        vector<unsigned int> face = *iterf;
        output << face.size();
        for (unsigned int i=0; i<face.size(); i++) {
            output << " " << face[i];
        }
        output << "\n";
    }
}

void ViewSurface::dataGeom(ofstream &output){
    
    //print each point
    for (unsigned int i=0; i<points.size(); i++) {
        output << points[i][0] << " ";
        output << points[i][1] << " ";
        output << points[i][2] << "\n";
    }
    
    output << endl;
    
    //print each face
    for (unsigned int i=0; i<faces.size(); i++) {
        output << faces[i].size();
        for (unsigned int j=0; j<faces[i].size(); j++) {
            output << " " << faces[i][j];
        }
        output << endl;
    }
}

bool ViewSurface::newTmpFile(string &stmp,ofstream &output){
    char tmp [16];
    int fd;
    
    strcpy(tmp,"geomviewXXXXXX\0");
    
    fd=mkstemp(tmp);
    
    if(fd==-1){
        cout << "  Error: Can't create a tmp file. Try saving before";
        cout << " viewing.\n";
        return false;
    }
    
    stmp+=tmp;
    output.open(tmp);
    return true;
}
