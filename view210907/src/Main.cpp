#include "Main.h"

using namespace std;

void emsg(){
	cout << "use: ./viewer -o file.m3d\n";
	cout << " where -o can be:\n";
	cout << "  -d: for mdl\n";
    cout << "  -o: for off\n";
	cout << "  -m: for m3d\n";
    cout << "  -v: for mvm\n"; //mixed volume mesh
    cout << "  -t: for vtk\n";
    cout << "  -e: for enabling inner edges display (volumen meshes)\n";
}

int main(int argc,char** argv){
	string op;
	
	if(argc!=3){
		emsg();
		return 1;
	}
	
    ViewVolume vv;
    ViewSurface vs;
    
    switch (argv[1][1]) {
        case 'm':
            if (strlen(argv[1])==3 && argv[1][2]=='e') {
                vv.execute(argv[2],true,0);
            }
            else{
                vv.execute(argv[2],false,0);
            }
            break;
        case 'v':
            if (strlen(argv[1])==3 && argv[1][2]=='e') {
                vv.execute(argv[2],true,1);
            }
            else{
                vv.execute(argv[2],false,1);
            }
            break;
        case 't':
            if (strlen(argv[1])==3 && argv[1][2]=='e') {
                vv.execute(argv[2],true,2);
            }
            else{
                vv.execute(argv[2],false,2);
            }
            break;
        case 'd':
            vs.execute(argv[2]);
            break;
        default:
            emsg();
            return 1;
    }
	
	return 0;
}

