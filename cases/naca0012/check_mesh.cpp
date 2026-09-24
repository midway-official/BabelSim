// Independent audit of the written file using BabelSim's actual mesh reader.
#include "babelsim/mesh_io.h"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
using namespace babelsim;
int main(int argc,char** argv) {
    if(argc!=2) return 2;
    const Mesh mesh=readMeshFile(argv[1]);
    mesh.validate();
    double vmin=1e99,vsum=0,maxangle=0,maxskew=0,maxclosure=0,dmin=1e99,dmax=0;
    int hex=0,prism=0;
    for(Index c=0;c<mesh.cellCount();++c){
        vmin=std::min(vmin,mesh.cellVolume(c));vsum+=mesh.cellVolume(c);
        Vec3 closure{}; double area=0;
        for(auto f:mesh.cellFaces(c)) {closure+=(mesh.owner(f)==c?1.0:-1.0)*mesh.faceAreaVector(f);area+=mesh.faceArea(f);}
        maxclosure=std::max(maxclosure,norm(closure)/area);
        hex+=mesh.cellVertices(c).size()==8; prism+=mesh.cellVertices(c).size()==6;
    }
    for(Index f=0;f<mesh.faceCount();++f){
        auto owner=mesh.owner(f); auto next=mesh.neighbour(f);
        if(mesh.boundaryFace(f)) {
            if(mesh.patchName(mesh.boundaryPatch(f))=="airfoil"){
                double d=dot(mesh.faceCentre(f)-mesh.cellCentre(owner),mesh.faceNormal(f));
                dmin=std::min(dmin,d); dmax=std::max(dmax,d);
            }
            continue;
        }
        auto delta=mesh.cellCentre(next)-mesh.cellCentre(owner);
        auto normal=mesh.faceNormal(f);
        maxangle=std::max(maxangle,std::acos(std::clamp(dot(delta,normal)/norm(delta),-1.,1.))*180/std::acos(-1.));
        double fraction=dot(mesh.faceCentre(f)-mesh.cellCentre(owner),normal)/dot(delta,normal);
        maxskew=std::max(maxskew,norm(mesh.faceCentre(f)-mesh.cellCentre(owner)-fraction*delta)/norm(delta));
    }
    std::cout<<std::setprecision(16)<<"{\n  \"cells\": "<<mesh.cellCount()<<",\n  \"hexahedra\": "<<hex<<",\n  \"prisms\": "<<prism
      <<",\n  \"min_volume\": "<<vmin<<",\n  \"total_volume\": "<<vsum<<",\n  \"max_nonorthogonality_deg\": "<<maxangle
      <<",\n  \"max_normalized_skewness\": "<<maxskew<<",\n  \"max_relative_cell_area_closure\": "<<maxclosure
      <<",\n  \"wall_cell_distance_min\": "<<dmin<<",\n  \"wall_cell_distance_max\": "<<dmax<<"\n}\n";
    return !(vmin>0 && maxangle<70 && maxskew<.5 && maxclosure<1e-10 && hex+prism==mesh.cellCount());
}
