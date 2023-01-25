#pragma once
#include "section.hpp"
//#include "CompGeo.hpp"
/*
typedef struct sectionStructType
{
    sectionStructType(void):sType(' '), sect(NULL){}
    sectionStructType(const sectionStructType &sIn):sType(sIn.sType), sect(sIn.sect){}
    sectionStructType(const char & st, void * v):sType(st), sect(v){}
    ~sectionStructType(void);

    char sType; //'c'=cylinder, 'p'=plane
    void * sect;

} sectionType, * psectionType;
*/


class Model
{
public:
    Model(void){}
    ~Model(void){v.clear(); t.clear(); 
        for (std::map<std::string, psection>::iterator mit = s.begin(); mit != s.end(); ++mit) {delete mit->second;}
        s.clear();}

    void writeDBG(const std::string &);
    void writePLY(const std::string &);
    void writeSTL(const std::string &);
    static void rotateT(const char &, const double &);
    static void translateT(CompGeo::XYZ &);
    void transformVerts(void);
    void finalizeModel(void);
    
    static std::vector<CompGeo::XYZ> v; // vertices
    static std::vector<unsigned int> t; // triangles as sets of 3 indices into v
    static std::map<std::string, psection> s; // key = id, value = ptr to section
    static Matrix T;
};