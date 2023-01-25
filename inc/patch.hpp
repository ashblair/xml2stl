#pragma once
#include "section.hpp"

class xmlnode;
typedef class xmlnode * pxmlnode;

typedef class patch * ppatch;

class patch: public section
{
public:
    patch(void);
    patch(const patch &);
    ~patch(void);

    char what(void) {return 'x';}
    void getSection(const pxmlnode &);
    std::vector<unsigned int> getPart(const std::string &);
    void buildSection(const pxmlnode &);
    std::string printSection(void);
    void add_arc(const psection &);
    void sweep2line(const psection &);
    void make_triangles(const std::vector<std::vector<unsigned int>> &, const std::vector<std::vector<double>> &);




    std::string id, cylinderID;
    char cylinderPART; // a = start_angle alpha, b = end_angle beta, f = first_circle, l = last_circle
    static std::map<char, std::string> fullPART;
    static std::map<std::string, char> idxPART;
    std::vector<unsigned int> fromVerts;
    std::vector<double> oFacts;  // could be angles for circular part or t for line part
    struct toType
    {
        toType(void):sFrom(""), sPart("") {}
        toType(const toType & t):sFrom(t.sFrom), sPart(t.sPart), verts(t.verts), orderFactors(t.orderFactors){}
        toType(const std::string & f, const std::string & p, const std::vector<unsigned int> & v):sFrom(f), sPart(p), verts(v){}
        std::string sFrom, sPart;
        std::vector<unsigned int> verts;
        std::vector<double> orderFactors; // 1 to 1 array w/ verts        
    };
    std::vector<toType> sectsTo;
    unsigned int start_triangles, total_triangles;
    int orientation; //1 for counterclockwise, -1 for clockwise
};