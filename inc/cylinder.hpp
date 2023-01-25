#pragma once
#include "section.hpp"
//#include "CompGeo.hpp"

class xmlnode;
typedef class xmlnode * pxmlnode;


typedef class cylinder * pcylinder;

class cylinder: public section
{
public:
    cylinder(void);
    cylinder(const cylinder &);
    ~cylinder(void);

    //void handleError(const std::string &);
    //bool getTextFromNode(const pxmlnode &, std::string &);
    //pxmlnode getNode(const std::string &, std::vector<pxmlelement> &);
    struct seamType
    {
        std::string from, part;
    };
    void getSection(const pxmlnode &);
    std::vector<unsigned int> getPart(const std::string &);
    void calculateOtherParameters(void);
    bool handleSeam(const int &, unsigned int &, const std::vector<unsigned int> &, const CompGeo::XYZ &, const std::string &);
    void calculateArcs(void); //(std::map<char, seamType>);
    void calculateTriangles(void);
    void buildSection(const pxmlnode &);
    //std::string printXYZ(const std::string &, const CompGeo::XYZ &);
    std::string printSection(void);
    //char get_type(void) {return 'c';}
    //std::string get_id(void) {return id;}
    double sqr(const double & a) {return a * a;}
    char what(void) {return 'c';}

    std::string id;
    double radius, error, alpha, beta, axle_length, axle_delta;
    int n2r, N, total_arcs, total_triangles, start_triangles;
    CompGeo::XYZ ctr_first, ctr_last, axleHAT, NHAT, IHAT, JHAT;
    std::vector<unsigned int> arcs;
    std::map<char, seamType> seams;
};