#pragma once
#include "section.hpp"
//#include "CompGeo.hpp"

class xmlnode;
typedef class xmlnode * pxmlnode;

typedef class plane * pplane;

class plane: public section
{
public:
    plane(void);
    plane(const plane &);
    ~plane(void);

    //void handleError(const std::string &);
    //bool getTextFromNode(const pxmlnode &, std::string &);
    //pxmlnode getNode(const std::string &, std::vector<pxmlelement> &);
    //std::vector<pxmlnode> getAllNodes(const std::string &, std::vector<pxmlelement> &);
    //pxmlnode get_next_node(std::vector<pxmlelement>::iterator &, const std::vector<pxmlelement> &);
    //pxmlnode get_next_node(const std::string &, std::vector<pxmlelement>::iterator &, const std::vector<pxmlelement> &);
    //CompGeo::XYZ getNamedVector(const std::string &);
    void getSection(const pxmlnode &);
    std::vector<unsigned int> getPart(const std::string &);
    void triangularize(void);
    void buildSection(const pxmlnode &);
    //std::string printXYZ(const std::string &, const CompGeo::XYZ &);
    std::string printSection(void);
    char what(void) {return 'p';}

    std::string id;
    CompGeo::XYZ norm;
    double multiplyer; // this times norm gives origin on plane
    unsigned int triStart, triCount;
    std::vector<unsigned int> verts, offsets, shapes;

};