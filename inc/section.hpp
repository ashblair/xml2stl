#pragma once

class xmlnode;
typedef class xmlnode * pxmlnode;

typedef class section * psection;
class section
{
public:
    virtual ~section(void) {}
    virtual std::vector<unsigned int> getPart(const std::string &) = 0;
    virtual void buildSection(const pxmlnode &) = 0;
    virtual std::string printSection(void) = 0;
    virtual char what(void) = 0;
    //char get_type(void) {return 'c';}
    //std::string get_id(void) {return id;}
    double sqr(const double & a) {return a * a;}
    static void handleError(const std::string &);
    CompGeo::XYZ getNamedVector(const std::string &);
    
};