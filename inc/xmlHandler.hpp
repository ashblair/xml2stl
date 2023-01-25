#pragma once

class stringops
{
public:
    stringops(void){}

    static std::string fullTrim(const std::string &);
    static std::string getNameNoExt(const std::string &);
    static std::string backOne(const std::string &);

};

typedef class xmlnode * pxmlnode;
typedef class xmlelement * pxmlelement;
class xmlelement
{
public:
    xmlelement(void):eType(' '), elm(NULL) {}
    xmlelement(const xmlelement &);
    xmlelement(const char t, void * p): eType(t), elm(p) {}
    ~xmlelement(void);

    std::string write_element(void);

    char eType; // 'n'=node, 'c'=comment, 'p'=prolog, 't'=text
    void * elm; 
};

class xmlnode
{
public:
    xmlnode(void) {}
    xmlnode(const xmlnode & x): nodeID(x.nodeID), attributes(x.attributes), elements(x.elements) {}
    ~xmlnode(void); 

    std::string write_node(void);
    bool findAttribute(std::string &);
    //bool is_section(string &, string &, string &);
    //bool has_child(void);
    //bool is_widget(string &, string &);
    //void match_name_to_id(const string &);

    std::string nodeID; // this is right after the tag-start <
    std::vector<std::string> attributes;
    std::vector<pxmlelement> elements;

};

class xmlHandler
{
public:
    xmlHandler(void){}
    ~xmlHandler(void); 

    void readXMLfile(std::string);
    size_t addElement(const size_t, std::vector<pxmlelement> & );
    std::vector<pxmlnode> getAllNodes(std::vector<pxmlelement> &);
    size_t addNode(const size_t, std::vector<pxmlelement> &);
    size_t getAttributes(const size_t, std::vector<std::string> &);
    void LoadXML(const std::string &);
    void getRoot(void);
    void getSections(void);
    void buildSections(void);
    //void handleError(const std::string &);

    static bool getTextFromNode(const pxmlnode &, std::string &);
    static pxmlnode getNode(const std::string &, std::vector<pxmlelement> &);
    static std::vector<pxmlnode> getAllNodes(const std::string &, std::vector<pxmlelement> &);
    static pxmlnode get_next_node(std::vector<pxmlelement>::iterator &, const std::vector<pxmlelement> &);
    static pxmlnode get_next_node(const std::string &, std::vector<pxmlelement>::iterator &, const std::vector<pxmlelement> &);

    void getTransform(void);
    void LoadAndParseXML(const std::string &);

    
    static std::string xmlIN;
    static int nodeCount, elementCount;
    stringops sOps;
    pxmlnode root;
    pxmlelement lstElm;
    std::vector<pxmlelement> XML;
    std::vector<pxmlnode> sections;
};
