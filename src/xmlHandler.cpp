#include "first.hpp"
#include "xmlHandler.hpp"
#include "section.hpp"
#include "Model.hpp"
#include "cylinder.hpp"
#include "plane.hpp"
#include "patch.hpp"

using namespace std;

//static initialization:
string xmlHandler::xmlIN = "";
int xmlHandler::nodeCount = 0, xmlHandler::elementCount = 0;

string stringops::fullTrim(const string & s_in)
{
    string r = s_in;
    int i = 0; 
    while (i < r.length())
    {
        char a = r.at(i);
        if ((a == ' ') || (a == '\n') || (a == '\r') || (a == '\t'))
        {
            r = r.substr(0, i) + r.substr(i + 1);
        }
        else ++i;
    }
    return r;

}

string stringops::getNameNoExt(const string & full)
{
    if (full.length() == 0) return "";
    string rStr = full;
    if (rStr.at(0) == '~') rStr = rStr.substr(1);
    if (rStr.length() == 0) return "";
    size_t period = rStr.find_last_of('.'),
        slash = rStr.find_last_of('/', period);
    if ((period == string::npos) && (slash == string::npos)) return rStr;
    if (period != string::npos) rStr = rStr.substr(0, period);
    if (slash != string::npos) rStr = rStr.substr(slash + 1);
    return rStr;
}

string stringops::backOne(const string & full)
{
    if (full.length() == 0) return "";
    string rStr = full;
    if (rStr.at(0) == '~') rStr = "";
    if (rStr.length() == 0) return "";
    size_t slash = rStr.find_last_of('/');
    if ((0 == slash) || (string::npos == slash)) rStr = "";
    else
    {
        rStr = rStr.substr(0, slash);
    }
    return rStr;

}


xmlelement::xmlelement(const xmlelement & x)
{
    eType = x.eType;

    switch(eType)
    {
    case 'n':
        elm = new xmlnode(*static_cast<pxmlnode>(x.elm));
        break;
    case 'c':
    case 'p':
    case 't':
        elm = new string(*static_cast<string *>(x.elm));
        break; 
    }
}

xmlelement::~xmlelement(void)
{
    //pxmlnode p_n = NULL;
    //string * p_s = NULL;

    //cout << xmlHandler::elementCount << ") destructing a " << eType << " element.\n";
    ++xmlHandler::elementCount;

    switch(eType)
    {
    case 'n':
        delete static_cast<pxmlnode>(elm);
        break;
    case 'c':
    case 'p':
    case 't':
        delete static_cast<string *>(elm);
        break;
    }

}

string xmlelement::write_element(void)
{
    pxmlnode p_n = NULL;
    string * p_s = NULL;

    switch(eType)
    {
    case 'n':
        p_n = static_cast<pxmlnode>(elm);
        return p_n->write_node();
        break;
    case 'c':
        p_s = static_cast<string *>(elm);
        return "<!-- " + *p_s + "-->";
        break;
    case 'p':
        p_s = static_cast<string *>(elm);
        return "<?xml " + *p_s + "?>";
        break;
    case 't':
        p_s = static_cast<string *>(elm);
        return *p_s;
        break;
    }
    return "Duh!";
}

xmlnode::~xmlnode(void)
{
    //cout << xmlHandler::nodeCount << ") destructing node " << nodeID << "\n";
    ++xmlHandler::nodeCount;

    attributes.clear(); 
    for (vector<pxmlelement>::iterator eit = elements.begin(); eit != elements.end(); ++eit)
    {
        pxmlelement & e = *eit;
        delete e;
        /*
        switch(e->eType)
        {
        case 'n':
            delete static_cast<pxmlnode>(e->elm); // recursion
            break;
        case 'c':
        case 'p':
        case 't':
            delete static_cast<string *>(e->elm);
            break;
        }
        */
    }

    elements.clear();

}

string xmlnode::write_node(void)
{
    string buff = "<" + nodeID;
    for (size_t i = 0; i < attributes.size(); ++i)
    {
        if ((i % 2) == 0) buff += " " + attributes[i];
        else buff += "=\"" + attributes[i] + "\"";
    }
    size_t e = elements.size();
    if (e == 0)
    {
        buff += "/>";
        return buff;
    }
    buff += ">";
    for (size_t i = 0; i < e; ++i)
    {
        buff += elements[i]->write_element();
    }
    buff += "</" + nodeID + ">";
    return buff;

}

bool xmlnode::findAttribute(string & kv)
{ // returns true if found, kv set to value

    string key = kv;
    kv = "";
    for (vector<string>::iterator ait = attributes.begin(); ait != attributes.end(); ++ait)
    {
        string aKey = *ait, vKey = *(++ait);
        if (aKey.compare(key) == 0) 
        {
            kv = vKey;
            return true;
        }
    }
    return false;
}


/*
bool xmlnode::is_section(string &id, string &typ, string & err)
{
    typ = "";
    id = "";
    err = "";
    static unsigned int sectNum = 0;

    if (nodeID.compare("section") != 0) return false;

    bool idFound = false, typeFound = false;
    map<string, string> a;
    string akey = "", aval = "";
    vector<string>::iterator ait = attributes.begin();
    while (ait != attributes.end())
    {
        akey = *ait;
        ++ait;
        aval = *ait;
        a[akey] = aval;
        ++ait;
        if (!idFound) idFound = akey.compare("id") == 0;
        if (!typeFound) typeFound = akey.compare("type") == 0;
    }
    if (!typeFound)
    {
        err = "xml error: no type in section node!\n";
        return true;
    }
    typ = a["type"];
    id = idFound? a["id"]: string("__sect_") + to_string(++sectNum);

    return true;

}
*/

xmlHandler::~xmlHandler(void)
{
    nodeCount = 0; elementCount = 0;

    for (vector<pxmlelement>::iterator xit = XML.begin(); xit != XML.end(); ++xit)
    {
        pxmlelement & e = *xit;
        delete e;
    }

    XML.clear(); 
    
    sections.clear();

    //cout << "xmlHandler destructor called. nodes deleted:" << nodeCount << " elements deleted: " << elementCount << "\n";
}

void xmlHandler::readXMLfile(string gld)
{
    ifstream is (gld, ifstream::binary);
    assert(is);
    // get length of file:
    is.seekg (0, is.end);
    int length = is.tellg();
    is.seekg (0, is.beg);

    char * buffer = new char [length];

    //cout << "Reading " << to_string(length) << " characters... \n";
    // read data as a block:
    is.read (buffer,length);
    assert(is);

    xmlIN = string(buffer, length);

    delete[] buffer;    

    is.close();

}

size_t xmlHandler::addElement(const size_t pos, vector<pxmlelement> & elmts)
{ // pos & return are positions in static class member xmlIN (a string)
  // whenever pos does not point to a tag-start (<) a new xmlelement will be formed and added to elmts with the 
  // string of all characters from pos to just before the next tag-start (<)
    if (xmlIN[pos] == '<') return addNode(pos, elmts);
    size_t nxtTag = xmlIN.find('<', pos); assert (nxtTag != string::npos);
    string txt = xmlIN.substr(pos, nxtTag - pos);
    string *pTxt = new string;
    *pTxt = txt;
    pxmlelement x = new xmlelement('t', pTxt);
    elmts.push_back(x);
    lstElm = x;
    return nxtTag;
}

size_t xmlHandler::addNode(const size_t pos, vector<pxmlelement> & elmts)
{
    assert (xmlIN[pos] == '<');
    assert (xmlIN[pos + 1] != '/');
    // get nodeID
    size_t current = pos + 1, nxtPos = xmlIN.find_first_of(" />", current); assert (nxtPos != string::npos);

    string n_id = xmlIN.substr(current, nxtPos - current);

    if (n_id.compare("?xml") == 0)
    { // xml prolog s/b at most 1 of these right at the beginning
        current = nxtPos + 1; nxtPos = xmlIN.find("?>", current); assert (nxtPos != string::npos);
        string * pTxt = new string;
        *pTxt = xmlIN.substr(current, nxtPos - current);
        pxmlelement x = new xmlelement('p', pTxt);
        elmts.push_back(x);
        lstElm = x;
        return nxtPos + 2;
    }

    if (n_id.compare("!--") == 0)
    { // comment
        current = nxtPos + 1; nxtPos = xmlIN.find("-->", current); assert (nxtPos != string::npos);
        string * pTxt = new string;
        *pTxt = xmlIN.substr(current, nxtPos - current);
        pxmlelement x = new xmlelement('c', pTxt);
        elmts.push_back(x);
        lstElm = x;
        return nxtPos + 3;

    }

    pxmlnode pxml = new xmlnode;
    pxml->nodeID = n_id;
    current = nxtPos;
    nxtPos = getAttributes(current, pxml->attributes);
    if (xmlIN[nxtPos] == '/')
    {
        pxmlelement x = new xmlelement('n', pxml);
        elmts.push_back(x);
        lstElm = x;
        return  nxtPos + 2;
    }
    current = nxtPos + 1;
    string c_tag = "</" + n_id + ">"; // the closing tag
    nxtPos = xmlIN.find(c_tag, current); assert (nxtPos != string::npos);
    if (current == nxtPos)
    { // empty element like <{n_id} [attributes]></{n_id}>
        string * pTxt = new string;
        *pTxt = "";
        pxmlelement w = new xmlelement('t', pTxt);
        pxml->elements.push_back(w);
        lstElm = w;
    }
    else
    {
        while (c_tag != xmlIN.substr(current, c_tag.size()))
        {
            nxtPos = addElement(current, pxml->elements);
            current = nxtPos;
        }
    }
    pxmlelement x = new xmlelement('n', pxml);
    elmts.push_back(x);
    lstElm = x;
    return current + c_tag.size();
}

size_t xmlHandler::getAttributes(const size_t pos, vector<string> & a)
{
    size_t current = pos, nxtPos = xmlIN.find_first_not_of(" \n\t\r", current), endTag = xmlIN.find('>', current);
    assert (endTag != string::npos);
    if (xmlIN[endTag - 1] == '/') --endTag;
    while (nxtPos < endTag)
    {
        current = nxtPos;
        nxtPos = xmlIN.find('=', current); assert ((nxtPos != string::npos) && (nxtPos < endTag));
        a.push_back(sOps.fullTrim(xmlIN.substr(current, nxtPos - current)));
        current = nxtPos + 1; nxtPos = xmlIN.find('"', current); assert((nxtPos != string::npos) && (nxtPos < endTag));
        current = nxtPos + 1; nxtPos = xmlIN.find('"', current); assert((nxtPos != string::npos) && (nxtPos < endTag));
        a.push_back(xmlIN.substr(current, nxtPos - current));
        current = nxtPos + 1; nxtPos = xmlIN.find_first_not_of(" \n\t\r", current);
    }
    return endTag;
}

void xmlHandler::LoadXML(const std::string &xmlFilePath)
{
    readXMLfile(xmlFilePath);
    size_t i_pos = 0, f_pos = 0;
    do
    {
        f_pos = addElement(i_pos, XML); // recursive
        i_pos = f_pos;
    } while (lstElm->eType != 'n');

    root = static_cast<pxmlnode>(lstElm->elm);

    //debug:
    //string xmlOUT;
    //for (vector<xmlelement>::iterator xIT = XML.begin(); xIT != XML.end(); ++xIT)
    //{
    //    xmlOUT += (*xIT).write_element(); // this is recursive
    //}
    //cout << xmlOUT << "\n";
    //end debug

}
/*
void xmlHandler::handleError(const string & err)
{
    cout << err << "\n";
    exit(EXIT_FAILURE);
}
*/
bool xmlHandler::getTextFromNode(const pxmlnode & pNode, string & node_text)
{
    node_text = "";
    if (pNode == NULL) return false;
    for (vector<pxmlelement>::iterator xit = pNode->elements.begin(); xit != pNode->elements.end(); ++xit)
    {
        pxmlelement & e = *xit;
        if (e->eType == 'n') return false;
        if (e->eType == 't')
        {
            string * pStr = static_cast<string *>(e->elm);
            node_text += stringops::fullTrim(*pStr);
        }
    }
    return true;
}

pxmlnode xmlHandler::getNode(const string & nID, vector<pxmlelement> & elms)
{
    for (vector<pxmlelement>::iterator xit = elms.begin(); xit != elms.end(); ++xit)
    {
        pxmlelement & e = *xit;
        if (e->eType == 'n')
        {
            pxmlnode pNode = static_cast<pxmlnode>(e->elm);
            if (pNode->nodeID.compare(nID) == 0) return pNode;
        }
    }
    return NULL;
}

vector<pxmlnode> xmlHandler::getAllNodes(const string & nID, vector<pxmlelement> & elms)
{
    vector<pxmlnode> allNodes;
    for (vector<pxmlelement>::iterator xit = elms.begin(); xit != elms.end(); ++xit)
    {
        pxmlelement & e = *xit;
        if (e->eType == 'n')
        {
            pxmlnode pNode = static_cast<pxmlnode>(e->elm);
            if (pNode->nodeID.compare(nID) == 0) allNodes.push_back(pNode);
        }
    }
    return allNodes;

}

pxmlnode xmlHandler::get_next_node(vector<pxmlelement>::iterator & nit, const vector<pxmlelement> & elms)
{
    for (; nit != elms.end(); ++nit)
    {
        pxmlelement & e = *nit;
        if (e->eType == 'n')
        {
            pxmlnode pNode = static_cast<pxmlnode>(e->elm);
            return pNode;
        }
    }
    return NULL;

}

pxmlnode xmlHandler::get_next_node(const string & nID, vector<pxmlelement>::iterator & nit, const vector<pxmlelement> & elms)
{
    for (; nit != elms.end(); ++nit)
    {
        pxmlelement & e = *nit;
        if (e->eType == 'n')
        {
            pxmlnode pNode = static_cast<pxmlnode>(e->elm);
            if (pNode->nodeID.compare(nID) == 0) return pNode;
        }
    }
    return NULL;

}
/*
void xmlHandler::handleError(const string & err)
{
    cout << err << "\n";
    exit(EXIT_FAILURE);
}
*/

vector<pxmlnode> xmlHandler::getAllNodes(vector<pxmlelement> & elms)
{
    vector<pxmlnode> allNodes;
    for (vector<pxmlelement>::iterator xit = elms.begin(); xit != elms.end(); ++xit)
    {
        pxmlelement & e = *xit;
        if (e->eType == 'n')
        {
            pxmlnode pNode = static_cast<pxmlnode>(e->elm);
            allNodes.push_back(pNode);
        }
    }
    return allNodes;

}

void xmlHandler::getTransform(void)
{
    string rtnVal = "";
    vector<pxmlelement> & elms = root->elements;
    for (vector<pxmlelement>::iterator xit = elms.begin(); xit != elms.end(); ++xit)
    {
        pxmlelement & e = *xit;
        if (e->eType == 'n')
        {
            pxmlnode pNode = static_cast<pxmlnode>(e->elm);
            if (pNode->nodeID.compare("transform") == 0) 
            {
                vector<pxmlnode> tNodes = getAllNodes(pNode->elements);
                for (vector<pxmlnode>::iterator tit = tNodes.begin(); tit != tNodes.end(); ++tit)
                {
                    pxmlnode & tNode = *tit;
                    if (tNode->nodeID.compare("translate") == 0)
                    {
                        CompGeo::XYZ txyz;
                        pxmlnode xNode = getNode("x", tNode->elements),
                            yNode = getNode("y", tNode->elements),
                            zNode = getNode("z", tNode->elements);

                        if (getTextFromNode(xNode, rtnVal)) txyz.x = atof(rtnVal.c_str());
                        else section::handleError(string("error getting x of translation part of ") + pNode->write_node());

                        if (getTextFromNode(yNode, rtnVal)) txyz.y = atof(rtnVal.c_str());
                        else section::handleError(string("error getting y of translation part of ") + pNode->write_node());

                        if (getTextFromNode(zNode, rtnVal)) txyz.z = atof(rtnVal.c_str());
                        else section::handleError(string("error getting z of translation part of ") + pNode->write_node());
                        
                        Model::translateT(txyz);
                    }
                    if (tNode->nodeID.compare("rotate") == 0)
                    {
                        char axis = ' ';
                        double angle = 0.0;
                        pxmlnode axisNode = getNode("axis", tNode->elements),
                            angleNode = getNode("angle", tNode->elements);

                        if (getTextFromNode(axisNode, rtnVal)) axis = rtnVal[0];
                        else section::handleError(string("error getting axis of rotation part of ") + pNode->write_node());

                        if (getTextFromNode(angleNode, rtnVal)) axis = atof(rtnVal.c_str());
                        else section::handleError(string("error getting angle of rotation part of ") + pNode->write_node());

                        Model::rotateT(axis, angle);
                    }
                }
            }
        }
    }

}


void xmlHandler::getSections(void)
{
    vector<pxmlelement> & elms = root->elements;
    for (vector<pxmlelement>::iterator xit = elms.begin(); xit != elms.end(); ++xit)
    {
        pxmlelement & e = *xit;
        if (e->eType == 'n')
        {
            pxmlnode pNode = static_cast<pxmlnode>(e->elm);
            if (pNode->nodeID.compare("section") == 0) sections.push_back(pNode);
        }
    }
    
}

void xmlHandler::buildSections(void)
{
    for (vector<pxmlnode>::iterator nit = sections.begin(); nit != sections.end(); ++nit)
    {
        pxmlnode & xNode = *nit;
        string n_type = "type";
        if (xNode->findAttribute(n_type))
        {
            if (n_type.compare("cylinder") == 0)  // nodeID is "section"
            {
                pcylinder pCyl = new cylinder;
                pCyl->buildSection(xNode);
                //psectionType psectTyp = new sectionType('c', pCyl);
                Model::s[pCyl->id] = pCyl;
            }
            if (n_type.compare("plane") == 0)
            {
                pplane pPl = new plane;
                pPl->buildSection(xNode);
                //psectionType psectTyp = new sectionType('p', pPl);
                Model::s[pPl->id] = pPl;

            }
            if (n_type.compare("patch") == 0)
            {
                ppatch pPtch = new patch;
                pPtch->buildSection(xNode);
                Model::s[pPtch->id] = pPtch;
            }
        }
    }
}

void xmlHandler::LoadAndParseXML(const string & pth)
{
    LoadXML(pth);
    getTransform();
    getSections();
    buildSections();

}

