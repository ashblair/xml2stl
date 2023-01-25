#include "first.hpp"
#include "xmlHandler.hpp"
#include "section.hpp"
#include "Model.hpp"
#include "cylinder.hpp"
#include "plane.hpp"

using namespace std;


plane::plane(void):multiplyer(0.0), triStart(0), triCount(0){}

plane::plane(const plane & pln):id(pln.id), norm(pln.norm), multiplyer(pln.multiplyer), 
    triStart(pln.triStart), triCount(pln.triCount), verts(pln.verts), offsets(pln.offsets), shapes(pln.shapes){}

plane::~plane(void){}

/*
bool plane::getTextFromNode(const pxmlnode & pNode, string & node_text)
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
            node_text += *pStr;
        }
    }
    return true;
}

pxmlnode plane::getNode(const string & nID, vector<pxmlelement> & elms)
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

vector<pxmlnode> plane::getAllNodes(const string & nID, vector<pxmlelement> & elms)
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

pxmlnode plane::get_next_node(vector<pxmlelement>::iterator & nit, const vector<pxmlelement> & elms)
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

pxmlnode plane::get_next_node(const string & nID, vector<pxmlelement>::iterator & nit, const vector<pxmlelement> & elms)
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

void plane::handleError(const string & err)
{
    cout << err << "\n";
    exit(EXIT_FAILURE);
}
*/
/*
CompGeo::XYZ plane::getNamedVector(const string & vect)
{
    string v = vect;
    assert (v.size() > 0);
    char first = v[0];
    if ((first == '+') || (first == '-')) v = v.substr(1);
    bool isNeg = first == '-';
    assert (v.size() == 1);
    char vId = v[0];
    CompGeo::XYZ our_norm;
    assert ((vId == 'i') || (vId == 'j') || (vId == 'k'));
    switch (vId)
    {
        case 'i': our_norm = CompGeo::XYZ(1.0, 0.0, 0.0); break;
        case 'j': our_norm = CompGeo::XYZ(0.0, 1.0, 0.0); break;
        case 'k': our_norm = CompGeo::XYZ(0.0, 0.0, 1.0); break;
    }
    if (isNeg) our_norm = -our_norm;
    return our_norm;
}
*/
void plane::getSection(const pxmlnode & sect)
{ // who wants to edit check?

    pxmlnode nNode = xmlHandler::getNode("norm", sect->elements),
        mNode = xmlHandler::getNode("multiplyer", sect->elements);

    vector<pxmlnode> vNodes = xmlHandler::getAllNodes("verts", sect->elements);

    if ((nNode == NULL) || (mNode == NULL) || (vNodes.size() == 0))
    {
        cout << "error(s) in plane section " << id << " missing";
        if (nNode == NULL) cout << " norm ";
        if (mNode == NULL) cout << " multiplyer ";
        if (vNodes.size() == 0) cout << " verts ";
        handleError(string("!\nHERE'S THE SECTION:\n") + sect->write_node());
    }

    string rtnVal = "";
    
    if (xmlHandler::getTextFromNode(nNode, rtnVal)) norm = getNamedVector(rtnVal);
    else 
    {
        pxmlnode xNode = xmlHandler::getNode("x", nNode->elements),
            yNode = xmlHandler::getNode("y", nNode->elements),
            zNode = xmlHandler::getNode("z", nNode->elements);
        
        if (xmlHandler::getTextFromNode(xNode, rtnVal)) norm.x = atof(rtnVal.c_str());
        else handleError(string("error getting x of norm in ") + id + "\n" + sect->write_node());

        if (xmlHandler::getTextFromNode(yNode, rtnVal)) norm.y = atof(rtnVal.c_str());
        else handleError(string("error getting y of norm in ") + id + "\n" + sect->write_node());

        if (xmlHandler::getTextFromNode(zNode, rtnVal)) norm.z = atof(rtnVal.c_str());
        else handleError(string("error getting z of norm in ") + id + "\n" + sect->write_node());

        if (fabs(norm.GetMagnitude() - 1.0) > MAX_FLT_PRECISION) handleError(string("error: norm not normalized in ") + 
            id + "\n" + sect->write_node());
    }


    if (xmlHandler::getTextFromNode(mNode, rtnVal)) multiplyer = atof(rtnVal.c_str());
    else handleError(string("error getting multiplyer in ") + id + "\n" + sect->write_node());

    int vertsCount = 0;
    shapes.push_back(0);
    for (vector<pxmlnode>::iterator vit = vNodes.begin(); vit != vNodes.end(); ++vit)
    {
        offsets.push_back(verts.size());
        pxmlnode & vNode = *vit;
        string aStr = "type";
        if (vNode->findAttribute(aStr))
        {
            if ((verts.size() > 0) && (aStr.compare("first") == 0))
            {
                shapes.push_back(verts.size());
            }
        }
        vector<pxmlelement> e = vNode->elements;
        vector<pxmlelement>::iterator eit = e.begin();
        pxmlnode eNode = xmlHandler::get_next_node(eit, e);
        while (eNode != NULL)
        {
            if (eNode->nodeID.compare("lookup") == 0)
            {
                vector<pxmlelement> luElms = eNode->elements;
                vector<pxmlelement>::iterator luit = luElms.begin();

                pxmlnode sfNode = xmlHandler::get_next_node("section_from", luit, luElms),
                    spNode = xmlHandler::get_next_node("section_part", luit, luElms);

                while ((sfNode != NULL) && (spNode != NULL))
                {
                    string sFrom = "", sPart = "";
                    
                    if (!xmlHandler::getTextFromNode(sfNode, sFrom)) 
                    handleError(string("error getting section_from in verts [") + to_string(vertsCount) + "] in plane " + 
                        id + "\n" + sect->write_node());  

                    if (!xmlHandler::getTextFromNode(spNode, sPart)) 
                    handleError(string("error getting section_part in verts [") + to_string(vertsCount) + "] in plane " + 
                        id + "\n" + sect->write_node()); 

                    psection ps = Model::s[sFrom];
                    if (ps == NULL) handleError(string("error in plane [") + id + "] tried to get a seam from unknown section [" + sFrom + "]");
                    //void * vdPtr = psT->sect;
                    //char s = psT->sType;
                    //pplane pPl = NULL;
                    //pcylinder pCyl = NULL;
                    vector<unsigned int> vtxIn = ps->getPart(sPart);
                    /*
                    switch(s)
                    {
                        case 'c': 
                            pCyl = static_cast<pcylinder>(vdPtr);
                            vtxIn = pCyl->getPart(sPart);
                            break;
                        case 'p':
                            pPl = static_cast<pplane>(vdPtr);
                            vtxIn = pPl->getPart(sPart);
                            break;
                    }
                    */
                    int skipFirst = verts.size() == 0? 0: verts[verts.size() - 1] == vtxIn[0]? 1: 0;
                    verts.insert(verts.end(), vtxIn.begin() + skipFirst, vtxIn.end()); 

                    ++luit;
                    sfNode = luit != luElms.end()? xmlHandler::get_next_node("section_from", luit, luElms): NULL;
                    spNode = luit != luElms.end()? xmlHandler::get_next_node("section_part", luit, luElms): NULL;

                }
            }
            if (eNode->nodeID.compare("vtx") == 0)
            {
                pxmlnode vtxNode = eNode;
                vector<pxmlelement> vtxElms = vtxNode->elements;
                vector<pxmlelement>::iterator vtxit = vtxElms.begin();
                pxmlnode xNode = xmlHandler::get_next_node("x", vtxit, vtxElms),
                    yNode = xmlHandler::get_next_node("y", vtxit, vtxElms),
                    zNode = xmlHandler::get_next_node("z", vtxit, vtxElms);
                
                while ((xNode != NULL) && (yNode != NULL) && (zNode != NULL))
                {
                    CompGeo::XYZ vtxIn;

                    if (xmlHandler::getTextFromNode(xNode, rtnVal)) vtxIn.x = atof(rtnVal.c_str());
                    else handleError(string("error getting x of vtx in verts[") + to_string(vertsCount) + "] in plane " + 
                        id + "\n" + sect->write_node());

                    if (xmlHandler::getTextFromNode(yNode, rtnVal)) vtxIn.y = atof(rtnVal.c_str());
                    else handleError(string("error getting y of vtx in verts[") + to_string(vertsCount) + "] in plane " + 
                        id + "\n" + sect->write_node());

                    if (xmlHandler::getTextFromNode(zNode, rtnVal)) vtxIn.z = atof(rtnVal.c_str());
                    else handleError(string("error getting z of vtx in verts[") + to_string(vertsCount) + "] in plane " + 
                        id + "\n" + sect->write_node());

                    unsigned int idx = Model::v.size();
                    Model::v.push_back(vtxIn);
                    verts.push_back(idx);
                    ++vtxit;
                    xNode = vtxit != vtxElms.end()? xmlHandler::get_next_node("x", vtxit, vtxElms): NULL;
                    yNode = vtxit != vtxElms.end()? xmlHandler::get_next_node("y", vtxit, vtxElms): NULL;
                    zNode = vtxit != vtxElms.end()? xmlHandler::get_next_node("z", vtxit, vtxElms): NULL;
                }
            }
            ++eit;
            eNode = eit != e.end()? xmlHandler::get_next_node(++eit, e): NULL;
        }
        ++vertsCount;
    }
    if (verts[0] == verts[verts.size() - 1]) verts.pop_back();
}

vector<unsigned int> plane::getPart(const string & prt)
{ // plane parts are verts nodes indexed by xml order, or id's subnodes 
    // prt should be that index as verts[index] or the id of the part

    assert (prt.size() > 0); 
    string p = prt;
    char first = p[0];
    bool indicatorPresent = (first == '+') || (first == '-'); 
    if (indicatorPresent) p = p.substr(1);
    bool reversePart = first == '-';


    int vIdx = atoi(prt.c_str()), numPrtsAvail = offsets.size();
    assert (numPrtsAvail > vIdx);
    vector<unsigned int>::iterator pit1 = verts.begin() + offsets[vIdx], 
        pit2 = vIdx == (numPrtsAvail - 1)? verts.end(): verts.begin() + offsets[vIdx + 1];
    
    vector<unsigned int> vPrt;
    vPrt.insert(vPrt.end(), pit1, pit2);

    if (reversePart) reverse(vPrt.begin(), vPrt.end());

    return vPrt;

}


void plane::triangularize(void)
{
    //cout << "in triangularize for plane " << id << " " << to_string(verts.size()) << " points:\n";

    unsigned int vSz = verts.size();
    triStart = Model::t.size();

    if (vSz == 3)
    {
        Model::t.insert(Model::t.end(), verts.begin(), verts.end());
        triCount = 1;
        return;
    }

    vector<CompGeo::XYZ> vertices;

    string vts = "";
    unsigned int sMax = 0, sIdx = 0, topIdx = verts[0], i = 0;
    while (i < vSz)
    {
        unsigned int idx = verts[i];
        CompGeo::XYZ P = Model::v[idx];
        bool addMe = true;
        if (i == sMax)
        {
            ++sIdx;
            if (shapes.size() == sIdx) sMax = vSz;
            else sMax = shapes[sIdx];
            topIdx = idx;
        }
        else
        {
            if (idx == topIdx)
            {
                addMe = false;
                vector<unsigned int>::iterator uit = verts.begin() + i;
                verts.erase(uit);
                --sMax;
                --vSz;
                for(unsigned int j = sIdx; j < shapes.size(); ++j)
                {
                    --shapes[j];
                }
            }
        }
        //cout << "vtx["<< to_string(idx) << "]: ";
        //P.PrintXYZ();
        //cout << "\n";
        if (addMe)
        {
            vertices.push_back(P);
            //vts += P.toStr("");
            ++i;
        }
    }
    NGons ng;
    ng.Translate3DPlanePoints(norm, vertices, shapes, vts);
    string eMsg = "";
    Triangulation tri = Triangulation(ng, false, eMsg);
    if (eMsg.size() != 0) handleError(string("error: triangularization for plane <") + id + "> failed - " + eMsg + "\nVERTICES: " + vts);
    vector<unsigned int> tri_idx = tri.indexFaces(tri.tLst.Faces);

    //cout << "triangles:";

    for (unsigned int i = 0; i < tri_idx.size(); ++i)
    {
        tri_idx[i] = verts[tri_idx[i]];

        //if ((i % 3) == 0) cout << "\n[" << to_string(i / 3) << "]\t";
        //cout << to_string(tri_idx[i]) << "\t";
    }
    
    //cout << "\n";
    triCount = tri_idx.size() / 3;
    Model::t.insert(Model::t.end(), tri_idx.begin(), tri_idx.end());
}

void plane::buildSection(const pxmlnode & sect)
{
    id = "id";
    if(!sect->findAttribute(id)) handleError(string("error: no id with section in plane\n") + sect->write_node());
    getSection(sect);
    triangularize();
}

//string plane::printXYZ(const string & label, const CompGeo::XYZ & xyz)
//{
//    return string((label.size() > 0?string("\t" + label + ":("): string(" (")) + 
//        to_string(xyz.x) + ", " + to_string(xyz.y) + ", " + to_string(xyz.z) + ")");
//}

string plane::printSection(void)
{
    //std::string id;
    //CompGeo::XYZ norm;
    //double multiplyer; // this times norm gives origin on plane
    //unsigned int triStart, triCount;
    //std::vector<unsigned int> verts, offsets;

    string r = "plane [" + id + "]\n" +
        norm.toStr("norm") + "\tmultiplyer:[" + to_string(multiplyer) + "]" +
        "\ttriStart:[" + to_string(triStart) + "]\ttriCount:[" + to_string(triCount) + "]" +
        "\nPLANE VERTEX INDICES:\n=====================\n";

    unsigned int numVerts = verts.size(), numOffs = offsets.size(), current_start = 0, current_size = numVerts;
    assert (numOffs > 0);
    assert (offsets[0] == 0);

    while (current_start < numOffs)
    {
        r += "\t<verts>[" + to_string(current_start) + "]:\t{";

        unsigned int vtx_acc = offsets[current_start], 
            vtx_acc_next = current_start == (numOffs - 1)? numVerts: offsets[current_start + 1];

        current_size = vtx_acc_next - vtx_acc;

        for (unsigned int i = 0; i < current_size; ++i)
        {
            r += to_string(verts[i + vtx_acc]);
            if (i < (current_size - 1)) r += ", ";
        }

        r += "}\n";

        ++current_start;
    }

    r += "END plane " + id + "\n==================\n\n";
    return r;

}