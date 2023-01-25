#include "first.hpp"
#include "xmlHandler.hpp"
#include "Model.hpp"
#include "cylinder.hpp"
#include "plane.hpp"
#include "patch.hpp"

using namespace std;

//    char cylinderPART; // a = start_angle alpha, b = end_angle beta, f = first_circle, l = last_circle
map<char, string> patch::fullPART = 
    {pair('a', "start_angle"), pair('b', "end_angle"), pair('f', "first_circle"), pair('l', "last_circle")};

map<string, char> patch::idxPART = 
    {pair("end_angle", 'b'), pair("first_circle", 'f'), pair("last_circle", 'l'), pair("start_angle", 'a')};

patch::patch(void):id(""), cylinderID(""), cylinderPART(' '), start_triangles(0), total_triangles(0), orientation(0){}

patch::patch(const patch & p):id(p.id), cylinderID(p.cylinderID), cylinderPART(p.cylinderPART), fromVerts(p.fromVerts),
    oFacts(p.oFacts), sectsTo(p.sectsTo),
    start_triangles(p.start_triangles), total_triangles(p.total_triangles), orientation(p.orientation){}

patch::~patch(void) {}

void patch::getSection(const pxmlnode & sect)
{

    pxmlnode cfNode = xmlHandler::getNode("cylinder_from", sect->elements),
        cpNode = xmlHandler::getNode("cylinder_part", sect->elements);


    if ((cfNode == NULL) || (cpNode == NULL))
    {
        cout << "error(s) in patch section " << id << " missing";
        if (cfNode == NULL) cout << " cylinder_from ";
        if (cpNode == NULL) cout << " cylinder_part ";
        handleError(string("!\nHERE'S THE SECTION:\n") + sect->write_node());
    }

    string rtnVal = "";
    
    if (!xmlHandler::getTextFromNode(cfNode, cylinderID))
        handleError(string("error getting cylinder_from in ") + id + "\n" + sect->write_node());

    if (!xmlHandler::getTextFromNode(cpNode, rtnVal))
        handleError(string("error getting cylinder_part in ") + id + "\n" + sect->write_node());
    map<string, char>::iterator xit = idxPART.find(rtnVal);
    if (xit == idxPART.end()) 
        handleError(string("error: unknown cylinder part ") + rtnVal + " in " + id + "\n" + sect->write_node());
    cylinderPART = xit->second;

    fromVerts = Model::s[cylinderID]->getPart(rtnVal);

    vector<pxmlelement>::iterator nit = sect->elements.begin();
    while (nit != sect->elements.end())
    {
        pxmlnode tsNode = xmlHandler::get_next_node("section_to", nit, sect->elements);
        if (tsNode != NULL)
        {
            if (!xmlHandler::getTextFromNode(tsNode, rtnVal))
                handleError(string("error getting section_to from ") + id + "\n" + sect->write_node());
            
            string tsStr = rtnVal;

            vector<pxmlelement>::iterator pit = nit;
            pxmlnode psNode = xmlHandler::get_next_node("section_part", pit, sect->elements);
            if (psNode == NULL) 
                handleError(string("error: missing part for section ") + tsStr + " in " + id + "\n" + sect->write_node());
            if (!xmlHandler::getTextFromNode(psNode, rtnVal))
                handleError(string("error getting section_part for section ") + tsStr + " in " + id + "\n" + sect->write_node());

            string psStr = rtnVal;

            vector<unsigned int> tVts = Model::s[tsStr]->getPart(psStr);
            sectsTo.push_back(toType(tsStr, psStr, tVts));
            ++nit;
        }
    }

}

void patch::make_triangles(const vector<vector<unsigned int>> & verts, const vector<vector<double>> & orders)
{
    // triangle vertex sequences 2 counterclockwise  then 2 clockwise
    // the FROMs are A and B, the TOs are a and b
    start_triangles = Model::t.size();
    total_triangles = 0;
    unsigned int AaB[] = {0, 2, 1}, Bab[] = {1, 2, 3}, ABa[] = {0, 1, 2}, Bba[] = {1, 3, 2},
        * tSq[] = {AaB, Bab, ABa, Bba}, vtx_mat[2][2],
        orientation_offset = 1 - orientation, // 0 for counterclockwise, 2 for clockwise
        fromRow = 0, toRow = 1,
        vIdx[] = {0, 0},
        vSz[] = {static_cast<unsigned int>(verts[fromRow].size()), static_cast<unsigned int>(verts[toRow].size())},
        shift = 2, // 0=FROM, 1=TO, 2=both ==> shifted the vtx_mtx for this row
                    // use tSq[orientation_offset + shift] to make triangle
        lo_from = verts[fromRow][0], hi_from = verts[fromRow][vSz[fromRow] - 1],
        lo_to = verts[toRow][0], hi_to = verts[toRow][vSz[toRow] - 1];
    
    double ord_mat[2][2];

    if ((vSz[0] + vSz[1]) < 3) return; //handleError(string("error in patch ") + id + " not enough vertices for triangles.\n");
    // initial load:
    vtx_mat[fromRow][0] = lo_from; vtx_mat[fromRow][1] = 0; ord_mat[fromRow][0] = orders[fromRow][0]; ord_mat[fromRow][1] = 0.0;
    vtx_mat[toRow][0] = lo_to; 
    vtx_mat[toRow][1] = 0; 
    ord_mat[toRow][0] = orders[toRow][0]; 
    ord_mat[toRow][1] = 0.0;
    bool OnEnd[] = {vIdx[fromRow] == (vSz[fromRow] - 1), vIdx[toRow] == (vSz[toRow] - 1)};
    for (int i = 0; i < 2; ++i)
    {
        if (!OnEnd[i])
        {
            ++(vIdx[i]);
            vtx_mat[i][1] = verts[i][vIdx[i]];
            ord_mat[i][1] = orders[i][vIdx[i]];
            OnEnd[i] = vIdx[i] == (vSz[i] - 1);
        }
        else shift = (i + 1) % 2;
    }

    bool isStarting = true;

    while (!(OnEnd[0] && OnEnd[1]))
    {
        if (!isStarting)
        {
            // getting the next vertices based on the order of the last vertices in vtx_mat:
            shift = OnEnd[fromRow]? toRow: OnEnd[toRow]? fromRow: ord_mat[fromRow][1] <= ord_mat[toRow][1]? fromRow: toRow;
            ++(vIdx[shift]);
            // shift back:
            vtx_mat[shift][0] = vtx_mat[shift][1];
            ord_mat[shift][0] = ord_mat[shift][1];
            // shift in new:
            vtx_mat[shift][1] = verts[shift][vIdx[shift]];
            ord_mat[shift][1] = orders[shift][vIdx[shift]];
            OnEnd[shift] = vIdx[shift] == (vSz[shift] - 1);

            unsigned int other = (shift + 1) % 2;
            // shift back:
            vtx_mat[other][0] = vtx_mat[other][1];
            ord_mat[other][0] = ord_mat[other][1];
            
            if (!OnEnd[other] && (ord_mat[other][1] <= ord_mat[shift][1]))
            { // both from and to shall advance and we'll make 2 triangles
                ++(vIdx[other]);
                // shift back:
                //vtx_mat[other][0] = vtx_mat[other][1];
                //ord_mat[other][0] = ord_mat[other][1];
                // shift in new:
                vtx_mat[other][1] = verts[other][vIdx[other]];
                ord_mat[other][1] = orders[other][vIdx[other]];
                OnEnd[other] = vIdx[other] == (vSz[other] - 1);
                shift = 2;
            }
        }
        else isStarting = false;

        if ((shift == 0) || (shift == 2)) 
        { // make triangle using 2 FROM vertices:
            for (int i = 0; i < 3; ++i) 
            {
                unsigned int vtxIdx = tSq[orientation_offset + fromRow][i];
                Model::t.push_back(vtx_mat[vtxIdx / 2][vtxIdx % 2]);
            }
            ++total_triangles;
        }
        if ((shift == 1) || (shift == 2))
        { // make triangle using 2 TO vertices:
            for (int i = 0; i < 3; ++i)
            {
                unsigned int vtxIdx = tSq[orientation_offset + toRow][i];
                Model::t.push_back(vtx_mat[vtxIdx / 2][vtxIdx % 2]);
            }
            ++total_triangles;
        }
    }


}


void patch::add_arc(const psection & sct)
{
    pcylinder pCyl = dynamic_cast<pcylinder>(sct);
    
    assert ((pCyl->n2r == -1) || (pCyl->n2r == 1));
    assert ((cylinderPART == 'f') || (cylinderPART == 'l'));

    bool isConvex = pCyl->n2r == 1, 
        isForward = cylinderPART == 'l';

    orientation = isConvex? isForward? 1: -1: isForward? -1: 1;

    int arc_vtx_count = pCyl->n2r == 1? pCyl->N + 1: pCyl->N + 2;

    double ALPH = pCyl->alpha, BETA = pCyl->beta, da = BETA - ALPH, d0 = da / pCyl->N, delta = 0.0, 
        theta = ALPH;

    bool negsPresent = ALPH < 0.0;

    

    // getting radians for the angle from +IHAT to vertices on the arc:
    for (int i = 0; i < arc_vtx_count; ++i)
    {
        if (i > 0) 
        {
            delta = d0;
            if (pCyl->n2r == -1)
            {
                delta = d0;
                if ((i == 1) || (i == (arc_vtx_count - 1)))
                {
                    delta /= 2.0;
                }
            }
        }
        theta += delta;
        oFacts.push_back(theta);
    }

    CompGeo::XYZ norm= pCyl->NHAT, IHAT = pCyl->IHAT, JHAT= pCyl->JHAT, 
        origin = cylinderPART == 'f'? pCyl->ctr_first: pCyl->ctr_last;

    vector<unsigned int> to_verts;
    vector<double>to_orders;

    double oldOrder = 0.0, oldZ = 0.0;
    bool isStarting = true;

    for (vector<toType>::iterator tit = sectsTo.begin(); tit != sectsTo.end(); ++tit)
    {
        toType & tt = *tit;
        //to_verts.insert(to_verts.end(), tt.verts.begin(), tt.verts.end());
        for (vector<unsigned int>::iterator vit = tt.verts.begin(); vit != tt.verts.end(); ++vit)
        {
            to_verts.push_back(*vit);
            CompGeo::XYZ vtx = Model::v[*vit];
            vtx -= origin;
            double z = vtx * norm, x = vtx * IHAT, y = vtx * JHAT;
            bool xZero = fabs(x) < MAX_FLT_PRECISION, yZero = fabs(y) < MAX_FLT_PRECISION, zZero = fabs(z) < MAX_FLT_PRECISION;
            if (zZero) handleError(string("error: no offset in patch ") + id + " part " + tt.sPart + 
                " in section " + tt.sFrom + " has a vertex on the arc " + fullPART[cylinderPART] +
                " in cylinder " + cylinderID + "\n");
            if (xZero && yZero) handleError(string("error in patch ") + id + " part " + tt.sPart + 
                " in section " + tt.sFrom + " has a vertex projected to the origin of the arc " + fullPART[cylinderPART] +
                " in cylinder " + cylinderID + "\n");
            if (xZero || yZero)
            {
                //if (xZero) theta = y > 0.0? 0: M_PI;
                //else theta = x > 0.0? M_PI_2: 3 * M_PI_2;
                if (yZero) theta = x > 0.0? (isStarting || negsPresent)? 0: 2 * M_PI: M_PI;
                else theta = y > 0.0? M_PI_2: negsPresent? -M_PI_2: 3 * M_PI_2;
            }
            else
            {
                theta = atan(fabs(y) / fabs(x));
                if (y < 0.0)
                {
                    if (x < 0.0) theta += M_PI; // Quad III
                    else theta = 2 * M_PI - theta; // Quad IV
                    if (negsPresent) theta -= 2 * M_PI;
                }
                else
                {
                    if (x < 0.0) theta = M_PI - theta; // Quad III
                }
            }

            tt.orderFactors.push_back(theta);
            to_orders.push_back(theta);
            
            if (isStarting) isStarting = false;
            else
            {
                if ((oldZ * z) < 0.0) handleError(string("error in patch ") + id + " part " + tt.sPart + 
                    " in section " + tt.sFrom + " has an edge that intersects the arc " + fullPART[cylinderPART] +
                    " in cylinder " + cylinderID + "\n");
                if (theta <= oldOrder) handleError(string("error in patch ") + id + " part " + tt.sPart + 
                    " in section " + tt.sFrom + " is not radially monotonic relative to the arc " + fullPART[cylinderPART] +
                    " in cylinder " + cylinderID + "\n");
            }
            oldOrder = theta;
            oldZ = z;
        }
        //to_orders.insert(to_orders.end(), tt.orderFactors.begin(), tt.orderFactors.end());
    }
    vector<vector<unsigned int>> all_verts = {fromVerts, to_verts};
    vector<vector<double>> all_orders = {oFacts, to_orders};

    make_triangles(all_verts, all_orders);
}

void patch::sweep2line(const psection & sct)
{
    pcylinder pCyl = dynamic_cast<pcylinder>(sct);

    assert ((pCyl->n2r == -1) || (pCyl->n2r == 1));
    assert ((cylinderPART == 'a') || (cylinderPART == 'b'));

    bool isConvex = pCyl->n2r == 1, 
        isAlpha = cylinderPART == 'a';

    orientation = isConvex? isAlpha? 1: -1: isAlpha? -1: 1;

    // line eqn R=R0 + tL R0 will be line point closest to the model origin
    // orders will be t

    CompGeo::XYZ L = pCyl->axleHAT, Ri = Model::v[fromVerts[0]], R0 = Ri - (Ri * L) * L;
    bool xlZero = fabs(L.x) < MAX_FLT_PRECISION, ylZero = fabs(L.y) < MAX_FLT_PRECISION, zlZero = fabs(L.z) < MAX_FLT_PRECISION;

    for (vector<unsigned int>::iterator vit = fromVerts.begin(); vit != fromVerts.end(); ++ vit)
    {
        CompGeo::XYZ vtx = Model::v[*vit];
        double t = xlZero? ylZero? (vtx.z - R0.z) / L.z: (vtx.y - R0.y) / L.y: (vtx.x - R0.x) / L.x;
        oFacts.push_back(t); 
    }

    vector<unsigned int> to_verts;
    vector<double> to_orders;

    double oldOrder = 0.0;
    bool isStarting = true;

    for (vector<toType>::iterator tit = sectsTo.begin(); tit != sectsTo.end(); ++tit)
    {
        toType & tt = *tit;
        for (vector<unsigned int>::iterator vit = tt.verts.begin(); vit != tt.verts.end(); ++vit)
        {
            to_verts.push_back(*vit);
            CompGeo::XYZ vtx = Model::v[*vit];
            vtx -= R0;
            double t = vtx * L;
            tt.orderFactors.push_back(t);
            to_orders.push_back(t);
            if (isStarting) isStarting = false;
            else
            {
                if (t <= oldOrder) handleError(string("error in patch ") + id + " part " + tt.sPart + 
                    " in section " + tt.sFrom + " is unaligned with line at " + fullPART[cylinderPART] +
                    " of cylinder " + cylinderID + "\n");
            }
            oldOrder = t;
        }
    }
    vector<vector<unsigned int>> all_verts = {fromVerts, to_verts};
    vector<vector<double>> all_orders = {oFacts, to_orders};

    make_triangles(all_verts, all_orders);

}

vector<unsigned int> patch::getPart(const std::string & prt)
{ // prt = "lo_from", "lo_to", "hi_to", "hi_from", "lo", "hi"
    vector<unsigned int> v;

    assert (prt.size() >= 2);

    string pre = prt.substr(0, 2), suf = "";

    bool isLo = pre.compare("lo") == 0, isHi = pre.compare("hi") == 0, hasSfx = prt.size() > 2;

    if (hasSfx) suf = prt.substr(2);

    assert(suf[0] == '_');
    suf = suf.substr(1);

    bool isFrom = suf.compare("from") == 0, isTo = suf.compare("to") == 0;

    assert (isLo || isHi);
    if (hasSfx) assert (isFrom || isTo);

    unsigned int loIdx = 0, hiIdx = 1, fromIdx = 0, ToIdx = 1, fromPts[] = {fromVerts[0], fromVerts[fromVerts.size() - 1]},
        toPts[] = {sectsTo[0].verts[0], sectsTo[sectsTo.size() - 1].verts[sectsTo[sectsTo.size() - 1].verts.size() - 1]},
        * oPts[] = {fromPts, toPts}, 
        heightIdx = isLo? loIdx: hiIdx, 
        srcIdx = isFrom? fromIdx: ToIdx, 
        startAt = hasSfx? 0: isLo? orientation == 1? 0: 1: orientation == 1? 1: 0;

    if (hasSfx) v.push_back(oPts[srcIdx][heightIdx]);
    else
    {
        v.push_back(oPts[srcIdx][(heightIdx + startAt) % 2]);
        v.push_back(oPts[srcIdx][(heightIdx + startAt + 1) % 2]);
    }


    return v;
}

void patch::buildSection(const pxmlnode & sect)
{
    id = "id";
    if(!sect->findAttribute(id)) handleError(string("error: no id with section in plane\n") + sect->write_node());
    getSection(sect);
    assert ((cylinderPART == 'a') || (cylinderPART == 'b') || (cylinderPART == 'f') || (cylinderPART == 'l'));
    if ((cylinderPART == 'a') || (cylinderPART == 'b')) sweep2line(Model::s[cylinderID]);
    if ((cylinderPART == 'f') || (cylinderPART == 'l')) add_arc(Model::s[cylinderID]);
}

string patch::printSection(void)
{

    string r = "patch [" + id + "]\n" +
        "\ttriangles start:[" + to_string(start_triangles) + "]\ttotal triangles:[" + to_string(total_triangles) + "]" +
        "\tfrom cylinder:[" + cylinderID + "]\tpart:[" + fullPART[cylinderPART] + "]" +
        "\nFROM VERTEX INDICES / ORDER:\n=====================\n";
    
    for (unsigned int i = 0; i < fromVerts.size(); ++i)
    {
        r+= "\t[" + to_string(i) + "]: v" + to_string(fromVerts[i]) + " / o" + to_string(oFacts[i]);
    } 

    r+= "\nTO VERTEX INDICES / ORDER:\n\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\\/\n";

    unsigned int acc = 0;
    for (unsigned int i = 0; i < sectsTo.size(); ++i)
    {
        toType tt = sectsTo[i];
        r+=  "\t[" + to_string(i) + "] section " + tt.sFrom + " part " + tt.sPart + "\n\t";
        for (unsigned int j = 0; j < tt.verts.size(); ++j)
        {
            r+= "\t[" + to_string(j + acc) + "]: v" + to_string(tt.verts[j]) + " / o" + to_string(tt.orderFactors[j]);
        }
        acc += tt.verts.size();
        r+= "\n";
    }

    r += "END patch " + id + "\n==================\n\n";
 


    return r;
}

