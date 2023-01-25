#include "first.hpp"
#include "xmlHandler.hpp"
#include "section.hpp"
#include "Model.hpp"
#include "cylinder.hpp"
#include "plane.hpp"


using namespace std;


cylinder::cylinder(void):radius(0.0), error(0.0), alpha(0.0), beta(2.0 * M_PI), axle_length(0.0), axle_delta(0.0), 
    n2r(0), N(0), total_arcs(0), total_triangles(0), start_triangles(0) {}

cylinder::cylinder(const cylinder & c):id(c.id),
    radius(c.radius), error(c.error), alpha(c.alpha), beta(c.beta), axle_length(c.axle_length), axle_delta(c.axle_delta),
    n2r(c.n2r), N(c.N), total_arcs(c.total_arcs), total_triangles(c.total_triangles), start_triangles(c.start_triangles),
    ctr_first(c.ctr_first), ctr_last(c.ctr_last), axleHAT(c.axleHAT), NHAT(c.NHAT), IHAT(c.IHAT), JHAT(c.JHAT), 
    arcs(c.arcs) {}
 
cylinder::~cylinder(void) {}

/*
bool cylinder::getTextFromNode(const pxmlnode & pNode, string & node_text)
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

pxmlnode cylinder::getNode(const string & nID, vector<pxmlelement> & elms)
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

void cylinder::handleError(const string & err)
{
    cout << err << "\n";
    exit(EXIT_FAILURE);
}
*/

void cylinder::getSection(const pxmlnode & sect)
{ // who wants to edit check?

    //map<char, seamType> seams;
    pxmlnode rNode = xmlHandler::getNode("radius", sect->elements),
        nNode = xmlHandler::getNode("norm", sect->elements),
        eNode = xmlHandler::getNode("error", sect->elements),
        fNode = xmlHandler::getNode("first_circle", sect->elements),
        lNode = xmlHandler::getNode("last_circle", sect->elements),
        sNode = xmlHandler::getNode("start_angle", sect->elements),
        dNode = xmlHandler::getNode("end_angle", sect->elements);
    
    if ((rNode == NULL) || (nNode == NULL) || (eNode == NULL) || (fNode == NULL) || (lNode == NULL))
    {
        cout << "error(s) in cylinder section " << id << " missing";
        if (rNode == NULL) cout << " radius ";
        if (nNode == NULL) cout << " norm ";
        if (eNode == NULL) cout << " error ";
        if (fNode == NULL) cout << " first_circle ";
        if (lNode == NULL) cout << " last_circle ";
        cout << "!\nHERE'S THE SECTION:\n" << sect->write_node();
        exit(EXIT_FAILURE);
    }

    string rtnVal = "";
    
    if (xmlHandler::getTextFromNode(rNode, rtnVal)) radius = atof(rtnVal.c_str());
    else handleError(string("error getting radius in ") + id + "\n" + sect->write_node());
    if (radius <= 0.0) handleError(string("error: radius must be > 0 in ") + id + "\n" + sect->write_node());

    if (xmlHandler::getTextFromNode(nNode, rtnVal)) n2r = atoi(rtnVal.c_str());
    else handleError(string("error getting norm in ") + id + "\n" + sect->write_node());
    if (!((n2r == 1) || (n2r == -1))) handleError(string("error: norm must be [+]1 or -1 in ") + id + "\n" + sect->write_node());

    if (xmlHandler::getTextFromNode(eNode, rtnVal)) error = atof(rtnVal.c_str());
    else handleError(string("error getting max error in ") + id + "\n" + sect->write_node());
    if (error <= 0.0) handleError(string("error: max error must be > 0 in ") + id + "\n" + sect->write_node());
    if (radius <= error) handleError(string("error: max error must be < radius in ") + id + "\n" + sect->write_node());

    if (sNode != NULL)
    { // start_angle -> alpha
        if (xmlHandler::getTextFromNode(sNode, rtnVal)) alpha = atof(rtnVal.c_str());
        else 
        {
            pxmlnode vNode = xmlHandler::getNode("value", sNode->elements),
                mNode = xmlHandler::getNode("seam", sNode->elements);
            
            if (xmlHandler::getTextFromNode(vNode, rtnVal)) alpha = atof(rtnVal.c_str());
            else handleError(string("error getting start_angle in ") + id + "\n" + sect->write_node());

            if (mNode != NULL)
            {
                seamType sT;

                pxmlnode sfNode = xmlHandler::getNode("section_from", mNode->elements),
                    spNode = xmlHandler::getNode("section_part", mNode->elements);

                if (xmlHandler::getTextFromNode(sfNode, rtnVal)) sT.from = rtnVal;
                else handleError(string("error getting start_angle > seam > section_from in ") + id + "\n" + sect->write_node());

                if (xmlHandler::getTextFromNode(spNode, rtnVal)) sT.part = rtnVal;
                else handleError(string("error getting start_angle > seam > section_part in ") + id + "\n" + sect->write_node());

                seams['a'] = sT;
            }
        }
    }
    if (dNode != NULL)
    { // end_angle -> beta
        if (xmlHandler::getTextFromNode(dNode, rtnVal)) beta = atof(rtnVal.c_str());
        else 
        {
            pxmlnode vNode = xmlHandler::getNode("value", dNode->elements),
                mNode = xmlHandler::getNode("seam", dNode->elements);
            
            if (xmlHandler::getTextFromNode(vNode, rtnVal)) beta = atof(rtnVal.c_str());
            else handleError(string("error getting end_angle in ") + id + "\n" + sect->write_node());

            if (mNode != NULL)
            {
                seamType sT;

                pxmlnode sfNode = xmlHandler::getNode("section_from", mNode->elements),
                    spNode = xmlHandler::getNode("section_part", mNode->elements);

                if (xmlHandler::getTextFromNode(sfNode, rtnVal)) sT.from = rtnVal;
                else handleError(string("error getting end_angle > seam > section_from in ") + id + "\n" + sect->write_node());

                if (xmlHandler::getTextFromNode(spNode, rtnVal)) sT.part = rtnVal;
                else handleError(string("error getting end_angle > seam > section_part in ") + id + "\n" + sect->write_node());

                seams['b'] = sT;
            }
        }
    }

    if (beta <= alpha) handleError(string("error: start_angle must be < end_angle in ") + id + "\n" + sect->write_node());
    if (alpha < -M_PI) handleError(string("error: start_angle must be > -pi in ") + id + "\n" + sect->write_node());
    if (beta > (2 * M_PI)) handleError(string("error: end_angle must be < 2pi in ") + id + "\n" + sect->write_node());
    if ((beta - alpha) > (2 * M_PI)) handleError(string("error: 1 rotation max in ") + id + "\n" + sect->write_node());
    //circles:

    pxmlnode fcNode = xmlHandler::getNode("ctr", fNode->elements);

    if (fcNode == NULL) handleError(string("error getting first_circle > ctr in ") + id + "\n" + sect->write_node());
    else
    {
        pxmlnode xNode = xmlHandler::getNode("x", fcNode->elements),
            yNode = xmlHandler::getNode("y", fcNode->elements),
            zNode = xmlHandler::getNode("z", fcNode->elements);

        if (xmlHandler::getTextFromNode(xNode, rtnVal)) ctr_first.x = atof(rtnVal.c_str());
        else handleError(string("error getting first_circle > ctr > x in ") + id + "\n" + sect->write_node());

        if (xmlHandler::getTextFromNode(yNode, rtnVal)) ctr_first.y = atof(rtnVal.c_str());
        else handleError(string("error getting first_circle > ctr > y in ") + id + "\n" + sect->write_node());

        if (xmlHandler::getTextFromNode(zNode, rtnVal)) ctr_first.z = atof(rtnVal.c_str());
        else handleError(string("error getting first_circle > ctr > z in ") + id + "\n" + sect->write_node());
        
    }

    pxmlnode fmNode = xmlHandler::getNode("seam", fNode->elements);

    if (fmNode != NULL)
    {
        seamType sT;

        pxmlnode sfNode = xmlHandler::getNode("section_from", fmNode->elements),
            spNode = xmlHandler::getNode("section_part", fmNode->elements);

        if (xmlHandler::getTextFromNode(sfNode, rtnVal)) sT.from = rtnVal;
        else handleError(string("error getting first_circle > seam > section_from in ") + id + "\n" + sect->write_node());

        if (xmlHandler::getTextFromNode(spNode, rtnVal)) sT.part = rtnVal;
        else handleError(string("error getting first_circle > seam > section_part in ") + id + "\n" + sect->write_node());

        seams['f'] = sT;

    }

    pxmlnode lcNode = xmlHandler::getNode("ctr", lNode->elements);

    if (lcNode == NULL) handleError(string("error getting last_circle > ctr in ") + id + "\n" + sect->write_node());
    else
    {
        pxmlnode xNode = xmlHandler::getNode("x", lcNode->elements),
            yNode = xmlHandler::getNode("y", lcNode->elements),
            zNode = xmlHandler::getNode("z", lcNode->elements);

        if (xmlHandler::getTextFromNode(xNode, rtnVal)) ctr_last.x = atof(rtnVal.c_str());
        else handleError(string("error getting last_circle > ctr > x in ") + id + "\n" + sect->write_node());

        if (xmlHandler::getTextFromNode(yNode, rtnVal)) ctr_last.y = atof(rtnVal.c_str());
        else handleError(string("error getting last_circle > ctr > y in ") + id + "\n" + sect->write_node());

        if (xmlHandler::getTextFromNode(zNode, rtnVal)) ctr_last.z = atof(rtnVal.c_str());
        else handleError(string("error getting last_circle > ctr > z in ") + id + "\n" + sect->write_node());
        
    }

    pxmlnode lmNode = xmlHandler::getNode("seam", lNode->elements);

    if (lmNode != NULL)
    {
        seamType sT;

        pxmlnode sfNode = xmlHandler::getNode("section_from", lmNode->elements),
            spNode = xmlHandler::getNode("section_part", lmNode->elements);

        if (xmlHandler::getTextFromNode(sfNode, rtnVal)) sT.from = rtnVal;
        else handleError(string("error getting last_circle > seam > section_from in ") + id + "\n" + sect->write_node());

        if (xmlHandler::getTextFromNode(spNode, rtnVal)) sT.part = rtnVal;
        else handleError(string("error getting last_circle > seam > section_part in ") + id + "\n" + sect->write_node());

        seams['l'] = sT;

    }

    if (ctr_first == ctr_last) handleError(string("error: circle centers must not be equal in ") + id + "\n" + sect->write_node());

    //return seams;
}

void cylinder::calculateOtherParameters(void)
{

    CompGeo::XYZ AXLE = ctr_last - ctr_first;
    axle_length = AXLE.GetMagnitude();
    axleHAT = AXLE / axle_length;
    NHAT = -axleHAT; // normal vector corresponding to k
    double w = sqrt(sqr(NHAT.x) + sqr(NHAT.z));
    bool wIsZero = fabs(w) < MAX_FLT_PRECISION;
    double X = wIsZero? 0.0: NHAT.x / w, Z = wIsZero? 1.0: NHAT.z / w;
    // basis vectors corresponding to i and j
    IHAT = CompGeo::XYZ(Z, 0.0, -X);
    JHAT = CompGeo::XYZ(-NHAT.y * X, w, -NHAT.y * Z);
    // the polygonal approximation to the convex circular arc from alpha to beta
    // will use equal length segments w/ endpoints on the circle to span an integral (N) number
    // of subdivivisions of (beta - alpha) on the even arcs {0, 2, 4, ...}
    // if N > 1 the odd arcs will be staggered so that the first & last segments
    // will be half length from peak to peak
    // for concave arcs the segment midpoints will be on the circle
    // except for first and last which will start or end on the circle 
    // for the odd arcs the first and last segments need to turn at (beta - alpha) / 4
    // odd arcs should have one more vertex than even arcs
    // convex even arcs should have N + 1 vertices
    // concave even arcs should have N + 2 vertices
    // calculating the minimum number of sides of a polygonal approximation to 
    // an arc within our error
    double da = beta - alpha;
    N = static_cast<int>(ceil(da / (2.0 * atan(n2r == 1? (sqrt(2.0 * error * radius - sqr(error)) / (radius - error)): 
        (sqrt(2.0 * error * radius + sqr(error)) / radius)))));
    if (N < 2) ++N; // need at least 2 to stagger
    total_arcs = 1 + static_cast<int>(ceil(axle_length / error));
    if ((total_arcs % 2) == 0) ++total_arcs; // guarantee an even final arc
    axle_delta = axle_length / (total_arcs - 1); // this is <= error
}

vector<unsigned int> cylinder::getPart(const string & prt)
{
    assert (prt.size() > 0); 
    string p = prt;
    char first = p[0];
    bool indicatorPresent = (first == '+') || (first == '-'); // || (first == '@') || (first == '^');
    if (indicatorPresent) p = p.substr(1);
    bool reversePart = first == '-', 
        //firstOnly = first == '@',
        //lastOnly = first == '^',
        at_alpha = p.compare("start_angle") == 0, 
        at_beta = p.compare("end_angle") == 0,
        at_first = p.compare("first_circle") == 0,
        at_last = p.compare("last_circle") == 0;
    assert (at_alpha || at_beta || at_first || at_last);
    vector<unsigned int> vtx;
    int even_vtx = n2r == 1? N + 1: N + 2,
        odd_vtx = n2r == 1? N + 2: N + 3;

    if (at_alpha || at_beta)
    {
        int acc = 0;
        for (int i = 0; i < total_arcs; ++i)
        {
            int delta_acc = (i % 2) == 0? even_vtx: odd_vtx;
            if (at_alpha) vtx.push_back(arcs[acc]);
            if (at_beta) vtx.push_back(arcs[acc + delta_acc - 1]);
            acc += delta_acc;
        }
    }
    if (at_first || at_last)
    {

        int delta_acc = at_first? even_vtx: (total_arcs % 2) == 0? odd_vtx: even_vtx, 
            acc = at_first? 0: arcs.size() - delta_acc; 
        // arcs.size() s/b (total_arcs / 2) * (even_vtx + odd_vtx) + ((total_arcs % 2) == 0? 0: even_vtx)
        for (int i = acc; i < (acc + delta_acc); ++i)
        {
            vtx.push_back(arcs[i]);
        }
    }
    if (reversePart) reverse(vtx.begin(), vtx.end());

    //cout << "cylinder " << id << " getPart(" << prt << ") returning:\n";
    //for (vector<unsigned int>::iterator vit = vtx.begin(); vit != vtx.end(); ++vit)
    //{
    //    cout << to_string(*vit) << " ";
    //}
    //cout << "\n";
    
    return vtx;
}

bool cylinder::handleSeam(const int & i, unsigned int & idx, const std::vector<unsigned int> & sm, 
    const CompGeo::XYZ & calcd, const std::string & descr)
{

    //CompGeo::XYZ P = origin + evens[j];

    //cout << "arc[" << to_string(i) << "], vtx[" << to_string(j) << "] calc'd: "; P.PrintXYZ();

    idx = sm[i];
    CompGeo::XYZ A = Model::v[idx], D = calcd - A;


    double d = D.GetMagnitude();
    if (error < d) 
    {
        cout << "\taseam point: "; A.PrintXYZ();
        cout << "\tcalculated point: "; calcd.PrintXYZ();
        cout << "\t too large difference: " << to_string(d) << "\n";

        handleError(string("error: ") + descr + " seam too far from arc in [" + id + "]");
    }
    return false; // indicating that index is not new
}


void cylinder::calculateArcs(void) //(map<char, seamType> sms)
{
    calculateOtherParameters();
    double da = beta - alpha, d0 = da / N, delta = 0.0, theta = alpha, r = radius,
        // these next radii are for concave (n2r == -1) arcs
        R = radius / cos(d0 /  2.0), R_I = radius / cos(d0 / 4.0);
    bool full_rotation = fabs(da - (2 * M_PI)) < MAX_FLT_PRECISION;
    int even_vtx = n2r == 1? N + 1: N + 2,
        odd_vtx = n2r == 1? N + 2: N + 3;

    vector<CompGeo::XYZ> evens, odds;

    // building templates for even & odd arcs:
    for (int i = 0; i < even_vtx; ++i)
    {
        if (i > 0) 
        {
            delta = d0;
            if (n2r == -1)
            {
                r = R;
                delta = d0;
                if ((i == 1) || (i == (even_vtx - 1)))
                {
                    delta /= 2.0;
                    if (i == (even_vtx - 1)) r = radius;
                }
            }
        }
        theta += delta;
        CompGeo::XYZ P = (r * cos(theta)) * IHAT + (r * sin(theta)) * JHAT;
        evens.push_back(P);
    }

    theta = alpha; r = radius; delta = 0.0;
    for (int i = 0; i < odd_vtx; ++i)
    { // this is the staggered (odd) arc
        if (i > 0) 
        {
            delta = d0;
            if (n2r == 1)
            { // convex arc
                if ((i == 1) || (i == (odd_vtx - 1))) delta /= 2.0;
            }
            else //if (n2r == -1)
            { // concave arc
                r = R;
                if ((i == 1) || (i == 2) || (i == (odd_vtx - 1)) || (i == (odd_vtx - 2)))
                {
                    if ((i == 1) || (i == (odd_vtx - 1))) delta /= 4.0;
                    else delta *= 0.75;

                    if (i == (odd_vtx - 1)) r = radius;
                    else r = R_I;

                    if ((i == 2) && ((odd_vtx - 2) > 2)) r = R;
                }

            }
        }
        theta += delta;
        CompGeo::XYZ P = (r * cos(theta)) * IHAT + (r * sin(theta)) * JHAT;
        odds.push_back(P);
    }

    // getting seams:
    vector<unsigned int> aseam, bseam, fseam, lseam;
    char smType[] = {'a', 'b', 'f', 'l'};

    for (int i = 0; i < sizeof(smType); ++i)
    {
        vector<unsigned int> current_seam;
        map<char, seamType>::const_iterator sit = seams.find(smType[i]);
        if (sit != seams.end())
        {
            seamType smT = sit->second;
            psection ps = Model::s[smT.from];
            if (ps == NULL) handleError(string("error in cylinder [") + id + "] tried to get a seam from unknown section [" + smT.from + "]");
            current_seam = ps->getPart(smT.part);
            /*
            char sTyp = ps->sType;
            void * v = ps->sect;
            pcylinder pcyl = NULL;
            pplane ppl = NULL;
            switch(sTyp)
            {
                case 'c':
                    pcyl = static_cast<pcylinder>(v);
                    current_seam = pcyl->getPart(smT.part);
                    break;
                case 'p':
                    ppl = static_cast<pplane>(v);
                    current_seam = ppl->getPart(smT.part);
                    break;
            }
            */
            if ((smType[i] == 'a') || (smType[i] == 'b')) assert (current_seam.size() == total_arcs);
            else assert (current_seam.size() == even_vtx);

            switch(smType[i])
            {
                case 'a': aseam = current_seam; break;
                case 'b': bseam = current_seam; break;
                case 'f': fseam = current_seam; break;
                case 'l': lseam = current_seam; break;
            }

        }

    }
    // debug:
    //cout << "calculating arcs in cylinder " << id << ":\n";


    CompGeo::XYZ origin = ctr_first;
    unsigned int idx0 = 0;
    for (int i = 0; i < total_arcs; ++i)
    {
        if ((i % 2) == 0)
        { // evens
            for (int j = 0; j < evens.size(); ++j)
            {
                CompGeo::XYZ P = origin + evens[j];
                //cout << "arc[" << to_string(i) << "], vtx[" << to_string(j) << "] calc'd: "; P.PrintXYZ();

                unsigned int idx = Model::v.size();
                bool idxNew = true;
                if ((j == 0) && (aseam.size() > 0)) idxNew = handleSeam(i, idx, aseam, P, "start_angle");
                if ((j == (evens.size() - 1)) && (bseam.size() > 0)) idxNew = handleSeam(i, idx, bseam, P, "end_angle");
                if ((i == 0) && (fseam.size() > 0)) idxNew = handleSeam(j, idx, fseam, P, "first_circle");
                if ((i == (total_arcs - 1)) && (lseam.size() > 0)) idxNew = handleSeam(j, idx, lseam, P, "last_circle");
                if (j == 0) idx0 = idx;
                if ((j == (evens.size() - 1)) && full_rotation && idxNew) 
                {
                    idxNew = false;
                    idx = idx0;
                }

                if (idxNew) Model::v.push_back(P);
                arcs.push_back(idx);

                //cout << "index in: " << to_string(idx) << "\n";
            }

            //cout << "\n";
        }
        else
        { // odds (cylinders begin and end w/ even arcs)
            for (int j = 0; j < odds.size(); ++j)
            {
                CompGeo::XYZ P = origin + odds[j];
                //cout << "arc[" << to_string(i) << "], vtx[" << to_string(j) << "] calc'd: "; P.PrintXYZ();

                unsigned int idx = Model::v.size();
                bool idxNew = true;
                if ((j == 0) && (aseam.size() > 0)) idxNew = handleSeam(i, idx, aseam, P, "start_angle");
                if ((j == (odds.size() - 1)) && (bseam.size() > 0)) idxNew = handleSeam(i, idx, bseam, P, "end_angle");
                if (j == 0) idx0 = idx;
                if ((j == (odds.size() - 1)) && full_rotation && idxNew) 
                {
                    idxNew = false;
                    idx = idx0;
                }

                if (idxNew) Model::v.push_back(P);
                arcs.push_back(idx);

                //cout << "index in: " << to_string(idx) << "\n";
            }

            //cout << "\n";
        }
        origin += axle_delta * axleHAT;
        //double mult_length = (i + 1) * error;
        //if (mult_length > axle_length) origin = ctr_last;
    }

}

void cylinder::calculateTriangles(void)
{
    int even_vtx = n2r == 1? N + 1: N + 2,
        odd_vtx = n2r == 1? N + 2: N + 3,
        acc = 0, triCount = 0;

    bool onFinal = false, evenFirst = true, startedTris = false;

    for (int i = 0; i < (total_arcs - 1); ++i)
    {
        int acc_delta = evenFirst? even_vtx: odd_vtx;
        for (int j = 0; j < even_vtx; ++j)
        { // uses fact that even_vtx = odd_vtx - 1
            onFinal = j == (even_vtx - 1);
            unsigned int A1 = arcs[acc + j], B1 = evenFirst && onFinal? A1: arcs[acc + j + 1],
                A2 = arcs[acc + acc_delta + j], B2 = !evenFirst && onFinal? A2: arcs[acc + acc_delta + j + 1];
            if (!startedTris)
            {
                start_triangles = Model::t.size();
                startedTris = true;
            }
            if (A1 != B1) {Model::t.push_back(A1); Model::t.push_back(n2r == 1? A2: B1); Model::t.push_back(n2r == 1? B1: A2); ++triCount;}
            if (A2 != B2) {Model::t.push_back(B1); Model::t.push_back(n2r == 1? A2: B2); Model::t.push_back(n2r == 1? B2: A2); ++triCount;}

        }

        evenFirst = !evenFirst;
        acc += acc_delta;
    }
    total_triangles = triCount;
}

void cylinder::buildSection(const pxmlnode & sect)
{
    id = "id";
    if(!sect->findAttribute(id)) handleError(string("error: no id with section in cylinder\n") + sect->write_node());
    //getSection(sect);
    //calculateOtherParameters(); called from calculateArcs()
    getSection(sect);
    calculateArcs();
    calculateTriangles();

}

//string cylinder::printXYZ(const string & label, const CompGeo::XYZ & xyz)
//{
//    return string("\t" + label + ":(" + to_string(xyz.x) + ", " + to_string(xyz.y) + ", " + to_string(xyz.z) + ")");
//}

string cylinder::printSection(void)
{
    string r = "cylinder [" + id + "]\n" +
        "\tradius:[" + to_string(radius) + "]\terror:[" + to_string(error) + "]"+
        "\talpha:[" + to_string(alpha) + "]\tbeta:[" + to_string(beta) + "]" +
        "\taxle_length:[" + to_string(axle_length) + "]\taxle_delta:[" + to_string(axle_delta) + "]\tn2r:[" + to_string(n2r) + "]" +
        "\tN:[" + to_string(N) + "]\ttotal_arcs:[" + to_string(total_arcs) + "]" +
        "\ttotal_triangles:[" + to_string(total_triangles) + "]\tstart_triangles:[" + to_string(start_triangles) + "]" +
        ctr_first.toStr("ctr_first") + ctr_last.toStr("ctr_last") + axleHAT.toStr("axleHAT") +
        NHAT.toStr("NHAT") + IHAT.toStr("IHAT") + JHAT.toStr("JHAT") + "\nARC VERTEX INDICES\n==============\n";

    int even_vtx = n2r == 1? N + 1: N + 2,
        odd_vtx = n2r == 1? N + 2: N + 3,
        acc = 0, acc_delta = even_vtx;
    
    for(int i = 0; i < total_arcs; ++i)
    {
        acc_delta = (i % 2) == 0? even_vtx: odd_vtx;
        r += "arc " + to_string(i) + "):{";
        for (int j = 0; j < acc_delta; ++j)
        {
            r += to_string(arcs[acc + j]);
            if (j < (acc_delta - 1)) r += ", ";
        }
        r += "}\n";
        acc += acc_delta;
    }
    if (seams.size() > 0)
    {
        map<char, string> seamOn = {pair('a', "start_angle"), pair('b', "end_angle"), pair('f', "first_circle"), pair('l', "last_circle")};

        r += "\nSEAMS INCORPORATED\n^^^^^^^^^^^^^^^^^^^^^^\n";
        int sCount = -1;
        for (map<char, seamType>::iterator sit = seams.begin(); sit != seams.end(); ++sit)
        {
            char sectionType = sit->first;
            seamType s = sit->second;
            r +=  "\t[" + to_string(++sCount) + "]\tfrom: " + s.from + "::" + s.part + 
                "\tapplied to this cylinder's " + seamOn[sectionType] + "\n";
        }
    }
    r += "END cylinder " + id + "\n==================\n";
   return r;
}