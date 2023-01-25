#include "first.hpp"
#include "section.hpp"
#include "Model.hpp"
#include "cylinder.hpp"
#include "plane.hpp"

using namespace std;

// static initializations:
vector<CompGeo::XYZ> Model::v = {}; // vertices
vector<unsigned int> Model::t = {}; // triangles as sets of 3 indices into v
map<std::string, psection> Model::s = {}; // key = id, value = ptr to class
Matrix Model::T = Matrix(4);  // Transform initialized to 4X4 identity matrix

/*
sectionStructType::sectionStructType(const sectionStructType & s)
{
    sType = s.sType;
    pcylinder pc = NULL;
    pplane pp = NULL;
    switch(sType)
    {
        case 'c':
            pc = static_cast<pcylinder>(s.sect);
            sect = new cylinder(*pc);
            break;
        case 'p':
            pp = static_cast<pplane>(s.sect);
            sect = new plane(*pp);
            break; 
    }
}

sectionStructType::~sectionStructType(void)
{
    switch(sType)
    {
        case 'c':
            delete static_cast<pcylinder>(sect);
            break;
        case 'p':
            delete static_cast<pplane>(sect);
            break;
    }
}
*/
void Model::writeDBG(const string & pth)
{
    ofstream os (pth, ofstream::out);

    os << "debug file for xml model \n"
       << "there are " << s.size() << " model sections\n"
       << "followed by " << v.size() << " vertices\n"
       << "followed by " << t.size() << " triangles\n"
       << "\n";

    for (map<string, psection>::iterator sit = s.begin(); sit != s.end(); ++sit)
    { // this will be in alphabetical order, not the order in the xml
        psection & pst = sit->second;
        os << pst->printSection();
        /*
        switch(pst->sType)
        {
            case 'c':
            {
                pcylinder pCyl = static_cast<pcylinder>(pst->sect);
                os << pCyl->printSection();
            }
            break;
            
            case 'p':
            {
                pplane pPl = static_cast<pplane>(pst->sect);
                os << pPl->printSection();
            }
            break;
        }
        */
    } 

    int vCount = -1;
    os << "\nVERTICES: (X, Y, Z)\n";
    for(vector<CompGeo::XYZ>::iterator vit = v.begin(); vit != v.end(); ++vit)
    {
        CompGeo::XYZ xyz = *vit;
        os << "\t[" << ++vCount << "] (" << xyz.x << ", " << xyz.y << ", " << xyz.z << ")\n";
    }

    int tCount = -1;
    os << "\nTRIANGLES:\n";
    for (vector<unsigned int>::iterator tit = t.begin(); tit != t.end(); ++tit)
    {
        unsigned int A = *tit, B = *(++tit), C = *(++tit);
        os << "\t[" << ++tCount << "] {" << A << ", " << B << ", " << C << "}\n";
    }

    os.close();
    cout << "output: " << pth << "\n";
}

void Model::writePLY(const string & pth)
{
    ofstream os (pth, ofstream::out);
    
    os << "ply\n"
       << "format ascii 1.0\n"
       << "comment xml model output\n"
       << "element vertex " << v.size() << "\n"
       << "property float x\n"
       << "property float y\n"
       << "property float z\n"
       << "element face " << t.size() << "\n"
       << "property list uchar int vertex_indices\n"
       << "end_header\n";

    for(vector<CompGeo::XYZ>::iterator vit = v.begin(); vit != v.end(); ++vit)
    {
        CompGeo::XYZ xyz = *vit;
        os << xyz.x << " " << xyz.y << " " << xyz.z << "\n";
    }

    for (vector<unsigned int>::iterator tit = t.begin(); tit != t.end(); ++tit)
    {
        unsigned int A = *tit, B = *(++tit), C = *(++tit);
        os << "3 " << A << " " << B << " " << C << "\n";
    }

    os.close();
    cout << "output: " << pth << "\n";
}

void Model::writeSTL(const string & pth)
{
    ofstream os (pth, ofstream::binary);
    char hdr[80] = "STL BINARY                                                                     ",
        b2[2] = {'\0', '\0'};
    os.write(hdr, 80);
    union b4
    {
        char write_buff[4];
        int i_buff;
        float f_buff;
    } w4;

    int numTri = t.size() / 3;
    w4.i_buff = numTri;
    os.write(w4.write_buff, 4);

    struct floatify
    {
        floatify(void):x(0.0f), y(0.0f), z(0.0f){}
        floatify(CompGeo::XYZ xyz): x(static_cast<float>(CompGeo::Rounding(xyz.x, 6))), 
            y(static_cast<float>(CompGeo::Rounding(xyz.y, 6))), 
            z(static_cast<float>(CompGeo::Rounding(xyz.z, 6))) {}
        float x, y, z;
    };

    union b48
    {
        b48(void){}
        char wBuff[48];
        floatify fIn[4];
    } w48;
    
    for (int i = 0; i < numTri; ++i)
    {
        CompGeo::XYZ A = v[t[3 * i]], B = v[t[3 * i + 1]], C = v[t[3 * i + 2]],
            AB = B - A, AC = C - A, N = CompGeo::Cross(AB, AC);
        N /= N.GetMagnitude();
        w48.fIn[0] = floatify(N);
        w48.fIn[1] = floatify(A);
        w48.fIn[2] = floatify(B);
        w48.fIn[3] = floatify(C);

        os.write(w48.wBuff, 48);
        os.write(b2, 2);
    }
    
    os.close();
    cout << "output: " << pth << "\n";
}

void Model::rotateT(const char & axis, const double & angle)
{

    double c = cos(angle), s = sin(angle);
    Matrix r = GetRotationMatrix(axis, c, s);
    T *= r;
}

void Model::translateT(CompGeo::XYZ & xyz)
{
    Matrix buff = Matrix(4, 4); // all zero matrix
    double **& a = buff.alpha;
    CompGeo::XYZPlus xp = CompGeo::XYZPlus(xyz);
    for (unsigned char i = 0; i < 3; ++i)
    {
        a[i][3] = *(xp.dims[i]); //xyz.dim(i);

    }
    T += buff;
}

void Model::transformVerts(void)
{
    //unsigned int total_chars = 0;
    //char * mElms = T.GetElements("transform matrix: ", 8, 3, total_chars);
    //string tmStr = string(mElms, total_chars);
    //cout << tmStr << "\n";

    unsigned int vCount = 0;

    double colbuff[4] = {0.0};

    for (vector<CompGeo::XYZ>::iterator vit = v.begin(); vit != v.end(); ++vit)
    {
        CompGeo::XYZ & xyz = *vit;
        xyz.x = CompGeo::Rounding(xyz.x, 8); 
        xyz.y = CompGeo::Rounding(xyz.y, 8);
        xyz.z = CompGeo::Rounding(xyz.z, 8);
        Matrix xyzM = Matrix(xyz);
        //mElms = xyzM.GetElements("xyz before matrix: ", 8, 3, total_chars); tmStr = string(mElms, total_chars); cout << tmStr << "\n";
        xyzM.ReverseTimesEquals(T, colbuff);
        //mElms = xyzM.GetElements("xyz after matrix: ", 8, 3, total_chars); tmStr = string(mElms, total_chars); cout << tmStr << "\n";
        //Matrix r = T;
        //mElms = r.GetElements("r (copy of T) matrix: ", 8, 3, total_chars); tmStr = string(mElms, total_chars); cout << tmStr << "\n";
        //r *= xyzM;
        //mElms = r.GetElements("r * xyz matrix: ", 8, 3, total_chars); tmStr = string(mElms, total_chars); cout << tmStr << "\n";
        //xyz = r.GetXYZ();
        xyz = xyzM.GetXYZ();
        ++vCount;
    }
}

void Model::finalizeModel(void)
{
    transformVerts();
}
