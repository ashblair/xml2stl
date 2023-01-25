#include "first.hpp"
#include "xmlHandler.hpp"
#include "section.hpp"
#include "Model.hpp"
#include "cylinder.hpp"
#include "plane.hpp"
#include "patch.hpp"

using namespace std;

int main (int argc, char *argv[])
{
    if (argc < 2) 
    {
        cout << "call with {/path/ + model_name}.xml on the command line to output a model .stl file on the same path!\n";
        exit(EXIT_FAILURE);
    }
    Model * m = new Model;
    string xmlPath = argv[1], fileNoExt = stringops::getNameNoExt(xmlPath), path = stringops::backOne(xmlPath);
    xmlHandler * xH = new xmlHandler;
    xH->LoadAndParseXML(xmlPath);
    m->finalizeModel();
    m->writeSTL(path + "/" + fileNoExt + ".stl");
    //m->writePLY(path + "/" + fileNoExt + ".ply");
    //m->writeDBG(path + "/" + fileNoExt + ".dbg");
    delete xH;
    delete m;
    exit(EXIT_SUCCESS);
}
