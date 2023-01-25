#include "first.hpp"
#include "xmlHandler.hpp"
#include "section.hpp"

using namespace std;

void section::handleError(const string & err)
{
    cout << err << "\n";
    exit(EXIT_FAILURE);
}

CompGeo::XYZ section::getNamedVector(const string & vect)
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
