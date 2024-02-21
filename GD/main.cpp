//.......................BiB............................................
#include <DGtal/base/Common.h>                                       //.
#include <DGtal/helpers/StdDefs.h>                                   //.
#include <DGtal/images/ImageSelector.h>                              //.
#include "DGtal/io/readers/PGMReader.h"                              //.
#include "DGtal/io/writers/GenericWriter.h"                          //.
#include <DGtal/images/imagesSetsUtils/SetFromImage.h>               //.
#include <DGtal/io/boards/Board2D.h>                                 //.
#include <DGtal/io/colormaps/ColorBrightnessColorMap.h>              //.
#include <DGtal/topology/SurfelAdjacency.h>                          //.
#include <DGtal/topology/helpers/Surfaces.h>                         //.
#include "DGtal/io/Color.h"                                          //.
#include "DGtal/io/readers/GenericReader.h"                          //.
#include <DGtal/helpers/StdDefs.h>                                   //.
#include "DGtal/helpers/StdDefs.h"                                   //.
#include </usr/include/cairo/cairo.h>                                //.
#include "DGtal/topology/KhalimskySpaceND.h"                         //.
#include <DGtal/topology/DigitalTopology.h>                          //.
#include"DGtal/io/colormaps/GradientColorMap.h"                      //.
#include "DGtal/geometry/curves/ArithmeticalDSSComputer.h"           //.
#include "DGtal/geometry/curves/GreedySegmentation.h"                //.
#include "DGtal/geometry/curves/estimation/DSSLengthEstimator.h"     //.
#include "DGtal/geometry/curves/FreemanChain.h"                      //.
#include <sstream>                                                   //.
//......................................................................
typedef DGtal::ImageContainerBySTLMap<DGtal::Z3i::Domain, unsigned int> Image3D;
typedef DGtal::ImageContainerBySTLMap<DGtal::Z2i::Domain, unsigned int> Image2D;
using namespace std;
using namespace DGtal;
using namespace Z2i;
double calculeAreaShoelace(const std::vector<LibBoard::Point>& Vertices) {//Calculer Are avec Shoelace
    double Area = 0.0;
    for (int i = 0; i < Vertices.size(); ++i) {
        const auto& c1 = Vertices[i];
        const auto& c2 = Vertices[(i + 1) % Vertices.size()];
        Area += (c1.x * c2.y) - (c2.x * c1.y);
    }
    return std::abs(Area) / 2.0;
}
double calculedistance(const LibBoard::Point& p1, const LibBoard::Point& p2) {  //calculer la distance entre deux points
    double X= p1.x - p2.x;
    double Y = p1.y - p2.y;
    return std::sqrt(X * X + Y * Y);
}
double calculePerimetre(const std::vector<LibBoard::Point>& Vertices) {  //calculer perimetre
    double Perimetre = 0.0;
    for (int i = 0; i < Vertices.size(); ++i) {
        const auto& c1 = Vertices[i];
        const auto& c2 = Vertices[(i + 1) % Vertices.size()];
        Perimetre += calculedistance(c1, c2);
    }
    return Perimetre;
}
template<class T>
std::vector<Z2i::SCell> getBoundary(T & object)
{
    //Khalimsky space
    typedef KhalimskySpaceND< 2, int > KSpace;
    KSpace K;
    // we need to add a margine to prevent situations such that an object touch the bourder of the domain
    K.init( object.domain().lowerBound() - Point(1,1),
                   object.domain().upperBound() + Point(1,1), true );
    
    // 1) Call Surfaces::findABel() to find a cell which belongs to the border
    Z2i::SCell cell = Surfaces<Z2i::KSpace>::findABel(K, object.pointSet(),10000); //documentation DGtal
    std::vector<Z2i::SCell> boundaryPoints; // boundary points are going to be stored here
    // 2) Call Surfece::track2DBoundaryPoints to extract the boundary of the object
    SurfelAdjacency<2> SAdj( true );
    Surfaces<Z2i::KSpace>::track2DBoundary( boundaryPoints, K, SAdj, object.pointSet(), cell );//documentation DGtal
    return boundaryPoints;
}
template<class T>
void sendToBoard(int i,Board2D & board, T & p_Object, DGtal::Color p_Color) {
    board << CustomStyle(p_Object.className(), new DGtal::CustomFillColor(p_Color));
    board << p_Object;
    std::ostringstream fileName;
    fileName << "Etape2 Digital Object image(" << i << ").eps";
    board.saveEPS(fileName.str().c_str());
}
std::string computeFCC(const Z2i::Curve::PointsRange& Range) {    //calculer freemanchain 
    std::string fcc;
    Z2i::Point beginPoint = *Range.begin();
    Z2i::Point prevPoint = *Range.begin();
    Z2i::Point LastD = beginPoint - *Range.end();
    std::vector<Z2i::Point> ListD = {{1, 0}, {0, 1}, {-1, 0}, {0, -1}};
    for (Z2i::Point R : Range) {
        Z2i::Point displacement = R - prevPoint;
        int direction = -1;       
        for (int i = 0; i < ListD.size(); ++i) {
            if (ListD[i] == displacement) {
                direction = i;
                break;
            }
        }
        if (direction != -1) {
            fcc+= std::to_string(direction);
        }
        prevPoint = R;
    }
    
    int direction = -1;
    for (int i = 0; i < ListD.size(); ++i) {
        if (ListD[i] == LastD) {
            direction = i;
            break;
        }
    }
    if (direction != -1) {
        fcc += std::to_string(direction);
    }
    return fcc;
}
typedef FreemanChain<int> Contour4;
typedef ArithmeticalDSSComputer<Contour4::ConstIterator, int, 4> DSS4;
typedef GreedySegmentation<DSS4> Decomposition4;
// Créer un objet FreemanChain (Contour4) à partir de freemanChainCode
Contour4 ContourDSS(const std::string& freemanChainCode,int numberOfImage) { 
    std::stringstream ss(std::stringstream::in | std::stringstream::out);
    ss << "10000 10000 " + freemanChainCode << std::endl;
    Contour4 theContour(ss);
    return theContour;
}
int main(int argc, char** argv)
{
    setlocale(LC_NUMERIC, "us_US"); //To prevent French local settings
    typedef ImageSelector<Domain, unsigned char >::Type Image; // type of image
    typedef DigitalSetSelector< Domain, BIG_DS+HIGH_BEL_DS >::Type DigitalSet; // Digital set type
    typedef Object<DT4_8, DigitalSet> ObjectType; 
//==============================Etape 2 =============================================
    std::vector<std::string> pathsImages = {
        "./RiceGrains/Rice_camargue_seg_bin.pgm",
        "./RiceGrains/Rice_japonais_seg_bin.pgm",
        "./RiceGrains/Rice_basmati_seg_bin.pgm"
                                            };
    int IndiceImage=1;
    for (const std::string &Patheimage : pathsImages) {
        // read an image
        //Image image = PGMReader<Image>::importPGM ("../RiceGrains/Rice_japonais_seg_bin.pgm"); // you have to provide a correct path as a parameter
    Image image = PGMReader<Image>::importPGM (Patheimage); 
    // 1) make a "digital set" of proper size
    Z2i::DigitalSet digital_set (image.domain());
    // 2) populate a digital set from the image using SetFromImage::append()
    SetFromImage<Z2i::DigitalSet>::append<Image>(digital_set, image, 1, 255);
    std::cout << "======================================================================="<< std::endl;
    std::cout << "|                                                                     |"<< std::endl;
    std::cout << "===================== Image  "<< IndiceImage << " ========================================" << endl;
    std::cout << "Path= "<< Patheimage<< endl;
    Board2D aBoard;
    aBoard << digital_set;
    std::ostringstream fileName;
    fileName <<"Etape 2 Digital Set image (" << IndiceImage << ").eps";
    aBoard.saveEPS(fileName.str().c_str());
    /////////////////////////////////// Etape 3 /////////////////////////////////////////////////
    typedef SpaceND< 2,int > Z2;
    typedef MetricAdjacency< Z2, 1 > Adj4;
    typedef MetricAdjacency< Z2, 2 > Adj8;
    typedef DigitalTopology< Adj4, Adj8 > DT4_8;
    Adj4 adj4;
    Adj8 adj8;
    DT4_8 dt4_8( adj4, adj8);
    // 3) Create a digital object from the digital set
    std::vector< ObjectType > objects;          // All connected components are going to be stored in it
    std::back_insert_iterator< std::vector< ObjectType > > inserter( objects ); // Iterator used to populated "objects".
   
    ObjectType objj( dt4_8, digital_set);
    // 4) Set the adjacency pair and obtain the connected components using "writeComponents"
    objj.writeComponents( inserter );
    std::cout << " number of components image  "<< IndiceImage << " = " << objects.size() << endl; // Right now size of "objects" is the number of conected components
    //=================Elimination des grains situées en dehors des bords=========================
    std::vector<int> ObjectsEliminate;
    const auto OutBounds = [&image](const Z2i::Point& point) {
        const auto imageSize = image.domain().upperBound();
        return point[0] <= 0 || point[0] >= imageSize[0] || point[1] <= 0 || point[1] >= imageSize[1];
    };
    for (int i = 0; i < objects.size(); ++i) {
        const ObjectType& obj = objects[i];
        if (std::any_of(obj.pointSet().begin(), obj.pointSet().end(), OutBounds)) {
            ObjectsEliminate.push_back(i);
        }
    }
    for (auto it = ObjectsEliminate.rbegin(); it != ObjectsEliminate.rend(); ++it) {
    int index = *it;
    if (index >= 0 && index < objects.size()) {
        objects.erase(objects.begin() + index);
        }
    }
    std::cout << "number of components after elimination image "<< IndiceImage << " = " << objects.size() << endl;
   //===============================================================================
    Board2D aaBoard;
    Board2D aaaboard;
    //Board2D board;
    //board.saveEPS(filename.str().c_str());
    sendToBoard( IndiceImage,aaBoard, objects[45], Color::Red);   
    std::vector<vector<Z2i::SCell>> boundryVector;
    for (int i = 0; i < objects.size(); i++)//calculer les bords des grains dans l'image
    {
       boundryVector.push_back(getBoundary(objects[i]));
    }   
    std::vector<Z2i::Curve > curveVector ;
    for (int i = 0; i < boundryVector.size(); i++)
    {
        Z2i::Curve curve;
        curve.initFromSCellsVector(boundryVector.at(i));
        curveVector.push_back(curve);
    } 
    for (int i = 0; i < boundryVector.at(45).size(); i++)
    {
       aaaboard<< CustomStyle((boundryVector.at(45).at(i)).className() ,new DGtal::CustomFillColor( Color::Black))<< boundryVector.at(45).at(i);
    }
    std::ostringstream fileNamee ;
    fileNamee << "Etape 3 (Digital Object Boundary)image " << IndiceImage << ".eps";
    aaaboard.saveEPS(fileNamee.str().c_str());
    std::vector<std::vector<LibBoard::Point >> polygonVertices;
     for (int i = 0; i < curveVector.size(); i++)// ploygonisation (Greedy segmentation)
    {
        typedef Z2i::Curve::PointsRange Range;
        Range range = curveVector.at(i).getPointsRange();
        std::string fcc = computeFCC(range);
        Contour4 theContour = ContourDSS(fcc,IndiceImage);
        Decomposition4 theDecomposition(theContour.begin(), theContour.end(), DSS4());
        std::vector<LibBoard::Point> vertices;
        LibBoard::Point vertex0(0, 0);
    for (Decomposition4::SegmentComputerIterator it = theDecomposition.begin(), itEnd = theDecomposition.end(); it != itEnd; ++it) {
        if (vertex0.x == 0) {
            vertex0.x = (*it).Uf()[0];
            vertex0.y = (*it).Uf()[1];
        }
        vertices.push_back(LibBoard::Point((*it).Uf()[0], (*it).Uf()[1]));
    }
        vertices.push_back(LibBoard::Point(vertex0.x, vertex0.y));
        polygonVertices.push_back(vertices);
    }
    Board2D aaaaBoard;
    aaaaBoard.setPenColor(Color::Blue);
    aaaaBoard.drawPolyline(polygonVertices.at(45));
    std::ostringstream fileName2;
    fileName2 << "Etape 4 (Polygonisation) image " << IndiceImage << ".eps";
    aaaaBoard.saveEPS(fileName2.str().c_str());
    
    // //===============Etape 5===========================
    //     //===============methode 1=======================================
    //     std::cout << "Areas (methode 1)image  "<< IndiceImage<< std::endl;
    //     for (int i=0;i<objects.size();i++) {
    //     DigitalSet objjj = objects.at(i).pointSet();
    //     std::cout << "Area (methode 1) Component " << i << "=" << objjj.size() << std::endl;
    //     //===============methode 2========================================
    //     double polygonArea = calculeAreaShoelace(polygonVertices.at(i));
    //     std::cout << "Area (methode 2) Component " << i << " = " << polygonArea << std::endl;
    //     }
        
    // //========================Etape 6===================================================
    //     //======methode 1 ==========================================
    //     std::cout << "Perimetres (methode 1/2 image : "<< IndiceImage<< std::endl;
    //     for (int i=0;i<polygonVertices.size();i++) {
    //     std::cout << "Perimetere Component "<< i <<"(méthode 1) = "  << boundryVector.at(i).size() << std::endl;
    //     //======methode 2 ==========================================
    //     double PerimetrePloyg = calculePerimetre(polygonVertices.at(i));
    //     std::cout << "Perimetere Component "<< i<< "(méthode 2) = " << PerimetrePloyg << std::endl;
    //     }

    // //=========================Etape 7==================================================
    //     std::cout << "CIRCULARITY Components "<< std::endl;
    //     for (int i=0;i<polygonVertices.size();i++) {
    //     double PerimetrePloyg = calculePerimetre(polygonVertices.at(i));
    //     double polygonArea = calculeAreaShoelace(polygonVertices.at(i));
    //     double CIR = 4*3.141*polygonArea/(PerimetrePloyg*PerimetrePloyg);
    //     std::cout << "CIRCULARITY Component "<< i << " = " << CIR << std::endl;
    // }
      
    IndiceImage++;
    }
        
  return 0;
}
//g++ -o main main.cpp -lDGtal
// ./main