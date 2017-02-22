//  ===================================================================
//
//  Program:   Convex Quadrilateral Mesh Generator (CQMesh)
//  Module:    compquad.h
//  Language:  C++
//  Date:      December 2002
//  Version:   1.0.0
//
//  Copyright (c) 2002, 2003
//
//  Marcelo Siqueira
//
//  Graduate Student
//  Department of Computer and Information Science
//  Moore School Building
//  200 South 33rd Street, Suite #62
//  Philadelphia, PA 19104-6389
//  Phone No.: (215) 898-5879
//
//  marcelos@seas.upenn.edu
//
// This program  is distributed  WITHOUT ANY WARRANTY,  and it  may be
// freely redistributed under the condition that the copyright notices
// are  not  removed,  and  no  compensation  is  received.   Private,
// research, and institutional use  is free. Distribution of this code
// as  part of  a  commercial  system is  permissible  ONLY BY  DIRECT
// ARRANGEMENT WITH THE AUTHOR.
//
//  ==================================================================

#ifndef __compquad_h
#define __compquad_h

#include <list>
#include <map>

#include <cmath>

#include "typedefs.h"
#include "err.h"

#include "vertex2.h"
#include "halfedge2.h"
#include "edge2.h"
#include "face2.h"
#include "complex2.h"
#include "sptree.h"
#include "polygon2.h"
#include "quad2.h"

const int MAX_EDGE_DISP = 3;
const double TOL = 1e-10;

typedef std::map<tCQMVertex2*,int> QuadVertexSet;

typedef enum {R1E24E34, R4E13E12, R1E24E23, R4E13E23, R1R2E24E34,
              R3R4E13E12, R1R2E24E23, R3R4E13E23, R1R4E24E34,
              R1R4E13E12, R1R4E24E23, R1R4E13E23, R1E13E23, R4E24E23,
              R1E13E12, R4E24E34} PentCase;

typedef enum {CCCCCC, RCCCCC, RCRCCC, RRCCCC, RCCRCC, RCRCRC, RRCRCC,
	      RRCCRC, RRRCCC} HexCase;

#define IS_EQUAL(a,b) ((fabs(a-b) < TOL))
#define IS_ZERO(x) ((fabs(x) < 1e-24))

// -------------------------------------------------------------------
// tCQMGeneralControl
// ------------------
class tCQMCompQuad {
   private:
     // pointer  to a  triangular mesh  that is  supposed to  be built
     // alreay.
     tCQMComplex2* m_mesh;

     // vertex  counter  used   to  generate  unique  identifiers  for
     // quadrilaterals vertices.
     int m_counter;

     // counter of Steiner points used by the quadrangulation
     int m_spoints;

     //
     // private methods
     //

     // inline
     int getNewVertexIndex();

     // not inline
     void computeRootVertices(std::list<tCQMHalfEdge2*>&) const;

     bool isPointOnLine(const Point2&,const Point2&,
       const Point2&);
     bool isSegInLine(const Point2&,const Point2&,
       const Point2&,const Point2&);
     bool segLineIntersection(const Point2&,const Point2&,
       const Point2&,const Point2& v2,Point2&);
     tCQMSpanningTreeVertex* getChildOfDegreeNParent(int, 
       const std::list<tCQMSpanningTreeVertex*>&);
     tCQMSpanningTree* getSpanningTree(tCQMDualGraphVertex*);
     tCQMIndVertex2 getQuadVertex(QuadVertexSet&, tCQMVertex2*);
     tCQMIndVertex2 createSteinerPoint(const Point2&);
     void createQuad(const tCQMIndVertex2&, const tCQMIndVertex2&,
       const tCQMIndVertex2&, const tCQMIndVertex2&, bool, bool,
       bool, bool, std::list<tCQMQuadrilateral2*>*)
       throw (tExceptionObject);
     void insertVertexInContour(tCQMSpanningTreeVertex*, 
       QuadVertexSet&, std::list<tCQMQuadrilateral2*>*);
     void halfSpaceIntersection(const std::list<Point2>&,
       const Point2&, const Point2&, std::list<Point2>&,
       CGAL::Oriented_side orient) throw (tExceptionObject);
     bool myAssign(Point2&, const Segment2&, const Line2&);
     tCQMSpanningTreeVertex* findOneQuadrilateral(
       const std::list<tCQMSpanningTreeVertex*>&);
     tCQMSpanningTreeVertex* findOnePentagon(
       const std::list<tCQMSpanningTreeVertex*>&);
     tCQMPolygon2* getPentagon(tCQMSpanningTreeVertex*,
       tCQMSpanningTreeVertex*, int&, QuadVertexSet&);
     tCQMPolygon2* getPentagon(tCQMSpanningTreeVertex*,
       tCQMSpanningTreeVertex*, tCQMSpanningTreeVertex*,
       QuadVertexSet&);
     tCQMPolygon2* getHexagon(tCQMSpanningTreeVertex*,
       tCQMSpanningTreeVertex*, tCQMSpanningTreeVertex*,
       QuadVertexSet&);
     tCQMPolygon2* getHexagon(tCQMSpanningTreeVertex*,
       tCQMSpanningTreeVertex*, tCQMSpanningTreeVertex*,
       tCQMSpanningTreeVertex*, QuadVertexSet&);
     PentCase classifyPentagon(tCQMPolygon2* p,
       tCQMSpanningTreeVertex*, tCQMSpanningTreeVertex*,
       tCQMSpanningTreeVertex*);
     HexCase classifyHexagon(tCQMPolygon2*);
     bool makeQuadIfConvex(tCQMSpanningTreeVertex*,
       tCQMSpanningTreeVertex*, QuadVertexSet& vs,
       std::list<tCQMQuadrilateral2*>*);
     tCQMIndVertex2 quadrangulateDegTriangle(const tCQMIndVertex2&,
       const tCQMIndVertex2&, const tCQMIndVertex2&, const tCQMIndVertex2&,
       bool, bool, std::list<tCQMQuadrilateral2*>*);
     tCQMIndVertex2 quadrangulateDegPentagon1(const tCQMIndVertex2&,
       const tCQMIndVertex2&, const tCQMIndVertex2&, const tCQMIndVertex2&,
       const tCQMIndVertex2&, bool, bool, bool, bool, bool,
       std::list<tCQMQuadrilateral2*>*);
     void quadrangulatePentR1E24E34(tCQMPolygon2*, 
       tCQMSpanningTreeVertex*, QuadVertexSet&, std::list<tCQMQuadrilateral2*>*);
     void quadrangulatePentR4E13E12(tCQMPolygon2*,
       tCQMSpanningTreeVertex*, QuadVertexSet&, std::list<tCQMQuadrilateral2*>*);
     void quadrangulatePentR1E24E23(tCQMPolygon2*,
       tCQMSpanningTreeVertex*, QuadVertexSet&, std::list<tCQMQuadrilateral2*>*);
     void quadrangulatePentR4E13E23(tCQMPolygon2*,
       tCQMSpanningTreeVertex*, QuadVertexSet&, std::list<tCQMQuadrilateral2*>*);
     void quadrangulatePentR1R2E24E34(tCQMPolygon2*,
       tCQMSpanningTreeVertex*, QuadVertexSet&, std::list<tCQMQuadrilateral2*>*)
       throw (tExceptionObject);
     void quadrangulatePentR3R4E13E12(tCQMPolygon2*,
       tCQMSpanningTreeVertex*, QuadVertexSet&, std::list<tCQMQuadrilateral2*>*)
       throw (tExceptionObject);
     void quadrangulatePentR1R2E24E23(tCQMPolygon2*,
       tCQMSpanningTreeVertex*, QuadVertexSet&, std::list<tCQMQuadrilateral2*>*)
       throw (tExceptionObject);
     void quadrangulatePentR3R4E13E23(tCQMPolygon2*,
       tCQMSpanningTreeVertex*, QuadVertexSet&, std::list<tCQMQuadrilateral2*>*)
       throw (tExceptionObject);
     void quadrangulatePentR1R4E24E34(tCQMPolygon2*,
       tCQMSpanningTreeVertex*, QuadVertexSet&, std::list<tCQMQuadrilateral2*>*);
     void quadrangulatePentR1R4E13E12(tCQMPolygon2*,
       tCQMSpanningTreeVertex*, QuadVertexSet&, std::list<tCQMQuadrilateral2*>*);
     void quadrangulatePentR1R4E24E23(tCQMPolygon2*,
       tCQMSpanningTreeVertex*, QuadVertexSet&, std::list<tCQMQuadrilateral2*>*);
     void quadrangulatePentR1R4E13E23(tCQMPolygon2*,
       tCQMSpanningTreeVertex*, QuadVertexSet&, std::list<tCQMQuadrilateral2*>*);
     void quadrangulatePentR1E13E23(tCQMPolygon2*,
       tCQMSpanningTreeVertex*, QuadVertexSet&, std::list<tCQMQuadrilateral2*>*);
     void quadrangulatePentR4E24E23(tCQMPolygon2*,
       tCQMSpanningTreeVertex*, QuadVertexSet&, std::list<tCQMQuadrilateral2*>*);
     void quadrangulatePentR1E13E12(tCQMPolygon2*,
       tCQMSpanningTreeVertex*, QuadVertexSet&, std::list<tCQMQuadrilateral2*>*);
     void quadrangulatePentR4E24E34(tCQMPolygon2*,
       tCQMSpanningTreeVertex*, QuadVertexSet&, std::list<tCQMQuadrilateral2*>*);
     void quadrangulateHexCCCCCC(tCQMPolygon2*,
       std::list<tCQMQuadrilateral2*>*);
     void quadrangulateHexRCCCCC1(const tCQMIndVertex2&,
       const tCQMIndVertex2&, const tCQMIndVertex2&,
       const tCQMIndVertex2&, const tCQMIndVertex2&,
       const tCQMIndVertex2&, bool, bool, bool, bool,
       bool, bool, std::list<tCQMQuadrilateral2*>*);
     void quadrangulateHexRCCCCC2(const tCQMIndVertex2&,
       const tCQMIndVertex2&, const tCQMIndVertex2&,
       const tCQMIndVertex2&, const tCQMIndVertex2&,
       const tCQMIndVertex2&, bool, bool, bool, bool,
       bool, bool, std::list<tCQMQuadrilateral2*>*) throw (tExceptionObject);
     void quadrangulateHexRCCCCC(tCQMPolygon2*,
       std::list<tCQMQuadrilateral2*>*);
     void quadrangulateHexRCRCCC1(const tCQMIndVertex2&,
       const tCQMIndVertex2&, const tCQMIndVertex2&,
       const tCQMIndVertex2&, const tCQMIndVertex2&,
       const tCQMIndVertex2&, bool, bool, bool, bool,
       bool, bool, std::list<tCQMQuadrilateral2*>*);
     void quadrangulateHexRCRCCC3(const tCQMIndVertex2&,
       const tCQMIndVertex2&, const tCQMIndVertex2&,
       const tCQMIndVertex2&, const tCQMIndVertex2&,
       const tCQMIndVertex2&, bool, bool, bool, bool,
       bool, bool, std::list<tCQMQuadrilateral2*>*);
     void quadrangulateHexRCRCCC(tCQMPolygon2*,
       std::list<tCQMQuadrilateral2*>*);
     void quadrangulateHexRRCCCC(tCQMPolygon2*,
       std::list<tCQMQuadrilateral2*>*);
     void quadrangulateHexRCCRCC(tCQMPolygon2*,
       std::list<tCQMQuadrilateral2*>*);
     void quadrangulateHexRCRCRC3(const tCQMIndVertex2&,
       const tCQMIndVertex2&, const tCQMIndVertex2&,
       const tCQMIndVertex2&, const tCQMIndVertex2&,
       const tCQMIndVertex2&, bool, bool, bool, bool,
       bool, bool, std::list<tCQMQuadrilateral2*>*);
     void quadrangulateHexRCRCRC(tCQMPolygon2*,
       std::list<tCQMQuadrilateral2*>*);
     void quadrangulateHexRRCRCC2(const tCQMIndVertex2&,
       const tCQMIndVertex2&, const tCQMIndVertex2&,
       const tCQMIndVertex2&, const tCQMIndVertex2&,
       const tCQMIndVertex2&, bool, bool, bool, bool,
       bool, bool, std::list<tCQMQuadrilateral2*>*);
     void quadrangulateHexRRCRCC3(const tCQMIndVertex2&,
       const tCQMIndVertex2&, const tCQMIndVertex2&,
       const tCQMIndVertex2&, const tCQMIndVertex2&,
       const tCQMIndVertex2&, bool, bool, bool, bool,
       bool, bool, std::list<tCQMQuadrilateral2*>*);
     void quadrangulateHexRRCRCC(tCQMPolygon2*,
       std::list<tCQMQuadrilateral2*>*);
     void quadrangulateHexRRCCRC2(const tCQMIndVertex2&,
       const tCQMIndVertex2&, const tCQMIndVertex2&,
       const tCQMIndVertex2&, const tCQMIndVertex2&,
       const tCQMIndVertex2&, bool, bool, bool, bool,
       bool, bool, std::list<tCQMQuadrilateral2*>*);
     void quadrangulateHexRRCCRC3(const tCQMIndVertex2&,
       const tCQMIndVertex2&, const tCQMIndVertex2&,
       const tCQMIndVertex2&, const tCQMIndVertex2&,
       const tCQMIndVertex2&, bool, bool, bool, bool,
       bool, bool, std::list<tCQMQuadrilateral2*>*);
     void quadrangulateHexRRCCRC(tCQMPolygon2*,
       std::list<tCQMQuadrilateral2*>*);
     void quadrangulateHexRRRCCC(tCQMPolygon2*,
       std::list<tCQMQuadrilateral2*>*);
     tCQMIndVertex2 quadrangulateSep1(const tCQMIndVertex2&,
       const tCQMIndVertex2&, const tCQMIndVertex2&,
       const tCQMIndVertex2&, const tCQMIndVertex2&,
       const tCQMIndVertex2&, const tCQMIndVertex2&,
       bool, bool, bool, bool, bool, bool, bool,
       std::list<tCQMQuadrilateral2*>*);
     tCQMIndVertex2 quadrangulateSep2(const tCQMIndVertex2&,
       const tCQMIndVertex2&, const tCQMIndVertex2&,
       const tCQMIndVertex2&, const tCQMIndVertex2&,
       const tCQMIndVertex2&, const tCQMIndVertex2&,
       bool, bool, bool, bool, bool, bool, bool,
       std::list<tCQMQuadrilateral2*>*);
     tCQMIndVertex2 quadrangulateSep3(const tCQMIndVertex2&,
       const tCQMIndVertex2&, const tCQMIndVertex2&,
       const tCQMIndVertex2&, const tCQMIndVertex2&,
       const tCQMIndVertex2&, const tCQMIndVertex2&,
       bool, bool, bool, bool, bool, bool, bool,
       std::list<tCQMQuadrilateral2*>*);
     tCQMIndVertex2 quadrangulateSep4(const tCQMIndVertex2&,
       const tCQMIndVertex2&, const tCQMIndVertex2&,
       const tCQMIndVertex2&, const tCQMIndVertex2&,
       const tCQMIndVertex2&, const tCQMIndVertex2&,
       bool, bool, bool, bool, bool, bool, bool,
       std::list<tCQMQuadrilateral2*>*);
     void makeQuad(tCQMSpanningTreeVertex*, QuadVertexSet&,
       std::list<tCQMQuadrilateral2*>*);
     void makeQuad(tCQMSpanningTreeVertex*, tCQMSpanningTreeVertex*,
       QuadVertexSet&, std::list<tCQMQuadrilateral2*>*);
     void makeQuad(tCQMSpanningTreeVertex*, tCQMSpanningTreeVertex*,
       tCQMSpanningTreeVertex*, QuadVertexSet&,
       std::list<tCQMQuadrilateral2*>*) throw (tExceptionObject);
     void makeQuadFromDegeneratePoly(tCQMSpanningTreeVertex*, 
       tCQMSpanningTreeVertex*, tCQMSpanningTreeVertex*,
       tCQMSpanningTreeVertex*, tCQMSpanningTreeVertex*,
       QuadVertexSet&, std::list<tCQMQuadrilateral2*>*);
     void makeQuadFromDegeneratePoly(tCQMSpanningTreeVertex*, 
       tCQMSpanningTreeVertex*, tCQMSpanningTreeVertex*,
       tCQMSpanningTreeVertex*, tCQMSpanningTreeVertex*,
       tCQMSpanningTreeVertex*, tCQMSpanningTreeVertex*,
       QuadVertexSet&, std::list<tCQMQuadrilateral2*>*);
     void makeQuad(tCQMSpanningTreeVertex*, tCQMSpanningTreeVertex*,
       tCQMSpanningTreeVertex*, tCQMSpanningTreeVertex*,
       tCQMSpanningTreeVertex*, QuadVertexSet&,
       std::list<tCQMQuadrilateral2*>*);
     void makeQuad(tCQMSpanningTreeVertex*, tCQMSpanningTreeVertex*,
       tCQMSpanningTreeVertex*, tCQMSpanningTreeVertex*,
       tCQMSpanningTreeVertex*, tCQMSpanningTreeVertex*,
       tCQMSpanningTreeVertex*, QuadVertexSet&,
       std::list<tCQMQuadrilateral2*>*);
     void makeConvexQuadFrom2Tri(tCQMSpanningTreeVertex*,
       tCQMSpanningTreeVertex*,QuadVertexSet&,
       std::list<tCQMQuadrilateral2*>*) throw (tExceptionObject);
     void makeConvexQuadFromQuad(tCQMPolygon2*,
       std::list<tCQMQuadrilateral2*>*);
     bool makeConvexQuadFromPent(tCQMSpanningTreeVertex*,
       QuadVertexSet&, std::list<tCQMQuadrilateral2*>*);
     bool makeConvexQuadFromPent(tCQMSpanningTreeVertex*,
       tCQMSpanningTreeVertex*, QuadVertexSet&,
       std::list<tCQMQuadrilateral2*>*);
     bool makeConvexQuadFromPent(tCQMSpanningTreeVertex*,
       tCQMSpanningTreeVertex*, tCQMSpanningTreeVertex*,
       QuadVertexSet&, std::list<tCQMQuadrilateral2*>*);
     bool makeConvexQuadFromHex(tCQMSpanningTreeVertex*,
       tCQMSpanningTreeVertex*, tCQMSpanningTreeVertex*,
       QuadVertexSet&, std::list<tCQMQuadrilateral2*>*);
     bool makeConvexQuadFromHex(tCQMSpanningTreeVertex*,
       tCQMSpanningTreeVertex*, tCQMSpanningTreeVertex*,
       tCQMSpanningTreeVertex*, QuadVertexSet&,
       std::list<tCQMQuadrilateral2*>*);
     bool makeConvexQuadFromHex(tCQMPolygon2*,
       QuadVertexSet&, std::list<tCQMQuadrilateral2*>*);
     void makeConvexQuadFromDegSep(tCQMSpanningTreeVertex*, 
       tCQMSpanningTreeVertex*, tCQMSpanningTreeVertex*,
       tCQMSpanningTreeVertex*, tCQMSpanningTreeVertex*,
       QuadVertexSet&, std::list<tCQMQuadrilateral2*>*);
     bool makeConvexQuadFromSep(tCQMSpanningTreeVertex*,
       tCQMSpanningTreeVertex*, tCQMSpanningTreeVertex*,
       tCQMSpanningTreeVertex*, tCQMSpanningTreeVertex*,
       QuadVertexSet&, std::list<tCQMQuadrilateral2*>*);

   public:
     // constructors
     tCQMCompQuad(tCQMComplex2*);

     // destructor
     ~tCQMCompQuad();

     // get method
     int getNumOfSteinerPoints();

     // modifiers
     std::list<tCQMQuadrilateral2*>* 
       computeConvexQuadrangulation(string fname);
};

// -------------------------------------------------------------------
// tCQMGeneralControl - inline implementation
// ------------------

// private method 
inline int
tCQMCompQuad::getNewVertexIndex()
{
  int x = 0;
  if (m_counter == 1141) {
    x += 1;
  }
  if (m_counter == 1142) {
    x += 2;
  }
  if (m_counter == 1143) {
    x += 3;
  }
  
  ++m_counter;
  return m_counter;
}

inline int
tCQMCompQuad::getNumOfSteinerPoints()
{
  return m_spoints;
}

#endif   // __gencontrol_h
