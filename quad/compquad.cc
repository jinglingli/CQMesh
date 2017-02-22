//  ===================================================================
//
//  Program:   Convex Quadrilateral Mesh Generator (CQMesh)
//  Module:    compquad.cc
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

#include "compquad.h"

#include "regiongraph.h"

#include <gmp.h>     // for all GMP and GMP++ functions
#include <gmpxx.h>
#include <execinfo.h>
#include <signal.h>

const char* const PentCaseStr[] = { "R1E24E34", "R4E13E12", "R1E24E23", "R4E13E23", "R1R2E24E34",
    "R3R4E13E12", "R1R2E24E23", "R3R4E13E23", "R1R4E24E34",
    "R1R4E13E12", "R1R4E24E23", "R1R4E13E23", "R1E13E23", "R4E24E23",
    "R1E13E12", "R4E24E34"};

int PentCaseCount[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

string spt_filename;

/*
void handler(int sig) {
    void *array[10];
    size_t size;
    
    // get void*'s for all entries on the stack
    size = backtrace(array, 10);
    
    // print out all the frames to stderr
    fprintf(stderr, "Error: signal %d:\n", sig);
    backtrace_symbols_fd(array, size, STDERR_FILENO);
    exit(1);
}
*/

// -------------------------------------------------------------------
// tCQMCompQuad - implementation
// ------------

// constructor
tCQMCompQuad::tCQMCompQuad(tCQMComplex2* mesh)
{
    m_mesh = mesh;
    m_counter = 0;
    m_spoints = 0;
}

// destructor
tCQMCompQuad::~tCQMCompQuad()
{}

//
// private methods
//

// -------------------------------------------------------------------
// Method getSpanningTree()
//
// It creates a spanning tree for  the dual graph of a mesh. The input
// is a vertex  of the dual graph,  which is chosen to be  the root of
// the spanning tree.
// -------------------------------------------------------------------
tCQMSpanningTree*
tCQMCompQuad::getSpanningTree(tCQMDualGraphVertex* v)
{
    // create a spanning tree vertex corresponding to v
    tCQMSpanningTreeVertex* par = (tCQMSpanningTreeVertex*)
    new tCQMSpanningTreeVertex(v);
    
    // create a spanning tree with root par
    tCQMSpanningTree* sp = (tCQMSpanningTree*) new tCQMSpanningTree(par);
    
    // create a working list of spanning tree vertices
    std::list<tCQMSpanningTreeVertex*> vlist;
    
    // insert the root vertex into vlist
    vlist.push_back(par);
    
    // mark vertex to indicate that it has been taken by
    // the spanning tree.
    v->setMark(true);
    
    // insert all vertices reachable from v into vlist
    while (!vlist.empty()) {
        // get first vertex in the working list of vertices
        // of the current spanning tree.
        par = vlist.front();
        
        // remove "par" from the working list of vertices
        vlist.pop_front();
        
        // get its corresponding vertex in the dual graph
        tCQMDualGraphVertex* vpar = par->getVertex();
        
        // get degree of "v" in the dual graph
        int dg = vpar->getDegree();
        
        // for each neighbor of vertex v in the dualgraph
        // create a corresponding vertex in the spanning
        // tree if this has not been done yet.
        for (int i=0; i<dg; i++) {
            // get i-th neighbor of v
            tCQMDualGraphVertex* vchild = vpar->getNeighbor(i);
            
            // has it been visited before?   is it the common edge of "vpar"
            // and "vchild" a constrained one?
            if (!vchild->isMarked() && !vpar->isDualofConstrainedEdge(i)) {
                // the common edge of vpar and vchild is not the dual edge
                // of a constained edge in the mesh.
                
                // mark vchild as taken by the current spanning tree
                vchild->setMark(true);
                
                // create a corresponding vertex  in the spanning tree and add
                // it to the spanning tree.
                tCQMSpanningTreeVertex* child = sp->addVert(par,vchild);
                
                // put it into the working list of vertices
                vlist.push_back(child);
            }
        }
    }
    
    // unmark vertices
    v->getDualGraph()->clearMark();
    
    // Debugging. -- A.C.
    //sp->dump();
    return sp;
}

tCQMSpanningTree*
tCQMCompQuad::getSpanningTree2(tCQMDualGraphVertex* v)
{
    // create a spanning tree vertex corresponding to v
    tCQMSpanningTreeVertex* par = (tCQMSpanningTreeVertex*)
    new tCQMSpanningTreeVertex(v);
    
    // create a spanning tree with root par
    tCQMSpanningTree* sp = (tCQMSpanningTree*) new tCQMSpanningTree(par);
    
    // create a working list of spanning tree vertices
    std::list<tCQMSpanningTreeVertex*> vlist;
    
    // insert the root vertex into vlist
    vlist.push_back(par);
    
    // mark vertex to indicate that it has been taken by
    // the spanning tree.
    v->setMark(true);
    
    // insert all vertices reachable from v into vlist
    while (!vlist.empty()) {
        // get first vertex in the working list of vertices
        // of the current spanning tree.
        par = vlist.front();
        
        // remove "par" from the working list of vertices
        vlist.pop_front();
        
        // get its corresponding vertex in the dual graph
        tCQMDualGraphVertex* vpar = par->getVertex();
        
        // get degree of "v" in the dual graph
        int dg = vpar->getDegree();
        
        // for each neighbor of vertex v in the dualgraph
        // create a corresponding vertex in the spanning
        // tree if this has not been done yet.
        for (int i=0; i<dg; i++) {
            // get i-th neighbor of v
            tCQMDualGraphVertex* vchild = vpar->getNeighbor(i);
            
            // has it been visited before?   is it the common edge of "vpar"
            // and "vchild" a constrained one?
            if (!vchild->isMarked() && !vpar->isDualofConstrainedEdge(i)) {
                // the common edge of vpar and vchild is not the dual edge
                // of a constained edge in the mesh.
                
                // mark vchild as taken by the current spanning tree
                vchild->setMark(true);
                
                // create a corresponding vertex  in the spanning tree and add
                // it to the spanning tree.
                tCQMSpanningTreeVertex* child = sp->addVert(par,vchild);
                
                // put it into the working list of vertices
                vlist.push_back(child);
            }
        }
    }
    
    // unmark vertices
    v->getDualGraph()->clearMark();
    
    // Debugging. -- A.C.
    //sp->dump();
    return sp;
}


// -------------------------------------------------------------------
// Method getChildOfDegreeNParent()
//
// It returns a vertex of a  spanning tree whose parent is a vertex of
// degree "n",  where "n" is  an input parameter. The  vertex returned
// belongs  to a  subset of  vertices of  the spanning  tree  given as
// another input parameter.
// -------------------------------------------------------------------
tCQMSpanningTreeVertex*
tCQMCompQuad::getChildOfDegreeNParent(
                                      int n,
                                      const std::list<tCQMSpanningTreeVertex*>& lv
                                      )
{
    std::list<tCQMSpanningTreeVertex*>::const_iterator vi;
    tCQMSpanningTreeVertex* vpar;
    
    for (vi= lv.begin(); vi != lv.end(); ++vi) {
        if ((*vi)->getParent() != 0) {
            vpar = (*vi)->getParent();
            if (vpar->getDegree() == n) return *vi;
        }
        else return 0;
    }
    
    return 0;
}

// -------------------------------------------------------------------
// Method getQuadVertex()
// -------------------------------------------------------------------
tCQMIndVertex2
tCQMCompQuad::getQuadVertex(
                            QuadVertexSet& vs,
                            tCQMVertex2* v
                            )
{
    // verify if vertex is already in the vertex set
    QuadVertexSet::iterator vpos = vs.find(v);
    
    int id;
    if (vpos == vs.end()) {
        // vertex "v" is not in the set of vertices, so create a new index
        // for it and insert it there.
        id = getNewVertexIndex();
        vs.insert(make_pair(v,id));
    }
    else {
        id = vpos->second;
    }
    
    // return vertex information
    return tCQMIndVertex2(id,v->getPoint());
}

// -------------------------------------------------------------------
// Method isPointOnLine()
// -------------------------------------------------------------------
bool
tCQMCompQuad::isPointOnLine(
                            const Point2& p1,
                            const Point2& p2,
                            const Point2& v
                            )
{
    if ((p1.x() == v.x()) && (p1.y() == v.y())) return true;
    if ((p2.x() == v.x()) && (p2.y() == v.y())) return true;
    
    Line2 l(p1,p2);
    return l.oriented_side(v) == CGAL::ON_ORIENTED_BOUNDARY;
}

// -------------------------------------------------------------------
// Method isSegInLine()
// -------------------------------------------------------------------
bool
tCQMCompQuad::isSegInLine(
                          const Point2& p1,
                          const Point2& p2,
                          const Point2& v1,
                          const Point2& v2
                          )
{
    bool res1 = isPointOnLine(p1,p2,v1);
    bool res2 = isPointOnLine(p1,p2,v2);
    
    return res1 && res2;
}


// -------------------------------------------------------------------
// Method segLineIntersection()
// -------------------------------------------------------------------
bool
tCQMCompQuad::segLineIntersection(
                                  const Point2& p1,
                                  const Point2& p2,
                                  const Point2& v1,
                                  const Point2& v2,
                                  Point2& pt
                                  )
{
    // if intersection is the segment v1v2 itself then return false
    if ( isSegInLine( p1 , p2 , v1 , v2 ) ) return false ;
    
    // test if intersection is the point "v1"
    if ( isPointOnLine( p1 , p2 , v1 ) ) {
        pt = v1 ;
        return true ;
    }
    
    // test if intersection is the point "v2"
    if ( isPointOnLine( p1 , p2 , v2 ) ) {
        pt = v2 ;
        return true ;
    }
    
    // find intersection point
    mpf_class a( p2.x() - p1.x() , sizeof( double) * 8 ) ;
    mpf_class b( v2.x() - v1.x() , sizeof( double) * 8 ) ;
    mpf_class c( v1.x() - p1.x() , sizeof( double) * 8 ) ;
    mpf_class d( p2.y() - p1.y() , sizeof( double) * 8 ) ;
    mpf_class e( v2.y() - v1.y() , sizeof( double) * 8 ) ;
    mpf_class f( v1.y() - p1.y() , sizeof( double) * 8 ) ;
    
    mpf_class g( 0 , sizeof( double) * 128 ) ;
    
    g = ( d * b ) - ( e * a ) ;
    
    mpf_class zero( 0 , sizeof( double) * 128 ) ;
    
    // if line and segment are parallel then there is no intersection
    if ( cmp( g , zero ) == 0 ) {
        return false ;
    }
    
    mpf_class t( 0 , sizeof( double) * 128 ) ;
    
    t = ( (f * a) - (d * c) ) / g ;
    
    mpf_class one( 1 , sizeof( double) * 128 ) ;
    
    if ( ( cmp(t , zero) < 0 ) || ( cmp(t , one ) > 0 ) ) {
        return false ;
    }
    
    double dt = t.get_d() ;
    double db = b.get_d() ;
    double de = e.get_d() ;
    
    pt = Point2( v1.x() + db * dt , v1.y() + de * dt ) ;
    
    if ( IS_EQUAL( pt.x() , v1.x() ) && IS_EQUAL( pt.y() , v1.y() ) ) {
        pt = v1 ;
        return true ;
    }
    
    if ( IS_EQUAL( pt.x() , v2.x() ) && IS_EQUAL( pt.y() , v2.y() ) ) {
        pt = v2 ;
        return true ;
    }
    
    return true ;
}


// -------------------------------------------------------------------
// Method myAssign()
// -------------------------------------------------------------------
bool
tCQMCompQuad::myAssign(
                       Point2& p ,
                       const Segment2& s ,
                       const Line2& l
                       )
{
    if ( !assign( p , CGAL::intersection( s , l ) ) ) {
        Point2 p1 = s.vertex( 0 ) ;
        Point2 p2 = s.vertex( 1 ) ;
        Point2 q1 = l.point( 0 ) ;
        Point2 q2 = l.point( 1 ) ;
        
        return segLineIntersection( q1 , q2 , p1 , p2 , p ) ;
    }
    
    return true;
}

// -------------------------------------------------------------------
//
// Method halfSpaceIntersection()
//
// This method is the source of all sort of numerical problems in this
// code. The only way of getting rid of these problems is to use exact
// arithmetic. We are looking forward  to doing this. In the meantime,
// be aware that robustness issues are  part of this code and they may
// cause the code to produce incorrect output.
//
// -------------------------------------------------------------------
void
tCQMCompQuad::halfSpaceIntersection(
                                    const std::list<Point2>& li,
                                    const Point2& p1,
                                    const Point2& p2,
                                    std::list<Point2>& lo,
                                    CGAL::Oriented_side orient
                                    )
throw (tExceptionObject)
{
    // create an oriented line "l" from "p1" to "p2"
    Line2 l(p1,p2);
    
    // create an iterator for list of vertices of the polygon
    std::list<Point2>::const_iterator pi = li.begin();
    
    // get the first vertex of the polygon and advance iterator
    Point2 v1 = *pi;
    Point2 vf = v1;
    ++pi;
    
    // find out position of "v1" with  respect to the line "l" from "p1"
    // to "p2". Variable  "side1" will hold true if "v1"  is on line "l"
    // or on its left.
    CGAL::Oriented_side result = l.oriented_side(v1);
    bool side1 = (result == CGAL::ON_ORIENTED_BOUNDARY) ||
    (result != orient) || (v1 == p1) || (v1 == p2);
    bool sidef = side1;
    
    // loop over all segments except for the last one
    for (; pi != li.end(); ++pi) {
        // get next vertex of the input polygon
        Point2 v2 = *pi;
        
        // find out  position of  "v2" with respect  to the line  "l" from
        // "p1" to  "p2". Variable  "side2" will hold  true if "v2"  is on
        // line "l" or on its left.
        result = l.oriented_side( v2 ) ;
        bool side2 = ( result == CGAL::ON_ORIENTED_BOUNDARY ) ||
        ( result != orient ) || ( v2 == p1 ) || ( v2 == p2 ) ;
        
        // if "v1"  and "v2" are on  distinct sides of "l",  then line "l"
        // intersects the segment "v1v2".
        Point2 p;
        if ( side1 != side2 ) {
            // compute  intersection  point  between  line "l"  and  segment
            // "v1v2".
            
            if ( !segLineIntersection( p1 , p2 , v1 , v2 , p ) ) {
                std::stringstream ss (std::stringstream::in | std::stringstream::out);
                ss << "halfSpaceIntersection(): cannot recover from numerical instability";
                throw tExceptionObject(__FILE__,__LINE__,ss.str());
            }
            else {
                if ((!IS_EQUAL(v1.x(),p.x()) || !IS_EQUAL(v1.y(),p.y())) &&
                    (!IS_EQUAL(v2.x(),p.x()) || !IS_EQUAL(v2.y(),p.y()))) {
                    lo.push_back(p);
                }
                else if ((v2.x() == p.x()) && (v2.y() == p.y()) && !side2) {
                    side2 = true;
                }
                else if (!side1 && side2 && (v1.x() == p.x()) && (v1.y() == p.y())) {
                    lo.push_back(p);
                }
            }
        }
        
        // insert "v2" inside output list of vertices if "v2" is on "l" or
        // on its left side.
        if (side2) lo.push_back(v2);
        
        // update "v1" and "side1" for the next loop iteration
        v1 = v2;
        side1 = side2;
    }
    
    // if the last and first vertices are on distinct sides of "l", then
    // line "l" intersects the segment defined by them.
    if (side1 != sidef) {
        Point2 p;
        
        if ( !segLineIntersection( p1 , p2 , v1 , vf , p ) ) {
            std::stringstream ss (std::stringstream::in | std::stringstream::out);
            ss << "halfSpaceIntersection(): cannot recover from numerical instability";
            throw tExceptionObject(__FILE__,__LINE__,ss.str());
        }
        else {
            if ((!IS_EQUAL(v1.x(),p.x()) || !IS_EQUAL(v1.y(),p.y())) &&
                (!IS_EQUAL(vf.x(),p.x()) || !IS_EQUAL(vf.y(),p.y()))) {
                lo.push_back(p);
            }
            else if ((vf.x() == p.x()) && (vf.y() == p.y()) && !sidef) {
                sidef = true;
            }
            else if (!side1 && sidef && (v1.x() == p.x()) && (v1.y() == p.y())) {
                lo.push_back(p);
            }
        }
    }
    
    // if first vertex is on line "l"  or on its left side, we insert it
    // into the output list of vertices.
    if (sidef) {
        lo.push_back(vf);
    }
}

// -------------------------------------------------------------------
// Method createSteinerPoint()
// -------------------------------------------------------------------
tCQMIndVertex2
tCQMCompQuad::createSteinerPoint(const Point2& p)
{
    // increment Steiner point counter
    ++m_spoints;
    
    // return new Steiner point
    return tCQMIndVertex2(getNewVertexIndex(),p);
}

// -------------------------------------------------------------------
// Method createQuad()
// -------------------------------------------------------------------
void
tCQMCompQuad::createQuad(
                         const tCQMIndVertex2& p1,
                         const tCQMIndVertex2& p2,
                         const tCQMIndVertex2& p3,
                         const tCQMIndVertex2& p4,
                         bool e1,
                         bool e2,
                         bool e3,
                         bool e4,
                         std::list<tCQMQuadrilateral2*>* lq
                         )
throw (tExceptionObject)
{
    // create a new quadrilateral
    tCQMQuadrilateral2* quad = (tCQMQuadrilateral2*)
    new tCQMQuadrilateral2();
    
    // insert points in it
    quad->insert(0,p1);
    quad->insert(1,p2);
    quad->insert(2,p3);
    quad->insert(3,p4);
    
    // set edge status
    quad->setConstraint(0,e1);
    quad->setConstraint(1,e2);
    quad->setConstraint(2,e3);
    quad->setConstraint(3,e4);
    
    if (!quad->isSimple()) {
        std::stringstream ss (std::stringstream::in |
                              std::stringstream::out);
        ss << "createQuad(): quadrilateral is not simple";
        throw tExceptionObject(__FILE__,__LINE__,ss.str());
    }
    
    if (!quad->isConvex()) {
        std::stringstream ss (std::stringstream::in |
                              std::stringstream::out);
        ss << "createQuad(): quadrilateral is not convex";
        throw tExceptionObject(__FILE__,__LINE__,ss.str());
    }
    
    // insert new quadrilateral into the list of quadrilaterals
    lq->push_back(quad);
}

// -------------------------------------------------------------------
// Method insertVertexInContour()
// -------------------------------------------------------------------
void
tCQMCompQuad::insertVertexInContour(
                                    tCQMSpanningTreeVertex* v,
                                    QuadVertexSet& vs,
                                    std::list<tCQMQuadrilateral2*>* lq
                                    )
{
    // define variables to store edge status
    bool ec1, ec2, ec3, ec4;
    
    // get  a constrained  edge in  the mesh  face corresponding  to the
    // vertex associated with "v" in the dual graph.
    tCQMFace2* f = v->getVertex()->getFace();
    
    // if face "f" is not adjacent  to any constrained edge, we throw an
    // exception.
    if (!f->isAdjToConstrainedEdge()) return;
    
    // look for a constrained edge in this face
    tCQMHalfEdge2* he = f->getHalfEdge();
    if (v->hasSteinerPoint()) {
        Point2 pt = v->getSteinerPointSuc();
        do {
            Point2 ptaux = he->getVertex()->getPoint();
            if ((pt.x() == ptaux.x()) && (pt.y() == ptaux.y())) {
                break;
            }
            he = he->getNext();
        } while (true);
        
        ec3 = false;
        ec4 = false;
    }
    else {
        while (!he->getEdge()->isConstrained()) {
            he = he->getNext();
        }
        
        ec3 = he->getNext()->getEdge()->isConstrained();
        ec4 = he->getPrev()->getEdge()->isConstrained();
    }
    
    // get edge status
    ec1 = true;
    ec2 = true;
    
    // get the two points defining the edge containing "he"
    tCQMIndVertex2 p1 = getQuadVertex(vs,he->getVertex());
    tCQMIndVertex2 p3 = getQuadVertex(vs,he->getNext()->getVertex());
    
    // get the remainder point in the triangle
    tCQMIndVertex2 p4;
    if (v->hasSteinerPoint()) {
        p4 = v->getSteinerPoint();
    }
    else {
        p4 = getQuadVertex(vs,he->getPrev()->getVertex());
    }
    
    // compute the new point and  split the face that shares the contour
    // edge with "f", if any.
    double x1 = p1.getPoint().x();
    double y1 = p1.getPoint().y();
    double x3 = p3.getPoint().x();
    double y3 = p3.getPoint().y();
    double x = (x1 + x3) / 2;
    double y = (y1 + y3) / 2;
    double dx = 0;
    double dy = 0;
    if (!he->getEdge()->isOnBoundary()) {
        Point2 c = he->getEdge()->getMate(he)->getFace()->getBarycenter();
        dx = 0.1 * (c.x() - x);
        dy = 0.1 * (c.y() - y);
    }
    
    if (fabs(dx) > MAX_EDGE_DISP) {
        dx = (dx < 0) ? -MAX_EDGE_DISP : MAX_EDGE_DISP;
    }
    
    if (fabs(dy) > MAX_EDGE_DISP) {
        dy = (dy < 0) ? -MAX_EDGE_DISP : MAX_EDGE_DISP;
    }
    
    // new point
    tCQMIndVertex2 p2 = createSteinerPoint(Point2(x + dx, y + dy));
    
    // split edge  containing point and subdivide face  that shares this
    // edge with face "f", if any.
    tCQMVertex2* newv = m_mesh->spFmkFE(he,p2.getPoint());
    
    // insert  new   vertex  in   the  set  of   the  vertices   of  the
    // quadrangulation.
    vs.insert(make_pair(newv,p2.getVertexId()));
    
    // create quadrilateral
    try {
        createQuad(p1,p2,p3,p4,ec1,ec2,ec3,ec4,lq);
    }
    catch (const tExceptionObject& xpt) {
        treatException(xpt);
        exit(0);
    }
}

// -------------------------------------------------------------------
// Method findOneQuadrilateral()
// -------------------------------------------------------------------
tCQMSpanningTreeVertex*
tCQMCompQuad::findOneQuadrilateral(
                                   const std::list<tCQMSpanningTreeVertex*>& lv)
{
    std::list<tCQMSpanningTreeVertex*>::const_iterator vi;
    for (vi= lv.begin(); vi != lv.end(); ++vi) {
        if ((*vi)->isQuadrilateral() != 0) {
            return *vi;
        }
    }
    
    return 0;
}

// -------------------------------------------------------------------
// Method findOnePentagon()
// -------------------------------------------------------------------
tCQMSpanningTreeVertex*
tCQMCompQuad::findOnePentagon(const std::list<tCQMSpanningTreeVertex*>& lv)
{
    std::list<tCQMSpanningTreeVertex*>::const_iterator vi;
    for (vi= lv.begin(); vi != lv.end(); ++vi) {
        if ((*vi)->isPentagon() != 0) {
            return *vi;
        }
    }
    
    return 0;
}

// -------------------------------------------------------------------
// Method getPentagon()
// -------------------------------------------------------------------
tCQMPolygon2*
tCQMCompQuad::getPentagon(
                          tCQMSpanningTreeVertex* v1,
                          tCQMSpanningTreeVertex* v2,
                          int& id,
                          QuadVertexSet& vs
                          )
{
    // get  the mesh  faces corresponding  to the  dual graph  vertex in
    // "v1" and "v2".
    tCQMFace2* f1 = v1->getVertex()->getFace();
    tCQMFace2* f2 = v2->getVertex()->getFace();
    
    // get the common edge of "f1" and "f2"
    tCQMEdge2* e = f1->getCommonEdge(f2);
    
    // get the half-edge of "e" in "f1"
    tCQMHalfEdge2* he1 = e->getHalfEdge();
    tCQMHalfEdge2* he2 = e->getMate(he1);
    if (he1->getFace() != f1) {
        he2 = he1;
        he1 = e->getMate(he1);
    }
    
    // compute pentagon vertices
    
    // "p1" and "p4" are the endpoints  of the edge "e", and "p5" is the
    // other vertex of the face corresponding to "v1".
    tCQMIndVertex2 p1 = getQuadVertex(vs, he1->getVertex());
    tCQMIndVertex2 p4 = getQuadVertex(vs, he2->getVertex());
    tCQMIndVertex2 p5;
    bool ec4, ec5;
    if (v1->hasSteinerPoint()) {
        p5 = v1->getSteinerPoint();
        ec4 = false;
        ec5 = false;
    }
    else {
        p5 = getQuadVertex(vs, he1->getPrev()->getVertex());
        ec4 = he1->getNext()->getEdge()->isConstrained();
        ec5 = he1->getPrev()->getEdge()->isConstrained();
    }
    
    // "p2" and  "p3" are  the other two  vertices of  the quadrilateral
    // corresponding to  "v2", and  "id" is the  index of one  vertex of
    // this quadrilateral that is  incident to the face corresponding to
    // the parent of "v2" in the spanning tree.
    tCQMIndVertex2 p2, p3;
    bool ec1, ec2, ec3;
    if (he2->getNext() == v2->getQuadHalfEdge()) {
        p2 = v2->getQuadSteinerPoint();
        p3 = getQuadVertex(vs, he2->getPrev()->getVertex());
        ec1 = false;
        ec2 = false;
        ec3 = he2->getPrev()->getEdge()->isConstrained();
        id = 3;
    }
    else {
        p2 = getQuadVertex(vs, he2->getPrev()->getVertex());
        p3 = v2->getQuadSteinerPoint();
        ec1 = he2->getNext()->getEdge()->isConstrained();
        ec2 = false;
        ec3 = false;
        id = 0;
    }
    
    // create the pentagon
    tCQMPolygon2* p = (tCQMPolygon2 *) new tCQMPolygon2(5);
    
    // insert the vertices in the pentagon in counterclockwise order
    p->insert(0,p1);
    p->insert(1,p2);
    p->insert(2,p3);
    p->insert(3,p4);
    p->insert(4,p5);
    
    p->setConstraint(0,ec1);
    p->setConstraint(1,ec2);
    p->setConstraint(2,ec3);
    p->setConstraint(3,ec4);
    p->setConstraint(4,ec5);
    
    return p;
}

// -------------------------------------------------------------------
// Method getPentagon()
// -------------------------------------------------------------------
tCQMPolygon2*
tCQMCompQuad::getPentagon(
                          tCQMSpanningTreeVertex* v1,
                          tCQMSpanningTreeVertex* v2,
                          tCQMSpanningTreeVertex* v3,
                          QuadVertexSet& vs
                          )
{
    // get  the mesh  faces corresponding  to the  dual graph  vertex in
    // "v1", "v2" and "v3".
    tCQMFace2* f1 = v1->getVertex()->getFace();
    tCQMFace2* f2 = v2->getVertex()->getFace();
    tCQMFace2* f3 = v3->getVertex()->getFace();
    
    // get the common edge of "f1" and "f2"
    tCQMEdge2* e1 = f1->getCommonEdge(f2);
    
    // get the common edge of "f2" and "f3"
    tCQMEdge2* e2 = f2->getCommonEdge(f3);
    
    // get the half-edge of "e1" in "f1"
    tCQMHalfEdge2* he1 = e1->getHalfEdge();
    tCQMHalfEdge2* he2 = e1->getMate(he1);
    if (he1->getFace() != f1) {
        he2 = he1;
        he1 = e1->getMate(he1);
    }
    
    // get the half-edge of "e2" in "f2"
    tCQMHalfEdge2* he3 = e2->getHalfEdge();
    tCQMHalfEdge2* he4 = e2->getMate(he3);
    if (he3->getFace() != f2) {
        he4 = he3;
        he3 = e2->getMate(he3);
    }
    
    // compute all vertices  of the pentagon.  Either "p1"  or "p4" is a
    // reflex  vertex.   Vertex  "v1"  corresponds to  triangle  p1p4p5,
    // vertex "v2" corresponds to  triangle p1p2p4 or p1p3p4, and vertex
    // "v3" corresponds to triangle p2p3p4 or p1p2p3.
    tCQMIndVertex2 p1 = getQuadVertex(vs, he1->getVertex());
    tCQMIndVertex2 p2, p3;
    bool ec1, ec2, ec3;
    if (he2->getNext() == he3) {
        // vertices  "v2"  and "v3"  correspond  to  triangles p1p3p4  and
        // p1p2p3.
        p2 = getQuadVertex(vs, he4->getPrev()->getVertex());
        p3 = getQuadVertex(vs, he4->getVertex());
        
        ec1 = he4->getNext()->getEdge()->isConstrained();
        ec2 = he4->getPrev()->getEdge()->isConstrained();
        ec3 = he2->getPrev()->getEdge()->isConstrained();
    }
    else {
        // vertices  "v2"  and "v3"  correspond  to  triangles p1p2p4  and
        // p2p3p4.
        p2 = getQuadVertex(vs, he3->getVertex());
        p3 = getQuadVertex(vs, he4->getPrev()->getVertex());
        
        ec1 = he2->getNext()->getEdge()->isConstrained();
        ec2 = he4->getNext()->getEdge()->isConstrained();
        ec3 = he4->getPrev()->getEdge()->isConstrained();
    }
    
    tCQMIndVertex2 p4 = getQuadVertex(vs, he2->getVertex());
    tCQMIndVertex2 p5;
    bool ec4, ec5;
    if (v1->hasSteinerPoint()) {
        p5 = v1->getSteinerPoint();
        ec4 = false;
        ec5 = false;
    }
    else {
        p5 = getQuadVertex(vs, he1->getPrev()->getVertex());
        ec4 = he1->getNext()->getEdge()->isConstrained();
        ec5 = he1->getPrev()->getEdge()->isConstrained();
    }
    
    // create the pentagon
    tCQMPolygon2* p = (tCQMPolygon2 *) new tCQMPolygon2(5);
    
    p->insert(0,p1);
    p->insert(1,p2);
    p->insert(2,p3);
    p->insert(3,p4);
    p->insert(4,p5);
    
    p->setConstraint(0,ec1);
    p->setConstraint(1,ec2);
    p->setConstraint(2,ec3);
    p->setConstraint(3,ec4);
    p->setConstraint(4,ec5);
    
    return p;
}

// -------------------------------------------------------------------
// Method getHexagon()
// -------------------------------------------------------------------
tCQMPolygon2*
tCQMCompQuad::getHexagon(
                         tCQMSpanningTreeVertex* v1,
                         tCQMSpanningTreeVertex* v2,
                         tCQMSpanningTreeVertex* v3,
                         QuadVertexSet& vs
                         )
{
    // get  the mesh  faces corresponding  to the  dual graph  vertex in
    // "v1", "v2" and "v3".
    tCQMFace2* f1 = v1->getVertex()->getFace();
    tCQMFace2* f2 = v2->getVertex()->getFace();
    tCQMFace2* f3 = v3->getVertex()->getFace();
    
    // get the common edge of "f1" and "f2"
    tCQMEdge2* e1 = f1->getCommonEdge(f2);
    
    // get the common edge of "f2" and "f3"
    tCQMEdge2* e2 = f2->getCommonEdge(f3);
    
    // get the half-edge of "e1" in "f1"
    tCQMHalfEdge2* he1 = e1->getHalfEdge();
    tCQMHalfEdge2* he2 = e1->getMate(he1);
    if (he1->getFace() != f1) {
        he2 = he1;
        he1 = e1->getMate(he1);
    }
    
    // get the half-edge of "e2" in "f2"
    tCQMHalfEdge2* he3 = e2->getHalfEdge();
    tCQMHalfEdge2* he4 = e2->getMate(he3);
    if (he3->getFace() != f2) {
        he4 = he3;
        he3 = e2->getMate(he3);
    }
    
    // create the hexagon
    tCQMPolygon2* p = (tCQMPolygon2 *) new tCQMPolygon2(6);
    
    // compute the hexagon vertices
    tCQMIndVertex2 p1, p2, p3, p4, p5, p6;
    bool ec1, ec2, ec3, ec4, ec5, ec6;
    if (v2->isQuadrilateral()) {
        p1 = getQuadVertex(vs, he1->getVertex());
        p5 = getQuadVertex(vs, he2->getVertex());
        if (v1->hasSteinerPoint()) {
            p6 = v1->getSteinerPoint();
            ec5 = false;
            ec6 = false;
        }
        else {
            p6 = getQuadVertex(vs, he1->getPrev()->getVertex());
            ec5 = he1->getNext()->getEdge()->isConstrained();
            ec6 = he1->getPrev()->getEdge()->isConstrained();
        }
        
        // check if  either "p2" or "p4"  will correspond to  the point in
        // the middle  of an  edge in the  face corresponding to  the dual
        // vertex associated with "v2".
        if (he2->getNext() == v2->getQuadHalfEdge()) {
            p2 = v2->getQuadSteinerPoint();
            ec1 = false;
            ec2 = false;
            p3 = getQuadVertex(vs, he3->getVertex());
            if (v3->hasSteinerPoint()) {
                p4 = v3->getSteinerPoint();
                ec3 = false;
                ec4 = false;
            }
            else {
                p4 = getQuadVertex(vs, he4->getPrev()->getVertex());
                ec3 = he4->getNext()->getEdge()->isConstrained();
                ec4 = he4->getPrev()->getEdge()->isConstrained();
            }
        }
        else {
            if (v3->hasSteinerPoint()) {
                p2 =  v3->getSteinerPoint();
                ec1 = false;
                ec2 = false;
            }
            else {
                p2 = getQuadVertex(vs, he4->getPrev()->getVertex());
                ec1 = he4->getNext()->getEdge()->isConstrained();
                ec2 = he4->getPrev()->getEdge()->isConstrained();
            }
            p3 = getQuadVertex(vs, he4->getVertex());
            p4 = v2->getQuadSteinerPoint();
            ec3 = false;
            ec4 = false;
        }
    }
    else {
        p1 = getQuadVertex(vs, he1->getVertex());
        p5 = getQuadVertex(vs, he1->getNext()->getVertex());
        if (v1->hasSteinerPoint()) {
            p6 =  v1->getSteinerPoint();
            ec5 = false;
            ec6 = false;
        }
        else {
            p6 = getQuadVertex(vs, he1->getPrev()->getVertex());
            ec5 = he1->getNext()->getEdge()->isConstrained();
            ec6 = he1->getPrev()->getEdge()->isConstrained();
        }
        
        // check if either "p1" or "p5"  is the common vertex of the faces
        // corresponding to  the dual vertices associated  with "v1", "v2"
        // and "v3".
        if (he1->getVertex() == he4->getNext()->getVertex()) {
            // "p1" is the  common vertex of the faces  corresponding to the
            // dual vertices associated with "v1", "v2" and "v3".
            p4 = getQuadVertex(vs, he2->getPrev()->getVertex());
            ec4 = he2->getPrev()->getEdge()->isConstrained();
            
            // check if either "p2" or  "p3" will correspond to the point in
            // the middle of  an edge in the face  corresponding to the dual
            // vertex associated with "v3".
            if (he4->getNext() == v3->getQuadHalfEdge()) {
                p2 = v3->getQuadSteinerPoint();
                p3 = getQuadVertex(vs, he4->getPrev()->getVertex());
                ec1 = false;
                ec2 = false;
                ec3 = he4->getPrev()->getEdge()->isConstrained();
            }
            else {
                p2 = getQuadVertex(vs, he4->getPrev()->getVertex());
                p3 = v3->getQuadSteinerPoint();
                ec1 = he4->getNext()->getEdge()->isConstrained();
                ec2 = false;
                ec3 = false;
            }
        }
        else {
            // "p5" is the  common vertex of the faces  corresponding to the
            // dual vertices associated with "v1", "v2" and "v3".
            p2 = getQuadVertex(vs, he2->getPrev()->getVertex());
            ec1 = he2->getNext()->getEdge()->isConstrained();
            
            // check if either "p3" or  "p4" will correspond to the point in
            // the middle of  an edge in the face  corresponding to the dual
            // vertex associated with "v3".
            if (he4->getNext() == v3->getQuadHalfEdge()) {
                p3 = v3->getQuadSteinerPoint();
                p4 = getQuadVertex(vs, he4->getPrev()->getVertex());
                ec2 = false;
                ec3 = false;
                ec4 = he4->getPrev()->getEdge()->isConstrained();
            }
            else {
                p3 = getQuadVertex(vs, he4->getPrev()->getVertex());
                p4 = v3->getQuadSteinerPoint();
                ec2 = he4->getNext()->getEdge()->isConstrained();
                ec3 = false;
                ec4 = false;
            }
        }
    }
    
    p->insert(0,p1);
    p->insert(1,p2);
    p->insert(2,p3);
    p->insert(3,p4);
    p->insert(4,p5);
    p->insert(5,p6);
    
    p->setConstraint(0,ec1);
    p->setConstraint(1,ec2);
    p->setConstraint(2,ec3);
    p->setConstraint(3,ec4);
    p->setConstraint(4,ec5);
    p->setConstraint(5,ec6);
    
    return p;
}

// -------------------------------------------------------------------
// Method getHexagon()
// -------------------------------------------------------------------
tCQMPolygon2*
tCQMCompQuad::getHexagon(
                         tCQMSpanningTreeVertex* v1,
                         tCQMSpanningTreeVertex* v2,
                         tCQMSpanningTreeVertex* v3,
                         tCQMSpanningTreeVertex* v4,
                         QuadVertexSet& vs
                         )
{
    
    // get  the mesh  faces corresponding  to the  dual graph  vertex in
    // "v1", "v2", "v3", and "v4".
    tCQMFace2* f1 = v1->getVertex()->getFace();
    tCQMFace2* f2 = v2->getVertex()->getFace();
    tCQMFace2* f3 = v3->getVertex()->getFace();
    tCQMFace2* f4 = v4->getVertex()->getFace();
    
    // get the common edge of "f1" and "f2"
    tCQMEdge2* e1 = f1->getCommonEdge(f2);
    
    // get the common edge of "f2" and "f4"
    tCQMEdge2* e2 = f2->getCommonEdge(f4);
    
    // get the half-edge of "e1" in "f1"
    tCQMHalfEdge2* he1 = e1->getHalfEdge();
    tCQMHalfEdge2* he2 = e1->getMate(he1);
    if (he1->getFace() != f1) {
        he2 = he1;
        he1 = e1->getMate(he1);
    }
    
    // get the half-edge of "e2" in "f2"
    tCQMHalfEdge2* he3 = e2->getHalfEdge();
    tCQMHalfEdge2* he4 = e2->getMate(he3);
    if (he3->getFace() != f2) {
        he4 = he3;
        he3 = e2->getMate(he3);
    }
    
    tCQMPolygon2* p;
    tCQMIndVertex2 p1, p2, p3, p4, p5, p6;
    bool ec1 = false,
    ec2 = false,
    ec3 = false,
    ec4 = false,
    ec5 = false,
    ec6 = false;
    
    // find out if "v3" is a child of either "v2" or "v4"
    if (v2->getNumChildren() == 1) {
        // "v3" is a child of "v4"
        
        // get the common edge of "f3" and "f4"
        tCQMEdge2* e3 = f3->getCommonEdge(f4);
        
        // get the half-edge of "e3" in "f4"
        tCQMHalfEdge2* he5 = e3->getHalfEdge();
        tCQMHalfEdge2* he6 = e3->getMate(he5);
        if (he5->getFace() != f4) {
            he6 = he5;
            he5 = e3->getMate(he5);
        }
        
        // if either  "f1" or  "f2" shares  an edge with  "f3", we  have a
        // quadrilateral instead of an hexagon.
        if ((f1->getCommonEdge(f3) != 0) || (f2->getCommonEdge(f3) != 0)) {
            
            if (he4->getNext() == he5) {
                // compute the hexagon vertices
                p1 = getQuadVertex(vs, he4->getPrev()->getVertex());
                ec1 = he4->getPrev()->getEdge()->isConstrained();
                
                if (f1->getCommonEdge(f3) != 0) {
                    p2 = getQuadVertex(vs, he3->getNext()->getVertex());
                    p3 = getQuadVertex(vs, he1->getNext()->getVertex());
                    ec2 = he3->getNext()->getEdge()->isConstrained();
                    ec3 = he1->getNext()->getEdge()->isConstrained();
                }
                else {
                    p2 = getQuadVertex(vs, he1->getNext()->getVertex());
                    p3 = getQuadVertex(vs, he1->getPrev()->getVertex());
                    ec2 = he1->getNext()->getEdge()->isConstrained();
                    ec3 = he1->getPrev()->getEdge()->isConstrained();
                }
                
                p4 = getQuadVertex(vs, he6->getPrev()->getVertex());
                ec4 = he6->getPrev()->getEdge()->isConstrained();
            }
            else {
                // compute the hexagon vertices
                p1 = getQuadVertex(vs, he4->getNext()->getVertex());
                ec1 = he4->getNext()->getEdge()->isConstrained();
                
                if (f1->getCommonEdge(f3) != 0) {
                    p2 = getQuadVertex(vs, he6->getNext()->getVertex());
                    p4 = getQuadVertex(vs, he3->getPrev()->getVertex());
                    ec2 = he6->getNext()->getEdge()->isConstrained();
                    ec4 = he3->getPrev()->getEdge()->isConstrained();
                }
                else {
                    p2 = getQuadVertex(vs, he6->getNext()->getVertex());
                    p4 = getQuadVertex(vs, he1->getPrev()->getVertex());
                    ec2 = he6->getNext()->getEdge()->isConstrained();
                    ec4 = he1->getPrev()->getEdge()->isConstrained();
                }
                
                p3 = getQuadVertex(vs, he1->getPrev()->getVertex());
                ec3 = he1->getPrev()->getEdge()->isConstrained();
            }
            
            // create a quadrilateral
            p = (tCQMPolygon2 *) new tCQMPolygon2(4);
            
            p->insert(0,p1);
            p->insert(1,p2);
            p->insert(2,p3);
            p->insert(3,p4);
            
            p->setConstraint(0,ec1);
            p->setConstraint(1,ec2);
            p->setConstraint(2,ec3);
            p->setConstraint(3,ec4);
        }
        else {
            // "f1", "f2", "f3" and "f4" define an hexagon
            
            // compute the hexagon vertices
            p1 = getQuadVertex(vs, he1->getVertex());
            p5 = getQuadVertex(vs, he2->getVertex());
            
            if (v1->hasSteinerPoint()) {
                p6 = v1->getSteinerPoint();
                ec5 = false;
                ec6 = false;
            }
            else {
                p6 = getQuadVertex(vs, he1->getPrev()->getVertex());
                ec5 = he1->getNext()->getEdge()->isConstrained();
                ec6 = he1->getPrev()->getEdge()->isConstrained();
            }
            
            if (he2->getNext() != he3) {
                p2 = getQuadVertex(vs, he2->getPrev()->getVertex());
                ec1 = he2->getNext()->getEdge()->isConstrained();
                if (he4->getNext() != he5) {
                    p3 = getQuadVertex(vs, he4->getPrev()->getVertex());
                    if (v3->hasSteinerPoint()) {
                        p4 = v3->getSteinerPoint();
                        ec3 = false;
                        ec4 = false;
                    }
                    else {
                        p4 = getQuadVertex(vs, he6->getPrev()->getVertex());
                        ec3 = he6->getNext()->getEdge()->isConstrained();
                        ec4 = he6->getPrev()->getEdge()->isConstrained();
                    }
                    ec2 = he4->getNext()->getEdge()->isConstrained();
                }
                else {
                    if (v3->hasSteinerPoint()) {
                        p3 = v3->getSteinerPoint();
                        ec2 = false;
                        ec3 = false;
                    }
                    else {
                        p3 = getQuadVertex(vs, he6->getPrev()->getVertex());
                        ec2 = he6->getNext()->getEdge()->isConstrained();
                        ec3 = he6->getPrev()->getEdge()->isConstrained();
                    }
                    p4 = getQuadVertex(vs, he4->getPrev()->getVertex());
                    ec4 = he4->getPrev()->getEdge()->isConstrained();
                }
            }
            else {
                p4 = getQuadVertex(vs, he2->getPrev()->getVertex());
                ec4 = he3->getNext()->getEdge()->isConstrained();
                if (he4->getNext() != he5) {
                    p2 = getQuadVertex(vs, he4->getPrev()->getVertex());
                    if (v3->hasSteinerPoint()) {
                        p3 = v3->getSteinerPoint();
                        ec2 = false;
                        ec3 = false;
                    }
                    else {
                        p3 = getQuadVertex(vs, he6->getPrev()->getVertex());
                        ec2 = he6->getNext()->getEdge()->isConstrained();
                        ec3 = he6->getPrev()->getEdge()->isConstrained();
                    }
                    ec1 = he4->getNext()->getEdge()->isConstrained();
                }
                else {
                    if (v3->hasSteinerPoint()) {
                        p2 = v3->getSteinerPoint();
                        ec1 = false;
                        ec2 = false;
                    }
                    else {
                        p2 = getQuadVertex(vs, he6->getPrev()->getVertex());
                        ec1 = he6->getNext()->getEdge()->isConstrained();
                        ec2 = he6->getPrev()->getEdge()->isConstrained();
                    }
                    
                    p3 = getQuadVertex(vs, he4->getPrev()->getVertex());
                    ec3 = he4->getPrev()->getEdge()->isConstrained();
                }
            }
            
            // create the hexagon
            p = (tCQMPolygon2 *) new tCQMPolygon2(6);
            
            p->insert(0,p1);
            p->insert(1,p2);
            p->insert(2,p3);
            p->insert(3,p4);
            p->insert(4,p5);
            p->insert(5,p6);
            
            p->setConstraint(0,ec1);
            p->setConstraint(1,ec2);
            p->setConstraint(2,ec3);
            p->setConstraint(3,ec4);
            p->setConstraint(4,ec5);
            p->setConstraint(5,ec6);
        }
    }
    else {
        // "v3" is a child of "v2"
        
        // get the common edge of "f2" and "f3"
        tCQMEdge2* e3 = f2->getCommonEdge(f3);
        
        // get the half-edge of "e3" in "f2"
        tCQMHalfEdge2* he5 = e3->getHalfEdge();
        tCQMHalfEdge2* he6 = e3->getMate(he5);
        if (he5->getFace() != f2) {
            he6 = he5;
            he5 = e3->getMate(he5);
        }
        
        // if either  "f1" or  "f2" shares  an edge with  "f3", we  have a
        // quadrilateral instead of an hexagon.
        if ((f1->getCommonEdge(f3) != 0) || (f1->getCommonEdge(f4) != 0)
            || (f3->getCommonEdge(f4) != 0)) {
            // "f1", "f2",  "f3" and "f4" define a  quadrilateral instead of
            // an hexagon.
            
            if (f1->getCommonEdge(f3) != 0) {
                // "f1" and "f3" share an edge
                if (he2->getNext() == he5) {
                    // "v1" is the left child of "v2"
                    p1 = getQuadVertex(vs, he6->getPrev()->getVertex());
                    p2 = getQuadVertex(vs, he4->getNext()->getVertex());
                    p3 = getQuadVertex(vs, he4->getPrev()->getVertex());
                    p4 = getQuadVertex(vs, he1->getNext()->getVertex());
                    
                    ec1 = he6->getPrev()->getEdge()->isConstrained();
                    ec2 = he4->getNext()->getEdge()->isConstrained();
                    ec3 = he4->getPrev()->getEdge()->isConstrained();
                    ec4 = he1->getNext()->getEdge()->isConstrained();
                }
                else {
                    // "v1" is the right child of "v2"
                    p1 = getQuadVertex(vs, he1->getPrev()->getVertex());
                    p2 = getQuadVertex(vs, he4->getNext()->getVertex());
                    p3 = getQuadVertex(vs, he4->getPrev()->getVertex());
                    p4 = getQuadVertex(vs, he6->getNext()->getVertex());
                    
                    ec1 = he1->getPrev()->getEdge()->isConstrained();
                    ec2 = he4->getNext()->getEdge()->isConstrained();
                    ec3 = he4->getPrev()->getEdge()->isConstrained();
                    ec4 = he6->getNext()->getEdge()->isConstrained();
                }
            }
            else if (f1->getCommonEdge(f4) != 0) {
                // "f1" and "f4" share an edge
                if (he2->getNext() == he5) {
                    // "v1" is the left child of "v2"
                    p1 = getQuadVertex(vs, he6->getNext()->getVertex());
                    if (v3->hasSteinerPoint()) {
                        p2 = v3->getSteinerPoint();
                        ec1 = false;
                        ec2 = false;
                    }
                    else {
                        p2 = getQuadVertex(vs, he6->getPrev()->getVertex());
                        ec1 = he6->getNext()->getEdge()->isConstrained();
                        ec2 = he6->getPrev()->getEdge()->isConstrained();
                    }
                    
                    p3 = getQuadVertex(vs, he4->getNext()->getVertex());
                    p4 = getQuadVertex(vs, he1->getPrev()->getVertex());
                    
                    ec3 = he4->getNext()->getEdge()->isConstrained();
                    ec4 = he1->getPrev()->getEdge()->isConstrained();
                }
                else {
                    // "v1" is the right child of "v2"
                    p1 = getQuadVertex(vs, he1->getNext()->getVertex());
                    p2 = getQuadVertex(vs, he4->getPrev()->getVertex());
                    p3 = getQuadVertex(vs, he6->getNext()->getVertex());
                    
                    ec1 = he1->getNext()->getEdge()->isConstrained();
                    ec2 = he2->getPrev()->getEdge()->isConstrained();
                    
                    if (v3->hasSteinerPoint()) {
                        p4 = v3->getSteinerPoint();
                        ec3 = false;
                        ec4 = false;
                    }
                    else {
                        p4 = getQuadVertex(vs, he6->getPrev()->getVertex());
                        ec3 = he6->getNext()->getEdge()->isConstrained();
                        ec4 = he6->getPrev()->getEdge()->isConstrained();
                    }
                }
            }
            else if (f3->getCommonEdge(f4) != 0) {
                if (he2->getNext() == he5) {
                    // "v1" is the left child of "v2"
                    p1 = getQuadVertex(vs, he6->getNext()->getVertex());
                    p2 = getQuadVertex(vs, he4->getPrev()->getVertex());
                    p3 = getQuadVertex(vs, he1->getNext()->getVertex());
                    
                    ec1 = he6->getNext()->getEdge()->isConstrained();
                    ec2 = he4->getPrev()->getEdge()->isConstrained();
                    
                    if (v1->hasSteinerPoint()) {
                        p4 = v1->getSteinerPoint();
                        ec3 = false;
                        ec4 = false;
                    }
                    else {
                        p4 = getQuadVertex(vs, he1->getPrev()->getVertex());
                        ec3 = he1->getNext()->getEdge()->isConstrained();
                        ec4 = he1->getPrev()->getEdge()->isConstrained();
                    }
                }
                else {
                    // "v1" is the right child of "v2"
                    p1 = getQuadVertex(vs, he1->getNext()->getVertex());
                    if (v1->hasSteinerPoint()) {
                        p2 = v1->getSteinerPoint();
                        ec1 = false;
                        ec2 = false;
                    }
                    else {
                        p2 = getQuadVertex(vs, he1->getPrev()->getVertex());
                        ec1 = he1->getNext()->getEdge()->isConstrained();
                        ec2 = he1->getPrev()->getEdge()->isConstrained();
                    }
                    
                    p3 = getQuadVertex(vs, he4->getNext()->getVertex());
                    p4 = getQuadVertex(vs, he6->getPrev()->getVertex());
                    
                    ec3 = he4->getNext()->getEdge()->isConstrained();
                    ec4 = he6->getPrev()->getEdge()->isConstrained();
                }
            }
            
            // create a quadrilateral
            p = (tCQMPolygon2 *) new tCQMPolygon2(4);
            
            p->insert(0,p1);
            p->insert(1,p2);
            p->insert(2,p3);
            p->insert(3,p4);
            
            p->setConstraint(0,ec1);
            p->setConstraint(1,ec2);
            p->setConstraint(2,ec3);
            p->setConstraint(3,ec4);
        }
        else {
            // "f1", "f2",  "f3" and "f4" define an hexagon
            
            // compute the hexagon vertices
            p1 = getQuadVertex(vs, he1->getVertex());
            p5 = getQuadVertex(vs, he2->getVertex());
            
            if (v1->hasSteinerPoint()) {
                p6 = v1->getSteinerPoint();
                ec5 = false;
                ec6 = false;
            }
            else {
                p6 = getQuadVertex(vs, he1->getPrev()->getVertex());
                ec5 = he1->getNext()->getEdge()->isConstrained();
                ec6 = he1->getPrev()->getEdge()->isConstrained();
            }
            
            if (he2->getNext() != he3) {
                if (v3->hasSteinerPoint()) {
                    p2 = v3->getSteinerPoint();
                    ec1 = false;
                    ec2 = false;
                }
                else {
                    p2 = getQuadVertex(vs, he6->getPrev()->getVertex());
                    ec1 = he6->getNext()->getEdge()->isConstrained();
                    ec2 = he6->getPrev()->getEdge()->isConstrained();
                }
                
                p3 = getQuadVertex(vs, he6->getVertex());
                p4 = getQuadVertex(vs, he4->getPrev()->getVertex());
                ec3 = he4->getNext()->getEdge()->isConstrained();
                ec4 = he4->getPrev()->getEdge()->isConstrained();
            }
            else {
                p2 = getQuadVertex(vs, he4->getPrev()->getVertex());
                p3 = getQuadVertex(vs, he5->getVertex());
                ec1 = he4->getNext()->getEdge()->isConstrained();
                ec2 = he4->getPrev()->getEdge()->isConstrained();
                if (v3->hasSteinerPoint()) {
                    p4 = v3->getSteinerPoint();
                    ec3 = false;
                    ec4 = false;
                }
                else {
                    p4 = getQuadVertex(vs, he6->getPrev()->getVertex());
                    ec3 = he6->getNext()->getEdge()->isConstrained();
                    ec4 = he6->getPrev()->getEdge()->isConstrained();
                }
            }
            
            // create the hexagon
            p = (tCQMPolygon2 *) new tCQMPolygon2(6);
            
            p->insert(0,p1);
            p->insert(1,p2);
            p->insert(2,p3);
            p->insert(3,p4);
            p->insert(4,p5);
            p->insert(5,p6);
            
            p->setConstraint(0,ec1);
            p->setConstraint(1,ec2);
            p->setConstraint(2,ec3);
            p->setConstraint(3,ec4);
            p->setConstraint(4,ec5);
            p->setConstraint(5,ec6);
        }
    }
    
    return p;
}

// -------------------------------------------------------------------
// Method classifyPentagon()
// -------------------------------------------------------------------
PentCase
tCQMCompQuad::classifyPentagon(
                               tCQMPolygon2* p,
                               tCQMSpanningTreeVertex* v1,
                               tCQMSpanningTreeVertex* v2,
                               tCQMSpanningTreeVertex* v3
                               )
{
    // get  the mesh  faces corresponding  to the  dual graph  vertex in
    // "v1", "v2" and "v3".
    tCQMFace2* f1 = v1->getVertex()->getFace();
    tCQMFace2* f2 = v2->getVertex()->getFace();
    tCQMFace2* f3 = v3->getVertex()->getFace();
    
    // get the common edge of "f1" and "f2"
    tCQMEdge2* e1 = f1->getCommonEdge(f2);
    
    // get the common edge of "f2" and "f3"
    tCQMEdge2* e2 = f2->getCommonEdge(f3);
    
    // get the half-edge of "e1" in "f1"
    tCQMHalfEdge2* he1 = e1->getHalfEdge();
    tCQMHalfEdge2* he2 = e1->getMate(he1);
    if (he1->getFace() != f1) {
        he2 = he1;
        he1 = e1->getMate(he1);
    }
    
    // get the half-edge of "e2" in "f2"
    tCQMHalfEdge2* he3 = e2->getHalfEdge();
    tCQMHalfEdge2* he4 = e2->getMate(he3);
    if (he3->getFace() != f2) {
        he4 = he3;
        he3 = e2->getMate(he3);
    }
    
    tCQMHalfEdge2* he5 = 0;
    tCQMSpanningTreeVertex* v3par = v3->getParent();
    if (v3par != 0) {
        // "v3" has a  parent, so we look for the common  edge of "v2" and
        // "v3" and the common edge of "v3" and its parent.
        tCQMFace2* f3par = v3par->getVertex()->getFace();
        tCQMEdge2* e4 = f3->getCommonEdge(f3par);
        
        // get the half-edge of "e4" in "f3"
        he5 = e4->getHalfEdge();
        if (he5->getFace() != f3) {
            he5 = e4->getMate(he5);
        }
    }
    
    // classify pentagon
    PentCase pc;
    
    if (p->isReflex(0)) {
        // vertex "p1" of the pentagon "p" is reflex
        if (p->isReflex(1)) {
            // vertices "p1" and "p2" of the pentagon "p" are reflex
            if ((he5 == 0) || (he5->getNext() != he4)) {
                // edge p2p3 is the common edge between "v3" and its parent
                pc = R1R2E24E23;
            }
            else {
                // edge p3p4 is the common edge between "v3" and its parent
                pc = R1R2E24E34;
            }
        }
        else if (p->isReflex(3)) {
            // vertices "p1" and "p4" of the pentagon "p" are reflex
            if (he2->getNext() != he3) {
                // diagonal p2p4 is the common edge of "v2" and "v3"
                if ((he5 == 0) || (he5->getNext() != he4)) {
                    // edge p2p3 is the common edge between "v3" and its parent
                    pc = R1R4E24E23;
                }
                else {
                    // edge p3p4 is the common edge between "v3" and its parent
                    pc = R1R4E24E34;
                }
            }
            else {
                // diagonal p1p3 is the common edge of "v2" and "v3"
                if ((he5 == 0) || (he5->getNext() == he4)) {
                    // edge p2p3 is the common edge between "v3" and its parent
                    pc = R1R4E13E23;
                }
                else {
                    // edge p1p2 is the common edge between "v3" and its parent
                    pc = R1R4E13E12;
                }
            }
        }
        else {
            // only vertex "p1" of the pentagon "p" is reflex
            if (he2->getNext() != he3) {
                // diagonal p2p4 is the common edge of "v2" and "v3"
                if ((he5 == 0) || (he5->getNext() == he4)) {
                    // edge p3p4 is the common edge between "v3" and its parent
                    pc = R1E24E34;
                }
                else {
                    // edge p2p3 is the common edge between "v3" and its parent
                    pc = R1E24E23;
                }
            }
            else {
                // diagonal p1p3 is the common edge of "v2" and "v3"
                if ((he5 == 0) || (he5->getNext() == he4)) {
                    // edge p2p3 is the common edge between "v3" and its parent
                    pc = R1E13E23;
                }
                else {
                    // edge p1p2 is the common edge between "v3" and its parent
                    pc = R1E13E12;
                }
            }
        }
    }
    else {
        // vertex "p4"  of the pentagon  "p" is reflex  and "p1" is  not a
        // reflex vertex.
        if (p->isReflex(2)) {
            // vertices "p4" and "p3" are reflex
            if ((he5 == 0) || (he5->getNext() == he4)) {
                // edge p2p3 is the common edge between "v3" and its parent
                pc = R3R4E13E23;
            }
            else {
                // edge p1p2 is the common edge between "v3" and its parent
                pc = R3R4E13E12;
            }
        }
        else {
            // only vertex "p4" of the pentagon "p" is a reflex vertex
            if (he2->getNext() == he3) {
                // diagonal p1p3 is the common edge of "v2" and "v3"
                if ((he5 == 0) || (he5->getNext() == he4)) {
                    // edge p2p3 is the common edge between "v3" and its parent
                    pc = R4E13E23;
                }
                else {
                    // edge p1p2 is the common edge between "v3" and its parent
                    pc = R4E13E12;
                }
            }
            else {
                // diagonal p2p4 is the common edge of "v2" and "v4"
                if ((he5 == 0) || (he5->getNext() != he4)) {
                    // edge p2p3 is the common edge between "v3" and its parent
                    pc = R4E24E23;
                }
                else {
                    // edge p3p4 is the common edge between "v3" and its parent
                    pc = R4E24E34;
                }
            }
        }
    }
    
    return pc;
}

// -------------------------------------------------------------------
// Method classifyHexagon()
// -------------------------------------------------------------------
HexCase
tCQMCompQuad::classifyHexagon(tCQMPolygon2* p)
{
    // compute number of reflex vertices of "p"
    int nre = p->getNumOfReflexPoints();
    
    if (nre == 1) {
        int rv = p->findReflexPoint();
        if (rv != 0) p->shiftLeft(rv);
        return RCCCCC;
    }
    else if (nre == 2) {
        int rv = p->findReflexPoint();
        if (rv != 0) p->shiftLeft(rv);
        if (p->isReflex(3)) return RCCRCC;
        if (p->isReflex(2)) return RCRCCC;
        else if (p->isReflex(4)) {
            p->shiftRight(2);
            return RCRCCC;
        }
        if (p->isReflex(1)) return RRCCCC;
        p->shiftRight(1);
        return RRCCCC;
    }
    else if (nre == 3) {
        int rv = p->findReflexPoint();
        if (rv != 0) p->shiftLeft(rv);
        if (p->isReflex(1)) {
            if (p->isReflex(2)) return RRRCCC;
            if (p->isReflex(3)) return RRCRCC;
            if (p->isReflex(4)) return RRCCRC;
            p->shiftRight(1);
            return RRRCCC;
        }
        else {
            if (p->isReflex(2)) {
                if (p->isReflex(4)) {
                    return RCRCRC;
                }
                if (p->isReflex(3)) {
                    p->shiftLeft(2);
                    return RRCCRC;
                }
                p->shiftRight(1);
                return RRCRCC;
            }
            if (p->isReflex(3)) {
                if (p->isReflex(4)) {
                    p->shiftLeft(3);
                    return RRCRCC;
                }
                p->shiftRight(1);
                return RRCCRC;
            }
            p->shiftRight(2);
            return RRRCCC;
        }
    }
    
    return CCCCCC;
}

// -------------------------------------------------------------------
// Method makeQuadIfConvex()
// -------------------------------------------------------------------
bool
tCQMCompQuad::makeQuadIfConvex(tCQMSpanningTreeVertex* v1,
                               tCQMSpanningTreeVertex* v2,
                               QuadVertexSet& vs,
                               std::list<tCQMQuadrilateral2*>* lq
                               )
{
    // get the mesh faces corresponding to the dual graph vertex in "v1"
    // and "v2".
    tCQMFace2* f1 = v1->getVertex()->getFace();
    tCQMFace2* f2 = v2->getVertex()->getFace();
    
    // get the common edge of "f1" and "f2"
    tCQMEdge2* e1 = f1->getCommonEdge(f2);
    
    // get the half-edge of e1 in f1
    tCQMHalfEdge2* he1 = e1->getHalfEdge();
    tCQMHalfEdge2* he2 = e1->getMate(he1);
    if (he1->getFace() != f1) {
        he2 = he1;
        he1 = e1->getMate(he1);
    }
    
    // get vertices of the quadrilateral
    tCQMIndVertex2 p1 = getQuadVertex(vs, he1->getVertex());
    bool ec1 = he2->getNext()->getEdge()->isConstrained();
    tCQMIndVertex2 p2 = getQuadVertex(vs, he2->getPrev()->getVertex());
    bool ec2 = he2->getPrev()->getEdge()->isConstrained();
    tCQMIndVertex2 p3 = getQuadVertex(vs, he2->getVertex());
    
    tCQMIndVertex2 p4;
    bool ec3, ec4;
    if (v1->hasSteinerPoint()) {
        p4 = v1->getSteinerPoint();
        ec3 = false;
        ec4 = false;
    }
    else {
        p4 = getQuadVertex(vs, he1->getPrev()->getVertex());
        ec3 = he1->getNext()->getEdge()->isConstrained();
        ec4 = he1->getPrev()->getEdge()->isConstrained();
    }
    
    tCQMQuadrilateral2* quad = (tCQMQuadrilateral2*)
    new tCQMQuadrilateral2();
    
    // insert points in it
    quad->insert(0,p1);
    quad->insert(1,p2);
    quad->insert(2,p3);
    quad->insert(3,p4);
    
    quad->setConstraint(0,ec1);
    quad->setConstraint(1,ec2);
    quad->setConstraint(2,ec3);
    quad->setConstraint(3,ec4);
    
    // if the  quadrilareal convex,  put it into  the list of  quads and
    // return true. Otherwise, return false.
    if (quad->isConvex()) {
        // insert new convex quadrilateral into the list of quadrilaterals
        lq->push_back(quad);
        return true;
    }
    
    delete quad;
    return false;
}

// -------------------------------------------------------------------
// Method quadrangulateDegTriangle()
//
// This method  convex quadrangulates  a triangle with  vertices "p1",
// "p2" and "p3"  by creating two quadrilaterals that  share the edges
// p1sp and  spp4, where "sp" is a  Steiner point and "p4"  is a given
// point inside triangle p1p2p3.
// -------------------------------------------------------------------
tCQMIndVertex2
tCQMCompQuad::quadrangulateDegTriangle(const tCQMIndVertex2& p1,
                                       const tCQMIndVertex2& p2,
                                       const tCQMIndVertex2& p3,
                                       const tCQMIndVertex2& p4,
                                       bool ec1,
                                       bool ec3,
                                       std::list<tCQMQuadrilateral2*>* lq
                                       )
{
    // compute the middle point of segment p1p4
    double x = (p1.getPoint().x() + p4.getPoint().x()) / 2;
    double y = (p1.getPoint().y() + p4.getPoint().y()) / 2;
    
    // create Steiner point
    tCQMIndVertex2 sp = createSteinerPoint(Point2(x,y));
    
    // create two quadrilaterals
    try {
        createQuad(p1,p2,p4,sp,ec1,false,false,false,lq);
        createQuad(p1,sp,p4,p3,false,false,false,ec3,lq);
    }
    catch (const tExceptionObject& xpt) {
        treatException(xpt);
        exit(0);
    }
    
    return p4;
}

// -------------------------------------------------------------------
// Method quadrangulateDegPentagon1()
//
// This method  convex quadrangulates  a pentagon with  vertices "p1",
// "p2", "p3", "p4" and "p5" such  that p3p4 is the edge shared by the
// root of the  subtree corresponding to this pentagon  and the parent
// of this root node in the spanning tree.
// -------------------------------------------------------------------
tCQMIndVertex2
tCQMCompQuad::quadrangulateDegPentagon1(const tCQMIndVertex2& p1,
                                        const tCQMIndVertex2& p2,
                                        const tCQMIndVertex2& p3,
                                        const tCQMIndVertex2& p4,
                                        const tCQMIndVertex2& p5,
                                        bool ec1,
                                        bool ec2,
                                        bool ec3,
                                        bool ec4,
                                        bool ec5,
                                        std::list<tCQMQuadrilateral2*>* lq
                                        )
{
    // we place only one Steiner point in the triangle defined by p1p3p4
    // such that this point is visible from "p5" and "p2".
    std::list<Point2> li;
    std::list<Point2> lo;
    
    // compute intersection of triangle p1p3p4 and the wedge of "p5"
    li.push_back(p1.getPoint()); li.push_back(p3.getPoint());
    li.push_back(p4.getPoint());
    
    try {
        halfSpaceIntersection(li,p5.getPoint(),p1.getPoint(),lo,
                              CGAL::ON_NEGATIVE_SIDE);
        li.clear();
        halfSpaceIntersection(lo,p5.getPoint(),p4.getPoint(),li,
                              CGAL::ON_POSITIVE_SIDE);
        lo.clear();
        
        // compute intersection of triangle p1p3p4 and the wedge of "p2"
        halfSpaceIntersection(li,p2.getPoint(),p1.getPoint(),lo,
                              CGAL::ON_POSITIVE_SIDE);
        li.clear();
        halfSpaceIntersection(lo,p2.getPoint(),p3.getPoint(),li,
                              CGAL::ON_NEGATIVE_SIDE);
        lo.clear();
    }
    catch (const tExceptionObject& xpt) {
        treatException(xpt);
        exit(0);
    }
    
    // compute the barycenter of the resulting polygon
    Point2 q;
    double x = 0;
    double y = 0;
    int num = 0;
    while (!li.empty()) {
        q = li.front();
        li.pop_front();
        x += q.x();
        y += q.y();
        ++num;
    }
    
    q = Point2(x / double(num), y / double(num));
    
    // create Steiner point
    tCQMIndVertex2 sp = createSteinerPoint(q);
    
    // create two quadrilaterals
    try {
        createQuad(p5,p1,sp,p4,ec5,false,false,ec4,lq);
        createQuad(p1,p2,p3,sp,ec1,ec2,false,false,lq);
    }
    catch (const tExceptionObject& xpt) {
        treatException(xpt);
        exit(0);
    }
    
    return sp;
}

// -------------------------------------------------------------------
// Method quadrangulatePentR1E24E34()
//
// This method  convex quadrangulates  a pentagon with  vertices "p1",
// "p2", "p3",  "p4" and "p5" such  that only "p1" is  a reflex vertex
// and p3p4  is the  common edge of  this pentagon with  some triangle
// corresponding to the parent node  of the subtree rooted at the node
// associated with the dual vertex of p2p3p4.
// -------------------------------------------------------------------
void
tCQMCompQuad::quadrangulatePentR1E24E34(
                                        tCQMPolygon2* p,
                                        tCQMSpanningTreeVertex* v,
                                        QuadVertexSet& vs,
                                        std::list<tCQMQuadrilateral2*>* lq
                                        )
{
    // get pentagon vertices
    tCQMIndVertex2 p1 = p->getIndVertex(0);
    tCQMIndVertex2 p2 = p->getIndVertex(1);
    tCQMIndVertex2 p3 = p->getIndVertex(2);
    tCQMIndVertex2 p4 = p->getIndVertex(3);
    tCQMIndVertex2 p5 = p->getIndVertex(4);
    
    // we need only  one Steiner point to convex  quadrangulate "p", and
    // this point is placed in  a region resulting from the intersection
    // of  the  triangle p1p2p4,  the  left  half-space  defined by  the
    // supporting  line p5p1,  and the  left half-space  defined  by the
    // supporting line p1p3.
    
    // compute  intersection of supporting  line p5p1  and edge  p2p4 of
    // triangle p1p2p4.
    
    Line2 l1(p5.getPoint(),p1.getPoint());
    Segment2 s1(p2.getPoint(),p4.getPoint());
    
    // compute intersection point of "l1" and "s1"
    Point2 pt1;
    if (!myAssign(pt1,s1,l1)) {
        pt1 = p2.getPoint();
    }
    
    // compute intersection of supporting line p1p3 and triangle p1pt1p4
    Line2 l2(p1.getPoint(),p3.getPoint());
    Segment2 s2(pt1,p4.getPoint());
    
    // compute intersection point of "l2" and "s2"
    Point2 pt2;
    if (!myAssign(pt2,s2,l2)) {
        pt2 = pt1;
    }
    
    // compute the barycenter of triangle p1pt2p4
    double x = (p4.getPoint().x() + p1.getPoint().x() + pt2.x()) / 3;
    double y = (p4.getPoint().y() + p1.getPoint().y() + pt2.y()) / 3;
    
    // create Steiner point
    tCQMIndVertex2 sp = createSteinerPoint(Point2(x,y));
    
    // get edge status
    bool ec1 = p->isConstrained(0);
    bool ec2 = p->isConstrained(1);
    bool ec4 = p->isConstrained(3);
    bool ec5 = p->isConstrained(4);
    
    // create two quadrilaterals
    try {
        createQuad(p5,p1,sp,p4,ec5,false,false,ec4,lq);
        createQuad(p1,p2,p3,sp,ec1,ec2,false,false,lq);
    }
    catch (const tExceptionObject& xpt) {
        treatException(xpt);
        exit(0);
    }
    
    v->setSteinerPoint(sp,p3.getPoint());
}

// -------------------------------------------------------------------
// Method quadrangulatePentR4E13E12()
//
// This method  convex quadrangulates  a pentagon with  vertices "p1",
// "p2", "p3",  "p4" and "p5" such  that only "p4" is  a reflex vertex
// and p1p2  is the  common edge of  this pentagon with  some triangle
// corresponding to the parent node  of the subtree rooted at the node
// associated with the dual vertex of p2p3p4.
// -------------------------------------------------------------------
void
tCQMCompQuad::quadrangulatePentR4E13E12(
                                        tCQMPolygon2* p,
                                        tCQMSpanningTreeVertex* v,
                                        QuadVertexSet& vs,
                                        std::list<tCQMQuadrilateral2*>* lq
                                        )
{
    // get pentagon vertices
    tCQMIndVertex2 p1 = p->getIndVertex(0);
    tCQMIndVertex2 p2 = p->getIndVertex(1);
    tCQMIndVertex2 p3 = p->getIndVertex(2);
    tCQMIndVertex2 p4 = p->getIndVertex(3);
    tCQMIndVertex2 p5 = p->getIndVertex(4);
    
    // we need only  one Steiner point to convex  quadrangulate "p", and
    // this point is placed in  a region resulting from the intersection
    // of  the triangle  p1p3p4,  the right  half-space  defined by  the
    // supporting  line p5p4, and  the right  half-space defined  by the
    // supporting line p4p2.
    
    // compute  intersection of supporting  line p5p4  and edge  p1p3 of
    // triangle p1p3p4.
    
    Line2 l1(p5.getPoint(),p4.getPoint());
    Segment2 s1(p3.getPoint(),p1.getPoint());
    
    // compute intersection point of "l1" and "s"
    Point2 pt1;
    if (!myAssign(pt1,s1,l1)) {
        pt1 = p3.getPoint();
    }
    
    // compute intersection of supporting line p4p2 and triangle p1pt1p4
    Line2 l2(p4.getPoint(),p2.getPoint());
    Segment2 s2(pt1,p1.getPoint());
    
    // compute intersection point of "l2" and "s2"
    Point2 pt2;
    if (!myAssign(pt2,s2,l2)) {
        pt2 = pt1;
    }
    
    // compute the barycenter of triangle p1pt2p4
    double x = (p1.getPoint().x() + p4.getPoint().x() + pt2.x()) / 3;
    double y = (p1.getPoint().y() + p4.getPoint().y() + pt2.y()) / 3;
    
    // create Steiner point
    tCQMIndVertex2 sp = createSteinerPoint(Point2(x,y));
    
    // get edge status
    bool ec2 = p->isConstrained(1);
    bool ec3 = p->isConstrained(2);
    bool ec4 = p->isConstrained(3);
    bool ec5 = p->isConstrained(4);
    
    // create two quadrilaterals
    try {
        createQuad(p5,p1,sp,p4,ec5,false,false,ec4,lq);
        createQuad(p2,p3,p4,sp,ec2,ec3,false,false,lq);
    }
    catch (const tExceptionObject& xpt) {
        treatException(xpt);
        exit(0);
    }
    
    v->setSteinerPoint(sp,p1.getPoint());
}

// -------------------------------------------------------------------
// Method quadrangulatePentR1E24E23()
//
// This method  convex quadrangulates  a pentagon with  vertices "p1",
// "p2", "p3",  "p4" and "p5" such  that only "p1" is  a reflex vertex
// and p2p3  is the  common edge of  this pentagon with  some triangle
// corresponding to the parent node  of the subtree rooted at the node
// associated with the dual vertex of p2p3p4.
// -------------------------------------------------------------------
void
tCQMCompQuad::quadrangulatePentR1E24E23(
                                        tCQMPolygon2* p,
                                        tCQMSpanningTreeVertex* v,
                                        QuadVertexSet& vs,
                                        std::list<tCQMQuadrilateral2*>* lq
                                        )
{
    
    // get pentagon vertices
    tCQMIndVertex2 p1 = p->getIndVertex(0);
    tCQMIndVertex2 p2 = p->getIndVertex(1);
    tCQMIndVertex2 p3 = p->getIndVertex(2);
    tCQMIndVertex2 p4 = p->getIndVertex(3);
    tCQMIndVertex2 p5 = p->getIndVertex(4);
    
    // we need two  Steiner points to convex quadrangulate  "p", and the
    // first one is  placed in a region resulting  from the intersection
    // of the  triangle p1p2p4  and the left  half-space defined  by the
    // supporting line p5p1.
    
    // compute  intersection of supporting  line p5p1  and edge  p2p4 of
    // triangle p1p2p4.
    
    Line2 l1(p5.getPoint(),p1.getPoint());
    Segment2 s(p2.getPoint(),p4.getPoint());
    
    // compute intersection point of "l1" and "s"
    Point2 pt1;
    if (!myAssign(pt1,s,l1)) {
        pt1 = p2.getPoint();
    }
    
    pt1 = Point2((pt1.x() + p4.getPoint().x() + p1.getPoint().x()) / 3,
                 (pt1.y() + p4.getPoint().y() + p1.getPoint().y()) / 3);
    
    // the second  Steiner point  is placed in  the intersection  of the
    // triangle pt1p2p3, the right  half-space defined by the supporting
    // line  p1pt1, and the  left half-space  defined by  the supporting
    // line p4pt1.
    std::list<Point2> li;
    std::list<Point2> lo;
    li.push_back(pt1); li.push_back(p2.getPoint());
    li.push_back(p3.getPoint());
    
    try {
        halfSpaceIntersection(li,p1.getPoint(),pt1,lo,
                              CGAL::ON_POSITIVE_SIDE);
        li.clear();
        halfSpaceIntersection(lo,p4.getPoint(),pt1,li,
                              CGAL::ON_NEGATIVE_SIDE);
        lo.clear();
    }
    catch (const tExceptionObject& xpt) {
        treatException(xpt);
        exit(0);
    }
    
    // compute the barycenter of the resulting polygon
    Point2 pt2;
    double x = 0;
    double y = 0;
    int num = 0;
    while (!li.empty()) {
        pt2 = li.front();
        li.pop_front();
        x += pt2.x();
        y += pt2.y();
        ++num;
    }
    
    pt2 = Point2(x / double(num), y / double(num));
    
    // create Steiner points
    tCQMIndVertex2 sp1 = createSteinerPoint(pt1);
    tCQMIndVertex2 sp2 = createSteinerPoint(pt2);
    
    // get edge status
    bool ec1 = p->isConstrained(0);
    bool ec3 = p->isConstrained(2);
    bool ec4 = p->isConstrained(3);
    bool ec5 = p->isConstrained(4);
    
    // create three quadrilaterals
    try {
        createQuad( p5, p1,sp1, p4,ec5,false,false,ec4,lq);
        createQuad( p1, p2,sp2,sp1,ec1,false,false,false,lq);
        createQuad(sp1,sp2, p3, p4,false,false,ec3,false,lq);
    }
    catch (const tExceptionObject& xpt) {
        treatException(xpt);
        exit(0);
    }
    
    v->setSteinerPoint(sp2,p2.getPoint());
}

// -------------------------------------------------------------------
// Method quadrangulatePentR4E13E23()
//
// This method  convex quadrangulates  a pentagon with  vertices "p1",
// "p2", "p3",  "p4" and "p5" such  that only "p4" is  a reflex vertex
// and p2p3  is the  common edge of  this pentagon with  some triangle
// corresponding to the parent node  of the subtree rooted at the node
// associated with the dual vertex of p1p3p4.
// -------------------------------------------------------------------
void
tCQMCompQuad::quadrangulatePentR4E13E23(
                                        tCQMPolygon2* p,
                                        tCQMSpanningTreeVertex* v,
                                        QuadVertexSet& vs,
                                        std::list<tCQMQuadrilateral2*>* lq
                                        )
{
    // get pentagon vertices
    tCQMIndVertex2 p1 = p->getIndVertex(0);
    tCQMIndVertex2 p2 = p->getIndVertex(1);
    tCQMIndVertex2 p3 = p->getIndVertex(2);
    tCQMIndVertex2 p4 = p->getIndVertex(3);
    tCQMIndVertex2 p5 = p->getIndVertex(4);
    
    // we need two  Steiner points to convex quadrangulate  "p", and the
    // first one is  placed in a region resulting  from the intersection
    // of the  triangle p1p3p4 and  the right half-space defined  by the
    // supporting line p5p4.
    
    // compute  intersection of supporting  line p5p4  and edge  p3p1 of
    // triangle p1p3p4.
    
    Line2 l1(p5.getPoint(),p4.getPoint());
    Segment2 s(p3.getPoint(),p1.getPoint());
    
    // compute intersection point of "l1" and "s"
    Point2 pt1;
    if (!myAssign(pt1,s,l1)) {
        pt1 = p3.getPoint();
    }
    
    pt1 = Point2((pt1.x() + p4.getPoint().x() + p1.getPoint().x()) / 3,
                 (pt1.y() + p4.getPoint().y() + p1.getPoint().y()) / 3);
    
    // the second  Steiner point  is placed in  the intersection  of the
    // triangle pt1p2p3,  the left half-space defined  by the supporting
    // line p4pt1,  and the right  half-space defined by  the supporting
    // line p1pt1.
    std::list<Point2> li;
    std::list<Point2> lo;
    li.push_back(pt1); li.push_back(p2.getPoint());
    li.push_back(p3.getPoint());
    
    try {
        halfSpaceIntersection(li,p4.getPoint(),pt1,lo,
                              CGAL::ON_NEGATIVE_SIDE);
        li.clear();
        halfSpaceIntersection(lo,p1.getPoint(),pt1,li,
                              CGAL::ON_POSITIVE_SIDE);
        lo.clear();
    }
    catch (const tExceptionObject& xpt) {
        treatException(xpt);
        exit(0);
    }
    
    // compute the barycenter of the resulting polygon
    Point2 pt2;
    double x = 0;
    double y = 0;
    int num = 0;
    while (!li.empty()) {
        pt2 = li.front();
        li.pop_front();
        x += pt2.x();
        y += pt2.y();
        ++num;
    }
    
    pt2 = Point2(x / double(num), y / double(num));
    
    // create Steiner points
    tCQMIndVertex2 sp1 = createSteinerPoint(pt1);
    tCQMIndVertex2 sp2 = createSteinerPoint(pt2);
    
    // get edge status
    bool ec1 = p->isConstrained(0);
    bool ec3 = p->isConstrained(2);
    bool ec4 = p->isConstrained(3);
    bool ec5 = p->isConstrained(4);
    
    // create three quadrilaterals
    try {
        createQuad( p5, p1,sp1, p4,ec5,false,false,ec4,lq);
        createQuad( p1, p2,sp2,sp1,ec1,false,false,false,lq);
        createQuad(sp1,sp2, p3, p4,false,false,ec3,false,lq);
    }
    catch (const tExceptionObject& xpt) {
        treatException(xpt);
        exit(0);
    }
    
    v->setSteinerPoint(sp2,p4.getPoint());
}

// -------------------------------------------------------------------
// Method quadrangulatePentR1R2E24E34()
//
// This method  convex quadrangulates  a pentagon with  vertices "p1",
// "p2",  "p3", "p4"  and  "p5" such  that  "p1" and  "p2" are  reflex
// vertices and  p3p4 is  the common edge  of this pentagon  with some
// triangle corresponding to the parent  node of the subtree rooted at
// the node associated with the dual vertex of p2p3p4.
// -------------------------------------------------------------------
void
tCQMCompQuad::quadrangulatePentR1R2E24E34(
                                          tCQMPolygon2* p,
                                          tCQMSpanningTreeVertex* v,
                                          QuadVertexSet& vs,
                                          std::list<tCQMQuadrilateral2*>* lq
                                          )
throw (tExceptionObject)
{
    // get pentagon vertices
    tCQMIndVertex2 p1 = p->getIndVertex(0);
    tCQMIndVertex2 p2 = p->getIndVertex(1);
    tCQMIndVertex2 p3 = p->getIndVertex(2);
    tCQMIndVertex2 p4 = p->getIndVertex(3);
    tCQMIndVertex2 p5 = p->getIndVertex(4);
    
    // we need two  Steiner points to convex quadrangulate  "p", and the
    // first one is  placed in a region resulting  from the intersection
    // of  the  triangle p1p2p4,  the  left  half-space  defined by  the
    // supporting  line p5p1, and  the right  half-space defined  by the
    // supporting line p3p2.
    
    // compute intersection of supporting line p5p1 and diagonal p2p4 of
    // triangle p1p2p4.
    
    Line2 l1(p5.getPoint(),p1.getPoint());
    Segment2 s1(p2.getPoint(),p4.getPoint());
    
    // compute intersection point of "l1" and "s1"
    double x;
    double y;
    Point2 pt1;
    if (!myAssign(pt1,s1,l1)) {
        pt1 = p2.getPoint();
    }
    
    if (IS_EQUAL(pt1.x(),p2.getPoint().x()) &&
        IS_EQUAL(pt1.y(),p2.getPoint().y())) {
        // compute intersection of supporting  line p3p2 and diagonal p1p4
        // of triangle p1p2p4.
        Line2 l2(p3.getPoint(),p2.getPoint());
        Segment2 s2(p1.getPoint(),p4.getPoint());
        
        Point2 pt2;
        if (!myAssign(pt2,s2,l2)) {
            pt2 = p1.getPoint();
        }
        
        x = (pt1.x() + p4.getPoint().x() + pt2.x()) / 3;
        y = (pt1.y() + p4.getPoint().y() + pt2.y()) / 3;
    }
    else {
        // compute intersection of supporting  line p3p2 and diagonal p1p4
        // of triangle p1p2p4.
        Line2 l2(p3.getPoint(),p2.getPoint());
        Segment2 s2(p1.getPoint(),p4.getPoint());
        
        Point2 pt2;
        if (!myAssign(pt2,s2,l2)) {
            pt2 = p1.getPoint();
        }
        
        if (IS_EQUAL(pt2.x(),p1.getPoint().x()) &&
            IS_EQUAL(pt2.y(),p1.getPoint().y())) {
            x = (pt1.x() + p4.getPoint().x() + pt2.x()) / 3;
            y = (pt1.y() + p4.getPoint().y() + pt2.y()) / 3;
        }
        else {
            // compute  intersection of  supporting line  p3p2  and diagonal
            // p1pt1.
            Segment2 s3(p1.getPoint(),pt1);
            Point2 pt3;
            if (!myAssign(pt3,s3,l2)) {
                if (!segLineIntersection(p3.getPoint(),p2.getPoint(),p1.getPoint(),pt1,pt3)) {
                    std::stringstream ss (std::stringstream::in | std::stringstream::out);
                    ss << "R1R2E24E34(): numerical instability on calculating line intersection";
                    throw tExceptionObject(__FILE__,__LINE__,ss.str());
                }
            }
            
            x = (pt1.x() + p4.getPoint().x() + pt2.x() + pt3.x()) / 4;
            y = (pt1.y() + p4.getPoint().y() + pt2.y() + pt3.y()) / 4;
        }
    }
    
    // create location of the first Steiner point
    Point2 q1(x,y);
    
    // Second Steiner point  is placed in the region  resulting from the
    // intersection of the triangle p2p3p4, the right half-space defined
    // by the  supporting line p5p1  and the left half-space  defined by
    // the supporting line p1p2.
    
    Segment2 s2(p3.getPoint(),p4.getPoint());
    
    // compute intersection point of "l1" and "s2"
    Point2 pt2;
    if (!myAssign(pt2,s2,l1)) {
        pt2 = p3.getPoint();
    }
    
    Point2 q2, q3;
    if (IS_EQUAL(pt2.x(),p3.getPoint().x()) &&
        IS_EQUAL(pt2.y(),p3.getPoint().y())) {
        // compute intersection of supporting line p1p2 and edge p3p4
        Line2 l2(p1.getPoint(),p2.getPoint());
        Point2 pt3;
        if (!myAssign(pt3,s2,l2)) {
            pt3 = p3.getPoint();
        }
        
        // compute intersection of supporting line q1pt3 and diagonal p2p4
        Line2 l3(q1,pt3);
        Point2 pt4;
        if (!myAssign(pt4,s1,l3)) {
            std::stringstream ss (std::stringstream::in | std::stringstream::out);
            ss << "R1R2E24E34(): numerical instability on calculating line intersection";
            throw tExceptionObject(__FILE__,__LINE__,ss.str());
        }
        
        x = (pt4.x() + pt3.x() + p2.getPoint().x()) / 3;
        y = (pt4.y() + pt3.y() + p2.getPoint().y()) / 3;
        
        q2 = Point2(x,y);
        
        // compute intersection of supporting line p2q2 and edge p3p4
        Line2 l4(p2.getPoint(),q2);
        Point2 pt5;
        if (!myAssign(pt5,s2,l4)) {
            std::stringstream ss (std::stringstream::in | std::stringstream::out);
            ss << "R1R2E24E34(): numerical instability on calculating line intersection";
            throw tExceptionObject(__FILE__,__LINE__,ss.str());
        }
        
        x = (pt3.x() + pt5.x()) / 2;
        y = (pt3.y() + pt5.y()) / 2;
        
        // third Steiner  point is  placed in the  segment defined  by the
        // point (x,y) and "q2".
        x = x + 0.3 * (q2.x() - x);
        y = y + 0.3 * (q2.y() - y);
        
        q3 = Point2(x,y);
    }
    else {
        // compute intersection point of "l2" and "s2"
        Line2 l2(p1.getPoint(),p2.getPoint());
        
        Point2 pt3;
        if (!myAssign(pt3,s2,l2)) {
            pt3 = p3.getPoint();
        }
        
        x = (p2.getPoint().x() + pt3.x() + pt2.x() + pt1.x()) / 4;
        y = (p2.getPoint().y() + pt3.y() + pt2.y() + pt1.y()) / 4;
        
        q2 = Point2(x,y);
        
        // compute intersection of supporting line p2q2 and edge p3p4
        Line2 l3(p2.getPoint(),q2);
        Point2 pt4;
        if (!myAssign(pt4,s2,l3)) {
            std::stringstream ss (std::stringstream::in | std::stringstream::out);
            ss << "R1R2E24E34(): numerical instability on calculating line intersection";
            throw tExceptionObject(__FILE__,__LINE__,ss.str());
        }
        
        x = (pt3.x() + pt4.x()) / 2;
        y = (pt3.y() + pt4.y()) / 2;
        
        // third Steiner  point is  placed in the  segment defined  by the
        // point (x,y) and "q2".
        x = x + 0.3 * (q2.x() - x);
        y = y + 0.3 * (q2.y() - y);
        
        q3 = Point2(x,y);
    }
    
    // create Steiner points
    tCQMIndVertex2 sp1 = createSteinerPoint(q1);
    tCQMIndVertex2 sp2 = createSteinerPoint(q2);
    tCQMIndVertex2 sp3 = createSteinerPoint(q3);
    
    // get edge status
    bool ec1 = p->isConstrained(0);
    bool ec2 = p->isConstrained(1);
    bool ec4 = p->isConstrained(3);
    bool ec5 = p->isConstrained(4);
    
    // create four quadrilaterals
    try {
        createQuad( p5, p1,sp1, p4,ec5,false,false,ec4,lq);
        createQuad( p1, p2,sp2,sp1,ec1,false,false,false,lq);
        createQuad(sp1,sp2,sp3, p4,false,false,false,false,lq);
        createQuad(sp3,sp2, p2, p3,false,false,ec2,false,lq);
    }
    catch (const tExceptionObject& xpt) {
        treatException(xpt);
        exit(0);
    }
    
    // get the mesh faces corresponding to the dual graph vertex in "v3"
    // and its parent.
    tCQMFace2* f1 = v->getVertex()->getFace();
    tCQMFace2* f2 = v->getParent()->getVertex()->getFace();
    
    // get the common edge of "f1" and "f2"
    tCQMEdge2* e = f1->getCommonEdge(f2);
    
    // get the half-edge of "e" in "f2"
    tCQMHalfEdge2* he = e->getHalfEdge();
    if (he->getFace() != f2) {
        he = e->getMate(he);
    }
    
    v->getParent()->setQuadSteinerPoint(sp3,he);
}

// -------------------------------------------------------------------
// Method quadrangulatePentR3R4E13E12()
//
// This method  convex quadrangulates  a pentagon with  vertices "p1",
// "p2",  "p3", "p4"  and  "p5" such  that  "p3" and  "p4" are  reflex
// vertices and  p1p2 is  the common edge  of this pentagon  with some
// triangle corresponding to the parent  node of the subtree rooted at
// the node associated with the dual vertex of p1p3p4.
// -------------------------------------------------------------------
void
tCQMCompQuad::quadrangulatePentR3R4E13E12(
                                          tCQMPolygon2* p,
                                          tCQMSpanningTreeVertex* v,
                                          QuadVertexSet& vs,
                                          std::list<tCQMQuadrilateral2*>* lq
                                          )
throw (tExceptionObject)
{
    // get pentagon vertices
    tCQMIndVertex2 p1 = p->getIndVertex(0);
    tCQMIndVertex2 p2 = p->getIndVertex(1);
    tCQMIndVertex2 p3 = p->getIndVertex(2);
    tCQMIndVertex2 p4 = p->getIndVertex(3);
    tCQMIndVertex2 p5 = p->getIndVertex(4);
    
    // we need two  Steiner points to convex quadrangulate  "p", and the
    // first one is  placed in a region resulting  from the intersection
    // of the  triangle p1p3p4 and  the right half-space defined  by the
    // supporting line p5p4.
    
    // compute  intersection of supporting  line p5p4  and edge  p3p1 of
    // triangle p1p3p4.
    
    Line2 l1(p5.getPoint(),p4.getPoint());
    Segment2 s1(p3.getPoint(),p1.getPoint());
    
    // compute intersection point of "l1" and "s1"
    double x;
    double y;
    Point2 pt1;
    if (!myAssign(pt1,s1,l1)) {
        pt1 = p3.getPoint();
    }
    
    if (IS_EQUAL(pt1.x(),p3.getPoint().x()) &&
        IS_EQUAL(pt1.y(),p3.getPoint().y())) {
        // compute intersection of supporting  line p2p3 and diagonal p4p1
        Line2 l2(p2.getPoint(),p3.getPoint());
        Segment2 s2(p4.getPoint(),p1.getPoint());
        
        Point2 pt2;
        if (!myAssign(pt2,s2,l2)) {
            pt2 = p4.getPoint();
        }
        
        x = (pt1.x() + p1.getPoint().x() + pt2.x()) / 3;
        y = (pt1.y() + p1.getPoint().y() + pt2.y()) / 3;
    }
    else {
        // compute intersection of supporting  line p2p3 and diagonal p4p1
        Line2 l2(p2.getPoint(),p3.getPoint());
        Segment2 s2(p4.getPoint(),p1.getPoint());
        
        Point2 pt2;
        if (!myAssign(pt2,s2,l2)) {
            pt2 = p4.getPoint();
        }
        
        if (IS_EQUAL(pt2.x(),p4.getPoint().x()) &&
            IS_EQUAL(pt2.y(),p4.getPoint().y())) {
            x = (pt1.x() + p1.getPoint().x() + pt2.x()) / 3;
            y = (pt1.y() + p1.getPoint().y() + pt2.y()) / 3;
        }
        else {
            // compute  intersection of  supporting line  p2p3  and diagonal
            // p4pt1.
            Segment2 s3(p4.getPoint(),pt1);
            Point2 pt3;
            if (!myAssign(pt3,s3,l2)) {
                std::stringstream ss (std::stringstream::in | std::stringstream::out);
                ss << "R3R4E13E12(): numerical instability on calculating line intersection";
                throw tExceptionObject(__FILE__,__LINE__,ss.str());
            }
            
            x = (pt1.x() + pt3.x() + pt2.x() + p1.getPoint().x()) / 4;
            y = (pt1.y() + pt3.y() + pt2.y() + p1.getPoint().y()) / 4;
        }
    }
    
    // create location of the first Steiner point
    Point2 q1(x,y);
    
    // Second Steiner point  is placed in the region  resulting from the
    // intersection of the triangle  p1p2p3, the left half-space defined
    // by the supporting  line p5p4 and the right  half-space defined by
    // the supporting line p4p3.
    
    Segment2 s2(p2.getPoint(),p1.getPoint());
    
    Point2 pt2;
    if (!myAssign(pt2,s2,l1)) {
        pt2 = p2.getPoint();
    }
    
    Point2 q2, q3;
    if (IS_EQUAL(pt2.x(),p2.getPoint().x()) &&
        IS_EQUAL(pt2.y(),p2.getPoint().y())) {
        // compute intersection of supporting line p4p3 and edge p2p1
        Line2 l2(p4.getPoint(),p3.getPoint());
        Point2 pt3;
        if (!myAssign(pt3,s2,l2)) {
            pt3 = p2.getPoint();
        }
        
        // compute intersection of supporting line q1pt3 and diagonal p3p1
        Line2 l3(q1,pt3);
        Point2 pt4;
        if (!myAssign(pt4,s1,l3)) {
            std::stringstream ss (std::stringstream::in | std::stringstream::out);
            ss << "R3R4E13E12(): numerical instability on calculating line intersection";
            throw tExceptionObject(__FILE__,__LINE__,ss.str());
        }
        
        x = (pt4.x() + pt3.x() + p3.getPoint().x()) / 3;
        y = (pt4.y() + pt3.y() + p3.getPoint().y()) / 3;
        
        q2 = Point2(x,y);
        
        // compute intersection of supporting line p3q2 and edge p2p1
        Line2 l4(p3.getPoint(),q2);
        Point2 pt5;
        if (!myAssign(pt5,s2,l4)) {
            std::stringstream ss (std::stringstream::in | std::stringstream::out);
            ss << "R3R4E13E12(): numerical instability on calculating line intersection";
            throw tExceptionObject(__FILE__,__LINE__,ss.str());
        }
        
        x = (pt3.x() + pt5.x()) / 2;
        y = (pt3.y() + pt5.y()) / 2;
        
        // third Steiner  point is  placed in the  segment defined  by the
        // point (x,y) and "q2".
        x = x + 0.3 * (q2.x() - x);
        y = y + 0.3 * (q2.y() - y);
        
        q3 = Point2(x,y);
    }
    else {
        // compute intersection point of "l2" and "s2"
        Line2 l2(p4.getPoint(),p3.getPoint());
        
        Point2 pt3;
        if (!myAssign(pt3,s2,l2)) {
            pt3 = p2.getPoint();
        }
        
        x = (p3.getPoint().x() + pt1.x() + pt2.x() + pt3.x()) / 4;
        y = (p3.getPoint().y() + pt1.y() + pt2.y() + pt3.y()) / 4;
        
        q2 = Point2(x,y);
        
        // compute intersection of supporting line p3q2 and edge p2p1
        Line2 l3(p3.getPoint(),q2);
        Point2 pt4;
        if (!myAssign(pt4,s2,l3)) {
            std::stringstream ss (std::stringstream::in | std::stringstream::out);
            ss << "R3R4E13E12(): numerical instability on calculating line intersection";
            throw tExceptionObject(__FILE__,__LINE__,ss.str());
        }
        
        x = (pt3.x() + pt4.x()) / 2;
        y = (pt3.y() + pt4.y()) / 2;
        
        // third Steiner  point is  placed in the  segment defined  by the
        // point (x,y) and "q2".
        x = x + 0.3 * (q2.x() - x);
        y = y + 0.3 * (q2.y() - y);
        
        q3 = Point2(x,y);
    }
    
    // create Steiner points
    tCQMIndVertex2 sp1 = createSteinerPoint(q1);
    tCQMIndVertex2 sp2 = createSteinerPoint(q2);
    tCQMIndVertex2 sp3 = createSteinerPoint(q3);
    
    // get edge status
    bool ec2 = p->isConstrained(1);
    bool ec3 = p->isConstrained(2);
    bool ec4 = p->isConstrained(3);
    bool ec5 = p->isConstrained(4);
    
    // create four quadrilaterals
    try {
        createQuad( p5, p1,sp1, p4,ec5,false,false,ec4,lq);
        createQuad( p4,sp1,sp2, p3,false,false,false,ec3,lq);
        createQuad(sp2,sp1, p1,sp3,false,false,false,false,lq);
        createQuad( p3,sp2,sp3, p2,false,false,false,ec2,lq);
    }
    catch (const tExceptionObject& xpt) {
        treatException(xpt);
        exit(0);
    }
    
    // get the mesh faces corresponding to the dual graph vertex in "v3"
    // and its parent.
    tCQMFace2* f1 = v->getVertex()->getFace();
    tCQMFace2* f2 = v->getParent()->getVertex()->getFace();
    
    // get the common edge of "f1" and "f2"
    tCQMEdge2* e = f1->getCommonEdge(f2);
    
    // get the half-edge of "e" in "f2"
    tCQMHalfEdge2* he = e->getHalfEdge();
    if (he->getFace() != f2) {
        he = e->getMate(he);
    }
    
    v->getParent()->setQuadSteinerPoint(sp3,he);
}

// -------------------------------------------------------------------
// Method quadrangulatePentR1R2E24E23()
//
// This method  convex quadrangulates  a pentagon with  vertices "p1",
// "p2",  "p3", "p4"  and  "p5" such  that  "p1" and  "p2" are  reflex
// vertices and  p2p3 is  the common edge  of this pentagon  with some
// triangle corresponding to the parent  node of the subtree rooted at
// the node associated with the dual vertex of p2p3p4.
// -------------------------------------------------------------------
void
tCQMCompQuad::quadrangulatePentR1R2E24E23(
                                          tCQMPolygon2* p,
                                          tCQMSpanningTreeVertex* v,
                                          QuadVertexSet& vs,
                                          std::list<tCQMQuadrilateral2*>* lq
                                          )
throw (tExceptionObject)
{
    // get pentagon vertices
    tCQMIndVertex2 p1 = p->getIndVertex(0);
    tCQMIndVertex2 p2 = p->getIndVertex(1);
    tCQMIndVertex2 p3 = p->getIndVertex(2);
    tCQMIndVertex2 p4 = p->getIndVertex(3);
    tCQMIndVertex2 p5 = p->getIndVertex(4);
    
    // we need two  Steiner points to convex quadrangulate  "p", and the
    // first one is  placed in a region resulting  from the intersection
    // of the triangle p1p2p4, the wedge of "p5" and the wedge of "p3".
    
    // compute intersection of supporting line p5p1 and diagonal p2p4 of
    // triangle p1p2p4.
    
    Line2 l1(p5.getPoint(),p1.getPoint());
    Segment2 s1(p2.getPoint(),p4.getPoint());
    
    // compute intersection point of "l1" and "s1"
    double x;
    double y;
    Point2 pt1;
    if (!myAssign(pt1,s1,l1)) {
        pt1 = p2.getPoint();
    }
    
    if (IS_EQUAL(pt1.x(),p2.getPoint().x()) &&
        IS_EQUAL(pt1.y(),p2.getPoint().y())) {
        // compute intersection of supporting  line p3p2 and diagonal p1p4
        // of triangle p1p2p4.
        Line2 l2(p3.getPoint(),p2.getPoint());
        Segment2 s2(p1.getPoint(),p4.getPoint());
        
        Point2 pt2;
        if (!myAssign(pt2,s2,l2)) {
            pt2 = p1.getPoint();
        }
        
        x = (pt1.x() + p4.getPoint().x() + pt2.x()) / 3;
        y = (pt1.y() + p4.getPoint().y() + pt2.y()) / 3;
    }
    else {
        // compute intersection of supporting  line p3p2 and diagonal p1p4
        // of triangle p1p2p4.
        Line2 l2(p3.getPoint(),p2.getPoint());
        Segment2 s2(p1.getPoint(),p4.getPoint());
        
        Point2 pt2;
        if (!myAssign(pt2,s2,l2)) {
            pt2 = p1.getPoint();
        }
        
        if (IS_EQUAL(pt2.x(),p1.getPoint().x()) &&
            IS_EQUAL(pt2.y(),p1.getPoint().y())) {
            x = (pt1.x() + p4.getPoint().x() + pt2.x()) / 3;
            y = (pt1.y() + p4.getPoint().y() + pt2.y()) / 3;
        }
        else {
            // compute  intersection of  supporting line  p3p2  and diagonal
            // p1pt1.
            Segment2 s3(p1.getPoint(),pt1);
            Point2 pt3;
            if (!myAssign(pt3,s3,l2)) {
                std::stringstream ss (std::stringstream::in | std::stringstream::out);
                ss << "R1R2E24E23(): numerical instability on calculating line intersection";
                throw tExceptionObject(__FILE__,__LINE__,ss.str());
            }
            
            x = (pt1.x() + p4.getPoint().x() + pt2.x() + pt3.x()) / 4;
            y = (pt1.y() + p4.getPoint().y() + pt2.y() + pt3.y()) / 4;
        }
    }
    
    // create location of the first Steiner point
    Point2 q1(x,y);
    
    // Second Steiner point  is placed inside the triangle q1p2p3
    Point2 q2((q1.x() + p2.getPoint().x() + p3.getPoint().x()) / 3,
              (q1.y() + p2.getPoint().y() + p3.getPoint().y()) / 3);
    
    // create Steiner points
    tCQMIndVertex2 sp1 = createSteinerPoint(q1);
    tCQMIndVertex2 sp2 = createSteinerPoint(q2);
    
    // get edge status
    bool ec1 = p->isConstrained(0);
    bool ec3 = p->isConstrained(2);
    bool ec4 = p->isConstrained(3);
    bool ec5 = p->isConstrained(4);
    
    // create three quadrilaterals
    try {
        createQuad( p5, p1,sp1, p4,ec5,false,false,ec4,lq);
        createQuad( p1, p2,sp2,sp1,ec1,false,false,false,lq);
        createQuad(sp1,sp2, p3, p4,false,false,ec3,false,lq);
    }
    catch (const tExceptionObject& xpt) {
        treatException(xpt);
        exit(0);
    }
    
    v->setSteinerPoint(sp2,p2.getPoint());
}

// -------------------------------------------------------------------
// Method quadrangulatePentR3R4E13E23()
//
// This method  convex quadrangulates  a pentagon with  vertices "p1",
// "p2",  "p3", "p4"  and  "p5" such  that  "p3" and  "p4" are  reflex
// vertices and  p2p3 is  the common edge  of this pentagon  with some
// triangle corresponding to the parent  node of the subtree rooted at
// the node associated with the dual vertex of p1p2p3.
// -------------------------------------------------------------------
void
tCQMCompQuad::quadrangulatePentR3R4E13E23(
                                          tCQMPolygon2* p,
                                          tCQMSpanningTreeVertex* v,
                                          QuadVertexSet& vs,
                                          std::list<tCQMQuadrilateral2*>* lq
                                          )
throw (tExceptionObject)
{
    // get pentagon vertices
    tCQMIndVertex2 p1 = p->getIndVertex(0);
    tCQMIndVertex2 p2 = p->getIndVertex(1);
    tCQMIndVertex2 p3 = p->getIndVertex(2);
    tCQMIndVertex2 p4 = p->getIndVertex(3);
    tCQMIndVertex2 p5 = p->getIndVertex(4);
    
    // we need two  Steiner points to convex quadrangulate  "p", and the
    // first one is  placed in a region resulting  from the intersection
    // of the triangle p1p3p4, the wedge of "p5" and the wedge of "p2".
    
    // compute  intersection of supporting  line p5p4  and edge  p3p1 of
    // triangle p1p3p4.
    
    Line2 l1(p5.getPoint(),p4.getPoint());
    Segment2 s1(p3.getPoint(),p1.getPoint());
    
    // compute intersection point of "l1" and "s1"
    double x;
    double y;
    Point2 pt1;
    if (!myAssign(pt1,s1,l1)) {
        pt1 = p3.getPoint();
    }
    
    if (IS_EQUAL(pt1.x(),p3.getPoint().x()) &&
        IS_EQUAL(pt1.y(),p3.getPoint().y())) {
        // compute intersection of supporting  line p2p3 and diagonal p4p1
        Line2 l2(p2.getPoint(),p3.getPoint());
        Segment2 s2(p4.getPoint(),p1.getPoint());
        
        Point2 pt2;
        if (!myAssign(pt2,s2,l2)) {
            pt2 = p4.getPoint();
        }
        
        x = (pt1.x() + p1.getPoint().x() + pt2.x()) / 3;
        y = (pt1.y() + p1.getPoint().y() + pt2.y()) / 3;
    }
    else {
        // compute intersection of supporting  line p2p3 and diagonal p4p1
        Line2 l2(p2.getPoint(),p3.getPoint());
        Segment2 s2(p4.getPoint(),p1.getPoint());
        
        Point2 pt2;
        if (!myAssign(pt2,s2,l2)) {
            pt2 = p4.getPoint();
        }
        
        if (IS_EQUAL(pt2.x(),p4.getPoint().x()) &&
            IS_EQUAL(pt2.y(),p4.getPoint().y())) {
            x = (pt1.x() + p1.getPoint().x() + pt2.x()) / 3;
            y = (pt1.y() + p1.getPoint().y() + pt2.y()) / 3;
        }
        else {
            // compute  intersection of  supporting line  p2p3  and diagonal
            // p4pt1.
            Segment2 s3(p4.getPoint(),pt1);
            Point2 pt3;
            if (!myAssign(pt3,s3,l2)) {
                std::stringstream ss (std::stringstream::in | std::stringstream::out);
                ss << "R3R4E13E23(): numerical instability on calculating line intersection";
                throw tExceptionObject(__FILE__,__LINE__,ss.str());
            }
            
            x = (pt1.x() + pt3.x() + pt2.x() + p1.getPoint().x()) / 4;
            y = (pt1.y() + pt3.y() + pt2.y() + p1.getPoint().y()) / 4;
        }
    }
    
    // create location of the first Steiner point
    Point2 q1(x,y);
    
    // second Steiner point  is placed inside the triangle q1p2p3
    Point2 q2((q1.x() + p2.getPoint().x() + p3.getPoint().x()) / 3,
              (q1.y() + p2.getPoint().y() + p3.getPoint().y()) / 3);
    
    // create Steiner points
    tCQMIndVertex2 sp1 = createSteinerPoint(q1);
    tCQMIndVertex2 sp2 = createSteinerPoint(q2);
    
    // get edge status
    bool ec1 = p->isConstrained(0);
    bool ec3 = p->isConstrained(2);
    bool ec4 = p->isConstrained(3);
    bool ec5 = p->isConstrained(4);
    
    // create three quadrilaterals
    try {
        createQuad( p5, p1,sp1,p4,ec5,false,false,ec4,lq);
        createQuad( p4,sp1,sp2,p3,false,false,false,ec3,lq);
        createQuad( p2,sp2,sp1,p1,false,false,false,ec1,lq);
    }
    catch (const tExceptionObject& xpt) {
        treatException(xpt);
        exit(0);
    }
    
    v->setSteinerPoint(sp2,p2.getPoint());
}

// -------------------------------------------------------------------
// Method quadrangulatePentR1R4E24E34()
//
// This method  convex quadrangulates  a pentagon with  vertices "p1",
// "p2",  "p3", "p4"  and  "p5" such  that  "p1" and  "p4" are  reflex
// vertices and  p3p4 is  the common edge  of this pentagon  with some
// triangle corresponding to the parent  node of the subtree rooted at
// the node associated with the dual vertex of p2p3p4.
// -------------------------------------------------------------------
void
tCQMCompQuad::quadrangulatePentR1R4E24E34(
                                          tCQMPolygon2* p,
                                          tCQMSpanningTreeVertex* v,
                                          QuadVertexSet& vs,
                                          std::list<tCQMQuadrilateral2*>* lq
                                          )
{
    // get pentagon vertices
    tCQMIndVertex2 p1 = p->getIndVertex(0);
    tCQMIndVertex2 p2 = p->getIndVertex(1);
    tCQMIndVertex2 p3 = p->getIndVertex(2);
    tCQMIndVertex2 p4 = p->getIndVertex(3);
    tCQMIndVertex2 p5 = p->getIndVertex(4);
    
    // we need three Steiner points to convex quadrangulate "p", and the
    // first one  is placed  on the diagonal  p2p4 between "p4"  and the
    // intersection of that diagonal with the supporting line p5p1.
    Line2 l1(p5.getPoint(),p1.getPoint());
    Segment2 s1(p2.getPoint(),p4.getPoint());
    
    // compute intersection point of "l1" and "s1"
    Point2 pt1;
    if (!myAssign(pt1,s1,l1)) {
        pt1 = p2.getPoint();
    }
    
    Point2 q1((p4.getPoint().x() + pt1.x()) / 2,
              (p4.getPoint().y() + pt1.y()) / 2);
    
    // compute the middle of edge p3p4
    Point2 pt2((p4.getPoint().x() + p3.getPoint().x()) / 2,
               (p4.getPoint().y() + p3.getPoint().y()) / 2);
    
    // the second steiner  point is placed in the  region resulting from
    // the  intersection   of  the  right  half-space   defined  by  the
    // supporting  line  p1q1,  the   left  half-space  defined  by  the
    // supporting line p5p1 and the triangle pt2p4p2.
    
    // compute intersection point of "l1" and "s2"
    Point2 pt3;
    Segment2 s2(p2.getPoint(),pt2);
    if (IS_EQUAL(pt1.x(),p2.getPoint().x()) &&
        IS_EQUAL(pt1.y(),p2.getPoint().y())) {
        pt3 = pt1;
    }
    else {
        if (!myAssign(pt3,s2,l1)) {
            pt3 = p2.getPoint();
        }
    }
    
    // compute  intersection  point  of  supporting line  p1q1  and  the
    // diagonal p2pt2.
    Line2 l2(p1.getPoint(),q1);
    Point2 pt4;
    if (!myAssign(pt4,s2,l2)) {
        pt4 = p3.getPoint();
    }
    
    double x = (pt1.x() + pt3.x() + pt4.x() + q1.x()) / 4;
    double y = (pt1.y() + pt3.y() + pt4.y() + q1.y()) / 4;
    
    Point2 q2(x,y);
    
    // the third Steiner point is placed off edge p3p4
    Line2 l3(p5.getPoint(),p4.getPoint());
    Segment2 s3(q2,pt2);
    
    // compute intersection point of "l2" and "s2"
    Point2 pt5;
    if (!myAssign(pt5,s3,l3)) {
        pt5 = q2;
    }
    
    x = pt2.x() + 0.3 * (pt5.x() - pt2.x());
    y = pt2.y() + 0.3 * (pt5.y() - pt5.y());
    
    Point2 q3(x,y);
    
    // create Steiner points
    tCQMIndVertex2 sp1 = createSteinerPoint(q1);
    tCQMIndVertex2 sp2 = createSteinerPoint(q2);
    tCQMIndVertex2 sp3 = createSteinerPoint(q3);
    
    // get edge status
    bool ec1 = p->isConstrained(0);
    bool ec2 = p->isConstrained(1);
    bool ec4 = p->isConstrained(3);
    bool ec5 = p->isConstrained(4);
    
    // create four quadrilaterals
    try {
        createQuad( p5, p1,sp1, p4,ec5,false,false,ec4,lq);
        createQuad( p1, p2,sp2,sp1,ec1,false,false,false,lq);
        createQuad(sp2, p2, p3,sp3,false,ec2,false,false,lq);
        createQuad(sp2,sp3, p4,sp1,false,false,false,false,lq);
    }
    catch (const tExceptionObject& xpt) {
        treatException(xpt);
        exit(0);
    }
    
    // get the mesh faces corresponding to the dual graph vertex in "v3"
    // and its parent.
    tCQMFace2* f1 = v->getVertex()->getFace();
    tCQMFace2* f2 = v->getParent()->getVertex()->getFace();
    
    // get the common edge of "f1" and "f2"
    tCQMEdge2* e = f1->getCommonEdge(f2);
    
    // get the half-edge of "e" in "f2"
    tCQMHalfEdge2* he = e->getHalfEdge();
    if (he->getFace() != f2) {
        he = e->getMate(he);
    }
    
    v->getParent()->setQuadSteinerPoint(sp3,he);
}

// -------------------------------------------------------------------
// Method quadrangulatePentR1R4E13E12()
//
// This method  convex quadrangulates  a pentagon with  vertices "p1",
// "p2",  "p3", "p4"  and  "p5" such  that  "p1" and  "p4" are  reflex
// vertices and  p1p2 is  the common edge  of this pentagon  with some
// triangle corresponding to the parent  node of the subtree rooted at
// the node associated with the dual vertex of p1p2p3.
// -------------------------------------------------------------------
void
tCQMCompQuad::quadrangulatePentR1R4E13E12(
                                          tCQMPolygon2* p,
                                          tCQMSpanningTreeVertex* v,
                                          QuadVertexSet& vs,
                                          std::list<tCQMQuadrilateral2*>* lq
                                          )
{
    // get pentagon vertices
    tCQMIndVertex2 p1 = p->getIndVertex(0);
    tCQMIndVertex2 p2 = p->getIndVertex(1);
    tCQMIndVertex2 p3 = p->getIndVertex(2);
    tCQMIndVertex2 p4 = p->getIndVertex(3);
    tCQMIndVertex2 p5 = p->getIndVertex(4);
    
    // we need three Steiner points to convex quadrangulate "p", and the
    // first one  is placed  on the diagonal  p3p1 between "p1"  and the
    // intersection of that diagonal with the supporting line p5p4.
    Line2 l1(p5.getPoint(),p4.getPoint());
    Segment2 s1(p3.getPoint(),p1.getPoint());
    
    // compute intersection point of "l1" and "s1"
    Point2 pt1;
    if (!myAssign(pt1,s1,l1)) {
        pt1 = p3.getPoint();
    }
    
    Point2 q1((p1.getPoint().x() + pt1.x()) / 2,
              (p1.getPoint().y() + pt1.y()) / 2);
    
    // compute the middle of edge p3p4
    Point2 pt2((p1.getPoint().x() + p2.getPoint().x()) / 2,
               (p1.getPoint().y() + p2.getPoint().y()) / 2);
    
    // the second steiner  point is placed in the  region resulting from
    // the intersection of the left half-space defined by the supporting
    // line p4q1,  the right half-space  defined by the  supporting line
    // p5p4 and the triangle pt2p3p1.
    
    // compute intersection point of "l1" and "s2"
    Point2 pt3;
    Segment2 s2(p3.getPoint(),pt2);
    if (IS_EQUAL(pt1.x(),p3.getPoint().x()) &&
        IS_EQUAL(pt1.y(),p3.getPoint().y())) {
        pt3 = pt1;
    }
    else {
        if (!myAssign(pt3,s2,l1)) {
            pt3 = p3.getPoint();
        }
    }
    
    // compute  intersection  point  of  supporting line  p4q1  and  the
    // diagonal p3pt2.
    Line2 l2(p4.getPoint(),q1);
    Point2 pt4;
    if (!myAssign(pt4,s2,l2)) {
        pt4 = p2.getPoint();
    }
    
    double x = (pt1.x() + q1.x() + pt4.x() + pt3.x()) / 4;
    double y = (pt1.y() + q1.y() + pt4.y() + pt3.y()) / 4;
    
    Point2 q2(x,y);
    
    // the third Steiner point is placed off edge p1p2
    Line2 l3(p5.getPoint(),p1.getPoint());
    Segment2 s3(q2,pt2);
    
    // compute intersection point of "l2" and "s2"
    Point2 pt5;
    if (!myAssign(pt5,s3,l3)) {
        pt3 = q2;
    }
    
    x = pt2.x() + 0.3 * (pt5.x() - pt2.x());
    y = pt2.y() + 0.3 * (pt5.y() - pt2.y());
    
    Point2 q3(x,y);
    
    // create Steiner points
    tCQMIndVertex2 sp1 = createSteinerPoint(q1);
    tCQMIndVertex2 sp2 = createSteinerPoint(q2);
    tCQMIndVertex2 sp3 = createSteinerPoint(q3);
    
    // get edge status
    bool ec2 = p->isConstrained(1);
    bool ec3 = p->isConstrained(2);
    bool ec4 = p->isConstrained(3);
    bool ec5 = p->isConstrained(4);
    
    // create four quadrilaterals
    try {
        createQuad( p5, p1,sp1, p4,ec5,false,false,ec4,lq);
        createQuad( p1,sp3,sp2,sp1,false,false,false,false,lq);
        createQuad(sp2,sp3, p2, p3,false,false,ec2,false,lq);
        createQuad(sp1,sp2, p3, p4,false,false,ec3,false,lq);
    }
    catch (const tExceptionObject& xpt) {
        treatException(xpt);
        exit(0);
    }
    
    // get the mesh faces corresponding to the dual graph vertex in "v3"
    // and its parent.
    tCQMFace2* f1 = v->getVertex()->getFace();
    tCQMFace2* f2 = v->getParent()->getVertex()->getFace();
    
    // get the common edge of "f1" and "f2"
    tCQMEdge2* e = f1->getCommonEdge(f2);
    
    // get the half-edge of "e" in "f2"
    tCQMHalfEdge2* he = e->getHalfEdge();
    if (he->getFace() != f2) {
        he = e->getMate(he);
    }
    
    v->getParent()->setQuadSteinerPoint(sp3,he);
}

// -------------------------------------------------------------------
// Method quadrangulatePentR1R4E24E23()
//
// This method  convex quadrangulates  a pentagon with  vertices "p1",
// "p2",  "p3", "p4"  and  "p5" such  that  "p1" and  "p4" are  reflex
// vertices and  p2p3 is  the common edge  of this pentagon  with some
// triangle corresponding to the parent  node of the subtree rooted at
// the node associated with the dual vertex of p2p3p4.
// -------------------------------------------------------------------
void
tCQMCompQuad::quadrangulatePentR1R4E24E23(
                                          tCQMPolygon2* p,
                                          tCQMSpanningTreeVertex* v,
                                          QuadVertexSet& vs,
                                          std::list<tCQMQuadrilateral2*>* lq
                                          )
{
    // get pentagon vertices
    tCQMIndVertex2 p1 = p->getIndVertex(0);
    tCQMIndVertex2 p2 = p->getIndVertex(1);
    tCQMIndVertex2 p3 = p->getIndVertex(2);
    tCQMIndVertex2 p4 = p->getIndVertex(3);
    tCQMIndVertex2 p5 = p->getIndVertex(4);
    
    // we need two  Steiner points to convex quadrangulate  "p", and the
    // first one  is placed  on the diagonal  p2p4 between "p4"  and the
    // intersection of that diagonal with the supporting line p5p1.
    Line2 l1(p5.getPoint(),p1.getPoint());
    Segment2 s1(p2.getPoint(),p4.getPoint());
    
    // compute intersection point of "l1" and "s1"
    Point2 pt1;
    if (!myAssign(pt1,s1,l1)) {
        pt1 = p2.getPoint();
    }
    
    Point2 q1((p4.getPoint().x() + pt1.x()) / 2,
              (p4.getPoint().y() + pt1.y()) / 2);
    
    // the second steiner  point is placed in the  region resulting from
    // the  intersection   of  the  right  half-space   defined  by  the
    // supporting  line  p1q1,  the   left  half-space  defined  by  the
    // supporting line p5p1 and the triangle p2p3p4.
    
    Segment2 s2(p2.getPoint(),p3.getPoint());
    
    // compute intersection point of "l1" and "s2"
    Point2 pt2;
    if (IS_EQUAL(pt1.x(),p2.getPoint().x()) &&
        IS_EQUAL(pt1.y(),p2.getPoint().y())) {
        pt2 = pt1;
    }
    else {
        if (!myAssign(pt2,s2,l1)) {
            pt2 = p2.getPoint();
        }
    }
    
    // compute intersection  point of supporting line p1q1  and the edge
    // p2p3.
    Line2 l2(p1.getPoint(),q1);
    Point2 pt3;
    if (!myAssign(pt3,s2,l2)) {
        pt3 = p3.getPoint();
    }
    
    double x = (pt1.x() + pt2.x() + pt3.x() + q1.x()) / 4;
    double y = (pt1.y() + pt2.y() + pt3.y() + q1.y()) / 4;
    
    Point2 q2(x,y);
    
    // create Steiner points
    tCQMIndVertex2 sp1 = createSteinerPoint(q1);
    tCQMIndVertex2 sp2 = createSteinerPoint(q2);
    
    // get edge status
    bool ec1 = p->isConstrained(0);
    bool ec3 = p->isConstrained(2);
    bool ec4 = p->isConstrained(3);
    bool ec5 = p->isConstrained(4);
    
    // create three quadrilaterals
    try {
        createQuad( p5, p1,sp1, p4,ec5,false,false,ec4,lq);
        createQuad( p1, p2,sp2,sp1,ec1,false,false,false,lq);
        createQuad( p3, p4,sp1,sp2,ec3,false,false,false,lq);
    }
    catch (const tExceptionObject& xpt) {
        treatException(xpt);
        exit(0);
    }
    
    v->setSteinerPoint(sp2,p2.getPoint());
}

// -------------------------------------------------------------------
// Method quadrangulatePentR1R4E13E23()
//
// This method  convex quadrangulates  a pentagon with  vertices "p1",
// "p2",  "p3", "p4"  and  "p5" such  that  "p1" and  "p4" are  reflex
// vertices and  p1p2 is  the common edge  of this pentagon  with some
// triangle corresponding to the parent  node of the subtree rooted at
// the node associated with the dual vertex of p1p2p3.
// -------------------------------------------------------------------
void
tCQMCompQuad::quadrangulatePentR1R4E13E23(
                                          tCQMPolygon2* p,
                                          tCQMSpanningTreeVertex* v,
                                          QuadVertexSet& vs,
                                          std::list<tCQMQuadrilateral2*>* lq
                                          )
{
    // get pentagon vertices
    tCQMIndVertex2 p1 = p->getIndVertex(0);
    tCQMIndVertex2 p2 = p->getIndVertex(1);
    tCQMIndVertex2 p3 = p->getIndVertex(2);
    tCQMIndVertex2 p4 = p->getIndVertex(3);
    tCQMIndVertex2 p5 = p->getIndVertex(4);
    
    // we need two  Steiner points to convex quadrangulate  "p", and the
    // first one  is placed  on the diagonal  p3p1 between "p1"  and the
    // intersection of that diagonal with the supporting line p5p4.
    Line2 l1(p5.getPoint(),p4.getPoint());
    Segment2 s1(p3.getPoint(),p1.getPoint());
    
    // compute intersection point of "l1" and "s1"
    Point2 pt1;
    if (!myAssign(pt1,s1,l1)) {
        pt1 = p3.getPoint();
    }
    
    Point2 q1((p1.getPoint().x() + pt1.x()) / 2,
              (p1.getPoint().y() + pt1.y()) / 2);
    
    // the second steiner  point is placed in the  region resulting from
    // the intersection of the left half-space defined by the supporting
    // line p4q1,  the right half-space  defined by the  supporting line
    // p5p4 and the triangle p1p2p3.
    
    // compute intersection point of "l1" and "s2"
    Segment2 s2(p3.getPoint(),p2.getPoint());
    
    Point2 pt2;
    if (IS_EQUAL(pt1.x(),p3.getPoint().x()) &&
        IS_EQUAL(pt1.y(),p3.getPoint().y())) {
        pt2 = pt1;
    }
    else {
        if (!myAssign(pt2,s2,l1)) {
            pt2 = p3.getPoint();
        }
    }
    
    // compute  intersection  point  of  supporting line  p4q1  and  the
    // diagonal p3pt2.
    Line2 l3(p4.getPoint(),q1);
    Point2 pt3;
    if (!myAssign(pt3,s2,l3)) {
        pt3 = p2.getPoint();
    }
    
    double x = (pt1.x() + q1.x() + pt3.x() + pt2.x()) / 4;
    double y = (pt1.y() + q1.y() + pt3.y() + pt2.y()) / 4;
    
    Point2 q2(x,y);
    
    // create Steiner points
    tCQMIndVertex2 sp1 = createSteinerPoint(q1);
    tCQMIndVertex2 sp2 = createSteinerPoint(q2);
    
    // get edge status
    bool ec1 = p->isConstrained(0);
    bool ec3 = p->isConstrained(2);
    bool ec4 = p->isConstrained(3);
    bool ec5 = p->isConstrained(4);
    
    // create three quadrilaterals
    try {
        createQuad( p5, p1,sp1, p4,ec5,false,false,ec4,lq);
        createQuad( p1, p2,sp2,sp1,ec1,false,false,false,lq);
        createQuad( p3, p4,sp1,sp2,ec3,false,false,false,lq);
    }
    catch (const tExceptionObject& xpt) {
        treatException(xpt);
        exit(0);
    }
    
    v->setSteinerPoint(sp2,p2.getPoint());
}

// -------------------------------------------------------------------
// Method quadrangulatePentR1E13E23()
//
// This method  convex quadrangulates  a pentagon with  vertices "p1",
// "p2", "p3",  "p4" and "p5" such  that only "p1" is  a reflex vertex
// and p2p3  is the  common edge of  this pentagon with  some triangle
// corresponding to the parent node  of the subtree rooted at the node
// associated with the dual vertex of p1p2p3.
// -------------------------------------------------------------------
void
tCQMCompQuad::quadrangulatePentR1E13E23(
                                        tCQMPolygon2* p,
                                        tCQMSpanningTreeVertex* v,
                                        QuadVertexSet& vs,
                                        std::list<tCQMQuadrilateral2*>* lq
                                        )
{
    // get pentagon vertices
    tCQMIndVertex2 p1 = p->getIndVertex(0);
    tCQMIndVertex2 p2 = p->getIndVertex(1);
    tCQMIndVertex2 p3 = p->getIndVertex(2);
    tCQMIndVertex2 p4 = p->getIndVertex(3);
    tCQMIndVertex2 p5 = p->getIndVertex(4);
    
    // we need two  Steiner points to convex quadrangulate  "p", and the
    // first  one  is  placed  inside  the  region  resulting  from  the
    // intersection of  the triangle p1p3p4,  the wedge of "p5"  and the
    // wedge of "p2".
    std::list<Point2> li;
    std::list<Point2> lo;
    li.push_back(p1.getPoint()); li.push_back(p3.getPoint());
    li.push_back(p4.getPoint());
    
    try {
        halfSpaceIntersection(li,p5.getPoint(),p1.getPoint(),lo,
                              CGAL::ON_NEGATIVE_SIDE);
        li.clear();
        halfSpaceIntersection(lo,p2.getPoint(),p1.getPoint(),li,
                              CGAL::ON_POSITIVE_SIDE);
        lo.clear();
    }
    catch (const tExceptionObject& xpt) {
        treatException(xpt);
        exit(0);
    }
    
    // compute the barycenter of the resulting polygon
    Point2 q1;
    double x = 0;
    double y = 0;
    int num = 0;
    while (!li.empty()) {
        q1 = li.front();
        li.pop_front();
        x += q1.x();
        y += q1.y();
        ++num;
    }
    
    q1 = Point2(x / double(num), y / double(num));
    
    // the  second Steiner  point is  placed  at the  barycenter of  the
    // region resulting from the intersection of the triangle q1p2p3 and
    // left half-space defined by the supporting line p4q1.
    li.push_back(q1); li.push_back(p2.getPoint());
    li.push_back(p3.getPoint());
    
    try {
        halfSpaceIntersection(li,p4.getPoint(),q1,lo,CGAL::ON_NEGATIVE_SIDE);
        li.clear();
    }
    catch (const tExceptionObject& xpt) {
        treatException(xpt);
        exit(0);
    }
    
    // compute the barycenter of the resulting polygon
    Point2 q2;
    x = 0;
    y = 0;
    num = 0;
    while (!lo.empty()) {
        q2 = lo.front();
        lo.pop_front();
        x += q2.x();
        y += q2.y();
        ++num;
    }
    
    q2 = Point2(x / double(num), y / double(num));
    
    // create Steiner points
    tCQMIndVertex2 sp1 = createSteinerPoint(q1);
    tCQMIndVertex2 sp2 = createSteinerPoint(q2);
    
    // get edge status
    bool ec1 = p->isConstrained(0);
    bool ec3 = p->isConstrained(2);
    bool ec4 = p->isConstrained(3);
    bool ec5 = p->isConstrained(4);
    
    // create three quadrilaterals
    try {
        createQuad( p5, p1,sp1, p4,ec5,false,false,ec4,lq);
        createQuad( p1, p2,sp2,sp1,ec1,false,false,false,lq);
        createQuad( p3, p4,sp1,sp2,ec3,false,false,false,lq);
    }
    catch (const tExceptionObject& xpt) {
        treatException(xpt);
        exit(0);
    }
    
    v->setSteinerPoint(sp2,p2.getPoint());
}

// -------------------------------------------------------------------
// Method quadrangulatePentR4E24E23()
//
// This method  convex quadrangulates  a pentagon with  vertices "p1",
// "p2", "p3",  "p4" and "p5" such  that only "p4" is  a reflex vertex
// and p2p3  is the  common edge of  this pentagon with  some triangle
// corresponding to the parent node  of the subtree rooted at the node
// associated with the dual vertex of p2p3p4.
// -------------------------------------------------------------------
void
tCQMCompQuad::quadrangulatePentR4E24E23(
                                        tCQMPolygon2* p,
                                        tCQMSpanningTreeVertex* v,
                                        QuadVertexSet& vs,
                                        std::list<tCQMQuadrilateral2*>* lq
                                        )
{
    // get pentagon vertices
    tCQMIndVertex2 p1 = p->getIndVertex(0);
    tCQMIndVertex2 p2 = p->getIndVertex(1);
    tCQMIndVertex2 p3 = p->getIndVertex(2);
    tCQMIndVertex2 p4 = p->getIndVertex(3);
    tCQMIndVertex2 p5 = p->getIndVertex(4);
    
    // we need two  Steiner points to convex quadrangulate  "p", and the
    // first  one  is  placed  inside  the  region  resulting  from  the
    // intersection of  the triangle p1p2p4,  the wedge of "p5"  and the
    // wedge of "p3".
    std::list<Point2> li;
    std::list<Point2> lo;
    li.push_back(p1.getPoint()); li.push_back(p2.getPoint());
    li.push_back(p4.getPoint());
    
    try {
        halfSpaceIntersection(li,p5.getPoint(),p4.getPoint(),lo,
                              CGAL::ON_POSITIVE_SIDE);
        li.clear();
        halfSpaceIntersection(lo,p3.getPoint(),p4.getPoint(),li,
                              CGAL::ON_NEGATIVE_SIDE);
        lo.clear();
    }
    catch (const tExceptionObject& xpt) {
        treatException(xpt);
        exit(0);
    }
    
    // compute the barycenter of the resulting polygon
    Point2 q1;
    double x = 0;
    double y = 0;
    int num = 0;
    while (!li.empty()) {
        q1 = li.front();
        li.pop_front();
        x += q1.x();
        y += q1.y();
        ++num;
    }
    
    q1 = Point2(x / double(num), y / double(num));
    
    // the  second Steiner  point is  placed  at the  barycenter of  the
    // region resulting from the intersection of the triangle q1p2p3 and
    // right half-space defined by the supporting line p1q1.
    li.push_back(q1); li.push_back(p2.getPoint());
    li.push_back(p3.getPoint());
    
    try {
        halfSpaceIntersection(li,p1.getPoint(),q1,lo,CGAL::ON_POSITIVE_SIDE);
        li.clear();
    }
    catch (const tExceptionObject& xpt) {
        treatException(xpt);
        exit(0);
    }
    
    // compute the barycenter of the resulting polygon
    Point2 q2;
    x = 0;
    y = 0;
    num = 0;
    while (!lo.empty()) {
        q2 = lo.front();
        lo.pop_front();
        x += q2.x();
        y += q2.y();
        ++num;
    }
    
    q2 = Point2(x / double(num), y / double(num));
    
    // create Steiner points
    tCQMIndVertex2 sp1 = createSteinerPoint(q1);
    tCQMIndVertex2 sp2 = createSteinerPoint(q2);
    
    // get edge status
    bool ec1 = p->isConstrained(0);
    bool ec3 = p->isConstrained(2);
    bool ec4 = p->isConstrained(3);
    bool ec5 = p->isConstrained(4);
    
    // create three quadrilaterals
    try {
        createQuad( p5, p1,sp1, p4,ec5,false,false,ec4,lq);
        createQuad( p1, p2,sp2,sp1,ec1,false,false,false,lq);
        createQuad( p3, p4,sp1,sp2,ec3,false,false,false,lq);
    }
    catch (const tExceptionObject& xpt) {
        treatException(xpt);
        exit(0);
    }
    
    v->setSteinerPoint(sp2,p2.getPoint());
}

// -------------------------------------------------------------------
// Method quadrangulatePentR1E13E12()
//
// This method  convex quadrangulates  a pentagon with  vertices "p1",
// "p2", "p3",  "p4" and "p5" such  that only "p1" is  a reflex vertex
// and p1p2  is the  common edge of  this pentagon with  some triangle
// corresponding to the parent node  of the subtree rooted at the node
// associated with the dual vertex of p1p2p3.
// -------------------------------------------------------------------
void
tCQMCompQuad::quadrangulatePentR1E13E12(
                                        tCQMPolygon2* p,
                                        tCQMSpanningTreeVertex* v,
                                        QuadVertexSet& vs,
                                        std::list<tCQMQuadrilateral2*>* lq
                                        )
{
    // get pentagon vertices
    tCQMIndVertex2 p1 = p->getIndVertex(0);
    tCQMIndVertex2 p2 = p->getIndVertex(1);
    tCQMIndVertex2 p3 = p->getIndVertex(2);
    tCQMIndVertex2 p4 = p->getIndVertex(3);
    tCQMIndVertex2 p5 = p->getIndVertex(4);
    
    // we need three Steiner points to convex quadrangulate "p", and the
    // first  one  is  placed  inside  the  region  resulting  from  the
    // intersection of  the triangle p1p3p4,  the wedge of "p5"  and the
    // wedge of "p2".
    std::list<Point2> li;
    std::list<Point2> lo;
    li.push_back(p1.getPoint()); li.push_back(p3.getPoint());
    li.push_back(p4.getPoint());
    
    try {
        halfSpaceIntersection(li,p5.getPoint(),p1.getPoint(),lo,
                              CGAL::ON_NEGATIVE_SIDE);
        li.clear();
        halfSpaceIntersection(lo,p2.getPoint(),p1.getPoint(),li,
                              CGAL::ON_POSITIVE_SIDE);
        lo.clear();
    }
    catch (const tExceptionObject& xpt) {
        treatException(xpt);
        exit(0);
    }
    
    // compute the barycenter of the resulting polygon
    Point2 q1;
    double x = 0;
    double y = 0;
    int num = 0;
    while (!li.empty()) {
        q1 = li.front();
        li.pop_front();
        x += q1.x();
        y += q1.y();
        ++num;
    }
    
    q1 = Point2(x / double(num), y / double(num));
    
    // the  second Steiner  point is  placed  at the  barycenter of  the
    // region resulting from the intersection of the triangle q1p2p3 and
    // left half-space defined by the supporting line p4q1.
    li.push_back(q1); li.push_back(p2.getPoint());
    li.push_back(p3.getPoint());
    
    try {
        halfSpaceIntersection(li,p4.getPoint(),q1,lo,CGAL::ON_NEGATIVE_SIDE);
        li.clear();
    }
    catch (const tExceptionObject& xpt) {
        treatException(xpt);
        exit(0);
    }
    
    // compute the barycenter of the resulting polygon
    Point2 q2;
    x = 0;
    y = 0;
    num = 0;
    while (!lo.empty()) {
        q2 = lo.front();
        lo.pop_front();
        x += q2.x();
        y += q2.y();
        ++num;
    }
    
    q2 = Point2(x / double(num), y / double(num));
    
    // the third  Steiner point  is off edge  p1p2 between "p2"  and the
    // intersection point of the supporting line p3q2 and edge p1p2.
    Line2 l(p3.getPoint(),q2);
    Segment2 s(p1.getPoint(),p2.getPoint());
    
    // compute intersection point of "l" and "s"
    Point2 pt1;
    if (!myAssign(pt1,s,l)) {
        // should never happen
        pt1 = p1.getPoint();
    }
    
    x = (p2.getPoint().x() + pt1.x()) / 2;
    y = (p2.getPoint().y() + pt1.y()) / 2;
    
    Point2 q3(x + 0.3 * (q2.x() - x), y + 0.3 * (q2.y() - y));
    
    // create Steiner points
    tCQMIndVertex2 sp1 = createSteinerPoint(q1);
    tCQMIndVertex2 sp2 = createSteinerPoint(q2);
    tCQMIndVertex2 sp3 = createSteinerPoint(q3);
    
    // get edge status
    bool ec2 = p->isConstrained(1);
    bool ec3 = p->isConstrained(2);
    bool ec4 = p->isConstrained(3);
    bool ec5 = p->isConstrained(4);
    
    // create four quadrilaterals
    try {
        createQuad( p5, p1,sp1, p4,ec5,false,false,ec4,lq);
        createQuad( p1,sp3,sp2,sp1,false,false,false,false,lq);
        createQuad( p3, p4,sp1,sp2,ec3,false,false,false,lq);
        createQuad( p3,sp2,sp3, p2,false,false,false,ec2,lq);
    }
    catch (const tExceptionObject& xpt) {
        treatException(xpt);
        exit(0);
    }
    
    // get the mesh faces corresponding to the dual graph vertex in "v3"
    // and its parent.
    tCQMFace2* f1 = v->getVertex()->getFace();
    tCQMFace2* f2 = v->getParent()->getVertex()->getFace();
    
    // get the common edge of "f1" and "f2"
    tCQMEdge2* e = f1->getCommonEdge(f2);
    
    // get the half-edge of "e" in "f2"
    tCQMHalfEdge2* he = e->getHalfEdge();
    if (he->getFace() != f2) {
        he = e->getMate(he);
    }
    
    v->getParent()->setQuadSteinerPoint(sp3,he);
}

// -------------------------------------------------------------------
// Method quadrangulatePentR4E24E34()
//
// This method  convex quadrangulates  a pentagon with  vertices "p1",
// "p2", "p3",  "p4" and "p5" such  that only "p4" is  a reflex vertex
// and p3p4  is the  common edge of  this pentagon with  some triangle
// corresponding to the parent node  of the subtree rooted at the node
// associated with the dual vertex of p2p3p4.
// -------------------------------------------------------------------
void
tCQMCompQuad::quadrangulatePentR4E24E34(
                                        tCQMPolygon2* p,
                                        tCQMSpanningTreeVertex* v,
                                        QuadVertexSet& vs,
                                        std::list<tCQMQuadrilateral2*>* lq
                                        )
{
    // get pentagon vertices
    tCQMIndVertex2 p1 = p->getIndVertex(0);
    tCQMIndVertex2 p2 = p->getIndVertex(1);
    tCQMIndVertex2 p3 = p->getIndVertex(2);
    tCQMIndVertex2 p4 = p->getIndVertex(3);
    tCQMIndVertex2 p5 = p->getIndVertex(4);
    
    // we need three Steiner points to convex quadrangulate "p", and the
    // first  one  is  placed  inside  the  region  resulting  from  the
    // intersection of  the triangle p1p2p4,  the wedge of "p5"  and the
    // wedge of "p3".
    std::list<Point2> li;
    std::list<Point2> lo;
    li.push_back(p1.getPoint()); li.push_back(p2.getPoint());
    li.push_back(p4.getPoint());
    
    try {
        halfSpaceIntersection(li,p5.getPoint(),p4.getPoint(),lo,
                              CGAL::ON_POSITIVE_SIDE);
        li.clear();
        halfSpaceIntersection(lo,p3.getPoint(),p4.getPoint(),li,
                              CGAL::ON_NEGATIVE_SIDE);
        lo.clear();
    }
    catch (const tExceptionObject& xpt) {
        treatException(xpt);
        exit(0);
    }
    
    // compute the barycenter of the resulting polygon
    Point2 q1;
    double x = 0;
    double y = 0;
    int num = 0;
    while (!li.empty()) {
        q1 = li.front();
        li.pop_front();
        x += q1.x();
        y += q1.y();
        ++num;
    }
    
    q1 = Point2(x / double(num), y / double(num));
    
    // the  second Steiner  point is  placed  at the  barycenter of  the
    // region resulting from the intersection of the triangle q1p2p3 and
    // the right half-space defined by the supporting line p1q1.
    li.push_back(q1); li.push_back(p2.getPoint());
    li.push_back(p3.getPoint());
    
    try {
        halfSpaceIntersection(li,p1.getPoint(),q1,lo,CGAL::ON_POSITIVE_SIDE);
        li.clear();
    }
    catch (const tExceptionObject& xpt) {
        treatException(xpt);
        exit(0);
    }
    
    // compute the barycenter of the resulting polygon
    Point2 q2;
    x = 0;
    y = 0;
    num = 0;
    while (!lo.empty()) {
        q2 = lo.front();
        lo.pop_front();
        x += q2.x();
        y += q2.y();
        ++num;
    }
    
    q2 = Point2(x / double(num), y / double(num));
    
    // the third  Steiner point  is off edge  p3p4 between "p2"  and the
    // intersection point of the supporting line p2q2 and edge p3p4.
    Line2 l(p2.getPoint(),q2);
    Segment2 s(p3.getPoint(),p4.getPoint());
    
    // compute intersection point of "l" and "s"
    Point2 pt1;
    if (!myAssign(pt1,s,l)) {
        // should never happen
        pt1 = p4.getPoint();
    }
    
    x = (p3.getPoint().x() + pt1.x()) / 2;
    y = (p3.getPoint().y() + pt1.y()) / 2;
    
    Point2 q3(x + 0.3 * (q2.x() - x), y + 0.3 * (q2.y() - y));
    
    // create Steiner points
    tCQMIndVertex2 sp1 = createSteinerPoint(q1);
    tCQMIndVertex2 sp2 = createSteinerPoint(q2);
    tCQMIndVertex2 sp3 = createSteinerPoint(q3);
    
    // get edge status
    bool ec1 = p->isConstrained(0);
    bool ec2 = p->isConstrained(1);
    bool ec4 = p->isConstrained(3);
    bool ec5 = p->isConstrained(4);
    
    // create four quadrilaterals
    try {
        createQuad( p5, p1,sp1, p4,ec5,false,false,ec4,lq);
        createQuad( p3,sp3,sp2, p2,false,false,false,ec2,lq);
        createQuad( p4,sp1,sp2,sp3,false,false,false,false,lq);
        createQuad( p1, p2,sp2,sp1,ec1,false,false,false,lq);
    }
    catch (const tExceptionObject& xpt) {
        treatException(xpt);
        exit(0);
    }
    
    // get the mesh faces corresponding to the dual graph vertex in "v3"
    // and its parent.
    tCQMFace2* f1 = v->getVertex()->getFace();
    tCQMFace2* f2 = v->getParent()->getVertex()->getFace();
    
    // get the common edge of "f1" and "f2"
    tCQMEdge2* e = f1->getCommonEdge(f2);
    
    // get the half-edge of "e" in "f2"
    tCQMHalfEdge2* he = e->getHalfEdge();
    if (he->getFace() != f2) {
        he = e->getMate(he);
    }
    
    v->getParent()->setQuadSteinerPoint(sp3,he);
}

// -------------------------------------------------------------------
// Method quadrangulateHexCCCCCC()
// -------------------------------------------------------------------
void
tCQMCompQuad::quadrangulateHexCCCCCC(
                                     tCQMPolygon2* p,
                                     std::list<tCQMQuadrilateral2*>* lq
                                     )
{
    // get all six vertices of the hexagon
    tCQMIndVertex2 p1 = p->getIndVertex(0);
    tCQMIndVertex2 p2 = p->getIndVertex(1);
    tCQMIndVertex2 p3 = p->getIndVertex(2);
    tCQMIndVertex2 p4 = p->getIndVertex(3);
    tCQMIndVertex2 p5 = p->getIndVertex(4);
    tCQMIndVertex2 p6 = p->getIndVertex(5);
    
    // get edge status
    bool ec1 = p->isConstrained(0);
    bool ec2 = p->isConstrained(1);
    bool ec3 = p->isConstrained(2);
    bool ec4 = p->isConstrained(3);
    bool ec5 = p->isConstrained(4);
    bool ec6 = p->isConstrained(5);
    
    // create quadrilaterals
    try {
        createQuad(p1,p2,p3,p4,ec1,ec2,ec3,false,lq);
        createQuad(p1,p4,p5,p6,false,ec4,ec5,ec6,lq);
    }
    catch (const tExceptionObject& xpt) {
        treatException(xpt);
        exit(0);
    }
}

// -------------------------------------------------------------------
// Method quadrangulateHexRCCCCC1()
// -------------------------------------------------------------------
void
tCQMCompQuad::quadrangulateHexRCCCCC1(
                                      const tCQMIndVertex2& p1,
                                      const tCQMIndVertex2& p2,
                                      const tCQMIndVertex2& p3,
                                      const tCQMIndVertex2& p4,
                                      const tCQMIndVertex2& p5,
                                      const tCQMIndVertex2& p6,
                                      bool ec1,
                                      bool ec2,
                                      bool ec3,
                                      bool ec4,
                                      bool ec5,
                                      bool ec6,
                                      std::list<tCQMQuadrilateral2*>* lq
                                      )
{
    // we assume that  "p1" is on the left of  the supporting line p3p5,
    // so that we can convex quadrangulate "p" by using only one Steiner
    // point.
    
    // Steiner point can be placed in the triangle p1pt1pt2, where "pt1"
    // and "pt2" are the intersection points of the segment p3p5 and the
    // supporting lines p2p1 and p6p1.
    Line2 l1(p2.getPoint(),p1.getPoint());
    Line2 l2(p6.getPoint(),p1.getPoint());
    Segment2 s(p3.getPoint(),p5.getPoint());
    Point2 pt1;
    Point2 pt2;
    if (!myAssign(pt1,s,l1)) pt1 = p5.getPoint();
    if (!myAssign(pt2,s,l2)) pt2 = p3.getPoint();
    
    // compute barycenter of the triangle p1pt1pt2
    Point2 q((pt1.x() + pt2.x() + p1.getPoint().x()) / 3,
             (pt1.y() + pt2.y() + p1.getPoint().y()) / 3);
    
    // create Steiner point
    tCQMIndVertex2 sp = createSteinerPoint(q);
    
    // create quadrilaterals
    try {
        createQuad(p1,p2,p3,sp,ec1,ec2,false,false,lq);
        createQuad(p1,sp,p5,p6,false,false,ec5,ec6,lq);
        createQuad(sp,p3,p4,p5,false,ec3,ec4,false,lq);
    }
    catch (const tExceptionObject& xpt) {
        treatException(xpt);
        exit(0);
    }
}

// -------------------------------------------------------------------
// Method quadrangulateHexRCCCCC2()
// -------------------------------------------------------------------
void
tCQMCompQuad::quadrangulateHexRCCCCC2(
                                      const tCQMIndVertex2& p1,
                                      const tCQMIndVertex2& p2,
                                      const tCQMIndVertex2& p3,
                                      const tCQMIndVertex2& p4,
                                      const tCQMIndVertex2& p5,
                                      const tCQMIndVertex2& p6,
                                      bool ec1,
                                      bool ec2,
                                      bool ec3,
                                      bool ec4,
                                      bool ec5,
                                      bool ec6,
                                      std::list<tCQMQuadrilateral2*>* lq
                                      )
throw (tExceptionObject)
{
    // we assume that  "p3" and "p5" does not see each other.
    
    // the first Steiner point is  placed in the triangle defined by the
    // intersection of the  wedge of "p1" and "p".   In order to compute
    // this  triangle, we  compute the  two intersection  points  of the
    // supporting lines  p2p1 and p6p1  with the edges p2p3,  p3p4, p4p5
    // and p5p6.
    
    // compute intersection of the hexagon and the wedge of "p1"
    Line2 l1(p2.getPoint(),p1.getPoint());
    Line2 l2(p6.getPoint(),p1.getPoint());
    if (l1.oriented_side(p4.getPoint()) != CGAL::ON_NEGATIVE_SIDE) {
        Point2 pt1;
        Point2 pt2;
        
        // line p2p1 intersects edge p3p4
        Segment2 s1(p3.getPoint(),p4.getPoint());
        if (!myAssign(pt1,s1,l1)) {
            std::stringstream ss (std::stringstream::in | std::stringstream::out);
            ss << "RCCCCC2(): cannot recover from numerical instability";
            throw tExceptionObject(__FILE__,__LINE__,ss.str());
        }
        
        // line p6p1 intersects either edge p3p4 or edge p2p3
        if (!assign(pt2,intersection(s1,l2))) {
            Segment2 s2(p2.getPoint(),p3.getPoint());
            if (!myAssign(pt2,s2,l2)) {
                std::stringstream ss (std::stringstream::in | std::stringstream::out);
                ss << "RCCCCC2(): cannot recover from numerical instability";
                throw tExceptionObject(__FILE__,__LINE__,ss.str());
            }
        }
        
        // compute first Steiner point
        Point2 q((pt1.x() + pt2.x() + p1.getPoint().x()) / 3,
                 (pt1.y() + pt2.y() + p1.getPoint().y()) / 3);
        
        // create the first Steiner point
        tCQMIndVertex2 sp = createSteinerPoint(q);
        
        // create one quadrilateral
        try {
            createQuad(p1,p2,p3,sp,ec1,ec2,false,false,lq);
        }
        catch (const tExceptionObject& xpt) {
            treatException(xpt);
            exit(0);
        }
        
        // we are left  with the hexagon qp3p4p5p6p1, which  can be convex
        // quadrangulate by inserting a  Steiner point in its intersection
        // with the wedge of q1.
        quadrangulateHexRCCCCC1(sp,p3,p4,p5,p6,p1,false,ec3,ec4,ec5,ec6,
                                false,lq);
    }
    else {
        Point2 pt1;
        Point2 pt2;
        
        // line p2p1 intersects either edge p4p5 or edge p5p6
        Segment2 s1(p4.getPoint(),p5.getPoint());
        if (!myAssign(pt1,s1,l1)) {
            Segment2 s2(p5.getPoint(),p6.getPoint());
            if (!myAssign(pt1,s2,l1)) {
                std::stringstream ss (std::stringstream::in | std::stringstream::out);
                ss << "RCCCCC2(): cannot recover from numerical instability";
                throw tExceptionObject(__FILE__,__LINE__,ss.str());
            }
        }
        
        // line p6p1 intersects edge p4p5
        if (!myAssign(pt2,s1,l2)) {
            std::stringstream ss (std::stringstream::in | std::stringstream::out);
            ss << "RCCCCC2(): cannot recover from numerical instability";
            throw tExceptionObject(__FILE__,__LINE__,ss.str());
        }
        
        // compute first Steiner point
        Point2 q((pt1.x() + pt2.x() + p1.getPoint().x()) / 3,
                 (pt1.y() + pt2.y() + p1.getPoint().y()) / 3);
        
        // create the first Steiner point
        tCQMIndVertex2 sp = createSteinerPoint(q);
        
        // create one quadrilateral
        try {
            createQuad(p1,sp,p5,p6,false,false,ec5,ec6,lq);
        }
        catch (const tExceptionObject& xpt) {
            treatException(xpt);
            exit(0);
        }
        
        // we are left  with the hexagon qp1p2p3p4p5, which  can be convex
        // quadrangulate by inserting a  Steiner point in its intersection
        // with the wedge of q1.
        quadrangulateHexRCCCCC1(sp,p1,p2,p3,p4,p5,false,ec1,ec2,ec3,ec4,
                                false,lq);
    }
}

// -------------------------------------------------------------------
// Method quadrangulateHexRCCCCC()
// -------------------------------------------------------------------
void
tCQMCompQuad::quadrangulateHexRCCCCC(
                                     tCQMPolygon2* p,
                                     std::list<tCQMQuadrilateral2*>* lq
                                     )
{
    // get all six vertices of the hexagon
    tCQMIndVertex2 p1 = p->getIndVertex(0);
    tCQMIndVertex2 p2 = p->getIndVertex(1);
    tCQMIndVertex2 p3 = p->getIndVertex(2);
    tCQMIndVertex2 p4 = p->getIndVertex(3);
    tCQMIndVertex2 p5 = p->getIndVertex(4);
    tCQMIndVertex2 p6 = p->getIndVertex(5);
    
    // get edge status
    bool ec1 = p->isConstrained(0);
    bool ec2 = p->isConstrained(1);
    bool ec3 = p->isConstrained(2);
    bool ec4 = p->isConstrained(3);
    bool ec5 = p->isConstrained(4);
    bool ec6 = p->isConstrained(5);
    
    // if "p4"  is on the right of  the supporting line p1p2  and on the
    // left   of  the  supporting   line  p6p1,   then  we   can  convex
    // quadrangulate "p" without using any Steiner point.
    Line2 l1(p2.getPoint(),p1.getPoint());
    Line2 l2(p6.getPoint(),p1.getPoint());
    if ((l1.oriented_side(p4.getPoint()) == CGAL::ON_NEGATIVE_SIDE) &&
        (l2.oriented_side(p4.getPoint()) == CGAL::ON_POSITIVE_SIDE)) {
        // create quadrilaterals
        try {
            createQuad(p1,p2,p3,p4,ec1,ec2,ec3,false,lq);
            createQuad(p1,p4,p5,p6,false,ec4,ec5,ec6,lq);
        }
        catch (const tExceptionObject& xpt) {
            treatException(xpt);
            exit(0);
        }
    }
    else {
        // if "p1" is on the left of the supporting line p3p5, then we can
        // convex quadrangulate "p" by using only one Steiner point.
        Line2 l3(p3.getPoint(),p5.getPoint());
        if (l3.oriented_side(p1.getPoint()) == CGAL::ON_POSITIVE_SIDE) {
            quadrangulateHexRCCCCC1(p1,p2,p3,p4,p5,p6,ec1,ec2,ec3,ec4,ec5,
                                    ec6,lq);
        }
        else {
            // "p1" is on  the right of the supporting line  p3p5 so that we
            // can convex quadrangulate "p" using two Steiner points.
            
            try {
                quadrangulateHexRCCCCC2(p1,p2,p3,p4,p5,p6,ec1,ec2,ec3,ec4,ec5,
                                        ec6,lq);
            }
            catch (const tExceptionObject& xpt) {
                treatException(xpt);
                exit(0);
            }
            
        }
    }
}

// -------------------------------------------------------------------
// Method quadrangulateHexRCRCCC1()
// -------------------------------------------------------------------
void
tCQMCompQuad::quadrangulateHexRCRCCC1(
                                      const tCQMIndVertex2& p1,
                                      const tCQMIndVertex2& p2,
                                      const tCQMIndVertex2& p3,
                                      const tCQMIndVertex2& p4,
                                      const tCQMIndVertex2& p5,
                                      const tCQMIndVertex2& p6,
                                      bool ec1,
                                      bool ec2,
                                      bool ec3,
                                      bool ec4,
                                      bool ec5,
                                      bool ec6,
                                      std::list<tCQMQuadrilateral2*>* lq
                                      )
{
    // a  pre-condition   to  use  this  method  is   that  the  polygon
    // p1p2p3p4p5p6  has two reflex  vertices, "p1"  and "p3",  and both
    // "p1" and "p3" see "p5".
    
    // we  can place  a Steiner  point in  the triangle  defined  by the
    // intersection of the  triangle p1p3p5, wedge of "p1"  and wedge of
    // "p3".
    
    // compute  intersection of  the triangle  p1p3p5 and  the  wedge of
    // "p1".
    std::list<Point2> li;
    std::list<Point2> lo;
    
    li.push_back(p1.getPoint()); li.push_back(p3.getPoint());
    li.push_back(p5.getPoint());
    
    // compute intersection of "p" and the wedge of "p1"
    try {
        halfSpaceIntersection(li,p2.getPoint(),p1.getPoint(),lo,
                              CGAL::ON_POSITIVE_SIDE);
        li.clear();
        halfSpaceIntersection(lo,p6.getPoint(),p1.getPoint(),li,
                              CGAL::ON_NEGATIVE_SIDE);
        lo.clear();
    }
    catch (const tExceptionObject& xpt) {
        treatException(xpt);
        exit(0);
    }
    
    // compute intersection of the resulting polygon above and the wedge
    // of "p3".
    try {
        halfSpaceIntersection(li,p2.getPoint(),p3.getPoint(),lo,
                              CGAL::ON_NEGATIVE_SIDE);
        li.clear();
        halfSpaceIntersection(lo,p4.getPoint(),p3.getPoint(),li,
                              CGAL::ON_POSITIVE_SIDE);
        lo.clear();
    }
    catch (const tExceptionObject& xpt) {
        treatException(xpt);
        exit(0);
    }
    
    // compute the barycenter of the resulting polygon
    Point2 q;
    double x = 0;
    double y = 0;
    int num = 0;
    while (!li.empty()) {
        q = li.front();
        li.pop_front();
        x += q.x();
        y += q.y();
        ++num;
    }
    
    q = Point2(x / double(num), y / double(num));
    
    // create Steiner point
    tCQMIndVertex2 sp = createSteinerPoint(q);
    
    // create quadrilaterals
    try {
        createQuad(p1,p2,p3,sp,ec1,ec2,false,false,lq);
        createQuad(sp,p3,p4,p5,false,ec3,ec4,false,lq);
        createQuad(p1,sp,p5,p6,false,false,ec5,ec6,lq);
    }
    catch (const tExceptionObject& xpt) {
        treatException(xpt);
        exit(0);
    }
}

// -------------------------------------------------------------------
// Method quadrangulateHexRCRCCC3()
// -------------------------------------------------------------------
void
tCQMCompQuad::quadrangulateHexRCRCCC3(
                                      const tCQMIndVertex2& p1,
                                      const tCQMIndVertex2& p2,
                                      const tCQMIndVertex2& p3,
                                      const tCQMIndVertex2& p4,
                                      const tCQMIndVertex2& p5,
                                      const tCQMIndVertex2& p6,
                                      bool ec1,
                                      bool ec2,
                                      bool ec3,
                                      bool ec4,
                                      bool ec5,
                                      bool ec6,
                                      std::list<tCQMQuadrilateral2*>* lq
                                      )
{
    // a  pre-condition   to  use  this  method  is   that  the  polygon
    // p1p2p3p4p5p6 has two reflex vertices, "p1" and "p3", and only one
    // of them can see "p5".
    
    // we  can place  a Steiner  point in  the triangle  defined  by the
    // intersection of the wedge of "p1", the wedge of "p3" and the left
    // half-space defined by the supporting line p1p3.
    
    // compute intersection of the polygon p1p2p3p4p5p6 and the wedge of
    // "p1".
    std::list<Point2> li;
    std::list<Point2> lo;
    
    li.push_back(p1.getPoint()); li.push_back(p2.getPoint());
    li.push_back(p3.getPoint()); li.push_back(p4.getPoint());
    li.push_back(p5.getPoint()); li.push_back(p6.getPoint());
    
    try {
        // compute intersection of "p" and the wedge of "p1".
        halfSpaceIntersection(li,p2.getPoint(),p1.getPoint(),lo,
                              CGAL::ON_POSITIVE_SIDE);
        li.clear();
        halfSpaceIntersection(lo,p6.getPoint(),p1.getPoint(),li,
                              CGAL::ON_NEGATIVE_SIDE);
        lo.clear();
        
        // compute  intersection of  the resulting  polygon above  and the
        // wedge of "p3".
        halfSpaceIntersection(li,p2.getPoint(),p3.getPoint(),lo,
                              CGAL::ON_NEGATIVE_SIDE);
        li.clear();
        halfSpaceIntersection(lo,p4.getPoint(),p3.getPoint(),li,
                              CGAL::ON_POSITIVE_SIDE);
        lo.clear();
        
        // compute  intersection of  the resulting  polygon above  and the
        // left half-space of the supporting line p1p3.
        halfSpaceIntersection(li,p1.getPoint(),p3.getPoint(),lo,
                              CGAL::ON_NEGATIVE_SIDE);
        li.clear();
    }
    catch (const tExceptionObject& xpt) {
        treatException(xpt);
        exit(0);
    }
    
    // compute the barycenter of the resulting polygon
    Point2 q;
    double x = 0;
    double y = 0;
    int num = 0;
    while (!lo.empty()) {
        q = lo.front();
        lo.pop_front();
        x += q.x();
        y += q.y();
        ++num;
    }
    
    q = Point2(x / double(num), y / double(num));
    
    // create Steiner point
    tCQMIndVertex2 sp = createSteinerPoint(q);
    
    // create one quadrilateral
    try {
        createQuad(p1,p2,p3,sp,ec1,ec2,false,false,lq);
    }
    catch (const tExceptionObject& xpt) {
        treatException(xpt);
        exit(0);
    }
    
    // create  a new hexagon,  spp3p4p5p6p1, which  has only  one reflex
    // vertex and therefore  can be convex quadrangulate with  at most 2
    // Steiner points.
    
    // create hexagon
    tCQMPolygon2* paux = (tCQMPolygon2 *) new tCQMPolygon2(6);
    
    paux->insert(0,sp);
    paux->insert(1,p3);
    paux->insert(2,p4);
    paux->insert(3,p5);
    paux->insert(4,p6);
    paux->insert(5,p1);
    
    paux->setConstraint(0,false);
    paux->setConstraint(1,ec3);
    paux->setConstraint(2,ec4);
    paux->setConstraint(3,ec5);
    paux->setConstraint(4,ec6);
    paux->setConstraint(5,false);
    
    // convex quadrangulate it
    quadrangulateHexRCCCCC(paux,lq);
    
    // delete hexagon
    delete paux;
}

// -------------------------------------------------------------------
// Method quadrangulateHexRCRCCC()
// -------------------------------------------------------------------
void
tCQMCompQuad::quadrangulateHexRCRCCC(
                                     tCQMPolygon2* p,
                                     std::list<tCQMQuadrilateral2*>* lq
                                     )
{
    // get all six vertices of the hexagon
    tCQMIndVertex2 p1 = p->getIndVertex(0);
    tCQMIndVertex2 p2 = p->getIndVertex(1);
    tCQMIndVertex2 p3 = p->getIndVertex(2);
    tCQMIndVertex2 p4 = p->getIndVertex(3);
    tCQMIndVertex2 p5 = p->getIndVertex(4);
    tCQMIndVertex2 p6 = p->getIndVertex(5);
    
    // get edge status
    bool ec1 = p->isConstrained(0);
    bool ec2 = p->isConstrained(1);
    bool ec3 = p->isConstrained(2);
    bool ec4 = p->isConstrained(3);
    bool ec5 = p->isConstrained(4);
    bool ec6 = p->isConstrained(5);
    
    // if both "p1"  and "p3" can see "p5"  then we convex quadrangulate
    // "p" by using  only one Steiner point. "p1" and  "p3" can see "p5"
    // if "p1" is on the left of the supporting line p3p5 and "p3" is on
    // the right of the supporting line p1p5.
    Line2 l1(p3.getPoint(),p5.getPoint());
    Line2 l2(p1.getPoint(),p5.getPoint());
    if ((l1.oriented_side(p1.getPoint()) == CGAL::ON_POSITIVE_SIDE) &&
        (l2.oriented_side(p3.getPoint()) == CGAL::ON_NEGATIVE_SIDE)) {
        quadrangulateHexRCRCCC1(p1,p2,p3,p4,p5,p6,ec1,ec2,ec3,ec4,ec5,
                                ec6,lq);
    }
    else {
        // "p3" and "p5" cannot see each other, so we convex quadrangulate
        // "p" by using at most 3 Steiner points.
        quadrangulateHexRCRCCC3(p1,p2,p3,p4,p5,p6,ec1,ec2,ec3,ec4,ec5,
                                ec6,lq);
    }
}

// -------------------------------------------------------------------
// Method quadrangulateHexRRCCCC()
// -------------------------------------------------------------------
void
tCQMCompQuad::quadrangulateHexRRCCCC(
                                     tCQMPolygon2* p,
                                     std::list<tCQMQuadrilateral2*>* lq
                                     )
{
    // get all six vertices of the hexagon
    tCQMIndVertex2 p1 = p->getIndVertex(0);
    tCQMIndVertex2 p2 = p->getIndVertex(1);
    tCQMIndVertex2 p3 = p->getIndVertex(2);
    tCQMIndVertex2 p4 = p->getIndVertex(3);
    tCQMIndVertex2 p5 = p->getIndVertex(4);
    tCQMIndVertex2 p6 = p->getIndVertex(5);
    
    // get edge status
    bool ec1 = p->isConstrained(0);
    bool ec2 = p->isConstrained(1);
    bool ec3 = p->isConstrained(2);
    bool ec4 = p->isConstrained(3);
    bool ec5 = p->isConstrained(4);
    bool ec6 = p->isConstrained(5);
    
    // "p1"   and  "p2"  are   reflex  vertices,   and  we   can  convex
    // quadrangulate  "p"  by placing  a  Steiner  point,  "sp", in  the
    // intersection region  of the wedge  of "p1", the  right half-space
    // defined  by the  supporting  line p1p5  and  the left  half-space
    // defined  by the  supporting line  p2p4. The  hexagon spp1p2p3p4p5
    // will have  2 reflex  vertices and it  will be of  type RCRCCCC-1,
    // which  means that  it can  be convex  quadrangulate by  using one
    // Steiner point.
    
    // compute intersection of the polygon p1p2p3p4p5p6 and the wedge of
    // "p1".
    std::list<Point2> li;
    std::list<Point2> lo;
    
    li.push_back(p1.getPoint()); li.push_back(p2.getPoint());
    li.push_back(p3.getPoint()); li.push_back(p4.getPoint());
    li.push_back(p5.getPoint()); li.push_back(p6.getPoint());
    
    try {
        // compute intersection of "p" and the wedge of "p1".
        halfSpaceIntersection(li,p2.getPoint(),p1.getPoint(),lo,
                              CGAL::ON_POSITIVE_SIDE);
        li.clear();
        halfSpaceIntersection(lo,p6.getPoint(),p1.getPoint(),li,
                              CGAL::ON_NEGATIVE_SIDE);
        lo.clear();
        
        // compute  intersection of  the resulting  polygon above  and the
        // right half-space defined by the supporting line p1p5.
        halfSpaceIntersection(li,p1.getPoint(),p5.getPoint(),lo,
                              CGAL::ON_POSITIVE_SIDE);
        li.clear();
        
        // compute  intersection of  the resulting  polygon above  and the
        // left half-space defined by the supporting line p2p4.
        halfSpaceIntersection(lo,p2.getPoint(),p4.getPoint(),li,
                              CGAL::ON_NEGATIVE_SIDE);
        lo.clear();
    }
    catch (const tExceptionObject& xpt) {
        treatException(xpt);
        exit(0);
    }
    
    // compute the barycenter of the resulting polygon
    Point2 q;
    double x = 0;
    double y = 0;
    int num = 0;
    while (!li.empty()) {
        q = li.front();
        li.pop_front();
        x += q.x();
        y += q.y();
        ++num;
    }
    
    q = Point2(x / double(num), y / double(num));
    
    // create Steiner point
    tCQMIndVertex2 sp = createSteinerPoint(q);
    
    // create one quadrilateral
    try {
        createQuad(p1,sp,p5,p6,false,false,ec5,ec6,lq);
    }
    catch (const tExceptionObject& xpt) {
        treatException(xpt);
        exit(0);
    }
    
    // convex quadrangulate the hexagon spp1p2p3p4p5
    quadrangulateHexRCRCCC1(sp,p1,p2,p3,p4,p5,false,ec1,ec2,ec3,ec4,
                            false,lq);
}

// -------------------------------------------------------------------
// Method quadrangulateHexRCCRCC()
// -------------------------------------------------------------------
void
tCQMCompQuad::quadrangulateHexRCCRCC(
                                     tCQMPolygon2* p,
                                     std::list<tCQMQuadrilateral2*>* lq
                                     )
{
    // get all six vertices of the hexagon
    tCQMIndVertex2 p1 = p->getIndVertex(0);
    tCQMIndVertex2 p2 = p->getIndVertex(1);
    tCQMIndVertex2 p3 = p->getIndVertex(2);
    tCQMIndVertex2 p4 = p->getIndVertex(3);
    tCQMIndVertex2 p5 = p->getIndVertex(4);
    tCQMIndVertex2 p6 = p->getIndVertex(5);
    
    // get edge status
    bool ec1 = p->isConstrained(0);
    bool ec2 = p->isConstrained(1);
    bool ec3 = p->isConstrained(2);
    bool ec4 = p->isConstrained(3);
    bool ec5 = p->isConstrained(4);
    bool ec6 = p->isConstrained(5);
    
    // "p1" and "p4" are reflex  vertices, and either the diagonals p1p5
    // and p2p4  are inside "p" or  the diagonals p1p3 and  p4p6 are. In
    // either case,  one Steiner point, "sp", suffices  to subdivide the
    // hexagon into  a convex quadrilateral  and another hexagon  of the
    // type  RCRCCC-1,  which  can  be convex  quadrangulate  using  one
    // Steiner point.
    
    std::list<Point2> li;
    std::list<Point2> lo;
    
    li.push_back(p1.getPoint()); li.push_back(p2.getPoint());
    li.push_back(p3.getPoint()); li.push_back(p4.getPoint());
    li.push_back(p5.getPoint()); li.push_back(p6.getPoint());
    
    // find  out which  pair of  diagonals, (p1p5,p2p4)  or (p1p3,p4p6),
    // contains internal diagonals.
    
    // check if  "p1p5" is internal  by finding out on  which half-space
    // defined by the supporting line p1p5 the point "p4" lies.
    Line2 l(p1.getPoint(),p5.getPoint());
    if (l.oriented_side(p4.getPoint()) != CGAL::ON_NEGATIVE_SIDE) {
        // "p4" lies on the left half-space defined by the supporting line
        // p1p5, so diagonals p1p3 and p4p6 are inside "p".
        
        try {
            // compute intersection of the wedge of "p1" and "p"
            halfSpaceIntersection(li,p2.getPoint(),p1.getPoint(),lo,
                                  CGAL::ON_POSITIVE_SIDE);
            li.clear();
            halfSpaceIntersection(lo,p6.getPoint(),p1.getPoint(),li,
                                  CGAL::ON_NEGATIVE_SIDE);
            lo.clear();
            
            // compute intersection  of the resulting polygon  above and the
            // left half-space defined by the supporting line p1p3.
            halfSpaceIntersection(li,p1.getPoint(),p3.getPoint(),lo,
                                  CGAL::ON_NEGATIVE_SIDE);
            li.clear();
            
            // compute intersection  of the resulting polygon  above and the
            // left half-space defined by the supporting line p4p6.
            halfSpaceIntersection(lo,p4.getPoint(),p6.getPoint(),li,
                                  CGAL::ON_NEGATIVE_SIDE);
            lo.clear();
        }
        catch (const tExceptionObject& xpt) {
            treatException(xpt);
            exit(0);
        }
        
        // compute the barycenter of the resulting polygon
        Point2 q;
        double x = 0;
        double y = 0;
        int num = 0;
        while (!li.empty()) {
            q = li.front();
            li.pop_front();
            x += q.x();
            y += q.y();
            ++num;
        }
        
        q = Point2(x / double(num), y / double(num));
        
        // create Steiner point
        tCQMIndVertex2 sp = createSteinerPoint(q);
        
        // get edge status
        bool ec1 = p->isConstrained(0);
        bool ec2 = p->isConstrained(1);
        
        // create one quadrilateral
        try {
            createQuad(p1,p2,p3,sp,ec1,ec2,false,false,lq);
        }
        catch (const tExceptionObject& xpt) {
            treatException(xpt);
            exit(0);
        }
        
        // convex quadrangulate the hexagon spp3p4p5p6p1
        quadrangulateHexRCRCCC1(sp,p3,p4,p5,p6,p1,false,ec3,ec4,ec5,ec6,
                                false,lq);
    }
    else {
        // "p4"  lies on the  right half-space  defined by  the supporting
        // line p1p5, so diagonals p1p5 and p2p4 are inside "p".
        
        try {
            // compute intersection of the wedge of "p1" and "p"
            halfSpaceIntersection(li,p2.getPoint(),p1.getPoint(),lo,
                                  CGAL::ON_POSITIVE_SIDE);
            li.clear();
            halfSpaceIntersection(lo,p6.getPoint(),p1.getPoint(),li,
                                  CGAL::ON_NEGATIVE_SIDE);
            lo.clear();
            
            // compute intersection  of the resulting polygon  above and the
            // left half-space defined by the supporting line p1p5.
            halfSpaceIntersection(li,p1.getPoint(),p5.getPoint(),lo,
                                  CGAL::ON_POSITIVE_SIDE);
            li.clear();
            
            // compute intersection  of the resulting polygon  above and the
            // left half-space defined by the supporting line p2p4.
            halfSpaceIntersection(lo,p2.getPoint(),p4.getPoint(),li,
                                  CGAL::ON_NEGATIVE_SIDE);
            lo.clear();
        }
        catch (const tExceptionObject& xpt) {
            treatException(xpt);
            exit(0);
        }
        
        // compute the barycenter of the resulting polygon
        Point2 q;
        double x = 0;
        double y = 0;
        int num = 0;
        while (!li.empty()) {
            q = li.front();
            li.pop_front();
            x += q.x();
            y += q.y();
            ++num;
        }
        
        q = Point2(x / double(num), y / double(num));
        
        // create Steiner point
        tCQMIndVertex2 sp = createSteinerPoint(q);
        
        // get edge status
        bool ec5 = p->isConstrained(4);
        bool ec6 = p->isConstrained(5);
        
        // create one quadrilateral
        try {
            createQuad(p1,sp,p5,p6,false,false,ec5,ec6,lq);
        }
        catch (const tExceptionObject& xpt) {
            treatException(xpt);
            exit(0);
        }
        
        // convex quadrangulate the hexagon p4p5spp1p2p3
        quadrangulateHexRCRCCC1(p4,p5,sp,p1,p2,p3,ec4,false,false,ec1,ec2,
                                ec3,lq);
    }
}

// -------------------------------------------------------------------
// Method quadrangulateHexRCRCRC3()
// -------------------------------------------------------------------
void
tCQMCompQuad::quadrangulateHexRCRCRC3(
                                      const tCQMIndVertex2& p1,
                                      const tCQMIndVertex2& p2,
                                      const tCQMIndVertex2& p3,
                                      const tCQMIndVertex2& p4,
                                      const tCQMIndVertex2& p5,
                                      const tCQMIndVertex2& p6,
                                      bool ec1,
                                      bool ec2,
                                      bool ec3,
                                      bool ec4,
                                      bool ec5,
                                      bool ec6,
                                      std::list<tCQMQuadrilateral2*>* lq
                                      )
{
    // we assume that  "p" is not star-shaped or  the triangle p1p3p5 is
    // not inside "p".
    
    // find out a diagonal of  the triangle p1p3p5 that is obstructed by
    // an edge of "p".
    
    tCQMIndVertex2 q1 = p1; tCQMIndVertex2 q2 = p2;
    tCQMIndVertex2 q3 = p3; tCQMIndVertex2 q4 = p4;
    tCQMIndVertex2 q5 = p5; tCQMIndVertex2 q6 = p6;
    
    bool eq1 = ec1; bool eq2 = ec2; bool eq3 = ec3;
    bool eq4 = ec4; bool eq5 = ec5; bool eq6 = ec6;
    
    // verify if diagonal ae is obstructed. If so, shift the vertices of
    // the polygon such that diagonal ae becomes ac.
    Line2 l1(q1.getPoint(),q5.getPoint());
    if ((l1.oriented_side(q3.getPoint()) != CGAL::ON_NEGATIVE_SIDE) &&
        (l1.oriented_side(q2.getPoint()) != CGAL::ON_POSITIVE_SIDE) &&
        (l1.oriented_side(q4.getPoint()) != CGAL::ON_POSITIVE_SIDE)) {
        tCQMIndVertex2 q1aux = q5;
        tCQMIndVertex2 q2aux = q6;
        q5 = q1;
        q6 = q2;
        q1 = q3;
        q2 = q4;
        q3 = q1aux;
        q4 = q2aux;
        
        eq5 = eq1;
        eq6 = eq2;
        eq1 = eq3;
        eq2 = eq4;
        eq3 = ec5;
        eq4 = ec6;
    }
    
    // diagonal ae is not obstructed by any edge in the polygon
    
    std::list<Point2> li;
    std::list<Point2> lo;
    
    // compute intersection of  "p", wedge of p1, wedge  of p5 and right
    // half-space defined by the supporting line p1p5.
    li.push_back(q1.getPoint()); li.push_back(q2.getPoint());
    li.push_back(q3.getPoint()); li.push_back(q4.getPoint());
    li.push_back(q5.getPoint()); li.push_back(q6.getPoint());
    
    
    try {
        // compute intersection of "p" and wedge of p1
        halfSpaceIntersection(li,q2.getPoint(),q1.getPoint(),lo,
                              CGAL::ON_POSITIVE_SIDE);
        li.clear();
        halfSpaceIntersection(lo,q6.getPoint(),q1.getPoint(),li,
                              CGAL::ON_NEGATIVE_SIDE);
        lo.clear();
        
        // compute intersection of the result of the previous intersection
        // and wedge of p5.
        halfSpaceIntersection(li,q6.getPoint(),q5.getPoint(),lo,
                              CGAL::ON_POSITIVE_SIDE);
        li.clear();
        halfSpaceIntersection(lo,q4.getPoint(),q5.getPoint(),li,
                              CGAL::ON_NEGATIVE_SIDE);
        lo.clear();
        
        // compute intersection of the result of the previous intersection
        // and the right half-space defined by the supporting line p1p5.
        halfSpaceIntersection(li,q1.getPoint(),q5.getPoint(),lo,
                              CGAL::ON_POSITIVE_SIDE);
        li.clear();
    }
    catch (const tExceptionObject& xpt) {
        treatException(xpt);
        exit(0);
    }
    
    // compute location of first Steiner point
    Point2 q;
    double x = 0;
    double y = 0;
    int num = 0;
    while (!lo.empty()) {
        q = lo.front();
        lo.pop_front();
        x += q.x();
        y += q.y();
        ++num;
    }
    
    q = Point2(x / double(num), y / double(num));
    
    // create Steiner point
    tCQMIndVertex2 sp = createSteinerPoint(q);
    
    // create one quadrilateral
    try {
        createQuad(q1,sp,q5,q6,false,false,eq5,eq6,lq);
    }
    catch (const tExceptionObject& xpt) {
        treatException(xpt);
        exit(0);
    }
    
    // convex quadrangulate the hexagon q1q2q3q4q5sp
    
    // create hexagon
    tCQMPolygon2* paux = (tCQMPolygon2 *) new tCQMPolygon2(6);
    
    paux->insert(0,sp);
    paux->insert(1,q1);
    paux->insert(2,q2);
    paux->insert(3,q3);
    paux->insert(4,q4);
    paux->insert(5,q5);
    
    paux->setConstraint(0,false);
    paux->setConstraint(1,eq1);
    paux->setConstraint(2,eq2);
    paux->setConstraint(3,eq3);
    paux->setConstraint(4,eq4);
    paux->setConstraint(5,false);
    
    // convex quadrangulate it
    quadrangulateHexRCCRCC(paux,lq);
    
    // delete hexagon
    delete paux;
}

// -------------------------------------------------------------------
// Method quadrangulateHexRCRCRC()
// -------------------------------------------------------------------
void
tCQMCompQuad::quadrangulateHexRCRCRC(
                                     tCQMPolygon2* p,
                                     std::list<tCQMQuadrilateral2*>* lq
                                     )
{
    // get all six vertices of the hexagon
    tCQMIndVertex2 p1 = p->getIndVertex(0);
    tCQMIndVertex2 p2 = p->getIndVertex(1);
    tCQMIndVertex2 p3 = p->getIndVertex(2);
    tCQMIndVertex2 p4 = p->getIndVertex(3);
    tCQMIndVertex2 p5 = p->getIndVertex(4);
    tCQMIndVertex2 p6 = p->getIndVertex(5);
    
    // get edge status
    bool ec1 = p->isConstrained(0);
    bool ec2 = p->isConstrained(1);
    bool ec3 = p->isConstrained(2);
    bool ec4 = p->isConstrained(3);
    bool ec5 = p->isConstrained(4);
    bool ec6 = p->isConstrained(5);
    
    // if the triangle p1p3p5 is inside "p" and "p" is star-shaped, then
    // one    Steiner   point    suffices   to    convex   quadrangulate
    // "p". Otherwise, we need 3 Steiner points.
    
    // check if p1p3p5 is inside "p"
    Line2 l(p1.getPoint(),p5.getPoint());
    if (l.oriented_side(p3.getPoint()) != CGAL::ON_NEGATIVE_SIDE) {
        // "p3" is on the left of or on the line p1p5, so p1p3p5 cannot be
        // inside "p".
        quadrangulateHexRCRCRC3(p1,p2,p3,p4,p5,p6,ec1,ec2,ec3,ec4,ec5,
                                ec6,lq);
    }
    else {
        // "p3" is on the right of the line p1p5, so p1p3p5 is inside "p"
        
        // let  us  find  out  if   "p"  is  star-shaped.  We  proceed  by
        // calculating the intersection of the wedge of "p2", the wedge of
        // "p4", the wedge of "p6" and the triangle p1p3p5.
        std::list<Point2> li;
        std::list<Point2> lo;
        
        li.push_back(p1.getPoint()); li.push_back(p3.getPoint());
        li.push_back(p5.getPoint());
        
        try {
            // compute intersection  of the wedge  of "p2" and  the triangle
            // p1p2p3.
            halfSpaceIntersection(li,p2.getPoint(),p1.getPoint(),lo,
                                  CGAL::ON_POSITIVE_SIDE);
            li.clear();
            halfSpaceIntersection(lo,p2.getPoint(),p3.getPoint(),li,
                                  CGAL::ON_NEGATIVE_SIDE);
            lo.clear();
            
            // compute intersection  of the wedge of "p4"  and the resulting
            // polygon from the previous intersection.
            halfSpaceIntersection(li,p4.getPoint(),p3.getPoint(),lo,
                                  CGAL::ON_POSITIVE_SIDE);
            li.clear();
            halfSpaceIntersection(lo,p4.getPoint(),p5.getPoint(),li,
                                  CGAL::ON_NEGATIVE_SIDE);
            lo.clear();
            
            // compute intersection  of the wedge of "p6"  and the resulting
            // polygon from the previous intersection.
            halfSpaceIntersection(li,p6.getPoint(),p5.getPoint(),lo,
                                  CGAL::ON_POSITIVE_SIDE);
            li.clear();
            halfSpaceIntersection(lo,p6.getPoint(),p1.getPoint(),li,
                                  CGAL::ON_NEGATIVE_SIDE);
            lo.clear();
        }
        catch (const tExceptionObject& xpt) {
            treatException(xpt);
            exit(0);
        }
        
        // if  "p"  is  not   star-shaped,  the  result  of  the  previous
        // intersections is the empty set.
        if (li.empty()) {
            quadrangulateHexRCRCRC3(p1,p2,p3,p4,p5,p6,ec1,ec2,ec3,ec4,ec5,
                                    ec6,lq);
        }
        else {
            Point2 q;
            double x = 0;
            double y = 0;
            int num = 0;
            while (!li.empty()) {
                q = li.front();
                li.pop_front();
                x += q.x();
                y += q.y();
                ++num;
            }
            
            q = Point2(x / double(num), y / double(num));
            
            // create Steiner point
            tCQMIndVertex2 sp = createSteinerPoint(q);
            
            // create three quadrilateral
            try {
                createQuad(p1,sp,p5,p6,false,false,ec5,ec6,lq);
                createQuad(p1,p2,p3,sp,ec1,ec2,false,false,lq);
                createQuad(sp,p3,p4,p5,false,ec3,ec4,false,lq);
            }
            catch (const tExceptionObject& xpt) {
                treatException(xpt);
                exit(0);
            }
        }
    }
}

// -------------------------------------------------------------------
// Method quadrangulateHexRRCRCC2()
// -------------------------------------------------------------------
void
tCQMCompQuad::quadrangulateHexRRCRCC2(
                                      const tCQMIndVertex2& p1,
                                      const tCQMIndVertex2& p2,
                                      const tCQMIndVertex2& p3,
                                      const tCQMIndVertex2& p4,
                                      const tCQMIndVertex2& p5,
                                      const tCQMIndVertex2& p6,
                                      bool ec1,
                                      bool ec2,
                                      bool ec3,
                                      bool ec4,
                                      bool ec5,
                                      bool ec6,
                                      std::list<tCQMQuadrilateral2*>* lq
                                      )
{
    // a pre-condition to use this method is that "p1" sees "p5.
    
    // we  can convex  quadrangulate p1p2p3p4p5p6  by using  two Steiner
    // points. The first Steiner point  is placed in the intersection of
    // the kernel  of p1p2p3p4p5p6, the right half-space  defined by the
    // supporting  line p1p5  and  the left  half-space  defined by  the
    // supporting line p2p4.
    
    // compute the kernel of p1p2p3p4p5p6
    std::list<Point2> li;
    std::list<Point2> lo;
    
    // hexagon p1p2p3p4p5p6
    li.push_back(p1.getPoint()); li.push_back(p2.getPoint());
    li.push_back(p3.getPoint()); li.push_back(p4.getPoint());
    li.push_back(p5.getPoint()); li.push_back(p6.getPoint());
    
    try {
        // intersection of the hexagon p1p2p3p4p5p6 and the wedge of p3
        halfSpaceIntersection(li,p3.getPoint(),p2.getPoint(),lo,
                              CGAL::ON_POSITIVE_SIDE);
        li.clear();
        halfSpaceIntersection(lo,p3.getPoint(),p4.getPoint(),li,
                              CGAL::ON_NEGATIVE_SIDE);
        lo.clear();
        
        // intersection of  the resulting set from  the previous operation
        // and the wedge of p6.
        halfSpaceIntersection(li,p6.getPoint(),p1.getPoint(),lo,
                              CGAL::ON_NEGATIVE_SIDE);
        li.clear();
        
        // intersection of  the resulting set from  the previous operation
        // and the left half-space defined by the supporting line p2p4.
        halfSpaceIntersection(lo,p2.getPoint(),p4.getPoint(),li,
                              CGAL::ON_NEGATIVE_SIDE);
        lo.clear();
        
        // intersection of  the resulting set from  the previous operation
        // and the right half-space defined by the supporting line p1p5.
        halfSpaceIntersection(li,p1.getPoint(),p5.getPoint(),lo,
                              CGAL::ON_POSITIVE_SIDE);
        li.clear();
    }
    catch (const tExceptionObject& xpt) {
        treatException(xpt);
        exit(0);
    }
    
    
    // compute location of first Steiner point
    Point2 q;
    double x = 0;
    double y = 0;
    int num = 0;
    while (!lo.empty()) {
        q = lo.front();
        lo.pop_front();
        x += q.x();
        y += q.y();
        ++num;
    }
    
    q = Point2(x / double(num), y / double(num));
    
    // create Steiner point
    tCQMIndVertex2 sp = createSteinerPoint(q);
    
    // create one quadrilateral
    try {
        createQuad(p1,sp,p5,p6,false,false,ec5,ec6,lq);
    }
    catch (const tExceptionObject& xpt) {
        treatException(xpt);
        exit(0);
    }
    
    // the remaining hexagon spp1p2p3p4p5 is of the type RCRCRC with the
    // triangle spp2p4 inside it, so it can be convex quadrangulate with
    // one Steiner point.
    
    // create hexagon
    tCQMPolygon2* paux = (tCQMPolygon2 *) new tCQMPolygon2(6);
    
    paux->insert(0,sp);
    paux->insert(1,p1);
    paux->insert(2,p2);
    paux->insert(3,p3);
    paux->insert(4,p4);
    paux->insert(5,p5);
    
    paux->setConstraint(0,false);
    paux->setConstraint(1,ec1);
    paux->setConstraint(2,ec2);
    paux->setConstraint(3,ec3);
    paux->setConstraint(4,ec4);
    paux->setConstraint(5,false);
    
    // convex quadrangulate it
    quadrangulateHexRCRCRC(paux,lq);
    
    // delete hexagon
    delete paux;
}

// -------------------------------------------------------------------
// Method quadrangulateHexRRCRCC3()
// -------------------------------------------------------------------
void
tCQMCompQuad::quadrangulateHexRRCRCC3(
                                      const tCQMIndVertex2& p1,
                                      const tCQMIndVertex2& p2,
                                      const tCQMIndVertex2& p3,
                                      const tCQMIndVertex2& p4,
                                      const tCQMIndVertex2& p5,
                                      const tCQMIndVertex2& p6,
                                      bool ec1,
                                      bool ec2,
                                      bool ec3,
                                      bool ec4,
                                      bool ec5,
                                      bool ec6,
                                      std::list<tCQMQuadrilateral2*>* lq
                                      )
{
    // a pre-condition for  using this method is that  "p5" does not see
    // "p1".
    
    // we insert  a Steiner  point in the  intersection of the  wedge of
    // "p5"  and the  right half-space  defined by  the  supporting line
    // p6p4.
    
    // compute wedge of "p5"
    std::list<Point2> li;
    std::list<Point2> lo;
    
    // hexagon p1p2p3p4p5p6
    li.push_back(p1.getPoint()); li.push_back(p2.getPoint());
    li.push_back(p3.getPoint()); li.push_back(p4.getPoint());
    li.push_back(p5.getPoint()); li.push_back(p6.getPoint());
    
    try {
        // intersection of the hexagon p1p2p3p4p5p6 and the wedge of "p5".
        halfSpaceIntersection(li,p5.getPoint(),p4.getPoint(),lo,
                              CGAL::ON_POSITIVE_SIDE);
        li.clear();
        
        // intersection of the wedge  of "p5" and right half-space defined
        // by the supporting line p6p4.
        halfSpaceIntersection(lo,p6.getPoint(),p4.getPoint(),li,
                              CGAL::ON_POSITIVE_SIDE);
        lo.clear();
    }
    catch (const tExceptionObject& xpt) {
        treatException(xpt);
        exit(0);
    }
    
    // compute location of first Steiner point
    Point2 q;
    double x = 0;
    double y = 0;
    int num = 0;
    while (!li.empty()) {
        q = li.front();
        li.pop_front();
        x += q.x();
        y += q.y();
        ++num;
    }
    
    q = Point2(x / double(num), y / double(num));
    
    // create Steiner point
    tCQMIndVertex2 sp = createSteinerPoint(q);
    
    // create one quadrilateral
    try {
        createQuad(sp,p4,p5,p6,false,ec4,ec5,false,lq);
    }
    catch (const tExceptionObject& xpt) {
        treatException(xpt);
        exit(0);
    }
    
    // the remaining hexagon  p1p2p3p4spp6 is of the type  RRCCRC, so it
    // can be convex quadrangulate with two Steiner points.
    
    // convex quadrangulate it
    quadrangulateHexRRCCRC2(p1,p2,p3,p4,sp,p6,ec1,ec2,ec3,false,false,
                            ec6,lq);
}

// -------------------------------------------------------------------
// Method quadrangulateHexRRCRCC()
// -------------------------------------------------------------------
void
tCQMCompQuad::quadrangulateHexRRCRCC(
                                     tCQMPolygon2* p,
                                     std::list<tCQMQuadrilateral2*>* lq
                                     )
{
    // get all six vertices of the hexagon
    tCQMIndVertex2 p1 = p->getIndVertex(0);
    tCQMIndVertex2 p2 = p->getIndVertex(1);
    tCQMIndVertex2 p3 = p->getIndVertex(2);
    tCQMIndVertex2 p4 = p->getIndVertex(3);
    tCQMIndVertex2 p5 = p->getIndVertex(4);
    tCQMIndVertex2 p6 = p->getIndVertex(5);
    
    // get edge status
    bool ec1 = p->isConstrained(0);
    bool ec2 = p->isConstrained(1);
    bool ec3 = p->isConstrained(2);
    bool ec4 = p->isConstrained(3);
    bool ec5 = p->isConstrained(4);
    bool ec6 = p->isConstrained(5);
    
    // if "p5" sees "p1", then we can convex quadrangulate "p" using two
    // Steiner points. Otherwise, 3 Steiner points suffice.
    
    // verify if "p5" sees "p1", which is done by finding out which side
    // of the supporting line p1p5 "p4" lies on.
    Line2 l(p1.getPoint(),p5.getPoint());
    if (l.oriented_side(p4.getPoint()) != CGAL::ON_NEGATIVE_SIDE) {
        // "p4" lies  on the supporting line  p1p5 or on its  left side so
        // that  we can convex  quadrangulate "p"  by using  three Steiner
        // points.
        quadrangulateHexRRCRCC3(p1,p2,p3,p4,p5,p6,ec1,ec2,ec3,ec4,ec5,
                                ec6,lq);
        
    }
    else {
        // "p4" lies  on the  right sidfe of  the supporting line  p1p5 so
        // that  we can  convex  quadrangulate "p"  by  using two  Steiner
        // points.
        quadrangulateHexRRCRCC2(p1,p2,p3,p4,p5,p6,ec1,ec2,ec3,ec4,ec5,
                                ec6,lq);
    }
}

// -------------------------------------------------------------------
// Method quadrangulateHexRRCCRC2()
// -------------------------------------------------------------------
void
tCQMCompQuad::quadrangulateHexRRCCRC2(
                                      const tCQMIndVertex2& p1,
                                      const tCQMIndVertex2& p2,
                                      const tCQMIndVertex2& p3,
                                      const tCQMIndVertex2& p4,
                                      const tCQMIndVertex2& p5,
                                      const tCQMIndVertex2& p6,
                                      bool ec1,
                                      bool ec2,
                                      bool ec3,
                                      bool ec4,
                                      bool ec5,
                                      bool ec6,
                                      std::list<tCQMQuadrilateral2*>* lq
                                      )
{
    // a pre-condition to use this method is that "p1" sees "p4".
    
    // we  can convex  quadrangulate p1p2p3p4p5p6  by using  two Steiner
    // points. The first Steiner point  is placed in the intersection of
    // the kernel  of p1p2p3p4p5p6, the  left half-space defined  by the
    // supporting  line p1p4  and the  right half-space  defined  by the
    // supporting line p6p1.
    
    // compute the kernel of p1p2p3p4p5p6
    std::list<Point2> li;
    std::list<Point2> lo;
    
    // hexagon p1p2p3p4p5p6
    li.push_back(p1.getPoint()); li.push_back(p2.getPoint());
    li.push_back(p3.getPoint()); li.push_back(p4.getPoint());
    li.push_back(p5.getPoint()); li.push_back(p6.getPoint());
    
    try {
        // intersection of the hexagon p1p2p3p4p5p6 and the wedge of p6
        halfSpaceIntersection(li,p6.getPoint(),p1.getPoint(),lo,
                              CGAL::ON_NEGATIVE_SIDE);
        li.clear();
        halfSpaceIntersection(lo,p6.getPoint(),p5.getPoint(),li,
                              CGAL::ON_POSITIVE_SIDE);
        lo.clear();
        
        // intersection of  the resulting set from  the previous operation
        // and the wedge of p3.
        halfSpaceIntersection(li,p3.getPoint(),p2.getPoint(),lo,
                              CGAL::ON_POSITIVE_SIDE);
        li.clear();
        
        // intersection of  the resulting set from  the previous operation
        // and the left half-space defined by the supporting line p2p4.
        halfSpaceIntersection(lo,p2.getPoint(),p4.getPoint(),li,
                              CGAL::ON_NEGATIVE_SIDE);
        lo.clear();
        
        // intersection of  the resulting set from  the previous operation
        // and the right half-space defined by the supporting line p1p5.
        halfSpaceIntersection(li,p1.getPoint(),p5.getPoint(),lo,
                              CGAL::ON_POSITIVE_SIDE);
        li.clear();
    }
    catch (const tExceptionObject& xpt) {
        treatException(xpt);
        exit(0);
    }
    
    // compute location of first Steiner point
    Point2 q;
    double x = 0;
    double y = 0;
    int num = 0;
    while (!lo.empty()) {
        q = lo.front();
        lo.pop_front();
        x += q.x();
        y += q.y();
        ++num;
    }
    
    q = Point2(x / double(num), y / double(num));
    
    // create Steiner point
    tCQMIndVertex2 sp = createSteinerPoint(q);
    
    // create one quadrilateral
    try {
        createQuad(p2,p3,p4,sp,ec2,ec3,false,false,lq);
    }
    catch (const tExceptionObject& xpt) {
        treatException(xpt);
        exit(0);
    }
    
    // the remaining hexagon p1p2spp4p5p6 is of the type RCRCRC with the
    // triangle p1spp5 inside it, so it can be convex quadrangulate with
    // one Steiner point.
    
    // create hexagon
    tCQMPolygon2* paux = (tCQMPolygon2 *) new tCQMPolygon2(6);
    
    paux->insert(0,p1);
    paux->insert(1,p2);
    paux->insert(2,sp);
    paux->insert(3,p4);
    paux->insert(4,p5);
    paux->insert(5,p6);
    
    paux->setConstraint(0,ec1);
    paux->setConstraint(1,false);
    paux->setConstraint(2,false);
    paux->setConstraint(3,ec4);
    paux->setConstraint(4,ec5);
    paux->setConstraint(5,ec6);
    
    // convex quadrangulate it
    quadrangulateHexRCRCRC(paux,lq);
    
    // delete hexagon
    delete paux;
}

// -------------------------------------------------------------------
// Method quadrangulateHexRRCCRC3()
// -------------------------------------------------------------------
void
tCQMCompQuad::quadrangulateHexRRCCRC3(
                                      const tCQMIndVertex2& p1,
                                      const tCQMIndVertex2& p2,
                                      const tCQMIndVertex2& p3,
                                      const tCQMIndVertex2& p4,
                                      const tCQMIndVertex2& p5,
                                      const tCQMIndVertex2& p6,
                                      bool ec1,
                                      bool ec2,
                                      bool ec3,
                                      bool ec4,
                                      bool ec5,
                                      bool ec6,
                                      std::list<tCQMQuadrilateral2*>* lq
                                      )
{
    // a pre-condition for  using this method is that  "p4" does not see
    // "p2".
    
    // we insert  a Steiner  point in the  intersection of the  wedge of
    // "p4" and the left half-space defined by the supporting line p3p5.
    
    // compute wedge of "p4"
    std::list<Point2> li;
    std::list<Point2> lo;
    
    // hexagon p1p2p3p4p5p6
    li.push_back(p1.getPoint()); li.push_back(p2.getPoint());
    li.push_back(p3.getPoint()); li.push_back(p4.getPoint());
    li.push_back(p5.getPoint()); li.push_back(p6.getPoint());
    
    try {
        // intersection of the hexagon p1p2p3p4p5p6 and the wedge of "p4"
        halfSpaceIntersection(li,p4.getPoint(),p5.getPoint(),lo,
                              CGAL::ON_NEGATIVE_SIDE);
        li.clear();
        
        // intersection of  the wedge of "p4" and  left half-space defined
        // by the supporting line p3p5.
        halfSpaceIntersection(lo,p3.getPoint(),p5.getPoint(),li,
                              CGAL::ON_NEGATIVE_SIDE);
        lo.clear();
    }
    catch (const tExceptionObject& xpt) {
        treatException(xpt);
        exit(0);
    }
    
    // compute location of first Steiner point
    Point2 q;
    double x = 0;
    double y = 0;
    int num = 0;
    while (!li.empty()) {
        q = li.front();
        li.pop_front();
        x += q.x();
        y += q.y();
        ++num;
    }
    
    q = Point2(x / double(num), y / double(num));
    
    // create Steiner point
    tCQMIndVertex2 sp = createSteinerPoint(q);
    
    // create one quadrilateral
    try {
        createQuad(sp,p3,p4,p5,false,ec3,ec4,false,lq);
    }
    catch (const tExceptionObject& xpt) {
        treatException(xpt);
        exit(0);
    }
    
    // the remaining hexagon  p1p2p3spp5p6 is of the type  RRCRCC, so it
    // can be convex quadrangulate with two Steiner points.
    
    // convex quadrangulate it
    quadrangulateHexRRCRCC2(p1,p2,p3,sp,p5,p6,ec1,ec2,false,false,ec5,
                            ec6,lq);
}

// -------------------------------------------------------------------
// Method quadrangulateHexRRCRCC()
// -------------------------------------------------------------------
void
tCQMCompQuad::quadrangulateHexRRCCRC(
                                     tCQMPolygon2* p,
                                     std::list<tCQMQuadrilateral2*>* lq
                                     )
{
    // get all six vertices of the hexagon
    tCQMIndVertex2 p1 = p->getIndVertex(0);
    tCQMIndVertex2 p2 = p->getIndVertex(1);
    tCQMIndVertex2 p3 = p->getIndVertex(2);
    tCQMIndVertex2 p4 = p->getIndVertex(3);
    tCQMIndVertex2 p5 = p->getIndVertex(4);
    tCQMIndVertex2 p6 = p->getIndVertex(5);
    
    // get edge status
    bool ec1 = p->isConstrained(0);
    bool ec2 = p->isConstrained(1);
    bool ec3 = p->isConstrained(2);
    bool ec4 = p->isConstrained(3);
    bool ec5 = p->isConstrained(4);
    bool ec6 = p->isConstrained(5);
    
    // if "p4" sees "p1", then we can convex quadrangulate "p" using two
    // Steiner points. Otherwise, 3 Steiner points suffice.
    
    // verify if "p4" sees "p1", which is done by finding out which side
    // of the supporting line p1p4 "p5" lies on.
    Line2 l(p1.getPoint(),p4.getPoint());
    if (l.oriented_side(p5.getPoint()) != CGAL::ON_POSITIVE_SIDE) {
        // "p5" lies  on the supporting line  p1p4 or on its  left side so
        // that  we can convex  quadrangulate "p"  by using  three Steiner
        // points.
        quadrangulateHexRRCCRC3(p1,p2,p3,p4,p5,p6,ec1,ec2,ec3,ec4,ec5,
                                ec6,lq);
        
    }
    else {
        // "p5" lies on the right side of the supporting line p1p4 so that
        // we can convex quadrangulate "p" by using two Steiner points.
        quadrangulateHexRRCCRC2(p1,p2,p3,p4,p5,p6,ec1,ec2,ec3,ec4,ec5,
                                ec6,lq);
    }
}

// -------------------------------------------------------------------
// Method quadrangulateHexRRRCCC()
// -------------------------------------------------------------------
void
tCQMCompQuad::quadrangulateHexRRRCCC(
                                     tCQMPolygon2* p,
                                     std::list<tCQMQuadrilateral2*>* lq
                                     )
{
    // "p1,  "p2"   and  "p3"  are   reflex  vertices,  and   we  convex
    // quadrangulate "p" using three Steiner points.
    
    // get all six vertices of the hexagon
    tCQMIndVertex2 p1 = p->getIndVertex(0);
    tCQMIndVertex2 p2 = p->getIndVertex(1);
    tCQMIndVertex2 p3 = p->getIndVertex(2);
    tCQMIndVertex2 p4 = p->getIndVertex(3);
    tCQMIndVertex2 p5 = p->getIndVertex(4);
    tCQMIndVertex2 p6 = p->getIndVertex(5);
    
    // get edge status
    bool ec1 = p->isConstrained(0);
    bool ec2 = p->isConstrained(1);
    bool ec3 = p->isConstrained(2);
    bool ec4 = p->isConstrained(3);
    bool ec5 = p->isConstrained(4);
    bool ec6 = p->isConstrained(5);
    
    // First Steiner point is placed in the intersection of the wedge of
    // "p1", the  right half-space defined  by the supporting  line p1p5
    // and the left half-space defined by the supporting line p2p5.
    
    // compute position of first Steiner point
    std::list<Point2> li;
    std::list<Point2> lo;
    
    // triangle p1p2p5
    li.push_back(p1.getPoint()); li.push_back(p2.getPoint());
    li.push_back(p5.getPoint());
    
    try {
        // intersection of triangle p1p2p5 and the supporting line p6p1
        halfSpaceIntersection(li,p6.getPoint(),p1.getPoint(),lo,
                              CGAL::ON_NEGATIVE_SIDE);
        li.clear();
    }
    catch (const tExceptionObject& xpt) {
        treatException(xpt);
        exit(0);
    }
    
    Point2 q;
    double x = 0;
    double y = 0;
    int num = 0;
    while (!lo.empty()) {
        q = lo.front();
        lo.pop_front();
        x += q.x();
        y += q.y();
        ++num;
    }
    
    q = Point2(x / double(num), y / double(num));
    
    // create Steiner point
    tCQMIndVertex2 sp = createSteinerPoint(q);
    
    // create one quadrilateral
    try {
        createQuad(p1,sp,p5,p6,false,false,ec5,ec6,lq);
    }
    catch (const tExceptionObject& xpt) {
        treatException(xpt);
        exit(0);
    }
    
    // the remaining hexagon  p2p3p4p5spp1 is of the type  RRCCRC, so it
    // can be convex quadrangulate with two Steiner points.
    
    // convex quadrangulate it
    quadrangulateHexRRCCRC2(p2,p3,p4,p5,sp,p1,ec2,ec3,ec4,false,false,
                            ec1,lq);
}

// -------------------------------------------------------------------
// Method quadrangulateSep1()
// -------------------------------------------------------------------
tCQMIndVertex2
tCQMCompQuad::quadrangulateSep1(
                                const tCQMIndVertex2& p1,
                                const tCQMIndVertex2& p2,
                                const tCQMIndVertex2& p3,
                                const tCQMIndVertex2& p4,
                                const tCQMIndVertex2& p5,
                                const tCQMIndVertex2& p6,
                                const tCQMIndVertex2& p7,
                                bool ec1,
                                bool ec2,
                                bool ec3,
                                bool ec4,
                                bool ec5,
                                bool ec6,
                                bool ec7,
                                std::list<tCQMQuadrilateral2*>* lq
                                )
{
    // let us compute the visible  region R1 of the triangle p6p1p5 from
    // triangles p1p4p5 and p1p6p7.
    std::list<Point2> li;
    std::list<Point2> lo;
    
    // triangle p6p1p5
    li.push_back(p6.getPoint()); li.push_back(p1.getPoint());
    li.push_back(p5.getPoint());
    
    try {
        // intersection of triangle p6p1p5 and the wedge of "p7"
        halfSpaceIntersection(li,p7.getPoint(),p1.getPoint(),lo,
                              CGAL::ON_NEGATIVE_SIDE);
        li.clear();
        halfSpaceIntersection(lo,p7.getPoint(),p6.getPoint(),li,
                              CGAL::ON_POSITIVE_SIDE);
        lo.clear();
        
        // intersection of the polygon resulting from the previous operation
        // and the wedge of "p4" with respect to the triangle p1p4p5.
        halfSpaceIntersection(li,p4.getPoint(),p5.getPoint(),lo,
                              CGAL::ON_NEGATIVE_SIDE);
        li.clear();
        halfSpaceIntersection(lo,p4.getPoint(),p1.getPoint(),li,
                              CGAL::ON_POSITIVE_SIDE);
        lo.clear();
    }
    catch (const tExceptionObject& xpt) {
        treatException(xpt);
        exit(0);
    }
    
    // compute location of the first and second Steiner points
    Point2 q1, q2;
    
    // if "p5" is in the right half-space defined by the supporting line
    // p7p1, then the line p3p1 goes through the region R1.
    Line2 l1(p7.getPoint(),p1.getPoint());
    if (l1.oriented_side(p5.getPoint()) == CGAL::ON_NEGATIVE_SIDE) {
        // compute  intersection of  the right  half-space defined  by the
        // supporting line p3p1 and the region R1.
        try {
            halfSpaceIntersection(li,p3.getPoint(),p1.getPoint(),lo,
                                  CGAL::ON_POSITIVE_SIDE);
            li.clear();
        }
        catch (const tExceptionObject& xpt) {
            treatException(xpt);
            exit(0);
        }
        
        // compute location of the first Steiner point
        double x = 0;
        double y = 0;
        int num = 0;
        while (!lo.empty()) {
            q1 = lo.front();
            lo.pop_front();
            x += q1.x();
            y += q1.y();
            ++num;
        }
        
        q1 = Point2(x / double(num), y / double(num));
        
        
        // compute location of the second Steiner point
        
        // let us  compute the  visible region R2  of the  triangle p1p3p4
        // from triangles p1p2p3 and p1p4p5.
        li.push_back(p1.getPoint()); li.push_back(p3.getPoint());
        li.push_back(p4.getPoint());
        
        try {
            // intersection of triangle p1p3p4 and the wedge of "p2"
            halfSpaceIntersection(li,p2.getPoint(),p3.getPoint(),lo,
                                  CGAL::ON_NEGATIVE_SIDE);
            li.clear();
            halfSpaceIntersection(lo,p2.getPoint(),p1.getPoint(),li,
                                  CGAL::ON_POSITIVE_SIDE);
            lo.clear();
            
            // intersection  of  the  polygon  resulting from  the  previous
            // operation and the wedge of "p5".
            halfSpaceIntersection(li,p5.getPoint(),p1.getPoint(),lo,
                                  CGAL::ON_NEGATIVE_SIDE);
            li.clear();
            halfSpaceIntersection(lo,p5.getPoint(),p4.getPoint(),li,
                                  CGAL::ON_POSITIVE_SIDE);
            lo.clear();
        }
        catch (const tExceptionObject& xpt) {
            treatException(xpt);
            exit(0);
        }
        
        // computes the barycenter of R2.
        x = 0;
        y = 0;
        num = 0;
        while (!li.empty()) {
            q2 = li.front();
            li.pop_front();
            x += q2.x();
            y += q2.y();
            ++num;
        }
        
        q2 = Point2(x / double(num), y / double(num));
    }
    else {
        // compute  location of  the  first Steiner  point,  which is  the
        // barycenter of R1.
        double x = 0;
        double y = 0;
        int num = 0;
        while (!li.empty()) {
            q1 = li.front();
            li.pop_front();
            x += q1.x();
            y += q1.y();
            ++num;
        }
        
        q1 = Point2(x / double(num), y / double(num));
        
        // compute location of the second Steiner point
        
        // let us  compute the  visible region R2  of the  triangle p1p3p4
        // from triangles p1p2p3 and p1p4p5.
        li.push_back(p1.getPoint()); li.push_back(p3.getPoint());
        li.push_back(p4.getPoint());
        
        try {
            // intersection of triangle p1p3p4 and the wedge of "p2".
            halfSpaceIntersection(li,p2.getPoint(),p3.getPoint(),lo,
                                  CGAL::ON_NEGATIVE_SIDE);
            li.clear();
            halfSpaceIntersection(lo,p2.getPoint(),p1.getPoint(),li,
                                  CGAL::ON_POSITIVE_SIDE);
            lo.clear();
            
            // intersection  of  the  polygon  resulting from  the  previous
            // operation and the wedge of "p5".
            halfSpaceIntersection(li,p5.getPoint(),p1.getPoint(),lo,
                                  CGAL::ON_NEGATIVE_SIDE);
            li.clear();
            halfSpaceIntersection(lo,p5.getPoint(),p4.getPoint(),li,
                                  CGAL::ON_POSITIVE_SIDE);
            lo.clear();
        }
        catch (const tExceptionObject& xpt) {
            treatException(xpt);
            exit(0);
        }
        
        // if "p4"  is in  the left half-space  defined by  the supporting
        // line p2p1, then the line p6p1 goes through the region R2.
        Line2 l2(p2.getPoint(),p1.getPoint());
        if (l2.oriented_side(p4.getPoint()) == CGAL::ON_POSITIVE_SIDE) {
            // compute intersection  of the  left half-space defined  by the
            // supporting  line  p6p1 and  the region R2.
            try {
                halfSpaceIntersection(li,p6.getPoint(),p1.getPoint(),lo,
                                      CGAL::ON_NEGATIVE_SIDE);
                li.clear();
            }
            catch (const tExceptionObject& xpt) {
                treatException(xpt);
                exit(0);
            }
        }
        else {
            // "v5" is not in the right half-space defined by the supporting
            // line v7v1 and  "v4" is not in the  left half-space defined by
            // the supporting line v2v1,  so the second Steiner point should
            // be  placed in  the intersection  of triangle  p1p3p4  and the
            // wedges of  "p2" and "p5"  and the left half-space  defined by
            // the supporting line q1v1.
            try {
                halfSpaceIntersection(li,q1,p1.getPoint(),lo,
                                      CGAL::ON_NEGATIVE_SIDE);
                li.clear();
            }
            catch (const tExceptionObject& xpt) {
                treatException(xpt);
                exit(0);
            }
        }
        
        x = 0;
        y = 0;
        num = 0;
        while (!lo.empty()) {
            q2 = lo.front();
            lo.pop_front();
            x += q2.x();
            y += q2.y();
            ++num;
        }
        
        q2 = Point2(x / double(num), y / double(num));
    }
    
    // compute location of the third and fourth Steiner points
    Line2 l2(p6.getPoint(),q1);
    Line2 l4(p4.getPoint(),q1);
    
    Line2 l3(p3.getPoint(),q2);
    Line2 l5(p5.getPoint(),q2);
    
    Segment2 s2(p1.getPoint(),p5.getPoint());
    Segment2 s3(p1.getPoint(),p4.getPoint());
    
    // compute  intersection of the  supporting line  p6q1 and  the edge
    // p1p5.
    Point2 q3, qaux;
    myAssign(q3,s2,l2);
    myAssign(qaux,s2,l4);
    
    // if "qaux" is in between "p5" and "q3", then let "q3 = quax".
    if (l4.oriented_side(q3) == CGAL::ON_POSITIVE_SIDE) {
        q3 = qaux;
    }
    
    // compute  intersection of the  supporting line  p3q2 and  the edge
    // p1p4.
    Point2 q4;
    myAssign(q4,s3,l3);
    myAssign(qaux,s3,l5);
    
    // if "qaux" is in between "p4" and "q4", then let "q4 = quax".
    if (l5.oriented_side(q4) == CGAL::ON_POSITIVE_SIDE) {
        q4 = qaux;
    }
    
    q3 = Point2((p5.getPoint().x() + q3.x()) / 2.0,
                (p5.getPoint().y() + q3.y()) / 2.0);
    q4 = Point2((p4.getPoint().x() + q4.x()) / 2.0,
                (p4.getPoint().y() + q4.y()) / 2.0);
    
    // create Steiner points
    tCQMIndVertex2 sp1 = createSteinerPoint(q1);
    tCQMIndVertex2 sp2 = createSteinerPoint(q2);
    tCQMIndVertex2 sp3 = createSteinerPoint(q3);
    tCQMIndVertex2 sp4 = createSteinerPoint(q4);
    
    // create quadrilaterals
    try {
        createQuad( p1,sp1, p6, p7,false,false,ec6,ec7,lq);
        createQuad(sp1,sp3, p5, p6,false,false,ec5,false,lq);
        createQuad( p1, p2, p3,sp2,ec1,ec2,false,false,lq);
        createQuad(sp2, p3, p4,sp4,false,ec3,false,false,lq);
        createQuad(sp1, p1,sp2,sp3,false,false,false,false,lq);
        createQuad(sp3,sp2,sp4, p5,false,false,false,false,lq);
    }
    catch (const tExceptionObject& xpt) {
        treatException(xpt);
        exit(0);
    }
    
    // set Steiner point
    return sp4;
}

// -------------------------------------------------------------------
// Method quadrangulateSep2()
// -------------------------------------------------------------------
tCQMIndVertex2
tCQMCompQuad::quadrangulateSep2(
                                const tCQMIndVertex2& p1,
                                const tCQMIndVertex2& p2,
                                const tCQMIndVertex2& p3,
                                const tCQMIndVertex2& p4,
                                const tCQMIndVertex2& p5,
                                const tCQMIndVertex2& p6,
                                const tCQMIndVertex2& p7,
                                bool ec1,
                                bool ec2,
                                bool ec3,
                                bool ec4,
                                bool ec5,
                                bool ec6,
                                bool ec7,
                                std::list<tCQMQuadrilateral2*>* lq
                                )
{
    // let us compute the visible  region R1 of the triangle p1p5p6 from
    // triangles p1p6p7 and p1p4p5.
    std::list<Point2> li;
    std::list<Point2> lo;
    
    // triangle p6p1p5
    li.push_back(p6.getPoint()); li.push_back(p1.getPoint());
    li.push_back(p5.getPoint());
    
    // intersection of triangle p6p1p5 and the wedge of "p7"
    try {
        halfSpaceIntersection(li,p7.getPoint(),p1.getPoint(),lo,
                              CGAL::ON_NEGATIVE_SIDE);
        li.clear();
        halfSpaceIntersection(lo,p7.getPoint(),p6.getPoint(),li,
                              CGAL::ON_POSITIVE_SIDE);
        lo.clear();
        
        // intersection  of  the   polygon  resulting  from  the  previous
        // operation and the wedge of "p4".
        halfSpaceIntersection(li,p4.getPoint(),p5.getPoint(),lo,
                              CGAL::ON_NEGATIVE_SIDE);
        li.clear();
        halfSpaceIntersection(lo,p4.getPoint(),p1.getPoint(),li,
                              CGAL::ON_POSITIVE_SIDE);
        lo.clear();
    }
    catch (const tExceptionObject& xpt) {
        treatException(xpt);
        exit(0);
    }
    
    // compute location of the first Steiner point
    Point2 q1;
    double x = 0;
    double y = 0;
    int num = 0;
    while (!li.empty()) {
        q1 = li.front();
        li.pop_front();
        x += q1.x();
        y += q1.y();
        ++num;
    }
    
    q1 = Point2(x / double(num), y / double(num));
    
    
    // compute location of the second Steiner point
    
    // let us compute the visible  region R2 of the triangle p1p2p4 from
    // triangles p2p3p4 and p1p4p5.
    li.push_back(p1.getPoint()); li.push_back(p2.getPoint());
    li.push_back(p4.getPoint());
    
    // intersection of triangle p1p2p4 and the wedge of "p3"
    try {
        halfSpaceIntersection(li,p3.getPoint(),p4.getPoint(),lo,
                              CGAL::ON_NEGATIVE_SIDE);
        li.clear();
        halfSpaceIntersection(lo,p3.getPoint(),p2.getPoint(),li,
                              CGAL::ON_POSITIVE_SIDE);
        lo.clear();
        
        // intersection  of  the   polygon  resulting  from  the  previous
        // operation and the wedge of "p5".
        halfSpaceIntersection(li,p5.getPoint(),p1.getPoint(),lo,
                              CGAL::ON_NEGATIVE_SIDE);
        li.clear();
        halfSpaceIntersection(lo,p5.getPoint(),p4.getPoint(),li,
                              CGAL::ON_POSITIVE_SIDE);
        lo.clear();
    }
    catch (const tExceptionObject& xpt) {
        treatException(xpt);
        exit(0);
    }
    
    Point2 q2;
    x = 0;
    y = 0;
    num = 0;
    while (!li.empty()) {
        q2 = li.front();
        li.pop_front();
        x += q2.x();
        y += q2.y();
        ++num;
    }
    
    q2 = Point2(x / double(num), y / double(num));
    
    // compute location of the third Steiner point
    Line2 l2(p6.getPoint(),q1);
    Segment2 s2(p1.getPoint(),p5.getPoint());
    Point2 q3;
    myAssign(q3,s2,l2);
    q3 = Point2((p5.getPoint().x() + q3.x()) / 2.0,
                (p5.getPoint().y() + q3.y()) / 2.0);
    
    // compute location of the fourth Steiner point
    Line2 l3(q3,q2);
    Segment2 s3(p1.getPoint(),p4.getPoint());
    Point2 q4;
    myAssign(q4,s3,l3);
    q4 = Point2((p1.getPoint().x() + q4.x()) / 2.0,
                (p1.getPoint().y() + q4.y()) / 2.0);
    
    // create Steiner points
    tCQMIndVertex2 sp1 = createSteinerPoint(q1);
    tCQMIndVertex2 sp2 = createSteinerPoint(q2);
    tCQMIndVertex2 sp3 = createSteinerPoint(q3);
    tCQMIndVertex2 sp4 = createSteinerPoint(q4);
    
    // create quadrilaterals
    try {
        createQuad( p1,sp1, p6, p7,false,false,ec6,ec7,lq);
        createQuad(sp1,sp3, p5, p6,false,false,ec5,false,lq);
        createQuad( p2, p3, p4,sp2,ec2,ec3,false,false,lq);
        createQuad(sp2,sp4, p1, p2,false,false,ec1,false,lq);
        createQuad(sp1, p1,sp4,sp3,false,false,false,false,lq);
        createQuad(sp3,sp4,sp2, p5,false,false,false,false,lq);
    }
    catch (const tExceptionObject& xpt) {
        treatException(xpt);
        exit(0);
    }
    
    // set Steiner point
    return sp2;
}

// -------------------------------------------------------------------
// Method quadrangulateSep3()
// -------------------------------------------------------------------
tCQMIndVertex2
tCQMCompQuad::quadrangulateSep3(
                                const tCQMIndVertex2& p1,
                                const tCQMIndVertex2& p2,
                                const tCQMIndVertex2& p3,
                                const tCQMIndVertex2& p4,
                                const tCQMIndVertex2& p5,
                                const tCQMIndVertex2& p6,
                                const tCQMIndVertex2& p7,
                                bool ec1,
                                bool ec2,
                                bool ec3,
                                bool ec4,
                                bool ec5,
                                bool ec6,
                                bool ec7,
                                std::list<tCQMQuadrilateral2*>* lq
                                )
{
    // let us compute the visible  region R1 of the triangle p2p6p1 from
    // triangles p1p6p7 and p2p5p6.
    std::list<Point2> li;
    std::list<Point2> lo;
    
    // triangle p2p6p1
    li.push_back(p2.getPoint()); li.push_back(p6.getPoint());
    li.push_back(p1.getPoint());
    
    try {
        // intersection of triangle p2p6p1 and the wedge of "p7".
        halfSpaceIntersection(li,p7.getPoint(),p1.getPoint(),lo,
                              CGAL::ON_NEGATIVE_SIDE);
        li.clear();
        halfSpaceIntersection(lo,p7.getPoint(),p6.getPoint(),li,
                              CGAL::ON_POSITIVE_SIDE);
        lo.clear();
        
        // intersection  of  the   polygon  resulting  from  the  previous
        // operation and the wedge of "p5".
        halfSpaceIntersection(li,p5.getPoint(),p6.getPoint(),lo,
                              CGAL::ON_NEGATIVE_SIDE);
        li.clear();
        halfSpaceIntersection(lo,p5.getPoint(),p2.getPoint(),li,
                              CGAL::ON_POSITIVE_SIDE);
        lo.clear();
    }
    catch (const tExceptionObject& xpt) {
        treatException(xpt);
        exit(0);
    }
    
    // compute location of the first Steiner point
    Point2 q1;
    double x = 0;
    double y = 0;
    int num = 0;
    while (!li.empty()) {
        q1 = li.front();
        li.pop_front();
        x += q1.x();
        y += q1.y();
        ++num;
    }
    
    q1 = Point2(x / double(num), y / double(num));
    
    
    // compute location of the second Steiner point
    
    // let us compute the visible  region R2 of the triangle p2p4p5 from
    // triangles p2p3p4 and p2p5p6.
    li.push_back(p2.getPoint()); li.push_back(p4.getPoint());
    li.push_back(p5.getPoint());
    
    try {
        // intersection of triangle p1p2p4 and the wedge of "p3".
        halfSpaceIntersection(li,p3.getPoint(),p4.getPoint(),lo,
                              CGAL::ON_NEGATIVE_SIDE);
        li.clear();
        halfSpaceIntersection(lo,p3.getPoint(),p2.getPoint(),li,
                              CGAL::ON_POSITIVE_SIDE);
        lo.clear();
        
        // intersection  of  the   polygon  resulting  from  the  previous
        // operation and the wedge of "p6".
        halfSpaceIntersection(li,p6.getPoint(),p2.getPoint(),lo,
                              CGAL::ON_NEGATIVE_SIDE);
        li.clear();
        halfSpaceIntersection(lo,p6.getPoint(),p5.getPoint(),li,
                              CGAL::ON_POSITIVE_SIDE);
        lo.clear();
    }
    catch (const tExceptionObject& xpt) {
        treatException(xpt);
        exit(0);
    }
    
    Point2 q2;
    x = 0;
    y = 0;
    num = 0;
    while (!li.empty()) {
        q2 = li.front();
        li.pop_front();
        x += q2.x();
        y += q2.y();
        ++num;
    }
    
    q2 = Point2(x / double(num), y / double(num));
    
    // compute location of the fourth Steiner point
    Line2 l2(p4.getPoint(),q2);
    Segment2 s2(p2.getPoint(),p5.getPoint());
    Point2 q4;
    myAssign(q4,s2,l2);
    q4 = Point2((p5.getPoint().x() + q4.x()) / 2.0,
                (p5.getPoint().y() + q4.y()) / 2.0);
    
    // compute location of the third Steiner point
    Line2 l3(q4,q1);
    Segment2 s3(p2.getPoint(),p6.getPoint());
    Point2 q3;
    myAssign(q3,s3,l3);
    q3 = Point2((p2.getPoint().x() + q3.x()) / 2.0,
                (p2.getPoint().y() + q3.y()) / 2.0);
    
    // create Steiner points
    tCQMIndVertex2 sp1 = createSteinerPoint(q1);
    tCQMIndVertex2 sp2 = createSteinerPoint(q2);
    tCQMIndVertex2 sp3 = createSteinerPoint(q3);
    tCQMIndVertex2 sp4 = createSteinerPoint(q4);
    
    // create quadrilaterals
    try {
        createQuad( p1,sp1, p6, p7,false,false,ec6,ec7,lq);
        createQuad(sp1, p1, p2,sp3,false,ec1,false,false,lq);
        createQuad( p2, p3, p4,sp2,ec2,ec3,false,false,lq);
        createQuad(sp2, p4, p5,sp4,false,ec4,false,false,lq);
        createQuad(sp4, p5,sp1,sp3,false,false,false,false,lq);
        createQuad( p2,sp2,sp4,sp3,false,false,false,false,lq);
    }
    catch (const tExceptionObject& xpt) {
        treatException(xpt);
        exit(0);
    }
    
    // set Steiner point
    return sp1;
}

// -------------------------------------------------------------------
// Method quadrangulateSep4()
// -------------------------------------------------------------------
tCQMIndVertex2
tCQMCompQuad::quadrangulateSep4(
                                const tCQMIndVertex2& p1,
                                const tCQMIndVertex2& p2,
                                const tCQMIndVertex2& p3,
                                const tCQMIndVertex2& p4,
                                const tCQMIndVertex2& p5,
                                const tCQMIndVertex2& p6,
                                const tCQMIndVertex2& p7,
                                bool ec1,
                                bool ec2,
                                bool ec3,
                                bool ec4,
                                bool ec5,
                                bool ec6,
                                bool ec7,
                                std::list<tCQMQuadrilateral2*>* lq
                                )
{
    // let us compute the visible  region R1 of the triangle p1p2p6 from
    // triangles p1p6p7 and p2p5p6.
    std::list<Point2> li;
    std::list<Point2> lo;
    
    // triangle p1p2p6
    li.push_back(p1.getPoint()); li.push_back(p2.getPoint());
    li.push_back(p6.getPoint());
    
    try {
        // intersection of triangle p1p2p6 and the wedge of "p7".
        halfSpaceIntersection(li,p7.getPoint(),p1.getPoint(),lo,
                              CGAL::ON_NEGATIVE_SIDE);
        li.clear();
        halfSpaceIntersection(lo,p7.getPoint(),p6.getPoint(),li,
                              CGAL::ON_POSITIVE_SIDE);
        lo.clear();
        
        // intersection  of  the   polygon  resulting  from  the  previous
        // operation and the wedge of "p5".
        halfSpaceIntersection(li,p5.getPoint(),p6.getPoint(),lo,
                              CGAL::ON_NEGATIVE_SIDE);
        li.clear();
        halfSpaceIntersection(lo,p5.getPoint(),p2.getPoint(),li,
                              CGAL::ON_POSITIVE_SIDE);
        lo.clear();
    }
    catch (const tExceptionObject& xpt) {
        treatException(xpt);
        exit(0);
    }
    
    // compute location of the first Steiner point
    Point2 q1;
    double x = 0;
    double y = 0;
    int num = 0;
    while (!li.empty()) {
        q1 = li.front();
        li.pop_front();
        x += q1.x();
        y += q1.y();
        ++num;
    }
    
    q1 = Point2(x / double(num), y / double(num));
    
    
    // compute location of the second Steiner point
    
    try {
        // let us  compute the  visible region R2  of the  triangle p2p3p5
        // from triangles p3p4p5 and p2p5p6.
        li.push_back(p2.getPoint()); li.push_back(p3.getPoint());
        li.push_back(p5.getPoint());
        
        // intersection of triangle p2p3p5 and the wedge of "p4"
        halfSpaceIntersection(li,p4.getPoint(),p5.getPoint(),lo,
                              CGAL::ON_NEGATIVE_SIDE);
        li.clear();
        halfSpaceIntersection(lo,p4.getPoint(),p3.getPoint(),li,
                              CGAL::ON_POSITIVE_SIDE);
        lo.clear();
        
        // intersection  of  the   polygon  resulting  from  the  previous
        // operation and the wedge of "p6".
        halfSpaceIntersection(li,p6.getPoint(),p2.getPoint(),lo,
                              CGAL::ON_NEGATIVE_SIDE);
        li.clear();
        halfSpaceIntersection(lo,p6.getPoint(),p5.getPoint(),li,
                              CGAL::ON_POSITIVE_SIDE);
        lo.clear();
    }
    catch (const tExceptionObject& xpt) {
        treatException(xpt);
        exit(0);
    }
    
    Point2 q2;
    x = 0;
    y = 0;
    num = 0;
    while (!li.empty()) {
        q2 = li.front();
        li.pop_front();
        x += q2.x();
        y += q2.y();
        ++num;
    }
    
    q2 = Point2(x / double(num), y / double(num));
    
    // compute location of the third Steiner point
    
    // compute the  polygonal region R3  defined by the  intersection of
    // triangle  p2p5p6,  right half-spaces  defined  by the  supporting
    // lines  p1q1  and  p2p3,  and  left  half-spaces  defined  by  the
    // supporting lines p3q2 and p1p2.
    
    li.push_back(p2.getPoint()); li.push_back(p5.getPoint());
    li.push_back(p6.getPoint());
    
    try {
        halfSpaceIntersection(li,p1.getPoint(),q1,lo,
                              CGAL::ON_POSITIVE_SIDE);
        li.clear();
        halfSpaceIntersection(lo,p3.getPoint(),p2.getPoint(),li,
                              CGAL::ON_POSITIVE_SIDE);
        lo.clear();
        
        halfSpaceIntersection(li,p3.getPoint(),q2,lo,
                              CGAL::ON_NEGATIVE_SIDE);
        li.clear();
        halfSpaceIntersection(lo,p1.getPoint(),p2.getPoint(),li,
                              CGAL::ON_NEGATIVE_SIDE);
        lo.clear();
    }
    catch (const tExceptionObject& xpt) {
        treatException(xpt);
        exit(0);
    }
    
    Point2 q3;
    x = 0;
    y = 0;
    num = 0;
    while (!li.empty()) {
        q3 = li.front();
        li.pop_front();
        x += q3.x();
        y += q3.y();
        ++num;
    }
    
    q3 = Point2(x / double(num), y / double(num));
    
    // computes the location of the fourth Steiner point
    
    // compute the barycenter of the  triangle defined by "q3", "p5" and
    // "p6"
    x = (q3.x() + p5.getPoint().x() + p6.getPoint().x()) / 3;
    y = (q3.y() + p5.getPoint().y() + p6.getPoint().y()) / 3;
    
    Point2 q4(x,y);
    
    // create Steiner points
    tCQMIndVertex2 sp1 = createSteinerPoint(q1);
    tCQMIndVertex2 sp2 = createSteinerPoint(q2);
    tCQMIndVertex2 sp3 = createSteinerPoint(q3);
    tCQMIndVertex2 sp4 = createSteinerPoint(q4);
    
    // create quadrilaterals
    try {
        createQuad( p1,sp1, p6, p7,false,false,ec6,ec7,lq);
        createQuad(sp1, p1, p2,sp3,false,ec1,false,false,lq);
        createQuad( p3, p4, p5,sp2,ec3,ec4,false,false,lq);
        createQuad(sp3, p2, p3,sp2,false,ec2,false,false,lq);
        createQuad(sp3,sp2, p5,sp4,false,false,false,false,lq);
        createQuad(sp3,sp4, p6,sp1,false,false,false,false,lq);
    }
    catch (const tExceptionObject& xpt) {
        treatException(xpt);
        exit(0);
    }
    
    // set Steiner point
    return sp4;
}

// -------------------------------------------------------------------
// Method makeQuad
// -------------------------------------------------------------------
void
tCQMCompQuad::makeQuad(
                       tCQMSpanningTreeVertex* v,
                       QuadVertexSet& vs,
                       std::list<tCQMQuadrilateral2*>* lq
                       )
{
    // "v" is supposed to be a quadrilateral
    
    if (!v->isQuadrilateral()) {
        return;
    }
    
    if (v->getQuadHalfEdge()->getFace() != v->getVertex()->getFace()) {
        return;
    }
    
    // get the half-edge of the edge containing the Steiner point
    tCQMHalfEdge2* hes = v->getQuadHalfEdge();
    
    // create quadrilateral vertices
    tCQMIndVertex2 p1 = getQuadVertex(vs, hes->getVertex());
    tCQMIndVertex2 p2 = v->getQuadSteinerPoint();
    tCQMIndVertex2 p3 = getQuadVertex(vs, hes->getNext()->getVertex());
    tCQMIndVertex2 p4 = getQuadVertex(vs, hes->getPrev()->getVertex());
    
    bool ec1 = false;
    bool ec2 = false;
    bool ec3 = hes->getNext()->getEdge()->isConstrained();
    bool ec4 = hes->getPrev()->getEdge()->isConstrained();
    
    try {
        createQuad(p1,p2,p3,p4,ec1,ec2,ec3,ec4,lq);
    }
    catch (const tExceptionObject& xpt) {
        treatException(xpt);
        exit(0);
    }
}

// -------------------------------------------------------------------
// Method makeQuad
// -------------------------------------------------------------------
void
tCQMCompQuad::makeQuad(
                       tCQMSpanningTreeVertex* v1,
                       tCQMSpanningTreeVertex* v2,
                       QuadVertexSet& vs,
                       std::list<tCQMQuadrilateral2*>* lq
                       )
{
    // get the mesh faces corresponding to the dual graph vertex in "v1"
    // and "v2".
    tCQMFace2* f1 = v1->getVertex()->getFace();
    tCQMFace2* f2 = v2->getVertex()->getFace();
    
    // get the common edge of "f1" and "f2"
    tCQMEdge2* e1 = f1->getCommonEdge(f2);
    
    // get the half-edge of e1 in f1
    tCQMHalfEdge2* he1 = e1->getHalfEdge();
    tCQMHalfEdge2* he2 = e1->getMate(he1);
    if (he1->getFace() != f1) {
        he2 = he1;
        he1 = e1->getMate(he1);
    }
    
    // get vertices of the quadrilateral
    tCQMIndVertex2 p1 = getQuadVertex(vs, he1->getVertex());
    tCQMIndVertex2 p3 = getQuadVertex(vs, he2->getVertex());
    
    // compute edge status
    bool ec1 = he2->getNext()->getEdge()->isConstrained();
    
    tCQMIndVertex2 p4;
    bool ec3, ec4;
    if (v1->hasSteinerPoint()) {
        p4 = v1->getSteinerPoint();
        ec3 = false;
        ec4 = false;
    }
    else {
        p4 = getQuadVertex(vs, he1->getPrev()->getVertex());
        ec3 = he1->getNext()->getEdge()->isConstrained();
        ec4 = he1->getPrev()->getEdge()->isConstrained();
    }
    
    tCQMIndVertex2 p2;
    bool ec2;
    if (v2->hasSteinerPoint()) {
        p2 = v2->getSteinerPoint();
        ec2 = false;
    }
    else {
        p2 = getQuadVertex(vs, he2->getPrev()->getVertex());
        ec2 = he2->getPrev()->getEdge()->isConstrained();
    }
    
    try {
        createQuad(p1,p2,p3,p4,ec1,ec2,ec3,ec4,lq);
    }
    catch (const tExceptionObject& xpt) {
        treatException(xpt);
        exit(0);
    }
}

// -------------------------------------------------------------------
// Method makeQuad
// -------------------------------------------------------------------
void
tCQMCompQuad::makeQuad(
                       tCQMSpanningTreeVertex* v1,
                       tCQMSpanningTreeVertex* v2,
                       tCQMSpanningTreeVertex* vp,
                       QuadVertexSet& vs,
                       std::list<tCQMQuadrilateral2*>* lq
                       )
throw (tExceptionObject)
{
    // get the mesh faces corresponding to the dual graph vertex in "v1"
    // "v2", and "vp".
    tCQMFace2* f1 = v1->getVertex()->getFace();
    tCQMFace2* f2 = v2->getVertex()->getFace();
    tCQMFace2* fp = vp->getVertex()->getFace();
    
    // get the common edge of "f1" and "fp"
    tCQMEdge2* e1 = fp->getCommonEdge(f1);
    
    // get the common edge of "f2" and "fp"
    tCQMEdge2* e2 = fp->getCommonEdge(f2);
    
    // get the half-edge of e1 in f1 and "fp"
    tCQMHalfEdge2* he1 = e1->getHalfEdge();
    tCQMHalfEdge2* he2 = e1->getMate(he1);
    if (he1->getFace() != f1) {
        he1 = he2;
        he2 = e1->getHalfEdge();
    }
    
    // get the half-edge of e2 in fp and "f2"
    tCQMHalfEdge2* he3 = e2->getHalfEdge();
    tCQMHalfEdge2* he4 = e2->getMate(he3);
    if (he3->getFace() != fp) {
        he3 = he4;
        he4 = e2->getHalfEdge();
    }
    
    // compute vertices of the pentagon defined by "v1", "v2" and "vp"
    tCQMIndVertex2 p1, p2, p3, p4, p5;
    bool ec1, ec2, ec3, ec4, ec5;
    tCQMHalfEdge2* hesp;
    
    // if "f1" and "f2" share an edge then we create two quadrilaterals
    if (f1->getCommonEdge(f2) != 0) {
        if (he2->getNext() == he3) {
            p1 = getQuadVertex(vs, he3->getNext()->getVertex());
            p2 = getQuadVertex(vs, he1->getNext()->getVertex());
            p3 = getQuadVertex(vs, he4->getPrev()->getVertex());
            p4 = getQuadVertex(vs, he3->getVertex());
            
            ec2 = he1->getNext()->getEdge()->isConstrained();
            ec3 = he4->getPrev()->getEdge()->isConstrained();
            
            hesp = he3->getNext();
        }
        else {
            p1 = getQuadVertex(vs, he2->getNext()->getVertex());
            p2 = getQuadVertex(vs, he4->getNext()->getVertex());
            p3 = getQuadVertex(vs, he1->getPrev()->getVertex());
            p4 = getQuadVertex(vs, he2->getVertex());
            
            ec2 = he4->getNext()->getEdge()->isConstrained();
            ec3 = he1->getPrev()->getEdge()->isConstrained();
            
            hesp = he2->getNext();
        }
        
        // compute intersection of supporting line p3p4 and edge p1p2
        Line2 l(p3.getPoint(),p4.getPoint());
        Segment2 s(p1.getPoint(),p2.getPoint());
        
        Point2 pt1;
        if (!myAssign(pt1,s,l)) {
            std::stringstream ss (std::stringstream::in | std::stringstream::out);
            ss << "makeQuad(): cannot recover from numerical instability";
            throw tExceptionObject(__FILE__,__LINE__,ss.str());
        }
        
        double x = (pt1.x() + p4.getPoint().x()) / 2;
        double y = (pt1.y() + p4.getPoint().y()) / 2;
        
        Point2 q(x,y);
        
        // create Steiner point
        tCQMIndVertex2 sp = createSteinerPoint(q);
        
        // store Steiner point in "vp"
        vp->setSteinerPoint(sp,hesp);
        
        // create quadrilaterals
        try {
            createQuad(p1,sp,p4,p3,false,false,false,ec3,lq);
            createQuad(p3,p4,sp,p2,false,false,false,ec2,lq);
        }
        catch (const tExceptionObject& xpt) {
            treatException(xpt);
            exit(0);
        }
    }
    else {
        if (he2->getNext() == he3) {
            p1 = getQuadVertex(vs, he1->getVertex());
            ec1 = he4->getNext()->getEdge()->isConstrained();
            
            if (v2->hasSteinerPoint()) {
                p2 = v2->getSteinerPoint();
                ec2 = false;
            }
            else {
                p2 = getQuadVertex(vs, he4->getPrev()->getVertex());
                ec2 = he4->getPrev()->getEdge()->isConstrained();
            }
            
            p3 = getQuadVertex(vs, he4->getVertex());
            p4 = getQuadVertex(vs, he2->getVertex());
            
            ec3 = he3->getNext()->getEdge()->isConstrained();
            ec4 = he1->getNext()->getEdge()->isConstrained();
            
            if (v1->hasSteinerPoint()) {
                p5 = v1->getSteinerPoint();
                ec5 = false;
            }
            else {
                p5 = getQuadVertex(vs, he1->getPrev()->getVertex());
                ec5 = he1->getPrev()->getEdge()->isConstrained();
            }
            
            // half-edge  corresponding to  the  vertices that  will not  be
            // replaced by the Steiner point.
            hesp = he2->getPrev();
        }
        else {
            p1 = getQuadVertex(vs, he2->getVertex());
            ec1 =  he1->getNext()->getEdge()->isConstrained();
            
            if (v1->hasSteinerPoint()) {
                p2 = v1->getSteinerPoint();
                ec2 = false;
            }
            else {
                p2 = getQuadVertex(vs, he1->getPrev()->getVertex());
                ec2 = he1->getPrev()->getEdge()->isConstrained();
            }
            
            p3 = getQuadVertex(vs, he1->getVertex());
            p4 = getQuadVertex(vs, he3->getVertex());
            
            ec3 = he2->getNext()->getEdge()->isConstrained();
            ec4 = he4->getNext()->getEdge()->isConstrained();
            
            if (v2->hasSteinerPoint()) {
                p5 = v2->getSteinerPoint();
                ec5 = false;
            }
            else {
                p5 = getQuadVertex(vs, he4->getPrev()->getVertex());
                ec5 = he4->getPrev()->getEdge()->isConstrained();
            }
            
            // half-edge  corresponding to  the  vertices that  will not  be
            // replaced by the Steiner point.
            hesp = he2->getNext();
        }
        
        std::list<Point2> li;
        std::list<Point2> lo;
        
        li.push_back(p1.getPoint()); li.push_back(p3.getPoint());
        li.push_back(p4.getPoint());
        
        try {
            // compute  intersection of  triangle  p1p3p4 and  the wedge  of
            // "p2".
            halfSpaceIntersection(li,p2.getPoint(),p1.getPoint(),lo,
                                  CGAL::ON_POSITIVE_SIDE);
            li.clear();
            halfSpaceIntersection(lo,p2.getPoint(),p3.getPoint(),li,
                                  CGAL::ON_NEGATIVE_SIDE);
            lo.clear();
            
            // compute intersection  of the resulting polygon  above and the
            // wedge of "p5".
            halfSpaceIntersection(li,p5.getPoint(),p1.getPoint(),lo,
                                  CGAL::ON_NEGATIVE_SIDE);
            li.clear();
            halfSpaceIntersection(lo,p5.getPoint(),p4.getPoint(),li,
                                  CGAL::ON_POSITIVE_SIDE);
            lo.clear();
        }
        catch (const tExceptionObject& xpt) {
            treatException(xpt);
            exit(0);
        }
        
        // compute the barycenter of the resulting polygon
        Point2 q;
        double x = 0;
        double y = 0;
        int num = 0;
        while (!li.empty()) {
            q = li.front();
            li.pop_front();
            x += q.x();
            y += q.y();
            ++num;
        }
        
        q = Point2(x / double(num), y / double(num));
        
        // create Steiner point
        tCQMIndVertex2 sp = createSteinerPoint(q);
        
        // store Steiner point in "vp"
        vp->setSteinerPoint(sp,hesp);
        
        // create quadrilaterals
        try {
            createQuad(p1,p2,p3,sp,ec1,ec2,false,false,lq);
            createQuad(p5,p1,sp,p4,ec5,false,false,ec4,lq);
        }
        catch (const tExceptionObject& xpt) {
            treatException(xpt);
            exit(0);
        }
    }
}

// -------------------------------------------------------------------
// Method makeQuadFromDegeneratePoly()
// -------------------------------------------------------------------
void
tCQMCompQuad::makeQuadFromDegeneratePoly(
                                         tCQMSpanningTreeVertex* v1,
                                         tCQMSpanningTreeVertex* v2,
                                         tCQMSpanningTreeVertex* v3,
                                         tCQMSpanningTreeVertex* v4,
                                         tCQMSpanningTreeVertex* v5,
                                         QuadVertexSet& vs,
                                         std::list<tCQMQuadrilateral2*>* lq
                                         )
{
    // v1 and v2 are  children of v3, and v3 and v4  are children of v5,
    // and either "f1"  and "f2" share an edge or "f1"  and "f4" or "f2"
    // and "f4" does.
    
    // get the mesh faces corresponding to the dual graph vertex in "v1"
    // "v2","v3", "v4" and "v5".
    tCQMFace2* f1 = v1->getVertex()->getFace();
    tCQMFace2* f2 = v2->getVertex()->getFace();
    tCQMFace2* f3 = v3->getVertex()->getFace();
    tCQMFace2* f4 = v4->getVertex()->getFace();
    tCQMFace2* f5 = v5->getVertex()->getFace();
    
    // get the common edge of "f1" and "f3"
    tCQMEdge2* e1 = f3->getCommonEdge(f1);
    
    // get the common edge of "f2" and "f3"
    tCQMEdge2* e2 = f3->getCommonEdge(f2);
    
    // get the common edge of "f3" and "f5"
    tCQMEdge2* e3 = f5->getCommonEdge(f3);
    
    // get the common edge of "f4" and "f5"
    tCQMEdge2* e4 = f5->getCommonEdge(f4);
    
    // get the half-edges of "e1" in "f1" and "f3"
    tCQMHalfEdge2* he1 = e1->getHalfEdge();
    tCQMHalfEdge2* he2 = e1->getMate(he1);
    if (he1->getFace() != f1) {
        he2 = he1;
        he1 = e1->getMate(he1);
    }
    
    // get the half-edges of "e2" in "f2" and "f3"
    tCQMHalfEdge2* he3 = e2->getHalfEdge();
    tCQMHalfEdge2* he4 = e2->getMate(he3);
    if (he3->getFace() != f3) {
        he4 = he3;
        he3 = e2->getMate(he3);
    }
    
    // get the half-edge of "e3" in "f3" and "f5"
    tCQMHalfEdge2* he5 = e3->getHalfEdge();
    tCQMHalfEdge2* he6 = e3->getMate(he5);
    if (he5->getFace() != f3) {
        he6 = he5;
        he5 = e3->getMate(he5);
    }
    
    // get the half-edge of "e4" in "f4" and "f5"
    tCQMHalfEdge2* he7 = e4->getHalfEdge();
    tCQMHalfEdge2* he8 = e4->getMate(he7);
    if (he7->getFace() != f4) {
        he8 = he7;
        he7 = e4->getMate(he7);
    }
    
    tCQMIndVertex2 p1, p2, p3, p4, p5;
    bool ec1, ec2, ec3, ec4, ec5;
    tCQMHalfEdge2* hesp;
    
    if (f1->getCommonEdge(f2) != 0) {
        // "f1" and "f2" share an edge
        if (he2->getNext() == he3) {
            // "v1" is the left child of "v3"
            if (he6->getNext() == he8) {
                // "v3" is the left child of "v5"
                p1 = getQuadVertex(vs, he7->getNext()->getVertex());
                
                if (v4->hasSteinerPoint()) {
                    p2 = v4->getSteinerPoint();
                    ec1 = false;
                    ec2 = false;
                }
                else {
                    p2 = getQuadVertex(vs, he7->getPrev()->getVertex());
                    ec1 = he7->getNext()->getEdge()->isConstrained();
                    ec2 = he7->getPrev()->getEdge()->isConstrained();
                }
                
                p3 = getQuadVertex(vs, he8->getNext()->getVertex());
                p4 = getQuadVertex(vs, he1->getNext()->getVertex());
                p5 = getQuadVertex(vs, he4->getPrev()->getVertex());
                
                ec3 = he8->getNext()->getEdge()->isConstrained();
                ec4 = he1->getNext()->getEdge()->isConstrained();
                ec5 = he4->getPrev()->getEdge()->isConstrained();
                
                hesp = he8->getNext();
            }
            else {
                // "v3" is the right child of "v5"
                p1 = getQuadVertex(vs, he1->getNext()->getVertex());
                p2 = getQuadVertex(vs, he4->getPrev()->getVertex());
                p3 = getQuadVertex(vs, he6->getNext()->getVertex());
                p4 = getQuadVertex(vs, he7->getNext()->getVertex());
                
                if (v4->hasSteinerPoint()) {
                    p5 = v4->getSteinerPoint();
                    ec4 = false;
                    ec5 = false;
                }
                else {
                    p5 = getQuadVertex(vs, he7->getPrev()->getVertex());
                    ec4 = he7->getNext()->getEdge()->isConstrained();
                    ec5 = he7->getPrev()->getEdge()->isConstrained();
                }
                
                ec1 = he1->getNext()->getEdge()->isConstrained();
                ec2 = he4->getPrev()->getEdge()->isConstrained();
                ec3 = he6->getNext()->getEdge()->isConstrained();
                
                hesp = he6->getNext();
            }
        }
        else {
            // "v1" is the right child of "v3"
            if (he6->getNext() == he8) {
                // "v3" is the left child of "v5"
                p1 = getQuadVertex(vs, he7->getNext()->getVertex());
                
                if (v4->hasSteinerPoint()) {
                    p2 = v4->getSteinerPoint();
                    ec1 = false;
                    ec2 = false;
                }
                else {
                    p2 = getQuadVertex(vs, he7->getPrev()->getVertex());
                    ec1 = he7->getNext()->getEdge()->isConstrained();
                    ec2 = he7->getPrev()->getEdge()->isConstrained();
                }
                
                p3 = getQuadVertex(vs, he8->getNext()->getVertex());
                p4 = getQuadVertex(vs, he4->getNext()->getVertex());
                p5 = getQuadVertex(vs, he1->getPrev()->getVertex());
                
                ec3 = he8->getNext()->getEdge()->isConstrained();
                ec4 = he4->getNext()->getEdge()->isConstrained();
                ec5 = he1->getPrev()->getEdge()->isConstrained();
                
                hesp = he8->getNext();
            }
            else {
                // "v3" is the right child of "v5"
                p1 = getQuadVertex(vs, he4->getNext()->getVertex());
                p2 = getQuadVertex(vs, he1->getPrev()->getVertex());
                p3 = getQuadVertex(vs, he6->getNext()->getVertex());
                p4 = getQuadVertex(vs, he7->getNext()->getVertex());
                
                if (v4->hasSteinerPoint()) {
                    p5 = v4->getSteinerPoint();
                    ec4 = false;
                    ec5 = false;
                }
                else {
                    p5 = getQuadVertex(vs, he7->getPrev()->getVertex());
                    ec4 = he7->getNext()->getEdge()->isConstrained();
                    ec5 = he7->getPrev()->getEdge()->isConstrained();
                }
                
                ec1 = he4->getNext()->getEdge()->isConstrained();
                ec2 = he1->getPrev()->getEdge()->isConstrained();
                ec3 = he6->getNext()->getEdge()->isConstrained();
                
                hesp = he6->getNext();
            }
        }
        
        // quadrangulate the pentagon p1p2p3p4p5
        tCQMIndVertex2 sp = quadrangulateDegPentagon1(p1,p2,p3,p4,
                                                      p5,ec1,ec2,ec3,ec4,ec5,lq);
        
        v5->setSteinerPoint(sp,hesp);
    }
    else if (f1->getCommonEdge(f4) != 0) {
        // "f1" and "f4" share an edge
        
        if (he2->getNext() == he3) {
            // "v1" is the left child of "v3"
            p1 = getQuadVertex(vs, he6->getNext()->getVertex());
            p2 = getQuadVertex(vs, he7->getNext()->getVertex());
            p3 = getQuadVertex(vs, he1->getPrev()->getVertex());
            p4 = getQuadVertex(vs, he4->getNext()->getVertex());
            
            if (v2->hasSteinerPoint()) {
                p5 = v2->getSteinerPoint();
                ec4 = false;
                ec5 = false;
            }
            else {
                p5 = getQuadVertex(vs, he4->getPrev()->getVertex());
                ec4 = he4->getNext()->getEdge()->isConstrained();
                ec5 = he4->getPrev()->getEdge()->isConstrained();
            }
            
            ec1 = he6->getNext()->getEdge()->isConstrained();
            ec2 = he7->getNext()->getEdge()->isConstrained();
            ec3 = he1->getPrev()->getEdge()->isConstrained();
            
            tCQMIndVertex2 p6 = getQuadVertex(vs, he7->getVertex());
            
            tCQMIndVertex2 sp = quadrangulateDegPentagon1(p4,p5,p1,
                                                          p6,p3,ec4,ec5,false,false,ec3,lq);
            
            sp = quadrangulateDegPentagon1(p6,sp,p1,p2,p3,
                                           false,false,ec1,ec2,false,lq);
            
            v5->setSteinerPoint(sp,he6->getNext());
        }
        else {
            // "v1" is the right child of "v3"
            p1 = getQuadVertex(vs, he8->getNext()->getVertex());
            p2 = getQuadVertex(vs, he4->getNext()->getVertex());
            
            if (v2->hasSteinerPoint()) {
                p3 = v2->getSteinerPoint();
                ec2 = false;
                ec3 = false;
            }
            else {
                p3 = getQuadVertex(vs, he4->getPrev()->getVertex());
                ec2 = he4->getNext()->getEdge()->isConstrained();
                ec3 = he4->getPrev()->getEdge()->isConstrained();
            }
            
            p4 = getQuadVertex(vs, he1->getNext()->getVertex());
            p5 = getQuadVertex(vs, he7->getPrev()->getVertex());
            
            ec1 = he8->getNext()->getEdge()->isConstrained();
            ec4 = he1->getNext()->getEdge()->isConstrained();
            ec5 = he7->getPrev()->getEdge()->isConstrained();
            
            tCQMIndVertex2 p6 = getQuadVertex(vs, he8->getVertex());
            
            tCQMIndVertex2 sp = quadrangulateDegPentagon1(p4,p5,p6,
                                                          p2,p3,ec4,false,false,ec2,ec3,lq);
            
            sp = quadrangulateDegPentagon1(p6,p5,p1,p2,sp,false,ec5,
                                           ec1,false,false,lq);
            
            v5->setSteinerPoint(sp,he8->getNext());
        }
    }
    else if (f2->getCommonEdge(f4) != 0) {
        // "f2" and "f4" share an edge
        
        if (he2->getNext() == he3) {
            // "v1" is the left child of "v3"
            p1 = getQuadVertex(vs, he8->getNext()->getVertex());
            p2 = getQuadVertex(vs, he1->getNext()->getVertex());
            
            if (v1->hasSteinerPoint()) {
                p3 = v1->getSteinerPoint();
                ec2 = false;
                ec3 = false;
            }
            else {
                p3 = getQuadVertex(vs, he1->getPrev()->getVertex());
                ec2 = he1->getNext()->getEdge()->isConstrained();
                ec3 = he1->getPrev()->getEdge()->isConstrained();
            }
            
            p4 = getQuadVertex(vs, he4->getNext()->getVertex());
            p5 = getQuadVertex(vs, he7->getPrev()->getVertex());
            
            ec1 = he8->getNext()->getEdge()->isConstrained();
            ec4 = he4->getNext()->getEdge()->isConstrained();
            ec5 = he7->getPrev()->getEdge()->isConstrained();
            
            tCQMIndVertex2 p6 = getQuadVertex(vs, he4->getVertex());
            
            tCQMIndVertex2 sp = quadrangulateDegPentagon1(p4,p5,p6,
                                                          p2,p3,ec4,false,false,ec2,ec3,lq);
            
            sp = quadrangulateDegPentagon1(p6,p5,p1,p2,sp,false,ec5,
                                           ec1,false,false,lq);
            
            v5->setSteinerPoint(sp,he8->getNext());
        }
        else {
            // "v1" is the right child of "v3"
            p1 = getQuadVertex(vs, he6->getNext()->getVertex());
            p2 = getQuadVertex(vs, he7->getNext()->getVertex());
            p3 = getQuadVertex(vs, he4->getPrev()->getVertex());
            p4 = getQuadVertex(vs, he1->getNext()->getVertex());
            
            if (v1->hasSteinerPoint()) {
                p5 = v1->getSteinerPoint();
                ec4 = false;
                ec5 = false;
            }
            else {
                p5 = getQuadVertex(vs, he1->getPrev()->getVertex());
                ec4 = he1->getNext()->getEdge()->isConstrained();
                ec5 = he1->getPrev()->getEdge()->isConstrained();
            }
            
            ec1 = he6->getNext()->getEdge()->isConstrained();
            ec2 = he7->getNext()->getEdge()->isConstrained();
            ec3 = he4->getPrev()->getEdge()->isConstrained();
            
            tCQMIndVertex2 p6 = getQuadVertex(vs, he3->getVertex());
            
            tCQMIndVertex2 sp = quadrangulateDegPentagon1(p4,p5,p1,
                                                          p6,p3,ec4,ec5,false,false,ec3,lq);
            
            sp = quadrangulateDegPentagon1(p6,sp,p1,p2,p3,
                                           false,false,ec1,ec2,false,lq);
            
            v5->setSteinerPoint(sp,he6->getNext());
        }
    }
}

// -------------------------------------------------------------------
// Method makeQuadFromDegeneratePoly()
// -------------------------------------------------------------------
void
tCQMCompQuad::makeQuadFromDegeneratePoly(
                                         tCQMSpanningTreeVertex* v1,
                                         tCQMSpanningTreeVertex* v2,
                                         tCQMSpanningTreeVertex* v3,
                                         tCQMSpanningTreeVertex* v4,
                                         tCQMSpanningTreeVertex* v5,
                                         tCQMSpanningTreeVertex* v6,
                                         tCQMSpanningTreeVertex* v7,
                                         QuadVertexSet& vs,
                                         std::list<tCQMQuadrilateral2*>* lq
                                         )
{
    // "v1" and "v2" are children of "v3", "v4" and "v5" are children of
    // "v6", and "v3" and "v6" are children of "v7".
    
    // two of "v1",  "v2", "v4" and "v5" correspond  to faces that share
    // one edge.
    
    // get the mesh faces corresponding to the dual graph vertex in "v1"
    // "v2","v3", "v4", "v5", "v6", and "v7".
    tCQMFace2* f1 = v1->getVertex()->getFace();
    tCQMFace2* f2 = v2->getVertex()->getFace();
    tCQMFace2* f3 = v3->getVertex()->getFace();
    tCQMFace2* f4 = v4->getVertex()->getFace();
    tCQMFace2* f5 = v5->getVertex()->getFace();
    tCQMFace2* f6 = v6->getVertex()->getFace();
    tCQMFace2* f7 = v7->getVertex()->getFace();
    
    // get the common edge of "f1" and "f3"
    tCQMEdge2* e1 = f3->getCommonEdge(f1);
    
    // get the half-edges of "e1" in "f1" and "f3"
    tCQMHalfEdge2* he1 = e1->getHalfEdge();
    tCQMHalfEdge2* he2 = e1->getMate(he1);
    if (he1->getFace() != f1) {
        he2 = he1;
        he1 = e1->getMate(he1);
    }
    
    // get the common edge of "f2" and "f3"
    tCQMEdge2* e2 = f3->getCommonEdge(f2);
    
    // get the half-edges of "e2" in "f2" and "f3"
    tCQMHalfEdge2* he3 = e2->getHalfEdge();
    tCQMHalfEdge2* he4 = e2->getMate(he3);
    if (he3->getFace() != f2) {
        he4 = he3;
        he3 = e2->getMate(he3);
    }
    
    // get the common edge of "f3" and "f7"
    tCQMEdge2* e3 = f3->getCommonEdge(f7);
    
    // get the half-edges of "e3" in "f3" and "f7"
    tCQMHalfEdge2* he5 = e3->getHalfEdge();
    tCQMHalfEdge2* he6 = e3->getMate(he5);
    if (he5->getFace() != f3) {
        he6 = he5;
        he5 = e3->getMate(he5);
    }
    
    // get the common edge of "f4" and "f6"
    tCQMEdge2* e4 = f4->getCommonEdge(f6);
    
    // get the half-edges of "e4" in "f4" and "f6"
    tCQMHalfEdge2* he7 = e4->getHalfEdge();
    tCQMHalfEdge2* he8 = e4->getMate(he7);
    if (he7->getFace() != f4) {
        he8 = he7;
        he7 = e4->getMate(he7);
    }
    
    // get the common edge of "f5" and "f6"
    tCQMEdge2* e5 = f5->getCommonEdge(f6);
    
    // get the half-edges of "e5" in "f5" and "f6"
    tCQMHalfEdge2* he9  = e5->getHalfEdge();
    tCQMHalfEdge2* he10 = e5->getMate(he9);
    if (he9->getFace() != f5) {
        he10 = he9;
        he9  = e5->getMate(he9);
    }
    
    // get the common edge of "f6" and "f7"
    tCQMEdge2* e6 = f6->getCommonEdge(f7);
    
    // get the half-edges of "e6" in "f6" and "f7"
    tCQMHalfEdge2* he11 = e6->getHalfEdge();
    tCQMHalfEdge2* he12 = e6->getMate(he11);
    if (he11->getFace() != f6) {
        he12 = he11;
        he11 = e6->getMate(he11);
    }
    
    // we form a polygon according to the kind of degeneracy
    tCQMEdge2* e12 = f1->getCommonEdge(f2);
    tCQMEdge2* e14 = f1->getCommonEdge(f4);
    tCQMEdge2* e15 = f1->getCommonEdge(f5);
    tCQMEdge2* e24 = f2->getCommonEdge(f4);
    tCQMEdge2* e25 = f2->getCommonEdge(f5);
    tCQMEdge2* e45 = f4->getCommonEdge(f5);
    
    tCQMIndVertex2 p1, p2, p3, p4, p5, p6, p7, sp;
    bool ec1, ec2, ec3, ec4, ec5, ec6, ec7;
    tCQMHalfEdge2* hesp;
    
    if (he6->getNext() == he12) {
        // "v3" is the left child of "v7"
        if (he2->getNext() == he4) {
            // "v1" is the left child of "v3"
            if (he8->getNext() == he10) {
                // "v4" is the left child of "v6"
                if ((e12 != 0) && (e24 != 0) && (e45 != 0)) {
                    // polygon degenerated to a triangle
                    p1 = getQuadVertex(vs, he9->getPrev()->getVertex());
                    p2 = getQuadVertex(vs, he6->getPrev()->getVertex());
                    p3 = getQuadVertex(vs, he1->getNext()->getVertex());
                    p4 = getQuadVertex(vs, he8->getVertex());
                    
                    ec1 = he9->getPrev()->getEdge()->isConstrained();
                    ec3 = he1->getNext()->getEdge()->isConstrained();
                    
                    sp = quadrangulateDegTriangle(p1,p2,p3,p4,ec1,
                                                  ec3,lq);
                }
                else if ((e12 != 0) && (e24 != 0)) {
                    // polygon degenerated to a pentagon
                    p1 = getQuadVertex(vs, he9->getNext()->getVertex());
                    
                    if (v5->hasSteinerPoint()) {
                        p2 = v5->getSteinerPoint();
                        ec1 = false;
                        ec2 = false;
                    }
                    else {
                        p2 = getQuadVertex(vs, he9->getPrev()->getVertex());
                        ec1 = he9->getNext()->getEdge()->isConstrained();
                        ec2 = he9->getPrev()->getEdge()->isConstrained();
                    }
                    
                    p3 = getQuadVertex(vs, he6->getPrev()->getVertex());
                    p4 = getQuadVertex(vs, he1->getNext()->getVertex());
                    p5 = getQuadVertex(vs, he7->getPrev()->getVertex());
                    
                    ec3 = he6->getPrev()->getEdge()->isConstrained();
                    ec4 = he1->getNext()->getEdge()->isConstrained();
                    ec5 = he7->getPrev()->getEdge()->isConstrained();
                    
                    sp = quadrangulateDegPentagon1(p1,p2,p3,p4,p5,ec1,ec2,
                                                   ec3,ec4,ec5,lq);
                }
                else if ((e12 != 0) && (e45 != 0)) {
                    // polygon degenerated to a pentagon
                    p1 = getQuadVertex(vs, he7->getNext()->getVertex());
                    p2 = getQuadVertex(vs, he9->getPrev()->getVertex());
                    p3 = getQuadVertex(vs, he6->getPrev()->getVertex());
                    p4 = getQuadVertex(vs, he1->getNext()->getVertex());
                    p5 = getQuadVertex(vs, he3->getPrev()->getVertex());
                    
                    ec1 = he7->getNext()->getEdge()->isConstrained();
                    ec2 = he9->getPrev()->getEdge()->isConstrained();
                    ec3 = he6->getPrev()->getEdge()->isConstrained();
                    ec4 = he1->getNext()->getEdge()->isConstrained();
                    ec5 = he3->getPrev()->getEdge()->isConstrained();
                    
                    sp = quadrangulateDegPentagon1(p1,p2,p3,p4,p5,ec1,ec2,
                                                   ec3,ec4,ec5,lq);
                }
                else if ((e24 != 0) && (e45 != 0)) {
                    // polygon degenerated to a pentagon
                    p1 = getQuadVertex(vs, he3->getNext()->getVertex());
                    p2 = getQuadVertex(vs, he9->getPrev()->getVertex());
                    p3 = getQuadVertex(vs, he6->getPrev()->getVertex());
                    p4 = getQuadVertex(vs, he1->getNext()->getVertex());
                    
                    if (v1->hasSteinerPoint()) {
                        p5 = v1->getSteinerPoint();
                        ec4 = false;
                        ec5 = false;
                    }
                    else {
                        p5 = getQuadVertex(vs, he1->getPrev()->getVertex());
                        ec4 = he1->getNext()->getEdge()->isConstrained();
                        ec5 = he1->getPrev()->getEdge()->isConstrained();
                    }
                    
                    ec1 = he3->getNext()->getEdge()->isConstrained();
                    ec2 = he9->getPrev()->getEdge()->isConstrained();
                    ec3 = he6->getPrev()->getEdge()->isConstrained();
                    
                    sp = quadrangulateDegPentagon1(p1,p2,p3,p4,p5,ec1,ec2,
                                                   ec3,ec4,ec5,lq);
                }
                else if (e12 != 0) {
                    // polygon degenerated to a septagon
                    p1 = getQuadVertex(vs, he9->getNext()->getVertex());
                    
                    if (v5->hasSteinerPoint()) {
                        p2 = v5->getSteinerPoint();
                        ec1 = false;
                        ec2 = false;
                    }
                    else {
                        p2 = getQuadVertex(vs, he9->getPrev()->getVertex());
                        ec1 = he9->getNext()->getEdge()->isConstrained();
                        ec2 = he9->getPrev()->getEdge()->isConstrained();
                    }
                    
                    p3 = getQuadVertex(vs, he6->getPrev()->getVertex());
                    p4 = getQuadVertex(vs, he1->getNext()->getVertex());
                    p5 = getQuadVertex(vs, he3->getPrev()->getVertex());
                    
                    ec3 = he6->getPrev()->getEdge()->isConstrained();
                    ec4 = he1->getNext()->getEdge()->isConstrained();
                    ec5 = he3->getPrev()->getEdge()->isConstrained();
                    
                    p6 = getQuadVertex(vs, he7->getNext()->getVertex());
                    if (v4->hasSteinerPoint()) {
                        p7 = v4->getSteinerPoint();
                        ec6 = false;
                        ec7 = false;
                    }
                    else {
                        p7 = getQuadVertex(vs, he7->getPrev()->getVertex());
                        ec6 = he7->getNext()->getEdge()->isConstrained();
                        ec7 = he7->getPrev()->getEdge()->isConstrained();
                    }
                    
                    sp = quadrangulateDegPentagon1(p1,p2,p3,p6,p7,ec1,ec2,
                                                   false,ec6,ec7,lq);
                    sp = quadrangulateDegPentagon1(p6,sp,p3,p4,p5,false,
                                                   false,ec3,ec4,ec5,lq);
                }
                else if (e24 != 0) {
                    // polygon degenerated to a septagon
                    p1 = getQuadVertex(vs, he3->getNext()->getVertex());
                    p2 = getQuadVertex(vs, he7->getPrev()->getVertex());
                    
                    ec1 = he3->getNext()->getEdge()->isConstrained();
                    ec2 = he7->getPrev()->getEdge()->isConstrained();
                    
                    p3 = getQuadVertex(vs, he9->getNext()->getVertex());
                    
                    if (v5->hasSteinerPoint()) {
                        p4 = v5->getSteinerPoint();
                        ec3 = false;
                        ec4 = false;
                    }
                    else {
                        p4 = getQuadVertex(vs, he9->getPrev()->getVertex());
                        ec3 = he9->getNext()->getEdge()->isConstrained();
                        ec4 = he9->getPrev()->getEdge()->isConstrained();
                    }
                    
                    p5 = getQuadVertex(vs, he6->getPrev()->getVertex());
                    
                    ec5 = he6->getPrev()->getEdge()->isConstrained();
                    
                    p6 = getQuadVertex(vs, he1->getNext()->getVertex());
                    
                    if (v1->hasSteinerPoint()) {
                        p7 = v1->getSteinerPoint();
                        ec6 = false;
                        ec7 = false;
                    }
                    else {
                        p7 = getQuadVertex(vs, he1->getPrev()->getVertex());
                        ec6 = he1->getNext()->getEdge()->isConstrained();
                        ec7 = he1->getPrev()->getEdge()->isConstrained();
                    }
                    
                    tCQMIndVertex2 spaux1, spaux2;
                    
                    sp = getQuadVertex(vs, he5->getVertex());
                    spaux1 = quadrangulateDegPentagon1(p1,p2,sp,p6,p7,ec1,false,
                                                       false,ec6,ec7,lq);
                    spaux2 = quadrangulateDegPentagon1(p3,p4,p5,sp,p2,ec3,ec4,
                                                       false,false,ec2,lq);
                    sp = quadrangulateDegPentagon1(sp,spaux2,p5,p6,spaux1,false,
                                                   false,ec5,false,false,lq);
                }
                else {  // e45 != 0
                    // polygon degenerated to a septagon
                    p1 = getQuadVertex(vs, he3->getNext()->getVertex());
                    
                    if (v2->hasSteinerPoint()) {
                        p2 = v2->getSteinerPoint();
                        ec1 = false;
                        ec2 = false;
                    }
                    else {
                        p2 = getQuadVertex(vs, he3->getPrev()->getVertex());
                        ec1 = he3->getNext()->getEdge()->isConstrained();
                        ec2 = he3->getPrev()->getEdge()->isConstrained();
                    }
                    
                    p3 = getQuadVertex(vs, he7->getNext()->getVertex());
                    p4 = getQuadVertex(vs, he9->getPrev()->getVertex());
                    p5 = getQuadVertex(vs, he6->getPrev()->getVertex());
                    
                    ec3 = he7->getNext()->getEdge()->isConstrained();
                    ec4 = he9->getPrev()->getEdge()->isConstrained();
                    ec5 = he6->getPrev()->getEdge()->isConstrained();
                    
                    p6 = getQuadVertex(vs, he1->getNext()->getVertex());
                    if (v1->hasSteinerPoint()) {
                        p7 = v1->getSteinerPoint();
                        ec6 = false;
                        ec7 = false;
                    }
                    else {
                        p7 = getQuadVertex(vs, he1->getPrev()->getVertex());
                        ec6 = he1->getNext()->getEdge()->isConstrained();
                        ec7 = he1->getPrev()->getEdge()->isConstrained();
                    }
                    
                    sp = quadrangulateDegPentagon1(p1,p2,p3,p6,p7,ec1,ec2,
                                                   false,ec6,ec7,lq);
                    sp = quadrangulateDegPentagon1(p3,p4,p5,p6,sp,ec3,ec4,
                                                   ec5,false,false,lq);
                }
            }
            else {
                // "v4" is the right child of "v6"
                if ((e12 != 0) && (e25 != 0) && (e45 != 0)) {
                    // polygon degenerated to a triangle
                    p1 = getQuadVertex(vs, he7->getPrev()->getVertex());
                    p2 = getQuadVertex(vs, he6->getPrev()->getVertex());
                    p3 = getQuadVertex(vs, he1->getNext()->getVertex());
                    p4 = getQuadVertex(vs, he10->getVertex());
                    
                    ec1 = he7->getPrev()->getEdge()->isConstrained();
                    ec3 = he1->getNext()->getEdge()->isConstrained();
                    
                    sp = quadrangulateDegTriangle(p1,p2,p3,p4,ec1,
                                                  ec3,lq);
                }
                else if ((e12 != 0) && (e25 != 0)) {
                    // polygon degenerated to a pentagon
                    p1 = getQuadVertex(vs, he7->getNext()->getVertex());
                    
                    if (v4->hasSteinerPoint()) {
                        p2 = v4->getSteinerPoint();
                        ec1 = false;
                        ec2 = false;
                    }
                    else {
                        p2 = getQuadVertex(vs, he7->getPrev()->getVertex());
                        ec1 = he7->getNext()->getEdge()->isConstrained();
                        ec2 = he7->getPrev()->getEdge()->isConstrained();
                    }
                    
                    p3 = getQuadVertex(vs, he6->getPrev()->getVertex());
                    p4 = getQuadVertex(vs, he1->getNext()->getVertex());
                    p5 = getQuadVertex(vs, he9->getPrev()->getVertex());
                    
                    ec3 = he6->getPrev()->getEdge()->isConstrained();
                    ec4 = he1->getNext()->getEdge()->isConstrained();
                    ec5 = he9->getPrev()->getEdge()->isConstrained();
                    
                    sp = quadrangulateDegPentagon1(p1,p2,p3,p4,p5,ec1,ec2,
                                                   ec3,ec4,ec5,lq);
                }
                else if ((e12 != 0) && (e45 != 0)) {
                    // polygon degenerated to a pentagon
                    p1 = getQuadVertex(vs, he9->getNext()->getVertex());
                    p2 = getQuadVertex(vs, he7->getPrev()->getVertex());
                    p3 = getQuadVertex(vs, he6->getPrev()->getVertex());
                    p4 = getQuadVertex(vs, he1->getNext()->getVertex());
                    p5 = getQuadVertex(vs, he3->getPrev()->getVertex());
                    
                    ec1 = he9->getNext()->getEdge()->isConstrained();
                    ec2 = he7->getPrev()->getEdge()->isConstrained();
                    ec3 = he6->getPrev()->getEdge()->isConstrained();
                    ec4 = he1->getNext()->getEdge()->isConstrained();
                    ec5 = he3->getPrev()->getEdge()->isConstrained();
                    
                    sp = quadrangulateDegPentagon1(p1,p2,p3,p4,p5,ec1,ec2,
                                                   ec3,ec4,ec5,lq);
                }
                else if ((e25 != 0) && (e45 != 0)) {
                    // polygon degenerated to a pentagon
                    p1 = getQuadVertex(vs, he3->getNext()->getVertex());
                    p2 = getQuadVertex(vs, he7->getPrev()->getVertex());
                    p3 = getQuadVertex(vs, he6->getPrev()->getVertex());
                    p4 = getQuadVertex(vs, he1->getNext()->getVertex());
                    
                    if (v1->hasSteinerPoint()) {
                        p5 = v1->getSteinerPoint();
                        ec4 = false;
                        ec5 = false;
                    }
                    else {
                        p5 = getQuadVertex(vs, he1->getPrev()->getVertex());
                        ec4 = he1->getNext()->getEdge()->isConstrained();
                        ec5 = he1->getPrev()->getEdge()->isConstrained();
                    }
                    
                    ec1 = he3->getNext()->getEdge()->isConstrained();
                    ec2 = he7->getPrev()->getEdge()->isConstrained();
                    ec3 = he6->getPrev()->getEdge()->isConstrained();
                    
                    sp = quadrangulateDegPentagon1(p1,p2,p3,p4,p5,ec1,ec2,
                                                   ec3,ec4,ec5,lq);
                }
                else if (e12 != 0) {
                    // polygon degenerated to a septagon
                    p1 = getQuadVertex(vs, he7->getNext()->getVertex());
                    
                    if (v4->hasSteinerPoint()) {
                        p2 = v4->getSteinerPoint();
                        ec1 = false;
                        ec2 = false;
                    }
                    else {
                        p2 = getQuadVertex(vs, he7->getPrev()->getVertex());
                        ec1 = he7->getNext()->getEdge()->isConstrained();
                        ec2 = he7->getPrev()->getEdge()->isConstrained();
                    }
                    
                    p3 = getQuadVertex(vs, he6->getPrev()->getVertex());
                    p4 = getQuadVertex(vs, he1->getNext()->getVertex());
                    p5 = getQuadVertex(vs, he3->getPrev()->getVertex());
                    
                    ec3 = he6->getPrev()->getEdge()->isConstrained();
                    ec4 = he1->getNext()->getEdge()->isConstrained();
                    ec5 = he3->getPrev()->getEdge()->isConstrained();
                    
                    p6 = getQuadVertex(vs, he9->getNext()->getVertex());
                    if (v5->hasSteinerPoint()) {
                        p7 = v5->getSteinerPoint();
                        ec6 = false;
                        ec7 = false;
                    }
                    else {
                        p7 = getQuadVertex(vs, he9->getPrev()->getVertex());
                        ec6 = he9->getNext()->getEdge()->isConstrained();
                        ec7 = he9->getPrev()->getEdge()->isConstrained();
                    }
                    
                    sp = quadrangulateDegPentagon1(p1,p2,p3,p6,p7,ec1,ec2,
                                                   false,ec6,ec7,lq);
                    sp = quadrangulateDegPentagon1(p6,sp,p3,p4,p5,false,
                                                   false,ec3,ec4,ec5,lq);
                }
                else if (e25 != 0) {
                    // polygon degenerated to a septagon
                    p1 = getQuadVertex(vs, he3->getNext()->getVertex());
                    p2 = getQuadVertex(vs, he9->getPrev()->getVertex());
                    
                    ec1 = he3->getNext()->getEdge()->isConstrained();
                    ec2 = he9->getPrev()->getEdge()->isConstrained();
                    
                    p3 = getQuadVertex(vs, he7->getNext()->getVertex());
                    
                    if (v4->hasSteinerPoint()) {
                        p4 = v4->getSteinerPoint();
                        ec3 = false;
                        ec4 = false;
                    }
                    else {
                        p4 = getQuadVertex(vs, he7->getPrev()->getVertex());
                        ec3 = he7->getNext()->getEdge()->isConstrained();
                        ec4 = he7->getPrev()->getEdge()->isConstrained();
                    }
                    
                    p5 = getQuadVertex(vs, he6->getPrev()->getVertex());
                    
                    ec5 = he6->getPrev()->getEdge()->isConstrained();
                    
                    p6 = getQuadVertex(vs, he1->getNext()->getVertex());
                    
                    if (v1->hasSteinerPoint()) {
                        p7 = v1->getSteinerPoint();
                        ec6 = false;
                        ec7 = false;
                    }
                    else {
                        p7 = getQuadVertex(vs, he1->getPrev()->getVertex());
                        ec6 = he1->getNext()->getEdge()->isConstrained();
                        ec7 = he1->getPrev()->getEdge()->isConstrained();
                    }
                    
                    tCQMIndVertex2 spaux1, spaux2;
                    
                    sp = getQuadVertex(vs, he5->getVertex());
                    spaux1 = quadrangulateDegPentagon1(p1,p2,sp,p6,p7,ec1,false,
                                                       false,ec6,ec7,lq);
                    spaux2 = quadrangulateDegPentagon1(p3,p4,p5,sp,p2,ec3,ec4,
                                                       false,false,ec2,lq);
                    sp = quadrangulateDegPentagon1(sp,spaux2,p5,p6,spaux1,false,
                                                   false,ec5,false,false,lq);
                }
                else {  // e45 != 0
                    // polygon degenerated to a septagon
                    p1 = getQuadVertex(vs, he3->getNext()->getVertex());
                    
                    if (v2->hasSteinerPoint()) {
                        p2 = v2->getSteinerPoint();
                        ec1 = false;
                        ec2 = false;
                    }
                    else {
                        p2 = getQuadVertex(vs, he3->getPrev()->getVertex());
                        ec1 = he3->getNext()->getEdge()->isConstrained();
                        ec2 = he3->getPrev()->getEdge()->isConstrained();
                    }
                    
                    p3 = getQuadVertex(vs, he9->getNext()->getVertex());
                    p4 = getQuadVertex(vs, he7->getPrev()->getVertex());
                    p5 = getQuadVertex(vs, he6->getPrev()->getVertex());
                    
                    ec3 = he9->getNext()->getEdge()->isConstrained();
                    ec4 = he7->getPrev()->getEdge()->isConstrained();
                    ec5 = he6->getPrev()->getEdge()->isConstrained();
                    
                    p6 = getQuadVertex(vs, he1->getNext()->getVertex());
                    if (v1->hasSteinerPoint()) {
                        p7 = v1->getSteinerPoint();
                        ec6 = false;
                        ec7 = false;
                    }
                    else {
                        p7 = getQuadVertex(vs, he1->getPrev()->getVertex());
                        ec6 = he1->getNext()->getEdge()->isConstrained();
                        ec7 = he1->getPrev()->getEdge()->isConstrained();
                    }
                    
                    sp = quadrangulateDegPentagon1(p1,p2,p3,p6,p7,ec1,ec2,
                                                   false,ec6,ec7,lq);
                    sp = quadrangulateDegPentagon1(p3,p4,p5,p6,sp,ec3,ec4,
                                                   ec5,false,false,lq);
                }
            }
        }
        else {
            // "v1" is the right child of "v3"
            if (he8->getNext() == he10) {
                // "v4" is the left child of "v6"
                
                if ((e12 != 0) && (e14 != 0) && (e45 != 0)) {
                    // polygon degenerated to a triangle
                    p1 = getQuadVertex(vs, he9->getPrev()->getVertex());
                    p2 = getQuadVertex(vs, he6->getPrev()->getVertex());
                    p3 = getQuadVertex(vs, he3->getNext()->getVertex());
                    p4 = getQuadVertex(vs, he8->getVertex());
                    
                    ec1 = he9->getPrev()->getEdge()->isConstrained();
                    ec3 = he3->getNext()->getEdge()->isConstrained();
                    
                    sp = quadrangulateDegTriangle(p1,p2,p3,p4,ec1,
                                                  ec3,lq);
                }
                else if ((e12 != 0) && (e14 != 0)) {
                    // polygon degenerated to a pentagon
                    p1 = getQuadVertex(vs, he9->getNext()->getVertex());
                    
                    if (v5->hasSteinerPoint()) {
                        p2 = v5->getSteinerPoint();
                        ec1 = false;
                        ec2 = false;
                    }
                    else {
                        p2 = getQuadVertex(vs, he9->getPrev()->getVertex());
                        ec1 = he9->getNext()->getEdge()->isConstrained();
                        ec2 = he9->getPrev()->getEdge()->isConstrained();
                    }
                    
                    p3 = getQuadVertex(vs, he6->getPrev()->getVertex());
                    p4 = getQuadVertex(vs, he3->getNext()->getVertex());
                    p5 = getQuadVertex(vs, he7->getPrev()->getVertex());
                    
                    ec3 = he6->getPrev()->getEdge()->isConstrained();
                    ec4 = he3->getNext()->getEdge()->isConstrained();
                    ec5 = he7->getPrev()->getEdge()->isConstrained();
                    
                    sp = quadrangulateDegPentagon1(p1,p2,p3,p4,p5,ec1,ec2,
                                                   ec3,ec4,ec5,lq);
                }
                else if ((e12 != 0) && (e45 != 0)) {
                    // polygon degenerated to a pentagon
                    p1 = getQuadVertex(vs, he7->getNext()->getVertex());
                    p2 = getQuadVertex(vs, he9->getPrev()->getVertex());
                    p3 = getQuadVertex(vs, he6->getPrev()->getVertex());
                    p4 = getQuadVertex(vs, he3->getNext()->getVertex());
                    p5 = getQuadVertex(vs, he1->getPrev()->getVertex());
                    
                    ec1 = he7->getNext()->getEdge()->isConstrained();
                    ec2 = he9->getPrev()->getEdge()->isConstrained();
                    ec3 = he6->getPrev()->getEdge()->isConstrained();
                    ec4 = he3->getNext()->getEdge()->isConstrained();
                    ec5 = he1->getPrev()->getEdge()->isConstrained();
                    
                    sp = quadrangulateDegPentagon1(p1,p2,p3,p4,p5,ec1,ec2,
                                                   ec3,ec4,ec5,lq);
                }
                else if ((e14 != 0) && (e45 != 0)) {
                    // polygon degenerated to a pentagon
                    p1 = getQuadVertex(vs, he1->getNext()->getVertex());
                    p2 = getQuadVertex(vs, he9->getPrev()->getVertex());
                    p3 = getQuadVertex(vs, he6->getPrev()->getVertex());
                    p4 = getQuadVertex(vs, he3->getNext()->getVertex());
                    
                    if (v2->hasSteinerPoint()) {
                        p5 = v2->getSteinerPoint();
                        ec4 = false;
                        ec5 = false;
                    }
                    else {
                        p5 = getQuadVertex(vs, he3->getPrev()->getVertex());
                        ec4 = he3->getNext()->getEdge()->isConstrained();
                        ec5 = he3->getPrev()->getEdge()->isConstrained();
                    }
                    
                    ec1 = he1->getNext()->getEdge()->isConstrained();
                    ec2 = he9->getPrev()->getEdge()->isConstrained();
                    ec3 = he6->getPrev()->getEdge()->isConstrained();
                    
                    sp = quadrangulateDegPentagon1(p1,p2,p3,p4,p5,ec1,ec2,
                                                   ec3,ec4,ec5,lq);
                }
                else if (e12 != 0) {
                    // polygon degenerated to a septagon
                    p1 = getQuadVertex(vs, he9->getNext()->getVertex());
                    
                    if (v5->hasSteinerPoint()) {
                        p2 = v5->getSteinerPoint();
                        ec1 = false;
                        ec2 = false;
                    }
                    else {
                        p2 = getQuadVertex(vs, he9->getPrev()->getVertex());
                        ec1 = he9->getNext()->getEdge()->isConstrained();
                        ec2 = he9->getPrev()->getEdge()->isConstrained();
                    }
                    
                    p3 = getQuadVertex(vs, he6->getPrev()->getVertex());
                    p4 = getQuadVertex(vs, he3->getNext()->getVertex());
                    p5 = getQuadVertex(vs, he1->getPrev()->getVertex());
                    
                    ec3 = he6->getPrev()->getEdge()->isConstrained();
                    ec4 = he3->getNext()->getEdge()->isConstrained();
                    ec5 = he1->getPrev()->getEdge()->isConstrained();
                    
                    p6 = getQuadVertex(vs, he7->getNext()->getVertex());
                    if (v4->hasSteinerPoint()) {
                        p7 = v4->getSteinerPoint();
                        ec6 = false;
                        ec7 = false;
                    }
                    else {
                        p7 = getQuadVertex(vs, he7->getPrev()->getVertex());
                        ec6 = he7->getNext()->getEdge()->isConstrained();
                        ec7 = he7->getPrev()->getEdge()->isConstrained();
                    }
                    
                    sp = quadrangulateDegPentagon1(p1,p2,p3,p6,p7,ec1,ec2,
                                                   false,ec6,ec7,lq);
                    sp = quadrangulateDegPentagon1(p6,sp,p3,p4,p5,false,
                                                   false,ec3,ec4,ec5,lq);
                }
                else if (e14 != 0) {
                    // polygon degenerated to a septagon
                    p1 = getQuadVertex(vs, he1->getNext()->getVertex());
                    p2 = getQuadVertex(vs, he7->getPrev()->getVertex());
                    
                    ec1 = he1->getNext()->getEdge()->isConstrained();
                    ec2 = he7->getPrev()->getEdge()->isConstrained();
                    
                    p3 = getQuadVertex(vs, he9->getNext()->getVertex());
                    
                    if (v5->hasSteinerPoint()) {
                        p4 = v5->getSteinerPoint();
                        ec3 = false;
                        ec4 = false;
                    }
                    else {
                        p4 = getQuadVertex(vs, he9->getPrev()->getVertex());
                        ec3 = he7->getNext()->getEdge()->isConstrained();
                        ec4 = he9->getPrev()->getEdge()->isConstrained();
                    }
                    
                    p5 = getQuadVertex(vs, he6->getPrev()->getVertex());
                    
                    ec5 = he6->getPrev()->getEdge()->isConstrained();
                    
                    p6 = getQuadVertex(vs, he3->getNext()->getVertex());
                    if (v2->hasSteinerPoint()) {
                        p7 = v2->getSteinerPoint();
                        ec6 = false;
                        ec7 = false;
                    }
                    else {
                        p7 = getQuadVertex(vs, he3->getPrev()->getVertex());
                        ec6 = he3->getNext()->getEdge()->isConstrained();
                        ec7 = he3->getPrev()->getEdge()->isConstrained();
                    }
                    
                    tCQMIndVertex2 spaux1, spaux2;
                    
                    sp = getQuadVertex(vs, he5->getVertex());
                    spaux1 = quadrangulateDegPentagon1(p1,p2,sp,p6,p7,ec1,false,
                                                       false,ec6,ec7,lq);
                    spaux2 = quadrangulateDegPentagon1(p3,p4,p5,sp,p2,ec3,ec4,
                                                       false,false,ec2,lq);
                    sp = quadrangulateDegPentagon1(sp,spaux2,p5,p6,spaux1,false,
                                                   false,ec5,false,false,lq);
                }
                else {  // e45 != 0
                    // polygon degenerated to a septagon
                    p1 = getQuadVertex(vs, he1->getNext()->getVertex());
                    
                    if (v1->hasSteinerPoint()) {
                        p2 = v1->getSteinerPoint();
                        ec1 = false;
                        ec2 = false;
                    }
                    else {
                        p2 = getQuadVertex(vs, he1->getPrev()->getVertex());
                        ec1 = he1->getNext()->getEdge()->isConstrained();
                        ec2 = he1->getPrev()->getEdge()->isConstrained();
                    }
                    
                    p3 = getQuadVertex(vs, he7->getNext()->getVertex());
                    p4 = getQuadVertex(vs, he9->getPrev()->getVertex());
                    p5 = getQuadVertex(vs, he6->getPrev()->getVertex());
                    
                    ec3 = he7->getNext()->getEdge()->isConstrained();
                    ec4 = he9->getPrev()->getEdge()->isConstrained();
                    ec5 = he6->getPrev()->getEdge()->isConstrained();
                    
                    p6 = getQuadVertex(vs, he3->getNext()->getVertex());
                    if (v2->hasSteinerPoint()) {
                        p7 = v2->getSteinerPoint();
                        ec6 = false;
                        ec7 = false;
                    }
                    else {
                        p7 = getQuadVertex(vs, he3->getPrev()->getVertex());
                        ec6 = he3->getNext()->getEdge()->isConstrained();
                        ec7 = he3->getPrev()->getEdge()->isConstrained();
                    }
                    
                    sp = quadrangulateDegPentagon1(p1,p2,p3,p6,p7,ec1,ec2,
                                                   false,ec6,ec7,lq);
                    sp = quadrangulateDegPentagon1(p3,p4,p5,p6,sp,ec3,ec4,
                                                   ec5,false,false,lq);
                }
            }
            else {
                // "v4" is the right child of "v6"
                if ((e12 != 0) && (e15 != 0) && (e45 != 0)) {
                    // polygon degenerated to a triangle
                    p1 = getQuadVertex(vs, he7->getPrev()->getVertex());
                    p2 = getQuadVertex(vs, he6->getPrev()->getVertex());
                    p3 = getQuadVertex(vs, he3->getNext()->getVertex());
                    p4 = getQuadVertex(vs, he10->getVertex());
                    
                    ec1 = he7->getPrev()->getEdge()->isConstrained();
                    ec3 = he3->getNext()->getEdge()->isConstrained();
                    
                    sp = quadrangulateDegTriangle(p1,p2,p3,p4,ec1,
                                                  ec3,lq);
                }
                else if ((e12 != 0) && (e15 != 0)) {
                    // polygon degenerated to a pentagon
                    p1 = getQuadVertex(vs, he7->getNext()->getVertex());
                    
                    if (v4->hasSteinerPoint()) {
                        p2 = v4->getSteinerPoint();
                        ec1 = false;
                        ec2 = false;
                    }
                    else {
                        p2 = getQuadVertex(vs, he7->getPrev()->getVertex());
                        ec1 = he7->getNext()->getEdge()->isConstrained();
                        ec2 = he7->getPrev()->getEdge()->isConstrained();
                    }
                    
                    p3 = getQuadVertex(vs, he6->getPrev()->getVertex());
                    p4 = getQuadVertex(vs, he3->getNext()->getVertex());
                    p5 = getQuadVertex(vs, he9->getPrev()->getVertex());
                    
                    ec3 = he6->getPrev()->getEdge()->isConstrained();
                    ec4 = he3->getNext()->getEdge()->isConstrained();
                    ec5 = he9->getPrev()->getEdge()->isConstrained();
                    
                    sp = quadrangulateDegPentagon1(p1,p2,p3,p4,p5,ec1,ec2,
                                                   ec3,ec4,ec5,lq);
                }
                else if ((e12 != 0) && (e45 != 0)) {
                    // polygon degenerated to a pentagon
                    p1 = getQuadVertex(vs, he9->getNext()->getVertex());
                    p2 = getQuadVertex(vs, he7->getPrev()->getVertex());
                    p3 = getQuadVertex(vs, he6->getPrev()->getVertex());
                    p4 = getQuadVertex(vs, he3->getNext()->getVertex());
                    p5 = getQuadVertex(vs, he1->getPrev()->getVertex());
                    
                    ec1 = he9->getNext()->getEdge()->isConstrained();
                    ec2 = he7->getPrev()->getEdge()->isConstrained();
                    ec3 = he6->getPrev()->getEdge()->isConstrained();
                    ec4 = he3->getNext()->getEdge()->isConstrained();
                    ec5 = he1->getPrev()->getEdge()->isConstrained();
                    
                    sp = quadrangulateDegPentagon1(p1,p2,p3,p4,p5,ec1,ec2,
                                                   ec3,ec4,ec5,lq);
                }
                else if ((e15 != 0) && (e45 != 0)) {
                    // polygon degenerated to a pentagon
                    p1 = getQuadVertex(vs, he1->getNext()->getVertex());
                    p2 = getQuadVertex(vs, he7->getPrev()->getVertex());
                    p3 = getQuadVertex(vs, he6->getPrev()->getVertex());
                    p4 = getQuadVertex(vs, he3->getNext()->getVertex());
                    
                    if (v2->hasSteinerPoint()) {
                        p5 = v2->getSteinerPoint();
                        ec4 = false;
                        ec5 = false;
                    }
                    else {
                        p5 = getQuadVertex(vs, he3->getPrev()->getVertex());
                        ec4 = he3->getNext()->getEdge()->isConstrained();
                        ec5 = he3->getPrev()->getEdge()->isConstrained();
                    }
                    
                    ec1 = he1->getNext()->getEdge()->isConstrained();
                    ec2 = he7->getPrev()->getEdge()->isConstrained();
                    ec3 = he6->getPrev()->getEdge()->isConstrained();
                    
                    sp = quadrangulateDegPentagon1(p1,p2,p3,p4,p5,ec1,ec2,
                                                   ec3,ec4,ec5,lq);
                }
                else if (e12 != 0) {
                    // polygon degenerated to a septagon
                    p1 = getQuadVertex(vs, he7->getNext()->getVertex());
                    
                    if (v4->hasSteinerPoint()) {
                        p2 = v4->getSteinerPoint();
                        ec1 = false;
                        ec2 = false;
                    }
                    else {
                        p2 = getQuadVertex(vs, he7->getPrev()->getVertex());
                        ec1 = he7->getNext()->getEdge()->isConstrained();
                        ec2 = he7->getPrev()->getEdge()->isConstrained();
                    }
                    
                    p3 = getQuadVertex(vs, he6->getPrev()->getVertex());
                    p4 = getQuadVertex(vs, he3->getNext()->getVertex());
                    p5 = getQuadVertex(vs, he1->getPrev()->getVertex());
                    
                    ec3 = he6->getPrev()->getEdge()->isConstrained();
                    ec4 = he3->getNext()->getEdge()->isConstrained();
                    ec5 = he1->getPrev()->getEdge()->isConstrained();
                    
                    p6 = getQuadVertex(vs, he9->getNext()->getVertex());
                    if (v5->hasSteinerPoint()) {
                        p7 = v5->getSteinerPoint();
                        ec6 = false;
                        ec7 = false;
                    }
                    else {
                        p7 = getQuadVertex(vs, he9->getPrev()->getVertex());
                        ec6 = he9->getNext()->getEdge()->isConstrained();
                        ec7 = he9->getPrev()->getEdge()->isConstrained();
                    }
                    
                    sp = quadrangulateDegPentagon1(p1,p2,p3,p6,p7,ec1,ec2,
                                                   false,ec6,ec7,lq);
                    sp = quadrangulateDegPentagon1(p6,sp,p3,p4,p5,false,
                                                   false,ec3,ec4,ec5,lq);
                }
                else if (e15 != 0) {
                    // polygon degenerated to a septagon
                    p1 = getQuadVertex(vs, he1->getNext()->getVertex());
                    p2 = getQuadVertex(vs, he9->getPrev()->getVertex());
                    
                    ec1 = he1->getNext()->getEdge()->isConstrained();
                    ec2 = he9->getPrev()->getEdge()->isConstrained();
                    
                    p3 = getQuadVertex(vs, he7->getNext()->getVertex());
                    
                    if (v4->hasSteinerPoint()) {
                        p4 = v4->getSteinerPoint();
                        ec3 = false;
                        ec4 = false;
                    }
                    else {
                        p4 = getQuadVertex(vs, he7->getPrev()->getVertex());
                        ec3 = he7->getNext()->getEdge()->isConstrained();
                        ec4 = he7->getPrev()->getEdge()->isConstrained();
                    }
                    
                    p5 = getQuadVertex(vs, he6->getPrev()->getVertex());
                    
                    ec5 = he6->getPrev()->getEdge()->isConstrained();
                    
                    p6 = getQuadVertex(vs, he3->getNext()->getVertex());
                    
                    if (v2->hasSteinerPoint()) {
                        p7 = v2->getSteinerPoint();
                        ec6 = false;
                        ec7 = false;
                    }
                    else {
                        p7 = getQuadVertex(vs, he3->getPrev()->getVertex());
                        ec6 = he6->getNext()->getEdge()->isConstrained();
                        ec7 = he3->getPrev()->getEdge()->isConstrained();
                    }
                    
                    tCQMIndVertex2 spaux1, spaux2;
                    
                    sp = getQuadVertex(vs, he5->getVertex());
                    spaux1 = quadrangulateDegPentagon1(p1,p2,sp,p6,p7,ec1,false,
                                                       false,ec6,ec7,lq);
                    spaux2 = quadrangulateDegPentagon1(p3,p4,p5,sp,p2,ec3,ec4,
                                                       false,false,ec2,lq);
                    sp = quadrangulateDegPentagon1(sp,spaux2,p5,p6,spaux1,false,
                                                   false,ec5,false,false,lq);
                }
                else {  // e45 != 0
                    // polygon degenerated to a septagon
                    p1 = getQuadVertex(vs, he1->getNext()->getVertex());
                    
                    if (v1->hasSteinerPoint()) {
                        p2 = v1->getSteinerPoint();
                        ec1 = false;
                        ec2 = false;
                    }
                    else {
                        p2 = getQuadVertex(vs, he1->getPrev()->getVertex());
                        ec1 = he1->getNext()->getEdge()->isConstrained();
                        ec2 = he1->getPrev()->getEdge()->isConstrained();
                    }
                    
                    p3 = getQuadVertex(vs, he9->getNext()->getVertex());
                    p4 = getQuadVertex(vs, he7->getPrev()->getVertex());
                    p5 = getQuadVertex(vs, he6->getPrev()->getVertex());
                    
                    ec3 = he9->getNext()->getEdge()->isConstrained();
                    ec4 = he7->getPrev()->getEdge()->isConstrained();
                    ec5 = he6->getPrev()->getEdge()->isConstrained();
                    
                    p6 = getQuadVertex(vs, he3->getNext()->getVertex());
                    if (v2->hasSteinerPoint()) {
                        p7 = v2->getSteinerPoint();
                        ec6 = false;
                        ec7 = false;
                    }
                    else {
                        p7 = getQuadVertex(vs, he3->getPrev()->getVertex());
                        ec6 = he3->getNext()->getEdge()->isConstrained();
                        ec7 = he3->getPrev()->getEdge()->isConstrained();
                    }
                    
                    sp = quadrangulateDegPentagon1(p1,p2,p3,p6,p7,ec1,ec2,
                                                   false,ec6,ec7,lq);
                    sp = quadrangulateDegPentagon1(p3,p4,p5,p6,sp,ec3,ec4,
                                                   ec5,false,false,lq);
                }
            }
        }
        
        hesp = he6->getPrev();
    }
    else {
        // "v6" is the left child of "v7"
        if (he2->getNext() == he4) {
            // "v1" is the left child of "v3"
            if (he8->getNext() == he10) {
                // "v4" is the left child of "v6"
                if ((e12 != 0) && (e15 != 0) && (e45 != 0)) {
                    // polygon degenerated to a triangle
                    p1 = getQuadVertex(vs, he3->getPrev()->getVertex());
                    p2 = getQuadVertex(vs, he6->getNext()->getVertex());
                    p3 = getQuadVertex(vs, he7->getNext()->getVertex());
                    p4 = getQuadVertex(vs, he2->getVertex());
                    
                    ec1 = he3->getPrev()->getEdge()->isConstrained();
                    ec3 = he7->getNext()->getEdge()->isConstrained();
                    
                    sp = quadrangulateDegTriangle(p1,p2,p3,p4,ec1,
                                                  ec3,lq);
                }
                else if ((e12 != 0) && (e15 != 0)) {
                    // polygon degenerated to a pentagon
                    p1 = getQuadVertex(vs, he9->getNext()->getVertex());
                    p2 = getQuadVertex(vs, he3->getPrev()->getVertex());
                    p3 = getQuadVertex(vs, he6->getNext()->getVertex());
                    p4 = getQuadVertex(vs, he7->getNext()->getVertex());
                    
                    ec1 = he9->getNext()->getEdge()->isConstrained();
                    ec2 = he3->getPrev()->getEdge()->isConstrained();
                    ec3 = he6->getNext()->getEdge()->isConstrained();
                    
                    if (v4->hasSteinerPoint()) {
                        p5 = v4->getSteinerPoint();
                        ec4 = false;
                        ec5 = false;
                    }
                    else {
                        p5 = getQuadVertex(vs, he7->getPrev()->getVertex());
                        ec4 = he7->getNext()->getEdge()->isConstrained();
                        ec5 = he7->getPrev()->getEdge()->isConstrained();
                    }
                    
                    sp = quadrangulateDegPentagon1(p1,p2,p3,p4,p5,ec1,ec2,
                                                   ec3,ec4,ec5,lq);
                }
                else if ((e12 != 0) && (e45 != 0)) {
                    // polygon degenerated to a pentagon
                    p1 = getQuadVertex(vs, he1->getNext()->getVertex());
                    p2 = getQuadVertex(vs, he3->getPrev()->getVertex());
                    p3 = getQuadVertex(vs, he6->getNext()->getVertex());
                    p4 = getQuadVertex(vs, he7->getNext()->getVertex());
                    p5 = getQuadVertex(vs, he9->getPrev()->getVertex());
                    
                    ec1 = he1->getNext()->getEdge()->isConstrained();
                    ec2 = he3->getPrev()->getEdge()->isConstrained();
                    ec3 = he6->getNext()->getEdge()->isConstrained();
                    ec4 = he7->getNext()->getEdge()->isConstrained();
                    ec5 = he9->getPrev()->getEdge()->isConstrained();
                    
                    sp = quadrangulateDegPentagon1(p1,p2,p3,p4,p5,ec1,ec2,
                                                   ec3,ec4,ec5,lq);
                }
                else if ((e15 != 0) && (e45 != 0)) {
                    // polygon degenerated to a pentagon
                    p1 = getQuadVertex(vs, he3->getNext()->getVertex());
                    
                    if (v2->hasSteinerPoint()) {
                        p2 = v2->getSteinerPoint();
                        ec1 = false;
                        ec2 = false;
                    }
                    else {
                        p2 = getQuadVertex(vs, he3->getPrev()->getVertex());
                        ec1 = he3->getNext()->getEdge()->isConstrained();
                        ec2 = he3->getPrev()->getEdge()->isConstrained();
                    }
                    
                    p3 = getQuadVertex(vs, he6->getNext()->getVertex());
                    p4 = getQuadVertex(vs, he7->getNext()->getVertex());
                    p5 = getQuadVertex(vs, he1->getPrev()->getVertex());
                    
                    ec3 = he6->getNext()->getEdge()->isConstrained();
                    ec4 = he7->getNext()->getEdge()->isConstrained();
                    ec5 = he1->getPrev()->getEdge()->isConstrained();
                    
                    sp = quadrangulateDegPentagon1(p1,p2,p3,p4,p5,ec1,ec2,
                                                   ec3,ec4,ec5,lq);
                }
                else if (e12 != 0) {
                    // polygon degenerated to a septagon
                    p1 = getQuadVertex(vs, he9->getNext()->getVertex());
                    
                    if (v5->hasSteinerPoint()) {
                        p2 = v5->getSteinerPoint();
                        ec1 = false;
                        ec2 = false;
                    }
                    else {
                        p2 = getQuadVertex(vs, he9->getPrev()->getVertex());
                        ec1 = he9->getNext()->getEdge()->isConstrained();
                        ec2 = he9->getPrev()->getEdge()->isConstrained();
                    }
                    
                    p3 = getQuadVertex(vs, he1->getNext()->getVertex());
                    p4 = getQuadVertex(vs, he3->getPrev()->getVertex());
                    p5 = getQuadVertex(vs, he6->getNext()->getVertex());
                    
                    ec3 = he1->getNext()->getEdge()->isConstrained();
                    ec4 = he3->getPrev()->getEdge()->isConstrained();
                    ec5 = he6->getNext()->getEdge()->isConstrained();
                    
                    p6 = getQuadVertex(vs, he7->getNext()->getVertex());
                    if (v4->hasSteinerPoint()) {
                        p7 = v4->getSteinerPoint();
                        ec6 = false;
                        ec7 = false;
                    }
                    else {
                        p7 = getQuadVertex(vs, he7->getPrev()->getVertex());
                        ec6 = he7->getNext()->getEdge()->isConstrained();
                        ec7 = he7->getPrev()->getEdge()->isConstrained();
                    }
                    
                    sp = quadrangulateDegPentagon1(p1,p2,p3,p6,p7,ec1,ec2,
                                                   false,ec6,ec7,lq);
                    sp = quadrangulateDegPentagon1(p3,p4,p5,p6,sp,ec3,ec4,
                                                   ec5,false,false,lq);
                }
                else if (e15 != 0) {
                    // polygon degenerated to a septagon
                    p1 = getQuadVertex(vs, he9->getNext()->getVertex());
                    p2 = getQuadVertex(vs, he1->getPrev()->getVertex());
                    p3 = getQuadVertex(vs, he3->getNext()->getVertex());
                    
                    ec1 = he9->getNext()->getEdge()->isConstrained();
                    ec2 = he1->getPrev()->getEdge()->isConstrained();
                    
                    if (v2->hasSteinerPoint()) {
                        p4 = v2->getSteinerPoint();
                        ec3 = false;
                        ec4 = false;
                    }
                    else {
                        p4 = getQuadVertex(vs, he3->getPrev()->getVertex());
                        ec3 = he3->getNext()->getEdge()->isConstrained();
                        ec4 = he3->getPrev()->getEdge()->isConstrained();
                    }
                    
                    p5 = getQuadVertex(vs, he6->getNext()->getVertex());
                    ec5 = he6->getNext()->getEdge()->isConstrained();
                    
                    p6 = getQuadVertex(vs, he7->getNext()->getVertex());
                    if (v4->hasSteinerPoint()) {
                        p7 = v4->getSteinerPoint();
                        ec6 = false;
                        ec7 = false;
                    }
                    else {
                        p7 = getQuadVertex(vs, he7->getPrev()->getVertex());
                        ec6 = he7->getNext()->getEdge()->isConstrained();
                        ec7 = he7->getPrev()->getEdge()->isConstrained();
                    }
                    
                    tCQMIndVertex2 spaux1, spaux2;
                    
                    sp = getQuadVertex(vs, he6->getVertex());
                    spaux1 = quadrangulateDegPentagon1(p1,p2,sp,p6,p7,ec1,false,
                                                       false,ec6,ec7,lq);
                    spaux2 = quadrangulateDegPentagon1(p3,p4,p5,sp,p2,ec3,ec4,
                                                       false,false,ec2,lq);
                    sp = quadrangulateDegPentagon1(sp,spaux2,p5,p6,spaux1,false,
                                                   false,ec5,false,false,lq);
                }
                else {  // e45 != 0
                    // polygon degenerated to a septagon
                    p1 = getQuadVertex(vs, he3->getNext()->getVertex());
                    
                    if (v2->hasSteinerPoint()) {
                        p2 = v2->getSteinerPoint();
                        ec1 = false;
                        ec2 = false;
                    }
                    else {
                        p2 = getQuadVertex(vs, he3->getPrev()->getVertex());
                        ec1 = he3->getNext()->getEdge()->isConstrained();
                        ec2 = he3->getPrev()->getEdge()->isConstrained();
                    }
                    
                    p3 = getQuadVertex(vs, he6->getNext()->getVertex());
                    p4 = getQuadVertex(vs, he7->getNext()->getVertex());
                    p5 = getQuadVertex(vs, he9->getPrev()->getVertex());
                    
                    ec3 = he6->getNext()->getEdge()->isConstrained();
                    ec4 = he7->getNext()->getEdge()->isConstrained();
                    ec5 = he9->getPrev()->getEdge()->isConstrained();
                    
                    p6 = getQuadVertex(vs, he1->getNext()->getVertex());
                    if (v1->hasSteinerPoint()) {
                        p7 = v1->getSteinerPoint();
                        ec6 = false;
                        ec7 = false;
                    }
                    else {
                        p7 = getQuadVertex(vs, he1->getPrev()->getVertex());
                        ec6 = he1->getNext()->getEdge()->isConstrained();
                        ec7 = he1->getPrev()->getEdge()->isConstrained();
                    }
                    
                    sp = quadrangulateDegPentagon1(p1,p2,p3,p6,p7,ec1,ec2,
                                                   false,ec6,ec7,lq);
                    sp = quadrangulateDegPentagon1(p6,sp,p3,p4,p5,false,
                                                   false,ec3,ec4,ec5,lq);
                }
            }
            else {
                // "v4" is the right child of "v3"
                if ((e12 != 0) && (e14 != 0) && (e45 != 0)) {
                    // polygon degenerated to a triangle
                    p1 = getQuadVertex(vs, he3->getPrev()->getVertex());
                    p2 = getQuadVertex(vs, he6->getNext()->getVertex());
                    p3 = getQuadVertex(vs, he9->getNext()->getVertex());
                    p4 = getQuadVertex(vs, he2->getVertex());
                    
                    ec1 = he3->getPrev()->getEdge()->isConstrained();
                    ec3 = he9->getNext()->getEdge()->isConstrained();
                    
                    sp = quadrangulateDegTriangle(p1,p2,p3,p4,ec1,
                                                  ec3,lq);
                }
                else if ((e12 != 0) && (e14 != 0)) {
                    // polygon degenerated to a pentagon
                    p1 = getQuadVertex(vs, he7->getNext()->getVertex());
                    p2 = getQuadVertex(vs, he3->getPrev()->getVertex());
                    p3 = getQuadVertex(vs, he6->getNext()->getVertex());
                    p4 = getQuadVertex(vs, he9->getNext()->getVertex());
                    
                    ec1 = he7->getNext()->getEdge()->isConstrained();
                    ec2 = he3->getPrev()->getEdge()->isConstrained();
                    ec3 = he6->getNext()->getEdge()->isConstrained();
                    
                    if (v5->hasSteinerPoint()) {
                        p5 = v5->getSteinerPoint();
                        ec4 = false;
                        ec5 = false;
                    }
                    else {
                        p5 = getQuadVertex(vs, he7->getPrev()->getVertex());
                        ec4 = he7->getNext()->getEdge()->isConstrained();
                        ec5 = he7->getPrev()->getEdge()->isConstrained();
                    }
                    
                    sp = quadrangulateDegPentagon1(p1,p2,p3,p4,p5,ec1,ec2,
                                                   ec3,ec4,ec5,lq);
                }
                else if ((e12 != 0) && (e45 != 0)) {
                    // polygon degenerated to a pentagon
                    p1 = getQuadVertex(vs, he1->getNext()->getVertex());
                    p2 = getQuadVertex(vs, he3->getPrev()->getVertex());
                    p3 = getQuadVertex(vs, he6->getNext()->getVertex());
                    p4 = getQuadVertex(vs, he9->getNext()->getVertex());
                    p5 = getQuadVertex(vs, he7->getPrev()->getVertex());
                    
                    ec1 = he1->getNext()->getEdge()->isConstrained();
                    ec2 = he3->getPrev()->getEdge()->isConstrained();
                    ec3 = he6->getNext()->getEdge()->isConstrained();
                    ec4 = he9->getNext()->getEdge()->isConstrained();
                    ec5 = he7->getPrev()->getEdge()->isConstrained();
                    
                    sp = quadrangulateDegPentagon1(p1,p2,p3,p4,p5,ec1,ec2,
                                                   ec3,ec4,ec5,lq);
                }
                else if ((e14 != 0) && (e45 != 0)) {
                    // polygon degenerated to a pentagon
                    p1 = getQuadVertex(vs, he3->getNext()->getVertex());
                    
                    if (v2->hasSteinerPoint()) {
                        p2 = v2->getSteinerPoint();
                        ec1 = false;
                        ec2 = false;
                    }
                    else {
                        p2 = getQuadVertex(vs, he3->getPrev()->getVertex());
                        ec1 = he3->getNext()->getEdge()->isConstrained();
                        ec2 = he3->getPrev()->getEdge()->isConstrained();
                    }
                    
                    p3 = getQuadVertex(vs, he6->getNext()->getVertex());
                    p4 = getQuadVertex(vs, he9->getNext()->getVertex());
                    p5 = getQuadVertex(vs, he1->getPrev()->getVertex());
                    
                    ec3 = he6->getNext()->getEdge()->isConstrained();
                    ec4 = he9->getNext()->getEdge()->isConstrained();
                    ec5 = he1->getPrev()->getEdge()->isConstrained();
                    
                    sp = quadrangulateDegPentagon1(p1,p2,p3,p4,p5,ec1,ec2,
                                                   ec3,ec4,ec5,lq);
                }
                else if (e12 != 0) {
                    // polygon degenerated to a septagon
                    p1 = getQuadVertex(vs, he7->getNext()->getVertex());
                    
                    if (v4->hasSteinerPoint()) {
                        p2 = v4->getSteinerPoint();
                        ec1 = false;
                        ec2 = false;
                    }
                    else {
                        p2 = getQuadVertex(vs, he7->getPrev()->getVertex());
                        ec1 = he7->getNext()->getEdge()->isConstrained();
                        ec2 = he7->getPrev()->getEdge()->isConstrained();
                    }
                    
                    p3 = getQuadVertex(vs, he1->getNext()->getVertex());
                    p4 = getQuadVertex(vs, he3->getPrev()->getVertex());
                    p5 = getQuadVertex(vs, he6->getNext()->getVertex());
                    
                    ec3 = he1->getNext()->getEdge()->isConstrained();
                    ec4 = he3->getPrev()->getEdge()->isConstrained();
                    ec5 = he6->getNext()->getEdge()->isConstrained();
                    
                    p6 = getQuadVertex(vs, he9->getNext()->getVertex());
                    if (v5->hasSteinerPoint()) {
                        p7 = v5->getSteinerPoint();
                        ec6 = false;
                        ec7 = false;
                    }
                    else {
                        p7 = getQuadVertex(vs, he9->getPrev()->getVertex());
                        ec6 = he9->getNext()->getEdge()->isConstrained();
                        ec7 = he9->getPrev()->getEdge()->isConstrained();
                    }
                    
                    sp = quadrangulateDegPentagon1(p1,p2,p3,p6,p7,ec1,ec2,
                                                   false,ec6,ec7,lq);
                    sp = quadrangulateDegPentagon1(p3,p4,p5,p6,sp,ec3,ec4,
                                                   ec5,false,false,lq);
                }
                else if (e14 != 0) {
                    // polygon degenerated to a septagon
                    p1 = getQuadVertex(vs, he7->getNext()->getVertex());
                    p2 = getQuadVertex(vs, he1->getPrev()->getVertex());
                    p3 = getQuadVertex(vs, he3->getNext()->getVertex());
                    
                    ec1 = he7->getNext()->getEdge()->isConstrained();
                    ec2 = he1->getPrev()->getEdge()->isConstrained();
                    
                    if (v2->hasSteinerPoint()) {
                        p4 = v2->getSteinerPoint();
                        ec3 = false;
                        ec4 = false;
                    }
                    else {
                        p4 = getQuadVertex(vs, he3->getPrev()->getVertex());
                        ec3 = he3->getNext()->getEdge()->isConstrained();
                        ec4 = he3->getPrev()->getEdge()->isConstrained();
                    }
                    
                    p5 = getQuadVertex(vs, he6->getNext()->getVertex());
                    ec5 = he6->getNext()->getEdge()->isConstrained();
                    p6 = getQuadVertex(vs, he9->getNext()->getVertex());
                    
                    if (v5->hasSteinerPoint()) {
                        p7 = v5->getSteinerPoint();
                        ec6 = false;
                        ec7 = false;
                    }
                    else {
                        p7 = getQuadVertex(vs, he9->getPrev()->getVertex());
                        ec6 = he9->getNext()->getEdge()->isConstrained();
                        ec7 = he9->getPrev()->getEdge()->isConstrained();
                    }
                    
                    tCQMIndVertex2 spaux1, spaux2;
                    
                    sp = getQuadVertex(vs, he6->getVertex());
                    spaux1 = quadrangulateDegPentagon1(p1,p2,sp,p6,p7,ec1,false,
                                                       false,ec6,ec7,lq);
                    spaux2 = quadrangulateDegPentagon1(p3,p4,p5,sp,p2,ec3,ec4,
                                                       false,false,ec2,lq);
                    sp = quadrangulateDegPentagon1(sp,spaux2,p5,p6,spaux1,false,
                                                   false,ec5,false,false,lq);
                }
                else {  // e45 != 0
                    // polygon degenerated to a septagon
                    p1 = getQuadVertex(vs, he3->getNext()->getVertex());
                    
                    if (v2->hasSteinerPoint()) {
                        p2 = v2->getSteinerPoint();
                        ec1 = false;
                        ec2 = false;
                    }
                    else {
                        p2 = getQuadVertex(vs, he3->getPrev()->getVertex());
                        ec1 = he3->getNext()->getEdge()->isConstrained();
                        ec2 = he3->getPrev()->getEdge()->isConstrained();
                    }
                    
                    p3 = getQuadVertex(vs, he6->getNext()->getVertex());
                    p4 = getQuadVertex(vs, he9->getNext()->getVertex());
                    p5 = getQuadVertex(vs, he7->getPrev()->getVertex());
                    
                    ec3 = he6->getNext()->getEdge()->isConstrained();
                    ec4 = he9->getNext()->getEdge()->isConstrained();
                    ec5 = he7->getPrev()->getEdge()->isConstrained();
                    
                    p6 = getQuadVertex(vs, he1->getNext()->getVertex());
                    if (v1->hasSteinerPoint()) {
                        p7 = v1->getSteinerPoint();
                        ec6 = false;
                        ec7 = false;
                    }
                    else {
                        p7 = getQuadVertex(vs, he1->getPrev()->getVertex());
                        ec6 = he1->getNext()->getEdge()->isConstrained();
                        ec7 = he1->getPrev()->getEdge()->isConstrained();
                    }
                    
                    sp = quadrangulateDegPentagon1(p1,p2,p3,p6,p7,ec1,ec2,
                                                   false,ec6,ec7,lq);
                    sp = quadrangulateDegPentagon1(p6,sp,p3,p4,p5,false,
                                                   false,ec3,ec4,ec5,lq);
                }
            }
        }
        else {
            // "v1" is the right child of "v3"
            if (he8->getNext() == he10) {
                // "v4" is the left child of "v6"
                if ((e12 != 0) && (e25 != 0) && (e45 != 0)) {
                    // polygon degenerated to a triangle
                    p1 = getQuadVertex(vs, he1->getPrev()->getVertex());
                    p2 = getQuadVertex(vs, he6->getNext()->getVertex());
                    p3 = getQuadVertex(vs, he7->getNext()->getVertex());
                    p4 = getQuadVertex(vs, he4->getVertex());
                    
                    ec1 = he1->getPrev()->getEdge()->isConstrained();
                    ec3 = he7->getNext()->getEdge()->isConstrained();
                    
                    sp = quadrangulateDegTriangle(p1,p2,p3,p4,ec1,
                                                  ec3,lq);
                }
                else if ((e12 != 0) && (e25 != 0)) {
                    // polygon degenerated to a pentagon
                    p1 = getQuadVertex(vs, he9->getNext()->getVertex());
                    p2 = getQuadVertex(vs, he1->getPrev()->getVertex());
                    p3 = getQuadVertex(vs, he6->getNext()->getVertex());
                    p4 = getQuadVertex(vs, he7->getNext()->getVertex());
                    
                    ec1 = he9->getNext()->getEdge()->isConstrained();
                    ec2 = he1->getPrev()->getEdge()->isConstrained();
                    ec3 = he6->getNext()->getEdge()->isConstrained();
                    
                    if (v4->hasSteinerPoint()) {
                        p5 = v4->getSteinerPoint();
                        ec4 = false;
                        ec5 = false;
                    }
                    else {
                        p5 = getQuadVertex(vs, he7->getPrev()->getVertex());
                        ec4 = he7->getNext()->getEdge()->isConstrained();
                        ec5 = he7->getPrev()->getEdge()->isConstrained();
                    }
                    
                    sp = quadrangulateDegPentagon1(p1,p2,p3,p4,p5,ec1,ec2,
                                                   ec3,ec4,ec5,lq);
                }
                else if ((e12 != 0) && (e45 != 0)) {
                    // polygon degenerated to a pentagon
                    p1 = getQuadVertex(vs, he3->getNext()->getVertex());
                    p2 = getQuadVertex(vs, he1->getPrev()->getVertex());
                    p3 = getQuadVertex(vs, he6->getNext()->getVertex());
                    p4 = getQuadVertex(vs, he7->getNext()->getVertex());
                    p5 = getQuadVertex(vs, he9->getPrev()->getVertex());
                    
                    ec1 = he3->getNext()->getEdge()->isConstrained();
                    ec2 = he1->getPrev()->getEdge()->isConstrained();
                    ec3 = he6->getNext()->getEdge()->isConstrained();
                    ec4 = he7->getNext()->getEdge()->isConstrained();
                    ec5 = he9->getPrev()->getEdge()->isConstrained();
                    
                    sp = quadrangulateDegPentagon1(p1,p2,p3,p4,p5,ec1,ec2,
                                                   ec3,ec4,ec5,lq);
                }
                else if ((e25 != 0) && (e45 != 0)) {
                    // polygon degenerated to a pentagon
                    p1 = getQuadVertex(vs, he1->getNext()->getVertex());
                    
                    if (v2->hasSteinerPoint()) {
                        p2 = v2->getSteinerPoint();
                        ec1 = false;
                        ec2 = false;
                    }
                    else {
                        p2 = getQuadVertex(vs, he1->getPrev()->getVertex());
                        ec1 = he1->getNext()->getEdge()->isConstrained();
                        ec2 = he1->getPrev()->getEdge()->isConstrained();
                    }
                    
                    p3 = getQuadVertex(vs, he6->getNext()->getVertex());
                    p4 = getQuadVertex(vs, he7->getNext()->getVertex());
                    p5 = getQuadVertex(vs, he3->getPrev()->getVertex());
                    
                    ec3 = he6->getNext()->getEdge()->isConstrained();
                    ec4 = he7->getNext()->getEdge()->isConstrained();
                    ec5 = he3->getPrev()->getEdge()->isConstrained();
                    
                    sp = quadrangulateDegPentagon1(p1,p2,p3,p4,p5,ec1,ec2,
                                                   ec3,ec4,ec5,lq);
                }
                else if (e12 != 0) {
                    // polygon degenerated to a septagon
                    p1 = getQuadVertex(vs, he9->getNext()->getVertex());
                    
                    if (v5->hasSteinerPoint()) {
                        p2 = v5->getSteinerPoint();
                        ec1 = false;
                        ec2 = false;
                    }
                    else {
                        p2 = getQuadVertex(vs, he9->getPrev()->getVertex());
                        ec1 = he9->getNext()->getEdge()->isConstrained();
                        ec2 = he9->getPrev()->getEdge()->isConstrained();
                    }
                    
                    p3 = getQuadVertex(vs, he3->getNext()->getVertex());
                    p4 = getQuadVertex(vs, he1->getPrev()->getVertex());
                    p5 = getQuadVertex(vs, he6->getNext()->getVertex());
                    
                    ec3 = he3->getNext()->getEdge()->isConstrained();
                    ec4 = he1->getPrev()->getEdge()->isConstrained();
                    ec5 = he6->getNext()->getEdge()->isConstrained();
                    
                    p6 = getQuadVertex(vs, he7->getNext()->getVertex());
                    if (v4->hasSteinerPoint()) {
                        p7 = v4->getSteinerPoint();
                        ec6 = false;
                        ec7 = false;
                    }
                    else {
                        p7 = getQuadVertex(vs, he7->getPrev()->getVertex());
                        ec6 = he7->getNext()->getEdge()->isConstrained();
                        ec7 = he7->getPrev()->getEdge()->isConstrained();
                    }
                    
                    sp = quadrangulateDegPentagon1(p1,p2,p3,p6,p7,ec1,ec2,
                                                   false,ec6,ec7,lq);
                    sp = quadrangulateDegPentagon1(p3,p4,p5,p6,sp,ec3,ec4,
                                                   ec5,false,false,lq);
                }
                else if (e25 != 0) {
                    // polygon degenerated to a septagon
                    p1 = getQuadVertex(vs, he9->getNext()->getVertex());
                    p2 = getQuadVertex(vs, he3->getPrev()->getVertex());
                    p3 = getQuadVertex(vs, he1->getNext()->getVertex());
                    
                    ec1 = he9->getNext()->getEdge()->isConstrained();
                    ec2 = he3->getPrev()->getEdge()->isConstrained();
                    
                    if (v1->hasSteinerPoint()) {
                        p4 = v1->getSteinerPoint();
                        ec3 = false;
                        ec4 = false;
                    }
                    else {
                        p4 = getQuadVertex(vs, he1->getPrev()->getVertex());
                        ec3 = he1->getNext()->getEdge()->isConstrained();
                        ec4 = he1->getPrev()->getEdge()->isConstrained();
                    }
                    
                    p5 = getQuadVertex(vs, he6->getNext()->getVertex());
                    ec5 = he6->getNext()->getEdge()->isConstrained();
                    p6 = getQuadVertex(vs, he7->getNext()->getVertex());
                    if (v4->hasSteinerPoint()) {
                        p7 = v4->getSteinerPoint();
                        ec6 = false;
                        ec7 = false;
                    }
                    else {
                        p7 = getQuadVertex(vs, he7->getPrev()->getVertex());
                        ec6 = he7->getNext()->getEdge()->isConstrained();
                        ec7 = he7->getPrev()->getEdge()->isConstrained();
                    }
                    
                    tCQMIndVertex2 spaux1, spaux2;
                    
                    sp = getQuadVertex(vs, he6->getVertex());
                    spaux1 = quadrangulateDegPentagon1(p1,p2,sp,p6,p7,ec1,false,
                                                       false,ec6,ec7,lq);
                    spaux2 = quadrangulateDegPentagon1(p3,p4,p5,sp,p2,ec3,ec4,
                                                       false,false,ec2,lq);
                    sp = quadrangulateDegPentagon1(sp,spaux2,p5,p6,spaux1,false,
                                                   false,ec5,false,false,lq);
                }
                else {  // e45 != 0
                    // polygon degenerated to a septagon
                    p1 = getQuadVertex(vs, he1->getNext()->getVertex());
                    
                    if (v1->hasSteinerPoint()) {
                        p2 = v1->getSteinerPoint();
                        ec1 = false;
                        ec2 = false;
                    }
                    else {
                        p2 = getQuadVertex(vs, he1->getPrev()->getVertex());
                        ec1 = he1->getNext()->getEdge()->isConstrained();
                        ec2 = he1->getPrev()->getEdge()->isConstrained();
                    }
                    
                    p3 = getQuadVertex(vs, he6->getNext()->getVertex());
                    p4 = getQuadVertex(vs, he7->getNext()->getVertex());
                    p5 = getQuadVertex(vs, he9->getPrev()->getVertex());
                    
                    ec3 = he6->getNext()->getEdge()->isConstrained();
                    ec4 = he7->getNext()->getEdge()->isConstrained();
                    ec5 = he9->getPrev()->getEdge()->isConstrained();
                    
                    p6 = getQuadVertex(vs, he3->getNext()->getVertex());
                    if (v2->hasSteinerPoint()) {
                        p7 = v2->getSteinerPoint();
                        ec6 = false;
                        ec7 = false;
                    }
                    else {
                        p7 = getQuadVertex(vs, he3->getPrev()->getVertex());
                        ec6 = he3->getNext()->getEdge()->isConstrained();
                        ec7 = he3->getPrev()->getEdge()->isConstrained();
                    }
                    
                    sp = quadrangulateDegPentagon1(p1,p2,p3,p6,p7,ec1,ec2,
                                                   false,ec6,ec7,lq);
                    sp = quadrangulateDegPentagon1(p6,sp,p3,p4,p5,false,
                                                   false,ec3,ec4,ec5,lq);
                }
            }
            else {
                // "v4" is the right child of "v3"
                if ((e12 != 0) && (e24 != 0) && (e45 != 0)) {
                    // polygon degenerated to a triangle
                    p1 = getQuadVertex(vs, he1->getPrev()->getVertex());
                    p2 = getQuadVertex(vs, he6->getNext()->getVertex());
                    p3 = getQuadVertex(vs, he9->getNext()->getVertex());
                    p4 = getQuadVertex(vs, he6->getVertex());
                    
                    ec1 = he1->getPrev()->getEdge()->isConstrained();
                    ec3 = he9->getNext()->getEdge()->isConstrained();
                    
                    sp = quadrangulateDegTriangle(p1,p2,p3,p4,ec1,
                                                  ec3,lq);
                }
                else if ((e12 != 0) && (e24 != 0)) {
                    // polygon degenerated to a pentagon
                    p1 = getQuadVertex(vs, he7->getNext()->getVertex());
                    p2 = getQuadVertex(vs, he1->getPrev()->getVertex());
                    p3 = getQuadVertex(vs, he6->getNext()->getVertex());
                    p4 = getQuadVertex(vs, he9->getNext()->getVertex());
                    
                    ec1 = he7->getNext()->getEdge()->isConstrained();
                    ec2 = he1->getPrev()->getEdge()->isConstrained();
                    ec3 = he6->getNext()->getEdge()->isConstrained();
                    
                    if (v5->hasSteinerPoint()) {
                        p5 = v5->getSteinerPoint();
                        ec4 = false;
                        ec5 = false;
                    }
                    else {
                        p5 = getQuadVertex(vs, he9->getPrev()->getVertex());
                        ec4 = he9->getNext()->getEdge()->isConstrained();
                        ec5 = he9->getPrev()->getEdge()->isConstrained();
                    }
                    
                    sp = quadrangulateDegPentagon1(p1,p2,p3,p4,p5,ec1,ec2,
                                                   ec3,ec4,ec5,lq);
                }
                else if ((e12 != 0) && (e45 != 0)) {
                    // polygon degenerated to a pentagon
                    p1 = getQuadVertex(vs, he3->getNext()->getVertex());
                    p2 = getQuadVertex(vs, he1->getPrev()->getVertex());
                    p3 = getQuadVertex(vs, he6->getNext()->getVertex());
                    p4 = getQuadVertex(vs, he9->getNext()->getVertex());
                    p5 = getQuadVertex(vs, he7->getPrev()->getVertex());
                    
                    ec1 = he3->getNext()->getEdge()->isConstrained();
                    ec2 = he1->getPrev()->getEdge()->isConstrained();
                    ec3 = he6->getNext()->getEdge()->isConstrained();
                    ec4 = he9->getNext()->getEdge()->isConstrained();
                    ec5 = he7->getPrev()->getEdge()->isConstrained();
                    
                    sp = quadrangulateDegPentagon1(p1,p2,p3,p4,p5,ec1,ec2,
                                                   ec3,ec4,ec5,lq);
                }
                else if ((e24 != 0) && (e45 != 0)) {
                    // polygon degenerated to a pentagon
                    p1 = getQuadVertex(vs, he1->getNext()->getVertex());
                    
                    if (v1->hasSteinerPoint()) {
                        p2 = v1->getSteinerPoint();
                        ec1 = false;
                        ec2 = false;
                    }
                    else {
                        p2 = getQuadVertex(vs, he1->getPrev()->getVertex());
                        ec1 = he1->getNext()->getEdge()->isConstrained();
                        ec2 = he1->getPrev()->getEdge()->isConstrained();
                    }
                    
                    p3 = getQuadVertex(vs, he6->getNext()->getVertex());
                    p4 = getQuadVertex(vs, he9->getNext()->getVertex());
                    p5 = getQuadVertex(vs, he3->getPrev()->getVertex());
                    
                    ec3 = he6->getNext()->getEdge()->isConstrained();
                    ec4 = he9->getNext()->getEdge()->isConstrained();
                    ec5 = he3->getPrev()->getEdge()->isConstrained();
                    
                    sp = quadrangulateDegPentagon1(p1,p2,p3,p4,p5,ec1,ec2,
                                                   ec3,ec4,ec5,lq);
                }
                else if (e12 != 0) {
                    // polygon degenerated to a septagon
                    p1 = getQuadVertex(vs, he7->getNext()->getVertex());
                    
                    if (v4->hasSteinerPoint()) {
                        p2 = v4->getSteinerPoint();
                        ec1 = false;
                        ec2 = false;
                    }
                    else {
                        p2 = getQuadVertex(vs, he7->getPrev()->getVertex());
                        ec1 = he7->getNext()->getEdge()->isConstrained();
                        ec2 = he7->getPrev()->getEdge()->isConstrained();
                    }
                    
                    p3 = getQuadVertex(vs, he3->getNext()->getVertex());
                    p4 = getQuadVertex(vs, he1->getPrev()->getVertex());
                    p5 = getQuadVertex(vs, he6->getNext()->getVertex());
                    
                    ec3 = he3->getNext()->getEdge()->isConstrained();
                    ec4 = he1->getPrev()->getEdge()->isConstrained();
                    ec5 = he6->getNext()->getEdge()->isConstrained();
                    
                    p6 = getQuadVertex(vs, he9->getNext()->getVertex());
                    if (v5->hasSteinerPoint()) {
                        p7 = v5->getSteinerPoint();
                        ec6 = false;
                        ec7 = false;
                    }
                    else {
                        p7 = getQuadVertex(vs, he9->getPrev()->getVertex());
                        ec6 = he9->getNext()->getEdge()->isConstrained();
                        ec7 = he9->getPrev()->getEdge()->isConstrained();
                    }
                    
                    sp = quadrangulateDegPentagon1(p1,p2,p3,p6,p7,ec1,ec2,
                                                   false,ec6,ec7,lq);
                    sp = quadrangulateDegPentagon1(p3,p4,p5,p6,sp,ec3,ec4,
                                                   ec5,false,false,lq);
                }
                else if (e24 != 0) {
                    // polygon degenerated to a septagon
                    p1 = getQuadVertex(vs, he7->getNext()->getVertex());
                    p2 = getQuadVertex(vs, he3->getPrev()->getVertex());
                    p3 = getQuadVertex(vs, he1->getNext()->getVertex());
                    
                    ec1 = he7->getNext()->getEdge()->isConstrained();
                    ec2 = he3->getPrev()->getEdge()->isConstrained();
                    
                    if (v1->hasSteinerPoint()) {
                        p4 = v1->getSteinerPoint();
                        ec3 = false;
                        ec4 = false;
                    }
                    else {
                        p4 = getQuadVertex(vs, he1->getPrev()->getVertex());
                        ec3 = he1->getNext()->getEdge()->isConstrained();
                        ec4 = he1->getPrev()->getEdge()->isConstrained();
                    }
                    
                    p5 = getQuadVertex(vs, he6->getNext()->getVertex());
                    ec5 = he6->getNext()->getEdge()->isConstrained();
                    
                    p6 = getQuadVertex(vs, he9->getNext()->getVertex());
                    if (v5->hasSteinerPoint()) {
                        p7 = v5->getSteinerPoint();
                        ec6 = false;
                        ec7 = false;
                    }
                    else {
                        p7 = getQuadVertex(vs, he9->getPrev()->getVertex());
                        ec6 = he9->getNext()->getEdge()->isConstrained();
                        ec7 = he9->getPrev()->getEdge()->isConstrained();
                    }
                    
                    tCQMIndVertex2 spaux1, spaux2;
                    
                    sp = getQuadVertex(vs, he6->getVertex());
                    spaux1 = quadrangulateDegPentagon1(p1,p2,sp,p6,p7,ec1,false,
                                                       false,ec6,ec7,lq);
                    spaux2 = quadrangulateDegPentagon1(p3,p4,p5,sp,p2,ec3,ec4,
                                                       false,false,ec2,lq);
                    sp = quadrangulateDegPentagon1(sp,spaux2,p5,p6,spaux1,false,
                                                   false,ec5,false,false,lq);
                }
                else {  // e45 != 0
                    // polygon degenerated to a septagon
                    p1 = getQuadVertex(vs, he1->getNext()->getVertex());
                    
                    if (v1->hasSteinerPoint()) {
                        p2 = v1->getSteinerPoint();
                        ec1 = false;
                        ec2 = false;
                    }
                    else {
                        p2 = getQuadVertex(vs, he1->getPrev()->getVertex());
                        ec1 = he1->getNext()->getEdge()->isConstrained();
                        ec2 = he1->getPrev()->getEdge()->isConstrained();
                    }
                    
                    p3 = getQuadVertex(vs, he6->getNext()->getVertex());
                    p4 = getQuadVertex(vs, he9->getNext()->getVertex());
                    p5 = getQuadVertex(vs, he7->getPrev()->getVertex());
                    
                    ec3 = he6->getNext()->getEdge()->isConstrained();
                    ec4 = he9->getNext()->getEdge()->isConstrained();
                    ec5 = he7->getPrev()->getEdge()->isConstrained();
                    
                    p6 = getQuadVertex(vs, he3->getNext()->getVertex());
                    if (v2->hasSteinerPoint()) {
                        p7 = v2->getSteinerPoint();
                        ec6 = false;
                        ec7 = false;
                    }
                    else {
                        p7 = getQuadVertex(vs, he3->getPrev()->getVertex());
                        ec6 = he3->getNext()->getEdge()->isConstrained();
                        ec7 = he3->getPrev()->getEdge()->isConstrained();
                    }
                    
                    sp = quadrangulateDegPentagon1(p1,p2,p3,p6,p7,ec1,ec2,
                                                   false,ec6,ec7,lq);
                    sp = quadrangulateDegPentagon1(p6,sp,p3,p4,p5,false,
                                                   false,ec3,ec4,ec5,lq);
                }
            }
        }
        
        hesp = he6->getNext();
    }
    
    v7->setSteinerPoint(sp,hesp);
}

// -------------------------------------------------------------------
// Method makeQuad
// -------------------------------------------------------------------
void
tCQMCompQuad::makeQuad(
                       tCQMSpanningTreeVertex* v1,
                       tCQMSpanningTreeVertex* v2,
                       tCQMSpanningTreeVertex* v3,
                       tCQMSpanningTreeVertex* v4,
                       tCQMSpanningTreeVertex* v5,
                       QuadVertexSet& vs,
                       std::list<tCQMQuadrilateral2*>* lq
                       )
{
    // v1 and v2 are children of v3, and v3 and v4 are children of v5.
    
    // get the mesh faces corresponding to the dual graph vertex in "v1"
    // "v2","v3", "v4" and "v5".
    tCQMFace2* f1 = v1->getVertex()->getFace();
    tCQMFace2* f2 = v2->getVertex()->getFace();
    tCQMFace2* f3 = v3->getVertex()->getFace();
    tCQMFace2* f4 = v4->getVertex()->getFace();
    tCQMFace2* f5 = v5->getVertex()->getFace();
    
    if ((f1->getCommonEdge(f2) != 0) || (f1->getCommonEdge(f4) != 0) ||
        (f2->getCommonEdge(f4) != 0)) {
        makeQuadFromDegeneratePoly(v1,v2,v3,v4,v5,vs,lq);
    }
    else {
        
        // get the common edge of "f1" and "f3"
        tCQMEdge2* e1 = f3->getCommonEdge(f1);
        
        // get the common edge of "f2" and "f3"
        tCQMEdge2* e2 = f3->getCommonEdge(f2);
        
        // get the common edge of "f4" and "f5"
        tCQMEdge2* e3 = f5->getCommonEdge(f4);
        
        // get the half-edge of e1 in f1
        tCQMHalfEdge2* he1 = e1->getHalfEdge();
        if (he1->getFace() != f1) {
            he1 = e1->getMate(he1);
        }
        
        // get the half-edge of e2 in f2
        tCQMHalfEdge2* he2 = e2->getHalfEdge();
        if (he2->getFace() != f2) {
            he2 = e2->getMate(he2);
        }
        
        // get the half-edge of "e3" in "f4" and "f5"
        tCQMHalfEdge2* he3 = e3->getHalfEdge();
        tCQMHalfEdge2* he4 = e3->getMate(he3);
        if (he3->getFace() != f4) {
            he3 = he4;
            he4 = e3->getMate(he3);
        }
        
        // find half-edge  of "f5" that does  not belong to  a common edge
        // with "f4" neither "f3".
        tCQMEdge2* e4 = f5->getCommonEdge(f3);
        tCQMHalfEdge2* hesp = e3->getHalfEdge();
        if (hesp->getFace() != f5) {
            hesp = e4->getMate(hesp);
        }
        if (hesp->getNext() == he4) {
            hesp = hesp->getPrev();
        }
        
        // get  vertices of  the  face corresponding  to  the dual  vertex
        // associated with "v".
        tCQMIndVertex2 p1 = getQuadVertex(vs, he1->getVertex());
        tCQMIndVertex2 p2 = getQuadVertex(vs, he1->getNext()->getVertex());
        
        // compute edge status
        bool ec2 = he1->getNext()->getEdge()->isConstrained();
        
        tCQMIndVertex2 p3;
        bool ec3;
        if (v1->hasSteinerPoint()) {
            p3 = v1->getSteinerPoint();
            ec3 = false;
        }
        else {
            p3 = getQuadVertex(vs, he1->getPrev()->getVertex());
            ec3 = he1->getPrev()->getEdge()->isConstrained();
        }
        
        // get  vertices of  the  face corresponding  to  the dual  vertex
        // associated with "vsib".
        tCQMIndVertex2 p4 = getQuadVertex(vs, he2->getVertex());
        tCQMIndVertex2 p5 = getQuadVertex(vs, he2->getNext()->getVertex());
        
        // compute edge status
        bool ec5 = he2->getNext()->getEdge()->isConstrained();
        
        tCQMIndVertex2 p6;
        bool ec6;
        if (v2->hasSteinerPoint()) {
            p6 = v2->getSteinerPoint();
            ec6 = false;
        }
        else {
            p6 = getQuadVertex(vs, he2->getPrev()->getVertex());
            ec6 = he2->getPrev()->getEdge()->isConstrained();
        }
        
        // get  vertices of  the  face corresponding  to  the dual  vertex
        // associated with "vparsib".
        tCQMIndVertex2 p7 = getQuadVertex(vs, he3->getVertex());
        tCQMIndVertex2 p8 = getQuadVertex(vs, he3->getNext()->getVertex());
        
        // compute edge status
        bool ec8 = he3->getNext()->getEdge()->isConstrained();
        
        tCQMIndVertex2 p9;
        bool ec9;
        if (v4->hasSteinerPoint()) {
            p9 = v4->getSteinerPoint();
            ec9 = false;
        }
        else {
            p9 = getQuadVertex(vs, he3->getPrev()->getVertex());
            ec9 = he3->getPrev()->getEdge()->isConstrained();
        }
        
        // check if  there is a point  in "f5" that is  visible from "f1",
        // "f2", "f3" and  "f4". If so, we need only  this point to convex
        // quadrangulate the  polygon defined by the union  of those faces
        // and  we are left  with one  triangle.  Otherwise,  we subdivide
        // this  polygon into  one  triangle, one  quadrilateral, and  one
        // hexagon,  then we  convex quadrangulate  the hexagon.   In this
        // case, we use at most 4 Steiner points.
        
        // get  vertices of  the  face corresponding  to  the dual  vertex
        // associated with "vparvpar".
        tCQMIndVertex2 q1 = getQuadVertex(vs, he4->getVertex());
        tCQMIndVertex2 q2 = getQuadVertex(vs, he4->getNext()->getVertex());
        tCQMIndVertex2 q3 = getQuadVertex(vs, he4->getPrev()->getVertex());
        
        // let us compute a visible point in face "f5", if any.
        std::list<Point2> li;
        std::list<Point2> lo;
        
        li.push_back(q1.getPoint()); li.push_back(q2.getPoint());
        li.push_back(q3.getPoint());
        
        try {
            // compute intersection of face "f5" and the wedge of "p9"
            halfSpaceIntersection(li,p9.getPoint(),p7.getPoint(),lo,
                                  CGAL::ON_NEGATIVE_SIDE);
            li.clear();
            halfSpaceIntersection(lo,p9.getPoint(),p8.getPoint(),li,
                                  CGAL::ON_POSITIVE_SIDE);
            lo.clear();
        }
        catch (const tExceptionObject& xpt) {
            treatException(xpt);
            exit(0);
        }
        
        // copy resulting polygon to an auxiliary list
        std::list<Point2>::iterator ilist;
        std::list<Point2> liaux;
        for (ilist = li.begin(); ilist != li.end(); ++ilist) {
            liaux.push_back(*ilist);
        }
        
        // compute intersection  of the resulting polygon  above and the
        // wedge of "p3".
        try {
            halfSpaceIntersection(li,p3.getPoint(),p1.getPoint(),lo,
                                  CGAL::ON_NEGATIVE_SIDE);
            li.clear();
            halfSpaceIntersection(lo,p3.getPoint(),p2.getPoint(),li,
                                  CGAL::ON_POSITIVE_SIDE);
            lo.clear();
        }
        catch (const tExceptionObject& xpt) {
            treatException(xpt);
            exit(0);
        }
        
        if (li.size() > 2) {
            // compute intersection of the resulting polygon above and the
            // wedge of "p6".
            try {
                halfSpaceIntersection(li,p6.getPoint(),p4.getPoint(),lo,
                                      CGAL::ON_NEGATIVE_SIDE);
                li.clear();
                halfSpaceIntersection(lo,p6.getPoint(),p5.getPoint(),li,
                                      CGAL::ON_POSITIVE_SIDE);
                lo.clear();
            }
            catch (const tExceptionObject& xpt) {
                treatException(xpt);
                exit(0);
            }
            
            if (li.size() > 2) {
                // compute intersection  of the resulting  polygon above and
                // the wedge of the common vertex of "f1" and "f3".
                if (p1.getPoint() == p5.getPoint()) {
                    try {
                        halfSpaceIntersection(li,p1.getPoint(),p4.getPoint(),lo,
                                              CGAL::ON_NEGATIVE_SIDE);
                        li.clear();
                        halfSpaceIntersection(lo,p1.getPoint(),p2.getPoint(),li,
                                              CGAL::ON_POSITIVE_SIDE);
                        lo.clear();
                    }
                    catch (const tExceptionObject& xpt) {
                        treatException(xpt);
                        exit(0);
                    }
                }
                else {
                    try {
                        halfSpaceIntersection(li,p2.getPoint(),p1.getPoint(),lo,
                                              CGAL::ON_NEGATIVE_SIDE);
                        li.clear();
                        halfSpaceIntersection(lo,p2.getPoint(),p5.getPoint(),li,
                                              CGAL::ON_POSITIVE_SIDE);
                        lo.clear();
                    }
                    catch (const tExceptionObject& xpt) {
                        treatException(xpt);
                        exit(0);
                    }
                }
                
                if (li.size() > 2) {
                    // compute intersection of the resulting polygon above and
                    // the wedge of "p9".
                    try {
                        halfSpaceIntersection(li,p9.getPoint(),p7.getPoint(),lo,
                                              CGAL::ON_NEGATIVE_SIDE);
                        li.clear();
                        halfSpaceIntersection(lo,p9.getPoint(),p8.getPoint(),li,
                                              CGAL::ON_POSITIVE_SIDE);
                        lo.clear();
                    }
                    catch (const tExceptionObject& xpt) {
                        treatException(xpt);
                        exit(0);
                    }
                }
            }
        }
        
        if (li.size() > 2) {
            // compute the  barycenter of  the resulting polygon  from the
            // previous intersections.
            Point2 q;
            double x = 0;
            double y = 0;
            int num = 0;
            while (!li.empty()) {
                q = li.front();
                li.pop_front();
                x += q.x();
                y += q.y();
                ++num;
            }
            
            q = Point2(x / double(num), y / double(num));
            
            // create a Steiner point inside face "f5"
            tCQMIndVertex2 sp = createSteinerPoint(q);
            
            // store Steiner point in "v5"
            v5->setSteinerPoint(sp,hesp);
            
            // create quadrilaterals
            try {
                createQuad(p1,sp,p2,p3,false,false,ec2,ec3,lq);
                createQuad(p6,p4,sp,p5,ec6,false,false,ec5,lq);
                createQuad(p9,p7,sp,p8,ec9,false,false,ec8,lq);
            }
            catch (const tExceptionObject& xpt) {
                treatException(xpt);
                exit(0);
            }
        }
        else {
            // compute  the barycenter  of  the resulting  polygon from  the
            // intersection of "f5" and the wedge of "p9".
            Point2 q;
            double x = 0;
            double y = 0;
            int num = 0;
            while (!liaux.empty()) {
                q = liaux.front();
                liaux.pop_front();
                x += q.x();
                y += q.y();
                ++num;
            }
            
            q = Point2(x / double(num), y / double(num));
            
            // create a Steiner point inside face "f5"
            tCQMIndVertex2 sp = createSteinerPoint(q);
            
            // store Steiner point in "v5"
            v5->setSteinerPoint(sp,hesp);
            
            // create one quadrilateral
            try {
                createQuad(p9,p7,sp,p8,ec9,false,false,ec8,lq);
            }
            catch (const tExceptionObject& xpt) {
                treatException(xpt);
                exit(0);
            }
            
            // create an hexagon
            tCQMPolygon2* p = (tCQMPolygon2 *) new tCQMPolygon2(6);
            
            // convex quadrangulate  the hexagon  p6p4spp2p3p1 if p1  is the
            // same as p5, or the hexagon  p3p1spp5p6p4 if p2 is the same as
            // p4.
            if (p1.getPoint() == p5.getPoint()) {
                // insert the vertices in the hexagon
                p->insert(0,p6);
                p->insert(1,p4);
                p->insert(2,sp);
                p->insert(3,p2);
                p->insert(4,p3);
                p->insert(5,p1);
                
                p->setConstraint(0,ec6);
                p->setConstraint(1,false);
                p->setConstraint(2,false);
                p->setConstraint(3,ec2);
                p->setConstraint(4,ec3);
                p->setConstraint(5,ec5);
                
                // convex quadrangulate the hexagon p6p4spp2p3p1
                makeConvexQuadFromHex(p,vs,lq);
            }
            else {
                // insert the vertices in the hexagon
                p->insert(0,p3);
                p->insert(1,p1);
                p->insert(2,sp);
                p->insert(3,p5);
                p->insert(4,p6);
                p->insert(5,p4);
                
                p->setConstraint(0,ec3);
                p->setConstraint(1,false);
                p->setConstraint(2,false);
                p->setConstraint(3,ec5);
                p->setConstraint(4,ec6);
                p->setConstraint(5,ec2);
                
                // convex quadrangulate the hexagon p6p4spp2p3p1
                makeConvexQuadFromHex(p,vs,lq);
            }
            
            delete p;
        }
    }
}

// -------------------------------------------------------------------
// Method makeQuad
// -------------------------------------------------------------------
void
tCQMCompQuad::makeQuad(
                       tCQMSpanningTreeVertex* v1,
                       tCQMSpanningTreeVertex* v2,
                       tCQMSpanningTreeVertex* v3,
                       tCQMSpanningTreeVertex* v4,
                       tCQMSpanningTreeVertex* v5,
                       tCQMSpanningTreeVertex* v6,
                       tCQMSpanningTreeVertex* v7,
                       QuadVertexSet& vs,
                       std::list<tCQMQuadrilateral2*>* lq
                       )
{
    // recall that v1 and v2 are  children of v3, v4 and v5 are children
    // of v6, and v3 and v6 are children of v7.
    
    // get the mesh faces corresponding to the dual graph vertex in "v1"
    // "v2","v3", "v4", "v5", "v6", and "v7".
    tCQMFace2* f1 = v1->getVertex()->getFace();
    tCQMFace2* f2 = v2->getVertex()->getFace();
    tCQMFace2* f3 = v3->getVertex()->getFace();
    tCQMFace2* f4 = v4->getVertex()->getFace();
    tCQMFace2* f5 = v5->getVertex()->getFace();
    tCQMFace2* f6 = v6->getVertex()->getFace();
    tCQMFace2* f7 = v7->getVertex()->getFace();
    
    // verify if  we have a degenerate  case in which  "f1", "f2", "f3",
    // "f4", "f5",  "f6" and  "f7" define a  triangle, a pentagon,  or a
    // septagon.  If  so,  we  convex  quadrangulate  the  corresponding
    // degenerate polygon using other method than this one.
    if ((f1->getCommonEdge(f2) != 0) || (f1->getCommonEdge(f4) != 0) ||
        (f1->getCommonEdge(f5) != 0) || (f2->getCommonEdge(f4) != 0) ||
        (f2->getCommonEdge(f5) != 0) || (f4->getCommonEdge(f5) != 0)) {
        makeQuadFromDegeneratePoly(v1,v2,v3,v4,v5,v6,v7,vs,lq);
    }
    else {
        // get the common edge of "f1" and "f3"
        tCQMEdge2* e1 = f3->getCommonEdge(f1);
        
        // get the common edge of "f2" and "f3"
        tCQMEdge2* e2 = f3->getCommonEdge(f2);
        
        // get the common edge of "f4" and "f6"
        tCQMEdge2* e3 = f6->getCommonEdge(f4);
        
        // get the common edge of "f5" and "f6"
        tCQMEdge2* e4 = f6->getCommonEdge(f5);
        
        // get the half-edge of e1 in f1
        tCQMHalfEdge2* he1 = e1->getHalfEdge();
        if (he1->getFace() != f1) {
            he1 = e1->getMate(he1);
        }
        
        // get the half-edge of e2 in f2
        tCQMHalfEdge2* he2 = e2->getHalfEdge();
        if (he2->getFace() != f2) {
            he2 = e2->getMate(he2);
        }
        
        // get the half-edge of e3 in f4
        tCQMHalfEdge2* he3 = e3->getHalfEdge();
        if (he3->getFace() != f4) {
            he3 = e3->getMate(he3);
        }
        
        // get the half-edge of e4 in f5
        tCQMHalfEdge2* he4 = e4->getHalfEdge();
        if (he4->getFace() != f5) {
            he4 = e4->getMate(he4);
        }
        
        // get the half-edge of "f7" that does not belong to a common edge
        // with "f3" nor "f6".
        tCQMEdge2* e5 = f7->getCommonEdge(f3);
        tCQMEdge2* e6 = f7->getCommonEdge(f6);
        
        tCQMHalfEdge2* hesp = e5->getHalfEdge();
        if (hesp->getFace() != f7) {
            hesp = e5->getMate(hesp);
        }
        
        tCQMHalfEdge2* hespaux = e6->getHalfEdge();
        if (hespaux->getFace() != f7) {
            hespaux = e6->getMate(hespaux);
        }
        
        if (hesp->getNext() == hespaux) {
            hesp = hesp->getPrev();
        }
        
        // get vertices of face "f1"
        tCQMIndVertex2 p1 = getQuadVertex(vs, he1->getVertex());
        tCQMIndVertex2 p2 = getQuadVertex(vs, he1->getNext()->getVertex());
        
        // edge status
        bool ep2 = he1->getNext()->getEdge()->isConstrained();
        
        tCQMIndVertex2 p3;
        bool ep3;
        if (v1->hasSteinerPoint()) {
            p3 = v1->getSteinerPoint();
            ep3 = false;
        }
        else {
            p3 = getQuadVertex(vs, he1->getPrev()->getVertex());
            ep3 = he1->getPrev()->getEdge()->isConstrained();
        }
        
        // get vertices of face "f2"
        tCQMIndVertex2 p4 = getQuadVertex(vs, he2->getVertex());
        tCQMIndVertex2 p5 = getQuadVertex(vs, he2->getNext()->getVertex());
        
        // edge status
        bool ep5 = he2->getNext()->getEdge()->isConstrained();
        
        tCQMIndVertex2 p6;
        bool ep6;
        if (v2->hasSteinerPoint()) {
            p6 = v2->getSteinerPoint();
            ep6 = false;
        }
        else {
            p6 = getQuadVertex(vs, he2->getPrev()->getVertex());
            ep6 = he2->getPrev()->getEdge()->isConstrained();
        }
        
        // get vertices of face "f3"
        tCQMIndVertex2 q1, q2, q3;
        if (p1.getPoint() == p5.getPoint()) {
            q1 = p1;
            q2 = p4;
            q3 = p2;
        }
        else {
            q1 = p2;
            q2 = p1;
            q3 = p5;
        }
        
        // get vertices of face "f4"
        tCQMIndVertex2 s1 = getQuadVertex(vs, he3->getVertex());
        tCQMIndVertex2 s2 = getQuadVertex(vs, he3->getNext()->getVertex());
        
        // edge status
        bool es2 = he3->getNext()->getEdge()->isConstrained();
        
        tCQMIndVertex2 s3;
        bool es3;
        if (v4->hasSteinerPoint()) {
            s3 = v4->getSteinerPoint();
            es3 = false;
        }
        else {
            s3 = getQuadVertex(vs, he3->getPrev()->getVertex());
            es3 = he3->getPrev()->getEdge()->isConstrained();
        }
        
        // get vertices of face "f5"
        tCQMIndVertex2 s4 = getQuadVertex(vs, he4->getVertex());
        tCQMIndVertex2 s5 = getQuadVertex(vs, he4->getNext()->getVertex());
        
        // edge status
        bool es5 = he4->getNext()->getEdge()->isConstrained();
        
        tCQMIndVertex2 s6;
        bool es6;
        if (v5->hasSteinerPoint()) {
            s6 = v5->getSteinerPoint();
            es6 = false;
        }
        else {
            s6 = getQuadVertex(vs, he4->getPrev()->getVertex());
            es6 = he4->getPrev()->getEdge()->isConstrained();
        }
        
        // get vertices of face "f6"
        tCQMIndVertex2 q4, q5, q6;
        if (s1.getPoint() == s5.getPoint()) {
            q4 = s1;
            q5 = s4;
            q6 = s2;
        }
        else {
            q4 = s2;
            q5 = s1;
            q6 = s5;
        }
        
        // get vertices of face "f7"
        tCQMIndVertex2 r1, r2, r3;
        if (q2.getPoint() == q6.getPoint()) {
            r1 = q2;
            r2 = q5;
            r3 = q3;
        }
        else {
            r1 = q3;
            r2 = q2;
            r3 = q6;
        }
        
        // check if  there is a point  in "f7" that is  visible from "f1",
        // "f2", "f3", "f4", "f5" and "f6". If so, we need only this point
        // to  convex quadrangulate the  polygon defined  by the  union of
        // those faces  and we are  left with one triangle.  Otherwise, we
        // subdivide this polygon into one triangle and two hexagons, then
        // we convex quadrangulate both hexagons.  In this case, we use at
        // most 7 Steiner points.
        
        // let us compute a visible point in face "f7", if any.
        std::list<Point2> li;
        std::list<Point2> lo;
        
        li.push_back(r1.getPoint()); li.push_back(r2.getPoint());
        li.push_back(r3.getPoint());
        
        try {
            // compute intersection of face "f7" and the wedge of "q1"
            halfSpaceIntersection(li,q1.getPoint(),q2.getPoint(),lo,
                                  CGAL::ON_NEGATIVE_SIDE);
            li.clear();
            halfSpaceIntersection(lo,q1.getPoint(),q3.getPoint(),li,
                                  CGAL::ON_POSITIVE_SIDE);
            lo.clear();
            
            // compute  intersection  of  the  polygon  resulting  from  the
            // previous operation and the wedge of "q6".
            halfSpaceIntersection(li,q4.getPoint(),q5.getPoint(),lo,
                                  CGAL::ON_NEGATIVE_SIDE);
            li.clear();
            halfSpaceIntersection(lo,q4.getPoint(),q6.getPoint(),li,
                                  CGAL::ON_POSITIVE_SIDE);
            lo.clear();
        }
        catch (const tExceptionObject& xpt) {
            treatException(xpt);
            exit(0);
        }
        
        // copy resulting polygon to an auxiliary list
        std::list<Point2>::iterator ilist;
        std::list<Point2> liaux;
        for (ilist = li.begin(); ilist != li.end(); ++ilist) {
            liaux.push_back(*ilist);
        }
        
        // compute  intersection of  the resulting  polygon above  and the
        // wedge of "p3".
        try {
            halfSpaceIntersection(li,p3.getPoint(),p1.getPoint(),lo,
                                  CGAL::ON_NEGATIVE_SIDE);
            li.clear();
            halfSpaceIntersection(lo,p3.getPoint(),p2.getPoint(),li,
                                  CGAL::ON_POSITIVE_SIDE);
            lo.clear();
        }
        catch (const tExceptionObject& xpt) {
            treatException(xpt);
            exit(0);
        }
        
        if (li.size() > 2) {
            // compute intersection  of the resulting polygon  above and the
            // wedge of "p6".
            try {
                halfSpaceIntersection(li,p6.getPoint(),p4.getPoint(),lo,
                                      CGAL::ON_NEGATIVE_SIDE);
                li.clear();
                halfSpaceIntersection(lo,p6.getPoint(),p5.getPoint(),li,
                                      CGAL::ON_POSITIVE_SIDE);
                lo.clear();
            }
            catch (const tExceptionObject& xpt) {
                treatException(xpt);
                exit(0);
            }
            
            if (li.size() > 2) {
                // compute intersection of the resulting polygon above and the
                // wedge of "s3".
                try {
                    halfSpaceIntersection(li,s3.getPoint(),s1.getPoint(),lo,
                                          CGAL::ON_NEGATIVE_SIDE);
                    li.clear();
                    halfSpaceIntersection(lo,s3.getPoint(),s2.getPoint(),li,
                                          CGAL::ON_POSITIVE_SIDE);
                    lo.clear();
                }
                catch (const tExceptionObject& xpt) {
                    treatException(xpt);
                    exit(0);
                }
                
                if (li.size() > 2) {
                    // compute intersection  of the resulting  polygon above and
                    // the wedge of "s6".
                    try {
                        halfSpaceIntersection(li,s6.getPoint(),s4.getPoint(),lo,
                                              CGAL::ON_NEGATIVE_SIDE);
                        li.clear();
                        halfSpaceIntersection(lo,s6.getPoint(),s5.getPoint(),li,
                                              CGAL::ON_POSITIVE_SIDE);
                        lo.clear();
                    }
                    catch (const tExceptionObject& xpt) {
                        treatException(xpt);
                        exit(0);
                    }
                }
            }
        }
        
        if (li.size() > 2) {
            // compute  the barycenter  of  the resulting  polygon from  the
            // previous intersections.
            Point2 q;
            double x = 0;
            double y = 0;
            int num = 0;
            while (!li.empty()) {
                q = li.front();
                li.pop_front();
                x += q.x();
                y += q.y();
                ++num;
            }
            
            q = Point2(x / double(num), y / double(num));
            
            // create a Steiner point inside face "f7"
            tCQMIndVertex2 sp = createSteinerPoint(q);
            
            // store Steiner point in "v7"
            v7->setSteinerPoint(sp,hesp);
            
            // create quadrilaterals
            try {
                createQuad(p1,sp,p2,p3,false,false,ep2,ep3,lq);
                createQuad(p6,p4,sp,p5,ep6,false,false,ep5,lq);
                createQuad(s1,sp,s2,s3,false,false,es2,es3,lq);
                createQuad(s6,s4,sp,s5,es6,false,false,es5,lq);
            }
            catch (const tExceptionObject& xpt) {
                treatException(xpt);
                exit(0);
            }
        }
        else {
            // compute  the barycenter  of  the polygon  resulting from  the
            // intersection of "f7" and the wedges of "q1" and "q6".
            Point2 q;
            double x = 0;
            double y = 0;
            int num = 0;
            while (!liaux.empty()) {
                q = liaux.front();
                liaux.pop_front();
                x += q.x();
                y += q.y();
                ++num;
            }
            
            q = Point2(x / double(num), y / double(num));
            
            // create a Steiner point inside face "f7"
            tCQMIndVertex2 sp = createSteinerPoint(q);
            
            // store Steiner point in "v7"
            v7->setSteinerPoint(sp,hesp);
            
            // create an hexagon
            tCQMPolygon2* h1 = (tCQMPolygon2 *) new tCQMPolygon2(6);
            tCQMPolygon2* h2 = (tCQMPolygon2 *) new tCQMPolygon2(6);
            
            // store Steiner point in "v7"
            v7->setSteinerPoint(sp,hesp);
            
            // convex quadrangulate the hexagon  defined by "f1", "f2", "f3"
            // and the  triangle consisting of  the common edge of  "f3" and
            // "f7" and the Steiner point "sp".
            
            // find out the  vertices of the hexagon defined  by "f1", "f2",
            // "f3" and the  triangle consisting of the common  edge of "f3"
            // and "f7" and the Steiner point "sp".
            if (p1.getPoint() == p5.getPoint()) {
                // insert the vertices in the hexagon
                h1->insert(0,p6);
                h1->insert(1,p4);
                h1->insert(2,sp);
                h1->insert(3,p2);
                h1->insert(4,p3);
                h1->insert(5,p1);
                
                h1->setConstraint(0,ep6);
                h1->setConstraint(1,false);
                h1->setConstraint(2,false);
                h1->setConstraint(3,ep2);
                h1->setConstraint(4,ep3);
                h1->setConstraint(5,ep5);
                
                // convex quadrangulate the hexagon p6p4spp2p3p1
                makeConvexQuadFromHex(h1,vs,lq);
            }
            else {
                // insert the vertices in the hexagon
                h1->insert(0,p3);
                h1->insert(1,p1);
                h1->insert(2,sp);
                h1->insert(3,p5);
                h1->insert(4,p6);
                h1->insert(5,p4);
                
                h1->setConstraint(0,ep3);
                h1->setConstraint(1,false);
                h1->setConstraint(2,false);
                h1->setConstraint(3,ep5);
                h1->setConstraint(4,ep6);
                h1->setConstraint(5,ep2);
                
                // convex quadrangulate the hexagon p3p1spp5p6p4
                makeConvexQuadFromHex(h1,vs,lq);
            }
            
            // convex quadrangulate the hexagon  defined by "f4", "f5", "f6"
            // and the  triangle consisting of  the common edge of  "f6" and
            // "f7" and the Steiner point "sp".
            
            // find out the  vertices of the hexagon defined  by "f4", "f5",
            // "f6" and the  triangle consisting of the common  edge of "f6"
            // and "f7" and the Steiner point "sp".
            if (s1.getPoint() == s5.getPoint()) {
                // insert the vertices in the hexagon
                h2->insert(0,s6);
                h2->insert(1,s4);
                h2->insert(2,sp);
                h2->insert(3,s2);
                h2->insert(4,s3);
                h2->insert(5,s1);
                
                h2->setConstraint(0,es6);
                h2->setConstraint(1,false);
                h2->setConstraint(2,false);
                h2->setConstraint(3,es2);
                h2->setConstraint(4,es3);
                h2->setConstraint(5,es5);
                
                // convex quadrangulate the hexagon s6s4sps2s3s1
                makeConvexQuadFromHex(h2,vs,lq);
            }
            else {
                // insert the vertices in the hexagon
                h2->insert(0,s3);
                h2->insert(1,s1);
                h2->insert(2,sp);
                h2->insert(3,s5);
                h2->insert(4,s6);
                h2->insert(5,s4);
                
                h2->setConstraint(0,es3);
                h2->setConstraint(1,false);
                h2->setConstraint(2,false);
                h2->setConstraint(3,es5);
                h2->setConstraint(4,es6);
                h2->setConstraint(5,es2);
                
                // convex quadrangulate the hexagon s6s4sps3s1sp
                makeConvexQuadFromHex(h2,vs,lq);
            }
            
            delete h1;
            delete h2;
        }
    }
}

// -------------------------------------------------------------------
// Method makeConvexQuadFrom2Tri()
// -------------------------------------------------------------------
void
tCQMCompQuad::makeConvexQuadFrom2Tri(
                                     tCQMSpanningTreeVertex* v1,
                                     tCQMSpanningTreeVertex* v2,
                                     QuadVertexSet& vs,
                                     std::list<tCQMQuadrilateral2*>* lq
                                     )
throw (tExceptionObject)
{
    // both "v1" and "v2" are triangles, and "v2" is the parent of "v1"
    
    if (v1->isQuadrilateral() || v1->isPentagon() ||
        v2->isQuadrilateral()) {
        std::stringstream ss (std::stringstream::in | std::stringstream::out);
        ss << "makeConvexQuadFrom2Tri(): two triangles are expected";
        throw tExceptionObject(__FILE__,__LINE__,ss.str());
    }
    
    // get faces corresponding to "v1" and "v2"
    tCQMFace2* f1 = v1->getVertex()->getFace();
    tCQMFace2* f2 = v2->getVertex()->getFace();
    
    // get the common edge of "f1" and "f2"
    tCQMEdge2* e = f1->getCommonEdge(f2);
    
    // get the half-edge of "e" in "f1"
    tCQMHalfEdge2* he1 = e->getHalfEdge();
    tCQMHalfEdge2* he2 = e->getMate(he1);
    if (he1->getFace() != f1) {
        he2 = he1;
        he1 = e->getMate(he1);
    }
    
    // define variables to store edge status
    bool ec1, ec2, ec3, ec4;
    tCQMIndVertex2 p1 = getQuadVertex(vs, he1->getVertex());
    tCQMIndVertex2 p2 = getQuadVertex(vs, he2->getPrev()->getVertex());
    tCQMIndVertex2 p3 = getQuadVertex(vs, he2->getVertex());
    tCQMIndVertex2 p4;
    if (v1->hasSteinerPoint()) {
        p4 = v1->getSteinerPoint();
        ec3 = false;
        ec4 = false;
    }
    else {
        p4 = getQuadVertex(vs, he1->getPrev()->getVertex());
        ec3 = he1->getNext()->getEdge()->isConstrained();
        ec4 = he1->getPrev()->getEdge()->isConstrained();
    }
    
    ec1 = he2->getNext()->getEdge()->isConstrained();
    ec2 = he2->getPrev()->getEdge()->isConstrained();
    
    // computes the first Steiner point,  which is the barycenter of the
    // intersection of the triangle p1p2p3 and the wedge of "p4".
    std::list<Point2> li;
    std::list<Point2> lo;
    
    li.push_back(p1.getPoint()); li.push_back(p2.getPoint());
    li.push_back(p3.getPoint());
    
    try {
        halfSpaceIntersection(li,p4.getPoint(),p1.getPoint(),lo,
                              CGAL::ON_NEGATIVE_SIDE);
        li.clear();
        halfSpaceIntersection(lo,p4.getPoint(),p3.getPoint(),li,
                              CGAL::ON_POSITIVE_SIDE);
        lo.clear();
    }
    catch (const tExceptionObject& xpt) {
        treatException(xpt);
        exit(0);
    }
    
    // compute the barycenter of the resulting polygon
    Point2 q1;
    double x = 0;
    double y = 0;
    int num = 0;
    while (!li.empty()) {
        q1 = li.front();
        li.pop_front();
        x += q1.x();
        y += q1.y();
        ++num;
    }
    
    q1 = Point2(x / double(num), y / double(num));
    
    // computes the second and third Steiner points
    Point2 q2, q3;
    if (ec2) {
        Line2 l(p1.getPoint(),q1);
        Segment2 s(p2.getPoint(),p3.getPoint());
        
        Point2 pt;
        if (!myAssign(pt,s,l)) {
            // should never happen
            pt = p3.getPoint();
        }
        
        x = (p2.getPoint().x() + pt.x()) / 2;
        y = (p2.getPoint().y() + pt.y()) / 2;
        
        q2 = Point2(x,y);
        
        x = (p3.getPoint().x() + x) / 2;
        y = (p3.getPoint().y() + y) / 2;
        
        if (!he2->getPrev()->getEdge()->isOnBoundary()) {
            tCQMHalfEdge2* heaux = he2->getPrev()->getEdge()->getMate(he2->getPrev());
            double xf = heaux->getPrev()->getVertex()->getPoint().x();
            double yf = heaux->getPrev()->getVertex()->getPoint().y();
            x = x + (xf - x) * 0.05;
            y = y + (yf - y) * 0.05;
        }
        else {
            double xf = p3.getPoint().y() - y;
            double yf = -(p3.getPoint().x() - x);
            double ll = sqrt((xf * xf) + (yf * yf));
            xf /= ll;
            yf /= ll;
            
            x = x + 0.01 * xf;
            y = y + 0.01 * yf;
        }
        
        q3 = Point2(x,y);
        
        // create Steiner point
        tCQMIndVertex2 sp1 = createSteinerPoint(q1);
        tCQMIndVertex2 sp2 = createSteinerPoint(q2);
        tCQMIndVertex2 sp3 = createSteinerPoint(q3);
        
        // create two quadrilaterals
        try {
            createQuad( p1,sp1, p3, p4,false,false,ec3,ec4,lq);
            createQuad( p2,sp2,sp1, p1,ec2,false,false,ec1,lq);
            createQuad(sp2,sp3, p3,sp1,ec2,ec2,false,false,lq);
        }
        catch (const tExceptionObject& xpt) {
            treatException(xpt);
            exit(0);
        }
        
        // split edge containing point and subdivide face that shares this
        // edge with face "f2", if any.
        tCQMVertex2* newv2 = m_mesh->spFmkFE(he2->getPrev(),sp2.getPoint());
        tCQMVertex2* newv3 = m_mesh->spFmkFE(he2->getPrev(),sp3.getPoint());
        
        // insert   new   vertex  in   the   set   of   vertices  of   the
        // quadrangulation.
        vs.insert(make_pair(newv2,sp2.getVertexId()));
        vs.insert(make_pair(newv3,sp3.getVertexId()));
    }
    else {
        Line2 l(p3.getPoint(),q1);
        Segment2 s(p1.getPoint(),p2.getPoint());
        
        Point2 pt;
        if (!myAssign(pt,s,l)) {
            // should never happen
            pt = p1.getPoint();
        }
        
        x = (p2.getPoint().x() + pt.x()) / 2;
        y = (p2.getPoint().y() + pt.y()) / 2;
        
        q2 = Point2(x,y);
        
        x = (p1.getPoint().x() + x) / 2;
        y = (p1.getPoint().y() + y) / 2;
        
        if (!he2->getNext()->getEdge()->isOnBoundary()) {
            tCQMHalfEdge2* heaux = he2->getNext()->getEdge()->getMate(he2->getNext());
            double xf = heaux->getPrev()->getVertex()->getPoint().x();
            double yf = heaux->getPrev()->getVertex()->getPoint().y();
            x = x + (xf - x) * 0.05;
            y = y + (yf - y) * 0.05;
        }
        else {
            double xf = y - p1.getPoint().y();
            double yf = -(x - p1.getPoint().x());
            double ll = sqrt((xf * xf) + (yf * yf));
            xf /= ll;
            yf /= ll;
            
            x = x + 0.01 * xf;
            y = y + 0.01 * yf;
        }
        
        q3 = Point2(x,y);
        
        // create Steiner point
        tCQMIndVertex2 sp1 = createSteinerPoint(q1);
        tCQMIndVertex2 sp2 = createSteinerPoint(q2);
        tCQMIndVertex2 sp3 = createSteinerPoint(q3);
        
        // create two quadrilaterals
        try {
            createQuad( p1,sp1, p3, p4,false,false,ec3,ec4,lq);
            createQuad( p2, p3,sp1,sp2,ec2,false,false,ec1,lq);
            createQuad(sp2,sp1, p1,sp3,false,false,ec1,ec1,lq);
        }
        catch (const tExceptionObject& xpt) {
            treatException(xpt);
            exit(0);
        }
        
        // split edge containing point and subdivide face that shares this
        // edge with face "f2", if any.
        tCQMVertex2* newv2 = m_mesh->spFmkFE(he2->getNext(),sp2.getPoint());
        tCQMVertex2* newv3 = m_mesh->spFmkFE(he2->getNext(),sp3.getPoint());
        
        // insert  new  vertex   in  the  set  of  the   vertices  of  the
        // quadrangulation.
        vs.insert(make_pair(newv2,sp2.getVertexId()));
        vs.insert(make_pair(newv3,sp3.getVertexId()));
    }
}

// -------------------------------------------------------------------
// Method makeConvexQuadFromQuad()
// -------------------------------------------------------------------
void
tCQMCompQuad::makeConvexQuadFromQuad(
                                     tCQMPolygon2* p,
                                     std::list<tCQMQuadrilateral2*>* lq
                                     )
{
    // get polygon points and the status of its edges
    tCQMIndVertex2 p1 = p->getIndVertex(0);
    tCQMIndVertex2 p2 = p->getIndVertex(1);
    tCQMIndVertex2 p3 = p->getIndVertex(2);
    tCQMIndVertex2 p4 = p->getIndVertex(3);
    
    bool ec1 = p->isConstrained(0);
    bool ec2 = p->isConstrained(1);
    bool ec3 = p->isConstrained(2);
    bool ec4 = p->isConstrained(3);
    
    // any quadrilateral can have at most one reflex vertex
    if (p->isReflex(0)) {
        double x = (p1.getPoint().x() + p3.getPoint().x()) / 2;
        double y = (p1.getPoint().y() + p3.getPoint().y()) / 2;
        
        // create Steiner point
        tCQMIndVertex2 sp = createSteinerPoint(Point2(x,y));
        
        // create quadrilaterals
        try {
            createQuad(p1,p2,p3,sp,ec1,ec2,false,false,lq);
            createQuad(p1,sp,p3,p4,false,false,ec3,ec4,lq);
        }
        catch (const tExceptionObject& xpt) {
            treatException(xpt);
            exit(0);
        }
    }
    else if (p->isReflex(1)) {
        double x = (p2.getPoint().x() + p4.getPoint().x()) / 2;
        double y = (p2.getPoint().y() + p4.getPoint().y()) / 2;
        
        // create Steiner point
        tCQMIndVertex2 sp = createSteinerPoint(Point2(x,y));
        
        // create quadrilaterals
        try {
            createQuad(p2,p3,p4,sp,ec2,ec3,false,false,lq);
            createQuad(p2,sp,p4,p1,false,false,ec4,ec1,lq);
        }
        catch (const tExceptionObject& xpt) {
            treatException(xpt);
            exit(0);
        }
    }
    else if (p->isReflex(2)) {
        double x = (p3.getPoint().x() + p1.getPoint().x()) / 2;
        double y = (p3.getPoint().y() + p1.getPoint().y()) / 2;
        
        // create Steiner point
        tCQMIndVertex2 sp = createSteinerPoint(Point2(x,y));
        
        // create quadrilaterals
        try {
            createQuad(p3,p4,p1,sp,ec3,ec4,false,false,lq);
            createQuad(p3,sp,p1,p2,false,false,ec1,ec2,lq);
        }
        catch (const tExceptionObject& xpt) {
            treatException(xpt);
            exit(0);
        }
    }
    else if (p->isReflex(3)) {
        double x = (p4.getPoint().x() + p2.getPoint().x()) / 2;
        double y = (p4.getPoint().y() + p2.getPoint().y()) / 2;
        
        // create Steiner point
        tCQMIndVertex2 sp = createSteinerPoint(Point2(x,y));
        
        // create quadrilaterals
        try {
            createQuad(p4,p1,p2,sp,ec4,ec1,false,false,lq);
            createQuad(p4,sp,p2,p3,false,false,ec2,ec3,lq);
        }
        catch (const tExceptionObject& xpt) {
            treatException(xpt);
            exit(0);
        }
    }
    else {
        try {
            createQuad(p1,p2,p3,p4,ec1,ec2,ec3,ec4,lq);
        }
        catch (const tExceptionObject& xpt) {
            treatException(xpt);
            exit(0);
        }
    }
}

// -------------------------------------------------------------------
// Method makeConvexQuadFromPent()
// -------------------------------------------------------------------
bool
tCQMCompQuad::makeConvexQuadFromPent(
                                     tCQMSpanningTreeVertex* v,
                                     QuadVertexSet& vs,
                                     std::list<tCQMQuadrilateral2*>* lq
                                     )
{
    // verify if "v" is a quadrilateral or a pentagon
    if (v->isQuadrilateral()) {
        // "v" is a quadrilateral.
        // get quadrilateral vertices
        tCQMHalfEdge2* he = v->getQuadHalfEdge();
        
        tCQMIndVertex2 p1 = getQuadVertex(vs,he->getVertex());
        tCQMIndVertex2 p2 = v->getQuadSteinerPoint();
        tCQMIndVertex2 p3 = getQuadVertex(vs,he->getNext()->getVertex());
        tCQMIndVertex2 p4 = getQuadVertex(vs,he->getPrev()->getVertex());
        
        bool ec3 = he->getNext()->getEdge()->isConstrained();
        bool ec4 = he->getPrev()->getEdge()->isConstrained();
        
        try {
            createQuad(p1,p2,p3,p4,false,false,ec3,ec4,lq);
        }
        catch (const tExceptionObject& xpt) {
            treatException(xpt);
            exit(0);
        }
    }
    else {
        // "v" is a pentagon.
        
        // get the half-edges corresponding to the split edges
        tCQMHalfEdge2* he1 = v->getQuadHalfEdge();
        tCQMHalfEdge2* he2 = v->getPentHalfEdge();
        
        tCQMIndVertex2 p1 = getQuadVertex(vs,he1->getVertex());
        tCQMIndVertex2 p2 = v->getQuadSteinerPoint();
        tCQMIndVertex2 p3 = getQuadVertex(vs,he1->getNext()->getVertex());
        
        tCQMIndVertex2 p4;
        tCQMIndVertex2 p5;
        if (he2 == he1->getNext()) {
            p4 = v->getPentSteinerPoint();
            p5 = getQuadVertex(vs,he1->getPrev()->getVertex());
        }
        else {
            p4 = getQuadVertex(vs,he2->getVertex());
            p5 = v->getPentSteinerPoint();
        }
        
        // create quadrilateral
        try {
            createQuad(p1,p2,p3,p4,false,false,false,false,lq);
        }
        catch (const tExceptionObject& xpt) {
            treatException(xpt);
            exit(0);
        }
        
        // update vertex "v"
        v->setSteinerPoint(p4,p5.getPoint());
        v->setQuadPentOff();
    }
    
    return false;
}

// -------------------------------------------------------------------
// Method makeConvexQuadFromPent()
// -------------------------------------------------------------------
bool
tCQMCompQuad::makeConvexQuadFromPent(
                                     tCQMSpanningTreeVertex* v1,
                                     tCQMSpanningTreeVertex* v2,
                                     QuadVertexSet& vs,
                                     std::list<tCQMQuadrilateral2*>* lq
                                     )
{
    // compute the  non-convex pentagon  defined by the  triangular face
    // and  the quadrilateral  corresponding to  the faces  of  the dual
    // vertices associated with "v1" and "v2", respectively.
    int id;
    tCQMPolygon2* p = getPentagon(v1,v2,id,vs);
    
    tCQMIndVertex2 p1 = p->getIndVertex(0);
    tCQMIndVertex2 p2 = p->getIndVertex(1);
    tCQMIndVertex2 p3 = p->getIndVertex(2);
    tCQMIndVertex2 p4 = p->getIndVertex(3);
    tCQMIndVertex2 p5 = p->getIndVertex(4);
    
    // get edge status
    bool ec1 = p->isConstrained(0);
    bool ec2 = p->isConstrained(1);
    bool ec3 = p->isConstrained(2);
    bool ec4 = p->isConstrained(3);
    bool ec5 = p->isConstrained(4);
    
    // if  either  "p1"  or "p4"  is  a  reflex  vertex then  we  convex
    // quadrangulate  "p" using  one Steiner  point.  Otherwise,  "p" is
    // already convex and we subdivide  it into two quadrilaterals and a
    // triangle by adding a Steiner point.
    if (p->isReflex(0) || p->isReflex(3)) {
        // pentagon has one reflex vertex:  either "p1" or "p4". we convex
        // quadrangulate the pentagon by using one Steiner point.
        
        Point2 pt;
        if (p->isReflex(0)) {
            // "p1" is the reflex vertex, so the Steiner point is located in
            // the triangle defined  by "p1", "p4", and "pt",  where "pt" is
            // the intersection  point of the  supporting line p5p1  and the
            // segment p2p4.
            Line2 l(p5.getPoint(),p1.getPoint());
            Segment2 s(p2.getPoint(),p4.getPoint());
            myAssign(pt,s,l);
        }
        else {
            // "p4" is the reflex vertex, so the Steiner point is located in
            // the triangle defined  by "p1", "p4", and "pt",  where "pt" is
            // the intersection  point of the  supporting line p5p4  and the
            // segment p1p3.
            Line2 l(p5.getPoint(),p4.getPoint());
            Segment2 s(p1.getPoint(),p3.getPoint());
            myAssign(pt,s,l);
        }
        
        double x = (pt.x() + p1.getPoint().x() + p4.getPoint().x()) / 3;
        double y = (pt.y() + p1.getPoint().y() + p4.getPoint().y()) / 3;
        
        // create Steiner point
        tCQMIndVertex2 sp = createSteinerPoint(Point2(x,y));
        
        // create quadrilaterals
        if (id == 0) {
            // common edge of the faces  correspoding to "v2" and its parent
            // is p1p2.
            try {
                createQuad(p1,sp,p4,p5,false,false,ec4,ec5,lq);
                createQuad(p4,sp,p2,p3,false,false,ec2,ec3,lq);
            }
            catch (const tExceptionObject& xpt) {
                treatException(xpt);
                exit(0);
            }
            
            // "v2" is now a triangle  whose one vertex is replaced with the
            // Steiner point.
            v2->setSteinerPoint(sp,p1.getPoint());
        }
        else {
            // common edge of the faces  correspoding to "v2" and its parent
            // is p3p4.
            try {
                createQuad(p1,sp,p4,p5,false,false,ec4,ec5,lq);
                createQuad(p1,p2,p3,sp,ec1,ec2,false,false,lq);
            }
            catch (const tExceptionObject& xpt) {
                treatException(xpt);
                exit(0);
            }
            
            // "v2" is now a triangle  whose one vertex is replaced with the
            // Steiner point.
            v2->setSteinerPoint(sp,p3.getPoint());
        }
    }
    else {
        // pentagon is convex and  we subdivide it into two quadrilaterals
        // and a triangle by inserting one Steiner point.
        
        // create Steiner point
        tCQMIndVertex2 sp;
        
        if (id == 0) {
            // common edge of the faces  correspoding to "v2" and its parent
            // is p1p2.
            
            // Steiner  point is placed  in the  barycenter of  the triangle
            // defined by "p1", "p2" and "p4".
            Point2 pt((p1.getPoint().x() + p2.getPoint().x() +
                       p4.getPoint().x()) / 3, (p1.getPoint().y() +
                                                p2.getPoint().y() + p4.getPoint().y()) / 3);
            
            // create Steiner point
            sp = createSteinerPoint(pt);
            
            // create quadrilaterals
            try {
                createQuad(p1,sp,p4,p5,false,false,ec4,ec5,lq);
                createQuad(p2,p3,p4,sp,ec2,ec3,false,false,lq);
            }
            catch (const tExceptionObject& xpt) {
                treatException(xpt);
                exit(0);
            }
            
            // "v2" is now a triangle  whose one vertex is replaced with the
            // Steiner point.
            v2->setSteinerPoint(sp,p1.getPoint());
        }
        else {
            // common edge of the faces  correspoding to "v2" and its parent
            // is p3p4.
            
            // Steiner  point is placed  in the  barycenter of  the triangle
            // defined by "p3", "p4" and "p1".
            Point2 pt((p1.getPoint().x() + p3.getPoint().x() +
                       p4.getPoint().x()) / 3, (p1.getPoint().y() +
                                                p3.getPoint().y() + p4.getPoint().y()) / 3);
            
            // create Steiner point
            sp = createSteinerPoint(pt);
            
            // create quadrilaterals
            try {
                createQuad(p1,sp,p4,p5,false,false,ec4,ec5,lq);
                createQuad(p1,p2,p3,sp,ec1,ec2,false,false,lq);
            }
            catch (const tExceptionObject& xpt) {
                treatException(xpt);
                exit(0);
            }
            
            // "v2" is now a triangle  whose one vertex is replaced with the
            // Steiner point.
            v2->setSteinerPoint(sp,p3.getPoint());
        }
    }
    
    // "v2" is no longer a quadrilateral but a triangle
    v2->setQuadPentOff();
    
    // delete pentagon
    delete p;
    
    return false;
}

// -------------------------------------------------------------------
// Method makeConvexQuadFromPent()
// -------------------------------------------------------------------
bool
tCQMCompQuad::makeConvexQuadFromPent(
                                     tCQMSpanningTreeVertex* v1,
                                     tCQMSpanningTreeVertex* v2,
                                     tCQMSpanningTreeVertex* v3,
                                     QuadVertexSet& vs,
                                     std::list<tCQMQuadrilateral2*>* lq
                                     )
{
    // compute pentagon vertices
    tCQMPolygon2* p = getPentagon(v1,v2,v3,vs);
    
    // based on the number and position of the reflex vertices, classify
    // the pentagon and quadrangulate it.
    bool res = false;
    
    PentCaseCount[classifyPentagon(p,v1,v2,v3)]++;
    std::ofstream myfile;
    myfile.open(spt_filename, std::ios_base::app);
    
    v1->dump(myfile, 1);
    v2->dump(myfile, 2);
    v3->dump(myfile, 3);
    myfile.close();
    
    switch (classifyPentagon(p,v1,v2,v3)) {
            case R1E24E34:
            // only vertex "p1" is reflex, diagonal p2p4 is the common edge of
            // "v2" and "v3", and edge p3p4 is the common edge of "v3" and its
            // parent, if any.
            quadrangulatePentR1E24E34(p,v3,vs,lq);
            break;
            case R4E13E12:
            // only vertex "p4" is reflex, diagonal p1p3 is the common edge of
            // "v2" and "v3", and edge p1p2 is the common edge of "v3" and its
            // parent, if any.
            quadrangulatePentR4E13E12(p,v3,vs,lq);
            break;
            case R1E24E23:
            // only vertex "p1" is reflex, diagonal p2p4 is the common edge of
            // "v2" and "v3", and edge p2p3 is the common edge of "v3" and its
            // parent, if any.
            quadrangulatePentR1E24E23(p,v3,vs,lq);
            break;
            case R4E13E23:
            // only vertex "p4" is reflex, diagonal p1p3 is the common edge of
            // "v2" and "v3", and edge p2p3 is the common edge of "v3" and its
            // parent, if any.
            quadrangulatePentR4E13E23(p,v3,vs,lq);
            break;
            case R1R2E24E34:
            // vertices "p1" and "p2" are  reflex, diagonal p2p4 is the common
            // edge of "v2" and "v3", and edge p3p4 is the common edge of "v3"
            // and its parent, if any.
            try {
                quadrangulatePentR1R2E24E34(p,v3,vs,lq);
                res = true;
            }
            catch (const tExceptionObject& xpt) {
                treatException(xpt);
                exit(0);
            }
            break;
            case R3R4E13E12:
            // vertices "p1" and "p2" are  reflex, diagonal p2p4 is the common
            // edge of "v2" and "v3", and edge p1p2 is the common edge of "v3"
            // and its parent, if any.
            try {
                quadrangulatePentR3R4E13E12(p,v3,vs,lq);
                res = true;
            }
            catch (const tExceptionObject& xpt) {
                treatException(xpt);
                exit(0);
            }
            break;
            case R1R2E24E23:
            // vertices "p1" and "p2" are  reflex, diagonal p2p4 is the common
            // edge of "v2" and "v3", and edge p2p3 is the common edge of "v3"
            // and its parent, if any.
            try {
                quadrangulatePentR1R2E24E23(p,v3,vs,lq);
            }
            catch (const tExceptionObject& xpt) {
                treatException(xpt);
                exit(0);
            }
            break;
            case R3R4E13E23:
            // vertices "p3" and "p4" are  reflex, diagonal p1p3 is the common
            // edge of "v2" and "v3", and edge p2p3 is the common edge of "v3"
            // and its parent, if any.
            try {
                quadrangulatePentR3R4E13E23(p,v3,vs,lq);
            }
            catch (const tExceptionObject& xpt) {
                treatException(xpt);
                exit(0);
            }
            break;
            case R1R4E24E34:
            // vertices "p1" and "p4" are  reflex, diagonal p2p4 is the common
            // edge of "v2" and "v3", and edge p3p4 is the common edge of "v3"
            // and its parent, if any.
            quadrangulatePentR1R4E24E34(p,v3,vs,lq);
            res = true;
            break;
            case R1R4E13E12:
            // vertices "p1" and "p3" are  reflex, diagonal p1p3 is the common
            // edge of "v2" and "v3", and edge p1p2 is the common edge of "v3"
            // and its parent, if any.
            quadrangulatePentR1R4E13E12(p,v3,vs,lq);
            res = true;
            break;
            case R1R4E24E23:
            // vertices "p1" and "p4" are  reflex, diagonal p2p4 is the common
            // edge of "v2" and "v3", and edge p2p3 is the common edge of "v3"
            // and its parent, if any.
            quadrangulatePentR1R4E24E23(p,v3,vs,lq);
            break;
            case R1R4E13E23:
            // vertices "p1" and "p4" are  reflex, diagonal p1p3 is the common
            // edge of "v2" and "v3", and edge p2p3 is the common edge of "v3"
            // and its parent, if any.
            quadrangulatePentR1R4E13E23(p,v3,vs,lq);
            break;
            case R1E13E23:
            // only vertex "p1" is reflex, diagonal p1p3 is the common edge of
            // "v2" and "v3", and edge p2p3 is the common edge of "v3" and its
            // parent, if any.
            quadrangulatePentR1E13E23(p,v3,vs,lq);
            break;
            case R4E24E23:
            // only vertex "p4" is reflex, diagonal p2p4 is the common edge of
            // "v2" and "v3", and edge p2p3 is the common edge of "v3" and its
            // parent, if any.
            quadrangulatePentR4E24E23(p,v3,vs,lq);
            break;
            case R1E13E12:
            // only vertex "p1" is reflex, diagonal p1p3 is the common edge of
            // "v2" and "v3", and edge p1p2 is the common edge of "v3" and its
            // parent, if any.
            quadrangulatePentR1E13E12(p,v3,vs,lq);
            res = true;
            break;
            case R4E24E34:
            // only vertex "p4" is reflex, diagonal p2p4 is the common edge of
            // "v2" and "v3", and edge p3p4 is the common edge of "v3" and its
            // parent, if any.
            quadrangulatePentR4E24E34(p,v3,vs,lq);
            res = true;
            break;
    }
    
    delete p;
    
    return res;
}

// -------------------------------------------------------------------
// Method makeConvexQuadFromHex()
// -------------------------------------------------------------------
bool
tCQMCompQuad::makeConvexQuadFromHex(
                                    tCQMSpanningTreeVertex* v1,
                                    tCQMSpanningTreeVertex* v2,
                                    tCQMSpanningTreeVertex* v3,
                                    QuadVertexSet& vs,
                                    std::list<tCQMQuadrilateral2*>* lq
                                    )
{
    // compute  the hexagon defined  by the  faces corresponding  to the
    // dual vertices  associated with "v1", "v2" and  "v3", where either
    // the face corresponding to the  vertex associated with "v2" or the
    // one  corresponding  to  the  vertex  associated with  "v3"  is  a
    // quadrilateral.
    tCQMPolygon2* p = getHexagon(v1,v2,v3,vs);
    
    makeConvexQuadFromHex(p,vs,lq);
    delete p;
    
    return false;
}

// -------------------------------------------------------------------
// Method makeConvexQuadFromHex()
// -------------------------------------------------------------------
bool
tCQMCompQuad::makeConvexQuadFromHex(
                                    tCQMSpanningTreeVertex* v1,
                                    tCQMSpanningTreeVertex* v2,
                                    tCQMSpanningTreeVertex* v3,
                                    tCQMSpanningTreeVertex* v4,
                                    QuadVertexSet& vs,
                                    std::list<tCQMQuadrilateral2*>* lq
                                    )
{
    // compute the hexagon defined by the triangular faces corresponding
    // to the dual vertices associated with "v1", "v2", "v3" and "v4".
    tCQMPolygon2* p = getHexagon(v1,v2,v3,v4,vs);
    
    if (p->getNumVerts() == 6) {
        makeConvexQuadFromHex(p,vs,lq);
    }
    else {
        makeConvexQuadFromQuad(p,lq);
    }
    
    delete p;
    
    return false;
}

// -------------------------------------------------------------------
// Method makeConvexQuadFromHex()
// -------------------------------------------------------------------
bool
tCQMCompQuad::makeConvexQuadFromHex(
                                    tCQMPolygon2* p,
                                    QuadVertexSet& vs,
                                    std::list<tCQMQuadrilateral2*>* lq
                                    )
{
    // classify the  hexagon according to its number  of reflex vertices
    // and the relative position of  its reflex vertices with respect to
    // each other.
    switch (classifyHexagon(p)) {
            case CCCCCC:
            // "p" has  no reflex vertex,  so it can be  quadrangulate without
            // any Steiner point.
            quadrangulateHexCCCCCC(p,lq);
            break;
            case RCCCCC:
            // "p" has one  reflex vertex, so it can  be quadrangulate with at
            // most two Steiner points.
            quadrangulateHexRCCCCC(p,lq);
            break;
            case RCRCCC:
            // "p" has two reflex vertices, so it can be quadrangulate with at
            // most 3 Steiner points.
            quadrangulateHexRCRCCC(p,lq);
            break;
            case RRCCCC:
            // "p"  has two  consecutive reflex  vertices  so that  it can  be
            // convex quadrangulate with 2 Steiner points.
            quadrangulateHexRRCCCC(p,lq);
            break;
            case RCCRCC:
            // "p" has  two reflex vertices between two  convex vertices, both
            // clockwise  and  counterclockwise,  so  that it  can  be  convex
            // quadrangulate with 2 Steiner points.
            quadrangulateHexRCCRCC(p,lq);
            break;
            case RCRCRC:
            // "p" has  three reflex  vertices and they  alternate so  that we
            // convex quadrangulate "p" by using at most 3 Steiner points.
            quadrangulateHexRCRCRC(p,lq);
            break;
            case RRCRCC:
            // "p" has three  reflex vertices and two of  them are consecutive
            // and  preceded by  two consecutive  convex vertices  so  that we
            // convex quadrangulate "p" by using at most 3 Steiner points.
            quadrangulateHexRRCRCC(p,lq);
            break;
            case RRCCRC:
            // "p" has three  reflex vertices and two of  them are consecutive
            // and  suceeded by  two consecutive  convex vertices  so  that we
            // convex quadrangulate "p" by using at most 3 Steiner points.
            quadrangulateHexRRCCRC(p,lq);
            break;
            case RRRCCC:
            // "p" has three consecutive reflex vertices
            quadrangulateHexRRRCCC(p,lq);
            break;
    }
    
    return false;
}

// -------------------------------------------------------------------
// Method makeConvexQuadFromDegSep()
// -------------------------------------------------------------------
void
tCQMCompQuad::makeConvexQuadFromDegSep(
                                       tCQMSpanningTreeVertex* v1,
                                       tCQMSpanningTreeVertex* v2,
                                       tCQMSpanningTreeVertex* v3,
                                       tCQMSpanningTreeVertex* v4,
                                       tCQMSpanningTreeVertex* v5,
                                       QuadVertexSet& vs,
                                       std::list<tCQMQuadrilateral2*>* lq
                                       )
{
    // v1 and v2 are  children of v3, and v3 and v4  are children of v5,
    // and either "f1"  and "f2" share an edge or "f1"  and "f4" or "f2"
    // and "f4" does.
    
    // get the mesh faces corresponding to the dual graph vertex in "v1"
    // "v2","v3", "v4" and "v5".
    tCQMFace2* f1 = v1->getVertex()->getFace();
    tCQMFace2* f2 = v2->getVertex()->getFace();
    tCQMFace2* f3 = v3->getVertex()->getFace();
    tCQMFace2* f4 = v4->getVertex()->getFace();
    tCQMFace2* f5 = v5->getVertex()->getFace();
    
    // get the common edge of "f1" and "f2"
    tCQMEdge2* e1 = f1->getCommonEdge(f2);
    
    // get the common edge of "f2" and "f3"
    tCQMEdge2* e2 = f2->getCommonEdge(f3);
    
    // get the common edge of "f5" and "f4"
    tCQMEdge2* e3 = f5->getCommonEdge(f4);
    
    // get the common edge of "f4" and "f3"
    tCQMEdge2* e4 = f4->getCommonEdge(f3);
    
    // get the half-edges of "e1" in "f1" and "f2"
    tCQMHalfEdge2* he1 = e1->getHalfEdge();
    tCQMHalfEdge2* he2 = e1->getMate(he1);
    if (he1->getFace() != f1) {
        he2 = he1;
        he1 = e1->getMate(he1);
    }
    
    // get the half-edges of "e2" in "f2" and "f3"
    tCQMHalfEdge2* he3 = e2->getHalfEdge();
    tCQMHalfEdge2* he4 = e2->getMate(he3);
    if (he3->getFace() != f2) {
        he4 = he3;
        he3 = e2->getMate(he3);
    }
    
    // get the half-edge of "e3" in "f5" and "f4"
    tCQMHalfEdge2* he5 = e3->getHalfEdge();
    tCQMHalfEdge2* he6 = e3->getMate(he5);
    if (he5->getFace() != f5) {
        he6 = he5;
        he5 = e3->getMate(he5);
    }
    
    // get the half-edge of "e4" in "f4" and "f3"
    tCQMHalfEdge2* he7 = e4->getHalfEdge();
    tCQMHalfEdge2* he8 = e4->getMate(he7);
    if (he7->getFace() != f4) {
        he8 = he7;
        he7 = e4->getMate(he7);
    }
    
    // get all possible combinations of shared edges
    tCQMEdge2* e14 = f1->getCommonEdge(f4);
    tCQMEdge2* e15 = f1->getCommonEdge(f5);
    tCQMEdge2* e24 = f2->getCommonEdge(f4);
    tCQMEdge2* e25 = f2->getCommonEdge(f5);
    
    // compute polygon vertices
    tCQMIndVertex2 p1, p2, p3, p4, p5, p6;
    bool ec1, ec2, ec3, ec4, ec5, ec6;
    if (he4->getNext() == he8) {
        // "v2" is the left child of "v3"
        if ((e24 != 0) && (e15 != 0)) {
            p1 = getQuadVertex(vs, he5->getPrev()->getVertex());
            p2 = getQuadVertex(vs, he4->getPrev()->getVertex());
            p3 = getQuadVertex(vs, he1->getNext()->getVertex());
            p4 = getQuadVertex(vs, he3->getVertex());
            
            ec1 = he5->getPrev()->getEdge()->isConstrained();
            ec3 = he1->getNext()->getEdge()->isConstrained();
            
            tCQMIndVertex2 sp(quadrangulateDegTriangle(p1,p2,p3,p4,ec1,ec3,lq));
            v3->setSteinerPoint(sp,he8->getNext());
            return;
        }
        else {
            if (e14 != 0) {
                p1 = getQuadVertex(vs, he4->getPrev()->getVertex());
                p2 = getQuadVertex(vs, he8->getVertex());
                p3 = getQuadVertex(vs, he3->getNext()->getVertex());
                p4 = getQuadVertex(vs, he1->getNext()->getVertex());
                p5 = getQuadVertex(vs, he5->getNext()->getVertex());
                
                ec1 = false;
                ec2 = false;
                ec3 = he3->getNext()->getEdge()->isConstrained();
                ec4 = he1->getNext()->getEdge()->isConstrained();
                
                if (v5->hasSteinerPoint()) {
                    p6 = v5->getSteinerPoint();
                    ec5 = false;
                    ec6 = false;
                }
                else {
                    p6 = getQuadVertex(vs, he5->getPrev()->getVertex());
                    ec5 = he5->getNext()->getEdge()->isConstrained();
                    ec6 = he5->getPrev()->getEdge()->isConstrained();
                }
            }
            else if (e15 != 0) {
                p1 = getQuadVertex(vs, he4->getPrev()->getVertex());
                p2 = getQuadVertex(vs, he8->getVertex());
                p3 = getQuadVertex(vs, he3->getNext()->getVertex());
                p4 = getQuadVertex(vs, he1->getNext()->getVertex());
                p5 = getQuadVertex(vs, he5->getPrev()->getVertex());
                p6 = getQuadVertex(vs, he7->getPrev()->getVertex());
                
                ec1 = false;
                ec2 = false;
                ec3 = he3->getNext()->getEdge()->isConstrained();
                ec4 = he1->getNext()->getEdge()->isConstrained();
                ec5 = he5->getPrev()->getEdge()->isConstrained();
                ec6 = he7->getPrev()->getEdge()->isConstrained();
            }
            else if (e25 != 0) {
                p1 = getQuadVertex(vs, he4->getPrev()->getVertex());
                p2 = getQuadVertex(vs, he8->getVertex());
                p3 = getQuadVertex(vs, he1->getNext()->getVertex());
                
                if (v1->hasSteinerPoint()) {
                    p4 = v1->getSteinerPoint();
                    ec3 = false;
                    ec4 = false;
                }
                else {
                    p4 = getQuadVertex(vs, he1->getPrev()->getVertex());
                    ec3 = he1->getNext()->getEdge()->isConstrained();
                    ec4 = he1->getPrev()->getEdge()->isConstrained();
                }
                
                p5 = getQuadVertex(vs, he5->getPrev()->getVertex());
                p6 = getQuadVertex(vs, he7->getPrev()->getVertex());
                
                ec1 = false;
                ec2 = false;
                ec5 = he5->getPrev()->getEdge()->isConstrained();
                ec6 = he7->getPrev()->getEdge()->isConstrained();
            }
            else {  // e24 != 0
                p1 = getQuadVertex(vs, he4->getPrev()->getVertex());
                p2 = getQuadVertex(vs, he8->getVertex());
                p3 = getQuadVertex(vs, he1->getNext()->getVertex());
                
                if (v1->hasSteinerPoint()) {
                    p4 = v1->getSteinerPoint();
                    ec3 = false;
                    ec4 = false;
                }
                else {
                    p4 = getQuadVertex(vs, he1->getPrev()->getVertex());
                    ec3 = he1->getNext()->getEdge()->isConstrained();
                    ec4 = he1->getPrev()->getEdge()->isConstrained();
                }
                
                p5 = getQuadVertex(vs, he5->getNext()->getVertex());
                
                if (v5->hasSteinerPoint()) {
                    p6 = v5->getSteinerPoint();
                    ec5 = false;
                    ec6 = false;
                }
                else {
                    p6 = getQuadVertex(vs, he5->getPrev()->getVertex());
                    ec5 = he5->getNext()->getEdge()->isConstrained();
                    ec6 = he5->getPrev()->getEdge()->isConstrained();
                }
                
                ec1 = false;
                ec2 = false;
            }
        }
    }
    else {
        // "v6" is the right child of "v3"
        if ((e24 != 0) && (e15 != 0)) {
            p1 = getQuadVertex(vs, he1->getPrev()->getVertex());
            p2 = getQuadVertex(vs, he8->getPrev()->getVertex());
            p3 = getQuadVertex(vs, he5->getNext()->getVertex());
            p4 = getQuadVertex(vs, he4->getVertex());
            
            ec1 = he1->getPrev()->getEdge()->isConstrained();
            ec3 = he5->getNext()->getEdge()->isConstrained();
            
            tCQMIndVertex2 sp(quadrangulateDegTriangle(p1,p2,p3,p4,ec1,ec3,lq));
            v3->setSteinerPoint(sp,he4->getNext());
            return;
        }
        else {
            if (e14 != 0) {
                p1 = getQuadVertex(vs, he8->getPrev()->getVertex());
                p2 = getQuadVertex(vs, he4->getVertex());
                p3 = getQuadVertex(vs, he5->getNext()->getVertex());
                
                if (v5->hasSteinerPoint()) {
                    p4 = v5->getSteinerPoint();
                    ec3 = false;
                    ec4 = false;
                }
                else {
                    p4 = getQuadVertex(vs, he5->getPrev()->getVertex());
                    ec3 = he5->getNext()->getEdge()->isConstrained();
                    ec4 = he5->getPrev()->getEdge()->isConstrained();
                }
                
                p5 = getQuadVertex(vs, he1->getPrev()->getVertex());
                p6 = getQuadVertex(vs, he3->getPrev()->getVertex());
                
                ec1 = false;
                ec2 = false;
                ec5 = he1->getPrev()->getEdge()->isConstrained();
                ec6 = he3->getPrev()->getEdge()->isConstrained();
            }
            else if (e15 != 0) {
                p1 = getQuadVertex(vs, he8->getPrev()->getVertex());
                p2 = getQuadVertex(vs, he4->getVertex());
                p3 = getQuadVertex(vs, he7->getNext()->getVertex());
                p4 = getQuadVertex(vs, he5->getNext()->getVertex());
                p5 = getQuadVertex(vs, he1->getPrev()->getVertex());
                p6 = getQuadVertex(vs, he3->getPrev()->getVertex());
                
                ec1 = false;
                ec2 = false;
                ec3 = he7->getNext()->getEdge()->isConstrained();
                ec4 = he5->getNext()->getEdge()->isConstrained();
                ec5 = he1->getPrev()->getEdge()->isConstrained();
                ec6 = he3->getPrev()->getEdge()->isConstrained();
            }
            else if (e25 != 0) {
                p1 = getQuadVertex(vs, he8->getPrev()->getVertex());
                p2 = getQuadVertex(vs, he4->getVertex());
                p3 = getQuadVertex(vs, he7->getNext()->getVertex());
                p4 = getQuadVertex(vs, he5->getNext()->getVertex());
                p5 = getQuadVertex(vs, he1->getNext()->getVertex());
                
                if (v1->hasSteinerPoint()) {
                    p6 = v1->getSteinerPoint();
                    ec5 = false;
                    ec6 = false;
                }
                else {
                    p6 = getQuadVertex(vs, he1->getPrev()->getVertex());
                    ec5 = he1->getNext()->getEdge()->isConstrained();
                    ec6 = he1->getPrev()->getEdge()->isConstrained();
                }
                
                ec1 = false;
                ec2 = false;
                ec3 = he7->getNext()->getEdge()->isConstrained();
                ec4 = he5->getNext()->getEdge()->isConstrained();
            }
            else {  // e24 != 0
                p1 = getQuadVertex(vs, he8->getPrev()->getVertex());
                p2 = getQuadVertex(vs, he4->getVertex());
                p3 = getQuadVertex(vs, he5->getNext()->getVertex());
                
                if (v5->hasSteinerPoint()) {
                    p4 = v5->getSteinerPoint();
                    ec3 = false;
                    ec4 = false;
                }
                else {
                    p4 = getQuadVertex(vs, he5->getPrev()->getVertex());
                    ec3 = he5->getNext()->getEdge()->isConstrained();
                    ec4 = he5->getPrev()->getEdge()->isConstrained();
                }
                
                p5 = getQuadVertex(vs, he1->getNext()->getVertex());
                
                if (v1->hasSteinerPoint()) {
                    p6 = v1->getSteinerPoint();
                    ec5 = false;
                    ec6 = false;
                }
                else {
                    p6 = getQuadVertex(vs, he1->getPrev()->getVertex());
                    ec5 = he5->getNext()->getEdge()->isConstrained();
                    ec6 = he5->getPrev()->getEdge()->isConstrained();
                }
                
                ec1 = false;
                ec2 = false;
            }
        }
    }
    
    // create an hexagon
    tCQMPolygon2* p = (tCQMPolygon2 *) new tCQMPolygon2(6);
    
    p->insert(0,p1);
    p->insert(1,p2);
    p->insert(2,p3);
    p->insert(3,p4);
    p->insert(4,p5);
    p->insert(5,p6);
    
    p->setConstraint(0,ec1);
    p->setConstraint(1,ec2);
    p->setConstraint(2,ec3);
    p->setConstraint(3,ec4);
    p->setConstraint(4,ec5);
    p->setConstraint(5,ec6);
    
    makeConvexQuadFromHex(p,vs,lq);
    
    delete p;
}

// -------------------------------------------------------------------
// Method makeConvexQuadFromSep()
// -------------------------------------------------------------------
bool
tCQMCompQuad::makeConvexQuadFromSep(
                                    tCQMSpanningTreeVertex* v1,
                                    tCQMSpanningTreeVertex* v2,
                                    tCQMSpanningTreeVertex* v3,
                                    tCQMSpanningTreeVertex* v4,
                                    tCQMSpanningTreeVertex* v5,
                                    QuadVertexSet& vs,
                                    std::list<tCQMQuadrilateral2*>* lq
                                    )
{
    // "v1" is the only child of "v2", which is a child of "v3". "v5" is
    // the only child of "v4", which is another child of "v3".
    
    
    // get  the mesh  faces corresponding  to the  dual graph  vertex in
    // "v1", "v2", "v3", "v4" and "v5".
    tCQMFace2* f1 = v1->getVertex()->getFace();
    tCQMFace2* f2 = v2->getVertex()->getFace();
    tCQMFace2* f3 = v3->getVertex()->getFace();
    tCQMFace2* f4 = v4->getVertex()->getFace();
    tCQMFace2* f5 = v5->getVertex()->getFace();
    
    // verify  if  "f1", "f2",  "f3",  "f4"  and  "f5" degenerate  to  a
    // triangle  or  a pentagon.  If  so,  we  convex quadrangulate  the
    // degenerate pentagon using other method than this.
    if ((f1->getCommonEdge(f4) != 0) || (f1->getCommonEdge(f5) != 0) ||
        (f2->getCommonEdge(f4) != 0) || (f2->getCommonEdge(f4) != 0)) {
        makeConvexQuadFromDegSep(v1,v2,v3,v4,v5,vs,lq);
    }
    else {
        
        // get the common edge of "f1" and "f2"
        tCQMEdge2* e1 = f1->getCommonEdge(f2);
        
        // get the common edge of "f2" and "f3"
        tCQMEdge2* e2 = f2->getCommonEdge(f3);
        
        // get the common edge of "f3" and "f4"
        tCQMEdge2* e3 = f3->getCommonEdge(f4);
        
        // get the common edge of "f4" and "f5"
        tCQMEdge2* e4 = f4->getCommonEdge(f5);
        
        // get the half-edges of "e1" in "f1" and "f2"
        tCQMHalfEdge2* he1 = e1->getHalfEdge();
        tCQMHalfEdge2* he2 = e1->getMate(he1);
        if (he1->getFace() != f1) {
            he2 = he1;
            he1 = e1->getMate(he1);
        }
        
        // get the half-edges of "e2" in "f2" and "f3"
        tCQMHalfEdge2* he3 = e2->getHalfEdge();
        tCQMHalfEdge2* he4 = e2->getMate(he3);
        if (he3->getFace() != f2) {
            he4 = he3;
            he3 = e2->getMate(he3);
        }
        
        // get the half-edge of "e3" in "f3" and "f4"
        tCQMHalfEdge2* he5 = e3->getHalfEdge();
        tCQMHalfEdge2* he6 = e3->getMate(he5);
        if (he5->getFace() != f3) {
            he6 = he5;
            he5 = e3->getMate(he5);
        }
        
        // get the half-edge of "e4" in "f4" and "f5"
        tCQMHalfEdge2* he7 = e4->getHalfEdge();
        tCQMHalfEdge2* he8 = e4->getMate(he7);
        if (he7->getFace() != f4) {
            he8 = he7;
            he7 = e4->getMate(he7);
        }
        
        // find out the  vertices of the septagon defined  by the union of
        // "f1", "f2", "f3", "f4" and "f5".
        tCQMIndVertex2 p1 = getQuadVertex(vs, he1->getVertex());
        tCQMIndVertex2 p6 = getQuadVertex(vs, he1->getNext()->getVertex());
        tCQMIndVertex2 p7;
        bool ec6, ec7;
        if (v1->hasSteinerPoint()) {
            p7 = v1->getSteinerPoint();
            ec6 = false;
            ec7 = false;
        }
        else {
            p7 = getQuadVertex(vs, he1->getPrev()->getVertex());
            ec6 = he1->getNext()->getEdge()->isConstrained();
            ec7 = he1->getPrev()->getEdge()->isConstrained();
        }
        
        // if "he2->next =  he3" then the first vertex  of the septagon is
        // the common vertex of "f1", "f2" and "f3".
        tCQMIndVertex2 p2, p3, p4, p5;
        bool ec1, ec2, ec3, ec4, ec5;
        if (he2->getNext() == he3) {
            if (he6->getNext() == he7) {
                // either "p1" or "p5" is the common vertex of the faces "f3",
                // "f4" and "f5".
                if (he3->getVertex() == he7->getVertex()) {
                    // "p1" is the common vertex  of all faces "f1", "f2", "f3",
                    // "f4" and "f5".
                    if (v5->hasSteinerPoint()) {
                        p2 = v5->getSteinerPoint();
                        ec1 = false;
                        ec2 = false;
                    }
                    else {
                        p2 = getQuadVertex(vs, he8->getPrev()->getVertex());
                        ec1 = he8->getNext()->getEdge()->isConstrained();
                        ec2 = he8->getPrev()->getEdge()->isConstrained();
                    }
                    
                    p3 = getQuadVertex(vs, he8->getVertex());
                    p4 = getQuadVertex(vs, he6->getVertex());
                    p5 = getQuadVertex(vs, he4->getVertex());
                    
                    ec3 = he7->getNext()->getEdge()->isConstrained();
                    ec4 = he5->getNext()->getEdge()->isConstrained();
                    ec5 = he3->getNext()->getEdge()->isConstrained();
                    
                    v3->setSteinerPoint(quadrangulateSep1(p1,p2,p3,p4,p5,p6,p7,
                                                          ec1,ec2,ec3,ec4,ec5,ec6,ec7,lq),p4.getPoint());
                }
                else {
                    // "p1" is  the common  vertex of faces  "f1" and  "f2", and
                    // "p2" is the common vertex of faces "f3", "f4" and "f5".
                    p2 = getQuadVertex(vs, he7->getVertex());
                    if (v5->hasSteinerPoint()) {
                        p3 = v5->getSteinerPoint();
                        ec2 = false;
                        ec3 = false;
                    }
                    else {
                        p3 = getQuadVertex(vs, he8->getPrev()->getVertex());
                        ec2 = he8->getNext()->getEdge()->isConstrained();
                        ec3 = he8->getPrev()->getEdge()->isConstrained();
                    }
                    
                    p4 = getQuadVertex(vs, he8->getVertex());
                    p5 = getQuadVertex(vs, he6->getVertex());
                    
                    ec1 = he4->getNext()->getEdge()->isConstrained();
                    ec4 = he6->getPrev()->getEdge()->isConstrained();
                    ec5 = he3->getNext()->getEdge()->isConstrained();
                    
                    v3->setSteinerPoint(quadrangulateSep4(p4,p5,p6,p7,p1,p2,p3,
                                                          ec4,ec5,ec6,ec7,ec1,ec2,ec3,lq),p1.getPoint());
                }
            }
            else {
                // either "p4" or "p5" is the common vertex of the faces "f3",
                // "f4" and "f5".
                if (he3->getVertex() == he5->getVertex()) {
                    // "p4" is  the common  vertex of the  faces "f3",  "f4" and
                    // "f5".
                    p2 = getQuadVertex(vs, he7->getVertex());
                    
                    if (v5->hasSteinerPoint()) {
                        p3 = v5->getSteinerPoint();
                        ec2 = false;
                        ec3 = false;
                    }
                    else {
                        p3 = getQuadVertex(vs, he8->getPrev()->getVertex());
                        ec2 = he8->getNext()->getEdge()->isConstrained();
                        ec3 = he8->getPrev()->getEdge()->isConstrained();
                    }
                    
                    p4 = getQuadVertex(vs, he6->getVertex());
                    p5 = getQuadVertex(vs, he4->getVertex());
                    
                    ec1 = he6->getNext()->getEdge()->isConstrained();
                    ec4 = he5->getNext()->getEdge()->isConstrained();
                    ec5 = he3->getNext()->getEdge()->isConstrained();
                    
                    v3->setSteinerPoint(quadrangulateSep2(p1,p2,p3,p4,p5,p6,p7,
                                                          ec1,ec2,ec3,ec4,ec5,ec6,ec7,lq),p4.getPoint());
                }
                else {
                    // "p5" is  the common  vertex of the  faces "f3",  "f4" and
                    // "f5".
                    p2 = getQuadVertex(vs, he5->getVertex());
                    p3 = getQuadVertex(vs, he7->getVertex());
                    
                    ec1 = he4->getNext()->getEdge()->isConstrained();
                    ec2 = he6->getNext()->getEdge()->isConstrained();
                    if (v5->hasSteinerPoint()) {
                        p4 = v5->getSteinerPoint();
                        ec3 = false;
                        ec4 = false;
                    }
                    else {
                        p4 = getQuadVertex(vs, he8->getPrev()->getVertex());
                        ec3 = he8->getNext()->getEdge()->isConstrained();
                        ec4 = he8->getPrev()->getEdge()->isConstrained();
                    }
                    
                    p5 = getQuadVertex(vs, he6->getVertex());
                    ec5 = he3->getNext()->getEdge()->isConstrained();
                    
                    v3->setSteinerPoint(quadrangulateSep2(p5,p6,p7,p1,p2,p3,p4,
                                                          ec5,ec6,ec7,ec1,ec2,ec3,ec4,lq),p1.getPoint());
                }
            }
        }
        else {
            // "p1" is  the common vertex of  faces "f1" and "f2"  and it is
            // not incident to "f3".
            if (he6->getNext() == he7) {
                // either "p2"  or "p6"  is the common  vertex of  faces "f3",
                // "f4" and "f5".
                if (he3->getVertex() == he7->getVertex()) {
                    // "p2" is the common vertex of faces "f3", "f4" and "f5".
                    p2 = getQuadVertex(vs, he3->getVertex());
                    
                    ec1 = he2->getNext()->getEdge()->isConstrained();
                    if (v5->hasSteinerPoint()) {
                        p3 = v5->getSteinerPoint();
                        ec2 = false;
                        ec3 = false;
                    }
                    else {
                        p3 = getQuadVertex(vs, he8->getPrev()->getVertex());
                        ec2 = he8->getNext()->getEdge()->isConstrained();	  
                        ec3 = he8->getPrev()->getEdge()->isConstrained();
                    }
                    
                    p4 = getQuadVertex(vs, he8->getVertex());
                    p5 = getQuadVertex(vs, he6->getVertex());
                    
                    ec4 = he7->getNext()->getEdge()->isConstrained();
                    ec5 = he5->getNext()->getEdge()->isConstrained();
                    
                    v3->setSteinerPoint(quadrangulateSep3(p1,p2,p3,p4,p5,p6,p7,
                                                          ec1,ec2,ec3,ec4,ec5,ec6,ec7,lq),p5.getPoint());
                }
                else {
                    // "p6" is the common vertex of faces "f3", "f4" and "f5".
                    p2 = getQuadVertex(vs, he3->getVertex());
                    p3 = getQuadVertex(vs, he5->getVertex());
                    
                    ec1 = he2->getNext()->getEdge()->isConstrained();
                    ec2 = he4->getNext()->getEdge()->isConstrained();
                    
                    if (v5->hasSteinerPoint()) {
                        p4 = v5->getSteinerPoint();
                        ec3 = false;
                        ec4 = false;
                    }
                    else {
                        p4 = getQuadVertex(vs, he8->getPrev()->getVertex());
                        ec3 = he8->getNext()->getEdge()->isConstrained();	  
                        ec4 = he8->getPrev()->getEdge()->isConstrained();
                    }
                    
                    p5 = getQuadVertex(vs, he8->getVertex());
                    ec5 = he7->getNext()->getEdge()->isConstrained();
                    
                    v3->setSteinerPoint(quadrangulateSep3(p5,p6,p7,p1,p2,p3,p4,
                                                          ec5,ec6,ec7,ec1,ec2,ec3,ec4,lq),p2.getPoint());
                }
            }
            else {
                // either "p2" or "p6" is the common vertex of all faces
                if (he3->getVertex() == he5->getVertex()) {
                    // "p2" is the common vertex of all faces
                    p2 = getQuadVertex(vs, he3->getVertex());
                    p3 = getQuadVertex(vs, he7->getVertex());
                    
                    ec1 = he2->getNext()->getEdge()->isConstrained();
                    ec2 = he6->getNext()->getEdge()->isConstrained();
                    
                    if (v5->hasSteinerPoint()) {
                        p4 = v5->getSteinerPoint();
                        ec3 = false;
                        ec4 = false;
                    }
                    else {
                        p4 = getQuadVertex(vs, he8->getPrev()->getVertex());
                        ec3 = he8->getNext()->getEdge()->isConstrained(); 
                        ec4 = he8->getPrev()->getEdge()->isConstrained();
                    }
                    
                    p5 = getQuadVertex(vs, he8->getVertex());
                    ec5 = he5->getNext()->getEdge()->isConstrained();
                    
                    v3->setSteinerPoint(quadrangulateSep4(p1,p2,p3,p4,p5,p6,p7,
                                                          ec1,ec2,ec3,ec4,ec5,ec6,ec7,lq),p5.getPoint());
                }
                else {
                    // "p6" is the common vertex of all faces
                    p2 = getQuadVertex(vs, he3->getVertex());
                    p3 = getQuadVertex(vs, he5->getVertex());
                    p4 = getQuadVertex(vs, he7->getVertex());
                    
                    ec1 = he2->getNext()->getEdge()->isConstrained();
                    ec2 = he4->getNext()->getEdge()->isConstrained();
                    ec3 = he6->getNext()->getEdge()->isConstrained();
                    
                    if (v5->hasSteinerPoint()) {
                        p5 = v5->getSteinerPoint();
                        ec4 = false;
                        ec5 = false;
                    }
                    else {
                        p5 = getQuadVertex(vs, he8->getPrev()->getVertex());
                        ec4 = he8->getNext()->getEdge()->isConstrained();
                        ec5 = he8->getPrev()->getEdge()->isConstrained();
                    }
                    
                    v3->setSteinerPoint(quadrangulateSep1(p6,p7,p1,p2,p3,p4,p5,
                                                          ec6,ec7,ec1,ec2,ec3,ec4,ec5,lq),p2.getPoint());
                }
            }
        }
    }
    
    return false;
}


// -------------------------------------------------------------------
// Method computeRootVertices()
// 
// It  computes  the  "root"  vertices  of  the  spanning  trees  that
// correspond  to each  "simple"  region  of the  input  domain. By  a
// "simple"  region,  we  mean   a  simple  polygon  with  or  without
// roles. The boundary of each region consists of a set of constrained
// edges of  the input triangulation,  and no constrained edge  of the
// triangulation  can  intersect  the  interior  of  the  region.  The
// triangulation  of  each simple  region  will  be  converted into  a
// quadrangulation.
//
// This method returns a list of "half-edges" whose faces are the dual
// of the  vertices of the  dual graph that  will be the roots  of the
// spanning trees.
//
// -------------------------------------------------------------------
void
tCQMCompQuad::computeRootVertices(std::list<tCQMHalfEdge2*>& lhedges)
const
{
    // create the region graph.
    tCQMRegionGraph* rg = new tCQMRegionGraph(m_mesh->getDualGraph());
    
    rg->getListOfHalfEdges(lhedges);
    
    delete rg;
}

// -------------------------------------------------------------------
// Method computeConvexQuadrangulation()
// 
// It computes a quadrangulation  by converting a triangular mesh into
// a quadrangular one consisting only of convex quadrilaterals.
// -------------------------------------------------------------------

std::list<tCQMQuadrilateral2*>*
tCQMCompQuad::computeConvexQuadrangulation(string fname)
{
    std::cout << "Begin computing quad" << endl;
    // initialize number of Steiner points
    m_spoints = 0;
    spt_filename = fname;
    // create empty list of quadrilaterals
    std::list<tCQMQuadrilateral2*>* lquads = (std::list<tCQMQuadrilateral2*>*)
    new std::list<tCQMQuadrilateral2*>();
    
    // create   set  of  pair   <vertex,id>  for   uniquely  identifying
    // quadrilateral vertices.
    QuadVertexSet vs;
    
    
    // compute the list of vertices of the dual graph that correspond to
    // the "root"  vertices of the spanning  trees that will  be used to
    // compute the quadrangulation of  each "simple" region of the input
    // domain.
    std::list<tCQMHalfEdge2*> lhedges;
    computeRootVertices(lhedges); 
    
    // for each "root" vertex "v", create a spanning tree rooted at "v",
    // and   then   generate  a   quadrangulation   by  converting   the
    // triangulation   represented   by  the   spanning   tree  into   a
    // quadrangulation.
    while (!lhedges.empty()) {
        
        // get one contour
        tCQMHalfEdge2* rvhe = lhedges.front();
        
        // remove it from the list of contours
        lhedges.pop_front();
        
        // get dual vertex of the face containing the half-edge.
        tCQMDualGraphVertex* v = rvhe->getFace()->getDualVertex();
        
        // compute spanning tree inside contour
        tCQMSpanningTree* sp = this->getSpanningTree(v);
        
        // Debugging. -- J.L.
        sp->dump(fname);
        // get set of vertices per level of the spanning tree
        std::list<std::list<tCQMSpanningTreeVertex*>*>* spvl = sp->getVertPerLevel();
        
        // create a vector of set of vertices per level
        int vsetsize = spvl->size();
        vector<std::list<tCQMSpanningTreeVertex*>*> vset(vsetsize);
        
        // put the several sets into the vector
        std::list<std::list<tCQMSpanningTreeVertex*>*>::iterator vl;
        int i = 0;
        for (vl = spvl->begin(); vl != spvl->end(); ++vl) {
            vset[i] = *vl;
            ++i;
        }
        
        // loop over each set of vertex per level
        for (i=0; i<vsetsize; i++) {
            std::cout << "level: " << i << endl;
            // if there  is at least  one vertex in  the set of  vertices at
            // level "i", process it.
            if (!vset[i]->empty()) {
                tCQMSpanningTreeVertex* v;
                
                // process vertices  at level "i" which  are quadrilaterals or
                // pentagons.  After   this  loop,  there  will   be  no  more
                // quadrilaterals and pentagons at level "i".
                while ((v = findOneQuadrilateral(*(vset[i]))) != 0) {
                    makeQuad(v,vs,lquads);
                    // remove "v" from the set of vertices at level "i" and from
                    // the spanning tree.
                    vset[i]->remove(v);
                    sp->delLeaf(v);
                }
                
                while ((v = findOnePentagon(*(vset[i]))) != 0) {
                    makeConvexQuadFromPent(v,vs,lquads);
                }
                
                // process  vertices whose  parent  has degree  2. After  this
                // loop, there will  be no more vertices of  degree 2 at level
                // "i+1".
                while ((v = getChildOfDegreeNParent(2,*(vset[i]))) != 0) {
                    // get the parent of "v", "vpar"
                    tCQMSpanningTreeVertex* vpar = v->getParent();
                    
                    // get parent of "vpar", "vparvpar"
                    tCQMSpanningTreeVertex* vparvpar = vpar->getParent();
                    
                    // auxiliary variable
                    bool rmvpar2 = false;
                    
                    // if "vpar" has a  parent, then the triangles corresponding
                    // to "v" and "vpar" make a quadrilateral. Otherwise, "vpar"
                    // has other child than "v" and it is the tree root.
                    if (vparvpar != 0) {
                        // if "vpar"  is a  quadrilateral already, then  check one
                        // level up.
                        
                        if (vpar->isQuadrilateral()) {
                            // check if  "vparvpar" is also a  quadrilateral. If so,
                            // ignore  it  and  convex  quadrangulate  the  pentagon
                            // defined by "v" and "vpar".
                            
                            if (vparvpar->isQuadrilateral()) {
                                // convex  quadrangulate the  pentagon defined  by "v"
                                // and "vpar".
                                makeConvexQuadFromPent(v,vpar,vs,lquads);
                                
                                // remove "v"  from the set  of vertices at  level "i"
                                // and from the spanning tree.
                                vset[i]->remove(v);
                                sp->delLeaf(v);
                            }
                            else {
                                // if  "vparvpar"  has  no  other child  than  "vpar",
                                // convex  quadrangulate the  hexagon defined  by "v",
                                // "vpar"  and  "vparvpar"  with  at  most  3  Steiner
                                // points.
                                if (vparvpar->getNumChildren() == 1) {
                                    makeConvexQuadFromHex(v,vpar,vparvpar,vs,lquads);
                                    
                                    // remove "v" from the  set of vertices at level "i"
                                    // and from the spanning tree.
                                    vset[i]->remove(v);
                                    sp->delLeaf(v);
                                    
                                    // remove "vpar"  from the set of  vertices at level
                                    // "i+1".
                                    vset[i+1]->remove(vpar);
                                    sp->delLeaf(vpar);
                                    
                                    // if "vparvpar" can be also removed, do it
                                    vset[i+2]->remove(vparvpar);
                                    sp->delLeaf(vparvpar);
                                }
                                else {
                                    // convex quadrangulate the  pentagon defined by "v"
                                    // and "vpar".
                                    makeConvexQuadFromPent(v,vpar,vs,lquads);
                                    
                                    // remove "v" from the  set of vertices at level "i"
                                    // and from the spanning tree.
                                    vset[i]->remove(v);
                                    sp->delLeaf(v);
                                }
                            }
                        }
                        else {
                            // "vpar" is  not a  quadrilateral, and we  consider the
                            // quadrilateral defined by "v" and "vpar".
                            
                            // if "v"  and "vpar" form a  convex quadrilateral, make
                            // the   quadrilateral.   Otherwise,  we   consider  the
                            // subtree rooted at the parent of "vpar".
                            if (!makeQuadIfConvex(v,vpar,vs,lquads)) {
                                
                                // if "vparvpar" is  a quadrilateral, then "v", "vpar"
                                // and  "vparvpar" define  an hexagon,  and  we convex
                                // quadrangulate this  hexagon with at  most 3 Steiner
                                // points.
                                if (vparvpar->isQuadrilateral()) {
                                    makeConvexQuadFromHex(v,vpar,vparvpar,vs,lquads);
                                    
                                    // remove "v" from the  set of vertices at level "i"
                                    // and from the spanning tree.
                                    vset[i]->remove(v);
                                    sp->delLeaf(v);
                                    
                                    // remove "vpar"  from the set of  vertices at level
                                    // "i+1".
                                    vset[i+1]->remove(vpar);
                                    sp->delLeaf(vpar);
                                    
                                    // if "vparvpar" can be also removed, do it
                                    vset[i+2]->remove(vparvpar);
                                    sp->delLeaf(vparvpar);
                                }
                                else {
                                    // "vpar" and "vparvpar" are not quadrilaterals, and
                                    // if  "vparvpar" has  no child  other  than "vpar",
                                    // then  we  can  convex  quadrangulate  the  region
                                    // corresponding to "v", "vpar", and "vparvpar" with
                                    // at most 3 Steiner points.
                                    if (vparvpar->getNumChildren() == 1) {
                                        // we convex quadrangulate the pentagon defined by
                                        // "v",  "vpar"  and "vparvpar"  using  at most  3
                                        // Steiner points.
                                        rmvpar2 = makeConvexQuadFromPent(v,vpar,vparvpar,
                                                                         vs,lquads);
                                        
                                        // remove "v"  from the  set of vertices  at level
                                        // "i" and from the spanning tree.
                                        vset[i]->remove(v);
                                        sp->delLeaf(v);
                                        
                                        // remove "vpar" from the set of vertices at level
                                        // "i+1".
                                        vset[i+1]->remove(vpar);
                                        sp->delLeaf(vpar);
                                        
                                        // if "vparvpar" can be also removed, do it
                                        if (rmvpar2) {
                                            vset[i+2]->remove(vparvpar);
                                            sp->delLeaf(vparvpar);
                                        }
                                    }
                                    else {
                                        // "vparvpar"   has  other   child   than  "vpar",
                                        // "vparsib".
                                        tCQMSpanningTreeVertex* vparsib = 
                                        vparvpar->getSibling(vpar);
                                        
                                        // if   "vparsib"   is   a  pentagon   we   convex
                                        // quadrangulate  it  and  we  are left  with  one
                                        // triangle.  Later,  we convex quabdrangulate the
                                        // hexagon defined by the remaining triangle, "v",
                                        // "vpar" and "vparvpar".
                                        if (vparsib->isPentagon()) {
                                            // convex quadrangulate a pentagon
                                            makeConvexQuadFromPent(vparsib,vs,lquads);
                                            
                                            // convex  quadrangulate the hexagon  defined by
                                            // "v",  "vpar",  "vparvpar"  and the  remaining
                                            // triangle that replaced "vparsib".
                                            makeConvexQuadFromHex(v,vpar,vparsib,vparvpar,
                                                                  vs,lquads);
                                            
                                            // remove "v" from the  set of vertices at level
                                            // "i" and from the spanning tree.
                                            vset[i]->remove(v);
                                            sp->delLeaf(v);
                                            
                                            // remove  "vpar" from  the set  of  vertices at
                                            // level "i+1".
                                            vset[i+1]->remove(vpar);
                                            sp->delLeaf(vpar);
                                            
                                            // remove "vparsib" from  the set of vertices at
                                            // level "i+1".
                                            vset[i+1]->remove(vparsib);
                                            sp->delLeaf(vparsib);
                                            
                                            // remove "vparvpar" from the set of vertices at
                                            // level "i+2".
                                            vset[i+2]->remove(vparvpar);
                                            sp->delLeaf(vparvpar);
                                        }
                                        else if (vparsib->isQuadrilateral()) {
                                            // "vparsib"  is a  quadrilateral
                                            if (vparsib->getDegree() == 2) {
                                                // "vparsib"  has  a   child,  and  we  convex
                                                // quadrangulate   the  pentagon   defined  by
                                                // "vparsib" and its child.
                                                // child of "vparsib" is a triangle
                                                tCQMSpanningTreeVertex* vparsibchild = 
                                                vparsib->getChild(0);
                                                makeConvexQuadFromPent(vparsibchild,vparsib,vs,
                                                                       lquads);
                                                
                                                // remove  "vparsibchild"   from  the  set  of
                                                // vertices at level "i" and from the spanning
                                                // tree.
                                                vset[i]->remove(vparsibchild);
                                                sp->delLeaf(vparsibchild);
                                                
                                                // convex quadrangulate the hexagon defined by
                                                // "v",  "vpar", "vparvpar" and  the remaining
                                                // triangle that replaced "vparsib".
                                                makeConvexQuadFromHex(v,vpar,vparsib,vparvpar,
                                                                      vs,lquads);
                                                
                                                // remove  "v"  from the  set  of vertices  at
                                                // level "i" and from the spanning tree.
                                                vset[i]->remove(v);
                                                sp->delLeaf(v);
                                                
                                                // remove "vpar"  from the set  of vertices at
                                                // level "i+1".
                                                vset[i+1]->remove(vpar);
                                                sp->delLeaf(vpar);
                                                
                                                // remove "vparsib"  from the set  of vertices
                                                // at level "i+1".
                                                vset[i+1]->remove(vparsib);
                                                sp->delLeaf(vparsib);
                                                
                                                // remove "vparvpar" from  the set of vertices
                                                // at level "i+2".
                                                vset[i+2]->remove(vparvpar);
                                                sp->delLeaf(vparvpar);
                                            }
                                            else {
                                                // "vparsib"  has no  child, and  we  create a
                                                // quadrilateral corresponding to it.
                                                makeQuad(vparsib,vs,lquads);
                                                
                                                // remove  "v"  from the  set  of vertices  at
                                                // level "i" and from the spanning tree.
                                                vset[i+1]->remove(vparsib);
                                                sp->delLeaf(vparsib);
                                                
                                                // convex  quadrangulate the  pentagon defined
                                                // by "v", "vpar" and "vparvpar".
                                                rmvpar2 = makeConvexQuadFromPent(v,vpar,
                                                                                 vparvpar,vs,lquads);
                                                
                                                // remove  "v"  from the  set  of vertices  at
                                                // level "i" and from the spanning tree.
                                                vset[i]->remove(v);
                                                sp->delLeaf(v);
                                                
                                                // remove "vpar"  from the set  of vertices at
                                                // level "i+1".
                                                vset[i+1]->remove(vpar);
                                                sp->delLeaf(vpar);
                                                
                                                // if "vparvpar" can be also removed, do it
                                                if (rmvpar2) {
                                                    vset[i+2]->remove(vparvpar);
                                                    sp->delLeaf(vparvpar);
                                                }
                                            }
                                        }
                                        else {
                                            // "vparsib" is a triangle and, if "vparsib" has
                                            // no child, we convex quadrangulate the hexagon
                                            // defined  by the  faces  corresponding to  the
                                            // dual  vertices associated  with  "v", "vpar",
                                            // "vparvpar" and "vparsib".
                                            if (vparsib->getDegree() == 1) {
                                                // convex quadrangulate the hexagon defined by
                                                // "v", "vpar", "vparsib" and "vparvpar".
                                                makeConvexQuadFromHex(v,vpar,vparsib,vparvpar,
                                                                      vs,lquads);
                                                
                                                // remove  "v"  from the  set  of vertices  at
                                                // level "i" and from the spanning tree.
                                                vset[i]->remove(v);
                                                sp->delLeaf(v);
                                                
                                                // remove "vpar"  from the set  of vertices at
                                                // level "i+1".
                                                vset[i+1]->remove(vpar);
                                                sp->delLeaf(vpar);
                                                
                                                // remove "vparsib"  from the set  of vertices
                                                // at level "i+1".
                                                vset[i+1]->remove(vparsib);
                                                sp->delLeaf(vparsib);
                                                
                                                // remove "vparvpar" from  the set of vertices
                                                // at level "i+2".
                                                vset[i+2]->remove(vparvpar);
                                                sp->delLeaf(vparvpar);
                                            }
                                            else if (vparsib->getDegree() == 2) {
                                                // child of "vparsib" is a triangle
                                                tCQMSpanningTreeVertex* vparsibchild = 
                                                vparsib->getChild(0);
                                                
                                                // convex  quadrangulate the  septagon defined
                                                // by  "v",   "vpar",  "vparvpar",  "vparsib",
                                                // "vparsibchild".
                                                makeConvexQuadFromSep(v,vpar,vparvpar,vparsib,
                                                                      vparsibchild,vs,lquads);
                                                
                                                // remove  "vparsibchild"   from  the  set  of
                                                // vertices at level "i".
                                                vset[i]->remove(vparsibchild);
                                                sp->delLeaf(vparsibchild);
                                                
                                                // remove  "v"  from the  set  of vertices  at
                                                // level "i".
                                                vset[i]->remove(v);
                                                sp->delLeaf(v);
                                                
                                                // remove "vpar"  from the set  of vertices at
                                                // level "i+1".
                                                vset[i+1]->remove(vpar);
                                                sp->delLeaf(vpar);
                                                
                                                // remove "vparsib"  from the set  of vertices
                                                // at level "i+1".
                                                vset[i+1]->remove(vparsib);
                                                sp->delLeaf(vparsib);
                                            }
                                            else {
                                                // convex  quadrangulate the  pentagon defined
                                                // by "vparsib" and  its children, then convex
                                                // quadrangulate  the hexagon defined  by "v",
                                                // "vpar",   "vparvpar"   and   the   triangle
                                                // resulting       from      the      pentagon
                                                // quadrangulation.
                                                tCQMSpanningTreeVertex* vparsibch1 = 
                                                vparsib->getChild(0);
                                                tCQMSpanningTreeVertex* vparsibch2 = 
                                                vparsib->getChild(1);
                                                
                                                // convex quadrangulate a pentagon
                                                try {
                                                    makeQuad(vparsibch1,vparsibch2,vparsib,vs,
                                                             lquads);
                                                }
                                                catch (const tExceptionObject& xpt) {
                                                    treatException(xpt);
                                                    exit(0);
                                                }
                                                
                                                // remove children  of "vparsib" from  the set
                                                // of vertices at level "i".
                                                vset[i]->remove(vparsibch1);
                                                sp->delLeaf(vparsibch1);
                                                vset[i]->remove(vparsibch2);
                                                sp->delLeaf(vparsibch2);
                                                
                                                // convex quadrangulate an hexagon
                                                makeConvexQuadFromHex(v,vpar,vparsib,vparvpar,
                                                                      vs,lquads);
                                                
                                                // remove  "v"  from the  set  of vertices  at
                                                // level "i" and from the spanning tree.
                                                vset[i]->remove(v);
                                                sp->delLeaf(v);
                                                
                                                // remove "vpar"  from the set  of vertices at
                                                // level "i+1".
                                                vset[i+1]->remove(vpar);
                                                sp->delLeaf(vpar);
                                                
                                                // remove "vparsib"  from the set  of vertices
                                                // at level "i+1".
                                                vset[i+1]->remove(vparsib);
                                                sp->delLeaf(vparsib);
                                                
                                                vset[i+2]->remove(vparvpar);
                                                sp->delLeaf(vparvpar);
                                            }
                                        }
                                    }
                                }
                            }
                            else {
                                // remove "v"  from the set  of vertices at  level "i"
                                // and from the spanning tree.
                                vset[i]->remove(v);
                                sp->delLeaf(v);
                                
                                // remove "vpar" from the set of vertices at level "i+1"
                                vset[i+1]->remove(vpar);
                                sp->delLeaf(vpar);
                            }
                        }
                    }
                    else {
                        // "vpar" is the  tree root and has other  child than "v",
                        // "vsib".
                        // get "vsib"
                        tCQMSpanningTreeVertex* vsib = vpar->getSibling(v);
                        
                        // create  a Steiner  point  inside "vpar",  and make  two
                        // quadrilaterals  using  "v",  "vsib", and  that  Steiner
                        // point.
                        try {
                            makeQuad(v,vsib,vpar,vs,lquads);
                        }
                        catch (const tExceptionObject& xpt) {
                            treatException(xpt);
                            exit(0);
                        }
                        
                        // remove "v"  from the set  of vertices at level  "i" and
                        // from the spanning tree.
                        vset[i]->remove(v);
                        sp->delLeaf(v);
                        
                        // remove sibling  from set of  vertices at level  "i" and
                        // from the spanning tree.
                        vset[i]->remove(vsib);
                        sp->delLeaf(vsib);
                    }
                }   // end  of  while to  process  vertices  whose parent  has
                // degree 2.
                
                // process  vertices whose  parent  has degree  3. After  this
                // loop, there will be no vertex of degree 3 at level "i+1".
                while ((v = getChildOfDegreeNParent(3, *(vset[i]))) != 0) {
                    // get the parent of "v", "vpar"
                    tCQMSpanningTreeVertex* vpar = v->getParent();
                    
                    // get the parent  of "vpar", "vparvpar".
                    tCQMSpanningTreeVertex* vparvpar = vpar->getParent();
                    
                    // if  "vparvpar" has  only one  child, it  can be  either a
                    // triangle or a quadrilateral.
                    if (vparvpar->getNumChildren() == 1) {
                        // get sibling of "v"
                        tCQMSpanningTreeVertex* vsib = vpar->getSibling(v);
                        
                        if (vparvpar->isQuadrilateral()) {
                            // convex  quadrangulate the  pentagon  defined by  "v",
                            // "vsib" and "vpar".
                            try {
                                makeQuad(v,vsib,vpar,vs,lquads);
                            }
                            catch (const tExceptionObject& xpt) {
                                treatException(xpt);
                                exit(0);
                            }
                            
                            // remove "v" from the set  of vertices at level "i" and
                            // from the spanning tree.
                            vset[i]->remove(v);
                            sp->delLeaf(v);
                            
                            // remove "vsib"  from the set of vertices  at level "i"
                            // and from the spanning tree.
                            vset[i]->remove(vsib);
                            sp->delLeaf(vsib);
                        }
                        else {
                            // convex quadrangulate the  hexagon defined by "v", the
                            // sibling of "v", "vpar"  and "vparvpar" with at most 3
                            // steiner points.
                            makeConvexQuadFromHex(v,vpar,vsib,vparvpar,vs,lquads);
                            
                            // remove "v" from the set  of vertices at level "i" and
                            // from the spanning tree.
                            vset[i]->remove(v);
                            sp->delLeaf(v);
                            
                            // remove "vsib"  from the set of vertices  at level "i"
                            // and from the spanning tree.
                            vset[i]->remove(vsib);
                            sp->delLeaf(vsib);
                            
                            // remove "vpar" from the set of vertices at level "i+1"
                            // and from the spanning tree.
                            vset[i+1]->remove(vpar);
                            sp->delLeaf(vpar);
                            
                            // remove "vparvpar"  from the set of  vertices at level
                            // "i+2" and from the spanning tree.
                            vset[i+2]->remove(vparvpar);
                            sp->delLeaf(vparvpar);
                        }
                    }
                    else {
                        // "vparvpar" has other child than "vpar", and we can make
                        // quadrilaterals by using the children of "vparvpar".
                        
                        // get the sibling of "vpar", "vparsib".
                        tCQMSpanningTreeVertex* vparsib = vparvpar->getSibling(vpar);
                        
                        // if   "vparsib"   has  degree   1,   then   we  make   3
                        // quadrilaterals  by  creating  a  Steiner  point  inside
                        // "vparvpar", and by using the triangles corresponding to
                        // "v", "vsib", "vpar" and "vparsib".
                        if (vparsib->getDegree() == 1) {
                            
                            // if "vparsib" is a  quadrilateral then put it into the
                            // output list of quadrilaterals  and make an hexagon by
                            // using "v", "vsib", "vpar" and "vparvpar".
                            if (vparsib->isQuadrilateral()) {
                                
                                // make a quadrilateral
                                makeQuad(vparsib,vs,lquads);
                                
                                // remove "vparsib" from the  set of vertices at level
                                // "i+1" and from the spanning tree.
                                vset[i+1]->remove(vparsib);
                                sp->delLeaf(vparsib);
                                
                                // get sibling of "v"
                                tCQMSpanningTreeVertex* vsib = vpar->getSibling(v);
                                
                                makeConvexQuadFromHex(v,vpar,vsib,vparvpar,vs,lquads);
                                
                                // remove "v"  from the set  of vertices at  level "i"
                                // and from the spanning tree.
                                vset[i]->remove(v);
                                sp->delLeaf(v);
                                
                                // remove "vsib" from the set of vertices at level "i"
                                // and from the spanning tree.
                                vset[i]->remove(vsib);
                                sp->delLeaf(vsib);
                                
                                // remove "vpar" from the set of vertices at level
                                // "i+1" and from the spanning tree.
                                vset[i+1]->remove(vpar);
                                sp->delLeaf(vpar);
                                
                                // remove "vparvpar" from the set of vertices at level
                                // "i+2" and from the spanning tree.
                                vset[i+2]->remove(vparvpar);
                                sp->delLeaf(vparvpar);
                            }
                            else {
                                // if   "vparsib"  is  a   pentagon,  we   can  convex
                                // quadrangulate it and be left with a triangle.
                                
                                if (vparsib->isPentagon()) {
                                    // "vparsib"  is a  pentagon so  that we  can convex
                                    // quadrangulate  it and  be left  with  a triangle.
                                    makeConvexQuadFromPent(vparsib,vs,lquads);
                                }
                                
                                // get the sibling of "v", "vsib"
                                tCQMSpanningTreeVertex* vsib = vpar->getSibling(v);
                                
                                // if it  is possible to  insert one Steiner  point in
                                // "vparvpar" such that this point is visible from any
                                // point in "v", "vsib", "vpar", and "vparsib", we use
                                // it to convex  quadrangulate the septagon defined by
                                // them.  Otherwise,  we can break  this septagon into
                                // three pentagons,  and convex quadrangulate  them by
                                // using at most 2 Steiner points per pentagon.
                                makeQuad(v,vsib,vpar,vparsib,vparvpar,vs,lquads);
                                
                                // remove "v"  from the set  of vertices at  level "i"
                                // and from the spanning tree.
                                vset[i]->remove(v);
                                sp->delLeaf(v);
                                
                                // remove "vsib" from  "vset[i]" and from the spanning
                                // tree.
                                vset[i]->remove(vsib);
                                sp->delLeaf(vsib);
                                
                                // remove  "vpar" and  "vparsib" from  "vset[i+1]" and
                                // from the spanning tree.
                                vset[i+1]->remove(vpar);
                                sp->delLeaf(vpar);
                                
                                vset[i+1]->remove(vparsib);
                                sp->delLeaf(vparsib);
                            }
                        }
                        else {
                            // the  sibling of  "vpar", "vparsib",  has degree  3 as
                            // there is no longer vertex of degree 2 at level "i+1".
                            
                            // get the sibling of  "v", "vsib1", and the children of
                            // "vparsib", "vsib2" and "vsib3".
                            tCQMSpanningTreeVertex* vsib1 = vpar->getSibling(v);
                            tCQMSpanningTreeVertex* vsib2 = vparsib->getChild(0);
                            tCQMSpanningTreeVertex* vsib3 = vparsib->getChild(1);
                            
                            // create a Steiner point  inside "vparvpar", and make 4
                            // quadrilaterals  by using this  Steiner point  and the
                            // triangles  corresponding  to  "v",  "vsib1",  "vpar",
                            // "vparsib", "vsib2", and "vsib3".
                            makeQuad(v,vsib1,vpar,vsib2,vsib3,vparsib,vparvpar,vs,lquads);
                            
                            // remove "v" from the set  of vertices at level "i" and
                            // from the spanning tree.
                            vset[i]->remove(v);
                            sp->delLeaf(v);
                            
                            // remove "vsib1",  "vsib2", and "vsib3"  from "vset[i]"
                            // and from the spanning tree.
                            vset[i]->remove(vsib1);
                            sp->delLeaf(vsib1);
                            vset[i]->remove(vsib2);
                            sp->delLeaf(vsib2);
                            vset[i]->remove(vsib3);
                            sp->delLeaf(vsib3);
                            
                            // remove "vpar" and "vparsib" from "vset[i+1]" and from
                            // the spanning tree.
                            vset[i+1]->remove(vpar);
                            sp->delLeaf(vpar);
                            vset[i+1]->remove(vparsib);
                            sp->delLeaf(vparsib);
                        }
                    }
                }  // end of while to process vertices whose parent has degree
                // 3.
                
                // verify if there is any vertex left inside "vset[i]". If so,
                // there  must be  only  one  vertex.  If  this  vertex has  a
                // parent, its  parent must be the tree  root. Otherwise, this
                // vertex is the tree root.
                if (!vset[i]->empty()) {
                    // get the only vertex in "vset[i]"
                    v = vset[i]->front();
                    
                    // remove it from "vset[i]"
                    vset[i]->pop_front();
                    
                    // verify if "v" has a parent
                    tCQMSpanningTreeVertex* vpar = v->getParent();
                    if (vpar != 0) {
                        // parent of  "v", "vpar",  is the tree  root
                        
                        // if "vpar"  is a quadrilateral,  we convex quadrangulate
                        // the pentagon defined by "v" and "vpar" and one triangle
                        // is left.
                        if (vpar->isQuadrilateral()) {
                            // convex quadrangulate the  pentagon defined by "v" and
                            // "vpar".
                            makeConvexQuadFromPent(v,vpar,vs,lquads);
                            
                            // remove "v" from the spanning tree
                            sp->delLeaf(v);
                        }
                        else {
                            // "v" and "vpar" define a quadrilateral
                            if (!makeQuadIfConvex(v,vpar,vs,lquads)) {
                                try {
                                    makeConvexQuadFrom2Tri(v,vpar,vs,lquads);
                                }
                                catch (const tExceptionObject& xpt) {
                                    treatException(xpt);
                                    exit(0);
                                }
                            }
                            
                            // remove "v" from the spanning tree
                            sp->delLeaf(v);
                            
                            // remove  "vpar" from "vset[i+1]"  and from  the spanning
                            // tree.
                            vset[i+1]->remove(vpar);
                            sp->delLeaf(vpar);
                        }
                    }
                    else {
                        // "v" is the tree root, which means that we have only one
                        // triangle left. In this case, we must modify the contour
                        // enclosing this triangle by inserting a Steiner point in
                        // common  segment of  the  triangle and  the contour.  By
                        // doing so, we make the triangle into a quadrilateral.
                        insertVertexInContour(v,vs,lquads);
                        
                        // remove "v" from the spanning tree
                        sp->delLeaf(v);
                    }
                }
            } // end of  the if-statement that processes "vset[i]"  if it is
            // not empty. By this point, "vset[i]" should be empty.
            
        } // end of  the for-statement that processes all  set of vertices
        // of a given spanning tree.  By this point, we have generated a
        // quadrangulation by converting the triangulation corresponding
        // to this spanning tree.
        
        
    } // end  of the while-statement  that processes  all contours  of a
    // triangular  mesh.    By  this   point,  we  have   generated  a
    // quadrangulation of the interior of each contour.
    
    for(int i = 0; i < 16; i++){
        printf("Case %s: %d\n", PentCaseStr[i], PentCaseCount[i]);
    }
    
    return lquads;
}

