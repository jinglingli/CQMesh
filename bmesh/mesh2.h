//  ===================================================================
//
//  Program:   Convex Quadrilateral Mesh Generator (CQMesh)
//  Module:    mesh2.h
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

#ifndef __mesh2_h
#define __mesh2_h

#include <string>
#include <exception>
#include <list>
#include <map>

#include "typedefs.h"
#include "err.h"

#include "edgekey.h"
#include "tri2.h"

#include "iomesh2.h"

#include "vertex2.h"
#include "halfedge2.h"
#include "edge2.h"
#include "face2.h"
#include "boundary2.h"
#include "complex2.h"

typedef std::multimap<int,tCQMVertex2*> VertexSet;
typedef std::multimap<tCQMEdgeKey,tCQMHalfEdge2*> EdgeSet;

// -----------------------------------------------------
// tCQMGlueEdge
// -----------------------------------------------------
class tCQMGlueEdge2 {
   public:
     tCQMHalfEdge2* m_he1;
     tCQMHalfEdge2* m_he2;

     tCQMGlueEdge2() :
       m_he1(0),
       m_he2(0)
     {}
};

// -----------------------------------------------------
// tCQMMesh2Builder
// ----------------
class tCQMMesh2Builder {
   private:
     VertexSet vs;
     EdgeSet es;
     tCQMComplex2* m_cp;

   public:
     // constructors
     tCQMMesh2Builder();

     // destructor
     ~tCQMMesh2Builder();

     // searchers
     bool findVertexInBound(int, const tCQMEdgeKey&) const;

     // query
     bool isConstrainedEdge(const tCQMEdgeKey&, const ConstrainedEdgeSet&) const;

     // modifiers
     void insertGlue(tCQMHalfEdge2*, tCQMHalfEdge2*, list<tCQMGlueEdge2>&);
     void makeTriFace(const tCQMTriangle2&, int, int, int, const Point2&, const Point2&,
		      const Point2&, const tCQMEdgeKey&, const tCQMEdgeKey&, const tCQMEdgeKey&,
		      const ConstrainedEdgeSet&);
     void makeTriFace(const tCQMTriangle2&, int, const Point2&, const tCQMEdgeKey&,
		      const tCQMEdgeKey&, const tCQMEdgeKey&, const ConstrainedEdgeSet&);
     void makeTriFaceGlue1(const tCQMTriangle2&, int, int, int, const Point2&, const Point2&,
			   const Point2&, const tCQMEdgeKey&, const tCQMEdgeKey&, const tCQMEdgeKey&,
			   list<tCQMGlueEdge2>&, const ConstrainedEdgeSet&);
     void makeTriFaceGlue2(const tCQMTriangle2&, int, int, int, const Point2&, const Point2&,
			   const Point2&, const tCQMEdgeKey&, const tCQMEdgeKey&, const tCQMEdgeKey&,
			   int, list<tCQMGlueEdge2>&, const ConstrainedEdgeSet&);
     void makeTriFaceGlue3(const tCQMTriangle2&, int, int, int, const Point2&, const Point2&,
			   const Point2&, const tCQMEdgeKey&, const tCQMEdgeKey&, const tCQMEdgeKey&,
			   list<tCQMGlueEdge2>&, const ConstrainedEdgeSet&);
     void glueFaces(tCQMHalfEdge2*, tCQMHalfEdge2*);

     tCQMComplex2* buildTriMesh(const std::string&);
};

#endif   // __mesh2_h
