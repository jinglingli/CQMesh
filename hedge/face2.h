//  ===================================================================
//
//  Program:   Convex Quadrilateral Mesh Generator (CQMesh)
//  Module:    face2.h
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

#ifndef __face2_h
#define __face2_h

#include "halfedge2.h"
#include "edge2.h"
#include "face2attrib.h"

// -----------------------------------------------------
// forward definition of class
// -----------------------------------------------------
class tCQMDualGraphVertex;
class tCQMComplex2;

// -----------------------------------------------------
// tCQMFace2
// ---------
class tCQMFace2 {
   private:
     // private attributes
     tCQMHalfEdge2* m_he;         // one of the halfedges on this face cycle
     tCQMComplex2* m_cp;          // complex component containing this face
     tCQMFace2Attrib* m_att;      // face attributes
     tCQMDualGraphVertex* m_dual; // pointer to vertex in the dual graph of the
                                  // mesh containing this face 

     // private methods
     // 
     // set methods
     // only methods in tCQMComplex2 will use them
     void setHalfEdge(tCQMHalfEdge2*);
     void setComplex(tCQMComplex2*);
     void setDualVertex(tCQMDualGraphVertex*);

   public:
     // constructors
     tCQMFace2();
     tCQMFace2(tCQMHalfEdge2*, tCQMComplex2*);

     // get methods
     tCQMHalfEdge2* getHalfEdge() const;
     tCQMComplex2* getComplex() const;
     tCQMFace2Attrib* getAttribute();
     tCQMEdge2* getCommonEdge(tCQMFace2*);
     tCQMDualGraphVertex* getDualVertex();
     int getNumVerts() const;
     Point2 getBarycenter() const;

     // set method
     void setAttribute(tCQMFace2Attrib*);

     // query method
     bool isAdjToBoundary() const;
     bool isAdjToConstrainedEdge() const;

     friend class tCQMComplex2;
};

// inline implementation of methods

// public method implementation
// get methods
inline tCQMHalfEdge2*
tCQMFace2::getHalfEdge() const
{
  return m_he;
}

inline tCQMComplex2*
tCQMFace2::getComplex() const
{
  return m_cp;
}

inline tCQMFace2Attrib*
tCQMFace2::getAttribute()
{
  return m_att;
}

inline tCQMDualGraphVertex*
tCQMFace2::getDualVertex()
{
  return m_dual;
}

// private method implementation
// set methods
inline void
tCQMFace2::setHalfEdge(tCQMHalfEdge2* he)
{
  m_he = he;
}

inline void
tCQMFace2::setComplex(tCQMComplex2* cp)
{
  m_cp = cp;
}

inline void
tCQMFace2::setAttribute(tCQMFace2Attrib* att)
{
  m_att = att;
}

inline void
tCQMFace2::setDualVertex(tCQMDualGraphVertex* dual)
{
  m_dual = dual;
}

#endif   // __face2_h
