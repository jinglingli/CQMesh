//  ===================================================================
//
//  Program:   Convex Quadrilateral Mesh Generator (CQMesh)
//  Module:    vertex2.h
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

#ifndef __vertex2_h
#define __vertex2_h

#include "typedefs.h"
#include "vertex2attrib.h"

// -----------------------------------------------------
// forward definition of classes
// -----------------------------------------------------
class tCQMComplex2;
class tCQMHalfEdge2;

// -----------------------------------------------------
// tCQMVertex2
// -----------
class tCQMVertex2 {
   private:
     // private attributes 
     Point2 m_point;             // location of the vertex
     tCQMHalfEdge2* m_shedge;    // pointer to some halfedge whose initial
                                 // vertex is this vertex.
     tCQMHalfEdge2* m_bhedge;    // pointer to the boundary half edge whose
                                 // initial vertex this vertex.
                                 // this pointer is null if this vertex is an
                                 // interior one.
     tCQMVertex2Attrib* m_att;   // vertex attributes

     // private methods
     //
     // set methods
     // only methods in tCQMComplex2 will use them 
     void setPoint(const Point2&);
     void setHalfEdge(tCQMHalfEdge2*);
     void setBoundHalfEdge(tCQMHalfEdge2*);
     
   public:
     // constructors
     tCQMVertex2();
     tCQMVertex2(const Point2&);
     tCQMVertex2(const Point2&,tCQMHalfEdge2*,tCQMHalfEdge2*);
    	
     // get methods
     Point2 getPoint() const;
     tCQMHalfEdge2* getHalfEdge() const;
     tCQMHalfEdge2* getBoundHalfEdge() const;
     tCQMVertex2Attrib* getAttribute();

     // set method
     void setAttribute(tCQMVertex2Attrib*);

     // query method
     bool isOnBoundary() const;

     // friend class
     friend class tCQMComplex2;
};

// inline implementation of methods

// public method implementation
// get methods
inline Point2
tCQMVertex2::getPoint() const
{
  return m_point;
}

inline tCQMHalfEdge2*
tCQMVertex2::getHalfEdge() const
{
  return m_shedge;
}

inline tCQMHalfEdge2*
tCQMVertex2::getBoundHalfEdge() const
{
  return m_bhedge;
}

inline tCQMVertex2Attrib*
tCQMVertex2::getAttribute()
{
  return m_att;
}

// private method implementation
// set methods
inline void
tCQMVertex2::setPoint(const Point2& p)
{
  m_point = p;
}

inline void
tCQMVertex2::setHalfEdge(tCQMHalfEdge2* she)
{
  m_shedge = she;
}

inline void
tCQMVertex2::setBoundHalfEdge(tCQMHalfEdge2* bhe)
{
  m_bhedge = bhe;
}

inline void
tCQMVertex2::setAttribute(tCQMVertex2Attrib* att)
{
  m_att = att;
}

// query method

inline bool
tCQMVertex2::isOnBoundary() const
{
  return (m_bhedge != 0);
}

#endif   // __vertex2_h
