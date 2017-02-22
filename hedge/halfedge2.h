//  ===================================================================
//
//  Program:   Convex Quadrilateral Mesh Generator (CQMesh)
//  Module:    halfedge2.h
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

#ifndef __halfedge2_h
#define __halfedge2_h

#include "vertex2.h"

// -----------------------------------------------------
// forward definition of classes
// -----------------------------------------------------
class tCQMComplex2;
class tCQMFace2;
class tCQMEdge2;


// -----------------------------------------------------
// tCQMHalfEdge2
// -------------
class tCQMHalfEdge2 {
   private:
     // private attributes
     tCQMVertex2* m_vert;    // intial vertex of this halfedge
     tCQMEdge2* m_edge;      // edge where it belongs to
     tCQMFace2* m_face;      // face where it belongs to
     tCQMHalfEdge2* m_next;  // next halfedge on the face cycle
     tCQMHalfEdge2* m_prev;  // previous halfedge on the face cycle

     // private methods
     // 
     // set methods
     // only methods in tCQMComplex2 will use them
     void setVertex(tCQMVertex2*);
     void setEdge(tCQMEdge2*);
     void setFace(tCQMFace2*);
     void setNext(tCQMHalfEdge2*);
     void setPrev(tCQMHalfEdge2*);

   public:
     // constructors
     tCQMHalfEdge2();
     tCQMHalfEdge2(tCQMVertex2*, tCQMEdge2*, tCQMFace2*, tCQMHalfEdge2*, tCQMHalfEdge2*);
    	
     // get methods
     tCQMVertex2* getVertex() const;
     tCQMEdge2* getEdge() const;
     tCQMFace2* getFace() const;
     tCQMHalfEdge2* getNext() const;
     tCQMHalfEdge2* getPrev() const;

     friend class tCQMComplex2;
};

// inline implementation of methods

// public method implementation
// get methods
inline tCQMVertex2*
tCQMHalfEdge2::getVertex() const
{
  return m_vert;
}

inline tCQMEdge2*
tCQMHalfEdge2::getEdge() const
{
  return m_edge;
}

inline tCQMFace2*
tCQMHalfEdge2::getFace() const
{
  return m_face;
}

inline tCQMHalfEdge2*
tCQMHalfEdge2::getNext() const
{
  return m_next;
}

inline tCQMHalfEdge2*
tCQMHalfEdge2::getPrev() const
{
  return m_prev;
}

// private method implementation
// set methods
inline void
tCQMHalfEdge2::setVertex(tCQMVertex2* v)
{
  m_vert = v;
}

inline void
tCQMHalfEdge2::setEdge(tCQMEdge2* e)
{
  m_edge = e;
}

inline void
tCQMHalfEdge2::setFace(tCQMFace2* f)
{
  m_face = f;
}

inline void
tCQMHalfEdge2::setNext(tCQMHalfEdge2* he)
{
  m_next = he;
}

inline void
tCQMHalfEdge2::setPrev(tCQMHalfEdge2* he)
{
  m_prev = he;
}

#endif   // __halfedge2_h
