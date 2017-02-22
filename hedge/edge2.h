//  ===================================================================
//
//  Program:   Convex Quadrilateral Mesh Generator (CQMesh)
//  Module:    edge2.h
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

#ifndef __edge2_h
#define __edge2_h

#include "vertex2.h"
#include "halfedge2.h"
#include "edge2attrib.h"

// -----------------------------------------------------
// forward definition of classes
// -----------------------------------------------------
class tCQMComplex2;
class tCQMBoundary2;

// -----------------------------------------------------
// tCQMEdge2
// -----------
class tCQMEdge2 {
   private:
     // private attributes
     tCQMHalfEdge2* m_he1;     // first halfedge of this edge
     tCQMHalfEdge2* m_he2;     // second halfedge of this edge
     tCQMBoundary2* m_bd;      // boundary curve to which this edge belongs
     bool m_constrained;       // flag to indicate if edge is or not constrained
     tCQMEdge2Attrib* m_att;   // edge attributes

     // private methods
     // 
     // set methods
     // only methods in tCQMComplex2 will use them
     void setHalfEdge1(tCQMHalfEdge2*);
     void setHalfEdge2(tCQMHalfEdge2*);
     void setBoundary(tCQMBoundary2*);

   public:
     // constructors
     tCQMEdge2();
     tCQMEdge2(tCQMHalfEdge2*, tCQMHalfEdge2*, tCQMBoundary2*);

     // get methods
     tCQMVertex2* getVertex1() const;
     tCQMVertex2* getVertex2() const;
     tCQMHalfEdge2* getHalfEdge() const;
     tCQMHalfEdge2* getMate(tCQMHalfEdge2*) const;
     tCQMBoundary2* getBoundary() const;
     tCQMEdge2Attrib* getAttribute();
     tCQMVertex2* getCommonVertex(tCQMEdge2*) const;

     // set method
     void setAttribute(tCQMEdge2Attrib*);
     void constraintOn();
     void constraintOff();

     // query methods
     bool isOnBoundary() const;
     bool isConstrained() const;

     friend class tCQMComplex2;
};

// inline implementation of methods

// public method implementation
// get methods
inline tCQMVertex2*
tCQMEdge2::getVertex1() const
{
  return m_he1->getVertex();
}

inline tCQMVertex2*
tCQMEdge2::getVertex2() const
{
  if (m_he2 != 0) return m_he2->getVertex();
  return 0;
}

inline tCQMHalfEdge2*
tCQMEdge2::getHalfEdge() const
{
  return m_he1;
}

inline tCQMHalfEdge2*
tCQMEdge2::getMate(tCQMHalfEdge2* he) const
{
  if (he == m_he1) return m_he2;

  return m_he1;
}

inline tCQMBoundary2*
tCQMEdge2::getBoundary() const
{
  return m_bd;
}

inline tCQMEdge2Attrib*
tCQMEdge2::getAttribute()
{
  return m_att;
}

// query method

inline bool
tCQMEdge2::isOnBoundary() const
{
  return (m_bd != 0);
}

inline bool
tCQMEdge2::isConstrained() const
{
  return m_constrained;
}

// set methods

inline void
tCQMEdge2::constraintOn()
{
  m_constrained = true;
}

inline void
tCQMEdge2::constraintOff()
{
  m_constrained = false;
}

// private method implementation
// set methods
inline void
tCQMEdge2::setHalfEdge1(tCQMHalfEdge2* he)
{
  m_he1 = he;
}

inline void
tCQMEdge2::setHalfEdge2(tCQMHalfEdge2* he)
{
  m_he2 = he;
}

inline void
tCQMEdge2::setBoundary(tCQMBoundary2* bd)
{
  m_bd = bd;
}

inline void
tCQMEdge2::setAttribute(tCQMEdge2Attrib* att)
{
  m_att = att;
}

#endif   // __edge2_h
