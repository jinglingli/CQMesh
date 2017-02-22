//  ===================================================================
//
//  Program:   Convex Quadrilateral Mesh Generator (CQMesh)
//  Module:    boundary.h
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

#ifndef __boundary2_h
#define __boundary2_h

#include "edge2.h"

// -----------------------------------------------------
// tCQMBoundary2
// -------------
class tCQMBoundary2 {
   private:
     // private attribute
     tCQMEdge2* m_edge;       // one of the boundary edges on this boundary

     // private method
     // 
     // set methods
     // only methods in tCQMComplex2 will use them
     void setBoundaryEdge(tCQMEdge2*);

   public:
     // constructors
     tCQMBoundary2();
     tCQMBoundary2(tCQMEdge2*);

     // get method
     tCQMEdge2* getBoundaryEdge() const;

     friend class tCQMComplex2;
};

// inline implementation of methods

// public method implementation
// get method
inline tCQMEdge2*
tCQMBoundary2::getBoundaryEdge() const
{
  return m_edge;
}

// private method implementation
// set method
inline void
tCQMBoundary2::setBoundaryEdge(tCQMEdge2* e)
{
  m_edge = e;
}

#endif   // __boundary2_h
