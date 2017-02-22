//  ===================================================================
//
//  Program:   Convex Quadrilateral Mesh Generator (CQMesh)
//  Module:    edge2.cc
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

#include "edge2.h"

// -----------------------------------------------------
// tCQMEdge2
// ---------

// constructors
tCQMEdge2::tCQMEdge2()
{
  m_he1 = 0;
  m_he2 = 0;
  m_bd = 0;
  m_att = 0;
  m_constrained = false;
}

tCQMEdge2::tCQMEdge2(
		     tCQMHalfEdge2* he1,
		     tCQMHalfEdge2* he2,
		     tCQMBoundary2* bd
		     )
{
  m_he1 = he1;
  m_he2 = he2;
  m_bd = bd;
  m_att = 0;
  m_constrained = false;
}

// get method
tCQMVertex2*
tCQMEdge2::getCommonVertex(tCQMEdge2* e) const
{
  tCQMHalfEdge2* he1 = this->getHalfEdge();
  
  tCQMVertex2* v1 = he1->getVertex();
  tCQMVertex2* v2 = he1->getNext()->getVertex();

  he1 = e->getHalfEdge();
  
  tCQMVertex2* v3 = he1->getVertex();
  tCQMVertex2* v4 = he1->getNext()->getVertex();

  if ((v1 == v3) || (v1 == v4)) return v1;
  if ((v2 == v3) || (v2 == v4)) return v2;

  return 0;
}
