//  ===================================================================
//
//  Program:   Convex Quadrilateral Mesh Generator (CQMesh)
//  Module:    halfedge2.cc
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

#include "halfedge2.h"

// -----------------------------------------------------
// tCQMHalfEdge2
// -------------

// constructors
tCQMHalfEdge2::tCQMHalfEdge2()
{
  m_vert = 0;
  m_edge = 0;
  m_face = 0;
  m_next = 0;
  m_prev = 0;
}

tCQMHalfEdge2::tCQMHalfEdge2(
			     tCQMVertex2* v,
			     tCQMEdge2* e,
			     tCQMFace2* f,
			     tCQMHalfEdge2* nxt,
			     tCQMHalfEdge2* prv
			     )
{
  m_vert = v;
  m_edge = e;
  m_face = f;
  m_next = nxt;
  m_prev = prv;
}
