//  ===================================================================
//
//  Program:   Convex Quadrilateral Mesh Generator (CQMesh)
//  Module:    quad2.cc
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
#include <iostream>

#include "quad2.h"

// -----------------------------------------------------
// tCQMQuadrilateral2
// ------------------

// constructors
tCQMQuadrilateral2::tCQMQuadrilateral2() : tCQMPolygon2(4)
{}

tCQMQuadrilateral2::tCQMQuadrilateral2(const tCQMQuadrilateral2& q):
  tCQMPolygon2(4)
{
  m_verts[0] = q.m_verts[0];
  m_verts[1] = q.m_verts[1];
  m_verts[2] = q.m_verts[2];
  m_verts[3] = q.m_verts[3];
}

// modifiers
void
tCQMQuadrilateral2::reverse()
{
  tCQMIndVertex2 v;

  v = m_verts[3];
  m_verts[3] = m_verts[0];
  m_verts[0] = v;

  v = m_verts[2];
  m_verts[2] = m_verts[1];
  m_verts[1] = v;

  bool aux = m_ec[0];
  m_ec[0] = m_ec[3];
  m_ec[3] = aux;

  aux = m_ec[1];
  m_ec[1] = m_ec[2];
  m_ec[2] = aux;
}
