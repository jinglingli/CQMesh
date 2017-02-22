//  ===================================================================
//
//  Program:   Convex Quadrilateral Mesh Generator (CQMesh)
//  Module:    tri2.cc
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

#include "tri2.h"

// -----------------------------------------------------
// tCQMTriangle2
// -------------

// constructors
tCQMTriangle2::tCQMTriangle2() : tCQMPolygon2(3)
{}

tCQMTriangle2::tCQMTriangle2(const tCQMTriangle2& t) :
  tCQMPolygon2(3)
{
  m_verts[0] = t.m_verts[0];
  m_verts[1] = t.m_verts[1];
  m_verts[2] = t.m_verts[2];
}

// query method

bool
tCQMTriangle2::isCCW() const
{
  Point2 p0 = m_verts[0].getPoint();
  Point2 p1 = m_verts[1].getPoint();
  Point2 p2 = m_verts[2].getPoint();

  double a1 = p1.x() - p0.x();
  double a2 = p1.y() - p0.y();
  double b1 = p2.x() - p0.x();
  double b2 = p2.y() - p0.y();

  return (a1*b2 - a2*b1) > 0;
}

// modifier
void
tCQMTriangle2::reverse()
{
  tCQMIndVertex2 v = m_verts[0];
  m_verts[0] = m_verts[2];
  m_verts[2] = v;

  bool e = m_ec[0];
  m_ec[0] = m_ec[2];
  m_ec[2] = e;
}
