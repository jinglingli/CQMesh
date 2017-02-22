//  ===================================================================
//
//  Program:   Convex Quadrilateral Mesh Generator (CQMesh)
//  Module:    bbox2.cc
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

#include "bbox2.h"

// -----------------------------------------------------
// tCQMBBox2
// ---------

// constructors
tCQMBBox2::tCQMBBox2()
{
  m_pmin = Point2(0,0);
  m_pmax = Point2(0,0);
}

tCQMBBox2::tCQMBBox2(const Point2& p1, const Point2& p2)
{
  m_pmin = p1;
  m_pmax = p2;
}

tCQMBBox2::tCQMBBox2(const tCQMBBox2& bb)
{
  m_pmin = bb.m_pmin;
  m_pmax = bb.m_pmax;
}
