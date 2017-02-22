//  ===================================================================
//
//  Program:   Convex Quadrilateral Mesh Generator (CQMesh)
//  Module:    vertex2.cc
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

#include "vertex2.h"

// -----------------------------------------------------
// tCQMVertex2
// -----------

// constructors
tCQMVertex2::tCQMVertex2()
{
  m_shedge = 0;
  m_bhedge = 0;
  m_att = 0;
}

tCQMVertex2::tCQMVertex2(const Point2& p)
{
  m_point = p;
  m_shedge = 0;
  m_bhedge = 0;
  m_att = 0;
}

tCQMVertex2::tCQMVertex2(
			 const Point2& p,
			 tCQMHalfEdge2* she,
			 tCQMHalfEdge2* bhe
			 )
{
  m_point = p;
  m_shedge = she;
  m_bhedge = bhe;
  m_att = 0;
}
