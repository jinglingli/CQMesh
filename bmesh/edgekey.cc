//  ===================================================================
//
//  Program:   Convex Quadrilateral Mesh Generator (CQMesh)
//  Module:    edgekey.cc
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

#include "edgekey.h"

// -----------------------------------------------------
// tCQMEdgeKey
// -----------

// constructors
tCQMEdgeKey::tCQMEdgeKey() :
  m_E1(-1),
  m_E2(-1)
{}

tCQMEdgeKey::tCQMEdgeKey(int e1, int e2) :
  m_E1(e1),
  m_E2(e2)
{}

tCQMEdgeKey::tCQMEdgeKey(const tCQMEdgeKey& k) :
  m_E1(k.m_E1),
  m_E2(k.m_E2)
{}
  
// comparator method 
bool
tCQMEdgeKey::operator<(const tCQMEdgeKey& k) const
{
  if (m_E1 == k.m_E1) {
    return m_E2 < k.m_E2;
  }
	
  return m_E1 < k.m_E1;
}
  	
// get methods
int
tCQMEdgeKey::getEdgeC1() const
{
  return m_E1;
}
  		
int
tCQMEdgeKey::getEdgeC2() const
{
  return m_E2;
}

// get methods
tCQMEdgeKey
tCQMEdgeKey::reverse() const
{
  return tCQMEdgeKey(m_E2,m_E1);
}
