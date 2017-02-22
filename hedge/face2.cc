//  ===================================================================
//
//  Program:   Convex Quadrilateral Mesh Generator (CQMesh)
//  Module:    face2.cc
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

#include "face2.h"

// -----------------------------------------------------
// tCQMFace2
// ---------

// constructors
tCQMFace2::tCQMFace2():
  m_he(0),
  m_cp(0),
  m_att(0),
  m_dual(0)
{}

tCQMFace2::tCQMFace2(tCQMHalfEdge2* he, tCQMComplex2* cp):
  m_he(he),
  m_cp(cp),
  m_att(0),
  m_dual(0)
{}

// get method
int
tCQMFace2::getNumVerts() const
{
  tCQMHalfEdge2* he = m_he->getNext();
  int count = 1;
  while (he != m_he) {
    ++count;
    he = he->getNext();
  }

  return count;
}

tCQMEdge2*
tCQMFace2::getCommonEdge(tCQMFace2* f)
{
  tCQMHalfEdge2 *he1 = m_he;
  tCQMEdge2* e = he1->getEdge();
  do {
    tCQMHalfEdge2 *he2 = e->getMate(he1);
    if ((he2 != 0) && (he2->getFace() == f)) {
      return e;
    }
    else {
      he1 = he1->getNext();
      e = he1->getEdge();
    }
  } while (he1 != m_he);

  return 0;
}

Point2
tCQMFace2::getBarycenter() const
{
  tCQMHalfEdge2* he = m_he->getNext();
  double xc = m_he->getVertex()->getPoint().x();
  double yc = m_he->getVertex()->getPoint().y();
  int count = 1;

  while (he != m_he) {
    xc += he->getVertex()->getPoint().x();
    yc += he->getVertex()->getPoint().y();
    ++count;
    he = he->getNext();
  }

  xc /= double(count);
  yc /= double(count);

  return Point2(xc,yc);
}

// query method
bool
tCQMFace2::isAdjToBoundary() const
{
  if (m_he->getEdge()->isOnBoundary()) return true;
  else {
    tCQMHalfEdge2* he = m_he->getNext();
    while (he != m_he) {
      if (he->getEdge()->isOnBoundary()) return true;
      else he = he->getNext();
    }
  }

  return false;
}

bool
tCQMFace2::isAdjToConstrainedEdge() const
{
  if (m_he->getEdge()->isConstrained()) return true;
  else {
    tCQMHalfEdge2* he = m_he->getNext();
    while (he != m_he) {
      if (he->getEdge()->isConstrained()) return true;
      else he = he->getNext();
    }
  }

  return false;
}
