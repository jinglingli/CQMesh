//  ===================================================================
//
//  Program:   Convex Quadrilateral Mesh Generator (CQMesh)
//  Module:    polygon2.h
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

#include "polygon2.h"

// -----------------------------------------------------
// tCQMIndVertex2
// --------------

// constructors
tCQMIndVertex2::tCQMIndVertex2()
{
  m_id = -1;
  m_v = Point2(0,0);
}

tCQMIndVertex2::tCQMIndVertex2(int id, const Point2& p)
{
  m_id = id;
  m_v = p;
}

tCQMIndVertex2::tCQMIndVertex2(const tCQMIndVertex2& v)
{
  m_id = v.m_id;
  m_v = v.m_v;
}

// get methods
int
tCQMIndVertex2::getVertexId() const
{
  return m_id;
}

Point2
tCQMIndVertex2::getPoint() const
{
  return m_v;
}

// set methods
void
tCQMIndVertex2::setVertexId(int id)
{
  m_id = id;
}

void
tCQMIndVertex2::setPoint(const Point2& p)
{
  m_v = p;
}

// -----------------------------------------------------
// tCQMPolygon2
// ------------

// constructors
tCQMPolygon2::tCQMPolygon2()
{
  m_size = 0;
}

tCQMPolygon2::tCQMPolygon2(int n)
{
  m_size = 0;
  m_verts.resize(n);
  m_ec.resize(n);
}

tCQMPolygon2::tCQMPolygon2(const tCQMPolygon2& q)
{
  m_size = q.m_size;
  m_verts.resize(m_size);
  m_ec.resize(m_size);

  for (int i=0; i<m_size; i++) {
    m_verts[i] = q.m_verts[i];
    m_ec[i] = q.m_ec[i];
  }
}


// -------------------------------------------------------------------
// Method: findReflexPoint()
// 
// Returns the index of a reflex vertex in the polygon, if any.
//
// -------------------------------------------------------------------
int
tCQMPolygon2::findReflexPoint() const throw (tExceptionObject)
{
   for (int i=0; i<m_size; i++) if (isReflex(i)) return i;

   std::stringstream ss(std::stringstream::in | 
			std::stringstream::out);
   ss << "Polygon has no reflex vertex";
   throw tExceptionObject(__FILE__,__LINE__,ss.str());
}

// -------------------------------------------------------------------
// Method: getIndexOfPoint()
// 
// Returns  the  index  of a  vertex  in  this  polygon with  a  given
// location.   If such there  is no  such vertex  in this  polygon, an
// exception is thrown.
//
// -------------------------------------------------------------------
int
tCQMPolygon2::getIndexOfPoint(const Point2& p) const
  throw (tExceptionObject)
{
  for (int i=0; i<m_size; i++) {
    if (p == m_verts[i].getPoint()) return i;
  }

  std::stringstream ss(std::stringstream::in | 
		       std::stringstream::out);
  ss << "Polygon has no point at the given location";
  throw tExceptionObject(__FILE__,__LINE__,ss.str());
}

int
tCQMPolygon2::getNumOfReflexPoints() const
{
  int counter = 0;
  for (int i=0; i<m_size; i++) {
    if (this->isReflex(i)) ++counter;
  }

  return counter;
}

double
tCQMPolygon2::area() const
{
  Polygon2 q;

  for (int i=0; i<m_size; i++) {
    q.push_back(m_verts[i].getPoint());
  }

  return q.area();
}

// query methods
bool
tCQMPolygon2::isCW() const
{
  return !this->isCCW();
}

bool
tCQMPolygon2::isCCW() const
{
  Polygon2 q;

  for (int i=0; i<m_size; i++) {
    q.push_back(m_verts[i].getPoint());
  }

  return q.is_counterclockwise_oriented();
}

bool
tCQMPolygon2::isConvex() const
{
  for (int i=0; i<m_size; i++) {
    if (isReflex(i)) return false;
  }

  return true;
}

bool
tCQMPolygon2::isReflex(int i) const
{
  CGAL::Oriented_side right = (this->isCW()) ? CGAL::ON_POSITIVE_SIDE :
    CGAL::ON_NEGATIVE_SIDE;

  Point2 p1 = m_verts[i].getPoint();
  Point2 p2 = m_verts[(i+1) % m_size].getPoint();
  Point2 p3 = m_verts[(i+m_size-1) % m_size].getPoint();
  Line2 l1(p1,p2);
  Line2 l2(p3,p1);

  // there  is a  redundancy  here, but  we  do so  to avoid  annoying
  // numerical problems.
  return ((l1.oriented_side(p3) == right) && 
	  (l2.oriented_side(p2) == right));
}

bool
tCQMPolygon2::isSimple() const
{
  Polygon2 q;

  for (int i=0; i<m_size; i++) {
    q.push_back(m_verts[i].getPoint());
  }

  return q.is_simple();
}

// modifiers
void
tCQMPolygon2::reverse()
{
  tCQMIndVertex2 aux1;
  bool aux2;
  int size2 = m_size >> 1;

  for (int i=0; i<size2; i++) {
    aux1  = m_verts[i];
    m_verts[i] = m_verts[m_size-i-1];
    m_verts[m_size-i-1] = aux1;

    aux2  = m_ec[i];
    m_ec[i] = m_ec[m_size-i-1];
    m_ec[m_size-i-1] = aux2;
  }
}

void 
tCQMPolygon2::insert(int i, int id, const Point2& p)
{
  m_verts[i].setVertexId(id);
  m_verts[i].setPoint(p);
  ++m_size;
}

void 
tCQMPolygon2::insert(int i, const tCQMIndVertex2& v)
{
  m_verts[i] = v;
  ++m_size;
}

void 
tCQMPolygon2::shiftLeft(int d)
{
  for (int i=0; i<d; i++) {
    tCQMIndVertex2 aux1 = m_verts[0];
    bool aux2 = m_ec[0];
    for (int j=0; j<m_size-1; j++) {
      m_verts[j] = m_verts[j+1];
      m_ec[j] = m_ec[j+1];
    }
    m_verts[m_size-1] = aux1;
    m_ec[m_size-1] = aux2;
  }
}

void 
tCQMPolygon2::shiftRight(int d)
{
  for (int i=0; i<d; i++) {
    tCQMIndVertex2 aux1 = m_verts[m_size-1];
    bool aux2 = m_ec[m_size-1];
    for (int j=m_size-1; j>0; j--) {
      m_verts[j] = m_verts[j-1];
      m_ec[j] = m_ec[j-1];
    }
    m_verts[0] = aux1;
    m_ec[0] = aux2;
  }
}
