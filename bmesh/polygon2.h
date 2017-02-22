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

#ifndef __polygon2_h
#define __polygon2_h

#include <sstream>
#include <vector>

#include "typedefs.h"
#include "err.h"

// -----------------------------------------------------
// tCQMIndVertex2
// --------------
class tCQMIndVertex2 {
   private:
     int m_id;
     Point2 m_v;

   public:
     // constructor
     tCQMIndVertex2();
     tCQMIndVertex2(int, const Point2&);
     tCQMIndVertex2(const tCQMIndVertex2&);

     // get methods
     int getVertexId() const;
     Point2 getPoint() const;

     // set methods
     void setVertexId(int);
     void setPoint(const Point2&);
};

// -----------------------------------------------------
// tCQMPolygon2
// ------------
class tCQMPolygon2 {
   protected:
     std::vector<tCQMIndVertex2> m_verts;
     std::vector<bool> m_ec;
     int m_size;

   public:
     // constructors
     tCQMPolygon2();
     tCQMPolygon2(int);
     tCQMPolygon2(const tCQMPolygon2&);

     virtual ~tCQMPolygon2() {};

     // get methods
     int getNumVerts() const;
     int getVertexId(int) const;
     Point2 getPoint(int) const;
     tCQMIndVertex2 getIndVertex(int) const;
     int findReflexPoint() const throw (tExceptionObject);
     int getIndexOfPoint(const Point2&) const 
       throw (tExceptionObject) ;
     int getNumOfReflexPoints() const;
     virtual double area() const;

     // set methods
     void setConstraint(int,bool);

     // query methods
     virtual bool isCW() const;
     virtual bool isCCW() const;
     virtual bool isConvex() const;
     virtual bool isReflex(int) const;
     virtual bool isConstrained(int) const;
     virtual bool isSimple() const;

     // modifiers
     virtual void reverse();
     virtual void insert(int, int, const Point2&);
     virtual void insert(int, const tCQMIndVertex2&);
     virtual void shiftLeft(int);
     virtual void shiftRight(int);
};

// -----------------------------------------------------
// tCQMPolygon2 - inline implementation
// ------------
inline int
tCQMPolygon2::getNumVerts() const
{
  return m_size;
}

inline int
tCQMPolygon2::getVertexId(int i) const
{
  return m_verts[i].getVertexId();
}

inline Point2
tCQMPolygon2::getPoint(int i) const
{
  return m_verts[i].getPoint();
}

inline tCQMIndVertex2
tCQMPolygon2::getIndVertex(int i) const
{
  return m_verts[i];
}

inline void
tCQMPolygon2::setConstraint(int i,bool e)
{
  m_ec[i] = e;
}

inline bool
tCQMPolygon2::isConstrained(int i) const
{
  return m_ec[i];
}

#endif   // __polygon2_h
