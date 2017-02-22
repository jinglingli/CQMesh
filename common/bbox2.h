//  ===================================================================
//
//  Program:   Convex Quadrilateral Mesh Generator (CQMesh)
//  Module:    bbox2.h
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

#ifndef __bbox2_h
#define __bbox2_h

#include "typedefs.h"

// -----------------------------------------------------
// tCQMBBox2
// ---------
class tCQMBBox2 {
   private:
     // private attributes 
     Point2 m_pmin;     // lower-left point of this bounding box          
     Point2 m_pmax;     // upper-right point of this bounding box
     
   public:
     // constructors
     tCQMBBox2();
     tCQMBBox2(const Point2&, const Point2&);
     tCQMBBox2(const tCQMBBox2&);
    	
     // get methods
     Point2 getPMin() const;
     Point2 getPMax() const;

     // set methods
     void setPMin(const Point2&);
     void setPMax(const Point2&);
};

// inline implementation of methods

// public method implementation
// get methods
inline Point2
tCQMBBox2::getPMin() const
{
  return m_pmin;
}

inline Point2
tCQMBBox2::getPMax() const
{
  return m_pmax;
}

inline void
tCQMBBox2::setPMin(const Point2& p)
{
  m_pmin = p;
}

inline void
tCQMBBox2::setPMax(const Point2& p)
{
  m_pmax = p;
}

#endif   // __bbox2_h
