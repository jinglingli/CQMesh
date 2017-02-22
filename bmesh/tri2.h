//  ===================================================================
//
//  Program:   Convex Quadrilateral Mesh Generator (CQMesh)
//  Module:    tri2.h
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

#ifndef __tri2_h
#define __tri2_h

#include "typedefs.h"
#include "polygon2.h"

// -----------------------------------------------------
// tCQMTriangle2
// -------------
class tCQMTriangle2 : public tCQMPolygon2 {
   public:
     // constructor
     tCQMTriangle2();
     tCQMTriangle2(const tCQMTriangle2&);

     ~tCQMTriangle2()
     {}

     // query methods
     bool isCCW() const;
     
     // modifiers
     void reverse();
};

#endif   // __tri2_h
