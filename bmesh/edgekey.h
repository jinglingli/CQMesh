//  ===================================================================
//
//  Program:   Convex Quadrilateral Mesh Generator (CQMesh)
//  Module:    edgekey.h
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

#ifndef __edgekey_h
#define __edgekey_h

// -----------------------------------------------------
// tCQMEdgeKey
// -----------
class tCQMEdgeKey {
   private:
     int m_E1;
     int m_E2;

   public:
     // constructors
     tCQMEdgeKey();
     tCQMEdgeKey(int, int);
     tCQMEdgeKey(const tCQMEdgeKey&);
  
     // comparator method 
     bool operator<(const tCQMEdgeKey&) const;
  	
     // get methods
     int getEdgeC1() const;
     int getEdgeC2() const;

     // compute inverse edge key
     tCQMEdgeKey reverse() const;
};

#endif   // __edgekey_h
