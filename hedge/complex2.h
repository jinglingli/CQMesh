//  ===================================================================
//
//  Program:   Convex Quadrilateral Mesh Generator (CQMesh)
//  Module:    complex2.h
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

#ifndef __complex2_h
#define __complex2_h

#include <list>

#include "typedefs.h"
#include "bbox2.h"

#include "vertex2.h"
#include "halfedge2.h"
#include "edge2.h"
#include "face2.h"
#include "boundary2.h"
#include "dualgraph.h"

// -----------------------------------------------------
// tCQMComplex2
// ------------
class tCQMComplex2 {
   private:
     list<tCQMVertex2*> m_lverts;       // list of all vertices in this complex
     list<tCQMEdge2*> m_ledges;         // list of all edges in this complex
     list<tCQMFace2*> m_lfaces;         // list of all faces in this complex
     list<tCQMBoundary2*> m_lbounds;    // list of all boundaries in this complex
     tCQMBBox2 m_bb;                    // bounding box enclosing this complex
     tCQMDualGraph* m_dual;             // dual graph of this complex

     // private methods
     void createDualVertex(tCQMFace2*);
     void createDualEdge(tCQMFace2*, tCQMFace2*);
     void deleteDualVertex(tCQMFace2*);

   public:
     // constructor
     tCQMComplex2();

     // destructor
     ~tCQMComplex2();

     // get methods
     list<tCQMVertex2*>* getListVerts();
     list<tCQMEdge2*>* getListEdges();
     list<tCQMFace2*>* getListFaces();
     list<tCQMBoundary2*>* getListBounds();

     int getNumVerts() const;
     int getNumEdges() const;
     int getNumFaces() const;
     int getNumHoles() const;
     int getNumBounds() const;

     tCQMBBox2 getBoundingBox() const;
     tCQMDualGraph* getDualGraph();
     tCQMHalfEdge2* getNextStarVert(tCQMHalfEdge2*);
     tCQMHalfEdge2* getPrevStarVert(tCQMHalfEdge2*);
     tCQMEdge2* getNextBoundEdge(tCQMEdge2*);
     tCQMEdge2* getPrevBoundEdge(tCQMEdge2*);

     // modifiers
     tCQMVertex2* createVertex(const Point2&);
     tCQMFace2* mkFBEV(tCQMVertex2*);
     tCQMEdge2* mkEV(tCQMVertex2*, tCQMHalfEdge2*);
     tCQMEdge2* spEmkEV(tCQMVertex2*, tCQMHalfEdge2*);
     tCQMFace2* mkFE(tCQMHalfEdge2*);
     tCQMFace2* mkFE(tCQMHalfEdge2*, tCQMHalfEdge2*);
     tCQMVertex2* spFmkFE(tCQMHalfEdge2*, const Point2&);
     void rmBEV(tCQMHalfEdge2*, tCQMHalfEdge2*);
     void mkGrmEV(tCQMHalfEdge2*, tCQMHalfEdge2*);
};

// inline implementation of methods

// public method implementation
// get methods
inline list<tCQMVertex2*>*
tCQMComplex2::getListVerts()
{
  return &m_lverts;
}

inline list<tCQMEdge2*>*
tCQMComplex2::getListEdges()
{
  return &m_ledges;
}

inline list<tCQMFace2*>*
tCQMComplex2::getListFaces()
{
  return &m_lfaces;
}

inline list<tCQMBoundary2*>*
tCQMComplex2::getListBounds()
{
  return &m_lbounds;
}

inline int
tCQMComplex2::getNumVerts() const
{ 
  return m_lverts.size(); 
}

inline int
tCQMComplex2::getNumEdges() const
{
  return m_ledges.size();
}

inline int
tCQMComplex2::getNumFaces() const
{
  return m_lfaces.size();
}

inline int
tCQMComplex2::getNumHoles() const
{
  return 1 - (this->getNumVerts()) + (this->getNumEdges()) - (this->getNumFaces());
}

inline int
tCQMComplex2::getNumBounds() const
{
  return m_lbounds.size();
}

inline tCQMBBox2
tCQMComplex2::getBoundingBox() const
{
  return m_bb;
}

inline tCQMDualGraph*
tCQMComplex2::getDualGraph()
{
  return m_dual;
}

#endif   // __complex2_h
