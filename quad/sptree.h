//  ===================================================================
//
//  Program:   Convex Quadrilateral Mesh Generator (CQMesh)
//  Module:    sptree.h
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

#ifndef _sptree_h
#define _sptree_h

#include <list>
#include <vector>
#include <iostream>
#include <fstream>

#include "typedefs.h"
#include "polygon2.h"
#include "halfedge2.h"
#include "dualgraph.h"

const int DEFAULT_TREE_DEGREE = 4;

// -------------------------------------------------------------------
// forward definition of class
// -------------------------------------------------------------------
class tCQMSpanningTree;

// -------------------------------------------------------------------
// tCQMSpanningTreeVert
// --------------------
class tCQMSpanningTreeVertex {
private:
  // pointer to the corresponding vertex in the dual graph
  tCQMDualGraphVertex* m_vert;
  
  // pointer to the parent of this vertex
  tCQMSpanningTreeVertex* m_parent;
  
  // array of pointers to adjacent vertices
  std::vector<tCQMSpanningTreeVertex*> m_child;
  
  // maximum number of neighbor vertices
  int m_maxdg;
  
  // number of children in this vertex
  int m_nchildren;
  
  // flag to indicate if this vertex store a Steiner point
  bool m_haspoint;
  
  // Steiner point. This  point replaces one of the  vertices of the
  // face  corresponding to  the graph  vertex associated  with this
  // tree vertex.
  tCQMIndVertex2 m_point1;
  
  // point following Steiner point in counterclockwise direction
  Point2 m_nextpt;
  
  // flag to indicate if this vertex corresponds to a quadrilateral
  bool m_isquad;
  
  // flag to indicate if this vertex corresponds to a pentagon
  bool m_ispent;
  
  // Steiner  point   that  makes   this  vertex  correspond   to  a
  // quadrilateral
  tCQMIndVertex2 m_point2;
  
  // Steiner point that makes this vertex correspond to a pentagon
  tCQMIndVertex2 m_point3;
  
  // pointer  to  the  half-edge  of  the edge  that  was  split  by
  // inserting the Steiner point 2.
  tCQMHalfEdge2* m_he1;
  
  // pointer  to  the  half-edge  of  the edge  that  was  split  by
  // inserting the Steiner point 2.
  tCQMHalfEdge2* m_he2;
  
  // modifiers
  void setParent(tCQMSpanningTreeVertex*);
  void addChild(tCQMSpanningTreeVertex*);
  void dettach();
  
public:
  // constructor
  tCQMSpanningTreeVertex(tCQMDualGraphVertex*);
  
  // destructor
  ~tCQMSpanningTreeVertex();
  
  // get methods
  int getDegree() const;
  int getMaxDegree() const;
  int getNumChildren() const;
  tCQMSpanningTreeVertex* getParent() const;
  tCQMSpanningTreeVertex* getChild(int) const;
  tCQMSpanningTreeVertex* getSibling(tCQMSpanningTreeVertex*) const;
  tCQMDualGraphVertex* getVertex() const;
  tCQMIndVertex2 getSteinerPoint() const;
  Point2 getSteinerPointSuc() const;
  tCQMIndVertex2 getQuadSteinerPoint() const;
  tCQMIndVertex2 getPentSteinerPoint() const;
  tCQMHalfEdge2* getQuadHalfEdge() const;
  tCQMHalfEdge2* getPentHalfEdge() const;
  
  // query method
  bool isInSubtree(tCQMSpanningTreeVertex*) const;
  bool hasSteinerPoint() const;
  bool isQuadrilateral() const;
  bool isPentagon() const;
  
  // modifiers
  void setSteinerPoint(const tCQMIndVertex2&, tCQMHalfEdge2*);
  void setSteinerPoint(const tCQMIndVertex2&, const Point2&);
  void setQuadSteinerPoint(const tCQMIndVertex2&, tCQMHalfEdge2*);
  void setPentSteinerPoint(const tCQMIndVertex2&, tCQMHalfEdge2*);
  void setQuadPentOff();
  
  // print
  void dump(std::ofstream& myfile, int type);
  friend class tCQMSpanningTree;
};

// -------------------------------------------------------------------
// tCQMSpanningTree
// ----------------
class tCQMSpanningTree {
private:
  tCQMSpanningTreeVertex* m_root;  // root
  
public:
  // constructors
  tCQMSpanningTree(tCQMSpanningTreeVertex*);
  
  // destructor
  ~tCQMSpanningTree();
  
  // get methods
  std::list<std::list<tCQMSpanningTreeVertex*>*>* getVertPerLevel();
  
  // modifiers
  tCQMSpanningTreeVertex* addVert(tCQMSpanningTreeVertex*,
                                  tCQMDualGraphVertex*);
  void delLeaf(tCQMSpanningTreeVertex*);
  
  // query
  bool isInTree(tCQMSpanningTreeVertex*) const;
  
  // dump sptree
  void dump(string fname);
};

// -------------------------------------------------------------------
// tCQMSpanningTreeVertex - inline implementation
// ----------------------

//
// public methods
//

// get methods
inline int
tCQMSpanningTreeVertex::getDegree() const
{
  return this->getNumChildren() + ((this->getParent() != 0) ? 1 : 0);
}

inline int
tCQMSpanningTreeVertex::getMaxDegree() const
{
  return m_maxdg;
}

inline int
tCQMSpanningTreeVertex::getNumChildren() const
{
  return m_nchildren;
}

inline tCQMSpanningTreeVertex*
tCQMSpanningTreeVertex::getParent() const
{
  return m_parent;
}

inline tCQMSpanningTreeVertex*
tCQMSpanningTreeVertex::getChild(int i) const
{
  return m_child[i];
}

inline tCQMDualGraphVertex*
tCQMSpanningTreeVertex::getVertex() const
{
  return m_vert;
}

inline tCQMIndVertex2
tCQMSpanningTreeVertex::getSteinerPoint() const
{
  return m_point1;
}

inline Point2
tCQMSpanningTreeVertex::getSteinerPointSuc() const
{
  return m_nextpt;
}

inline tCQMIndVertex2
tCQMSpanningTreeVertex::getQuadSteinerPoint() const
{
  return m_point2;
}

inline tCQMIndVertex2
tCQMSpanningTreeVertex::getPentSteinerPoint() const
{
  return m_point3;
}

inline tCQMHalfEdge2*
tCQMSpanningTreeVertex::getQuadHalfEdge() const
{
  return m_he1;
}

inline tCQMHalfEdge2*
tCQMSpanningTreeVertex::getPentHalfEdge() const
{
  return m_he2;
}

inline bool
tCQMSpanningTreeVertex::hasSteinerPoint() const
{
  return m_haspoint;
}

inline bool
tCQMSpanningTreeVertex::isPentagon() const
{
  return m_ispent;
}

inline bool
tCQMSpanningTreeVertex::isQuadrilateral() const
{
  return m_isquad;
}

// modifiers

inline void
tCQMSpanningTreeVertex:: setSteinerPoint(
                                         const tCQMIndVertex2& p,
                                         const Point2& nxt
                                         )
{
  m_haspoint = true;
  m_point1 = p;
  m_nextpt = nxt;
}

inline void
tCQMSpanningTreeVertex:: setSteinerPoint(
                                         const tCQMIndVertex2& p,
                                         tCQMHalfEdge2* he
                                         )
{
  setSteinerPoint(p,he->getVertex()->getPoint());
}

inline void
tCQMSpanningTreeVertex:: setQuadSteinerPoint(
                                             const tCQMIndVertex2& p,
                                             tCQMHalfEdge2* he
                                             )
{
  m_isquad = true;
  m_point2 = p;
  m_he1 = he;
}

inline void
tCQMSpanningTreeVertex:: setPentSteinerPoint(
                                             const tCQMIndVertex2& p,
                                             tCQMHalfEdge2* he
                                             )
{
  m_ispent = true;
  m_point3 = p;
  m_he2 = he;
}

inline void
tCQMSpanningTreeVertex::setQuadPentOff()
{
  m_isquad = false;
  m_ispent = false;
}

inline void
tCQMSpanningTreeVertex::setParent(tCQMSpanningTreeVertex* par)
{
  m_parent = par;
}

// -------------------------------------------------------------------
// tCQMSpanningTree - inline implementation
// ----------------

inline tCQMSpanningTreeVertex*
tCQMSpanningTree::addVert(
                          tCQMSpanningTreeVertex* par,
                          tCQMDualGraphVertex* v
                          )
{
  if (isInTree(par)) {
    tCQMSpanningTreeVertex* child = (tCQMSpanningTreeVertex*) 
    new tCQMSpanningTreeVertex(v);
    child->setParent(par);
    par->addChild(child);
    return child;
  }
  
  return 0;
}

// query method
inline bool
tCQMSpanningTree::isInTree(tCQMSpanningTreeVertex* v) const {
  return m_root->isInSubtree(v);
} 

#endif   // __sptree_h
