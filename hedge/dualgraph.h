//  ===================================================================
//
//  Program:   Convex Quadrilateral Mesh Generator (CQMesh)
//  Module:    dualgraph.h
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

#ifndef _dualgraph_h
#define _dualgraph_h

#include <list>
#include <vector>

#include "edge2.h"
#include "face2.h"
#include "dgvattrib.h"


using namespace std;


const int DEFAULT_GRAPH_DEGREE = 4;


// -----------------------------------------------------
// forward definition of class
// -----------------------------------------------------
class tCQMComplex2;
class tCQMDualGraph;

// -----------------------------------------------------
// tCQMDualGraphVertex
// -------------------
class tCQMDualGraphVertex {
  private:
    // pointer to the dual graph this vertex belongs to.
    tCQMDualGraph* m_dual;

    // pointer to the dual face of this vertex in the mesh.
    tCQMFace2* m_face;

    // array of pointers to adjacent vertices.
    vector<tCQMDualGraphVertex*> m_neigh; 

    // maximum number of neighbors.
    int m_maxdg; 

    // current vertex degree.                       
    int m_degree;

    // flag used by graph algorithms.
    bool m_marked;

    // pointer to the  vertex of the region graph  this vertex belongs
    // to.
    tCQMDualGraphVertexAttrib* m_dgvattrib;

    // private methods
    void addNeighbor(tCQMDualGraphVertex*);
    void delNeighbor(tCQMDualGraphVertex*);
    void setFace(tCQMFace2*);

 public:
    // constructors
    tCQMDualGraphVertex(tCQMDualGraph*, tCQMFace2*);

    // destructor
    ~tCQMDualGraphVertex() {}

    // get methods
    int getDegree() const;
    int getMaxDegree() const;
    tCQMDualGraphVertex* getNeighbor(int) const;
    tCQMFace2* getFace() const;
    tCQMDualGraph* getDualGraph() const;
    tCQMDualGraphVertexAttrib* getDualGraphVertexAttrib() const;

    // query method
    bool isMarked() const;
    bool isDualofConstrainedEdge(int);

    // set method
    void setMark(bool);
    void setDualGraphVertexAttrib(tCQMDualGraphVertexAttrib*);

    friend class tCQMDualGraph;
};

// -----------------------------------------------------
// tCQMDualGraph
// -------------
class tCQMDualGraph {
  private:
    list<tCQMDualGraphVertex*> m_lverts;    // list of vertices of this graph
    tCQMComplex2* m_mesh;

    // modifiers
    void addVert(tCQMDualGraphVertex*);
    void delVert(tCQMDualGraphVertex*);
    void addEdge(tCQMDualGraphVertex*, tCQMDualGraphVertex*);
    void delEdge(tCQMDualGraphVertex*, tCQMDualGraphVertex*);

  public:

    // constructors
    tCQMDualGraph(tCQMComplex2*);
    
    // destructor
    ~tCQMDualGraph();

    // get methods
    int getNumVerts();
    list<tCQMDualGraphVertex*>* getVertices();

    // modifier
    void clearMark();

    friend class tCQMComplex2;
};

// -----------------------------------------------------
// tCQMDualGraphVertex - inline implementation
// -------------------

//
// public methods
//

// get methods
inline int
tCQMDualGraphVertex::getDegree() const
{
  return m_degree;
}

inline int
tCQMDualGraphVertex::getMaxDegree() const
{
  return m_maxdg;
}

inline tCQMDualGraphVertex*
tCQMDualGraphVertex::getNeighbor(int i) const
{
  return m_neigh[i];
}

inline tCQMFace2*
tCQMDualGraphVertex::getFace() const
{
  return m_face;
}

inline tCQMDualGraph*
tCQMDualGraphVertex::getDualGraph() const
{
  return m_dual;
}

inline tCQMDualGraphVertexAttrib*
tCQMDualGraphVertex::getDualGraphVertexAttrib() const
{
  return m_dgvattrib;
}

// query method
inline bool
tCQMDualGraphVertex::isMarked() const
{
  return m_marked;
}

inline bool
tCQMDualGraphVertex::isDualofConstrainedEdge(int i)
{
  return m_face->getCommonEdge(m_neigh[i]->getFace())->isConstrained();
}

// set method
inline void
tCQMDualGraphVertex::setMark(bool mark)
{
  m_marked = mark;
}

inline void
tCQMDualGraphVertex::setDualGraphVertexAttrib(tCQMDualGraphVertexAttrib* dgva)
{
  m_dgvattrib = dgva;
}

//
// private method
//

// set method
inline void
tCQMDualGraphVertex::setFace(tCQMFace2* f)
{
  m_face = f;
}

// -----------------------------------------------------
// tCQMDualGraph - inline implementation
// -------------

//
// public methods
//

// get methods
inline int
tCQMDualGraph::getNumVerts()
{
  return m_lverts.size();
}

inline list<tCQMDualGraphVertex*>*
 tCQMDualGraph::getVertices()
{
  return &m_lverts;
}

//
// private method
//

// modifiers
inline void
tCQMDualGraph::addVert(tCQMDualGraphVertex* v)
{
  m_lverts.push_back(v);
}

inline void
tCQMDualGraph::addEdge(tCQMDualGraphVertex* v1, tCQMDualGraphVertex* v2)
{
  v1->addNeighbor(v2);
  v2->addNeighbor(v1);
}

inline void
tCQMDualGraph::delEdge(tCQMDualGraphVertex* v1, tCQMDualGraphVertex* v2)
{
  v1->delNeighbor(v2);
  v2->delNeighbor(v1);
}

#endif   // __dualgraph_h










