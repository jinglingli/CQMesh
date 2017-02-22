//  ===================================================================
//
//  Program:   Convex Quadrilateral Mesh Generator (CQMesh)
//  Module:    regiongraph.h
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

#ifndef _regiongraph_h
#define _regiongraph_h

#include <list>
#include <vector>

#include "dgvattrib.h"
#include "dualgraph.h"
#include "halfedge2.h"

using namespace std;

const int DEFAULT_REGION_GRAPH_DEGREE = 4;

// -----------------------------------------------------
// forward definition of class
// -----------------------------------------------------
class tCQMRegionGraph;

// -----------------------------------------------------
// tCQMRegionGraphVertex
// -------------------
class tCQMRegionGraphVertex : public tCQMDualGraphVertexAttrib {
  private:
    // pointer to the region graph this vertex belongs to.
    tCQMRegionGraph* m_rg;

    // array of pointers to adjacent vertices.
    vector<tCQMRegionGraphVertex*> m_neigh;

    // array of pointers to constrained half-edges.
    vector<tCQMHalfEdge2*> m_hedge; 

    // maximum number of neighbors.
    int m_maxdg;

    // current vertex degree.
    int m_degree;

    // flag to indicate whether the region has a boundary edge.
    bool m_isouter;

    // pointer to a boundary half-edge (if  any) of the face dual to a
    // vertex of this region.
    tCQMHalfEdge2* m_he;

    // flag used by graph algorithms.
    bool m_marked;

    // working list of dual graph vertices that are dual of faces with
    // at least one constrained edge.
    list<tCQMDualGraphVertex*> m_lbdverts;

    // private methods
    void addNeighbor(tCQMRegionGraphVertex*, tCQMHalfEdge2*);
    void delNeighbor(tCQMRegionGraphVertex*);

 public:
    // constructors
    tCQMRegionGraphVertex(tCQMRegionGraph*);

    // destructor
    ~tCQMRegionGraphVertex() {}

    // get methods
    int getDegree() const;
    int getMaxDegree() const;
    tCQMRegionGraphVertex* getNeighbor(int) const;
    tCQMHalfEdge2* getDualEdge(int) const;
    tCQMRegionGraph* getRegionGraph() const;
    tCQMHalfEdge2* getBoundHalfEdge() const;

    // query method
    bool isOuterRegion() const;
    bool isMarked() const;

    // set method
    void setMark(bool);
    void setOuter(bool);
    void setBoundHalfEdge(tCQMHalfEdge2*);

    friend class tCQMRegionGraph;
};

// -----------------------------------------------------
// tCQMRegionGraph
// ---------------
class tCQMRegionGraph {
  private:
    // list of vertices of this graph
    list<tCQMRegionGraphVertex*> m_lverts;    
    tCQMDualGraph* m_dual;

    // modifiers
    void addVert(tCQMRegionGraphVertex*);
    void delVert(tCQMRegionGraphVertex*);

    void addEdge(tCQMRegionGraphVertex*, tCQMDualGraphVertex*,
                 tCQMRegionGraphVertex*, tCQMDualGraphVertex*);
    void delEdge(tCQMRegionGraphVertex*, tCQMRegionGraphVertex*);
    void create_vertices();
    void create_edges();

  public:
    // constructors
    tCQMRegionGraph(tCQMDualGraph*);
    
    // destructor
    ~tCQMRegionGraph();

    // modifiers
    void getListOfHalfEdges(list<tCQMHalfEdge2*>&);

    // get methods
    int getNumVerts();
    list<tCQMRegionGraphVertex*>* getVertices();

    // modifier
    void clearMark();
};


// -----------------------------------------------------
// tCQMRegionGraphVertex - inline implementation
// ---------------------

//
// public methods
//

// get methods
inline int
tCQMRegionGraphVertex::getDegree() const
{
  return m_degree;
}

inline int
tCQMRegionGraphVertex::getMaxDegree() const
{
  return m_maxdg;
}

inline tCQMRegionGraphVertex*
tCQMRegionGraphVertex::getNeighbor(int i) const
{
  return m_neigh[i];
}

inline tCQMHalfEdge2*
tCQMRegionGraphVertex::getDualEdge(int i) const
{
  return m_hedge[i];
}

inline tCQMRegionGraph*
tCQMRegionGraphVertex::getRegionGraph() const
{
  return m_rg;
}

inline tCQMHalfEdge2*
tCQMRegionGraphVertex::getBoundHalfEdge() const
{
  return m_he;
}

// query method
inline bool
tCQMRegionGraphVertex::isOuterRegion() const
{
  return m_isouter;
}

inline bool
tCQMRegionGraphVertex::isMarked() const
{
  return m_marked;
}


// set method
inline void
tCQMRegionGraphVertex::setMark(bool mark)
{
  m_marked = mark;
}

inline void
tCQMRegionGraphVertex::setOuter(bool outer)
{
  m_isouter = outer;
}

inline void
tCQMRegionGraphVertex::setBoundHalfEdge(tCQMHalfEdge2* he)
{
  m_he = he;
}


// -----------------------------------------------------
// tCQMRegionGraph - inline implementation
// ---------------

//
// public methods
//

// get methods
inline int
tCQMRegionGraph::getNumVerts()
{
  return m_lverts.size();
}

inline list<tCQMRegionGraphVertex*>*
tCQMRegionGraph::getVertices()
{
  return &m_lverts;
}

//
// private method
//

// modifiers
inline void
tCQMRegionGraph::addVert(tCQMRegionGraphVertex* v)
{
  m_lverts.push_back(v);
}

inline void
tCQMRegionGraph::delEdge(tCQMRegionGraphVertex* v1, tCQMRegionGraphVertex* v2)
{
  v1->delNeighbor(v2);
  v2->delNeighbor(v1);
}

#endif   // __regiongraph_h










