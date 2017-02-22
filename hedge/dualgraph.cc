//  ===================================================================
//
//  Program:   Convex Quadrilateral Mesh Generator (CQMesh)
//  Module:    dualgraph.cc
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

#include "dualgraph.h"

// -----------------------------------------------------
// tCQMDualGraphVertex - implementation
// -------------------

// constructor
tCQMDualGraphVertex::tCQMDualGraphVertex(
					 tCQMDualGraph* g, 
					 tCQMFace2* f
					 )
{
  m_dual = g;
  m_face = f;
  m_degree = 0;
  m_maxdg = DEFAULT_GRAPH_DEGREE;
  m_neigh.resize(m_maxdg);
  m_marked = false;
  m_dgvattrib = 0;
}

//
// private method
//

// modifier
void
tCQMDualGraphVertex::addNeighbor(tCQMDualGraphVertex* neigh)
{
  if (m_degree == m_maxdg) m_neigh.resize(m_neigh.capacity() << 1);
 
  m_neigh[m_degree] = neigh;
  ++m_degree;
}

void
tCQMDualGraphVertex::delNeighbor(tCQMDualGraphVertex* v)
{
  int i=0;
  while ((i < m_degree) && (m_neigh[i] != v)) {
    i++;
  }

  if (i < m_degree) {
    for (int j=i; j<m_degree-1; j++) {
      m_neigh[j] = m_neigh[j+1];
    }

    --m_degree; --v->m_degree;
  }
}

// -----------------------------------------------------
// tCQMDualGraph - implementation
// -------------

// constructor
tCQMDualGraph::tCQMDualGraph(tCQMComplex2* mesh):
  m_mesh(mesh)
{}

// destructor
tCQMDualGraph::~tCQMDualGraph()
{
  list<tCQMDualGraphVertex*>::iterator vi;

  vi = m_lverts.begin();
  while (vi != m_lverts.end()) {
    delete *vi;
    vi = m_lverts.erase(vi);
  }
}

// modifiers
void
tCQMDualGraph::delVert(tCQMDualGraphVertex* v)
{
  // for each neighbor, remove adjacency information
  int dg = v->getDegree();
  for (int i=0; i<dg; i++) {
    v->getNeighbor(i)->delNeighbor(v);
  }

  m_lverts.remove(v);
}

void 
tCQMDualGraph::clearMark() {
  list<tCQMDualGraphVertex*>::iterator vi;
  for (vi = m_lverts.begin(); vi != m_lverts.end(); ++vi) {
    (*vi)->setMark(false);
  }
}
