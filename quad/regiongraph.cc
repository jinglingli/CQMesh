//  ===================================================================
//
//  Program:   Convex Quadrilateral Mesh Generator (CQMesh)
//  Module:    regiongraph.cc
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

#include "regiongraph.h"

#include "edge2.h"
#include "face2.h"


#include <iostream>


// -----------------------------------------------------
// tCQMDualGraphVertex - implementation
// -------------------

// constructor
tCQMRegionGraphVertex::tCQMRegionGraphVertex(tCQMRegionGraph* g) :
  m_rg(g),
  m_maxdg(DEFAULT_REGION_GRAPH_DEGREE),
  m_degree(0),
  m_isouter(false),
  m_he(0),
  m_marked(false)
{
  m_neigh.resize(m_maxdg);
  m_hedge.resize(m_maxdg);

}

//
// private method
//

// modifier
void
tCQMRegionGraphVertex::addNeighbor(
                                   tCQMRegionGraphVertex* neigh,
                                   tCQMHalfEdge2* he
                                  )
{
  if (m_degree == m_maxdg) {
    m_neigh.resize(m_neigh.capacity() << 1);
    m_hedge.resize(m_hedge.capacity() << 1);
  }
 
  m_neigh[m_degree] = neigh;
  m_hedge[m_degree] = he;
  
  ++m_degree;
}

void
tCQMRegionGraphVertex::delNeighbor(tCQMRegionGraphVertex* v)
{
  int i=0;
  while ((i < m_degree) && (m_neigh[i] != v)) {
    i++;
  }

  if (i < m_degree) {
    for (int j=i; j<m_degree-1; j++) {
      m_neigh[j] = m_neigh[j+1];
      m_hedge[j] = m_hedge[j+1];
    }

    --m_degree; --v->m_degree;
  }
}


// -----------------------------------------------------
// tCQMRegionGraph - implementation
// ---------------

// constructor
tCQMRegionGraph::tCQMRegionGraph(tCQMDualGraph* dual):
  m_dual(dual)
{
  // create the vertices of the region graph.
  create_vertices();

  // create the edges of the region graph.
  create_edges();
}


// destructor
tCQMRegionGraph::~tCQMRegionGraph()
{
  list<tCQMRegionGraphVertex*>::iterator vi;

  vi = m_lverts.begin();
  while (vi != m_lverts.end()) {
    delete *vi;
    vi = m_lverts.erase(vi);
  }
}


// modifiers
void
tCQMRegionGraph::create_vertices()
{
  // get the list of all vertices of the dual graph.
  list<tCQMDualGraphVertex*>* vlist = m_dual->getVertices();

  // loop over  the list of  vertices of the  dual graph and  for each
  // vertex "v", find all vertices that are reachable from "v" through
  // a  path of  unconstrained edges  ("v" and  these vertices  form a
  // "region").
  list<tCQMDualGraphVertex*>::iterator vit;
  for (vit = vlist->begin(); vit != vlist->end(); vit++) {
    tCQMDualGraphVertex* v = *vit;

    // if a dual graph vertex is unmarked then it has not been visited
    // yet and hence  it can be the first vertex  of a distinct region
    // of the dual graph.
    if (!v->isMarked()) {
      // create a new region graph vertex.
      tCQMRegionGraphVertex* rgv = new tCQMRegionGraphVertex(this);

      // insert vertex into the list of vertices of the region graph.
      this->addVert(rgv);

      // mark  "v" and make  it a  vertex of  the newly  created graph
      // region.
      v->setMark(true);
      v->setDualGraphVertexAttrib(rgv);

      // insert "v"  into a temporary  list of all  vertices reachable
      // from "v" through a path of unconstrained edges ("v" and these
      // vertices form a "region").
      list<tCQMDualGraphVertex*> rgvlist;

      rgvlist.push_back(v);

      // perform  a breadth first  search (BFS)  to find  all vertices
      // reachable from  "v" through a path of  unconstrained edges of
      // the dual graph.
      while (!rgvlist.empty()) {
        // remove the first vertex of the list.
        tCQMDualGraphVertex* v1 = rgvlist.back();
        rgvlist.pop_back();

        // if  the  dual graph  vertex  is the  dual  of  a face  that
        // contains  a boundary  edge, then  the  corresponding region
        // graph vertex is considered an "outer" vertex; so, the outer
        // flag is set to true.
        if (!rgv->isOuterRegion()) {
          tCQMHalfEdge2* he1 = v1->getFace()->getHalfEdge();
          if (he1->getEdge()->isOnBoundary()) {
            rgv->setOuter(true);
            rgv->setBoundHalfEdge(he1);
          }
          else {
            he1 = he1->getNext();
            if (he1->getEdge()->isOnBoundary()) {
              rgv->setOuter(true);
              rgv->setBoundHalfEdge(he1);
            }
            else {
              he1 = he1->getNext();
              if (he1->getEdge()->isOnBoundary()) {
                rgv->setOuter(true);
                rgv->setBoundHalfEdge(he1);
              }
            }
          }
        }

        // find  all  neighbors of  the  vertex  that  do not  form  a
        // constrained edge with the vertex, and then insert them into
        // the list.
        int i;
        for (i=0; i<v1->getDegree(); i++) {
          // if the  "i"-th neighbor does not form  a constrained edge
          // of the  dual graph  with "v1", insert  it into  the list;
          // otherwise, "v1" is a boundary vertex of the region and we
          // insert "v1"  into a list of boundary  vertices for future
          // use.
          if (!v1->isDualofConstrainedEdge(i)) {
            tCQMDualGraphVertex* v2 = v1->getNeighbor(i);
            // the  vertex may  have been  visited before;  if  so, we
            // ignore it.
            if (!v2->isMarked()) {
              v2->setMark(true);
              v2->setDualGraphVertexAttrib(rgv);
              rgvlist.push_back(v2);
            }
          }
          else {
            // "v1"  is a  vertex of  a constrained  edge of  the dual
            // graph, so we insert it into a list of boundary vertices
            // for future use.
            if (rgv->m_lbdverts.front() != v1) {
              rgv->m_lbdverts.push_front(v1);
            }
          }
        }
      }
    }
  }

  // clear the "marks" of all vertices of the dual graph.
  m_dual->clearMark();
}

void
tCQMRegionGraph::create_edges()
{
  // for each region  of the dual graph (i.e., for  each vertex "v" of
  // the  region  graph)"v", find  all  neighboring  regions and  then
  // create the graph edges.
  list<tCQMRegionGraphVertex*>::iterator vit;
  for (vit = m_lverts.begin(); vit != m_lverts.end(); vit++) {
    tCQMRegionGraphVertex* v = *vit;

    // get the list of boundary vertices of the region associated with
    // "v".
    while (!v->m_lbdverts.empty()) {
      tCQMDualGraphVertex* v1 = v->m_lbdverts.front();
      v->m_lbdverts.pop_front();

      // for each  neighbor "v2" of "v1",  if "v1" and  "v2" belong to
      // distinct  regions then  create an  edge in  the  region graph
      // (unless such an edge has already been created).
      int i;
      for (i=0; i<v1->getDegree(); i++) {
        if (v1->isDualofConstrainedEdge(i)) {
          tCQMDualGraphVertex* v2 = v1->getNeighbor(i);
          tCQMRegionGraphVertex* vaux = dynamic_cast<tCQMRegionGraphVertex*>
            (v2->getDualGraphVertexAttrib());

          // check if the regions "v"  and "vaux" are distinct (as two
          // vertices  of  the  same  region  may be  connected  by  a
          // constrained edge).
          if (v != vaux) {
            // create an edge of the region graph corresponding to the
            // pair  (v,  vaux)  unless  this edge  has  already  been
            // created.
            bool found = false;
            int j = 0;
            while ((j < v->getDegree()) && !found) {
              if (v->getNeighbor(j) == vaux) {
                found = true;
              }
              else {
                j++;
              }
            }
          
            if (!found) {
              this->addEdge(v,v1,vaux,v2);
            }
          }
        }
      }
    }
  }
}

void
tCQMRegionGraph::addEdge(
                         tCQMRegionGraphVertex* v1,
                         tCQMDualGraphVertex* v3,
                         tCQMRegionGraphVertex* v2, 
                         tCQMDualGraphVertex* v4
                        )
{
  tCQMEdge2* e = v3->getFace()->getCommonEdge(v4->getFace());

  tCQMHalfEdge2* he1 = e->getHalfEdge();
  tCQMHalfEdge2* he2;
  if (he1->getFace() == v3->getFace()) {
    he2 = e->getMate(he1);
  }
  else {
    he2 = he1;
    he1 = e->getMate(he1);
  }

  v1->addNeighbor(v2,he1);
  v2->addNeighbor(v1,he2);
}

void
tCQMRegionGraph::delVert(tCQMRegionGraphVertex* v)
{
  // for each neighbor, remove adjacency information
  int dg = v->getDegree();
  for (int i=0; i<dg; i++) {
    v->getNeighbor(i)->delNeighbor(v);
  }

  m_lverts.remove(v);
}

void
tCQMRegionGraph::getListOfHalfEdges(list<tCQMHalfEdge2*>& lhe)
{
  // for each "outer" vertex "v" of the region graph, start a BFS from
  // "v" to create a sorted list of regions.
  list<tCQMRegionGraphVertex*>::const_iterator vit;
  for (vit = m_lverts.begin(); vit != m_lverts.end(); vit++) {
    tCQMRegionGraphVertex* v = *vit;

    // if the region is an outer region, then perform a BFS from it to
    // find  all regions  path-connected  regions that  are not  outer
    // regions.
    if (v->isOuterRegion()) {
      // create a working list to perform the BFS.
      list<tCQMRegionGraphVertex*> rgvlist;
      rgvlist.push_back(v);

      lhe.push_front(v->getBoundHalfEdge());

      while (!rgvlist.empty()) {
        tCQMRegionGraphVertex* v1 = rgvlist.front();
        rgvlist.pop_front();

        int i;
        for (i=0; i<v1->getDegree(); i++) {
          tCQMRegionGraphVertex* v2 = v1->getNeighbor(i);
          if (!v2->isOuterRegion() && !v2->isMarked()) {
            v2->setMark(true);
            int j = 0;
            while (v2->getNeighbor(j) != v1) {
              ++j;
            }

            lhe.push_front(v2->getDualEdge(j));

            rgvlist.push_back(v2);
          }
        }
      }
    }
  }

  this->clearMark();
}

void 
tCQMRegionGraph::clearMark() {
  list<tCQMRegionGraphVertex*>::iterator vi;
  for (vi = m_lverts.begin(); vi != m_lverts.end(); ++vi) {
    (*vi)->setMark(false);
  }
}
