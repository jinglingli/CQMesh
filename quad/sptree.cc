//  ===================================================================
//
//  Program:   Convex Quadrilateral Mesh Generator (CQMesh)
//  Module:    sptree.cc
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

#include "sptree.h"

// -------------------------------------------------------------------
// tCQMSpanningTreeVertex - implementation
// ----------------------

// constructors
tCQMSpanningTreeVertex::tCQMSpanningTreeVertex(tCQMDualGraphVertex* v)
{
  m_vert = v;
  m_parent = 0;
  m_nchildren = 0;
  m_maxdg = DEFAULT_TREE_DEGREE;
  m_child.resize(m_maxdg-1);
  m_haspoint = false;
  m_point1 = tCQMIndVertex2(-1,Point2(0,0));
  m_nextpt = Point2(0,0);
  m_isquad = false;
  m_ispent = false;
  m_point2 = tCQMIndVertex2(-1,Point2(0,0));
  m_point3 = tCQMIndVertex2(-1,Point2(0,0));
  m_he1 = 0;
  m_he2 = 0;
}

// destructor
tCQMSpanningTreeVertex::~tCQMSpanningTreeVertex()
{
  for (int i=0; i<m_nchildren; i++) {
    delete m_child[i];
  }
}

// get method
tCQMSpanningTreeVertex* 
tCQMSpanningTreeVertex::getSibling(tCQMSpanningTreeVertex* v) const
{
  for (int i=0; i<m_nchildren; i++) {
    if (m_child[i] == v) {
      return m_child[(i+1) % m_nchildren];
    }
  }

  return 0;
}

// query method
bool
tCQMSpanningTreeVertex::isInSubtree(tCQMSpanningTreeVertex* v) const
{
  if (this == v) return true;

  for (int i=0; i<getNumChildren(); i++) {
    if (getChild(i)->isInSubtree(v)) return true;
  }

  return false;
}

//
// private methods
//

// modifier
void
tCQMSpanningTreeVertex::addChild(tCQMSpanningTreeVertex* child)
{
  if (m_nchildren == m_maxdg-1) {
    m_child.resize(m_child.capacity() << 1);
  }
 
  m_child[m_nchildren] = child;
  ++m_nchildren;
}

void
tCQMSpanningTreeVertex::dettach()
{
  if (m_parent == 0) return;

  int i=0;
  while (m_parent->m_child[i] != this) {
    i++;
  }
  m_child[i] = 0;

  for (int j=i; j<m_parent->m_nchildren-1; j++) {
    m_parent->m_child[j] = m_parent->m_child[j+1];
  }

  --m_parent->m_nchildren;
}

// -------------------------------------------------------------------
// tCQMSpanningTree - implementation
// ----------------

tCQMSpanningTree::tCQMSpanningTree(tCQMSpanningTreeVertex* root)
{
  m_root = root;
}

// destructor
tCQMSpanningTree::~tCQMSpanningTree()
{
  if (m_root) delete m_root;
}

// modifier
void
tCQMSpanningTree::delLeaf(tCQMSpanningTreeVertex* v)
{
  if (isInTree(v) && (v->getNumChildren() == 0)) {
    v->dettach();
    if (m_root == v) m_root = 0;
    delete v;
  }
}

// get methods
std::list<std::list<tCQMSpanningTreeVertex*>*>*
tCQMSpanningTree::getVertPerLevel()
{
  if (m_root == 0) return 0;

  std::list<std::list<tCQMSpanningTreeVertex*>*>* llist = 
    (std::list<std::list<tCQMSpanningTreeVertex*>*>*) 
    new std::list<std::list<tCQMSpanningTreeVertex*>*>();
  std::list<tCQMSpanningTreeVertex*>* vlist = 
    (std::list<tCQMSpanningTreeVertex*>*) 
    new std::list<tCQMSpanningTreeVertex*>();

  // insert list of vertices at level  0 into list of list of vertices
  // per level.
  vlist->push_back(m_root);
  llist->push_front(vlist);

  while (!vlist->empty()) 
  {
    std::list<tCQMSpanningTreeVertex*>* vtemp = 
      (std::list<tCQMSpanningTreeVertex*>*) 
      new std::list<tCQMSpanningTreeVertex*>();

    std::list<tCQMSpanningTreeVertex*>::iterator vi;
    for (vi= vlist->begin(); vi != vlist->end(); ++vi) {
      int nc = (*vi)->getNumChildren();
      for (int i=0; i<nc; i++) {
	vtemp->push_back((*vi)->getChild(i));
      }
    }

    if (!vtemp->empty()) {
      llist->push_front(vtemp);
    }

    vlist = vtemp;
  }

  delete vlist;

  return llist;
}


void tCQMSpanningTreeVertex::dump(std::ofstream& myfile, int type) {
  tCQMDualGraphVertex* dual_vertex = getVertex();
  tCQMFace2* face = dual_vertex->getFace();
  Point2 barycenter = face->getBarycenter();
  //std::cout << barycenter.x() << "," << barycenter.y();
    if(type == 0){
        myfile << barycenter.x() << "," << barycenter.y();
    } else {
        myfile << "sptree, ";
        myfile << barycenter.x() << "," << barycenter.y() << "," << type;
        myfile << endl;
    }
}

void tCQMSpanningTree::dump(string fname) {
   std::ofstream myfile;
   myfile.open(fname);
   //myfile << "Writing this to a file.\n";
    
  std::vector<tCQMSpanningTreeVertex*> *to_process = new std::vector<tCQMSpanningTreeVertex*>();
  to_process->push_back(m_root);
  while (!to_process->empty()) {
    tCQMSpanningTreeVertex* cur = to_process->back();
    to_process->pop_back();
    int num = cur->getNumChildren();
    for (int i = 0; i < num; i++) {
      to_process->push_back(cur->getChild(i));
    }
    // print current node ..
    //std::cout << "sptree, "; // easier for filter output.
    myfile << "sptree, ";
    cur->dump(myfile, 0);
    // and its parent.
    if (cur->m_parent) {
      //std::cout << ", \t";
      myfile << ", \t";
      cur->m_parent->dump(myfile, 0);
    }
    //std::cout << endl;
      myfile << endl;
  }
    myfile.close();
}
