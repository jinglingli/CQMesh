//  ===================================================================
//
//  Program:   Convex Quadrilateral Mesh Generator (CQMesh)
//  Module:    complex2.cc
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

#include "complex2.h"

// -----------------------------------------------------
// tCQMComplex2
// ------------

// constructor
tCQMComplex2::tCQMComplex2() 
{
  m_dual = (tCQMDualGraph*) new tCQMDualGraph(this);
}

// destructor
tCQMComplex2::~tCQMComplex2() 
{
  // delete all vertices
  list<tCQMVertex2*>::iterator vi = m_lverts.begin();
  while (vi != m_lverts.end()) {
    tCQMVertex2* vaux = *vi;
    vi = m_lverts.erase(vi);
    delete vaux;
  }

  // delete half-edges and edges
  list<tCQMEdge2*>::iterator ei = m_ledges.begin();
  while (ei != m_ledges.end()) {
    tCQMEdge2* eaux = *ei;
    tCQMHalfEdge2* he1 = eaux->getHalfEdge();
    tCQMHalfEdge2* he2 = eaux->getMate(he1);

    delete he1;
    if (he2 != 0) delete he2;

    ei = m_ledges.erase(ei);

    delete eaux;
  }

  // delete all faces
  list<tCQMFace2*>::iterator fi = m_lfaces.begin();
  while (fi != m_lfaces.end()) {
    tCQMFace2* faux = *fi;
    fi = m_lfaces.erase(fi);
    delete faux;
  }

  // delete all boundaries
  list<tCQMBoundary2*>::iterator bi = m_lbounds.begin();
  while (bi != m_lbounds.end()) {
    tCQMBoundary2* baux = *bi;
    bi = m_lbounds.erase(bi);
    delete baux;
  }

  // delete dual graph
  delete m_dual;
}

// public methods

// ---------------------------------------------------------------
// Method createVertex
// ---------------------------------------------------------------
tCQMVertex2* 
tCQMComplex2::createVertex(const Point2& p) 
{
  // create vertex
  tCQMVertex2* newv = (tCQMVertex2*) new tCQMVertex2(p);

  // insert it into list of vertices of this complex
  m_lverts.push_back(newv);

  // update bounding box
  Point2 p1 = m_bb.getPMin();
  Point2 p2 = m_bb.getPMax();

  double xmin = p1.x();
  double ymin = p1.y();
  double xmax = p2.x();
  double ymax = p2.y();

  if (xmin > p.x()) xmin = p.x();
  if (ymin > p.y()) ymin = p.y();
  if (xmax < p.x()) xmax = p.x();
  if (ymax < p.y()) ymax = p.y();

  m_bb.setPMin(Point2(xmin,ymin));
  m_bb.setPMax(Point2(xmax,ymax));

  return newv;
}

// ---------------------------------------------------------------
// Method createDualVertex
// ---------------------------------------------------------------
void
tCQMComplex2::createDualVertex(tCQMFace2* f) 
{
  // create dual vertex
  tCQMDualGraphVertex* newv = (tCQMDualGraphVertex*) 
    new tCQMDualGraphVertex(this->getDualGraph(),f);

  // make face f point to its dual vertex in the dual graph
  f->setDualVertex(newv);

  // add vertex to dual graph
  m_dual->addVert(newv);
}

// ---------------------------------------------------------------
// Method deleteDualVertex()
// ---------------------------------------------------------------
void
tCQMComplex2::deleteDualVertex(tCQMFace2* f) 
{
  // get dual vertex
  tCQMDualGraphVertex* v = f->getDualVertex();

  // remove vertex from dual graph
  m_dual->delVert(v);

  // face "f" is not pointing to any vertex
  f->setDualVertex(0);

  // delete vertex
  delete v;
}

// ---------------------------------------------------------------
// Method createDualEdge
// ---------------------------------------------------------------
void
tCQMComplex2::createDualEdge(tCQMFace2* f1, tCQMFace2* f2) 
{
  // get dual vertex of f1 and f2
  tCQMDualGraphVertex* v1 = f1->getDualVertex();
  tCQMDualGraphVertex* v2 = f2->getDualVertex();

  // add edge to dual graph
  m_dual->addEdge(v1,v2);
}

// ---------------------------------------------------------------
// Method getNextStarVert
// ---------------------------------------------------------------
tCQMHalfEdge2*
tCQMComplex2::getNextStarVert(tCQMHalfEdge2* he)
{
  tCQMHalfEdge2* heaux = he->getEdge()->getMate(he);

  if (heaux != 0) return heaux->getNext();
  return heaux;
}

// ---------------------------------------------------------------
// Method getPrevStarVert
// ---------------------------------------------------------------
tCQMHalfEdge2*
tCQMComplex2::getPrevStarVert(tCQMHalfEdge2* he)
{
  tCQMHalfEdge2* heaux = he->getPrev();

  if (heaux != 0) return heaux->getEdge()->getMate(heaux);
  return heaux;
}

// ---------------------------------------------------------------
// Method getNextBoundEdge
// ---------------------------------------------------------------
tCQMEdge2*
tCQMComplex2::getNextBoundEdge(tCQMEdge2* e)
{
  tCQMHalfEdge2* he1 = e->getHalfEdge()->getNext();
  tCQMHalfEdge2* he2 = he1->getEdge()->getMate(he1);
  while (he2 != 0) {
    he1 = he2->getNext();
    he2 = he1->getEdge()->getMate(he1); 
  }

  return he1->getEdge();
}

// ---------------------------------------------------------------
// Method getPrevBoundEdge
// ---------------------------------------------------------------
tCQMEdge2*
tCQMComplex2::getPrevBoundEdge(tCQMEdge2* e)
{
  tCQMHalfEdge2* he1 = e->getHalfEdge()->getPrev();
  tCQMHalfEdge2* he2 = he1->getEdge()->getMate(he1);
  while (he2 != 0) {
    he1 = he2->getPrev();
    he2 = he1->getEdge()->getMate(he1); 
  }

  return he1->getEdge();
}

// ---------------------------------------------------------------
// Method mkFBEV
//
// It creates a new face with only one vertex, edge, and boundary, 
// which is used to create an isolated face.
// ---------------------------------------------------------------
tCQMFace2*
tCQMComplex2::mkFBEV(tCQMVertex2* v)
{
  // create boundary, face and edge
  tCQMBoundary2* newbd = (tCQMBoundary2*) new tCQMBoundary2();
  tCQMFace2* newf = (tCQMFace2*) new tCQMFace2();
  tCQMEdge2* newe = (tCQMEdge2*) new tCQMEdge2();

  // insert the new elements into the corresponding lists
  m_lbounds.push_back(newbd);
  m_lfaces.push_back(newf);
  m_ledges.push_back(newe);
  
  // create the first half-edge of the new edge
  tCQMHalfEdge2* newhe = (tCQMHalfEdge2*) new tCQMHalfEdge2();

  // initialize first half-edge
  newhe->setVertex(v);
  newhe->setEdge(newe);
  newhe->setFace(newf);
  newhe->setNext(newhe);
  newhe->setPrev(newhe);
  
  // update vertex
  v->setHalfEdge(newhe);
  v->setBoundHalfEdge(newhe);

  // initialize new edge
  newe->setHalfEdge1(newhe);
  newe->setBoundary(newbd);

  // initialize new face
  newf->setHalfEdge(newhe);
  newf->setComplex(this);

  // initialize new boundary
  newbd->setBoundaryEdge(newe);

  // create a corresponding vertex in the dual graph
  createDualVertex(newf);

  // returns the new face
  return newf;
}

// -------------------------------------------------------------------
// Method mkEV
//
// It adds a  new vertex and a  new edge to an existent  face. The new
// edge follows  a given  edge of  the face with  respect to  the face
// cycle.
// -------------------------------------------------------------------
tCQMEdge2*
tCQMComplex2::mkEV(
		   tCQMVertex2* v,
		   tCQMHalfEdge2* he
		   )
{
  // create edge
  tCQMEdge2* newe = (tCQMEdge2*) new tCQMEdge2();

  // insert the new edge into its corresponding list
  m_ledges.push_back(newe);

  // create the first half-edge of the new edge
  tCQMHalfEdge2* newhe = (tCQMHalfEdge2*) new tCQMHalfEdge2();

  // initialize first half-edge
  newhe->setVertex(v);
  newhe->setEdge(newe);
  newhe->setFace(he->getFace());

  // insert new half-edge in the face cycle
  tCQMHalfEdge2* henext = he->getNext();
  henext->setPrev(newhe);
  he->setNext(newhe);
  newhe->setNext(henext);
  newhe->setPrev(he);

  // update vertex
  v->setHalfEdge(newhe);
  v->setBoundHalfEdge(newhe);

  // initialize new edge
  newe->setHalfEdge1(newhe);

  // boundary is the same as that one of he
  newe->setBoundary(he->getEdge()->getBoundary());

  // return the new edge
  return newe;
}

// -------------------------------------------------------------------
// Method spEmkEV()
//
// It splits an edge by adding a new vertex inbetween its endpoints.
//
// -------------------------------------------------------------------
tCQMEdge2*
tCQMComplex2::spEmkEV(
		      tCQMVertex2* v,
		      tCQMHalfEdge2* he
		      )
{
  // create edge
  tCQMEdge2* newe = (tCQMEdge2*) new tCQMEdge2();

  // insert the new edge into its corresponding list
  m_ledges.push_back(newe);

  // create half-edges of the new edge
  tCQMHalfEdge2* newhe1 = (tCQMHalfEdge2*) new tCQMHalfEdge2();
  tCQMHalfEdge2* newhe2 = (tCQMHalfEdge2*) new tCQMHalfEdge2();

  // initialize first half-edge
  newhe1->setVertex(v);
  newhe1->setEdge(newe);
  newhe1->setFace(he->getFace());

  // initialize second half-edge
  tCQMHalfEdge2* hemate = he->getEdge()->getMate(he);
  newhe2->setVertex(hemate->getVertex());
  newhe2->setEdge(newe);
  newhe2->setFace(hemate->getFace());

  // insert first half-edge in the face cycle
  tCQMHalfEdge2* henext = he->getNext();
  henext->setPrev(newhe1);
  he->setNext(newhe1);
  newhe1->setNext(henext);
  newhe1->setPrev(he);

  // insert second half-edge in the face cycle
  tCQMHalfEdge2* heprev = hemate->getPrev();
  heprev->setNext(newhe2);
  hemate->setPrev(newhe2);
  newhe2->setNext(hemate);
  newhe2->setPrev(heprev);

  // update vertices
  hemate->setVertex(v);
  newhe2->getVertex()->setHalfEdge(newhe2);
  v->setHalfEdge(newhe1);
  v->setBoundHalfEdge(0);

  // initialize new edge
  newe->setHalfEdge1(newhe1);
  newe->setHalfEdge2(newhe2);

  // boundary is the same as that one of he
  newe->setBoundary(0);

  // update dual gaph to include one more edge
  createDualEdge(newhe1->getFace(),newhe2->getFace());

  // return the new edge
  return newe;
}

// ---------------------------------------------------------------
// Method mkFE
//
// It creates a new face and a new edge. The new face is incident
// to an edge of an existent face. The half-edge of the commom
// edge in the adjacent face is given as an input parameter.
// ---------------------------------------------------------------
tCQMFace2*
tCQMComplex2::mkFE(tCQMHalfEdge2* he)
{
  // create new edge and new face
  tCQMEdge2* newe = (tCQMEdge2*) new tCQMEdge2();
  tCQMFace2* newf = (tCQMFace2*) new tCQMFace2();

  // insert the new elements into the corresponding lists
  m_ledges.push_back(newe);
  m_lfaces.push_back(newf);

  // get common face
  tCQMFace2* face = he->getFace();

  // get commom edge
  tCQMEdge2* edge = he->getEdge();

  // get commom vertices
  tCQMVertex2* v1 = he->getNext()->getVertex();
  tCQMVertex2* v2 = he->getVertex();

  // create the first half-edge of the new edge
  tCQMHalfEdge2* newhe1 = (tCQMHalfEdge2*) new tCQMHalfEdge2();

  // create the second half-edge of the commom edge
  tCQMHalfEdge2* newhe2 = (tCQMHalfEdge2*) new tCQMHalfEdge2();

  // initialize first half-edge of the new edge
  newhe1->setVertex(v2);
  newhe1->setEdge(newe);
  newhe1->setFace(newf);
  newhe1->setNext(newhe2);
  newhe1->setPrev(newhe2);

  // initialize second half-edge of the commom edge
  newhe2->setVertex(v1);
  newhe2->setEdge(edge);
  newhe2->setFace(newf);
  newhe2->setNext(newhe1);
  newhe2->setPrev(newhe1);
  
  // initialize new edge
  newe->setHalfEdge1(newhe1);

  // it is on the same boundary as edge
  newe->setBoundary(edge->getBoundary());

  // update vertex
  v2->setBoundHalfEdge(newhe1);

  // initialize new face
  newf->setHalfEdge(newhe1);
  newf->setComplex(this);

  // update the boundary
  edge->getBoundary()->setBoundaryEdge(newe);

  // update commom edge
  edge->setHalfEdge2(newhe2);
  edge->setBoundary(0);

  // create a corresponding vertex in the dual graph
  createDualVertex(newf);
  
  // create an edge in the dual graph corresponding to
  // the common edge between newf and  face
  createDualEdge(face,newf);

  // return the new face
  return newf;
}

// -------------------------------------------------------------------
// Method mkFE
//
// It creates a new face and a new edge by splitting an existent face.
//
// -------------------------------------------------------------------
tCQMFace2*
tCQMComplex2::mkFE(tCQMHalfEdge2* he1, tCQMHalfEdge2* he2)
{
  // create new edge and new face
  tCQMEdge2* newe = (tCQMEdge2*) new tCQMEdge2();
  tCQMFace2* newf = (tCQMFace2*) new tCQMFace2();

  // insert the new elements into the corresponding lists
  m_ledges.push_back(newe);
  m_lfaces.push_back(newf);

  // get common face
  tCQMFace2* face = he1->getFace();

  // remove dual vertex of "face" from the dual graph
  deleteDualVertex(face);

  // get commom vertices
  tCQMVertex2* v3 = he1->getVertex();
  tCQMVertex2* v4 = he2->getVertex();

  // create the first half-edge of the new edge
  tCQMHalfEdge2* newhe1 = (tCQMHalfEdge2*) new tCQMHalfEdge2();

  // create the second half-edge of the commom edge
  tCQMHalfEdge2* newhe2 = (tCQMHalfEdge2*) new tCQMHalfEdge2();

  // initialize first half-edge of the new edge
  newhe1->setVertex(v4);
  newhe1->setEdge(newe);
  newhe1->setFace(face);

  // initialize second half-edge of the commom edge
  newhe2->setVertex(v3);
  newhe2->setEdge(newe);
  newhe2->setFace(newf);

  // split cycle of half-edges
  newhe1->setNext(he1);
  newhe1->setPrev(he2->getPrev());
  newhe2->setNext(he2);
  newhe2->setPrev(he1->getPrev());
  he2->getPrev()->setNext(newhe1);
  he1->getPrev()->setNext(newhe2);
  he2->setPrev(newhe2);
  he1->setPrev(newhe1);
  
  // initialize new edge
  newe->setHalfEdge1(newhe1);
  newe->setHalfEdge2(newhe2);
  newe->setBoundary(0);

  // initialize new face
  newf->setHalfEdge(newhe2);
  newf->setComplex(this);

  // update old face
  face->setHalfEdge(he1);

  // create corresponding vertices in the dual graph
  createDualVertex(face);
  createDualVertex(newf);
  
  // create all  edges corresponding to  the adjacent faces  to "face"
  // and "newf" in the dual  graph. Also, update half-edges in the new
  // face.
  tCQMHalfEdge2* heaux1 = face->getHalfEdge();
  tCQMHalfEdge2* heaux2 = heaux1;
  do {
    tCQMHalfEdge2* heaux3 = heaux2->getEdge()->getMate(heaux2);
    if (heaux3 != 0) {
      createDualEdge(face,heaux3->getFace());
    }
    heaux2 = heaux2->getNext();
  } while (heaux1 != heaux2);

  heaux2 = he2;
  do {
    // update pointer to new face  of all half-edges in "newf", except
    // for "newhe2" that is already set.
    heaux2->setFace(newf);

    heaux1 = heaux2->getEdge()->getMate(heaux2);
    if (heaux1 != 0) {
      createDualEdge(newf,heaux1->getFace());
    }
    heaux2 = heaux2->getNext();
  } while (heaux2 != newhe2);

  // return the new face
  return newf;
}

// -------------------------------------------------------------------
// Method spFmkFE()
// -------------------------------------------------------------------
tCQMVertex2*
tCQMComplex2::spFmkFE(tCQMHalfEdge2* he, const Point2& p) 
{
  // create a  new vertice and  a new edge  in the face  contained the
  // half-edge "he". The  new vertex has location at  "p", and the new
  // edge has endpoints  "p" and "q", where "q"  is the initial vertex
  // of the half-edge following "he" in its face.

  // if new edge is not on the boundary, split the face that shares it
  // with the face containing "he".
  tCQMEdge2* e;
  tCQMVertex2 *v = createVertex(p);
  if (!he->getEdge()->isOnBoundary()) {
    // create new vertex and new edge
    e = this->spEmkEV(v,he);

    // get the mate of he
    tCQMHalfEdge2* he1 = he->getEdge()->getMate(he);

    // get a half-edge  that has an initial vertex that  is not in the
    // common edge.
    tCQMHalfEdge2* he2 = he1->getNext()->getNext();
    
    // split the face containing "he1"  and "he2" by inserting an edge
    // linking their initial vertex.
    this->mkFE(he1,he2);
  }
  else {
    // create new vertex and new edge
    e = this->mkEV(v,he);
  }

  e->constraintOn();

  return v;
}

// ---------------------------------------------------------------
// Method rmBEV
//
// It glues two connected components of a complex by removing an
// edge, two vertice and one boundary.
// ---------------------------------------------------------------
void
tCQMComplex2::rmBEV(
		    tCQMHalfEdge2* he1, 
		    tCQMHalfEdge2* he2
		    )
{
  // get initial vertex of he1, he2, he1->next and he2->next
  tCQMVertex2* v1 = he1->getVertex();
  tCQMVertex2* v2 = he1->getNext()->getVertex();
  tCQMVertex2* v3 = he2->getVertex();
  tCQMVertex2* v4 = he2->getNext()->getVertex();

  // get the edge containing he1 and the edge containing he2
  tCQMEdge2* e1 = he1->getEdge();
  tCQMEdge2* e2 = he2->getEdge();

  // get the face containing he1 and the face containing he2
  tCQMFace2* f1 = he1->getFace();
  tCQMFace2* f2 = he2->getFace();

  // get the boundary containing he1 and the boundary containing he2
  tCQMBoundary2* bd1 = e1->getBoundary();
  tCQMBoundary2* bd2 = e2->getBoundary();

  // merge boundary
  //
  // set boundary edge pointer in face containing he2 to
  // point to the boundary containing he1
  tCQMEdge2* eaux = this->getNextBoundEdge(e2);
  while (eaux != e2) {
    eaux->setBoundary(bd1);
    eaux = this->getNextBoundEdge(eaux);
  }

  // if boundary containing he1 points to the edge containing
  // he1, make it point to other edge in the same boundary.
  //
  // get next edge in the boundary
  if (e1 == bd1->getBoundaryEdge()) bd1->setBoundaryEdge(this->getNextBoundEdge(e1));

  // update every half-edge with initial vertex v3 and v4
  tCQMHalfEdge2* heaux = this->getNextStarVert(v3->getHalfEdge());
  while (heaux != 0) {
    heaux->setVertex(v2);
    heaux = this->getNextStarVert(heaux);
  }

  heaux = this->getPrevStarVert(v3->getHalfEdge());
  while (heaux != 0)  {
    heaux->setVertex(v2);
    heaux = this->getPrevStarVert(heaux);
  }
  
  v3->getHalfEdge()->setVertex(v2);

  heaux = this->getNextStarVert(v4->getHalfEdge());
  while (heaux != 0) {
    heaux->setVertex(v1);
    heaux = this->getNextStarVert(heaux);
  }

  heaux = this->getPrevStarVert(v4->getHalfEdge());
  while (heaux != 0)  {
    heaux->setVertex(v1);
    heaux = this->getPrevStarVert(heaux);
  }

  v4->getHalfEdge()->setVertex(v1);

  // make e1 the edge of he2
  he2->setEdge(e1);

  // make he2 the mate of he1
  e1->setHalfEdge2(he2);

  // e1 is no longer a boundary edge
  e1->setBoundary(0);

  // remove vertices
  m_lverts.remove(v3);
  m_lverts.remove(v4);
  delete v3;
  delete v4;

  // remove edge that contained he2
  m_ledges.remove(e2);
  delete e2;

  // remove boundary that contained he2
  m_lbounds.remove(bd2);
  delete bd2;

  // create an edge in the dual graph corresponding to
  // the common edge between newf and  face
  createDualEdge(f1,f2);
}

// ---------------------------------------------------------------
// Method mkGrmEV
//
// It glues two faces of the same connected component, and removes
// an edge and at most two vertices. It can also make a genus by
// spliting one boundary.
// ---------------------------------------------------------------
void
tCQMComplex2::mkGrmEV(
		      tCQMHalfEdge2* he1, 
		      tCQMHalfEdge2* he2
		      )
{
  // get initial vertex of he1, he2, he1->next and he2->next
  tCQMVertex2* v1 = he1->getVertex();
  tCQMVertex2* v2 = he1->getNext()->getVertex();
  tCQMVertex2* v3 = he2->getVertex();
  tCQMVertex2* v4 = he2->getNext()->getVertex();

  // get the edge containing he1 and the edge containing he2
  tCQMEdge2* e1 = he1->getEdge();
  tCQMEdge2* e2 = he2->getEdge();

  // get the face containing he1 and the face containing he2
  tCQMFace2* f1 = he1->getFace();
  tCQMFace2* f2 = he2->getFace();

  // get the boundary containing he1 and he2
  tCQMBoundary2* bd = e1->getBoundary();

  // get edges preceding and succeeding e1 on the boundary bd
  tCQMEdge2* eprev = this->getPrevBoundEdge(e1);
  tCQMEdge2* enext = this->getNextBoundEdge(e1);

  if (eprev == enext) {
    // boundary has only two edges: e1 and e2
    // we must remove the boundary

    m_lbounds.remove(bd);
    delete bd;

    // update initial vertex of he1 and he2
    he1->getVertex()->setBoundHalfEdge(0);
    he1->getNext()->getVertex()->setBoundHalfEdge(0);
  }
  else if (eprev == e2) {
    // he1 follows he2 on the boundary bd

    // update edge pointer of boundary bd
    bd->setBoundaryEdge(enext);

    // update initial vertex of he1
    he1->getVertex()->setBoundHalfEdge(0);
  }
  else if (enext == e2) {
    // he1 precedes he2 on the boundary bd

    // update edge pointer of boundary bd
    bd->setBoundaryEdge(eprev);

    // update initial vertex of he1->next
    he1->getNext()->getVertex()->setBoundHalfEdge(0);
  }
  else {
    // boundary bd must be split into two ones

    // create another boundary
    tCQMBoundary2* bdaux = (tCQMBoundary2*) new tCQMBoundary2();

    // insert it into list of boundaries of the complex
    m_lbounds.push_back(bdaux);

    // set first edge of the new boundary
    bdaux->setBoundaryEdge(enext);

    // update edges of the new boundary
    tCQMEdge2* eaux = enext;
    while (eaux != e2) {
      eaux->setBoundary(bdaux);
      eaux = this->getNextBoundEdge(eaux);
    }

    // update edge pointer of the old boundary
    bd->setBoundaryEdge(eprev);
  }

  // update every half-edge with initial vertex v3 and v4
  tCQMHalfEdge2* heaux = this->getNextStarVert(v3->getHalfEdge());
  while (heaux != 0) {
    heaux->setVertex(v2);
    heaux = this->getNextStarVert(heaux);
  }

  heaux = this->getPrevStarVert(v3->getHalfEdge());
  while (heaux != 0)  {
    heaux->setVertex(v2);
    heaux = this->getPrevStarVert(heaux);
  }

  v3->getHalfEdge()->setVertex(v2);

  heaux = this->getNextStarVert(v4->getHalfEdge());
  while (heaux != 0) {
    heaux->setVertex(v1);
    heaux = this->getNextStarVert(heaux);
  }

  heaux = this->getPrevStarVert(v4->getHalfEdge());
  while (heaux != 0)  {
    heaux->setVertex(v1);
    heaux = this->getPrevStarVert(heaux);
  }

  v4->getHalfEdge()->setVertex(v1);

  // make e1 the edge of he2
  he2->setEdge(e1);

  // make he2 the mate of he1
  e1->setHalfEdge2(he2);

  // e1 is no longer a boundary edge
  e1->setBoundary(0);

  // remove vertices if needed
  if (v3 != v2) {
    m_lverts.remove(v3);
    delete v3;
  }

  if (v4 != v1) {
    m_lverts.remove(v4);
    delete v4;
  }

  // remove edge that contained he2
  m_ledges.remove(e2);
  delete e2;

  // create an edge in the dual graph corresponding to
  // the common edge between newf and  face
  createDualEdge(f1,f2);
}
