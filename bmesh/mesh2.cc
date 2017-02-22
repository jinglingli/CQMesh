//  ===================================================================
//
//  Program:   Convex Quadrilateral Mesh Generator (CQMesh)
//  Module:    mesh2.cc
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

#include "mesh2.h"

// -------------------------------------------------------------------
// tCQMMesh2Builder
// -------------------------------------------------------------------

// constructors
tCQMMesh2Builder::tCQMMesh2Builder()
{
  m_cp = (tCQMComplex2*) new tCQMComplex2();
}

// destructor
tCQMMesh2Builder::~tCQMMesh2Builder()
{
  m_cp = 0;
}

//
// searchers
//

// -------------------------------------------------------------------
//
// Method findVertexInBound()
//
// It returns true if a given vertex belongs to the same boundary as a
// given edge. Otherwise, it returns false.
//
// -------------------------------------------------------------------
bool
tCQMMesh2Builder::findVertexInBound(
				    int cv, 
				    const tCQMEdgeKey& ce
				    )
                                    const
{
  VertexSet::const_iterator vspos;

  vspos = vs.find(cv);
  if (vspos != vs.end()) {
    tCQMHalfEdge2* he1 = (es.find(ce))->second;
    tCQMHalfEdge2* he2;
    while ((vspos != vs.end()) && (vspos->first == cv)) {
      he2 = (vspos->second)->getBoundHalfEdge();
      if (he1->getEdge()->getBoundary() == 
	  he2->getEdge()->getBoundary()) {
	return true;
      }
      else {
	++vspos;
      }
    }
  }

  return false;
}

//
// query
//

// -------------------------------------------------------------------
//
// Method IsConstrainedEdge()
//
// It returns true  if a given edge is  constrained, and returns false
// if it is not.
//
// -------------------------------------------------------------------
bool
tCQMMesh2Builder::isConstrainedEdge(
				    const tCQMEdgeKey& ce,
				    const ConstrainedEdgeSet& le
				    )
                                    const
{
  int cv1 = ce.getEdgeC1();
  int cv2 = ce.getEdgeC2();

  ConstrainedEdgeSet::const_iterator epos;
  epos = le.find(cv1);
  while ((epos != le.end()) && (epos->first == cv1)) {
    if (epos->second == cv2) return true;
    epos++;
  }

  epos = le.find(cv2);
  while ((epos != le.end()) && (epos->first == cv2)) {
    if (epos->second == cv1) return true;
    epos++;
  }

  return false;
}

//
// modifiers
//

// -------------------------------------------------------------------
//
// Method insertGlue()
//
// It  inserts  two half-edges  into  the  list  of half-edges  to  be
// identified  when the  faces containing  them are  attached  to each
// other.
//
// -------------------------------------------------------------------
void
tCQMMesh2Builder::insertGlue(
			     tCQMHalfEdge2* he1, 
			     tCQMHalfEdge2* he2, 
			     std::list<tCQMGlueEdge2>& glist
			     ) 
{
  tCQMGlueEdge2 ge;
  ge.m_he1 = he1;
  ge.m_he2 = he2;
  glist.push_back(ge);
}

// -------------------------------------------------------------------
//
// Method makeTriFace()
//
// It creates  a new triangular face  as a new  connected component of
// the mesh.
//
// -------------------------------------------------------------------
void
tCQMMesh2Builder::makeTriFace(
			      const tCQMTriangle2& f,
			      int cv1,
			      int cv2,
			      int cv3,
			      const Point2& p1,
			      const Point2& p2,
			      const Point2& p3,
			      const tCQMEdgeKey& ce1,
			      const tCQMEdgeKey& ce2,
			      const tCQMEdgeKey& ce3,
			      const ConstrainedEdgeSet& le
			      )
{
  // create all three vertices
  tCQMVertex2* v1 = m_cp->createVertex(p1);
  tCQMVertex2* v2 = m_cp->createVertex(p2);
  tCQMVertex2* v3 = m_cp->createVertex(p3);

  // make the face and keep track of its edges
  tCQMFace2* newf = m_cp->mkFBEV(v1);
  tCQMHalfEdge2* he1 = newf->getHalfEdge();
  tCQMHalfEdge2* he2 = (m_cp->mkEV(v2,he1))->getHalfEdge();
  tCQMHalfEdge2* he3 = (m_cp->mkEV(v3,he2))->getHalfEdge();

  // insert new vertice codes into boundary vertex set
  vs.insert(make_pair(cv1,v1));
  vs.insert(make_pair(cv2,v2));
  vs.insert(make_pair(cv3,v3));

  // insert new edge codes into boundary edge set
  es.insert(make_pair(ce1.reverse(),he1));
  es.insert(make_pair(ce2.reverse(),he2));
  es.insert(make_pair(ce3.reverse(),he3));

  // set status of constrained edges
  if (isConstrainedEdge(ce1,le)) he1->getEdge()->constraintOn();
  if (isConstrainedEdge(ce2,le)) he2->getEdge()->constraintOn();
  if (isConstrainedEdge(ce3,le)) he3->getEdge()->constraintOn();
}

// -------------------------------------------------------------------
//
// Method makeTriFace()
//
// It creates a new triangular face and attaches it to another face in
// the mesh.
//
// -------------------------------------------------------------------
void
tCQMMesh2Builder::makeTriFace(
			      const tCQMTriangle2& f,
			      int cv,
			      const Point2& p,
			      const tCQMEdgeKey& ce,
			      const tCQMEdgeKey& ce2,
			      const tCQMEdgeKey& ce3,
			      const ConstrainedEdgeSet& le
			      )
{
  // create all three vertices
  tCQMVertex2* v = m_cp->createVertex(p);

  // get first half-edge of the commom edge
  tCQMHalfEdge2* he = (es.find(ce))->second;

  // make a new face with a commom edge and a new edge
  tCQMEdge2* e = ((m_cp->mkFE(he))->getHalfEdge())->getEdge();

  // get first half-edge of the new edge
  tCQMHalfEdge2* he2 = e->getHalfEdge();

  // create another edge and get its first half-edge
  tCQMHalfEdge2* he3 = (m_cp->mkEV(v,he2))->getHalfEdge();

  // insert new vertice code into boundary vertex set
  vs.insert(make_pair(cv,v));

  // insert new edge codes into boundary edge set
  es.insert(make_pair(ce2.reverse(),he2));
  es.insert(make_pair(ce3.reverse(),he3));

  // remove commom edge from the boundary edge set
  es.erase(es.find(ce));

  // set status of constrained edges
  if (isConstrainedEdge(ce2,le)) he2->getEdge()->constraintOn();
  if (isConstrainedEdge(ce3,le)) he3->getEdge()->constraintOn();
}

// -------------------------------------------------------------------
//
// Method makeTriFaceGlue1()
//
// It creates  a new triangular face  as a new  connected component of
// the mesh  and keeps track of  its common edge with  another face of
// the mesh, so that we can glue this connected component later.
//
// -------------------------------------------------------------------
void
tCQMMesh2Builder:: makeTriFaceGlue1(
				    const tCQMTriangle2& f,
				    int cv1,
				    int cv2,
				    int cv3,
				    const Point2& p1,
				    const Point2& p2,
				    const Point2& p3,
				    const tCQMEdgeKey& ce,
				    const tCQMEdgeKey& ce2,
				    const tCQMEdgeKey& ce3,
				    std::list<tCQMGlueEdge2>& glist,
				    const ConstrainedEdgeSet& le
				    )
{
  // create all three vertices
  tCQMVertex2* v1 = m_cp->createVertex(p1);
  tCQMVertex2* v2 = m_cp->createVertex(p2);
  tCQMVertex2* v3 = m_cp->createVertex(p3);

  // make the face and keep track of its edges
  tCQMFace2* newf = m_cp->mkFBEV(v1);
  tCQMHalfEdge2* he1 = newf->getHalfEdge();
  tCQMHalfEdge2* he2 = (m_cp->mkEV(v2,he1))->getHalfEdge();
  tCQMHalfEdge2* he3 = (m_cp->mkEV(v3,he2))->getHalfEdge();

  // insert new vertice codes into boundary vertex set
  vs.insert(make_pair(cv1,v1));
  vs.insert(make_pair(cv2,v2));
  vs.insert(make_pair(cv3,v3));

  // insert new edge codes into boundary edge set
  es.insert(make_pair(ce2.reverse(),he2));
  es.insert(make_pair(ce3.reverse(),he3));

  // keep track of edges to be glued at second stage
  EdgeSet::iterator hepos = es.find(ce);
  insertGlue(hepos->second,he1,glist);
  
  // remove commom edge from the boundary edge set
  es.erase(hepos);

  // set status of constrained edges
  if (isConstrainedEdge(ce,le)) he1->getEdge()->constraintOn();
  if (isConstrainedEdge(ce2,le)) he2->getEdge()->constraintOn();
  if (isConstrainedEdge(ce3,le)) he3->getEdge()->constraintOn();
}

// -------------------------------------------------------------------
//
// Method makeTriFaceGlue2()
//
// It creates  a new triangular face  as a new  connected component of
// the mesh and  keeps track of its two common  edges with other faces
// of the mesh, so that we can glue this connected component later.
//
// -------------------------------------------------------------------
void
tCQMMesh2Builder:: makeTriFaceGlue2(
				    const tCQMTriangle2& f,
				    int cv1,
				    int cv2,
				    int cv3,
				    const Point2& p1,
				    const Point2& p2,
				    const Point2& p3,
				    const tCQMEdgeKey& ce1,
				    const tCQMEdgeKey& ce2,
				    const tCQMEdgeKey& ce3,
				    int ne,
				    std::list<tCQMGlueEdge2>& glist,
				    const ConstrainedEdgeSet& le
				    )
{
  // create all three vertices
  tCQMVertex2* v1 = m_cp->createVertex(p1);
  tCQMVertex2* v2 = m_cp->createVertex(p2);
  tCQMVertex2* v3 = m_cp->createVertex(p3);

  // make the face and keep track of its edges
  tCQMFace2* newf = m_cp->mkFBEV(v1);
  tCQMHalfEdge2* he1 = newf->getHalfEdge();
  tCQMHalfEdge2* he2 = (m_cp->mkEV(v2,he1))->getHalfEdge();
  tCQMHalfEdge2* he3 = (m_cp->mkEV(v3,he2))->getHalfEdge();

  // insert new vertice codes into boundary vertex set
  vs.insert(make_pair(cv1,v1));
  vs.insert(make_pair(cv2,v2));
  vs.insert(make_pair(cv3,v3));

  // insert new edge code into boundary edge set and
  // keep track of the two edges to be glued in the
  // second stage.
  EdgeSet::iterator hepos1, hepos2;
  if (ne == 1) {
    hepos1 = es.find(ce2);
    hepos2 = es.find(ce3);
    es.insert(make_pair(ce1.reverse(),he1));
    insertGlue(hepos1->second,he2,glist);
    insertGlue(hepos2->second,he3,glist);
  }
  else if (ne == 2) {
    hepos1 = es.find(ce1);
    hepos2 = es.find(ce3);
    es.insert(make_pair(ce2.reverse(),he2));
    insertGlue(hepos1->second,he1,glist);
    insertGlue(hepos2->second,he3,glist);
  }
  else {
    hepos1 = es.find(ce1);
    hepos2 = es.find(ce2);
    es.insert(make_pair(ce3.reverse(),he3));
    insertGlue(hepos1->second,he1,glist);
    insertGlue(hepos2->second,he2,glist);
  }

  // remove commom edges from the boundary edge set
  es.erase(hepos1);
  es.erase(hepos2);

  // set status of constrained edges
  if (isConstrainedEdge(ce1,le)) he1->getEdge()->constraintOn();
  if (isConstrainedEdge(ce2,le)) he2->getEdge()->constraintOn();
  if (isConstrainedEdge(ce3,le)) he3->getEdge()->constraintOn();
}

// -------------------------------------------------------------------
//
// Method makeTriFaceGlue3()
//
// It creates  a new triangular face  as a new  connected component of
// the mesh and keeps track of its three common edges with other faces
// of the mesh, so that we can glue this connected component later.
//
// -------------------------------------------------------------------
void
tCQMMesh2Builder:: makeTriFaceGlue3(
				    const tCQMTriangle2& f,
				    int cv1,
				    int cv2,
				    int cv3,
				    const Point2& p1,
				    const Point2& p2,
				    const Point2& p3,
				    const tCQMEdgeKey& ce1,
				    const tCQMEdgeKey& ce2,
				    const tCQMEdgeKey& ce3,
				    std::list<tCQMGlueEdge2>& glist,
				    const ConstrainedEdgeSet& le
				    )
{
  // create all three vertices
  tCQMVertex2* v1 = m_cp->createVertex(p1);
  tCQMVertex2* v2 = m_cp->createVertex(p2);
  tCQMVertex2* v3 = m_cp->createVertex(p3);

  // make the face and keep track of its edges
  tCQMFace2* newf = m_cp->mkFBEV(v1);
  tCQMHalfEdge2* he1 = newf->getHalfEdge();
  tCQMHalfEdge2* he2 = (m_cp->mkEV(v2,he1))->getHalfEdge();
  tCQMHalfEdge2* he3 = (m_cp->mkEV(v3,he2))->getHalfEdge();

  // insert new vertice codes into boundary vertex set
  vs.insert(make_pair(cv1,v1));
  vs.insert(make_pair(cv2,v2));
  vs.insert(make_pair(cv3,v3));

  // insert new edge code into boundary edge set
  EdgeSet::iterator hepos1, hepos2, hepos3;
 
  // keep track of edges to be glued in the second stage
  hepos1 = es.find(ce1);
  hepos2 = es.find(ce2);
  hepos3 = es.find(ce3);

  insertGlue(hepos1->second,he1,glist);
  insertGlue(hepos2->second,he2,glist);
  insertGlue(hepos3->second,he3,glist);

  // remove commom edges from the boundary edge set
  es.erase(hepos1);
  es.erase(hepos2);
  es.erase(hepos3);

  // set status of constrained edges
  if (isConstrainedEdge(ce1,le)) he1->getEdge()->constraintOn();
  if (isConstrainedEdge(ce2,le)) he2->getEdge()->constraintOn();
  if (isConstrainedEdge(ce3,le)) he3->getEdge()->constraintOn();
}

// -------------------------------------------------------------------
//
// Method glueFaces()
//
// It attaches one  face to another along a common  edge. If the edges
// are in the same boundary, the glueing might create a genus. If they
// belong to distinct connected components, the componentes are merged
// together.
//
// -------------------------------------------------------------------
void
tCQMMesh2Builder::glueFaces(
			    tCQMHalfEdge2* he1,
			    tCQMHalfEdge2* he2
			    )
{
  // get the boundaries containing the half-edges
  tCQMBoundary2* bd1 = he1->getEdge()->getBoundary();
  tCQMBoundary2* bd2 = he2->getEdge()->getBoundary();

  // if the boundaries are distinct, we merge them.
  // otherwise, we will make a genus.
  if (bd1 != bd2) {
    m_cp->rmBEV(he1,he2);
  }
  else {
    m_cp->mkGrmEV(he1,he2);
  }
}

// -------------------------------------------------------------------
//
// Method buildTriMesh()
//
// It  stores  a  triangular  mesh  into  the  data  structure  handle
// edge. All  it needs is the  collection of triangles of  the mesh in
// any  given  order.   Triangles  are  expected  to  be  consistently
// oriented.
//
// -------------------------------------------------------------------
tCQMComplex2*
tCQMMesh2Builder::buildTriMesh(const string& name)
{
  // read files with information about vertex and triangles
  string infile1 = string(name) + string(".node");
  string infile2 = string(name) + string(".edge");
  string infile3 = string(name) + string(".ele");

  // create file reader
  tCQMMeshReader reader(infile1, infile2, infile3);

  // read files
  Point2* lv;
  ConstrainedEdgeSet le;
  std::list<tCQMTriangle2> lf;

  try {
    printf("\nstart reading vertex file\n");
    reader.readVertexFile(lv);
    printf("start reading edge file\n");
    reader.readEdgeFile(le);
    printf("start reading elem file\n");
    reader.readElemFile(lv,lf);
  }
  catch (const tExceptionObject& xpt) {
    treatException(xpt);
    exit(0);
  }

  delete [] lv;

  typedef std::list<tCQMTriangle2>::iterator LFI;

  std::list<tCQMGlueEdge2> gluelist;

  // first stage: making faces
 
  for (LFI fiter = lf.begin(); fiter != lf.end(); ++fiter) {
    // get face from iterator
    tCQMTriangle2 f(*fiter);

    // check orientation
    if (!f.isCCW()) {
      f.reverse();
    }

    // create vertex codes for every vertex of face f
    int cv1 = f.getVertexId(0);
    int cv2 = f.getVertexId(1);
    int cv3 = f.getVertexId(2);

    // get vertex locations
    Point2 lv1 = f.getPoint(0);
    Point2 lv2 = f.getPoint(1);
    Point2 lv3 = f.getPoint(2);

    // create edge codes for every edge of face f
    tCQMEdgeKey ce1 = tCQMEdgeKey(cv1,cv2);
    tCQMEdgeKey ce2 = tCQMEdgeKey(cv2,cv3);
    tCQMEdgeKey ce3 = tCQMEdgeKey(cv3,cv1);

    // search for existing edges in the boundary edge set
    int foundedge = 0;
    if (es.find(ce1) != es.end()) foundedge += 1;
    if (es.find(ce2) != es.end()) foundedge += 2;
    if (es.find(ce3) != es.end()) foundedge += 4;

    // if only one of the edges was found, new face should
    // be created by attaching commom edge if it does not
    // create a singular vertex.
    if ((foundedge == 1) || (foundedge == 2) || (foundedge == 4)) {
      int cvaux;
      Point2 lvaux;
      tCQMEdgeKey ceaux1, ceaux2, ceaux3;

      if (foundedge == 1) {
	cvaux = cv3;
	lvaux = lv3;
	ceaux1 = ce2;
        ceaux2 = ce3; 
	ceaux3 = ce1;
      }
      else if (foundedge == 2) {
	cvaux = cv1;
	lvaux = lv1;
	ceaux1 = ce3;
	ceaux2 = ce1;
	ceaux3 = ce2;
      }
      else {  // foundedge == 4
	cvaux = cv2;
	lvaux = lv2;
	ceaux1 = ce1;
	ceaux2 = ce2;
	ceaux3 = ce3;
      }

      // if vertex cvaux is not in the same boundary as edge
      // ceaux3.reverse(), we can attach face f to the boundary
      // containing ceaux3.reverse().
      // if (vs.find(cvaux) == vs.end()) {
      if (!findVertexInBound(cvaux,ceaux3)) {
	// either cvaux has not been created yet or it is not
	// on the boundary of ceaux3.reverse().

	this->makeTriFace(f,cvaux,lvaux,ceaux3,ceaux1,ceaux2,le);
      }
      else {
	// create an isolated face, which should be attached to f
	// in the second stage.

        this->makeTriFaceGlue1(f,cv1,cv2,cv3,lv1,lv2,lv3,ceaux3,
			       ceaux1,ceaux2,gluelist,le);
      }
    }
    else if (foundedge != 0) {   // foundedge == 3, 5, 6 or 7

      // if new face has two common edges with existent faces,
      // we must attach it to one of the faces, and glue it to
      // the other in the second stage.

      if (foundedge == 3) {
	// edges ce1 and ce2 are common edges
	//
	this->makeTriFaceGlue2(f,cv1,cv2,cv3,lv1,lv2,lv3,ce1,ce2,
			       ce3,3,gluelist,le);
      }
      else if (foundedge == 5) {
	// edges ce1 and ce3 are common edges
	//
	this->makeTriFaceGlue2(f,cv1,cv2,cv3,lv1,lv2,lv3,ce1,ce2,
			       ce3,2,gluelist,le);
      }
      else if (foundedge == 6) {
	// edges ce2 and ce3 are common edges
	//
	this->makeTriFaceGlue2(f,cv1,cv2,cv3,lv1,lv2,lv3,ce1,ce2,
			       ce3,1,gluelist,le);
      }
      else {
	// all edges in common
	//
	this->makeTriFaceGlue3(f,cv1,cv2,cv3,lv1,lv2,lv3,ce1,ce2,
			       ce3,gluelist,le);
      }
    }
    else { // foundedge == 0

      // new face has no vertex in common with existent faces
      // or one, two or three vertices in common with existent
      // faces.

      // we will create isolated faces in any case

      this->makeTriFace(f,cv1,cv2,cv3,lv1,lv2,lv3,ce1,ce2,ce3,le);
    }
  }

  // removes all elements from vertex and edge boundary sets
  vs.clear();
  es.clear();

  // remove all elements from constrained edge set

  // second stage: glueing faces
  typedef std::list<tCQMGlueEdge2>::iterator LGI;
  LGI giter = gluelist.begin();
  while (giter != gluelist.end()) {
    // get the half-edges to be identified in the faces
    // to be glued.
    tCQMHalfEdge2* he1 = (*giter).m_he1;
    tCQMHalfEdge2* he2 = (*giter).m_he2;

    // glue faces
    glueFaces(he1,he2);

    // remove the half-edges from the glue list
    giter = gluelist.erase(giter);
  }

  // return complex
  return m_cp;
}
