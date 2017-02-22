//  ===================================================================
//
//  Program:   Convex Quadrilateral Mesh Generator (CQMesh)
//  Module:    gencontrol.cc
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

#include <sstream>
#include <fstream>
#include <cctype>
#include <sys/time.h>

#include "gencontrol.h"

// -------------------------------------------------------------------
// tCQMGeneralControl - implementation
// ------------------

// constructors
tCQMGeneralControl::tCQMGeneralControl(string fn)
{
  m_fname = fn;
  m_mesh = 0;
  m_lcquads = 0;
  m_ledges = 0;
  m_lverts = 0;
  m_triverts = 0;
  m_triedges = 0;
  m_trieles = 0;
  m_spoints = 0;
}

// destructor
tCQMGeneralControl::~tCQMGeneralControl()
{
  if (m_mesh != 0) delete m_mesh;

  if (m_lcquads != 0) {
    while (!m_lcquads->empty()) {
      tCQMQuadrilateral2* q = m_lcquads->front();
       m_lcquads->pop_front();
       delete q;
    }

    delete m_lcquads;
  }

  if (m_ledges != 0) {
    m_ledges->clear();
    delete m_ledges;
  }

  if (m_lverts != 0) {
    m_lverts->clear();
    delete m_lverts;
  }
}

// private methods
void
tCQMGeneralControl::getListOfVerts()
{
  if ((m_lcquads == 0) || (m_lverts != 0)) return;

  m_lverts = (std::map<int,Point2>*) new  std::map<int,Point2>();

  std::list<tCQMQuadrilateral2*>::iterator li = m_lcquads->begin();
  for (; li != m_lcquads->end(); ++li) {
    tCQMQuadrilateral2* p = *li;
    for (int i=0; i<p->getNumVerts(); i++) {
      int vc = p->getVertexId(i);
      if (m_lverts->find(vc) == m_lverts->end()) { 
	m_lverts->insert(make_pair(vc,p->getPoint(i)));
      }
    }
  }
}

void
tCQMGeneralControl::getListOfEdges()
{
  if ((m_lcquads == 0) || (m_ledges != 0)) return;

  m_ledges = (std::map<tCQMEdgeKey,bool>*) 
    new std::map<tCQMEdgeKey,bool>();

  std::list<tCQMQuadrilateral2*>::iterator li = m_lcquads->begin();
  for (; li != m_lcquads->end(); ++li) {
    tCQMQuadrilateral2* p = *li;
    for (int i=0; i<p->getNumVerts(); i++) {
      int vc1 = p->getVertexId(i);
      int vc2 = p->getVertexId((i + 1) % p->getNumVerts());
      tCQMEdgeKey ec(vc2,vc1);
      if (m_ledges->find(ec) == m_ledges->end()) { 
	m_ledges->insert(make_pair(ec.reverse(),p->isConstrained(i)));
      }
    }
  }
}

// modifiers
void
tCQMGeneralControl::computeTriangulation()
{
  if (m_mesh == 0) {
    m_mesh = tCQMMesh2Builder().buildTriMesh(m_fname);

    m_triverts = m_mesh->getNumVerts();
    m_triedges = m_mesh->getNumEdges();
    m_trieles = m_mesh->getNumFaces();
  }
}

void
tCQMGeneralControl::computeConvexQuadrangulation()
{
  if (m_mesh == 0) return;

  if (m_lcquads == 0) {
    tCQMCompQuad compquad(m_mesh);

    struct timeval stime, etime;

    gettimeofday(&stime,0);
    m_lcquads = compquad.computeConvexQuadrangulation();
    gettimeofday(&etime,0);

    long tt = 1000l * (etime.tv_sec - stime.tv_sec) 
      + (etime.tv_usec - stime.tv_usec) / 1000l;

    std::cout << "Quadrangulation generated in " << tt 
	      << " miliseconds." << std::endl;

    m_spoints = compquad.getNumOfSteinerPoints();
    getListOfVerts();
    getListOfEdges();
  }
}

void
tCQMGeneralControl::saveQuadrangulation()
{
  if (m_lcquads == 0) return;

  if (m_lcquads->empty()) return;


  string fn1 = string(m_fname) + string("-quad.node");
  string fn2 = string(m_fname) + string("-quad.edge");
  string fn3 = string(m_fname) + string("-quad.ele");
  string fn4 = string(m_fname) + string("-quad.poly");

  tCQMMeshWriter mw(fn1,fn2,fn3,fn4);

  try {
    mw.writeVertexFile(*m_lverts);
    mw.writeEdgeFile(*m_ledges);
    mw.writeElemFile(*m_lcquads);
    mw.writePolyFile(*m_lverts,*m_ledges);
  }
  catch (const tExceptionObject& xpt) {
    treatException(xpt);
    exit(0);
  }
}

void
tCQMGeneralControl::generateStatistics() throw (tExceptionObject)
{
  if (m_lcquads == 0) return;

  if (m_lcquads->empty()) return;

  // compute minimum, maximum, and average area of quadrilaterals
  double minarea = 0;
  double maxarea = 0;
  double avgarea = 0;

  std::list<tCQMQuadrilateral2*>::iterator li = m_lcquads->begin();
  if (li != m_lcquads->end()) {
    tCQMQuadrilateral2* p = *li;
    double area = p->area();
    minarea = maxarea = avgarea = area;

    if (area < 0) {
      std::cout << "area: " << area << std::endl; 
      std::cout << "(" << p->getPoint(0).x() << "," 
		<< p->getPoint(0).y() << "), "
		<< "(" << p->getPoint(1).x() 
		<< "," << p->getPoint(1).y() << "), "
		<< "(" << p->getPoint(2).x() << "," 
		<< p->getPoint(2).y() << "), "
		<< "(" << p->getPoint(3).x() << "," 
		<< p->getPoint(3).y() << ")"
		<< std::endl;

      if (!p->isSimple()) {
	std::cout << "polygon is not simple!" << std::endl;
      }
      else if (p->isCW()) {
	std::cout << "polygon is CW-oriented" << std::endl;
      }

      area = -area;
    }

    ++li;
    for (; li != m_lcquads->end(); ++li) {
      p = *li;
      area = p->area();

      if (area < 0) {
	std::cout << "area: " << area << std::endl;
	std::cout << "(" << p->getPoint(0).x() << "," 
		  << p->getPoint(0).y() << "), "
		  << "(" << p->getPoint(1).x() << "," 
		  << p->getPoint(1).y() << "), "
		  << "(" << p->getPoint(2).x() << "," 
		  << p->getPoint(2).y() << "), "
		  << "(" << p->getPoint(3).x() << "," 
		  << p->getPoint(3).y() << ")"
		  << std::endl;

	if (!p->isSimple()) {
	  std::cout << "polygon is not simple!" << std::endl;
	}
	else if (p->isCW()) {
	  std::cout << "polygon is CW-oriented" << std::endl;
	}

	area = -area;
      }

      if (area > maxarea) maxarea = area;
      if (area < minarea) minarea = area;
      avgarea += area;
    }

    avgarea /= m_lcquads->size();
  }
  else {
    return;
  }

  // create output file
  std::string fn = std::string(m_fname) + std::string(".stat");

  // write statistics
  std::filebuf fb;
  if (fb.open(fn.c_str(), std::ios::out) == 0) {
    std::stringstream ss (std::stringstream::in | 
			  std::stringstream::out);
    ss << "Unable to create output file";
    fb.close();
    throw tExceptionObject(__FILE__,__LINE__,ss.str());
  }

  std::ostream ostr(&fb);

  ostr << "Triangulation" << std::endl;
  ostr << "-------------" << std::endl << std::endl;
  ostr << "Number of triangles: " << m_trieles  << std::endl;
  ostr << "Number of vertices : " << m_triverts << std::endl;
  ostr << "Number of edges    : " << m_triedges << std::endl;
  ostr << std::endl << std::endl << std::endl;
  ostr << "Quadrangulation" << std::endl;
  ostr << "---------------" << std::endl << std::endl;
  ostr << "Number of quadrilaterals  : " 
	    << m_lcquads->size() << std::endl;
  ostr << "Number of vertices        : " 
	    << m_lverts->size() << std::endl;
  ostr << "Number of edges           : " 
	    << m_ledges->size() << std::endl;
  ostr << "Number of Steiner points  : " 
	    << m_spoints << std::endl;
  ostr << "Maximum quadrilateral area: " 
	    << maxarea << std::endl;
  ostr << "Minimum quadrilateral area: " 
	    << minarea << std::endl;
  ostr << "Average quadrilateral area: " 
	    << avgarea << std::endl;
  ostr << std::endl << std::endl;

  fb.close();
}
