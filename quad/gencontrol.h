//  ===================================================================
//
//  Program:   Convex Quadrilateral Mesh Generator (CQMesh)
//  Module:    gencontrol.h
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

#ifndef __gencontrol_h
#define __gencontrol_h

#include <iostream>
#include <string>
#include <exception>
#include <list>
#include <map>

#include "typedefs.h"
#include "err.h"

#include "edgekey.h"
#include "vertex2.h"
#include "halfedge2.h"
#include "edge2.h"
#include "complex2.h"
#include "mesh2.h"

#include "tri2.h"
#include "quad2.h"
#include "compquad.h"
#include "iomesh2.h"

// -------------------------------------------------------------------
// tCQMGeneralControl
// ------------------
class tCQMGeneralControl {
   private:
     // file  name of  the input  files describing  a  triangular mesh
     // (*.node, *.edge, *.ele files).
     std::string m_fname;

     // pointer to the triangular mesh data structure
     tCQMComplex2* m_mesh;

     // list   of  quadrilaterals   resulting   from  converting   the
     // triangulation into a convex quadrangulation.                                       
     std::list<tCQMQuadrilateral2*>* m_lcquads;

     // list  of the edges  of all  quadrilaterals resulting  from the
     // convex quadrangulation.
     std::map<tCQMEdgeKey,bool>* m_ledges;

     // list of the vertices  of all quadrilaterals resulting from the
     // convex quadrangulation.
     std::map<int,Point2>* m_lverts;

     // number of vertices in the triangulation
     int m_triverts;

     // number of edges in the triangulation
     int m_triedges;

     // number of triangles in the triangulation
     int m_trieles;

     // number of Steiner points in the quadrangulation
     int m_spoints;

     // private methods
     void getListOfEdges();
     void getListOfVerts();

   public:
     // constructors
     tCQMGeneralControl(string);

     // destructor
     ~tCQMGeneralControl();

     // get methods
     tCQMComplex2* getTriMesh();

     // modifiers
     void computeTriangulation();
     void computeConvexQuadrangulation();

     // output mesh
     void saveQuadrangulation();

     // output statistics
     void generateStatistics() throw (tExceptionObject);
};

// -------------------------------------------------------------------
// tCQMGeneralControl - inline implementation
// ------------------
inline tCQMComplex2*
tCQMGeneralControl::getTriMesh()
{
  return m_mesh;
}

#endif   // __gencontrol_h
