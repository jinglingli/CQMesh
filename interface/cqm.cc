//  ===================================================================
//
//  Program:   Convex Quadrilateral Mesh Generator (CQMesh)
//  Module:    cqm.cc
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

#include "cqm.h"

int main(int argc, char** argv)
{
  // check whether  user provided  a filename for  the files  with the
  // vertices and elements of a mesh.
  if (argc != 2) {
    std::cout << "Usage: cqm <mesh filename>" << std::endl;
    exit(0);
  }

  // create an object to control the process of reading the input mesh
  // and creating the quadrangulation.
  tCQMGeneralControl* control = (tCQMGeneralControl*) 
    new tCQMGeneralControl(std::string(argv[1]));

  // read triangulation
  std::cout << std::endl << "reading triangulation..." << std::endl;
  control->computeTriangulation();

  // compute quadrangulation
  std::cout << "computing quadrangulation..." << std::endl;
  control->computeConvexQuadrangulation();

  // saving quadrangulation
  std::cout << "writing quadrangulation to a file..." << std::endl;
  control->saveQuadrangulation();

  // generating statistics
  std::cout << "generating statistics about quadrangulation..." << std::endl;
  try {
    control->generateStatistics();
  }
  catch (const tExceptionObject& xpt) {
    treatException(xpt);
    exit(0);
  }

  delete control;

  std::cout << "done!" << std::endl;
}
