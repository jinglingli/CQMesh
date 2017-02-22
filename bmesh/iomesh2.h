//  ===================================================================
//
//  Program:   Convex Quadrilateral Mesh Generator (CQMesh)
//  Module:    iomesh2.h
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

#ifndef __iomesh2_h
#define __iomesh2_h

#include <iostream>
#include <string>
#include <list>
#include <map>
#include <exception>

#include "typedefs.h"
#include "err.h"
#include "edgekey.h"
#include "tri2.h"
#include "quad2.h"

typedef std::multimap<int,int> ConstrainedEdgeSet;

// -------------------------------------------------------------------
// Class tCQMLexer: a lexer for "node" and "ele" mesh files
// ---------------

class tCQMLexer {
  private:
    std::istream& real_stream_;
    int linecount_;

    // private methods
    void throw_error(const std::string& s) throw (tExceptionObject);
    void increment_line_();
    void skip_space_();

  public:
    // constructor
    tCQMLexer (std::istream& is, int init_linecount = -1) : 
      real_stream_(is), 
      linecount_(init_linecount) 
    {
      increment_line_(); 
    }

    // public methods
    bool get_int(int& i);
    bool get_double(double& d);
    int linecount() const {return linecount_;}
};

// -------------------------------------------------------------------
// tCQMMeshReader: a reader for "node", "edge" and "ele" mesh file
// --------------
class tCQMMeshReader {
  private:
    std::string filename1;
    std::string filename2;
    std::string filename3;
    int linecount_;
  
   public:
    // constructor
    tCQMMeshReader(const std::string& fname1,
		   const std::string& fname2,
		   const std::string& fname3) :
      filename1(fname1),
      filename2(fname2),
      filename3(fname3),
      linecount_(0)
    {}

    int linecount() {return linecount_;}

    void readVertexFile(Point2*&) throw (tExceptionObject);
    void readEdgeFile(ConstrainedEdgeSet&) throw (tExceptionObject); 
    void readElemFile(Point2*, std::list<tCQMTriangle2>&) 
      throw (tExceptionObject);
};

// -------------------------------------------------------------------
// tCQMMeshWriter: a writer for "node" and "ele" mesh file
// --------------
class tCQMMeshWriter {
  private:
    std::string filename1;
    std::string filename2;
    std::string filename3;
    std::string filename4;

  public:
    // constructor
    tCQMMeshWriter(const std::string& fname1, 
		   const std::string& fname2, 
		   const std::string& fname3, 
		   const std::string& fname4) :
      filename1(fname1),
      filename2(fname2),
      filename3(fname3),
      filename4(fname4)
    {}

    void writeVertexFile(const std::map<int,Point2>&) 
      throw (tExceptionObject);
    void writeEdgeFile(const std::map<tCQMEdgeKey,bool>&) 
      throw (tExceptionObject);
    void writeElemFile(const std::list<tCQMQuadrilateral2*>&) 
      throw (tExceptionObject);
    void writePolyFile(const std::map<int,Point2>&,
		       const std::map<tCQMEdgeKey,bool>&) 
      throw (tExceptionObject);
};

#endif    // #define __iomesh2_h
