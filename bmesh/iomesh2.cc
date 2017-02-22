//  ===================================================================
//
//  Program:   Convex Quadrilateral Mesh Generator (CQMesh)
//  Module:    iomesh2.cc
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

#include "iomesh2.h"

// -------------------------------------------------------------------
// Class tCQMLexer - implementation
// -------------------------------------------------------------------

//
// private methods
//

void
tCQMLexer::throw_error(const std::string& s) throw (tExceptionObject) {
  std::stringstream ss (std::stringstream::in | 
			std::stringstream::out);
  ss << "File line " << linecount_ << ": " << s;
  throw tExceptionObject(__FILE__,__LINE__,ss.str());
}

void
tCQMLexer::increment_line_() {
  if (linecount_ >= 0) {
    ++linecount_;
  }
}

void
tCQMLexer::skip_space_() {
  char c;
  for (;;) {
    if (!real_stream_.get(c)) {
      return;
    }

    if (c == '\n') {
      increment_line_();
    }

    if (!isspace(c)) {
      real_stream_.putback(c);
      return;
    }
  }
}


//
// public methods
//

bool
tCQMLexer::get_int(int& i) {
  skip_space_();
  char c;
  if (!real_stream_.get(c)) {
    throw_error("Unexpected eof in the middle of reading an int");
  }

  real_stream_.putback(c);
  if (!(real_stream_ >> i)) { 
    throw_error("Error while reading integer");
  }

  return true;
}
  
bool
tCQMLexer::get_double(double& d) {
  skip_space_();
  char c;
  if (!real_stream_.get(c)) {
    throw_error("Unexpected eof in the middle of reading a double");
  }
   
  real_stream_.putback(c);
  if (!(real_stream_ >> d)) {
    throw_error("Error while reading double");
  }

  return true;
}


// -----------------------------------------------------
// Class tCQMMeshReader - implementation
// -----------------------------------------------------

//
// public methods
//
void
tCQMMeshReader::readVertexFile(Point2*& lverts) throw (tExceptionObject) 
{
  // opens the file
  std::filebuf fb;
  if (fb.open(filename1.c_str(), std::ios::in) == 0) {
    std::stringstream ss (std::stringstream::in | 
			  std::stringstream::out);
    ss << "Unable to open input file\n";
    ss << filename1.c_str();
    fb.close();
    throw tExceptionObject(__FILE__,__LINE__,ss.str());
  }

  std::istream istr(&fb);
  tCQMLexer is(istr, linecount_);

  // reads first line:
  //
  // <# of vertices><dimension><# of attributes><boundary markers (0 or 1)> 
  //
  int nVerts;
  int dim;
  int nAttr;       // it has not been used so far
  int bMarker;     // it has not been used so far

  is.get_int(nVerts);
  if (nVerts < 3) {
    std::stringstream ss (std::stringstream::in | 
			  std::stringstream::out);
    ss << "Line " << linecount_ << ": " 
       << filename1.c_str() << " contains " << nVerts << " points";
    fb.close();
    throw tExceptionObject(__FILE__,__LINE__,ss.str());
  }

  is.get_int(dim);
  if (dim != 2) {
    std::stringstream ss (std::stringstream::in | 
			  std::stringstream::out);
    ss << "Line " << linecount_ << ": " 
       << "I can only deal with dimension 2";
    fb.close();
    throw tExceptionObject(__FILE__,__LINE__,ss.str());
  }

  is.get_int(nAttr);
  if (nAttr != 0) {
    std::stringstream ss (std::stringstream::in | 
			  std::stringstream::out);
    ss << "Line " << linecount_ << ": " 
       << "Sorry, I do not support attributes yet.";
    fb.close();
    throw tExceptionObject(__FILE__,__LINE__,ss.str());
  }

  is.get_int(bMarker);
  if (bMarker < 0) {
    std::stringstream ss (std::stringstream::in | 
			  std::stringstream::out);
    ss << "Line " << linecount_ << ": " 
       << "Negative value for number of point markers";
    fb.close();
    throw tExceptionObject(__FILE__,__LINE__,ss.str());
  }

  if (bMarker > 1) {
    std::stringstream ss (std::stringstream::in | 
			  std::stringstream::out);
    ss << "Line " << linecount_ << ": " 
       << "Value greater than one for number of point markers";
    fb.close();
    throw tExceptionObject(__FILE__,__LINE__,ss.str());
  }

  // allocate two vectors of dimension nVerts each
  // for storing X,Y vertex coordinates
  lverts = new Point2[nVerts+1];
  double cx, cy;

  //
  // reads vertex coordinates
  //
  int pNumber = -1;
  int nNumber;
  for (int i=0; i<nVerts; i++) {

    // reads vertex number
    is.get_int(nNumber);

    // verifies if vertex numbering starts with 0 or 1
    if ((i == 0) && (nNumber == 1)) {
      pNumber = 0;
    }
    
    // verifies vertex number sequence
    if (nNumber != (pNumber + 1)) {
      std::stringstream ss (std::stringstream::in | 
			    std::stringstream::out);
      ss << "Line " << linecount_ << ": " 
	 << "Vertices are not numbered correctly";
      fb.close();
      throw tExceptionObject(__FILE__,__LINE__,ss.str());
    }

    // updates pNumber
    pNumber = nNumber;

    // reads X coordinate of vertex nNumber
    if (!is.get_double(cx)) {
      std::stringstream ss (std::stringstream::in | 
			    std::stringstream::out);
      ss << "Line " << linecount_ << ": " 
	 << "Problem reading X coordinate of vertex " << nNumber;
      fb.close();
      throw tExceptionObject(__FILE__,__LINE__,ss.str());
    }

    // reads Y coordinate of vertex nNumber
    if (!is.get_double(cy)) {
      std::stringstream ss (std::stringstream::in | 
			    std::stringstream::out);
      ss << "Line " << linecount_ << ": " 
	 << "Problem reading Y coordinate of vertex " << nNumber;
      fb.close();
      throw tExceptionObject(__FILE__,__LINE__,ss.str());
    }

    // create vertex with location and id information
    lverts[nNumber] = Point2(cx,cy);

    if (bMarker == 1) {
      int marker;
      if (!is.get_int(marker)) {
	std::stringstream ss (std::stringstream::in | 
			      std::stringstream::out);
	ss << "Line " << linecount_ << ": " 
	   << "Problem reading boundary marker";
	fb.close();
	throw tExceptionObject(__FILE__,__LINE__,ss.str());
      }
    }
  }

  // close file buffer
  fb.close();
}

void
tCQMMeshReader::readEdgeFile(ConstrainedEdgeSet& le) 
  throw (tExceptionObject) 
{
  // opens the file
  std::filebuf fb;
  if (fb.open(filename2.c_str(), std::ios::in) == 0) {
    std::stringstream ss (std::stringstream::in | 
			  std::stringstream::out);
    ss << "Unable to open input file";
    fb.close();
    throw tExceptionObject(__FILE__,__LINE__,ss.str());
  }

  std::istream istr(&fb);
  tCQMLexer is(istr, linecount_);

  // reads first line:
  //
  // <# of edges><boundary markers (0 or 1)> 
  //
  int nEdges;
  int bMarker;

  is.get_int(nEdges);
  if (nEdges < 3) {
    std::stringstream ss (std::stringstream::in | 
			  std::stringstream::out);
    ss << "Line " << linecount_ << ": " 
       << filename1.c_str() << " contains " << nEdges << " edges";
    fb.close();
    throw tExceptionObject(__FILE__,__LINE__,ss.str());
  }

  is.get_int(bMarker);
  if (bMarker < 0) {
    std::stringstream ss (std::stringstream::in | std::stringstream::out);
    ss << "Line " << linecount_ << ": " 
       << "Negative value for number of point markers";
    fb.close();
    throw tExceptionObject(__FILE__,__LINE__,ss.str());
  }

  if (bMarker > 1) {
    std::stringstream ss (std::stringstream::in | 
			  std::stringstream::out);
    ss << "Line " << linecount_ << ": " 
       << "Value greater than one for number of point markers";
    fb.close();
    throw tExceptionObject(__FILE__,__LINE__,ss.str());
  }

  //
  // read code of edge endpoints
  //
  int pNumber = -1;
  int nNumber;
  int cv1, cv2;
  for (int i=0; i<nEdges; i++) {

    // read edge number
    is.get_int(nNumber);

    // verifies if edge numbering starts with 0 or 1
    if ((i == 0) && (nNumber == 1)) {
      pNumber = 0;
    }
    
    // verifies edge number sequence
    if (nNumber != (pNumber + 1)) {
      std::stringstream ss (std::stringstream::in | 
			    std::stringstream::out);
      ss << "Line " << linecount_ << ": " 
         << "Edges are not numbered correctly";
      fb.close();
      throw tExceptionObject(__FILE__,__LINE__,ss.str());
    }

    // updates pNumber
    pNumber = nNumber;

    // read code of first endpoint of the edge
    if (!is.get_int(cv1)) {
      std::stringstream ss (std::stringstream::in | 
			    std::stringstream::out);
      ss << "Line " << linecount_ << ": " 
         << "Problem reading vertex code of edge first endpoint " << nNumber;
      fb.close();
      throw tExceptionObject(__FILE__,__LINE__,ss.str());
    }

    // read code of first endpoint of the edge
    if (!is.get_int(cv2)) {
      std::stringstream ss (std::stringstream::in | 
			    std::stringstream::out);
      ss << "Line " << linecount_ << ": " 
         << "Problem reading vertex code of edge second endpoint " << nNumber;
      fb.close();
      throw tExceptionObject(__FILE__,__LINE__,ss.str());
    }

    // reads marker
    int marker;
    if (bMarker == 1) {
      if (!is.get_int(marker)) {
	std::stringstream ss (std::stringstream::in | 
			      std::stringstream::out);
	ss << "Line " << linecount_ << ": " 
	   << "Problem reading first marker";
	fb.close();
	throw tExceptionObject(__FILE__,__LINE__,ss.str());
      }

      if (marker == 1) {
	le.insert(std::make_pair(cv1,cv2));
	le.insert(std::make_pair(cv2,cv1));
      }
    }
  }

  // close file buffer
  fb.close();
}

void
tCQMMeshReader::readElemFile(
			     Point2* lverts, 
			     std::list<tCQMTriangle2>& lfaces) 
  throw (tExceptionObject)
{
  // opens the file
  std::filebuf fb;
  if (fb.open(filename3.c_str(), std::ios::in) == 0) {
    std::stringstream ss (std::stringstream::in | 
			  std::stringstream::out);
    ss << "Unable to open input file";
    fb.close();
    throw tExceptionObject(__FILE__,__LINE__,ss.str());
  }

  std::istream istr(&fb);
  tCQMLexer is(istr, linecount_);

  // reads first line:
  //
  // <# of elements><# of vertices per element><# of attributes>
  //
  int nElems;
  int nVerts;      // it should match with argument nv
  int nAttrs;      // it has not been used so far

  is.get_int(nElems);
  if (nElems < 1) {
    std::stringstream ss (std::stringstream::in | 
			  std::stringstream::out);
    ss << "Line " << linecount_ << ": " 
       << filename1.c_str() << " contains " << nElems << " elements";
    fb.close();
    tExceptionObject(__FILE__,__LINE__,ss.str());
  }

  is.get_int(nVerts);
  if (nVerts != 3) {
    std::stringstream ss (std::stringstream::in | 
			  std::stringstream::out);
    ss << "Line " << linecount_ << ": " 
       << "Mesh element has to have 3 vertices";
    fb.close();
    tExceptionObject(__FILE__,__LINE__,ss.str());
  }

  is.get_int(nAttrs);
  if (nAttrs < 0) {
    std::stringstream ss (std::stringstream::in | 
			  std::stringstream::out);
    ss << "Line " << linecount_ << ": " 
       << "Negative value for number of attributes";
    fb.close();
    tExceptionObject(__FILE__,__LINE__,ss.str());
  }

  // allocate space to store node indices
  //
  // i'm using a dummy variable instead
  int index;

  // allocate space to store attributes
  // we're not using this information, so
  // i'm using a dummy variable here as well
  int attrib;

  // auxiliary variable to store face information
  tCQMTriangle2 tri;

  //
  // reads node indices
  //
  int pNumber = -1;
  int eNumber;
  int i, j;
  for (i=0; i<nElems; i++) {

    // reads element number
    is.get_int(eNumber);

    // verifies if element numbering starts with 0 or 1
    if ((i == 0) && (eNumber == 1)) {
      pNumber = 0;
    }
    
    // verifies element number sequence
    if (eNumber != (pNumber + 1)) {
      std::stringstream ss (std::stringstream::in | 
			    std::stringstream::out);
      ss << "Line " << linecount_ << ": " 
         << "Elements are not numbered correctly";
      fb.close();
      tExceptionObject(__FILE__,__LINE__,ss.str());
    }

    // updates pNumber
    pNumber = eNumber;

    // reads vertex indices
    for (j=0; j<nVerts; j++) {
      if (!is.get_int(index)) {
        std::stringstream ss (std::stringstream::in | 
			      std::stringstream::out);
        ss << "Line " << linecount_ << ": " 
           << "Problem reading " << j+1 << "-th node index of element " << eNumber;
        fb.close();
	tExceptionObject(__FILE__,__LINE__,ss.str());
      }

      // insert vertex in face
      tri.insert(j,index,lverts[index]);
    }

    // insert face into list of faces
    lfaces.push_back(tri);

    // reads attributes of eNumber
    // you will not use this information, so keep this loop as it is
    for (j=0; j<nAttrs; j++) {
      if (!is.get_int(attrib)) {
        std::stringstream ss (std::stringstream::in | 
			      std::stringstream::out);
        ss << "Line " << linecount_ << ": " 
           << "Problem reading " << j+1 << "-th attribute of element " << eNumber;
        fb.close();
	tExceptionObject(__FILE__,__LINE__,ss.str());
      }
    }
  }

  // close file buffer
  fb.close();
}

// -------------------------------------------------------------------
// class  tCQMMeshWriter: implementation
// -------------------------------------
void
tCQMMeshWriter::writeVertexFile(const std::map<int,Point2>& lv) 
  throw (tExceptionObject)
{
  // create an output file
  std::filebuf fb;
  if (fb.open(filename1.c_str(), std::ios::out) == 0) {
    std::stringstream ss (std::stringstream::in | 
			  std::stringstream::out);
    ss << "Unable to create output file";
    fb.close();
    throw tExceptionObject(__FILE__,__LINE__,ss.str());
  }

  std::ostream ostr(&fb);

  // write first line:
  //
  // <# of vertices> <dimension> <# of attributes> <boundary markers (0 or 1)> 
  //
  int nVert = lv.size();
  int dim = 2;
  int nAttr = 0;
  int bMarker = 0;

  ostr << nVert << "\t" << dim << "\t" << nAttr << "\t" 
       << bMarker << std::endl;

  std::map<int,Point2>::const_iterator pos;
  for (pos = lv.begin(); pos != lv.end(); ++pos) {
    int vn = pos->first;
    double x = pos->second.x();
    double y = pos->second.y();

    // write vertex information
    //
    // <# vertex> <x coordinate> <y coordinate>
    //
    ostr << vn << "\t" << x << "\t" << y << std::endl;
  }

  fb.close();
}

void
tCQMMeshWriter::writeEdgeFile(const std::map<tCQMEdgeKey,bool>& le)
  throw (tExceptionObject)
{
  // create an output file
  std::filebuf fb;
  if (fb.open(filename2.c_str(), std::ios::out) == 0) {
    std::stringstream ss (std::stringstream::in | 
			  std::stringstream::out);
    ss << "Unable to create output file";
    fb.close();
    throw tExceptionObject(__FILE__,__LINE__,ss.str());
  }

  std::ostream ostr(&fb);

  // write first line:
  //
  // <# of edges><boundary markers (0 or 1)> 
  //
  int nEdge = le.size();
  int bMarker = 1;

  ostr << nEdge << "\t" << bMarker << std::endl;

  std::map<tCQMEdgeKey,bool>::const_iterator pos;
  int count = 0;
  for (pos = le.begin(); pos != le.end(); ++pos) {
    int e1 = pos->first.getEdgeC1();
    int e2 = pos->first.getEdgeC2();
    ++count;
    if (pos->second) {
      ostr << count << "\t" << e1 << "\t" << e2 << "\t" 
	   << "1" << std::endl;
    }
    else {
      ostr << count << "\t" << e1 << "\t" << e2 << "\t" 
	   << "0" << std::endl;
    }
  }

  fb.close();
}

void
tCQMMeshWriter::writeElemFile(const std::list<tCQMQuadrilateral2*>& lp)
  throw (tExceptionObject)
{
  // create an output file
  std::filebuf fb;
  if (fb.open(filename3.c_str(), std::ios::out) == 0) {
    std::stringstream ss (std::stringstream::in | 
			  std::stringstream::out);
    ss << "Unable to create output file";
    fb.close();
    throw tExceptionObject(__FILE__,__LINE__,ss.str());
  }

  std::ostream ostr(&fb);

  // write first line:
  //
  // <# of elements><# of vertices per element><# of attributes>
  //
  int nElem = lp.size();
  if (nElem == 0) {
    std::stringstream ss (std::stringstream::in | 
			  std::stringstream::out);
    ss << "List of polygons is empty";
    fb.close();
    throw tExceptionObject(__FILE__,__LINE__,ss.str());
  }

  int nVert = 4;
  int nAttr = 0;

  ostr << nElem << "\t" << nVert << "\t" << nAttr << std::endl;

  // write the elements
  std::list<tCQMQuadrilateral2*>::const_iterator li = lp.begin();
  unsigned int count = 0;
  for (; li != lp.end(); ++li) {
    tCQMQuadrilateral2* p = *li;

    if (p->getNumVerts() != nVert) {
      std::stringstream ss (std::stringstream::in | 
			    std::stringstream::out);
      ss << "Polygons have distinct number of vertices";
      fb.close();
      throw tExceptionObject(__FILE__,__LINE__,ss.str());
      fb.close();
    }

    // write element index
    ++count; ostr << count;

    // find the (lowest) leftmost vertex
    Point2 qll = p->getPoint(0);
    int ill = 0;
    int i;
    for (i=1; i<nVert; i++) {
      Point2 qllaux = p->getPoint(i);
      if ((qll.x() > qllaux.x()) || ((qll.x() == qllaux.x()) && (qll.y() > qllaux.y()))) {
	qll = qllaux;
	ill = i;
      }
    }
    
    // write the index of the element vertices
    for (int i=0; i<nVert; i++) {
      ostr << "\t" << p->getVertexId((ill + i) % nVert);
    }

    ostr << std::endl;
  }

  fb.close();
}

void
tCQMMeshWriter::writePolyFile(
			      const std::map<int,Point2>& lv,
			      const std::map<tCQMEdgeKey,bool>& le
			      )
  throw (tExceptionObject)
{
  // create an output file
  std::filebuf fb;
  if (fb.open(filename4.c_str(), std::ios::out) == 0) {
    std::stringstream ss (std::stringstream::in | 
			  std::stringstream::out);
    ss << "Unable to create output file";
    fb.close();
    throw tExceptionObject(__FILE__,__LINE__,ss.str());
  }

  std::ostream ostr(&fb);

  // write first line:
  //
  // <# of vertices> <dimension> <# of attributes> <# of boundary markers (0 or 1)>

  int nVert = lv.size();
  int dim = 2;
  int nAttr = 0;
  int bMarker = 1;

  ostr << nVert << "\t" << dim << "\t" << nAttr << "\t" 
       << bMarker << std::endl;

  std::map<int,Point2>::const_iterator pos1;
  for (pos1 = lv.begin(); pos1 != lv.end(); ++pos1) {
    int vn = pos1->first;
    double x = pos1->second.x();
    double y = pos1->second.y();

    // write vertex information
    //
    // <# vertex> <x coordinate> <y coordinate>
    //
    ostr << vn << "\t" << x << "\t" << y << "\t" << 1 << std::endl;
  }

  // create a list of constrained edges
  std::list<tCQMEdgeKey> lcedges;

  std::map<tCQMEdgeKey,bool>::const_iterator pos2;
  for (pos2 = le.begin(); pos2 != le.end(); ++pos2) {
    if (pos2->second) {
      lcedges.push_back(pos2->first);
    }
  }

  // write next line
  //
  // <# of segments> <# of boundary markers (0 or 1)>
  //
  ostr << lcedges.size() << "\t" << "1" << std::endl;

  // write following lines
  //
  // <segment #> <endpoint> <endpoint> [boundary marker]
  //

  std::list<tCQMEdgeKey>::const_iterator pos3;
  int count = 0;
  for (pos3 = lcedges.begin(); pos3 != lcedges.end(); ++pos3) {
    int e1 = pos3->getEdgeC1();
    int e2 = pos3->getEdgeC2();
    ++count;
    ostr << count << "\t" << e1 << "\t" << e2 << "\t" 
	 << "1" << std::endl;
  }

  // write number of holes
  ostr << "0" << std::endl;

  fb.close();
}
