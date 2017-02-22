//  ===================================================================
//
//  Program:   Convex Quadrilateral Mesh Generator (CQMesh)
//  Module:    err.h
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

#ifndef __err_h
#define __err_h

#include <string>
#include <stdexcept>

#define treatException(e) \
  cout << "Exception: " << e.getDescription() << endl \
       << "File: " << e.getFile() << endl \
       << "Line: " << e.getLine() << endl;

// -------------------------------------------------------------------
// class tExceptionObject
//
// This class extends the class  exception of STL and provides us with
// a customized way of dealing with exceptions.
//
// -------------------------------------------------------------------

class tExceptionObject : public std::exception {
public:

  // Constructors
  tExceptionObject(const char *file="Unknown", unsigned int ln=0,
                   const char *desc="None", const char *loc="Unknown")
  {
    m_loca = loc;
    m_desc = desc;
    m_file = file;
    m_line = ln;
  }

  tExceptionObject(const std::string& file, unsigned int ln,
                   const std::string& desc="None", 
		   const std::string& loc="Unknown")
  {
    m_loca = loc;
    m_desc = desc;
    m_file = file;
    m_line = ln;
  }

  tExceptionObject(const tExceptionObject &xpt) : exception()
  {
    m_loca = xpt.m_loca;
    m_desc = xpt.m_desc;
    m_file = xpt.m_file;
    m_line = xpt.m_line;
  }
  
  // Virtual destructor
  virtual ~tExceptionObject() throw() {}

  // Assignment operator
  tExceptionObject &operator=(const tExceptionObject &xpt)
  {
    m_loca = xpt.m_loca;
    m_desc = xpt.m_desc;
    m_file = xpt.m_file;
    m_line = xpt.m_line;

    return *this;
  }
          
  virtual const char *GetNameOfClass() const 
    {return "tExceptionObject";}

  // Set methods  to set location and  description
  virtual void setLocation(const std::string& s)    
    { m_loca = s; }

  virtual void setLocation(const char * s)          
    { m_loca = s; }

  virtual void setDescription(const std::string& s) 
    { m_desc = s; }

  virtual void setDescription(const char *s)       
    { m_desc = s; }

  // Get methods to obtain location, description, file and line number
  // where the error occured.
  virtual const char *getLocation() const 
    { return m_loca.c_str(); }

  virtual const char *getDescription() const 
    { return m_desc.c_str(); }
  
  virtual const char *getFile() const 
    { return m_file.c_str(); }

  virtual unsigned int getLine() const 
    { return m_line; }
  
  virtual const char* what() const throw()
    { return m_desc.c_str(); }
  
private:
  std::string  m_loca;     // location of the error
  std::string  m_desc;     // description of the error
  std::string  m_file;     // file where the error occured
  unsigned int m_line;     // line of the file where the error occured
 
};

#endif   // __err_h
