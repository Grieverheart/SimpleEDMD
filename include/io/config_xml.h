#ifndef __IO_XML_H
#define __IO_XML_H

#include "configuration.h"

//WARNING: Temporary code, until we figure out how we want to layout
//the data in xml format.

bool xml_load_config(const char* filename, Configuration&);
bool xml_save_config(const char* filename, const Configuration&);

#endif
