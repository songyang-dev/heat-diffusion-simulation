#if !defined(__utils_h__)
#define __utils_h__

#include "trimesh_types.h"
#include <map>

// maps

// From edge to cotangents
typedef std::map<trimesh::index_t, std::pair<double, double>> EdgeToCotan;

#endif // __utils_h__
