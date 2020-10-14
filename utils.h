#if !defined(__utils_h__)
#define __utils_h__

#include "trimesh_types.h"
#include <vector>

// maps

// From edge index to cotangents
typedef std::vector<std::pair<double, double>> EdgeToCotan;

// Covert cosine of an angle to cotangent
double cosToCot(double cosine);

#endif // __utils_h__
