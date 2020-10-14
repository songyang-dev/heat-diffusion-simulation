#include "utils.h"

#include <cmath>

// Turn cos to cotan
double cosToCot(double cosine)
{
    // double sine = std::sqrt(1 - std::pow(cosine, 2));
    // return cosine/sine;

    return 1/tan(acos(cosine));
}