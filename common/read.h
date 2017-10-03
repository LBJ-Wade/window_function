#ifndef READ_H
#define READ_H 1

#include <vector>
#include "particle.h"

//void read_fof_ascii_file(const char filename[], std::vector<Particle>& v);
//void read_mock_text(const char filename[], std::vector<Particle>& v);
//void read_pos_ascii_file(const char filename[], std::vector<Particle>& v);
//void read_random_catalogue(std::vector<Particle>& v, const size_t np, const float boxsize);
//void read_mock_catalogue(const char filename[], std::vector<Particle>& v, const bool centrals_only);
void read_ascii_fkp(const char filename[], const double Pest,
		    const double subsample_rate, const unsigned long seed,
		    std::vector<Particle>& v);

#endif
