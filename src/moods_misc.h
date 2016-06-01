#ifndef MOODS_MISC_H
#define MOODS_MISC_H

#include "moods.h"

namespace MOODS { namespace misc{
    
    size_t shift(size_t a);
    bits_t mask(size_t a);
    size_t q_gram_size(size_t rows, size_t a);
    bits_t rc_tuple(bits_t CODE, size_t a, size_t q);

    std::vector<size_t> preprocess_seq(const std::string& s, size_t a, const std::vector<unsigned char>& alphabet_map); 

}}



#endif