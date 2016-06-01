#include "moods.h"
#include "moods_misc.h"


using std::vector;
using std::size_t;

namespace MOODS { namespace misc{

    // basically base-2 logarithm of a, rounded up
    size_t shift(size_t a)
    {
        size_t s = 0;
        size_t b = 1;
        
        while (b < a){
            s += 1;
            b = b << 1;
        }
        return s;
    }
        
    bits_t mask(size_t a){
        bits_t b = 1;
        
        while (b < a){
            b = b << 1;
        }
        return b-1;
    }

    size_t q_gram_size(size_t rows, size_t a){
        size_t q = 0;
        size_t s = 1;
        

        while (s < rows){
            q += 1;
            s *= a;
        }
        return q;   
    }

    bits_t rc_tuple(bits_t CODE, size_t a, size_t q){
        size_t SHIFT = shift(a);
        bits_t A_MASK = (1 << SHIFT) - 1;

        bits_t RC = 0;
        for (size_t i = 0; i < q; ++i){
            bits_t c = (CODE >> ((q - i - 1) * SHIFT)) & A_MASK;
            RC = RC | ((a - c - 1) << (i * SHIFT));
        }

        return RC;
    }

    // checks a sequence for non-scan regions and returns the corresponding bounds
    std::vector<size_t> preprocess_seq(const std::string& s, size_t a, const std::vector<unsigned char>& alphabet_map){
        
        vector<size_t> bounds;
        
        bool scannable = false;
        unsigned char c;
        
        for (size_t i = 0; i < s.size(); ++i){
            c = alphabet_map[(unsigned char)s[i]];
            
            if (c < a){
                if (!scannable){
                    scannable = true;
                    bounds.push_back(i);
                }
            }
            else {
                if (scannable){
                    scannable = false;
                    bounds.push_back(i);
                }
            }
        }
        if (scannable){
            bounds.push_back(s.size());
        }
        
        return bounds;
        
    }
}}