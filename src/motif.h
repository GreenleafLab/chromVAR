#ifndef MOTIF_H
#define MOTIF_H


#include "moods.h"

#include <utility>

using std::vector;
using std::string;
using std::size_t;

namespace MOODS { namespace scan{

class Motif {
public:
    virtual std::pair<bool, double> window_match(bits_t seq, bits_t shift) = 0;
    virtual std::pair<bool, double> check_hit(const std::string& s,
                                              const std::vector<unsigned char>& alphabet_map,
                                              const std::size_t window_match_pos,
                                              double score) = 0;
    virtual unsigned int size()=0;
    virtual unsigned int alphabet_size() = 0;
    virtual unsigned int window_pos() = 0;
    virtual double threshold() = 0;
};


// standard 0-order PWM
class Motif0 : public Motif {
private:
    score_matrix mat;
    std::vector<unsigned int> lookahead_order;
    std::vector<double> lookahead_scores;

    unsigned int l; // window size
    unsigned int m; // length
    unsigned int a; // alphabet size

    unsigned int wp; // window position
    double T;
public:
    Motif0 (const score_matrix& matrix, const vector<double>& bg, unsigned int window_size, double threshold);

    std::pair<bool, double> window_match(bits_t seq, bits_t shift);
    std::pair<bool, double> check_hit(const std::string& s, const std::vector<unsigned char>& alphabet_map, const std::size_t window_match_pos, double score);

    unsigned int size() { return m; }
    unsigned int alphabet_size() { return a; }
    unsigned int window_pos() { return wp; }
    double threshold() { return T; }

};


}}

#endif
