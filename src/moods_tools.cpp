// Copyright (C) 2007-2015  Pasi Rastas, Janne H. Korhonen, Petri Martinmäki
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version, or under the terms of the Biopython
// License.


#include "moods.h"
#include "moods_tools.h"
#include "moods_misc.h"
#include "match_types.h"

#include <limits>

using std::vector;
using std::size_t;

namespace MOODS { namespace tools{

const double PVAL_DP_MULTIPLIER = 2000.0;

// Generates a flat background distribution
vector<double> flat_bg(const unsigned int alphabet_size)
{

    vector<double> bg(alphabet_size, 1.0 / alphabet_size);
    return bg;
}

// Calculates the actual background distribution from a DNA sequence
vector<double> bg_from_sequence_dna(const std::string &seq, const double ps)
{
    unsigned int alphabet_size = 4;
    vector<double> bg(alphabet_size);
    vector<unsigned int> counts(5, 0);

    char c;
    // char value;

    vector<unsigned char> m(256, 4);

    m[(unsigned char)'a'] = 0;
    m[(unsigned char)'A'] = 0;

    m[(unsigned char)'c'] = 1;
    m[(unsigned char)'C'] = 1;

    m[(unsigned char)'g'] = 2;
    m[(unsigned char)'G'] = 2;

    m[(unsigned char)'t'] = 3;
    m[(unsigned char)'T'] = 3;

    for (size_t i = 0; i < seq.size(); ++i){
        c = (unsigned char)seq[i];
        counts[m[c]]++;
    }

    for (unsigned int j = 0; j < alphabet_size; ++j){
        bg[j] = (((double)counts[j] + ps)/ ((double)seq.size() + alphabet_size * ps));
    }

    return bg;
}

score_matrix reverse_complement(const score_matrix &mat)
{
    size_t a = mat.size();
    size_t n = mat[0].size();

    score_matrix ret(a, vector<double>(n));

    for(size_t i = 0; i < a; i++) {
        for (size_t j = 0; j < n; j++) {
            ret[i][j] = mat[a - i - 1][n - j - 1];
        }
    }
    return ret;
}

// Transforms a weight matrix into a PSSM
score_matrix log_odds(const score_matrix &mat, const vector<double> &bg, const double ps)
{
    size_t a = mat.size();
    size_t n = mat[0].size();

    score_matrix ret(a, vector<double>(n));

    for (size_t i = 0; i < n; ++i)
    {
        double column_sum = 0;
        for (size_t j = 0; j < a; ++j){
            column_sum += mat[j][i] + ps*bg[j];
        }
        for (size_t j = 0; j < a; ++j) {
            ret[j][i] = log((mat[j][i] + ps*bg[j])/column_sum) - log(bg[j]);
        }
    }
    return ret;
}

// Transforms a weight matrix into a PSSM with non-natural logarithm
score_matrix log_odds(const score_matrix &mat, const vector<double> &bg, const double ps, const double log_base)
{

    size_t a = mat.size();
    size_t n = mat[0].size();

    score_matrix ret = log_odds(mat, bg, ps);

    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < a; ++j) {
            ret[j][i] = ret[j][i] / log(log_base);
        }
    }

    return ret;

}

// // Calculates a threshold for a scoring matrix from a given p value
double threshold_from_p(const score_matrix &pssm, const vector<double> &bg, const double &p)
{

    // Approximate the scoring matrix with integer matrix for DP
    long a = pssm.size();
    long n = pssm[0].size();


    vector<vector<long> > mat(a, vector<long>(n));

    long maxT = 0;
    long minV = std::numeric_limits<long>::max();

    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < a; ++j)
        {
            if (pssm[j][i] > 0.0){
                mat[j][i] = (long) ( PVAL_DP_MULTIPLIER * pssm[j][i] + 0.5 );
            }
            else {
                mat[j][i] = (long) ( PVAL_DP_MULTIPLIER * pssm[j][i] - 0.5 );
            }
        }
    }

    for (size_t i = 0; i < n; ++i)
    {
        long max = mat[0][i];
        long min = max;
        for (size_t j = 1; j < a; ++j)
        {
            long v = mat[j][i];
            if (max < v)
                max = v;
            else if (min > v)
                min = v;
        }
        maxT += max;
        if (minV > min)
            minV = min;
    }

    long R = maxT - n * minV;

    vector<double> table0(R + 1, 0.0);
    vector<double> table1(R + 1, 0.0);

    for (size_t j = 0; j < a; ++j)
        table0[mat[j][0] - minV] += bg[j];

    for (size_t i = 1; i < n; ++i)
    {
        for (size_t j = 0; j < a; ++j)
        {
            long s = mat[j][i] - minV;
            for (long r = s; r <= R; ++r)
                table1[r] += bg[j] * table0[r - s];
        }
        for (long r = 0; r <= R; ++r)
        {
            table0[r] = table1[r];
            table1[r] = 0.0;
        }
    }



    double sum = table0[R];

    if (sum > p){
        return max_score(pssm) - min_delta(pssm)/2;
    }

    for (long r = R-1; r >= 0; --r)
    {
        sum += table0[r];
        if (sum > p)
        {
            return (double) ((r + n * minV + 1) / PVAL_DP_MULTIPLIER);
        }
    }

    return min_score(pssm)-1.0;
}


// Calculates maximum possible score for a PSSM
double max_score(const score_matrix &mat)
{
    size_t a = mat.size();
    size_t n = mat[0].size();

    double ret = 0;
    for (size_t i = 0; i < n; ++i)
    {
        double max = -std::numeric_limits<double>::infinity();
        for (size_t j = 0; j < a; ++j)
        {
            max = std::max(mat[j][i], max);
        }
        ret += max;
    }
    return ret;
}

double min_score(const score_matrix &mat)
{
    size_t a = mat.size();
    size_t n = mat[0].size();

    double ret = 0;
    for (size_t i = 0; i < n; ++i)
    {
        double min = std::numeric_limits<double>::infinity();
        for (size_t j = 0; j < a; ++j)
        {
            min = std::min(mat[j][i], min);
        }
        ret += min;
    }
    return ret;
}

double min_delta(const score_matrix &mat)
{
    size_t a = mat.size();
    size_t n = mat[0].size();

    double min_delta = std::numeric_limits<double>::infinity();
    for (size_t i = 0; i < n; ++i)
    {
        double max = -std::numeric_limits<double>::infinity();
        double second_max = -std::numeric_limits<double>::infinity();
        for (size_t j = 0; j < a; ++j)
        {
            if (mat[j][i] > max){
                second_max = max;
                max = mat[j][i];
            }
            else if (mat[j][i] < max && mat[j][i] > second_max) {
                second_max = mat[j][i];
            }
        }
        min_delta = std::min(min_delta, max - second_max);
    }
    return min_delta;
}


} // namespace tools
} // namespace MOODS
