/*  $Id: compat.cpp 5413 2024-01-06 19:27:10Z jcherry $
 * ===========================================================================
 *
 *                            PUBLIC DOMAIN NOTICE
 *               National Center for Biotechnology Information
 *
 *  This software/database is a "United States Government Work" under the
 *  terms of the United States Copyright Act.  It was written as part of
 *  the author's official duties as a United States Government employee and
 *  thus cannot be copyrighted.  This software/database is freely available
 *  to the public for use. The National Library of Medicine and the U.S.
 *  Government have not placed any restriction on its use or reproduction.
 *
 *  Although all reasonable efforts have been taken to ensure the accuracy
 *  and reliability of the software and data, the NLM and the U.S.
 *  Government do not and cannot warrant the performance or results that
 *  may be obtained by using this software or data. The NLM and the U.S.
 *  Government disclaim all warranties, express or implied, including
 *  warranties of performance, merchantability or fitness for any particular
 *  purpose.
 *
 *  Please cite the author in any work or product based on this material.
 *
 * ===========================================================================
 *
 * Authors:  Joshua L. Cherry
 *
 * File Description:  Maximum compatibility phylogentic reconstruction
 *
 */


#include <vector>
#include <set>
#include <map>
#include <iostream>
#include <algorithm>
#include <utility>
#include <stdexcept>

#include "compat.hpp"

using namespace std;

#define COUNT 1

// Comparison function for sorting in "reverse" order,
// i.e., from largest to smallest.
inline bool my_cmp(unsigned int lhs, unsigned int rhs) {
    return lhs > rhs;
}


// Does intersection and returns total weight of removed nodes
unsigned int do_intersection(const vector<unsigned int>& s,
                             const vector<char>& conn,
                             const vector<unsigned int>& weights,
                             vector<unsigned int>& result) {
    unsigned int removed_weight = 0;
    unsigned int sz = s.size();
    const unsigned int *ps = &s[0];
    result.resize(sz);
    unsigned int *presult = &result[0];
    const char *pconn = &conn[0];
    for (unsigned int i = 0; i < sz; ++i) {
        if (pconn[ps[i]]) {
            *presult = ps[i];
            presult++;
        } else{
            removed_weight += weights[ps[i]];
        }
    }
    result.resize(presult - &result[0]);
    return removed_weight;
}


template<class T>
unsigned int Weight(const T& nodes,
                    const vector<unsigned int>& weights) {
    unsigned int sum = 0;
    for (unsigned i : nodes) {
        sum += weights[i];
    }
    return sum;
}


MaxClique::MaxClique() {
    m_Report_maxes = true;
    m_Report_progress = true;
    m_Upper_bound = e_None;
    ub_calls = 0;
    intersections = 0;
    pointless_intersections = 0;
}


void MaxClique::SetColorings(const vector<vector<unsigned> >& colorings) {
    m_Colorings = colorings;
    m_NColors.clear();
    for (auto colors : colorings) {
        m_NColors.push_back(*max_element(colors.begin(), colors.end() ) + 1);
    }
};


vector<unsigned int>*
MaxClique::Recurse(vector<unsigned int>& nodes,
                   const vector<unsigned int>& weights,
                   const vector<vector<char> >& conn,
                   unsigned int weight,
                   unsigned int nodes_weight)
{
    if (nodes.size() == 0) {
        if (weight >= mx) {
            if (weight > mx) {
                mx = weight;
            }
            if (m_Report_maxes) {
                cerr << "max of " << mx << endl;
            }
            return new vector<unsigned int>;
        }
        return 0;
    }
    vector<unsigned int> *my_best = 0;
    unsigned int i;
    vector<unsigned int> inter;
    unsigned int removed_weight;

    bool done = false;
    while (nodes.size() != 0 && !done) {
        if (weight + nodes_weight < mx) {
            break;
        }
        i = nodes.back();
        nodes_weight -= weights[i];
        nodes.pop_back();
        if (weight == 0) {
            if (m_Report_progress) {
                cerr << i << endl;
            }
        }
        
        removed_weight = do_intersection(nodes, conn[i], weights, inter);

#ifdef COUNT
        intersections++;
#endif
        if (inter.size() == nodes.size()) {
            done = true;
#ifdef COUNT
            pointless_intersections++;
#endif
        }
        vector<unsigned int> *res = Recurse(inter,
                                            weights,
                                            conn,
                                            weight + weights[i],
                                            nodes_weight - removed_weight);
        if (res) {
            if (my_best) {
                delete my_best;
            }
            my_best = res;
            res->push_back(i);
        }
    }
    return my_best;
}


bool
MaxClique::RecurseHCoW(vector<unsigned int>& nodes,
                       const vector<unsigned int>& weights,
                       const vector<vector<char> >& conn,
                       unsigned int weight,
                       unsigned int nodes_weight)
{
    if (nodes.size() == 0) {
        if (weight >= mx) {
            if (weight > mx) {
                mx = weight;
            }
            if (m_Report_maxes) {
                cerr << "max of " << mx << endl;
            }
            return true;
        }
        return false;
    }
    unsigned int i;
    vector<unsigned int> inter;
    unsigned int removed_weight;

    bool done = false;
    while (nodes.size() != 0 && !done) {
        if (weight + nodes_weight < mx) {
            break;
        }
        i = nodes.back();
        nodes_weight -= weights[i];
        nodes.pop_back();
        if (weight == 0) {
            if (m_Report_progress) {
                cerr << i << endl;
            }
        }
        
        removed_weight = do_intersection(nodes, conn[i], weights, inter);

#ifdef COUNT
        intersections++;
#endif
        if (inter.size() == nodes.size()) {
            done = true;
#ifdef COUNT
            pointless_intersections++;
#endif
        }
        bool res = RecurseHCoW(inter,
                               weights,
                               conn,
                               weight + weights[i],
                               nodes_weight - removed_weight);
        if (res) {
            return true;
        }
    }
    return false;
}


vector<unsigned int>
MaxClique::Run(const vector<unsigned int>& nodes,
               const vector<unsigned int>& weights,
               const vector<vector<char> >& conn,
               unsigned int lower_bound)
{
    mx = lower_bound;

    vector<unsigned int> sorted_nodes = nodes;
    sort(sorted_nodes.begin(), sorted_nodes.end(), my_cmp);

    unsigned int nodes_weight = Weight(nodes, weights);

    vector<unsigned int> *res = Recurse(sorted_nodes, weights, conn,
                                        0, nodes_weight);
    if (res == 0) {
        return vector<unsigned int>();
    }
    return *res;
}


vector<unsigned int>
MaxClique::Run2(const vector<unsigned int>& nodes,
                const vector<unsigned int>& weights,
                const vector<vector<char> >& conn,
                const vector<unsigned int>& ranges,
                unsigned int lower_bound,
                const vector<unsigned int>& regroup_at)
{
    mx = lower_bound;

    vector<unsigned int> sorted_nodes = nodes;
    sort(sorted_nodes.begin(), sorted_nodes.end(), my_cmp);

    m_RegroupAt = regroup_at;
    sort(m_RegroupAt.begin(), m_RegroupAt.end(), my_cmp);
    m_Groups.clear();

    unsigned int nodes_weight = Weight(nodes, weights);

    vector<unsigned int> weights_copy(weights);  // So we can pass non-const &

    vector<unsigned int> *res = Recurse2(sorted_nodes, weights_copy,
                                         conn, ranges, 0, nodes_weight);
    if (res == 0) {
        return vector<unsigned int>();
    }
    return *res;
}


bool
MaxClique::HasCliqueOfWeight(const vector<unsigned int>& nodes,
                             const vector<unsigned int>& weights,
                             const vector<vector<char> >& conn,
                             unsigned int lower_bound)
{
    mx = lower_bound;

    vector<unsigned int> sorted_nodes = nodes;
    sort(sorted_nodes.begin(), sorted_nodes.end(), my_cmp);

    unsigned int nodes_weight = Weight(nodes, weights);

    return RecurseHCoW(sorted_nodes, weights, conn, 0, nodes_weight);
}


bool Equal(const vector<char>& c1, const vector<char>& c2, unsigned int start) {
    for (unsigned int i = start; i < c1.size(); ++i) {
        if (c1[i] != c2[i]) {
            return false;
        }
    }
    return true;
}


void
Regroup(vector<unsigned int>& weights,
        const vector<unsigned int>& ranges,
        unsigned int start,
        const vector<vector<char> >& conn,
        vector<unsigned int>& nodes,
        vector<vector<unsigned int> >& groups) {

    // Figure out where ranges start
    vector<unsigned int> range_starts;
    range_starts.push_back(start);
    for (unsigned int i = start + 1; i < ranges.size(); ++i) {
        if (ranges[i] != range_starts.back()) {
            range_starts.push_back(i);
        }
    }
    
    // Crude hashing to avoid lots of expensive comparisons
    vector<unsigned int> nconns(conn.size());
    vector<unsigned int> sum_confs(conn.size());
    for (unsigned int i = start; i < conn.size(); ++i) {
        if (weights[i] == 0) {
            continue;
        }
        for (unsigned int j = start; j < conn.size(); ++j) {
            nconns[i] += conn[i][j];
            if (!conn[i][j]) {
                sum_confs[i] += j;
            }
        }
    }

    // Iterate over ranges
    unsigned int begin, end;
    for (unsigned int i = 0; i < range_starts.size(); ++i) {
        begin = range_starts[i];
        end = ranges[begin];
        for (unsigned int j = begin; j < end; ++j) {
            if (weights[j] == 0) {
                continue;
            }
            for (unsigned int k = j + 1; k <= end; ++k) {
                if (weights[k] == 0) {
                    continue;
                }
                if (sum_confs[j] == sum_confs[k]
                    && nconns[j] == nconns[k]
                    && Equal(conn[j], conn[k], start)) {
                    weights[j] += weights[k];
                    weights[k] = 0;
                    for (unsigned int n = 0; n < groups[k].size(); ++n) {
                        groups[j].push_back(groups[k][n]);
                    }
                }
            }
        }
    }

    // Remove nodes as appropriate
    vector<unsigned int> tmp;
    tmp.reserve(nodes.size());
    for (unsigned int i = 0; i < nodes.size(); ++i) {
        if (weights[nodes[i]]) {
            tmp.push_back(nodes[i]);
        }
    }
    nodes.swap(tmp);
}


vector<unsigned int>*
MaxClique::Recurse2(vector<unsigned int>& nodes,
                    vector<unsigned int>& weights,
                    const vector<vector<char> >& conn,
                    const vector<unsigned int>& ranges,
                    unsigned int weight,
                    unsigned int nodes_weight)
{
    if (nodes.size() == 0) {
        if (weight >= mx) {
            if (weight > mx) {
                mx = weight;
            }
            if (m_Report_maxes) {
                cerr << "max of " << mx << endl;
            }
            return new vector<unsigned int>;
        }
        return 0;
    }
    vector<unsigned int> *my_best = 0;
    unsigned int i;
    vector<unsigned int> inter;
    unsigned int removed_weight;
    unsigned int included_weight;
    bool done = false;
    while (nodes.size() != 0 && !done) {
        if (weight + nodes_weight < mx) {
            break;
        }

        if (weight == 0) {
            if (m_RegroupAt.size() && nodes.back() >= m_RegroupAt.back()) {
                if (m_Groups.size() == 0) {
                    m_Groups.resize(conn.size());
                    for (unsigned int i = 0; i < m_Groups.size(); ++i) {
                        m_Groups[i].push_back(i);
                    }
                }
                Regroup(weights, ranges, nodes.back(), conn, nodes, m_Groups);
                while (m_RegroupAt.size()
                       && nodes.back() >= m_RegroupAt.back()) {
                    m_RegroupAt.pop_back();
                }
            }
        }



        vector<unsigned int> included;
        included_weight = 0;
        i = nodes.back();
        while (nodes.size() > 0 && nodes.back() <= ranges[i]) {        
            nodes_weight -= weights[nodes.back()];
            included_weight += weights[nodes.back()];
            included.push_back(nodes.back());
            nodes.pop_back();
        }
        if (weight == 0) {
            if (m_Report_progress) {
                cerr << i << endl;
            }
        }
        
        removed_weight = do_intersection(nodes, conn[i], weights, inter);

#ifdef COUNT
        intersections++;
#endif
        if (inter.size() == nodes.size()) {
            done = true;
#ifdef COUNT
            pointless_intersections++;
#endif
        }
        vector<unsigned int> *res = Recurse2(inter,
                                             weights,
                                             conn,
                                             ranges,
                                             weight + included_weight,
                                             nodes_weight - removed_weight);
        if (res) {
            if (my_best) {
                delete my_best;
            }
            my_best = res;
            for (unsigned int n = 0; n < included.size(); ++n) {
                res->push_back(included[n]);
            }
            if (weight == 0 && m_Groups.size()) {
                // need to expand groups
                vector<unsigned int> tmp;
                for (unsigned int i = 0; i < res->size(); ++i) {
                    for (unsigned int j = 0;
                         j < m_Groups[(*res)[i]].size();
                         ++j) {
                        tmp.push_back(m_Groups[(*res)[i]][j]);
                    }
                }
                res->swap(tmp);
            }
        }

    }
    return my_best;
}


vector<vector<unsigned int> >
MaxClique::FindAll(const vector<unsigned int>& nodes,
                   const vector<unsigned int>& weights,
                   const vector<vector<char> >& conn,
                   const vector<unsigned int>& ranges,
                   unsigned int lower_bound,
                   const vector<unsigned int>& regroup_at)
{
    mx = lower_bound;

    vector<unsigned int> sorted_nodes = nodes;
    sort(sorted_nodes.begin(), sorted_nodes.end(), my_cmp);

    m_RegroupAt = regroup_at;
    sort(m_RegroupAt.begin(), m_RegroupAt.end(), my_cmp);
    m_Groups.clear();

    unsigned int nodes_weight = Weight(nodes, weights);

    vector<unsigned int> weights_copy(weights);  // So we can pass non-const &

    vector<vector<unsigned int> > *res;
    if (m_Upper_bound == e_None) {
        res = RecurseFA(sorted_nodes, weights_copy,
                        conn, ranges, 0, nodes_weight);
    } else {
        res = RecurseFAColors(sorted_nodes, weights_copy,
                              conn, ranges, 0, nodes_weight);
    }
    if (res == 0) {
        return vector<vector<unsigned int> >(1);
    }
    return *res;
}


static void s_ExpandGroups(vector<unsigned>& nodes,
                           const vector<vector<unsigned int> >& groups) {
    vector<unsigned int> tmp;
    for (unsigned int i = 0; i < nodes.size(); ++i) {
        for (unsigned int j = 0;
             j < groups[nodes[i]].size();
             ++j) {
            tmp.push_back(groups[nodes[i]][j]);
        }
    }
    nodes.swap(tmp);
}


static void s_ExpandGroups(vector<vector<unsigned> >& node_sets,
                           const vector<vector<unsigned int> >& groups) {
    for (unsigned int i = 0; i < node_sets.size(); ++i) {
        s_ExpandGroups(node_sets[i], groups);
    }
}


// Recursion for finding all max. weight cliques
vector<vector<unsigned int> >*
MaxClique::RecurseFA(vector<unsigned int>& nodes,
                     vector<unsigned int>& weights,
                     const vector<vector<char> >& conn,
                     const vector<unsigned int>& ranges,
                     unsigned int weight,
                     unsigned int nodes_weight)
{
    if (nodes.size() == 0) {
        if (weight >= mx) {
            if (weight > mx) {
                mx = weight;
            }
            if (m_Report_maxes) {
                cerr << "max of " << mx << endl;
            }
            return new vector<vector<unsigned int> >(1);
        }
        return 0;
    }
    vector<vector<unsigned int> > *my_best = 0;
    unsigned int last_mx = mx;
    unsigned int i;
    vector<unsigned int> inter;
    unsigned int removed_weight;
    unsigned int included_weight;
    bool done = false;
    while (nodes.size() != 0 && !done) {
        if (weight + nodes_weight < mx) {
            break;
        }

        if (weight == 0) {
            if (m_RegroupAt.size() && nodes.back() >= m_RegroupAt.back()) {
                if (m_Groups.size() == 0) {
                    m_Groups.resize(conn.size());
                    for (unsigned int i = 0; i < m_Groups.size(); ++i) {
                        m_Groups[i].push_back(i);
                    }
                }
                Regroup(weights, ranges, nodes.back(), conn, nodes, m_Groups);
                while (m_RegroupAt.size()
                       && nodes.back() >= m_RegroupAt.back()) {
                    m_RegroupAt.pop_back();
                }
            }
        }



        vector<unsigned int> included;
        included_weight = 0;
        i = nodes.back();
        while (nodes.size() > 0 && nodes.back() <= ranges[i]) {        
            nodes_weight -= weights[nodes.back()];
            included_weight += weights[nodes.back()];
            included.push_back(nodes.back());
            nodes.pop_back();
        }
        if (weight == 0) {
            if (m_Report_progress) {
                cerr << i << endl;
            }
        }
        
        removed_weight = do_intersection(nodes, conn[i], weights, inter);

#ifdef COUNT
        intersections++;
#endif
        if (inter.size() == nodes.size()) {
            done = true;
#ifdef COUNT
            pointless_intersections++;
#endif
        }
        vector<vector<unsigned int> > *res;
        res = RecurseFA(inter,
                        weights,
                        conn,
                        ranges,
                        weight + included_weight,
                        nodes_weight - removed_weight);
        if (res) {
            for (unsigned int m = 0; m < res->size(); ++m) {
                for (unsigned int n = 0; n < included.size(); ++n) {
                    (*res)[m].push_back(included[n]);
                }
            }

            if (weight == 0 && m_Groups.size()) {
                s_ExpandGroups(*res, m_Groups);
            }

            if (my_best && (mx == last_mx)) {
                // Combine results
                my_best->resize(my_best->size() + res->size());
                for (unsigned n = 0; n < res->size(); ++n) {
                    (*my_best)[n + my_best->size() - res->size()].swap((*res)[n]);
                }
                delete res;
            } else {
                if (my_best) {
                    delete my_best;
                }
                my_best = res;
                last_mx = mx;
            }

        }

    }
    return my_best;
}


vector<vector<unsigned int> >*
MaxClique::RecurseFAColors(vector<unsigned int>& nodes,
                           vector<unsigned int>& weights,
                           const vector<vector<char> >& conn,
                           const vector<unsigned int>& ranges,
                           unsigned int weight,
                           unsigned int nodes_weight)
{
    if (nodes.size() == 0) {
        if (weight >= mx) {
            if (weight > mx) {
                mx = weight;
            }
            if (m_Report_maxes) {
                cerr << "max of " << mx << endl;
            }
            return new vector<vector<unsigned int> >(1);
        }
        return 0;
    }

    if (false) {
        unsigned color_ub = ColorUpperBound(nodes, m_Colorings[nodes.back()], weights);
        if (color_ub + weight < mx) {
            return 0;
        }
    }
    vector<vector<unsigned int> > *my_best = 0;
    unsigned int last_mx = mx;
    unsigned int i;
    vector<unsigned int> inter;
    unsigned int nodes_weight_last_cub = nodes_weight + 1;  // Force first ub
    unsigned int cub_excess = 0; 
    unsigned int removed_weight;
    unsigned int included_weight;
    bool done = false;
    while (nodes.size() != 0 && !done) {
        if (weight + nodes_weight < mx) {
            break;
        }

        if (weight == 0) {
            if (m_RegroupAt.size() && nodes.back() >= m_RegroupAt.back()) {
                if (m_Groups.size() == 0) {
                    m_Groups.resize(conn.size());
                    for (unsigned int i = 0; i < m_Groups.size(); ++i) {
                        m_Groups[i].push_back(i);
                    }
                }
                Regroup(weights, ranges, nodes.back(), conn, nodes, m_Groups);
                while (m_RegroupAt.size()
                       && nodes.back() >= m_RegroupAt.back()) {
                    m_RegroupAt.pop_back();
                }
            }
        }

        if ((nodes_weight_last_cub - nodes_weight) > cub_excess) {
            unsigned color_ub;
            if (m_Upper_bound == e_Coloring) {
                color_ub = ColorUpperBound(nodes,
                                           m_Colorings[nodes.back()],
                                           weights,
                                           m_NColors[nodes.back()]);
            } else {
                color_ub = DynProgUpperBound(nodes,
                                             m_Dyn_progs[nodes.back()],
                                             conn,
                                             weights);
            }
# if COUNT
            ub_calls++;
# endif
            if (color_ub + weight < mx) {
                break;
            }
            nodes_weight_last_cub = nodes_weight;
            cub_excess = color_ub + weight - mx;
        }


        vector<unsigned int> included;
        included_weight = 0;
        i = nodes.back();
        while (nodes.size() > 0 && nodes.back() <= ranges[i]) {        
            nodes_weight -= weights[nodes.back()];
            included_weight += weights[nodes.back()];
            included.push_back(nodes.back());
            nodes.pop_back();
        }
        if (weight == 0) {
            if (m_Report_progress) {
                cerr << i << endl;
            }
        }
        
        removed_weight = do_intersection(nodes, conn[i], weights, inter);

#ifdef COUNT
        intersections++;
        if (intersections_hist.size() < nodes.size() + 1) {
            intersections_hist.resize(nodes.size() + 1);
        }
        intersections_hist[nodes.size()]++;
#endif
        if (inter.size() == nodes.size()) {
            done = true;
#ifdef COUNT
            pointless_intersections++;
#endif
        }
        if (removed_weight < weights[i]) {
            // any maximum clique of this subset must contain node i
            done = true;
        }
        vector<vector<unsigned int> > *res;
        res = RecurseFAColors(inter,
                              weights,
                              conn,
                              ranges,
                              weight + included_weight,
                              nodes_weight - removed_weight);
        if (res) {
            for (unsigned int m = 0; m < res->size(); ++m) {
                for (unsigned int n = 0; n < included.size(); ++n) {
                    (*res)[m].push_back(included[n]);
                }
            }

            if (weight == 0 && m_Groups.size()) {
                s_ExpandGroups(*res, m_Groups);
            }

            if (my_best && (mx == last_mx)) {
                // Combine results
                my_best->resize(my_best->size() + res->size());
                for (unsigned n = 0; n < res->size(); ++n) {
                    (*my_best)[n + my_best->size() - res->size()].swap((*res)[n]);
                }
                delete res;
            } else {
                if (my_best) {
                    delete my_best;
                }
                my_best = res;
                last_mx = mx;
            }

        }

    }
    return my_best;
}


bool IsSubset(const vector<unsigned int>& s1,
              const vector<unsigned int>& s2){
    if (s1.size() > s2.size()) {
        return false;
    }
    const unsigned int *p1 = &s1[0];
    const unsigned int *p2 = &s2[0];
    unsigned int i1 = 0, i2 = 0;
    while (i1 < s1.size() && i2 < s2.size()) {
        if (*p1 < *p2) {
            return false;
        }
        if (*p1 == *p2) {
            ++p1, ++i1, ++p2, ++i2;
        } else {  // *p1 > *p2
            ++p2, ++i2;
        }
    }
    return i1 == s1.size();
}


bool Disjoint(const vector<unsigned int>& s1,
              const vector<unsigned int>& s2) {
    if (s1.front() > s2.back() || s2.front() > s1.back()) {
        return true;  // Large performance gain with tree-like renumbering
    }
    const unsigned int *p1 = &s1[0];
    const unsigned int *p2 = &s2[0];
    unsigned int i1 = 0, i2 = 0;
    while (i1 < s1.size() && i2 < s2.size()) {
        if (*p1 == *p2) {
            return false;
        }
        if (*p1 < *p2) {
            ++p1, ++i1;
        } else {  // *p1 > *p2
            ++p2, ++i2;
        }
    }
    return true;
}


bool JointlyExhaustive(const vector<unsigned int>& s1,
                       const vector<unsigned int>& s2,
                       unsigned int nleaves) {
    if (!(s1.front() == 0 || s2.front() == 0)) {
        return false;
    }
    if (s1.size() + s2.size() < nleaves) {
        return false;
    }
    const unsigned int *p1 = &s1[0];
    const unsigned int *p2 = &s2[0];
    unsigned int i1 = 0, i2 = 0;
    unsigned int sought = 0;
    bool found1 = false, found2 = false;

    while (i1 < s1.size() || i2 < s2.size()) {

        found1 = i1 < s1.size() && *p1 == sought;
        found2 = i2 < s2.size() && *p2 == sought;

        if (!(found1 || found2)) {
            return false;
        }

        if (found1) {
            ++p1, ++i1;
        }

        if (found2) {
            ++p2, ++i2;
        }

        ++sought;
        if (sought == nleaves) {
            return true;
        }
    }
    return false;
}


bool Compatible(const pair<vector<unsigned int>, vector<unsigned int> >& s1,
                const pair<vector<unsigned int>, vector<unsigned int> >& s2,
                unsigned int nleaves) {
    return s1.first.size() == 1
        || s2.first.size() == 1
        || IsSubset(s1.first, s2.second)
        || IsSubset(s2.first, s1.second)
        || Disjoint(s1.first, s2.first)
        || JointlyExhaustive(s1.second, s2.second, nleaves);
}


bool CannotConflict(const pair<vector<unsigned int>, vector<unsigned int> >& s1,
                    const pair<vector<unsigned int>, vector<unsigned int> >& s2,
                    unsigned int nleaves) {
    return s1.second.size() == 1
        || s2.second.size() == 1
        || Disjoint(s1.second, s2.second)
        || IsSubset(s1.second, s2.first)
        || IsSubset(s2.second, s1.first);
        //|| JointlyExhaustive(s1.first, s2.first, nleaves);
}


// Identify all unamb that cannot conflict with any ambiguous, however resolved
vector<unsigned int>
CannotConflict(vector<pair<vector<unsigned int>, vector<unsigned int> > > sets,
               unsigned int nleaves) {
    vector<unsigned int> rv;
    bool cannot_conflict;
    for (unsigned int i = 0; i < sets.size(); ++i) {
        pair<vector<unsigned int>, vector<unsigned int>> &s = sets[i];
        if (s.first.size() != s.second.size()) {
            // ambiguous
            continue;
        }
        cannot_conflict = true;
        for (unsigned int j = 0; j < sets.size(); ++j) {
            if (sets[j].first.size() == sets[j].second.size()) {
                // unambiguous
                continue;
            }
            if (!CannotConflict(s, sets[j], nleaves)) {
                cannot_conflict = false;
                break;
            }
        }
        if (cannot_conflict) {
            rv.push_back(i);
        }
    }
    return rv;
}


typedef vector<pair<vector<unsigned int>, vector<unsigned int> > > TSetPair;

// Comparison object for sorting
class CCmp {
public:
    CCmp(const TSetPair& sets) : m_Sets(sets) {}
    bool operator () (unsigned int lhs, unsigned int rhs) {
        if (m_Sets[lhs].first.size() != m_Sets[rhs].first.size()) {
            return m_Sets[lhs].first.size() < m_Sets[rhs].first.size();
        }
        return m_Sets[lhs].first < m_Sets[rhs].first;
    }
private:
    const TSetPair m_Sets;
};

void
Conflicts(vector<pair<vector<unsigned int>, vector<unsigned int> > > sets,
          unsigned int nleaves,
          vector<vector<unsigned int> >& conflicts) {

    conflicts.resize(sets.size());

    for (unsigned int i = 0; i < sets.size(); ++i) {
        sort(sets[i].first.begin(), sets[i].first.end());
        sort(sets[i].second.begin(), sets[i].second.end());
    }

    // Conflicts are possible only when smaller set has size > 1
    
    vector<unsigned int> candidates;
    for (unsigned int i = 0; i < sets.size(); ++i) {
        if (sets[i].first.size() > 1) {
            candidates.push_back(i);
        }
    }
    
    // Group things based on the value of the small set

    vector<vector<unsigned int> > groups;
    unsigned int c;
    for (unsigned int i = 0; i < candidates.size(); ++i) {
        c = candidates[i];
        const vector<unsigned int>& s = sets[c].first;
        unsigned int j;
        for (j = 0; j < groups.size(); ++j) {
            if (s == sets[groups[j].front()].first) {
                // found existing group with this set
                groups[j].push_back(c);
                break;
            }
        }
        if (j == groups.size()) {
            // no existing group; create a new one
            groups.resize(groups.size() + 1);
            groups.back().push_back(c);
        }
    }


    // Do the necessary comparisons
    unsigned int c1, c2;
    vector<unsigned int> possible_confs;
    for (unsigned int i = 0; i < groups.size(); ++i) {
        const vector<unsigned int>& s1 = sets[groups[i].front()].first;
        for (unsigned int j = i + 1; j < groups.size(); ++j) {
            const vector<unsigned int>& s2 = sets[groups[j].front()].first;
            if (Disjoint(s1, s2) || IsSubset(s1, s2) || IsSubset(s2, s1)) {
                // all members of the two groups must be consistent
                continue;
            }
            // test the between-group pairs
            possible_confs.clear();
            for (unsigned int n = 0; n < groups[i].size(); ++n) {
                c1 = groups[i][n];
                if (!IsSubset(sets[groups[j].front()].first, sets[c1].second)) {
                    // c1 not necessarily compatible with all from other group
                    possible_confs.push_back(c1);
                }
            }
            for (unsigned int m = 0; m < groups[j].size(); ++m) {
                c2 = groups[j][m];
                if (!IsSubset(sets[groups[i].front()].first, sets[c2].second)) {
                    for (unsigned int n = 0; n < possible_confs.size(); ++n) {
                        c1 = possible_confs[n];
                        if (!JointlyExhaustive(sets[c1].second,
                                               sets[c2].second,
                                               nleaves)){
                            conflicts[c1].push_back(c2);
                            conflicts[c2].push_back(c1);
                        }
                    }
                }
            }
        }
    }

    for (unsigned int i = 0; i < conflicts.size(); ++i) {
        sort(conflicts[i].begin(), conflicts[i].end());
    }

}

bool did_je = false;

// This assumes that first argument has ambiguities,
// though if it doesn't the penalty is just in
// performance.
void MutualDisambiguate(pair<vector<unsigned int>, vector<unsigned int> >& s1,
                        pair<vector<unsigned int>, vector<unsigned int> >& s2,
                        unsigned int nleaves) {
    if (Disjoint(s1.second, s2.second)) {
        // cannot conflict, however refined
        return;
    }
    bool subset = IsSubset(s1.first, s2.second);
    bool superset = IsSubset(s2.first, s1.second);
    bool disjoint = Disjoint(s1.first, s2.first);
    bool joint_exh = JointlyExhaustive(s1.second, s2.second, nleaves);
    if (subset + superset + disjoint + joint_exh > 1) {
        // cannot disambiguate anything if more than one of these holds
        return;
    }
    if (!(subset || superset || disjoint || joint_exh)) {
        // not compatible
        throw CInconsisExcept(s1, s2);
    }

    if (subset) {
        if (!IsSubset(s1.second, s2.second)) {
            vector<unsigned int> tmp;
            tmp.reserve(s1.second.size());
            set_intersection(s1.second.begin(), s1.second.end(),
                             s2.second.begin(), s2.second.end(),
                             back_inserter(tmp));
            s1.second.swap(tmp);
        }
        if (s2.second.size() != s2.first.size()
            && !IsSubset(s1.first, s2.first)) {
            vector<unsigned int> tmp;
            tmp.reserve(s2.second.size());
            set_union(s2.first.begin(), s2.first.end(),
                      s1.first.begin(), s1.first.end(),
                      back_inserter(tmp));
            s2.first.swap(tmp);
        }
    }
    if (superset) {
        if (!IsSubset(s2.first, s1.first)) {
            vector<unsigned int> tmp;
            tmp.reserve(s1.second.size());
            set_union(s1.first.begin(), s1.first.end(),
                      s2.first.begin(), s2.first.end(),
                      back_inserter(tmp));
            s1.first.swap(tmp);
        }
        if (s2.second.size() != s2.first.size()
            && !IsSubset(s2.second, s1.second)) {
            vector<unsigned int> tmp;
            tmp.reserve(s2.second.size());
            set_intersection(s2.second.begin(), s2.second.end(),
                             s1.second.begin(), s1.second.end(),
                             back_inserter(tmp));
            s2.second.swap(tmp);
        }
    }
    if (disjoint) {
        if (!Disjoint(s1.second, s2.first)) {
            vector<unsigned int> tmp;
            tmp.reserve(s1.second.size());
            set_difference(s1.second.begin(), s1.second.end(),
                           s2.first.begin(), s2.first.end(),
                           back_inserter(tmp));
            s1.second.swap(tmp);
        }
        if (s2.second.size() != s2.first.size()
            && !Disjoint(s2.second, s1.first)) {
            vector<unsigned int> tmp;
            tmp.reserve(s2.second.size());
            set_difference(s2.second.begin(), s2.second.end(),
                           s1.first.begin(), s1.first.end(),
                           back_inserter(tmp));
            s2.second.swap(tmp);
        }
    }
    if (joint_exh) {
        if (!JointlyExhaustive(s1.first, s2.second, nleaves)) {
            if (!did_je) {
                did_je = true;
            }
            vector<unsigned int> tmp1;
            vector<unsigned int> tmp2;
            tmp1.reserve(s1.first.size());
            set_difference(s1.second.begin(), s1.second.end(),
                           s2.second.begin(), s2.second.end(),
                           back_inserter(tmp1));
            set_union(s1.first.begin(), s1.first.end(),
                      tmp1.begin(), tmp1.end(),
                      back_inserter(tmp2));
            s1.first.swap(tmp2);
        }
        if (s2.second.size() != s2.first.size()
            && !JointlyExhaustive(s2.first, s1.second, nleaves)){
            vector<unsigned int> tmp1;
            vector<unsigned int> tmp2;
            tmp1.reserve(s1.first.size());
            set_difference(s2.second.begin(), s2.second.end(),
                           s1.second.begin(), s1.second.end(),
                           back_inserter(tmp1));
            set_union(s2.first.begin(), s2.first.end(),
                      tmp1.begin(), tmp1.end(),
                      back_inserter(tmp2));
            s2.first.swap(tmp2);
        }
    }
}


vector<pair<vector<unsigned int>, vector<unsigned int> > >
Disambiguate(vector<pair<vector<unsigned int>, vector<unsigned int> > > sets,
             unsigned int nleaves) {

    for (unsigned int i = 0; i < sets.size(); ++i) {
        sort(sets[i].first.begin(), sets[i].first.end());
        sort(sets[i].second.begin(), sets[i].second.end());
    }

    int old_diff_count = 0;  // For keeping track of whether
    int new_diff_count = 0;  // anything has changed during a pass
    for (unsigned int pass = 0;
         pass == 0 || new_diff_count != old_diff_count;
         ++pass) {
        for (unsigned int i = 0; i < sets.size(); ++i) {
            if (sets[i].first.size() == sets[i].second.size()) {
                continue;
            }
            for (unsigned int j = 0; j < sets.size(); ++j) {
                try {
                    MutualDisambiguate(sets[i], sets[j], nleaves);
                } catch (CInconsisExcept &e) {
                    cerr << "disambiguating incompatible sets" << endl;
                    cerr << i << " " << j << endl;
                    e.SetIdx1(i);
                    e.SetIdx2(j);
                    throw;
                }
                if (sets[i].first.size() == sets[i].second.size()) {
                    // It's now completely unambiguous
                    break;
                }
            }
        }
        // Compute the total count of ambiguous bits.
        // A pass should either decrease this or
        // leave it unchanged.  In the latter case,
        // nothing has changed and we're done.
        old_diff_count = new_diff_count;
        new_diff_count = 0;
        for (unsigned int i = 0; i < sets.size(); ++i) {
            new_diff_count += sets[i].second.size() - sets[i].first.size();
        }
    }

    // If leaf 0 disambiguated for a set, complement if appropriate
    vector<unsigned int> everything(nleaves);
    for (unsigned int i = 0; i < nleaves; ++i) {
        everything[i] = i;
    }
    vector<unsigned int> tmp;
    for (unsigned int i = 0; i < sets.size(); ++i) {
        if (sets[i].first.front() == 0) {
            // complement lower bound
            tmp.clear();
            tmp.reserve(nleaves - sets[i].first.size());
            set_difference(everything.begin(), everything.end(),
                           sets[i].first.begin(), sets[i].first.end(),
                           back_inserter(tmp));
            sets[i].first.swap(tmp);
            // and upper bound
            tmp.clear();
            tmp.reserve(nleaves - sets[i].second.size());
            set_difference(everything.begin(), everything.end(),
                           sets[i].second.begin(), sets[i].second.end(),
                           back_inserter(tmp));
            sets[i].second.swap(tmp);
            // interchange upper and lower bounds
            sets[i].second.swap(sets[i].first);
        }
    }

    // Check that all pairs are still consistent
    for (unsigned int i = 0; i < sets.size(); ++i) {
        if (sets[i].first.size() == 1) {
            continue;  // Singletons are consistent with everything
        }
        for (unsigned int j = i + 1; j < sets.size(); ++j) {
            if (!Compatible(sets[i], sets[j], nleaves)) {
                throw CInconsisExcept(i, j, sets[i], sets[j]);
            }
        }
    }

    return sets;
}


// For cases where all but a few are already disambiguated
vector<pair<vector<unsigned int>, vector<unsigned int> > >
Disambiguate(vector<pair<vector<unsigned int>, vector<unsigned int> > > sets,
             unsigned int nleaves,
             const vector<unsigned>& not_disambig) {

    for (unsigned int i = 0; i < sets.size(); ++i) {
        sort(sets[i].first.begin(), sets[i].first.end());
        sort(sets[i].second.begin(), sets[i].second.end());
    }

    vector<unsigned> to_process(not_disambig);

    // Size diff between UB and LB for detection of changes
    vector<unsigned> old_diff_count;

    while (true) {
        vector<unsigned> changed;
        old_diff_count.resize(to_process.size());
        for (unsigned int i = 0; i < to_process.size(); ++i) {
            old_diff_count[i] = sets[to_process[i]].second.size()
                - sets[to_process[i]].first.size();
        }
        for (unsigned int i = 0; i < sets.size(); ++i) {
            unsigned diff = sets[i].second.size() - sets[i].first.size();
            for (auto j : to_process) {
                try {
                    MutualDisambiguate(sets[j], sets[i], nleaves);
                } catch (CInconsisExcept &e) {
                    cerr << "disambiguating incompatible sets" << endl;
                    cerr << j << " " << i << endl;
                    e.SetIdx1(j);
                    e.SetIdx2(i);
                    throw;
                }
            }
            if (sets[i].second.size() - sets[i].first.size() != diff) {
                changed.push_back(i);
            }
        }

        for (unsigned int i = 0; i < to_process.size(); ++i) {
            if (sets[to_process[i]].second.size()
                - sets[to_process[i]].first.size()
                != old_diff_count[i]) {
                if (find(begin(changed), end(changed),
                         to_process[i]) == end(changed)) {
                    changed.push_back(to_process[i]);
                }
            }
        }

        // Are we done?
        if (changed.empty()) {
            break;
        }

        to_process = changed;
    }

    // If leaf 0 disambiguated for a set, complement if appropriate
    vector<unsigned int> everything(nleaves);
    for (unsigned int i = 0; i < nleaves; ++i) {
        everything[i] = i;
    }
    vector<unsigned int> tmp;
    for (unsigned int i = 0; i < sets.size(); ++i) {
        if (sets[i].first.front() == 0) {
            // complement lower bound
            tmp.clear();
            tmp.reserve(nleaves - sets[i].first.size());
            set_difference(everything.begin(), everything.end(),
                           sets[i].first.begin(), sets[i].first.end(),
                           back_inserter(tmp));
            sets[i].first.swap(tmp);
            // and upper bound
            tmp.clear();
            tmp.reserve(nleaves - sets[i].second.size());
            set_difference(everything.begin(), everything.end(),
                           sets[i].second.begin(), sets[i].second.end(),
                           back_inserter(tmp));
            sets[i].second.swap(tmp);
            // interchange upper and lower bounds
            sets[i].second.swap(sets[i].first);
        }
    }

    return sets;
}


vector<pair<vector<unsigned int>, vector<unsigned int> > >
Disambiguate(vector<pair<vector<unsigned int>, vector<unsigned int> > > sets,
             unsigned int nleaves,
             const vector<unsigned>& ambig_nodes,
             const vector<vector<unsigned> >& possible_conflicts) {
    vector<pair<unsigned int, unsigned int> > conflicts;
    for (unsigned int i = 0; i < sets.size(); ++i) {
        sort(sets[i].first.begin(), sets[i].first.end());
        sort(sets[i].second.begin(), sets[i].second.end());
    }

    int old_diff_count = 0;  // For keeping track of whether
    int new_diff_count = 0;  // anything has changed during a pass
    for (unsigned int pass = 0;
         pass == 0 || new_diff_count != old_diff_count;
         ++pass) {
        for (unsigned n = 0; n < ambig_nodes.size(); ++n) {
            unsigned i = ambig_nodes[n];
            if (sets[i].first.size() == sets[i].second.size()) {
                continue;
            }
            for (auto j : possible_conflicts[n]) {
                try {
                    MutualDisambiguate(sets[i], sets[j], nleaves);
                } catch (CInconsisExcept &e) {
                    cerr << "disambiguating incompatible sets" << endl;
                    cerr << i << " " << j << endl;
                    e.SetIdx1(i);
                    e.SetIdx2(j);
                    throw;
                }
                if (sets[i].first.size() == sets[i].second.size()) {
                    // It's now completely unambiguous
                    break;
                }
            }
        }
        // Compute the total count of ambiguous bits.
        // A pass should either decrease this or
        // leave it unchanged.  In the latter case,
        // nothing has changed and we're done.
        old_diff_count = new_diff_count;
        new_diff_count = 0;
        for (auto i : ambig_nodes) {
            new_diff_count += sets[i].second.size() - sets[i].first.size();
        }
    }

    // If leaf 0 disambiguated for a set, complement if appropriate
    vector<unsigned int> everything(nleaves);
    for (unsigned int i = 0; i < nleaves; ++i) {
        everything[i] = i;
    }
    vector<unsigned int> tmp;
    for (auto i : ambig_nodes) {
        if (sets[i].first.front() == 0) {
            // complement lower bound
            tmp.clear();
            tmp.reserve(nleaves - sets[i].first.size());
            set_difference(everything.begin(), everything.end(),
                           sets[i].first.begin(), sets[i].first.end(),
                           back_inserter(tmp));
            sets[i].first.swap(tmp);
            // and upper bound
            tmp.clear();
            tmp.reserve(nleaves - sets[i].second.size());
            set_difference(everything.begin(), everything.end(),
                           sets[i].second.begin(), sets[i].second.end(),
                           back_inserter(tmp));
            sets[i].second.swap(tmp);
            // interchange upper and lower bounds
            sets[i].second.swap(sets[i].first);
        }
    }

    // Check that all pairs are still consistent
    for (unsigned n = 0; n < ambig_nodes.size(); ++n) {
        unsigned i = ambig_nodes[n];
        if (sets[i].first.size() == 1) {
            continue;  // Singletons are consistent with everything
        }
        for (auto j : possible_conflicts[n]) {
            if (!Compatible(sets[i], sets[j], nleaves)) {
                throw CInconsisExcept(i, j, sets[i], sets[j]);
            }
        }
    }

    return sets;
}


vector<pair<vector<unsigned int>, vector<unsigned int> > >
ForceDisambiguation(vector<pair<vector<unsigned int>, vector<unsigned int> > > sets,
                    unsigned int nleaves) {
    // Force and disambiguate
    unsigned int idx;
    vector<unsigned int> to_process(1);
    for (unsigned int i = 0; i < sets.size(); ++i) {
        while (sets[i].first.size() != sets[i].second.size()) {
            // Remove one not in LB from UB
            if (sets[i].second.size() == sets[i].first.size() + 1) {
                sets[i].second = sets[i].first;
            } else {
                // Find first in UB but not in LB
                for (idx = 0; ; ++idx) {
                    if (sets[i].second[idx] != sets[i].first[idx]) {
                        break;
                    }
                }
                // Remove it
                for (unsigned int j = idx; j < sets[i].second.size() - 1; ++j) {
                    sets[i].second[j] = sets[i].second[j+1];
                }
                sets[i].second.resize(sets[i].second.size() - 1);
            }
            to_process[0] = i;
            sets = Disambiguate(sets, nleaves, to_process);
        }
    }
    return sets;
}


TRanges PairwiseDisambiguateMore(const TRanges& sets1,
                                 const TRanges& sets2,
                                 unsigned int nleaves) {
    TRanges rv;
    pair<vector<unsigned int>, vector<unsigned int> > tmp;
    for (auto s1 : sets1) {
        for (auto s2 : sets2) {
            bool subset = IsSubset(s1.first, s2.second);
            bool superset = IsSubset(s2.first, s1.second);
            bool disjoint = Disjoint(s1.first, s2.first);
            bool joint_exh = JointlyExhaustive(s1.second, s2.second, nleaves);
            if (subset) {
                tmp.first = s1.first;
                if (!IsSubset(s1.second, s2.second)) {
                    tmp.second.clear();
                    tmp.second.reserve(s1.second.size());
                    set_intersection(s1.second.begin(), s1.second.end(),
                                     s2.second.begin(), s2.second.end(),
                                     back_inserter(tmp.second));
                } else {
                    tmp.second = s1.second;
                }
                rv.push_back(tmp);
            }
            if (superset) {
                tmp.second = s1.second;
                if (!IsSubset(s2.first, s1.first)) {
                    tmp.first.clear();
                    tmp.first.reserve(s1.first.size());
                    set_union(s1.first.begin(), s1.first.end(),
                              s2.first.begin(), s2.first.end(),
                              back_inserter(tmp.first));
                } else {
                    tmp.first = s1.first;
                }
                rv.push_back(tmp);
            }
            if (disjoint) {
                tmp.first = s1.first;
                if (!Disjoint(s1.second, s2.first)) {
                    tmp.second.clear();
                    tmp.second.reserve(s1.second.size());
                    set_difference(s1.second.begin(), s1.second.end(),
                                   s2.first.begin(), s2.first.end(),
                                   back_inserter(tmp.second));
                } else {
                    tmp.second = s1.second;
                }
                rv.push_back(tmp);
            }
            if (joint_exh) {
                tmp.second = s1.second;
                if (!JointlyExhaustive(s1.first, s2.second, nleaves)) {
                    if (!did_je) {
                        did_je = true;
                    }
                    vector<unsigned int> tmp1;
                    tmp1.reserve(s1.first.size());
                    set_difference(s1.second.begin(), s1.second.end(),
                                   s2.second.begin(), s2.second.end(),
                                   back_inserter(tmp1));
                    tmp.first.clear();
                    set_union(s1.first.begin(), s1.first.end(),
                              tmp1.begin(), tmp1.end(),
                              back_inserter(tmp.first));
                } else {
                    tmp.first = s1.first;
                }
                rv.push_back(tmp);
            }
        }
    }
    return rv;
}


unsigned int ColorUpperBound(const vector<unsigned int>& nodes,
                             const vector<unsigned int>& colors,
                             const vector<unsigned int>& weights) {

    vector<unsigned int> maxes(colors.size());

    unsigned int sum = 0;

    for (unsigned int i = 0; i < nodes.size(); ++i) {
        if (weights[nodes[i]] > maxes[colors[nodes[i]]]) {
            sum += weights[nodes[i]] - maxes[colors[nodes[i]]];
            maxes[colors[nodes[i]]] = weights[nodes[i]];
        }
    }

    return sum;

}


unsigned int ColorUpperBound(const vector<unsigned int>& nodes,
                             const vector<unsigned int>& colors,
                             const vector<unsigned int>& weights,
                             unsigned ncolors) {

    vector<unsigned int> maxes(ncolors);

    unsigned int sum = 0;

    for (auto n : nodes) {
        if (weights[n] > maxes[colors[n]]) {
            sum += weights[n] - maxes[colors[n]];
            maxes[colors[n]] = weights[n];
        }
    }

    return sum;

}


class CCmpCounts {
public:
    CCmpCounts(const vector<unsigned>& counts) : m_Counts(counts) {};
    bool operator ()(unsigned lhs, unsigned rhs) {
        return m_Counts[lhs] > m_Counts[rhs];
    };
private:
    const vector<unsigned>& m_Counts;
};

vector<unsigned> Color(const vector<vector<char> >& conn) {
    
    // Count compatibilities (unweighted)
    vector<unsigned> compats(conn.size());
    unsigned count;
    for (unsigned i = 0; i < conn.size(); ++i) {
        count = 0;
        for (unsigned j = 0; j < conn.size(); ++j) {
            count += conn[i][j];
        }
        compats[i] = count;
    }
    
    // Sort by number of compatibilities, largest first
    vector<unsigned> nodes(conn.size());
    for (unsigned i = 0; i < conn.size(); ++i) {
        nodes[i] = i;
    }
    stable_sort(nodes.begin(), nodes.end(), CCmpCounts(compats));

    // Assign colors in order
    vector<vector<unsigned> > by_color(1);
    vector<unsigned> colors(conn.size());
    by_color[0].push_back(nodes[0]);  // First in order gets first color
    unsigned n, c, j;
    for (unsigned i = 1; i < nodes.size(); ++i) {
        n = nodes[i];
        for (c = 0; c < by_color.size(); ++c) {
            for (j = 0; j < by_color[c].size(); ++j) {
                if (conn[n][by_color[c][j]]) {
                    break;  // adjacent, so can't get same color
                }
            }
            if (j == by_color[c].size()) {
                // Didn't break, so not adjacent to anything with this color
                by_color[c].push_back(n);
                colors[n] = c;
                break;
            }
        }
        if (c == by_color.size()) {
            // Cannot be given any existing color; give it a new one
            by_color.resize(by_color.size() + 1);
            by_color.back().push_back(n);
            colors[n] = c;
        }
    }

    return colors;
    
}


vector<unsigned> Color2(const vector<vector<char> >& conn,
                        const vector<unsigned int>& weights) {
    
    // Count compatibilities (unweighted)
    vector<unsigned> compats(conn.size());
    unsigned count;
    for (unsigned i = 0; i < conn.size(); ++i) {
        count = 0;
        for (unsigned j = 0; j < conn.size(); ++j) {
            count += conn[i][j];
        }
        compats[i] = count;
    }
    
    // Sort by number of compatibilities, largest first
    vector<unsigned> nodes(conn.size());
    for (unsigned i = 0; i < conn.size(); ++i) {
        nodes[i] = i;
    }
    stable_sort(nodes.begin(), nodes.end(), CCmpCounts(compats));

    // Sort by weight, largest first
    stable_sort(nodes.begin(), nodes.end(), CCmpCounts(weights));

    // Assign colors in order
    vector<vector<unsigned> > by_color(1);
    vector<unsigned> colors(conn.size());
    by_color[0].push_back(nodes[0]);  // First in order gets first color
    unsigned n, c, j;
    for (unsigned i = 1; i < nodes.size(); ++i) {
        n = nodes[i];
        for (c = 0; c < by_color.size(); ++c) {
            for (j = 0; j < by_color[c].size(); ++j) {
                if (conn[n][by_color[c][j]]) {
                    break;  // adjacent, so can't get same color
                }
            }
            if (j == by_color[c].size()) {
                // Didn't break, so not adjacent to anything with this color
                by_color[c].push_back(n);
                colors[n] = c;
                break;
            }
        }
        if (c == by_color.size()) {
            // Cannot be given any existing color; give it a new one
            by_color.resize(by_color.size() + 1);
            by_color.back().push_back(n);
            colors[n] = c;
        }
    }

    return colors;
    
}


vector<unsigned> Color2(vector<unsigned> nodes,
                        const vector<vector<char> >& conn,
                        const vector<unsigned int>& weights) {
    
    // Count compatibilities (unweighted)
    vector<unsigned> compats(conn.size());
    unsigned count;
    for (auto i : nodes) {
        count = 0;
        for (auto j : nodes) {
            count += conn[i][j];
        }
        compats[i] = count;
    }
    
    // Sort by number of compatibilities, largest first
    stable_sort(nodes.begin(), nodes.end(), CCmpCounts(compats));

    // Sort by weight, largest first
    stable_sort(nodes.begin(), nodes.end(), CCmpCounts(weights));

    // Assign colors in order
    vector<vector<unsigned> > by_color(1);
    vector<unsigned> colors(conn.size());
    by_color[0].push_back(nodes[0]);  // First in order gets first color
    unsigned n, c, j;
    for (unsigned i = 1; i < nodes.size(); ++i) {
        n = nodes[i];
        for (c = 0; c < by_color.size(); ++c) {
            for (j = 0; j < by_color[c].size(); ++j) {
                if (conn[n][by_color[c][j]]) {
                    break;  // adjacent, so can't get same color
                }
            }
            if (j == by_color[c].size()) {
                // Didn't break, so not adjacent to anything with this color
                by_color[c].push_back(n);
                colors[n] = c;
                break;
            }
        }
        if (c == by_color.size()) {
            // Cannot be given any existing color; give it a new one
            by_color.resize(by_color.size() + 1);
            by_color.back().push_back(n);
            colors[n] = c;
        }
    }

    return colors;
    
}


vector<unsigned> ColorReduce(const vector<unsigned>& s,
                             const vector<unsigned>& colors,
                             const vector<unsigned>& weights,
                             const vector<vector<char> >& conn,
                             unsigned min_weight) {
    vector<unsigned> good;
    good.reserve(s.size());
    vector<unsigned> intersection;
    for (unsigned i = 0; i < s.size(); ++i) {
        intersection.clear();
        for (unsigned j = 0; j < s.size(); ++j) {
            if (conn[s[i]][s[j]]) {
                intersection.push_back(s[j]);
            }
        }
        if (ColorUpperBound(intersection, colors, weights) >= min_weight) {
            good.push_back(s[i]);
        }
    }
    return good;
}


vector<unsigned> Greedy(set<unsigned> nodes,
                        const vector<unsigned>& weights,
                        const vector<set<unsigned> >& conflicts) {

    vector<unsigned> total_weights(conflicts.size());
    unsigned i, conf;
    for (set<unsigned>::iterator pos = nodes.begin();
         pos != nodes.end();
         ++pos) {
        i = *pos;
        unsigned sum = 0;
        for (set<unsigned>::iterator cpos = conflicts[i].begin();
             cpos != conflicts[i].end();
             ++cpos) {
            conf = *cpos;
            if (nodes.count(conf)) {  // only count conflicts in our node set
                sum += weights[conf];
            }
        }
        total_weights[i] = sum;
    }

    vector<unsigned> removed_nodes;
    unsigned mx;
    while (true) {
        if (nodes.size() == 0) {
            break;
        }

        // Find the max of total conflict weight among current nodes
        vector<unsigned> wts;
        for (set<unsigned>::iterator pos = nodes.begin();
             pos != nodes.end();
             ++pos){
            wts.push_back(total_weights[*pos]);
        }
        mx = *max_element(wts.begin(), wts.end());
        if (mx == 0) {
            // We're down to a clique
            break;
        }

        // Find all current nodes with this max
        vector<unsigned> possibly_remove;
        for (set<unsigned>::iterator pos = nodes.begin(); pos != nodes.end(); ++pos){
            if (total_weights[*pos] == mx) {
                possibly_remove.push_back(*pos);
            }
        }

        // Remove all with identical pattern of conflicts in current set.
        // This makes them adjacent in the ordering.
        vector<unsigned> s;
        set_intersection(nodes.begin(),
                         nodes.end(),
                         conflicts[possibly_remove[0]].begin(),
                         conflicts[possibly_remove[0]].end(),
                         back_inserter(s));
        vector<unsigned> to_remove;
        to_remove.reserve(possibly_remove.size());
        for (vector<unsigned>::iterator pos = possibly_remove.begin();
             pos != possibly_remove.end();
             ++pos) {
            if (includes(conflicts[*pos].begin(),  // Check that s is a subset,
                         conflicts[*pos].end(),    // which is sufficient given
                         s.begin(),                // equality of wts
                         s.end())) {
                to_remove.push_back(*pos);
            }
        }
        unsigned w = 0;  // Total weight of removed nodes
        for (vector<unsigned>::iterator pos = to_remove.begin();
             pos != to_remove.end();
             ++pos) {
            nodes.erase(*pos);
            w += weights[*pos];
            removed_nodes.push_back(*pos);
        }
        // Adjust total weights of conflicts accordingly
        for (set<unsigned>::iterator conf = conflicts[to_remove[0]].begin();
             conf != conflicts[to_remove[0]].end();
             ++conf) {
            total_weights[*conf] -= w;
        }
    }

    return removed_nodes;

}


template<class TValue, class TKey> void GroupBy(const vector<TValue>& values,
                                                const vector<TKey>& keys,
                                                map<TKey, vector<TValue> >& groups) {
    for (unsigned i = 0; i < values.size(); ++i) {
        groups[keys[values[i]]].push_back(values[i]);
    }
}


// Comparison function for sorting by first element
template<class T>
inline bool CompFront(const T& lhs, const T& rhs) {
    return lhs.front() < rhs.front();
}


unsigned Preprocess(const vector<unsigned>& nodes,
                    const vector<unsigned>& weights,
                    const vector<vector<unsigned> >& conflicts,
                    vector<unsigned>& included,
                    vector<vector<unsigned> >& groups_new,
                    vector<unsigned>& weights_new,
                    vector<vector<char> >& conn,
                    vector<unsigned>& ranges) {

    vector<unsigned> inc_nodes;  // nodes with incompatibilities
    for (auto n : nodes) {
        if (!conflicts[n].empty()) {
            inc_nodes.push_back(n);
        }
    }

    map<vector<unsigned>, vector<unsigned> > group_map;
    GroupBy(inc_nodes, conflicts, group_map);

    vector<vector<unsigned> > groups;
    // map node to group num
    vector<unsigned> group_for_node(conflicts.size(), conflicts.size());
    for (const auto& pos : group_map) {
        groups.push_back(pos.second);
        for (unsigned n : pos.second) {
            group_for_node[n] = groups.size() - 1;
        }
    }

    vector<set<unsigned> > group_conflicts(groups.size());
    for (unsigned i = 0; i < groups.size(); ++i) {
        const auto& confs = conflicts[groups[i][0]];
        for (auto conf : confs) {
            if (group_for_node[conf] != conflicts.size()) {
                group_conflicts[i].insert(group_for_node[conf]);
            }
        }
    }

    set<unsigned> group_indices;
    for (unsigned i = 0; i < groups.size(); ++i) {
        group_indices.insert(i);
    }

    vector<unsigned> group_weights(groups.size());
    for (unsigned i = 0; i < groups.size(); ++i) {
        for (unsigned j = 0; j < groups[i].size(); ++j) {
            group_weights[i] += weights[groups[i][j]];
        }
    }
    cerr << group_weights.size() << " groups" << endl;
    vector<unsigned> removed_groups = 
        Greedy(group_indices, group_weights, group_conflicts);

    vector<unsigned> removed_nodes;
    for (auto i : removed_groups) {
        copy(groups[i].begin(), groups[i].end(), back_inserter(removed_nodes));
    }

    // Reorder accordingly
    set<unsigned> sremoved_nodes;
    for (auto i : removed_nodes) {
        sremoved_nodes.insert(i);
    }
    vector<unsigned> tmp(removed_nodes);
    for (auto i : inc_nodes) {
        if (!sremoved_nodes.count(i)) {
            tmp.push_back(i);
        }
    }
    inc_nodes.swap(tmp);

    // Remove those with too many incompatibilities
    unsigned max_incompat = 0;
    for (auto i : removed_nodes) {
        max_incompat += weights[i];
    }
    vector<unsigned> bad_groups;
    set<unsigned> bad_nodes;
    for (auto g : removed_groups) {
        if (Weight(group_conflicts[g], group_weights) > max_incompat) {
            bad_groups.push_back(g);
            for (auto n : groups[g]) {
                bad_nodes.insert(n);
            }
        }
    }
    tmp.clear();
    for (auto i : inc_nodes) {
        if (!bad_nodes.count(i)) {
            tmp.push_back(i);
        }
    }
    inc_nodes.swap(tmp);

    // Set aside those with no incompatibilities in smaller set
    vector<unsigned> candidates;
    set<unsigned> scandidates;
    for (auto i : inc_nodes) {
        if (!includes(bad_nodes.begin(), bad_nodes.end(),
                      conflicts[i].begin(), conflicts[i].end())) {
            candidates.push_back(i);
            scandidates.insert(i);
        }
    }

    // These are definitely included in any maximum clique
    for (auto i : nodes) {
        if (!scandidates.count(i) && !bad_nodes.count(i)) {
            included.push_back(i);
        }
    }

    // Group candidates with same conflict patterns among candidates

    map<unsigned, unsigned> mapping;
    for (unsigned i = 0; i < candidates.size(); ++i) {
        mapping[candidates[i]] = i;
    }

    vector<vector<unsigned> > conflicts2;
    conflicts2.reserve(candidates.size());
    for (auto i : candidates) {
        conflicts2.resize(conflicts2.size() + 1);
        set_intersection(conflicts[i].begin(), conflicts[i].end(),
                         scandidates.begin(), scandidates.end(),
                         back_inserter(conflicts2.back()));
    }

    map<vector<unsigned>, vector<unsigned> > d;
    vector<unsigned> r;
    r.reserve(conflicts2.size());
    for (unsigned i = 0; i < conflicts2.size(); ++i) {
        r.push_back(i);
    }
    GroupBy(r, conflicts2, d);

    groups.clear();
    for (auto& pair : d) {
        groups.push_back(pair.second);
    }
    sort(groups.begin(),
         groups.end(),
         CompFront<vector<unsigned> >);     // Order as in original candidates

    vector<unsigned> reps;                  // a representative of each group
    set<unsigned> sreps;
    weights_new.reserve(groups.size());
    for (auto& group : groups) {
        reps.push_back(candidates[group.front()]);
        sreps.insert(candidates[group.front()]);
        groups_new.resize(groups_new.size() + 1);
        for (auto i : group) {
            groups_new.back().push_back(candidates[i]);
        }
        weights_new.push_back(Weight(groups_new.back(), weights));
    }
    map<unsigned, unsigned> pos_for_rep;
    for (unsigned i = 0; i < reps.size(); ++i) {
        pos_for_rep[reps[i]] = i;
    }

    // Connection matrix
    
    conn.clear();
    conn.resize(groups.size(), vector<char>(groups.size(), 1));

    for (unsigned i = 0; i < groups.size(); ++i) {
        for (auto j : conflicts2[mapping[reps[i]]]) {
            if (sreps.count(j)) {
                conn[i][pos_for_rep[j]] = 0;
            }
        }
    }

    // Ranges with same conn pattern among remaining nodes
    vector<unsigned> same;
    if (conn.size() > 0) {
        for (unsigned i = 0; i < conn.size() - 1; ++i) {
            if (equal(conn[i].begin() + i,
                      conn[i].end(),
                      conn[i + 1].begin() + i)) {
                same.push_back(i);
            }
        }
    }

    ranges.reserve(groups.size());
    for (unsigned i = 0; i < conn.size(); ++i) {
        ranges.push_back(i);
    }

    for (auto iter = same.rbegin(); iter != same.rend(); ++iter) {
        ranges[*iter] = ranges[*iter + 1];
    }

    unsigned min_weight = Weight(nodes, weights)
        - Weight(removed_nodes, weights)
        - Weight(included, weights);

    return min_weight;
}


// Faster version of Python zip
vector<vector<char> > Zip(const vector<string >& rows) {
    unsigned nrows = rows.size();
    unsigned ncols = rows[0].size();
    vector<vector<char> > cols(ncols);
    for (unsigned i = 0; i < ncols; ++i) {
        cols[i].resize(nrows);
    }
    for (unsigned c = 0; c < ncols; ++c) {
        for (unsigned r = 0; r < nrows; ++r) {
            cols[c][r] = rows[r][c];
        }
    }
    return cols;
}


// Faster version of Python zip
void Zip(const vector<string >& rows, vector<vector<char> >& cols) {
    unsigned nrows = rows.size();
    unsigned ncols = rows[0].size();
    cols.resize(ncols);
    for (unsigned c = 0; c < ncols; ++c) {
        cols[c].resize(nrows);
        for (unsigned r = 0; r < nrows; ++r) {
            cols[c][r] = rows[r][c];
        }
    }
}


SDynProgEl::SDynProgEl(const vector<unsigned>& the_nodes) {
    nodes = the_nodes;
    scores.resize(the_nodes.size() + 1);
}


unsigned int DynProgUpperBound(const vector<unsigned int>& nodes,
                               SDynProg& dp,
                               const vector<vector<char> >& conn,
                               const vector<unsigned int>& weights) {
    // Represent set of nodes differently (allows fast look-up)
    vector<char> node_set(dp.max_node + 1);
    for (auto n : nodes) {
        node_set[n] = true;
    }
    // Do the recursion
    DynProgRecurse(node_set, dp.root, conn, weights);
    // Compute return value from result
    return *max_element(dp.root->scores.begin(), dp.root->scores.end());
}


void DynProgRecurse(const vector<char>& node_set,
                    SDynProgEl *dp,
                    const vector<vector<char> >& conn,
                    const vector<unsigned int>& weights) {
    
    for (auto ch : dp->children) {
        DynProgRecurse(node_set, ch, conn, weights);
    }
    
    unsigned mx = 0;
    unsigned n = 0;
    unsigned total;

    // Scores for taking each node from this category
    for (unsigned i = 0; i < dp->nodes.size(); ++i) {
        n = dp->nodes[i];
        if (!node_set[n]) {
            dp->scores[i] = 0;
        } else {
            total = weights[n];
            for (auto ch : dp->children) {
                mx = 0;  // max of allowed scores from this child
                for (unsigned j = 0; j < ch->nodes.size(); ++j) {
                    if (ch->scores[j] > mx && conn[n][ch->nodes[j]]) {
                        mx = ch->scores[j];
                    }
                }
                if (ch->scores.back() > mx) {
                    mx = ch->scores.back();
                }
                total += mx;
            }
            dp->scores[i] = total;
        }
    }

    // Score if we take nothing from this category
    total = 0;
    for (auto ch : dp->children) {
        mx = 0;  // max of allowed scores from this child
        for (unsigned j = 0; j < ch->nodes.size(); ++j) {
            if (ch->scores[j] > mx) {
                mx = ch->scores[j];
            }
        }
        if (ch->scores.back() > mx) {
            mx = ch->scores.back();
        }
        total += mx;
    }
    dp->scores.back() = total;

}


void DoNothing(vector<pair<vector<unsigned int>, vector<unsigned int> > > sets) {}
