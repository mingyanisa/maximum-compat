/*  $Id: compat.hpp 5413 2024-01-06 19:27:10Z jcherry $
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
#include <utility>
#include <stdexcept>
#include <algorithm>

using namespace std;


struct SDynProgEl {
    SDynProgEl(const vector<unsigned>& nodes);
    vector<unsigned> nodes;
    vector<unsigned> scores;
    vector<SDynProgEl *> children;
};

struct SDynProg {
    unsigned max_node;  // Largest node number
    SDynProgEl *root;   // Root of tree that represents dp
};


class MaxClique {
public:
    enum EUpperBound {
        e_None,
        e_Coloring,
        e_Dyn_prog
    };
    MaxClique();
    vector<unsigned int> Run(const vector<unsigned int>& nodes,
                             const vector<unsigned int>& weights,
                             const vector<vector<char> >& conn,
                             unsigned int lower_bound=0);
    vector<unsigned int>* Recurse(vector<unsigned int>& nodes,
                                  const vector<unsigned int>& weights,
                                  const vector<vector<char> >& conn,
                                  unsigned int weight,
                                  unsigned int nodes_weight);
    vector<unsigned int> Run2(const vector<unsigned int>& nodes,
                              const vector<unsigned int>& weights,
                              const vector<vector<char> >& conn,
                              const vector<unsigned int>& ranges,
                              unsigned int lower_bound=0,
                              const vector<unsigned int>& regroup_at
                              = vector<unsigned int>());
    vector<unsigned int>* Recurse2(vector<unsigned int>& nodes,
                                   vector<unsigned int>& weights,
                                   const vector<vector<char> >& conn,
                                   const vector<unsigned int>& ranges,
                                   unsigned int weight,
                                   unsigned int nodes_weight);
    vector<vector<unsigned int> > FindAll(const vector<unsigned int>& nodes,
                                          const vector<unsigned int>& weights,
                                          const vector<vector<char> >& conn,
                                          const vector<unsigned int>& ranges,
                                          unsigned int lower_bound=0,
                                          const vector<unsigned int>& regroup_at
                                          = vector<unsigned int>());
    vector<vector<unsigned int> >*
    RecurseFA(vector<unsigned int>& nodes,
              vector<unsigned int>& weights,
              const vector<vector<char> >& conn,
              const vector<unsigned int>& ranges,
              unsigned int weight,
              unsigned int nodes_weight);
    vector<vector<unsigned int> >*
    RecurseFAColors(vector<unsigned int>& nodes,
                    vector<unsigned int>& weights,
                    const vector<vector<char> >& conn,
                    const vector<unsigned int>& ranges,
                    unsigned int weight,
                    unsigned int nodes_weight);
    bool
    HasCliqueOfWeight(const vector<unsigned int>& nodes,
                      const vector<unsigned int>& weights,
                      const vector<vector<char> >& conn,
                      unsigned int lower_bound);
    bool
    RecurseHCoW(vector<unsigned int>& nodes,
                const vector<unsigned int>& weights,
                const vector<vector<char> >& conn,
                unsigned int weight,
                unsigned int nodes_weight);

    unsigned int GetMax() {return mx;};
    void SetReportMaxes(bool b) {m_Report_maxes = b;};
    void SetReportProgress(bool b) {m_Report_progress = b;};
    void SetUpperBoundType(EUpperBound b) {m_Upper_bound = b;};
    void SetColorings(const vector<vector<unsigned> >& colorings);
    void SetDynProgs(const vector<SDynProg>& dps) {m_Dyn_progs = dps;};
    unsigned int intersections;
    vector<unsigned> intersections_hist;
    unsigned int pointless_intersections;
    unsigned int ub_calls;
private:
    unsigned int mx;
    vector<vector<unsigned int> > m_Groups;
    vector<unsigned int> m_RegroupAt;
    bool m_Report_maxes;
    bool m_Report_progress;
    vector<vector<unsigned> > m_Colorings;
    vector<unsigned> m_NColors;
    EUpperBound m_Upper_bound;
    vector<SDynProg> m_Dyn_progs;
};


bool Compatible(const pair<vector<unsigned int>, vector<unsigned int> >& s1,
                const pair<vector<unsigned int>, vector<unsigned int> >& s2,
                unsigned int nleaves);

void
Conflicts(vector<pair<vector<unsigned int>, vector<unsigned int> > > sets,
          unsigned int nleaves,
          vector<vector<unsigned int> >& result);

vector<pair<vector<unsigned int>, vector<unsigned int> > >
Disambiguate(vector<pair<vector<unsigned int>, vector<unsigned int> > > sets,
             unsigned int nleaves);

vector<pair<vector<unsigned int>, vector<unsigned int> > >
Disambiguate(vector<pair<vector<unsigned int>, vector<unsigned int> > > sets,
             unsigned int nleaves,
             const vector<unsigned>& ambig_nodes,
             const vector<vector<unsigned> >& possible_conflicts);

// For cases where all but a few are already disambiguated
vector<pair<vector<unsigned int>, vector<unsigned int> > >
Disambiguate(vector<pair<vector<unsigned int>, vector<unsigned int> > > sets,
             unsigned int nleaves,
             const vector<unsigned>& not_disambig);

void MutualDisambiguate(pair<vector<unsigned int>, vector<unsigned int> >& s1,
                        pair<vector<unsigned int>, vector<unsigned int> >& s2,
                        unsigned int nleaves);

bool IsSubset(const vector<unsigned int>& s1,
              const vector<unsigned int>& s2);
bool Disjoint(const vector<unsigned int>& s1,
              const vector<unsigned int>& s2);
bool JointlyExhaustive(const vector<unsigned int>& s1,
                       const vector<unsigned int>& s2,
                       unsigned int nleaves);

bool CannotConflict(const pair<vector<unsigned int>, vector<unsigned int> >& s1,
                    const pair<vector<unsigned int>, vector<unsigned int> >& s2,
                    unsigned int nleaves);

// Identify all that cannot conflict with any ambiguous, however resolved
vector<unsigned int>
CannotConflict(vector<pair<vector<unsigned int>, vector<unsigned int> > > sets,
               unsigned int nleaves);

vector<pair<vector<unsigned int>, vector<unsigned int> > >
ForceDisambiguation(vector<pair<vector<unsigned int>, vector<unsigned int> > > sets,
                    unsigned int nleaves);

unsigned int ColorUpperBound(const vector<unsigned int>& nodes,
                             const vector<unsigned int>& colors,
                             const vector<unsigned int>& weights);

// This version takes the number of colors as argument to save a little time
unsigned int ColorUpperBound(const vector<unsigned int>& nodes,
                             const vector<unsigned int>& colors,
                             const vector<unsigned int>& weights,
                             unsigned ncolors);

vector<unsigned> ColorReduce(const vector<unsigned>& s,
                             const vector<unsigned>& colors,
                             const vector<unsigned>& weights,
                             const vector<vector<char> >& conn,
                             unsigned min_weight);

vector<unsigned> Color(const vector<vector<char> >& conn);

vector<unsigned> Color2(const vector<vector<char> >& conn,
                        const vector<unsigned int>& weights);

vector<unsigned> Color2(vector<unsigned> nodes,
                        const vector<vector<char> >& conn,
                        const vector<unsigned int>& weights);

unsigned int DynProgUpperBound(const vector<unsigned int>& nodes,
                               SDynProg& dp,
                               const vector<vector<char> >& conn,
                               const vector<unsigned int>& weights);

void DynProgRecurse(const vector<char>& nodes,
                    SDynProgEl *dp,
                    const vector<vector<char> >& conn,
                    const vector<unsigned int>& weights);

vector<unsigned> Greedy(set<unsigned> nodes,
                        const vector<unsigned>& weights,
                        const vector<set<unsigned> >& conflicts);

unsigned Preprocess(const vector<unsigned>& nodes,
                    const vector<unsigned>& weights,
                    const vector<vector<unsigned> >& conflicts,
                    vector<unsigned>& included,
                    vector<vector<unsigned> >& groups_new,
                    vector<unsigned>& weights_new,
                    vector<vector<char> >& conn,
                    vector<unsigned>& ranges);



void DoNothing(vector<pair<vector<unsigned int>, vector<unsigned int> > > sets);

// Union of two "sets" of unsigned int, represented as sorted vectors
inline vector<unsigned int> Union(const vector<unsigned int>& first,
                                  const vector<unsigned int>& second) {
    vector<unsigned int> rv;
    set_union(first.begin(), first.end(),
              second.begin(), second.end(),
              back_inserter(rv));
    return rv;
}


// Add conflicts for all pairs of elements of indices
inline void AddConflicts(vector<vector<unsigned int> >& conflicts,
                         const vector<unsigned int>& indices) {
    for (auto i : indices) {
        vector<unsigned int> tmp = Union(conflicts[i], indices);
        conflicts[i].swap(tmp);
    }
}

inline
vector<unsigned int>
ReturnVector(unsigned int sz) {
    vector<unsigned int> rv;
    for (unsigned int i = 0; i < sz; ++i) {
        rv.push_back(i);
    }
    return rv;
}

typedef vector<pair<vector<unsigned int>, vector<unsigned int> > > TRanges;
TRanges PairwiseDisambiguateMore(const TRanges& sets1,
                                 const TRanges& sets2,
                                 unsigned int nleaves);

class CInconsisExcept : public std::invalid_argument
{
public:
    CInconsisExcept(const pair<vector<unsigned int>,
                               vector<unsigned int> >& arg1,
                    const pair<vector<unsigned int>,
                               vector<unsigned int> >& arg2) 
        : invalid_argument("disambiguating inconsistent sets"),
          m_Arg1(arg1),
          m_Arg2(arg2) {};
    CInconsisExcept(unsigned int idx1,
                    unsigned int idx2,
                    const pair<vector<unsigned int>,
                               vector<unsigned int> >& arg1,
                    const pair<vector<unsigned int>,
                               vector<unsigned int> >& arg2) 
        : invalid_argument("disambiguating inconsistent sets"),
          m_Arg1(arg1),
          m_Arg2(arg2),
          m_Idx1(idx1),
          m_Idx2(idx2)
    {};
    void SetIdx1(unsigned int idx) {m_Idx1 = idx;}
    void SetIdx2(unsigned int idx) {m_Idx2 = idx;}
    unsigned int GetIdx1(void) const {return m_Idx1;}
    unsigned int GetIdx2(void) const {return m_Idx2;}
    const pair<vector<unsigned int>, vector<unsigned int> >&
    GetArg1() const {return m_Arg1;}
    const pair<vector<unsigned int>, vector<unsigned int> >&
    GetArg2() const {return m_Arg2;}

private:
    pair<vector<unsigned int>, vector<unsigned int> > m_Arg1;
    pair<vector<unsigned int>, vector<unsigned int> > m_Arg2;
    unsigned int m_Idx1;
    unsigned int m_Idx2;
};

vector<vector<char> > Zip(const vector<string >& rows);
void Zip(const vector<string >& rows, vector<vector<char> >& cols);
