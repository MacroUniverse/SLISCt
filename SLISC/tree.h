// tree utility
#pragma once
#include "search.h"

namespace slisc {

struct Node
{
    vector<Long> last; // last nodes
    vector<Long> next; // next nodes
};

// links[2n] -> links[2n+1]
void tree_gen(vector<Node> &tree, const vector<Long> &links)
{
    Long Nlink = links.size();
    tree.resize(Nlink);
    for (Long i = 0; i < Nlink; i += 2) {
        Long ind1 = links[i], ind2 = links[i + 1];
        tree[ind1].next.push_back(ind2);
        tree[ind2].last.push_back(ind1);
    }
}

// iterative implementation of tree_all_dep()
void tree_all_dep_imp(vector<Long> &deps, const vector<Node> &tree, Long_I ind)
{
    for (Long i = 0; i < Size(tree[ind].last); ++i) {
        Long ind0 = tree[ind].last[i];
        deps.push_back(ind0);
        tree_all_dep_imp(deps, tree, ind0);
    }
}

// find all upstream nodes of a tree, and the distances
Long tree_all_dep(vector<Long> &deps, const vector<Node> &tree, Long_I ind)
{
    deps.clear();
    tree_all_dep_imp(deps, tree, ind);
    while (true) {
        Long ret = find_repeat(deps);
        if (ret < 0) {
            return Size(deps);
        }
        deps.erase(deps.begin() + ret);
    }
}

// find redundant i.e. if A->B->...C, then A->C is redundent
Long tree_redundant(vector<Long> &links, const vector<Node> &tree)
{
    vector<Long> deps;
    for (Long i = 0; i < Size(tree); ++i) {
        deps.clear();
        tree_all_dep_imp(deps, tree, i);
        for (Long j = 0; j < Size(tree[i].last); ++j) {
            Long ind = search(tree[i].last[j], deps);
            if (ind >= 0) {
                deps.erase(deps.begin() + ind);
            }
        }
        for (Long j = 0; j < Size(tree[i].last); ++j) {
            Long ind = search(tree[i].last[j], deps);
            if (ind >= 0) {
                links.push_back(deps[ind]);
                links.push_back(i);
            }
        }
    }
    return Size(links) / 2;
}

} // namespace slisc
