#pragma once
#include "../SLISC/tree.h"

void test_tree()
{
    using namespace slisc;
    vector<Node> tree;
    vector<Long> links = {1, 4, 2, 4, 3, 5, 3, 10, 4, 7, 6, 8, 7, 8, 10, 9, 8, 9, 5, 6, 4, 6,
        1, 7, 3, 9, 2, 9, 7, 9, 5, 9};
    cout << "tree_gen()..." << endl;
    tree_gen(tree, links);
    vector<Long> deps;
    cout << "tree_all_dep()..." << endl;
    tree_all_dep(deps, tree, 9);
    /*for (Long i = 0; i < Size(deps); ++i)
        cout << deps[i] << endl;*/
    vector<Long> links_redundant;
    tree_redundant(links_redundant, tree);
    /*for (Long i = 0; i < Size(links_redundant); i += 2)
        cout << links_redundant[i] << " -> " << links_redundant[i+1] << endl;*/
}
