#!/bin/env python3
import random
import sys
from ete3 import Tree

def generate_random_tree(leaves):
    """
    Generate a random tree with the given leaves.

    Parameters:
        leaves (list): A list of leaf node names.

    Returns:
        Tree: A Tree object representing the generated random tree.
    """
    # Shuffle the leaves to randomize the tree
    random.shuffle(leaves)

    def build_tree(nodes):
        if len(nodes) == 1:
            return nodes[0]  # Base case: single node is returned as leaf

        # Randomly divide nodes into two groups
        split = random.randint(1, len(nodes) - 1)
        left_group = nodes[:split]
        right_group = nodes[split:]

        # Recursively build subtrees
        left_subtree = build_tree(left_group)
        right_subtree = build_tree(right_group)

        # Return the current tree as a nested structure
        return f"({left_subtree},{right_subtree})"

    # Generate the random tree structure in Newick format
    newick_tree = build_tree(leaves)
    
    # Create the Tree object from the Newick string
    return Tree(newick_tree+";")

if __name__ == "__main__":
    # Get leaves from command line arguments
    if len(sys.argv) < 2:
        print("Usage: python script.py leaf1 leaf2 leaf3 ...")
        sys.exit(1)

    leaves = sys.argv[1:]

    # Generate the random tree
    random_tree = generate_random_tree(leaves)

    # Save the tree to a Newick file
    random_tree.write(outfile="random_tree.nw", format=0)

    print("Generated Random Tree saved as random_tree.nw")
