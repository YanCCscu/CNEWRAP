#!/usr/bin/env python3
from ete3 import Tree
import sys

def load_tree(treefile):
	mytree=open(treefile)
	while True:
		treestr=mytree.readline()
		if treestr != "":
			break
	return(Tree(treestr,format=1))
	
def mark_tree_node(treefile):
	t=load_tree(treefile)
	#print(t)
	for node in t.traverse("postorder"):
		if not node.is_leaf():
			#print(node.get_children())
			node.name=node.get_children()[0].name.split("-")[0]+\
			"-"+node.get_children()[1].name.split("-")[0]
		else:
				node.name=node.name.split("_")[0]
	return(t.write(format=1))

def get_tree_node(treefile):
	t=load_tree(treefile)
	allnodes=[]
	for node in t.traverse("postorder"):
		if node.is_leaf():
			allnodes.append(node.name.split("_")[0])
	return(allnodes)
	
def traverse_from_leaf(treefile, leaf_name):
	tree=load_tree(treefile)
	leaf_node = tree.search_nodes(name=leaf_name)[0]
	
	ancestors = leaf_node.get_ancestors()
	ancestors.insert(0, leaf_node)  
	visited_leaves = []  
	for ancestor in ancestors:
		leaves = ancestor.get_leaves()
		for leaf in leaves:
			if leaf.name not in visited_leaves: 
				visited_leaves.append(leaf.name)
	return(visited_leaves)

def dist_by_tree(treefile, leaf_name):
	tree=load_tree(treefile)
	leaf_node = tree.search_nodes(name=leaf_name)[0]
	ancestors = leaf_node.get_ancestors()
	ancestors.insert(0, leaf_node)  
	visited_leaves = [leaf_name]  
	with open('lastz.dist','w') as LDIST:
		for ancestor in ancestors:
			leaves = ancestor.get_leaves()
			for leaf in leaves:
				if leaf.name not in visited_leaves: 
					visited_leaves.append(leaf.name)
					print("{}\t{}\tmedium".format(leaf_name,leaf.name),file=LDIST)

if __name__=="__main__":
	allnodes=get_tree_node(sys.argv[1])
	print(allnodes)
	print(mark_tree_node(sys.argv[1]))
	nodename='Tele'
	nodename=sys.argv[2]
	print(traverse_from_leaf(sys.argv[1],nodename))
	dist_by_tree(sys.argv[1],nodename)

