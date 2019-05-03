"""
Input: a dir of newick trees that end with '.tre' and
	bootstrap replicates that end with '.trees'
Output: newick trees with with nodes collapsed that
	fall below a support threshold
"""

import os,sys
import subprocess

def collapse_nodes(DIR,best_tree,threshold):
	assert best_tree.endswith(".raxml_bs.tre"),\
		"best tree infile "+best_tree+" does not end with .tre"
	clusterID = best_tree.split(".raxml_bs")[0]
	bs_trees = DIR+clusterID+".raxml_bs.trees"
	collapsed_tree = DIR+clusterID+".collapsed"
	if not os.path.exists(collapsed_tree):
		target_tree = best_tree if DIR == "./" else DIR+best_tree
		cmd = ["sumtrees.py","-d0","-p","-F","newick",\
			"-f",str(threshold),"--no-annotations",\
			"-t",target_tree,"-o",collapsed_tree,\
			bs_trees]
		print "".join(cmd)
		p = subprocess.Popen(cmd,stdout=subprocess.PIPE)
		out = p.communicate()
		assert p.returncode == 0,"Error sumtrees"+out[0]
	return collapsed_tree
	
def main(DIR,tree_ending,threshold):
	if DIR[-1] != "/": DIR += "/"
	filecount = 0
	for i in os.listdir(DIR):
		if i.endswith(".tre"):
			filecount += 1
			collapse_nodes(DIR,i,threshold)
	assert filecount > 0, "No file ends with .tre found in "+DIR
	
if __name__ == "__main__":
	if len(sys.argv) != 4:
		print "python collapse_nodes_below_threshold.py DIR tree_ending threshold"
		print "make sure that sumtrees.py is in the path"
		sys.exit(0)
	
	DIR,tree_ending,threshold = sys.argv[1:]
	main(DIR,tree_ending,threshold)
		
