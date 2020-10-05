from treetime import TreeAnc
from Bio import Phylo

T = Phylo.read('tree_raw.nwk', 'newick')
T.root_at_midpoint()

tt = TreeAnc(tree=T, aln='subsampled_alignment.fasta')

tt.optimize_tree(prune_short=True)

to_prune = []
for n in tt.tree.get_terminals():
    if len([x for x in n.mutations if x[2]!='N']) > 5:
        to_prune.append(n)

for n in to_prune:
    tt.tree.prune(n)


Phylo.write(tt.tree, 'tree_raw.nwk', 'newick')

