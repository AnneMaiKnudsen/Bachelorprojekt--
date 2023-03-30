# from Bio.Phylo.PAML import codeml

# cml = codeml.Codeml()
# cml.alignment = "/home/annemai/anne_ane_bsc/data/output_test/chr1/SCYL3/SCYL3.phylip"
# cml.tree = "/home/annemai/anne_ane_bsc/data/output_test/chr1/SCYL3/SCYL3.nw"
# cml.out_file = "/home/annemai/anne_ane_bsc/data/output_test/results.out"
# cml.working_dir = "/home/annemai/anne_ane_bsc/data/output_test"

# cml.set_options(noisy=4)
# cml.set_options(verbose=True)
# cml.set_options(NSsites=[0])
# cml.set_options(runmode=0)
# cml.set_options(seqtype=1)
# cml.set_options(omega=0.4)
# cml.set_options(kappa=2)
# cml.set_options(fix_alpha=1)
# cml.set_options(CodonFreq=2)

# cml.run("/home/annemai/anne_ane_bsc/codeml")

from Bio.Phylo.PAML import codeml

# Set the input and output file names
tree_file = "/home/annemai/anne_ane_bsc/data/output_test/chr1/SCYL3/SCYL3.nw"
alignment_file ="/home/annemai/anne_ane_bsc/data/output_test/chr1/SCYL3/SCYL3.phylip"
output_file = "/home/annemai/anne_ane_bsc/data/output_test"
control_file = "codeml.ctl"

# Set up the codeml object with the file names and working directory
cml = codeml.Codeml(alignment=alignment_file, tree=tree_file, out_file=output_file, working_dir="./")

# Set the options for the control file
cml.set_options(seqtype=1, CodonFreq=2, clock=0, aaDist=0, model=0, NSsites=[0], codeml="/home/annemai/anne_ane_bsc/codeml")

# Write the control file
with open(control_file, "w") as ctl:
    ctl.write(str(cml))

# Run codeml
cml.run()

cml.run()






