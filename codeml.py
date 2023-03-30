from Bio.Phylo.PAML import codeml
import os

cml = codeml.Codeml()
cml.alignment = "data/output_test/chr1/SCYL3/SCYL3.phylip"
cml.tree = "data/output_test/chr1/SCYL3/SCYL3.nw"
cml.out_file = "SCYL3_codeml_test.out"
cml.working_dir = "/home/annemai/anne_ane_bsc"

cml.set_options(clock=0)
cml.set_options(NSsites=[1, 2, 7, 8])
cml.set_options(verbose=1)
cml.set_options(seqtype=1)
cml.set_options(CodonFreq = 2, cleandata = 1, fix_blength = None, fix_omega = 0, ncatG = 1, runmode = 0, fix_kappa = 0, fix_alpha = 1, Small_Diff = .5e-06, method = 0, Malpha = 0, aaDist = 0, RateAncestor = 0, icode = 0, alpha = 0, omega = 1, getSE = 0, noisy = 0, Mgene = 0, kappa = 2, model = 0, ndata = 1)
cml.print_options()

cml.ctl_file = "SCYL3_control.ctl"
cml.write_ctl_file()
cml.read_ctl_file("SCYL3_control.ctl")


#cml.run()
#cml.run(command="/home/annemai/anne_ane_bsc/codeml")
cml.run(verbose = True, command=os.path.abspath('codeml'))

# from Bio.Phylo.PAML import codeml
# from Bio.Phylo.PAML.chi2 import cdf_chi2
# import argparse
# import sys
# import os
# from pprint import pprint

# parser = argparse.ArgumentParser(description='''
# Run codeml using biopython's wrapper.

# python codeml.py ~/scratch/genes/DYNLT3.phylib ~/scratch/genes/DYNLT3.nw DYNLT3.txt DYNLT3.ctl ~/scratch

# ''')

# parser.add_argument('alignment_file', type=str,
#                     help='Alignment file.')
# parser.add_argument('tree_file', type=str,
#                     help='Tree file.')
# parser.add_argument('results_file', type=str,
#                     help='Results file.')
# parser.add_argument('control_file', type=str,
#                     help='Control file to write model options.')
# parser.add_argument('working_directory', type=str,
#                     help='Working directory for temporary files.')

# args = parser.parse_args()

# ## SITE MODEL.
# cml = codeml.Codeml(alignment = args.alignment_file, tree = args.tree_file,
#                 out_file =  args.results_file,
#                 working_dir = args.working_directory)

# cml.set_options(verbose=1)
# cml.set_options(seqtype=1) # * 1:codon models; 2: amino acid models
# #cml.set_options(runmode=0) # 0: user tree; 1: semi-automatic; 2: automatic 3: StepwiseAddition; (4,5):PerturbationNNI
# cml.set_options(model=0) #  0:one-w for all branches; 2: w’s for branches
# cml.set_options(NSsites=[0, 1, 2, 7, 8]) #  0:one-ratio; 1:neutral; 2:selection; 3:discrete; 7:beta; 8:beta&w

# cml.ctl_file = args.control_file
# #cml.write_ctl_file()

# cml.run(verbose = True, command=os.path.abspath('paml4.9j/src/codeml'))


# python scripts/codeml.py 
# codeml/steps/cds_data/chrX/IL13RA2/IL13RA2.phylib
# codeml/steps/cds_data/chrX/IL13RA2/IL13RA2.nw
# steps/codeml/IL13RA2/IL13RA2.tmp
# IL13RA2.ctl steps/codeml/IL13RA2 && mv steps/codeml/IL13RA2/IL13RA2.tmp steps/codeml/IL13RA2/IL13RA2.txt