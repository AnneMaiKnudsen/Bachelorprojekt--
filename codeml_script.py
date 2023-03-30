from Bio.Phylo.PAML import codeml

cml = codeml.Codeml()
cml.alignment = "/home/anenaur/anne_ane_bsc/data/output_test/chr1/SCYL3/SCYL3.phylip"
cml.tree = "/home/anenaur/anne_ane_bsc/data/output_test/chr1/SCYL3/SCYL3.nw"
cml.out_file = "/home/anenaur/anne_ane_bsc/data/output_test/chr1/SCYL3/SCYL3_codeml_test.out"
cml.working_dir = "/home/anenaur/anne_ane_bsc"

cml.set_options(clock=0)
cml.set_options(NSsites=[1, 2, 7, 8])
cml.set_options(verbose=1)
cml.set_options(seqtype=1)
cml.set_options(CodonFreq = 2, cleandata = 1, fix_blength = None, fix_omega = 0, ncatG = 1, runmode = 0, fix_kappa = 0, fix_alpha = 1, Small_Diff = .5e-06, method = 0, Malpha = 0, aaDist = 0, RateAncestor = 0, icode = 0, alpha = 0, omega = 1, getSE = 0, noisy = 0, Mgene = 0, kappa = 2, model = 0, ndata = 1)
cml.print_options()

cml.ctl_file = "SCYL3_control.ctl"
cml.write_ctl_file()
cml.read_ctl_file("SCYL3_control.ctl")


cml.run()
#cml.run(command="/home/anenaur/anne_ane_bsc")

#hvordan kører man codeml i et script?
#hvordan får man den til at køre det rigtige codeml?

