{\rtf1\ansi\ansicpg1252\cocoartf2708
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fmodern\fcharset0 CourierNewPSMT;\f1\fswiss\fcharset0 ArialMT;}
{\colortbl;\red255\green255\blue255;\red0\green0\blue255;\red255\green255\blue255;\red26\green26\blue26;
}
{\*\expandedcolortbl;;\cssrgb\c0\c0\c100000;\cssrgb\c100000\c100000\c100000;\cssrgb\c13333\c13333\c13333;
}
\paperw11900\paperh16840\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\deftab720
\pard\pardeftab720\partightenfactor0

\f0\fs20 \cf2 \cb3 \expnd0\expndtw0\kerning0
seqfile = stewart.aa * sequence data file name\cb1 \
\cb3 outfile = mlc * main result file name\cb1 \
\cb3 treefile = stewart.trees * tree structure file name\cb1 \
\cb3 noisy = 9 * 0,1,2,3,9: how much rubbish on the screen\cb1 \
\cb3 verbose = 0 * 1: detailed output, 0: concise output\cb1 \
\cb3 runmode = 0 * 0: user tree; 1: semi-automatic; 2: automatic\cb1 \
\cb3 * 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise\cb1 \
\cb3 seqtype = 2 * 1:codons; 2:AAs; 3:codons-->AAs\cb1 \
\cb3 CodonFreq = 2 * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table\cb1 \
\cb3 * ndata = 10\cb1 \
\cb3 clock = 0 * 0:no clock, 1:clock; 2:local clock; 3:TipDate\cb1 \
\cb3 aaDist = 0 * 0:equal, +:geometric; -:linear, 1-6:G1974,Miyata,c,p,v,a\cb1 \
\cb3 * 7:AAClasses\cb1 \
\cb3 aaRatefile = wag.dat * only used for aa seqs with model=empirical(_F)\cb1 \
\cb3 * dayhoff.dat, jones.dat, wag.dat, mtmam.dat, or your own\cb1 \
\cb3 model = 2\cb1 \
\cb3 * models for codons:\cb1 \
\cb3 * 0:one, 1:b, 2:2 or more dN/dS ratios for branches\cb1 \
\cb3 * models for AAs or codon-translated AAs:\cb1 \
\cb3 * 0:poisson, 1:proportional,2:Empirical,3:Empirical+F\cb1 \
\cb3 * 6:FromCodon, 8:REVaa_0, 9:REVaa(nr=189)\cb1 \
\cb3 NSsites = 0 * 0:one w;1:neutral;2:selection; 3:discrete;4:freqs;\cb1 \
\cb3 * 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;\cb1 \
\cb3 * 10:beta&gamma+1; 11:beta&normal>1; 12:0&2normal>1;\cb1 \
\cb3 * 13:3normal>0\cb1 \
\cb3 icode = 0 * 0:universal code; 1:mammalian mt; 2-11:see below\cb1 \
\cb3 Mgene = 0 * 0:rates, 1:separate;\cb1 \
\cb3 fix_kappa = 0 * 1: kappa fixed, 0: kappa to be estimated\cb1 \
\cb3 kappa = 2 * initial or fixed kappa\cb1 \
\cb3 fix_omega = 0 * 1: omega or omega_1 fixed, 0: estimate\cb1 \
\cb3 omega = .4 * initial or fixed omega, for codons or codon-based AAs\cb1 \
\cb3 fix_alpha = 1 * 0: estimate gamma shape parameter; 1: fix it at alpha\cb1 \
\cb3 alpha = 0. * initial or fixed alpha, 0:infinity (constant rate)\cb1 \
\cb3 Malpha = 0 * different alphas for genes\cb1 \
\cb3 ncatG = 3 * # of categories in dG of NSsites models\cb1 \
\cb3 fix_rho = 1 * 0: estimate rho; 1: fix it at rho\cb1 \
\cb3 rho = 0. * initial or fixed rho, 0:no correlation\cb1 \
\cb3 getSE = 0 * 0: don't want them, 1: want S.E.s of estimates\cb1 \
\cb3 RateAncestor = 0 * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)
\f1\fs26\fsmilli13200 \cf4 \cb1 \

\f0\fs20 \cf2 \cb3 Small_Diff = .5e-6\cb1 \
\cb3 * cleandata = 0 * remove sites with ambiguity data (1:yes, 0:no)?\cb1 \
\cb3 * fix_blength = 0 * 0: ignore, -1: random, 1: initial, 2: fixed\cb1 \
\cb3 method = 0 * 0: simultaneous; 1: one branch at a time
\f1\fs26\fsmilli13200 \cf4 \cb1 \
}