#!/usr/bin/env python
# coding=utf-8
import sys
import os
import glob
from optparse import OptionParser
import re

#提供工作路径、模型、分组名(需要和存放树文件的目录名一致)、比对好的fasta存放目录
parser = OptionParser()
parser.add_option('-w', '--work_dir', type='string', help=('required, Absolute path of work directory'))
parser.add_option('-m', '--model', type='string', help=('required, models used: branch, branch_site, clade or site'))
parser.add_option('-g', '--group', type='string', help=('required, name of forebackground (Consistent with the tree files storage directory)'))
parser.add_option('-a', '--alignment_path', default=None, type='string', help=('optional, if provide alignment path, the script will convert aligned fasta to PAML phylip format'))

usage = "Usage:python %prog -w [work_dir] -m [model] -g [group] -a [alignment_path]"

(options, args) = parser.parse_args()

numOpts = len(sys.argv)
if numOpts < 2:
    parser.print_help()
    sys.exit()

def fa2phy(fas_files,output) :
    files = os.listdir(fas_files)
    for fas in files :
        if not fas.endswith('fasta') and not fas.endswith('fas') and not fas.endswith('fa'):
            continue
        phy_name = fas.split(".")[0] + ".phy"
        with open("%s/%s"%(fas_files,fas), 'r') as fin:
            sequences = [(m.group(1), ''.join(m.group(2).split()))
            for m in re.finditer(r'(?m)^>([^ \n]+)[^\n]*([^>]*)', fin.read())]
        with open("%s/%s"%(output,phy_name), 'w') as fout:
            fout.write('%d %d\n' % (len(sequences), len(sequences[0][1])))
            for item in sequences:
                fout.write('%-20s %s\n' % item)

def codeml(path,model,phy,tree):
    models = {}
    ###branch_model
    models['Model1ratio'] = [0, 0, 4, 0, 1]
    models['Model2Ratio'] = [2, 0, 4, 0, 1]
    models['ModelFreeRatio'] = [1, 0, 4, 0, 1]
    ###site_model
    models['Model1Neutral'] = [0, 1, 10, 0, 0.4]
    models['Model2Selection'] = [0, 2, 10, 0, 0.4]
    models['Model3Discrtk2']  = [0, 3, 2, 0, 0.4]
    models['Model3Discrtk3'] = [0, 3, 3, 0, 0.4]
    models['Model7beta'] = [0, 7, 10, 0, 0.4]
    models['Model8beta'] = [0, 8, 10, 0, 0.4]
    ###branch_site_model
    #models['ModelA'] = [2, 2, 3, 0, 0.4]
    #models['ModelB'] = [2, 3, 3, 0, 0.4]
    #models['ModelC'] = [3, 2, 3, 0, 0.4]
    #models['ModelD'] = [3, 3, 3, 0, 0.4]
    models['NullbranchSite'] = [2, 2, 4, 1, 1]
    models['AltbranchSite'] = [2, 2, 4, 0, 1.5]
    ###Clade_model
    models['CladeModelC'] = [3, 2, 3, 0, 1.5]
    #models['CladeModelD'] = [3, 3, 3, 0, 1.5]
    models['NullCladeModelC'] = [0, 22, 10, 0, 1.5]
    ctl = open(path+'/'+model+'.ctl','w')
    ctl.write('''
    seqfile  = %s
    treefile = %s
    outfile  = mlc.out
    noisy = 0           * 0,1,2,3,9: how much rubbish on the screen
    verbose = 1         * 0: concise; 1: detailed, 2: too much
    runmode = 0         * 0: user tree;  1: semi-automatic;  2: automatic
                        * 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise
    seqtype = 1         * 1:codons; 2:AAs; 3:codons-->AAs
    CodonFreq = 2       * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
    aaDist = 0          * 0:equal, +:geometric; -:linear, 1-6:G1974,Miyata,c,p,v,a
    aaRatefile = wag.dat           * only used for aa seqs with model=empirical(_F)
                                    * dayhoff.dat, jones.dat, wag.dat, mtmam.dat, or your own
    model = %s          * models for codons:
                        * 0:one, 1:b, 2:2 or more dN/dS ratios for branches
                        * models for AAs or codon-translated AAs:
                        * 0:poisson, 1:proportional, 2:Empirical, 3:Empirical+F
                        * 6:FromCodon, 7:AAClasses, 8:REVaa_0, 9:REVaa(nr=189)
    NSsites = %s        * 0:one w;1:neutral;2:selection; 3:discrete;4:freqs;
                        * 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;
                        * 10:beta&gamma+1; 11:beta&normal>1; 12:0&2normal>1;
                        * 13:3normal>0
    icode = 0           * 0:universal code; 1:mammalian mt; 2-10:see below
    fix_kappa = 0       * 1: kappa fixed, 0: kappa to be estimated
    kappa = 3           * initial or fixed kappa
    fix_omega = %s          * 1: omega or omega_1 fixed, 0: estimate
    omega = %s          * initial or fixed omega, for codons or codon-based AAs
    fix_alpha = 1           * 0: estimate gamma shape parameter; 1: fix it at alpha
    alpha = 0           * initial or fixed alpha, 0:infinity (constant rate)
    Malpha = 0          * different alphas for genes
    ncatG = %s          * # of categories in dG of NSsites models
    clock = 0           * 0:no clock, 1:global clock; 2:local clock; 3:TipDate
    getSE = 0           * 0: don't want them, 1: want S.E.s of estimates
    RateAncestor = 0            * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)
    Small_Diff = .5e-6
    cleandata = 1\n
    ''' %(phy, tree, models[model][0], models[model][1], models[model][3], models[model][4], models[model][2]))
    ctl.close()

def ModelSelect(model_select):
    branch_model = ['Model1ratio','Model2Ratio','ModelFreeRatio']
    site_model = ['Model1Neutral','Model2Selection','Model3Discrtk2','Model3Discrtk3','Model7beta','Model8beta']
    branch_site_model = ['NullbranchSite','AltbranchSite']
    clade_model = ['CladeModelC','NullCladeModelC']
    if model_select == "branch" :
        return branch_model
    if model_select == "site" :
        return site_model
    if model_select == "branch_site" :
        return branch_site_model
    if model_select == "clade" :
        return clade_model

work_path = options.work_dir
model_select = options.model
input_files = work_path + "/input"
fasta_dir = options.alignment_path
group = options.group
alignment_path = input_files+"/alignment_files"

fout = open("run.PAML.%s.%s.sh"%(group, model_select), "w")

if fasta_dir != None :
    alignment_path = input_files+"/alignment_files"
    if os.path.exists(alignment_path) :
        pass
    else :
        os.makedirs(alignment_path)
    fa2phy(fasta_dir,alignment_path)
    if os.path.exists(group + "/run_paml") :
        pass
    else :
        os.makedirs(group + "/run_paml")
    alignments = os.listdir(alignment_path)
    for phy in alignments :
        phy_name = phy.strip()
        phy_path = input_files+ "/alignment_files/" + phy_name
        phy_id = phy_name.split(".")[0]
        tree_name = phy_id+".nwk"
        if model_select == "clade" :
            for model in ModelSelect(model_select) :
                work_dir = group + "/run_paml/" + phy_id + "/" + model_select + "_model/" + model
                tree_path = input_files+"/tree_files/Phyllostomidae_clade/" + phy_id + ".nwk"
                os.system('mkdir -p %s && cp %s %s && cp %s %s' % (work_dir,tree_path,work_dir,phy_path,work_dir))
                codeml(work_dir,model,phy_name,tree_name)
                fout.write("cd %s/%s && codeml %s.ctl"%(work_path, work_dir, model)+"\n")
        else :
            for model in ModelSelect(model_select) :
                work_dir = group + "/run_paml/" + phy_id + "/" + model_select + "_model/" + model
                tree_path = input_files + "/tree_files/" + group + "/" + phy_id + ".nwk"
                os.system('mkdir -p %s && cp %s %s && cp %s %s' % (work_dir,tree_path,work_dir,phy_path,work_dir))
                codeml(work_dir,model,phy_name,tree_name)
                fout.write("cd %s/%s && codeml %s.ctl"%(work_path, work_dir, model)+"\n")
else :
    alignment_path = input_files+"/alignment_files"
    if os.path.exists(group + "/run_paml") :
        pass
    else :
        os.makedirs(group + "/run_paml")
    alignments = os.listdir(alignment_path)
    for phy in alignments :
        phy_name = phy.strip()
        phy_path = input_files+ "/alignment_files/" + phy_name
        phy_id = phy_name.split(".")[0]
        tree_name = phy_id+".nwk"
        if model_select == "clade" :
            for model in ModelSelect(model_select) :
                work_dir = group + "/run_paml/" + phy_id + "/" + model_select + "_model/" + model
                tree_path = input_files+"/tree_files/Phyllostomidae_clade/" + phy_id + ".nwk"
                os.system('mkdir -p %s && cp %s %s && cp %s %s' % (work_dir,tree_path,work_dir,phy_path,work_dir))
                codeml(work_dir,model,phy_name,tree_name)
                fout.write("cd %s/%s && codeml %s.ctl"%(work_path, work_dir, model)+"\n")
        else :
            for model in ModelSelect(model_select) :
                work_dir = group + "/run_paml/" + phy_id + "/" + model_select + "_model/" + model
                tree_path = input_files + "/tree_files/" + group + "/" + phy_id + ".nwk"
                os.system('mkdir -p %s && cp %s %s && cp %s %s' % (work_dir,tree_path,work_dir,phy_path,work_dir))
                codeml(work_dir,model,phy_name,tree_name)
                fout.write("cd %s/%s && codeml %s.ctl"%(work_path, work_dir, model)+"\n")

print("run command next:\nnohup ParaFly -c run.PAML.%s.%s.sh -CPU 60 &"%(group, model_select))
