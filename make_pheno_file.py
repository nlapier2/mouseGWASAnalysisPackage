"""
This script makes a phenotype file in pylmm forarray, which pylmm then uses to
  read in the phenotype values for each mouse for GWAS analysis. Essentially,
  we are simply extracting relevant information from the clinical traits file.
One of the plink files called a tfam is required input, as it specifies which
  mice are present and genotyped for this study.
This script transforms the phenotype by using a quantile transformation.
  This sorts the phenotypes and gives them the corresponding proporiton
  (1/n,2/n,...,n/n) and maps to the z-score. Alternatively, the user can normalizes
  phenotypes by subtracting their mean and dividing by their standard deviation.
  The user can specify no transformations, but generally shouldn't. The (default)
  quantile transformation is preferred/is more conservative.
Eventually this script will allow the user to specify covariates to
  regress out, but we currently do not have that option implemented.
"""


import argparse, math, sys, glob
import numpy as np
from sklearn import preprocessing
from scipy import stats
#copied all imports from stacked-manhattan.py may not need all
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure

def parseargs():    # handle user arguments
  parser = argparse.ArgumentParser(description="Given phenotype file," +
        " get into proper format for pylmm.")
  parser.add_argument('--clinical', required = True,
    help = 'Name of existing pheno file.')
  parser.add_argument('--target', required=True,
    help = 'EXACT name of target trait to study.')
  parser.add_argument('--tfam', required=True,
    help = 'Plink tfam file; exclude mice not in it.')
  parser.add_argument('--output', default='pheno_file_pylmm.txt',
    help = 'Name of output file.')
  parser.add_argument('--transform', default='standardize',
    choices=['quantile', 'standardize','boxcox', 'none'],
    help = 'Phenotype transformation. Options: quantile, standardize (default), boxcox, and none. The boxcox transformation also needs the --lambda argument')
  parser.add_argument('--lmbda', '--Lambda', default='0', type=float,
    help = 'Lambda value for boxcox transformation, 0 (default) which returns the log-transformed data.')
  parser.add_argument('--plot', action='store_true',
    help = 'Will return qqplot of the standardized, quantile, and boxcox for lambda= -2, -1, -0.5, 0, 0.5, and 2 in a stacked plot. This will be writen to the --output file prefix with .png ending. We will ignore  --transform and --lambda')
    #parser.add_argument('--regress', action='store_true',
  #  help = 'Use this to regress target phenotype on other phenotypes.')
  args = parser.parse_args()
  return args


def parse_tfam(args):
  """
  Reads the mice for this study from a tfam file. Used to ensure that mice
      that were excluded from analysis are also excluded from the pheno file.
  Argument: args are the user arguments parsed with argparse.
  Returns: tfams dict that maps strain to mouse number
  """
  tfams = {}
  if args.tfam != 'NONE':
    with(open(args.tfam, 'r')) as tfamfile:
      for line in tfamfile:
        splits = line.split(' ')  # splits[0, 1] == [strain, number]
        if splits[0] not in tfams:
          tfams[splits[0]] = [splits[1]]
        else:
          tfams[splits[0]].append(splits[1])
  return tfams


def parse_clinical_file(args, tfams):
  """
  Given a target phenotype, read the values for that phenotype.
  Arguments:
  -- args are the user arguments parsed with argparse
  -- tfams are the strain to mouse number mappings from the tfam file
  Returns:
  -- mouse2pheno: maps mouse to the value of the target phenotype and the
      the non-target phenos (other), which can be corrected for
  -- target_pheno: just the raw target phenotype values
  -- other_phenos: just the raw non-target (other) phenotype values
  """
  mouse2pheno, target_pheno, other_phenos = [], [], []
  with(open(args.clinical, 'r')) as infile:
    header = infile.readline().strip().split('\t')
    targetcol = header.index(args.target)  # find column of target phenotype
    mousecol, straincol = header.index('mouse_number'), header.index('Strain')
    for line in infile:
      splits = line.strip().split('\t')
      if len(splits) < 2:
        break
      if tfams != {} and (splits[straincol] not in tfams
	      or splits[mousecol] not in tfams[splits[straincol]]):
        continue  # this mouse not present in the tfam file
      mouse_fid_iid = splits[straincol] + ' ' + splits[mousecol]  # family & indiv. ID
      # grab target phenotype and other phenotyes in this line
      tgt = splits[targetcol]
      if tgt != 'NA':
        tgt = float(tgt)
      others = [splits[i] for i in range(len(splits))
        if i not in [straincol, mousecol, targetcol]]
      #others = splits[2 : targetcol] + splits[targetcol + 1 : ]
      # pack info together with mouse and also store target/other phenos
      mouse2pheno.append([mouse_fid_iid, tgt, others])
      if tgt != 'NA':
        target_pheno.append(tgt)
      other_phenos.append(others)
  return mouse2pheno, target_pheno, other_phenos


# NOT YET IMPLEMENTED
#def regress_pheno(args, target_pheno, other_phenos):
#  return target_pheno


def normalize_and_write(args, mouse2pheno, target_pheno):
  """
    By default, we preform a quantile transformation by sorting the traits and taking the
  relative location (1/n, 2/n, .., n/n) and map this to the z-scores.
  The user has the opportunity to standardize the target phenotype values by subtracting
    the mean and dividing by standard deviation
    The phenotype output is then writes results in pylmm format to args.output
  Arguments:
  -- mouse2pheno: maps mouse to the value of the target phenotype and the
      the non-target phenos (other), which can be corrected for
  -- target_pheno: just the raw target phenotype values
  -- other_phenos: just the raw non-target (other) phenotype values
  """
  if args.plot:
    #sklearn wants a matrix so here's a matrix
    mat = np.matrix(target_pheno)
    #there may be a less obnoxious way of getting everything to work but whatever
    mat = mat.reshape((len(target_pheno),1))
    #stats.boxcox wants a column vector, so here's a column vector
    array = np.array(target_pheno)
    quantile = preprocessing.quantile_transform(mat,output_distribution='normal',axis=0)
    #this is stupid but necessary
    quantile = quantile.flatten()
    saveloc = args.output.split('.')[0]
    saveloc = saveloc+'.png'
    # calculate mean and standard deviation of target phenotype values
    mean = sum(target_pheno) / float(len(target_pheno))
    sq_diffs = [(i - mean)**2 for i in target_pheno]
    sample_sd = math.sqrt(sum(sq_diffs) / float(len(target_pheno) - 1))
    standardize = [(i-mean)/sample_sd for i in target_pheno]
    boxneg2 = stats.boxcox(array, -2.0)
    boxneg1 = stats.boxcox(array, -1.0)
    boxneg05 = stats.boxcox(array,-0.5)
    box0 = stats.boxcox(array, 0.0)
    boxpos05 = stats.boxcox(array,0.5)
    boxpos2 = stats.boxcox(array,2.0)
    figure(num=None, figsize=(16,32), dpi=80, facecolor='w', edgecolor='k')
    plt.title('QQPlots of Different Transformations')
    plt.subplots_adjust(hspace=0.5)
    ax = plt.subplot(4, 2, 1)
    stats.probplot(standardize,plot=plt)
    ax.set_title('Standardized Phenotype')
    ax = plt.subplot(4, 2, 2)
    stats.probplot(quantile,plot=plt)
    ax.set_title('Quantile Transformed Phenotype')
    ax = plt.subplot(4, 2, 3)
    stats.probplot(boxneg2,plot=plt)
    ax.set_title('Box-Cox Lambda=-2 Transformed Phenotype')
    ax = plt.subplot(4, 2, 4)
    stats.probplot(boxneg1,plot=plt)
    ax.set_title('Box-Cox Lambda=-1 Transformed Phenotype')
    ax = plt.subplot(4, 2, 5)
    stats.probplot(boxneg05,plot=plt)
    ax.set_title('Box-Cox Lambda=-0.5 Transformed Phenotype')
    ax = plt.subplot(4, 2, 6)
    stats.probplot(box0,plot=plt)
    ax.set_title('Box-Cox Lambda=0 Transformed Phenotype')
    ax = plt.subplot(4, 2, 7)
    stats.probplot(boxpos05,plot=plt)
    ax.set_title('Box-Cox Lambda=0.5 Transformed Phenotype')
    ax = plt.subplot(4, 2, 8)
    stats.probplot(boxpos2,plot=plt)
    ax.set_title('Box-Cox Lambda=2 Transformed Phenotype')
    plt.savefig(saveloc)
    #just exit
    exit()

  if args.transform == 'boxcox':
    #stats.boxcox wants a column vector, so here's a column vector
    array = np.array(target_pheno)
    transform = stats.boxcox(array, float(args.lmbda))
    transform = transform.tolist()
    #target_pheno has no NAs so the loc of values don't necessarily map exactly
    targ_loc = 0
    for i in range(len(mouse2pheno)):
      if mouse2pheno[i][1] != 'NA':
        mouse2pheno[i][1] = transform[targ_loc]
        targ_loc += 1

  elif args.transform == 'quantile':
    #required to be a matrix and not an array
    mat = np.matrix(target_pheno)
    #the sklearn magic
    transform = preprocessing.quantile_transform(mat,output_distribution='normal',axis=1)
    #takes matrix and puts it into list
    transform = transform[0].tolist()
    #target_pheno has no NAs so the loc of values don't necessarily map exactly
    targ_loc = 0
    for i in range(len(mouse2pheno)):
      if mouse2pheno[i][1] != 'NA':
        mouse2pheno[i][1] = transform[targ_loc]
        targ_loc += 1

  elif args.transform == 'standardize':
    # calculate mean and standard deviation of target phenotype values
    mean = sum(target_pheno) / float(len(target_pheno))
    sq_diffs = [(i - mean)**2 for i in target_pheno]
    sample_sd = math.sqrt(sum(sq_diffs) / float(len(target_pheno) - 1))
    # normalize target phenotypes, stored in mouse2pheno[i][1]
    for i in range(len(mouse2pheno)):
      if mouse2pheno[i][1] != 'NA':
        mouse2pheno[i][1] = ((mouse2pheno[i][1] - mean) / sample_sd)

  with(open(args.output, 'w')) as outfile:
    for m in mouse2pheno:
      # m[0] is mouse strain and number, m[1] is target phenotype value
      outfile.write(m[0] + ' ' + str(m[1]) + '\n')


def main():
  args = parseargs()  # handle user arguments
  tfams = parse_tfam(args)  # read mice in this study from tfam file
  # read the phenotype values from the clinical trait file
  mouse2pheno, target_pheno, other_phenos = parse_clinical_file(args, tfams)
  #if args.regress:
  #  target_pheno = regress_pheno(args, target_pheno, other_phenos)
  # normalize phenotype values and write to args.output
  normalize_and_write(args, mouse2pheno, target_pheno)


if __name__ == '__main__':
  main()
#
