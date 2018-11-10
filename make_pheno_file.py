import argparse, math, sys
import numpy as np


def parseargs():    # handle user arguments
	parser = argparse.ArgumentParser(description="Given phenotype file," +
				" get into proper format for pylmm.")
	parser.add_argument('input', help = 'Name of existing pheno file.')
	parser.add_argument('output', help = 'Name of output file.')
	parser.add_argument('--no_normalization', action='store_true',
		help = 'Do not normalize target phenotype.')
	parser.add_argument('--regress', action='store_true',
		help = 'Use this to regress target phenotype on other phenotypes.')
	parser.add_argument('--target', type=int, default=3,
		help = 'Column number of target phenotype.')
	parser.add_argument('--tfam', default='NONE',
		help = 'Optionally provide a tfam file & exclude mice not in it.')
	args = parser.parse_args()
	return args


def parse_tfam(args):
	if args.tfam != 'NONE':
		tfams = {}
		with(open(args.tfam, 'r')) as tfamfile:
			for line in tfamfile:
				splits = line.split(' ')
				if splits[0] not in tfams:
					tfams[splits[0]] = [splits[1]]
				else:
					tfams[splits[0]].append(splits[1])
	else:
		tfams = None
	return tfams


def parse_input_file(args, tfams=None):
	mouse2pheno, target_pheno, other_phenos = [], [], []
	with(open(args.input, 'r')) as infile:
		infile.readline() ; infile.readline()  # skip header lines
		for line in infile:
			splits = line.split('  ')
			splits = [i.strip().replace(' ', '_') for i in splits if len(i) > 0]
			if len(splits) < 2:
				break
			splits[1] = splits[1].replace('/', '.')  # for plink
			if tfams != None and (splits[1] not in tfams or splits[0] not in tfams[splits[1]]):
				continue  # this mouse not present in the tfam file
			mouse_fid_iid = splits[1] + ' ' + splits[0]  # family & indiv. ID
			try:  # target phenotype must be number
				tgt = float(splits[args.target])
			except:
				print('Could not read all target phenotypes as floats.')
				print(target_pheno)
				print(splits[args.target])
				sys.exit()
			others = splits[2 : args.target] + splits[args.target + 1 : ]
			mouse2pheno.append([mouse_fid_iid, tgt, others])
			target_pheno.append(tgt)
			other_phenos.append(others)
	return mouse2pheno, target_pheno, other_phenos


def regress_pheno(args, target_pheno, other_phenos):
	return target_pheno


def normalize_and_write(args, mouse2pheno, target_pheno):
	# normalize via subtracting mean then dividing by standard deviation
	if not args.no_normalization:
		mean, sd = np.mean(target_pheno), math.sqrt(np.var(target_pheno))
		for i in range(len(mouse2pheno)):
			mouse2pheno[i][1] = ((mouse2pheno[i][1] - mean) / sd)
	with(open(args.output, 'w')) as outfile:
		for m in mouse2pheno:
			outfile.write(m[0] + ' ' + str(m[1]) + '\n')


def main():
	args = parseargs()
	tfams = parse_tfam(args)
	mouse2pheno, target_pheno, other_phenos = parse_input_file(args, tfams)
	if args.regress:
		target_pheno = regress_pheno(args, target_pheno, other_phenos)
	normalize_and_write(args, mouse2pheno, target_pheno)


if __name__ == '__main__':
	main()
#

