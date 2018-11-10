import os, subprocess, sys

if len(sys.argv) != 4:
	print('Usage: python loco_kinship.py plink_basename output_dir/ path/to/pylmmKinship.py')
	print('Given plink tped/tfam/map files, generates LOCO tpeds/maps and kinship matrices.')
	sys.exit()
plink_basename, outdir, pylmmKinship = sys.argv[1:]
tpedname, tfamfile, mapname = plink_basename+'.tped', plink_basename+'.tfam', plink_basename+'.map'

# figure out number of chromosomes in tped file
num_chroms = 0
with(open(tpedname, 'r')) as tped:
	for line in tped:
		chrom = int(line.split()[0])
		if chrom > num_chroms:
			num_chroms = chrom

# make directories and open file handlers
if not(outdir.endswith('/')):
	outdir += '/'
if not os.path.isdir(outdir):
	os.makedirs(outdir)
# LOCO files for LOCO kinship analysis, one_chr files for subsequent GWAS
loco_fnames = [outdir + 'loco_chr' + str(i+1) for i in range(num_chroms)]
onechr_fnames = [outdir + 'one_chr' + str(i+1) for i in range(num_chroms)]
tped_locos = [open(loco_fnames[i] + '.tped', 'w') for i in range(num_chroms)]
map_locos = [open(loco_fnames[i] + '.map', 'w') for i in range(num_chroms)]
tped_onechrs = [open(onechr_fnames[i] + '.tped', 'w') for i in range(num_chroms)]
map_onechrs = [open(onechr_fnames[i] + '.map', 'w') for i in range(num_chroms)]

# write full tped file to tped LOCO files, leaving appropriate chromosome out
with(open(tpedname, 'r')) as tped:
	with(open(mapname, 'r')) as map:
		for line in tped:
			mapline = map.readline()
			chrom = int(line.split()[0])
			for i in range(len(tped_locos)):
				if i+1 != chrom:
					tped_locos[i].write(line)
					map_locos[i].write(mapline)
				else:
					tped_onechrs[i].write(line)
					map_onechrs[i].write(mapline)
for i in range(len(tped_locos)):
	tped_locos[i].close()
	map_locos[i].close()
	tped_onechrs[i].close()
	map_onechrs[i].close()
for fname in loco_fnames:  # tfam file doesn't change by chromosome
	subprocess.Popen(['cp', tfamfile, fname + '.tfam']).wait()
for fname in onechr_fnames:
	subprocess.Popen(['cp', tfamfile, fname + '.tfam']).wait()

# run pylmmKinship for each LOCO file
for i in range(num_chroms):
	print('Running pylmm on ' + loco_fnames[i])
	subprocess.Popen(['python', pylmmKinship, '--tfile', loco_fnames[i],  loco_fnames[i] + '.kin']).wait()
#

