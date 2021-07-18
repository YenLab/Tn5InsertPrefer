# -*- coding: utf-8 -*-
Description = '''
This scripts can calculate the single-base sites enrichment in given genomic features.
NOTE: This is initially developed for Tn5 cut sites enrichment.

For v3, the backgroud of each location was decided by MPE-seq on naked DNA
'''.lstrip()

Copyright = '''
***
Modified by Houyu Zhang.
Issue report on Hughiez047@gmail.com
Copyright (c) 2021 __YenLab@SCUT__. All rights reserved.
'''.lstrip()

import sys,os
import glob
import argparse
import pybedtools
import scipy.stats as stats
import statistics
import numpy as np

def run(options):

	#read sites file and calculate average cut frequency
	site_base = os.path.basename(options.Sites)
	Sites = pybedtools.BedTool(options.Sites)
	total_sites = Sites.count()

	#open result file and write header
	OUT = open(options.OUT,"w+")
	OUT.write("Feature\tStatistics\tP-value\n")

	#handle each feature file
	for feature_file in glob.glob(os.path.join(options.Feature,"*bed")):
		Feature = pybedtools.BedTool(feature_file)
		feature_base = os.path.basename(feature_file)

		print("Calculating Overlap between {} sites of {} and {} features of {} ...".format(total_sites,site_base,Feature.count(),feature_base))
		Feature_with_Sites = Feature.intersect(Sites, c=True, wa=True, sorted=False)

		observed, expected = [], []
		with open(Feature_with_Sites.fn) as feature:
			for i in feature:
				line = i.strip().split("\t")
				#length = int(line[2]) - int(line[1])
				expected.append(float(line[6]) * total_sites/1000000)
				observed.append(int(line[7])) #bed format 7 cols


		#Remove outliers
		observed, expected = zip(*sorted(zip(observed, expected)))
		L = len(observed)
		outliers_lower = int(L * options.remove)
		outliers_upper = L - int(L * options.remove)
		observed = observed[outliers_lower:outliers_upper]
		expected = expected[outliers_lower:outliers_upper]

		ratio = round(statistics.mean(observed)/statistics.mean(expected),4)
		ii = np.where(np.array(expected) == 0)[0]
		expected_r0 = np.delete(np.array(expected),ii)
		observed_r0 = np.delete(np.array(observed),ii)

		print("Starting chi-square test...")
		res = stats.chisquare(f_obs = observed_r0, f_exp = expected_r0)
		pvalue = "{0:.2e}".format(res.pvalue)

		OUT.write(feature_base.rstrip(".bed") + "\t" + str(ratio) + "\t" + pvalue + "\n")

	OUT.close()

def main():
	parser = argparse.ArgumentParser(prog = os.path.basename(sys.argv[0]), usage = 'python %(prog)s [options]',
									 description = Description, epilog = Copyright,
									 formatter_class = argparse.RawDescriptionHelpFormatter)

	parser.add_argument('-g', '--GenomeSize',action = 'store', type=int, dest = 'GenomeSize',required=True,
						help = 'Total genomesize in bp [ce11 = 100286401;\n mm10 = 2725537669; hg38 = 3088286401]')
	parser.add_argument('-f', '--Feature',action = 'store', type=str, dest = 'Feature', required=True,
						help = 'Directory of bed formatted feature files')
	parser.add_argument('-s', '--Sites',action = 'store', type=str, dest = 'Sites', required=True,
						help = 'Sites file in bed format')
	parser.add_argument('-r', '--remove', action = 'store', type=float, dest = 'remove',
						help = 'Remove outliers from bins')
	parser.add_argument('-o', '--OUT',action = 'store', type=str, dest = 'OUT',
						help = 'Destination of result file')

	options = parser.parse_args()

	if not options:
		parser.print_help()
		sys.exit(3)

	run(options)

if __name__ == '__main__':
	main()
