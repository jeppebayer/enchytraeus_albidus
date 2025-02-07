#!/bin/env python3
import sys, os.path, matplotlib.pyplot as plt, pandas as pd, numpy as np
from statsmodels.stats.weightstats import DescrStatsW

fullPath = os.path.abspath(sys.argv[1])
columnNameFile = os.path.abspath(sys.argv[2])
minThres = int(sys.argv[3])
plotDir = os.path.abspath(sys.argv[4])
plotName = os.path.basename(fullPath)
columnNames = pd.read_csv(columnNameFile).columns

nColumns = pd.read_table(fullPath, header=None, nrows=1).shape[1]
depthDistribution = pd.DataFrame()

for i, name in zip(range(2, nColumns), columnNames):
	sampleDepthDistribution = pd.DataFrame()
	depthFileIterator = pd.read_table(fullPath, header=None, iterator=True, chunksize=10000, usecols=[i])
	for chunk in depthFileIterator:
		sampleDepthDistribution = pd.concat([sampleDepthDistribution, chunk.value_counts()], axis=1, join='outer').sum(axis=1)

	depthDistribution = pd.concat([depthDistribution, sampleDepthDistribution.rename(name)], axis=1, join='outer')

# Create observation column from index
depthDistribution = depthDistribution.reset_index().convert_dtypes()
depthDistribution.rename(columns={"level_0": "Depth"}, inplace=True)
depthDistribution = depthDistribution.sort_values(by='Depth').reset_index(drop=True)
# Remove depth 0 observations
depthDistribution = depthDistribution[depthDistribution['Depth'] > 0]

# Create tsv file of summary values
with open(f'{plotDir}/{plotName}.tsv', 'w') as outfile:
	outfile.write(f'sample_name\tmean\tstd\tmax_coverage_value\tpeak\tmax_threshold\n')
	plt.figure(figsize = (40, 40))
	for index, sample in enumerate(columnNames, start=1):
		subset = depthDistribution[['Depth', sample]]
		subset = subset[subset[sample].notna()]
		subset.sort_values(by=sample, inplace=True, ascending=False)

		subsetWeightedStatsBessel = DescrStatsW(subset['Depth'], weights=subset[sample], ddof=1)
		maxCov = subset.max(axis=0)['Depth']
		mean = subsetWeightedStatsBessel.mean
		stdDev = subsetWeightedStatsBessel.std

		for row in subset.itertuples():
			if row.Depth > minThres:
				peak = row.Depth
				maxThres = peak + 2 * stdDev
				break
		subsetMaxThres = subset[subset['Depth'] <= maxThres]

		outfile.write(f'{sample}\t{mean}\t{stdDev}\t{maxCov}\t{peak}\t{maxThres}\n')

		plt.subplot(6, 4, index)
		plt.bar(x=subsetMaxThres['Depth'], height=subsetMaxThres[sample], color='r', width=1)
		plt.ylabel("Count", fontweight = 'bold')
		plt.xlabel("Coverage", fontweight = 'bold')
		plt.xticks(np.arange(0, maxThres, 50.0), rotation='vertical')
		plt.axvline(mean, color = 'black', linestyle = 'dashed', label = 'Mean')
		plt.axvline(peak, color = 'blue', linestyle = 'dashed', label = 'Peak')
		plt.text(round(mean + stdDev), plt.ylim()[1]/2, f'Maximum: {round(maxThres, ndigits=2)}x\nMean: {round(mean, ndigits=2)}x\nStd: {round(stdDev, ndigits=2)}\nPeak: {peak}x\nUpper threshold: {round(maxThres, ndigits=2)}x', va = 'center', color = 'black', bbox=dict(facecolor='white'))
		plt.legend()
		plt.grid(True)
		plt.title(f'Coverage distribution {sample}', fontweight = 'bold')

plt.tight_layout()
plt.savefig(f'{plotDir}/{plotName}.png')
plt.close()

for index, sample in enumerate(columnNames, start=1):
	subset = depthDistribution[['Depth', sample]]
	subset = subset[subset[sample].notna()]

	subset.sort_values(by=sample, inplace=True, ascending=False)

	subsetWeightedStatsBessel = DescrStatsW(subset['Depth'], weights=subset[sample], ddof=1)
	maxCov = subset.max(axis=0)['Depth']
	mean = subsetWeightedStatsBessel.mean
	stdDev = subsetWeightedStatsBessel.std

	for row in subset.itertuples():
		if row.Depth > minThres:
			peak = row.Depth
			maxThres = peak + 2 * stdDev
			break
	subsetMaxThres = subset[subset['Depth'] <= maxThres]

	plt.figure(figsize = (15, 12))
	plt.bar(x=subsetMaxThres['Depth'], height=subsetMaxThres[sample], color='r', width=1)
	plt.ylabel("Count", fontweight = 'bold')
	plt.xlabel("Coverage", fontweight = 'bold')
	plt.xticks(np.arange(0, maxThres, 50.0), rotation='vertical')
	plt.axvline(mean, color = 'black', linestyle = 'dashed', label = 'Mean')
	plt.axvline(peak, color = 'blue', linestyle = 'dashed', label = 'Peak')
	plt.text(round(mean + stdDev), plt.ylim()[1]/2, f'Maximum: {maxThres}x\nMean: {round(mean, ndigits=2)}x\nStd: {round(stdDev, ndigits=2)}\nPeak: {peak}x\nUpper threshold: {maxThres}x', va = 'center', color = 'black', bbox=dict(facecolor='white'))
	plt.legend()
	plt.grid(True)
	plt.title(f'Coverage distribution {sample}', fontweight = 'bold')

	plt.tight_layout()
	plt.savefig(f'{plotDir}/{sample}.png')
	plt.close()
