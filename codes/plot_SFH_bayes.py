import matplotlib
matplotlib.use('agg')
import numpy as np
from scipy.io import readsav
import matplotlib.pyplot as plt
import cPickle

type_galaxy='real' #Choose from real or fake

#name_graphs='SSP_1Gyr_Calzetti'
#name_folder='SSP_1Gyr_Calzetti/'

name_graphs='Lgalaxies_test_750_MGS'
name_folder='Lgalaxies_test_750/'

#results: string of .sav dir and filename that contains the results of VESPA's run on original spectrum
#results: string of .sav dir and filename that contains the results of VESPA's run on hp spectrum
#inputsfh: optional input of file containing input sfh and z fraction that mocks were created from, string containing dir and filename

sampler_chain=cPickle.load(open('SavedData/'+name_folder+'sampler_chain_'+name_graphs+'.pkl','rb'))

log_liklihoods=cPickle.load(open('SavedData/'+name_folder+'sampler_logliklihoods_'+name_graphs+'.pkl','rb'))

index_max_liklihood=np.unravel_index(log_liklihoods.argmax(),log_liklihoods.shape)

best_fit_params=np.exp(sampler_chain[index_max_liklihood[0],index_max_liklihood[1],:])
if type_galaxy in ['fake']:
	inputsfh=np.loadtxt('SavedData/'+name_folder+'input_arr'+name_graphs+'.out')

bin_edges=np.concatenate((np.array([0]),10**(np.linspace(np.log10(0.02),np.log10(14),16))))
bin_edges_log=np.log10(bin_edges[1:])
av_bin_width=[]
for i in range(1,len(bin_edges_log)):
	av_bin_width.append((bin_edges_log[i]-bin_edges_log[i-1]))
first_bin_edge=np.zeros((1))
first_bin_edge[0]=bin_edges_log[0]-np.average(av_bin_width)

bin_centres=bin_edges_log-(np.average(av_bin_width)/2.)
mass_weighted_age=10.**(sum(bin_centres*best_fit_params))

print mass_weighted_age

bin_edges_log=np.concatenate([first_bin_edge,bin_edges_log]) #VESPA's bin edges

sffrecon=np.concatenate(([0.],best_fit_params))
sffrecon=sffrecon/np.average(av_bin_width)

if type_galaxy in ['fake']:
	inputsfh=np.concatenate(([0.],inputsfh))
	inputsfh=inputsfh/np.average(av_bin_width)

my_xticks=[]
for i in range(0,len(bin_edges)):
	my_xticks.append(str('%.2f' % bin_edges[i]))

plt.figure(figsize=(10,5))
#ax3=plt.subplot(211)
plt.ylabel('SFF/dex')
if type_galaxy in ['fake']:
	plt.step(bin_edges_log,inputsfh,'r',label='Input SFF/dex')
	for i in range(0,len(bin_edges_log)):
		plt.plot([bin_edges_log[i],bin_edges_log[i]],[0,inputsfh[i]],'r')
plt.step(bin_edges_log,sffrecon,'k',ls='--',label='Recovered SFF/dex')
for i in range(0,len(bin_edges_log)):
	plt.plot([bin_edges_log[i],bin_edges_log[i]],[0,sffrecon[i]],'k',ls='--')
plt.xticks(bin_edges_log,my_xticks)
plt.xlim(min(bin_edges_log),max(bin_edges_log))
plt.legend(loc='best')
#plt.ylim(0.,1.1)

'''

ax4=plt.subplot(212)
plt.ylabel('Z')
plt.step(bin_edges_log,zrecon,'r',ls='--')
for i in range(0,len(bin_edges_log)):
	plt.plot([bin_edges_log[i],bin_edges_log[i]],[0,hpzrecon[i]],'r',ls='--')
plt.step(bin_edges_log,hpzrecon,'b',ls='--')
for i in range(0,len(bin_edges_log)):
	plt.plot([bin_edges_log[i],bin_edges_log[i]],[0,hpzrecon[i]],'b',ls='--')
plt.xticks(bin_edges_log,my_xticks)
plt.xlim(min(bin_edges_log),max(bin_edges_log))
if max(zrecon)>max(hpzrecon):
	plt.ylim(0.,max(zrecon)+.1)
else:
	plt.ylim(0.,max(hpzrecon)+.1)
	
'''
plt.xlabel('Lookback Time (Gyr)')
plt.savefig('Images/'+name_folder+'SFH_pDex_Reconstruction'+name_graphs+'.pdf')
plt.show()