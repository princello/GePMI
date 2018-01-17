import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import beta
import sys
import os
from sklearn import metrics
from sklearn.metrics import average_precision_score
from sklearn.metrics import precision_recall_curve
import statsmodels.sandbox.stats.multicomp as sssm

def welcome():
	print "\n============================"
	print "Thank you for using GePMI :D\n============================\n"
	print "usage: GePMI -i input.csv -p 0.001 -q 0.01 -s 0 -o outputDir -t\n"
	print "-i your input sourmash csv file\nPlease name it as follows:\n[prefix]-[base number(millions)]-[k-mer length]-[hashes used in sourmash(thousands)].csv\nfor eaxamle:\ntest-100-18-10.csv\n"
	print "-p threshold of p-value (default 0.001)\n"
	print "-q threshold of q-value (default  0.01)\n"
	print "-s threshold of similarity (default  0)\n"
	print "-o path of output dirctory (default ./output)\n"
	print "-f save temple p/q value files\n"
	print "-t save temple p/q value figures\n"

args = sys.argv
check_point = 0
save_label_files = 0
save_label_figures = 0
input_csv=''
if len(args) == 1:
	print "\nWelcome! -h for help\n"
elif '-h' in args:
	welcome()
else:
	if '-i' in args:
		try:
			input_csv = args[args.index('-i')+1]
			if input_csv[-4:] != '.csv':
				print "ERROR: please check the input file!"
				welcome()
			else:
				check_point = 1
		except:
			print "ERROR: Please use the standard input!"
			welcome()
	if '-p' in args:
		try:
			threshold_p = float(args[args.index('-p')+1])
		except:
			print "ERROR: Please use the correct p-value threshold value!"
			check_point = 0
			welcome()
	else:
		threshold_p = 0.001
	if '-q' in args:
		try:
			threshold_q = float(args[args.index('-q')+1])
		except:
			print "ERROR: Please use the correct q-value threshold value!"
			check_point = 0
			welcome()
	else:
		threshold_q = 0.01
	if '-s' in args:
		try:
			threshold_s = float(args[args.index('-s')+1])
		except:
			print "ERROR: Please use the correct similarity threshold value!"
			check_point = 0
			welcome()
	else:
		threshold_s = 0.01
	if '-o' in args:
		try:
			path = args[args.index('-o')+1]
			if os.path.exists(path):
				while True:
					print 'ERROR: output dirctory exists! continue? [y/n]\n'
					y = raw_input()
					if y == 'y':
						check_point = 1
						break
					elif y == 'n':
						check_point = 0
						break
					else:
						print "ERROR: output dirctory exists! continue [y/n]?\n"

			else:
				os.system('mkdir '+path)
		except:
			print "ERROR: Please check the output dirctory!"
			check_point = 0
	else:
		path='output'
		x = os.system('mkdir '+path)
		if x == 256:
			while True:
				print "ERROR: output dirctory exists! continue [y/n]?\n"
				y = raw_input()
				if y == 'y':
					check_point = 1
					break
				elif y == 'n':
					check_point = 0
					break
				else:
					print "ERROR: output dirctory exists! continue [y/n]?\n"

	if '-f' in args:
		save_label_files = 1
	if '-t' in args:
		save_label_figures = 1

def rename(sample_name,min_s):
	original = sample_name
	try:
		if sample_name.index('.') < sample_name.index('-'):
			sample_name = sample_name.replace('.','',1)
			if sample_name.index('.')<sample_name.index('-'):
				rename(sample_name)
			else:
				#print original+' replaced with '+sample_name
				return sample_name.split('.')[0]
		else:
			return sample_name.split('.')[0]
	except:
		if sample_name.count('.')<=min_s:
			return sample_name.split('.')[0]
		else:
			sample_name = sample_name.replace('.','',sample_name.count('.')-min_s)
			return sample_name.split('.')[0]

def read(sn):
	f = open(sn,'r')
	name = f.readline().split(',')
	name[-1] = name[-1][:-2]
	suffix_length = []
	for i in range(len(name)):
		name[i] =name[i].split('/')[-1]
		suffix_length.append(name[i].count('.'))
	min_s = min(suffix_length)
	for i in range(len(name)):
		name[i] = rename(name[i],min_s)
	data = []
	while True:
		line = f.readline()
		if not line:
			break
		else:
			l = line.split(',')
			l[-1] = l[-1][:-1]
			data.append(l)
	return name,np.array(data).astype(dtype='float')

def inwhat(name,data):
	inner = []
	intra = []
	for i in range(len(data)):
		for j in range(len(data)):
			if name[i].split('-')[0] == name[j].split('-')[0] and name[i]!=name[j]:
				inner.append(data[i][j])
			elif name[i]!=name[j]:
				intra.append(data[i][j])
	return inner,intra

def pvalue(name,data):
	p_value = []
	total_number = len(name)
	for i in range(len(name)):
		sys.stdout.write('Sample tested: %d/'%(i)+str(total_number)+' \r')
		sys.stdout.flush()
		ppvp = []
		d1 = np.where(data[i]>0,data[i],0.0000000001)
		d1 = np.where(d1<1,d1,0.9999999999)
		rl = []
		for j in range(len(d1)):
			if name[j].split('-')[0] == name[i].split('-')[0]:
				rl.append(j)
		d2 = []
		for k in range(len(d1)):
			if k in rl:
				continue
			else:
				d2.append(d1[k])
		try:
			param = beta.fit(d2,floc=0, fscale=1)
		except:
			print d2
		rv = beta(param[0],param[1],0,1)
		for j in range(len(d1)):
			if j!=i:
				ppvp.append(1-rv.cdf(d1[j]))
			else:
				ppvp.append(1)
		p_value.append(ppvp)
	if save_label_files == 1:
		np.savetxt(path+'/p_values.txt',np.array(p_value))
	return np.array(p_value)

def fig(name,data,p_value,q_value,sm,kl,sk):
	sns.set(style='white')
	plt.clf()
	w = open(path+'/detail.txt','w')
	w.write('Target\tTest\tMinHash\tp-value\tq-value\ttag(known individual samples=1)\n')
	d_true_score = []
	p_true_score = []
	q_true_score = []
	d_false_score = []
	p_false_score = []
	q_false_score = []
	for i in range(len(data)):
		sys.stdout.write("sample: %d \r" %(i))
		sys.stdout.flush()
		d1 = data[i]
		p1 = p_value[i]
		q1 = q_value[i]
		rl = []
		for j in range(len(d1)):
			if name[j].split('-')[0] == name[i].split('-')[0]:
				rl.append(j)
		sl = 'abc'
		for j in range(len(d1)):
			if j in rl and j!=i and name[j].split('_')[0][:3]!=sl:
				d_true_score.append(d1[j])
				p_true_score.append(1-p1[j])
				q_true_score.append(1-q1[j])
				w.write(name[i]+'\t'+name[j]+'\t'+str(d1[j])+'\t'+str(p1[j])+'\t'+str(q1[j])+'\t'+str(1)+'\n')
			elif j!=i and name[j].split('_')[0][:3]!=sl:
				d_false_score.append(d1[j])
				p_false_score.append(1-p1[j])
				q_false_score.append(1-q1[j])
				w.write(name[i]+'\t'+name[j]+'\t'+str(d1[j])+'\t'+str(p1[j])+'\t'+str(q1[j])+'\t'+str(0)+'\n')
	w.close()
	result = []
	f = open(path+'/detail.txt','r')
	w = open(path+'/result.txt','w')
	w.write('Target sample\tsignificant similar Test sample\n')
	line = f.readline()
	sample_now = ''
	tmp = []
	while True:
		line = f.readline()
		if not line:
			break
		else:
			ll = line.split('\t')
			if ll[0] == sample_now:
				if float(ll[3]) <= threshold_p and float(ll[4]) <= threshold_q and float(ll[2]) >= threshold_s:
					tmp.append(ll[1])
			else:
				if sample_now!='':
					w.write(sample_now+'\t')
					for test in tmp:
						w.write(test+' ')
					w.write('\n')
				tmp = []
				if float(ll[3]) <= threshold_p and float(ll[4]) <= threshold_q and float(ll[2]) >= threshold_s:
					tmp.append(ll[1])
			sample_now = ll[0]
	w.close()
	if save_label_figures == 1:
		true_len = len(d_true_score)
		false_len = len(d_false_score)
		d_scores = []
		d_scores.extend(d_true_score)
		d_scores.extend(d_false_score)
		d_scores = np.array(d_scores)

		d_y = []
		d_y.extend([1 for n in range(true_len)])
		d_y.extend([0 for m in range(false_len)])
		d_y = np.array(d_y)
		fpr, tpr, thresholds = metrics.roc_curve(d_y, d_scores, pos_label=1)
		roc_auc = metrics.auc(fpr, tpr)
		d, = plt.plot(fpr,tpr,color='#E91E63',lw=7, label='MinHash (auROC = %0.4f)' % roc_auc)

		p_scores = []
		p_scores.extend(p_true_score)
		p_scores.extend(p_false_score)
		p_scores = np.array(p_scores)

		p_y = []
		p_y.extend([1 for n in range(true_len)])
		p_y.extend([0 for m in range(false_len)])
		p_y = np.array(p_y)
		fpr, tpr, thresholds = metrics.roc_curve(p_y, p_scores, pos_label=1)
		roc_auc = metrics.auc(fpr, tpr)
		p, = plt.plot(fpr,tpr,color='#9C27B0',lw=7,label='p_value (auROC = %0.4f)' % roc_auc)

		q_scores = []
		q_scores.extend(q_true_score)
		q_scores.extend(q_false_score)
		q_scores = np.array(q_scores)

		q_y = []
		q_y.extend([1 for n in range(true_len)])
		q_y.extend([0 for m in range(false_len)])
		q_y = np.array(q_y)
		fpr, tpr, thresholds = metrics.roc_curve(q_y, q_scores, pos_label=1)
		roc_auc = metrics.auc(fpr, tpr)
		q, = plt.plot(fpr,tpr,color='#3F51B5',lw=7,label='q_value (auROC = %0.4f)' % roc_auc)

		plt.xlim([0, 1.05])
		plt.ylim([0.0, 1.05])
		plt.xlabel('False Positive Rate',fontsize=30)
		plt.ylabel('True Positive Rate',fontsize=30)

		#plt.title('ROC curve',fontsize=36)
		#plt.plot([],[],'*',color='red', label='k-mer number: '+sm+',000,000,000' )

		spac, = plt.plot([],[],'*',color='white', label=' ')
		kmer, = plt.plot([],[],'*',color='red', label='k-mer size: '+kl)
		sket, = plt.plot([],[],'*',color='red', label='number of hashes: '+sk )
		plt.legend(loc="lower right",fontsize=18)
		axisx = plt.gca().xaxis
		for label in axisx.get_ticklabels():
			label.set_fontsize(24)
			label.set_rotation(0)

		axisy = plt.gca().yaxis
		for label in axisy.get_ticklabels():
			label.set_fontsize(24)
			label.set_rotation(0)
		plt.tight_layout()
		#plt.show()
		plt.savefig(path+'/roc.'+sm+'-'+kl+'-'+sk+'.pdf',dpi=300)

		plt.clf()

		precision, recall, thresholds = precision_recall_curve(d_y, d_scores,pos_label=1)
		average_precision=average_precision_score(d_y,d_scores)
		plt.plot(recall, precision, lw=7, color='#E91E63',label='MinHash (auPRC={0:0.4f})'.format(average_precision))

		precision, recall, thresholds = precision_recall_curve(p_y, p_scores,pos_label=1)
		average_precision=average_precision_score(p_y,p_scores)
		plt.plot(recall, precision, '#9C27B0',lw=7, label='p_value (auPRC={0:0.4f})'.format(average_precision))

		precision, recall, thresholds = precision_recall_curve(q_y, q_scores,pos_label=1)
		average_precision=average_precision_score(q_y,q_scores)
		plt.plot(recall, precision, '#3F51B5',lw=7, label='q_value (auPRC={0:0.4f})'.format(average_precision))


		spac, = plt.plot([],[],'*',color='white', label=' ')
		kmer, = plt.plot([],[],'*',color='red', label='k-mer size: '+kl)
		sket, = plt.plot([],[],'*',color='red', label='number of hashes: '+sk)

		plt.xlabel('Recall',fontsize=30)
		plt.ylabel('Precision',fontsize=30)
		plt.ylim([0.0, 1.05])
		plt.xlim([0.0, 1.05])
		#plt.title('Precision-Recall curve',fontsize=36)
		plt.legend(loc="lower left",fontsize=18)
		axisx = plt.gca().xaxis
		for label in axisx.get_ticklabels():
			label.set_fontsize(24)
			label.set_rotation(0)

		axisy = plt.gca().yaxis
		for label in axisy.get_ticklabels():
			label.set_fontsize(24)
			label.set_rotation(0)
		plt.tight_layout()
		#plt.show()
		plt.savefig(path+'/prc.'+sm+'-'+kl+'-'+sk+'.pdf',dpi=300)

		plt.clf()
		precision_p = []
		precision_q = []
		q = np.arange(0.00001,1,0.01)
	 	for i in q:
			p1p = np.sum(1-np.array(p_true_score)<=i)
			p1q = np.sum(1-np.array(q_true_score)<=i)
			pcp = p1p + np.sum(1-np.array(p_false_score)<=i)
			pcq = p1q + np.sum(1-np.array(q_false_score)<=i)
			precision_p.append(1-p1p/float(pcp))
			precision_q.append(1-p1q/float(pcq))
		#plt.plot([0, 1], [0, 1], color='#673AB7', lw=3, linestyle='--')
		plt.plot(q, precision_p, color='#9C27B0',lw=6,label='p-value')
		plt.plot(q, precision_q, color='#3F51B5',lw=6,label='q-value')
		plt.xlabel('values cutoff',fontsize=30)
		plt.ylabel('FDR',fontsize=30)
		plt.ylim([0.0, 1.05])
		plt.xlim([0.0, 1.01])
		#plt.title('FDR control',fontsize=36)
		plt.legend(loc="center right",fontsize=18)
		axisx = plt.gca().xaxis
		for label in axisx.get_ticklabels():
			label.set_fontsize(24)
			label.set_rotation(0)

		axisy = plt.gca().yaxis
		for label in axisy.get_ticklabels():
			label.set_fontsize(24)
			label.set_rotation(0)
		plt.tight_layout()
		plt.savefig(path+'/fdr.'+sm+'-'+kl+'-'+sk+'.pdf',dpi=300)

		plt.clf()
		precision_p = []
		precision_q = []
		q = np.arange(0.00001,0.1,0.001)
	 	for i in q:
			p1p = np.sum(1-np.array(p_true_score)<=i)
			p1q = np.sum(1-np.array(q_true_score)<=i)
			pcp = p1p + np.sum(1-np.array(p_false_score)<=i)
			pcq = p1q + np.sum(1-np.array(q_false_score)<=i)
			precision_p.append(1-p1p/float(pcp))
			precision_q.append(1-p1q/float(pcq))
		plt.plot([0, 0.1], [0, 0.1], color='#673AB7', lw=3, linestyle='--')
		plt.plot(q, precision_p, color='#9C27B0',lw=6,label='p-value')
		plt.plot(q, precision_q, color='#3F51B5',lw=6,label='q-value')
		plt.xlabel('values cutoff',fontsize=30)
		plt.ylabel('FDR',fontsize=30)
		plt.ylim([0.0, 0.105])
		plt.xlim([0.0, 0.101])
		#plt.title('FDR control',fontsize=36)
		plt.legend(loc="center right",fontsize=18)
		axisx = plt.gca().xaxis
		for label in axisx.get_ticklabels():
			label.set_fontsize(24)
			label.set_rotation(0)

		axisy = plt.gca().yaxis
		for label in axisy.get_ticklabels():
			label.set_fontsize(24)
			label.set_rotation(0)
		plt.tight_layout()
		plt.savefig(path+'/fdr.'+sm+'-'+kl+'-'+sk+'.d.pdf',dpi=300)


def allin(sn,check_point):
	if check_point==1:
		sl = sn.split('/')[-1].split('.')[0].split('-')
		kl = sl[2]
		sm = sl[1]
		sk = sl[3]
		name,data = read(sn)
		p_value = pvalue(name,data)
		q_value = []
		for i in range(len(p_value[0])):
			qv = sssm.multipletests(p_value[:,i],method='fdr_by')[1]
			q_value.append(qv)
		q_value = np.transpose(np.array(q_value))
		if save_label_files == 1:
			np.savetxt(path+'/q_values.txt',q_value)
		try:
			fig(name,data,p_value,q_value,sm,kl,sk)
		except:
			welcome()
	else:
		print 'please check!'

allin(input_csv,check_point)




