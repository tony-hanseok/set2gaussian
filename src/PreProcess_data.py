import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")
import numpy as np
import collections
import operator
import sys
import random
import gc
import os
import matplotlib
import cPickle as pickle
matplotlib.use('agg')
import matplotlib.pyplot as plt
from Set2Gaussian.BioNetwork import BioNetwork
from Set2Gaussian.RandomWalkRestart import RandomWalkRestart, DCA_vector
from Set2Gaussian.GenerateDiffusion import GenerateDiffusion
from Set2Gaussian.utils import evaluate_pathway_member
from Set2Gaussian.model import Graph2Gauss
from sklearn.model_selection import train_test_split
import scipy.spatial as sp


def read_node_embedding(DCA_dim,network_file, DCA_rst=0.8):
	net_file_l = []
	net_file_l.append(network_file)
	Net_obj = BioNetwork(net_file_l)
	network = Net_obj.sparse_network.toarray()
	print ('#node:',np.shape(network)[0])
	i2g = Net_obj.i2g
	g2i = Net_obj.g2i
	nnode = len(i2g)
	Node_RWR = RandomWalkRestart(network, DCA_rst)
	node_emb, _,_,_,node_context = DCA_vector(Node_RWR,DCA_dim)
	print ('preprocess node embedding is finished')
	return Net_obj, Node_RWR, node_emb, node_context

def create_matrix(Node_RWR, Net_obj, p_train, gene_set_file):
	GR_obj = GenerateDiffusion(Node_RWR, gene_set_file, Net_obj=Net_obj)
	p2i = GR_obj.p2i
	_, Path_mat_train_all, Path_mat_test_all, _, _ = GR_obj.RunDiffusion(p_train=p_train, random_state=0,all_gene_cv=False)

	npath,nnode = np.shape(Path_mat_train_all)
	print ('#node set', npath)
	nsmooth = max(npath,nnode)
	alpha = 1./(nsmooth*nsmooth)
	node_alpha = 1./(nnode*nnode)
	log_node_RWR =  np.log(Node_RWR +node_alpha) - np.log(node_alpha)
	Path_RWR = np.dot(Path_mat_train_all, Node_RWR)
	log_Path_RWR = -1 * np.log(Path_RWR +alpha)
	train_ind, test_ind = train_test_split(range(nnode), test_size=0.01)
	train_ind = np.array(range(nnode))
	return Path_RWR, log_Path_RWR, log_node_RWR, train_ind, test_ind, Path_mat_train_all, Path_mat_test_all, p2i

def save_embedding(p2g, path_mu, path_cov, Grep_node_emb, p2i, output_file):
	np.save(output_file + 'p2g.out', p2g)
	np.save(output_file + 'path_mu.out', path_mu)
	np.save(output_file + 'path_cov.out', path_cov)
	np.save(output_file + 'g2g_node_emb.out', Grep_node_emb)
	with open(output_file + 'g2i', "wb") as f:
		pickle.dump(p2i, f)

def run_embedding_method(method,log_Path_RWR, log_node_RWR, Path_RWR, node_emb, node_context,train_ind,test_ind, Path_mat_train, para_dict, metric= 'cosine'):
	npath, nnode = np.shape(Path_RWR)
	if method == 'Set2Gaussian':
		Grep_obj = Graph2Gauss(log_Path_RWR, log_node_RWR, Path_RWR, node_emb, node_context,
		path_batch_size = 20, node_batch_size = 5000,lr = para_dict['lr'],
		L=para_dict['node_emb_dim'],optimize_diag_path=para_dict['optimize_diag_path'],optimize_path_mean = para_dict['optimize_path_mean'],n_hidden = [para_dict['nhidden']],early_stopping=para_dict['early_stopping'],gene_loss_lambda=para_dict['gene_loss_lambda'],max_iter=para_dict['max_iter'],seed=0,train_ind=train_ind,test_ind = test_ind)#change 200 to 20
		path_mu, path_cov, Grep_node_emb, p2g = Grep_obj.train()
		return path_mu, path_cov, Grep_node_emb, p2g
	if method == 'Network_smoothed_mean':
		Path_emb = np.dot( Path_RWR, node_emb)
		Path_avg_emb = sp.distance.cdist(Path_emb, node_emb, metric)
		return [],[],[],Path_avg_emb
	if method == 'Sum':
		Path_emb = np.dot(Path_mat_train, node_emb)
		Path_avg_emb = sp.distance.cdist(Path_emb,node_emb, metric)
		return [],[],[],Path_avg_emb
	if method == 'Mean':
		Path_emb = np.dot(Path_mat_train, node_emb)
		for i in range(npath):
			Path_emb[i,:] /= np.sum(Path_mat_train[i,:])
		Path_avg_emb = sp.distance.cdist(Path_emb,node_emb, metric)
		return [],[],[],Path_avg_emb
	if method == 'Max':
		for i in range(npath):
			path_gene = np.where(Path_mat_train[i,:]>0)[0]
			Path_emb[i,:] = np.max(node_emb[path_gene,:],axis=0)
		Path_avg_emb = sp.distance.cdist(Path_emb, node_emb, metric)
		return [],[],[],Path_avg_emb

def evaluate_embedding(flog_file, p2g, Path_mat_test, Path_mat_train):
	auroc_d, auroc_l,prec_d,prec_l = evaluate_pathway_member(p2g, Path_mat_test, Path_mat_train, low_b=[3,11,31],up_b=[10,30,1000])
	flog = open(flog_file,'w')
	for d in auroc_d:
		aup,l = auroc_d[d]
		prec,l = prec_d[d]
		flog.write(str(aup)+'\t'+str(prec)+'\n')
		for rs in auroc_l[d]:
			flog.write('part'+'\t'+str(d)+'\t'+str(rs)+'\n')
	auroc_d, auroc_l,prec_d,prec_l = evaluate_pathway_member(p2g, Path_mat_test, Path_mat_train, low_b=[3],up_b=[1000])
	for d in auroc_d:
		aup,l = auroc_d[d]
		prec,l = prec_d[d]
		flog.write(str(aup)+'\t'+str(prec)+'\n')
		for rs in auroc_l[d]:
			flog.write('all'+'\t'+str(d)+'\t'+str(rs)+'\n')
	flog.write('\n')
	flog.close()
