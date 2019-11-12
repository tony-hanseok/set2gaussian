import numpy as np
import tensorflow as tf
import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
from .utils import *
#from tensorflow_probability import distributions as tfd
import sys
from scipy.optimize import minimize
from scipy import stats
import GPUtil
import time
import psutil


class Graph2Gauss:
	def __init__(self, log_Path_RWR, log_node_RWR, Path_RWR, Path_mat_train, Path_mat_test, node_emb, node_context,
	 L=50, path_batch_size = 100,node_batch_size = 100,evalute_every_iter=False,early_stopping=50,
	 optimize_gene_vec=False,optimize_path_mean=True,optimize_diag_path=0,lr = 1e-2, K=1, n_hidden=[50],use_piecewise_loss=False,gene_loss_lambda=1.,
				 max_iter=200,  seed=0, verbose=True,eval_obj={}, train_ind=[],
			 	test_ind = [], flog='',EVALUATE=True):
		"""
		Parameters
		----------
		A : scipy.sparse.spmatrix
			Sparse unweighted adjacency matrix
		X : scipy.sparse.spmatrix
			Sparse attribute matirx
		L : int
			Dimensionality of the node embeddings
		K : int
			Maximum distance to consider
		p_val : float
			Percent of edges in the validation set, 0 <= p_val < 1
		p_test : float
			Percent of edges in the test set, 0 <= p_test < 1
		p_nodes : float
			Percent of nodes to hide (inductive learning), 0 <= p_nodes < 1
		n_hidden : list(int)
			A list specifying the size of each hidden layer, default n_hidden=[512]
		max_iter :  int
			Maximum number of epoch for which to run gradient descent
		tolerance : int
			Used for early stopping. Number of epoch to wait for the score to improve on the validation set
		seed : int
			Random seed used to split the edges into train-val-test set
		verbose : bool
			Verbosity.
		"""

		tf.reset_default_graph()
		tf.set_random_seed(seed)
		np.random.seed(seed)
		self.evalute_every_iter = evalute_every_iter
		self.EVALUATE = EVALUATE
		self.lr = lr
		self.flog = flog
		if self.flog == '':
			self.flog = open('test','w')
		self.optimize_gene_vec = optimize_gene_vec
		self.early_stopping = early_stopping
		self.optimize_path_mean = optimize_path_mean
		self.optimize_diag_path = optimize_diag_path
		self.gene_reg = gene_loss_lambda

		self.path_batch_size = path_batch_size
		self.node_batch_size = node_batch_size
		self.log_Path_RWR = log_Path_RWR.astype(np.float32)
		self.log_node_RWR = log_node_RWR.astype(np.float32)
		self.Path_RWR = Path_RWR.astype(np.float32)
		self.Path_mat_train = Path_mat_train.astype(np.float32)
		self.Path_mat_test = Path_mat_test.astype(np.float32)
		self.use_piecewise_loss = use_piecewise_loss
		self.npath, self.nnode = self.log_Path_RWR.shape
		#self.log_Path_RWR = tf.convert_to_tensor(self.log_Path_RWR)
		self.log_node_RWR = tf.convert_to_tensor(self.log_node_RWR)
		node_emb = node_emb.astype(np.float32)
		self.node_emb = node_emb
		node_context = node_context.astype(np.float32)
		self.D = node_emb.shape[1]
		self.L = L
		self.max_iter = max_iter
		self.verbose = verbose
		self.eval_obj = eval_obj

		self.train_ind = train_ind
		self.test_ind = test_ind

		self.ntest = len(self.test_ind)
		self.ntrain = len(self.train_ind)

		if n_hidden is None:
			n_hidden = [512]
		self.n_hidden = n_hidden

		self.__build( node_emb, node_context)
		sys.stdout.flush()
		#self.__dataset_generator(hops, scale_terms)
		self.__build_loss()
		sys.stdout.flush()
		# setup the validation set for easy evaluation

	def __build(self, node_mu, node_context):#0:diag path cov, 1: full path cov, 2:one side optimize, 3:three way dot
		w_init = tf.contrib.layers.xavier_initializer
		self.node_mu = tf.get_variable(name='node_mu', initializer=node_mu)
		self.node_context = tf.get_variable(name='node_context', initializer= node_context)
		self.path_cov_w = tf.get_variable(name='Path_sigma_w', dtype=tf.float32, shape=[self.npath, self.L,  self.n_hidden[0]],initializer= w_init())
		self.path_cov_x = tf.get_variable(name='Path_sigma_x', dtype=tf.float32, shape=[self.npath, self.L , self.n_hidden[0]],initializer= w_init())
		Path_emb = np.dot( self.Path_RWR,node_mu)
		self.path_mu = tf.get_variable(name='Path_mu',initializer= Path_emb)


	def __build_loss(self):
		self.path_slice_st = tf.placeholder(tf.int32)
		self.path_slice_ed = tf.placeholder(tf.int32)
		self.gene_slice_st = tf.placeholder(tf.int32)
		self.gene_slice_ed = tf.placeholder(tf.int32)
		self.node_ind_old = tf.placeholder(shape=[None], dtype=tf.int32)
		self.node_ind = tf.convert_to_tensor(self.node_ind_old)
		node_ind_L = tf.shape(self.node_ind)[0]
		#node_ind = self.train_ind[self.gene_slice_st:self.gene_slice_ed]
		self.Y_target = tf.placeholder(shape=[None, None], dtype=tf.float32)

		path_cov_w_batch = self.path_cov_w[self.path_slice_st:self.path_slice_ed,:,:]
		path_cov_x_batch = self.path_cov_x[self.path_slice_st:self.path_slice_ed,:,:]
		path_cov_x_trans = tf.transpose(path_cov_x_batch, perm=[0,2,1])
		C_k = tf.matmul(path_cov_w_batch, path_cov_x_trans)
		C_k_trans = tf.transpose(C_k, perm=[0,2,1])
		self.cov_inv_k =  tf.matmul(C_k, C_k_trans)

		path_mu_batch = self.path_mu[self.path_slice_st:self.path_slice_ed,:]
		path_mu_batch = tf.expand_dims(path_mu_batch,1)
		path_bar_batch_all_node = tf.tile(path_mu_batch,[1,node_ind_L,1])
		#node_mu = self.node_mu[self.node_ind,:]
		node_mu = tf.gather(self.node_mu,self.node_ind)

		node_mu = tf.expand_dims(node_mu,0)
		node_bar_batch_all_node = tf.tile(node_mu,[self.path_slice_ed-self.path_slice_st,1,1])
		pgd = path_bar_batch_all_node - node_bar_batch_all_node
		pdd = self.cov_inv_k
		pgd_pdd = tf.matmul(pgd, pdd)
		self.p2g = tf.reduce_sum( tf.multiply( pgd, pgd_pdd ), 2, keep_dims=False )


		if self.use_piecewise_loss:
			self.pathway_loss = tf.losses.mean_pairwise_squared_error(self.p2g , self.Y_target)
		else:
			self.pathway_loss = tf.losses.mean_squared_error(self.p2g, self.Y_target)

		node_mu = tf.gather(self.node_mu,self.node_ind)
		g2g_diff = tf.matmul(node_mu, tf.transpose(self.node_context))
		log_node_RWR_batch = tf.gather(self.log_node_RWR, self.node_ind)
		self.gene_loss = self.gene_reg * tf.losses.mean_squared_error(g2g_diff, log_node_RWR_batch)

		self.loss = self.pathway_loss + self.gene_loss
		#self.loss = self.gene_loss

	def get_pathway_loss(self, sess, node_ind):
		pathway_loss = 0.
		nnode = len(node_ind)
		for pi in range(0, self.npath, self.path_batch_size):
			pst = pi
			ped = pi + self.path_batch_size
			ped = min(ped,self.npath)
			for bi in range(0, nnode, self.node_batch_size):
				gst = bi
				ged = bi + self.node_batch_size
				ged = min(ged, nnode)
				pathway_loss+= sess.run(self.pathway_loss,feed_dict={self.path_slice_st:pst,self.path_slice_ed:ped,self.node_ind:node_ind[gst:ged],
				self.Y_target:self.log_Path_RWR[pst:ped, node_ind[gst:ged]]})
		return pathway_loss



	def train(self, gpu_list='0'):
		"""
		Trains the model.

		Parameters
		----------
		gpu_list : string
			A list of available GPU devices.

		Returns
		-------
		sess : tf.Session
			Tensorflow session that can be used to obtain the trained embeddings

		"""
		final_path_cov = {}
		final_g2g_node_emb = {}
		final_p2g = {}
		final_path_mu = {}

		early_stopping_score_max = -1.0
		cost_val = []
		train_op = tf.train.AdamOptimizer(learning_rate=self.lr).minimize(self.loss)

		sess = tf.Session(config=tf.ConfigProto(gpu_options=tf.GPUOptions(allow_growth=True)))
		sess.run(tf.global_variables_initializer())

		for epoch in range(self.max_iter):
			start_time = time.time()

			test_pathway_loss = self.get_pathway_loss( sess, self.test_ind)

			test_gene_loss = sess.run(self.gene_loss,feed_dict={self.path_slice_st:0,self.path_slice_ed:2,self.node_ind:self.test_ind,self.Y_target:self.log_Path_RWR[0:2,self.test_ind]})
			test_loss_all = test_gene_loss + test_pathway_loss
			cost_val.append(test_loss_all)

			train_pathway_loss = self.get_pathway_loss( sess, self.train_ind)
			train_gene_loss = sess.run(self.gene_loss,feed_dict={self.path_slice_st:0,self.path_slice_ed:2,self.node_ind:self.train_ind,self.Y_target:self.log_Path_RWR[0:2,self.train_ind]})
			train_loss_all = train_gene_loss + train_pathway_loss

			loss = 0.
			for pi in range(0,self.npath,self.path_batch_size):
				pst, ped = pi, min(pi + self.path_batch_size,self.npath)
				for bi in range(0,self.ntrain,self.node_batch_size):
					gst, ged = bi, min(bi + self.node_batch_size,self.ntrain)
					tmp_loss, _ = sess.run([self.loss, train_op],feed_dict={self.path_slice_st:pst,self.path_slice_ed:ped,self.node_ind:self.train_ind[gst:ged],self.Y_target:self.log_Path_RWR[pst:ped,self.train_ind[gst:ged]]})
					loss += tmp_loss
			#print epoch,train_loss_all,test_loss_all,test_gene_loss,test_pathway_loss

			if epoch > self.early_stopping and cost_val[-1] >  np.mean(cost_val[-(self.early_stopping+1):-1]):
				print 'Early stopping ..',cost_val[-1],np.mean(cost_val[-(self.early_stopping+1):-1]),len(cost_val[-(self.early_stopping+1):-1])
				break
			sys.stdout.flush()
			if epoch >= self.early_stopping-1 or epoch >= self.max_iter-1:
				final_p2g = np.zeros((self.npath, self.nnode))
				for pi in range(0,self.npath,self.path_batch_size):
					pst = pi
					ped = pi + self.path_batch_size
					ped = min(ped,self.npath)
					final_p2g[pst:ped,:] = sess.run(self.p2g,feed_dict={self.path_slice_st:pst,self.path_slice_ed:ped,self.node_ind:np.array(range(self.nnode)),self.Y_target:self.log_Path_RWR[pst:ped,:]})#npath*L

		tf.reset_default_graph()
		sess.close()
		return final_path_mu, final_path_cov, final_g2g_node_emb, final_p2g
