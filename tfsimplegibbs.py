import tensorflow as tf
import numpy as np
import pdb
sweeps = 10000

p = 10
X0 = np.float32(np.random.randn(100,p))
w0 = np.float32(np.random.randn(p,1))
y0 = X0.dot(w0)

X = tf.Variable(X0)#np.atleast_2d(np.linspace(-1,1,100,dtype=np.float32)).T)#tf.placeholder("float")
y0 = tf.Variable(y0)#np.atleast_2d(3*np.linspace(-1,1,100,dtype=np.float32)).T)#tf.placeholder("float")
y = y0+tf.random_normal(y0.get_shape())
D = tf.Variable(np.eye(p,dtype=np.float32))#tf.placeholder("float")
alpha = tf.Variable(1.0)
beta = tf.Variable(np.zeros([p,1],dtype=np.float32))
sigvar = tf.Variable(1.0)

#feed_dict={D: np.eye(p), X: , y: 

Sinv = alpha*D+tf.matmul(tf.transpose(X),X)
S = tf.matrix_inverse(Sinv)
mu = tf.reshape(tf.matmul(S,tf.matmul(tf.transpose(X),y)),(p,1))

def gibbs(beta, mu, S):
	for param in xrange(p):
		vark = 1./S[param,param]
		span = tf.reshape(S[param,:],(1,p))
		muk = mu[param,:] - vark*(tf.reshape(tf.matmul(span,beta-mu),[1])-S[param,param]*(beta[param,:]-mu[param,:]))
		beta = tf.scatter_update(beta,param, tf.random_normal([1],mean=muk, stddev=tf.sqrt(vark)))
	return beta
sweepstore = []

gibbssweep = gibbs(beta, mu, S) 
init = tf.initialize_all_variables()

with tf.Session() as sess:
	#pdb.set_trace()
	X.initializer.run()
	y0.initializer.run()
	D.initializer.run()
	alpha.initializer.run()
	beta.initializer.run()
	sigvar.initializer.run()
	#init.run()
	for sweep in xrange(sweeps):
		sweepstore.append(sess.run(gibbssweep))
	print S.eval()
	print mu.eval()
sess.close()
allbeta = np.column_stack(sweepstore)

runningmean = np.cumsum(allbeta,axis=1)/(1+np.arange(allbeta.shape[1])).T[None,:]

import matplotlib.pyplot as plt
plt.plot(runningmean.T)
plt.show()

