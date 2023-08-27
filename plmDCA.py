import numpy as np
from numpy.linalg import norm
import pickle as pkl
from scipy.spatial.distance import squareform, pdist
from scipy.optimize import minimize
import sys

# Global variable
q = 22

def plmDCA_main(fastafile,outputfile,reweighting_threshold):
    # maxcor: Defined in processInputOptions as 100
    # ftol: Not sure what it is equivalent to, assuming c1 from processInputOptions
    # gtol: Not sure what it is equivalent to, assuming c1 from processInputOptions
    # eps: Defined in main as 1e-9
    # maxFun: Defined in processInputOptions as 1000
    # maxIter: Defined in processInputOptions as 500
    
    options = {'disp':None,'maxcor':100,'ftol':1e-4,'gtol':1e-4,'eps':1e-9,'maxfun':1000,'maxiter':500,'iprint':-1,'maxls':100}     
    
    N,B_with_id_seq,q,Y = return_alignment(fastafile)
    Y = np.unique(Y,axis=0)
    B,N = Y.shape
    weights = np.ones(B)
    if reweighting_threshold > 0.0:
#        print('Starting to calculate weights \n...')
        weights = 1./(1+np.sum(squareform(pdist(Y,"hamming")),axis=0))
#        print('Finished calculating weights \n')
    B_eff = np.sum(weights)
    
    print('### N = {} B_with_id_seq = {} B = {} B_eff = {} q = {}\n'.format(N,B_with_id_seq,B,B_eff,q))
    
    if B_eff>500:
        lambda_J=0.01
    else:
        lambda_J=0.1-(0.1-0.01)*B_eff/500
        
    lambda_h=lambda_J
    scaled_lambda_h=lambda_h*B_eff
    scaled_lambda_J=lambda_J*B_eff/2
    
    w = np.zeros((q+q*q*(N-1),N))
    for r in range(N):
        print('Minimizing g_r for node r=' + str(r))
        wr = min_g_r(Y,weights,N,q,scaled_lambda_h,scaled_lambda_J,r,options)
        w[:,r] = wr
        
    #Extracting the coupling estimates from w.
    JJ=np.reshape(w[q:,:],(q,q,N-1,N))
    Jtemp1=np.zeros((q,q,int(N*(N-1)/2)))
    Jtemp2=np.zeros((q,q,int(N*(N-1)/2)))
    l=0
    for i in range(N-2):
        for j in range(i,N-1):
            Jtemp1[:,:,l] = JJ[:,:,j-1,i]
            Jtemp2[:,:,l] = JJ[:,:,i,j].transpose()
            l=l+1

    #Shift the coupling estimates into the Ising gauge.
    J1=np.zeros((q,q,int(N*(N-1)/2)))
    J2=np.zeros((q,q,int(N*(N-1)/2)))
    for l in range(int(N*(N-1)/2)):
        J1[:,:,l] = Jtemp1[:,:,l] - np.mean(Jtemp1[:,:,l],axis=0)[None,:]-np.mean(Jtemp2[:,:,l],axis=1)[:,None] + np.mean(Jtemp1[:,:,l])
        J2[:,:,l] = Jtemp2[:,:,l] - np.mean(Jtemp1[:,:,l],axis=0)[None,:]- np.mean(Jtemp2[:,:,l],axis=1)[:,None] + np.mean(Jtemp1[:,:,l])

    #Take J_ij as the average of the estimates from g_i and g_j.
    J=0.5*(J1+J2)
    
    #Calculate frob. norms FN_ij.
    NORMS=np.zeros((N,N))
    l=1
    for i in range(N-1):
        for j in range(i+1,N):
            NORMS[i,j] = norm(J[1:,1:,l])
            NORMS[j,i] = NORMS[i,j]

    #Calculate scores CN_ij=FN_ij-(FN_i-)(FN_-j)/(FN_--), where '-'
    #denotes average
    norm_means = np.mean(NORMS,axis=0)*N/(N-1)
    norm_means_all = np.mean(NORMS)*N/(N-1)
    CORRNORMS=NORMS-norm_means.transpose()*norm_means/norm_means_all

    
    data = {'J': J,
           'h': None,
           'frobenius_norm': NORMS,
           'corrected_norm': CORRNORMS}

    values = []
    for k in range(data['J'].shape[-1]):
        cur_value = 0
        for i in range(22):
            for j in range(22):
                cur_value += data['J'][i][j][k] * data['J'][i][j][k]
        cur_value = np.sqrt(cur_value)
        values.append(cur_value)
    
    save(values,outputfile)

    
    
    
weights=0
N=0

    
    

def save(arr,fileName):
    fileObject = open(fileName, 'wb')
    pkl.dump(arr, fileObject)
    fileObject.close()

def return_alignment(fastafile):
    with open(fastafile) as f:
        lines = f.read().split()
    B=len(lines)
    N=len(lines[0])
    Y=np.zeros((B,N),dtype=np.dtype(int))
    for i in range(B):
        counter = 0
        for j in range(N):
            Y[i,j] = letter2number(lines[i][j])
    q=len(aa_order)    
    return N,B,q,Y

aa_order = '-ACDEFGHIKLMNPQRSTVWYX'

def letter2number(a):
    if a in aa_order:
        return aa_order.find(a)
    else:
        return 0
    
def min_g_r(Y,weights,N,q,scaled_lambda_h,scaled_lambda_J,r,options):
    funObj= lambda wr: g_r(wr,Y,weights,q,scaled_lambda_h,scaled_lambda_J,r)
    funGrad=lambda wr: g_r_grad(wr,Y,weights,q,scaled_lambda_h,scaled_lambda_J,r)
    wr0 = np.zeros(q+q*q*(N-1))
    wr = minFunc(funObj,funGrad,wr0,options)
    return wr

# I'm calculating g_r and g_r_grad seperately.
# In the MATLAB implementation, they are calculated
# together since it is computationally more efficient.
# However, the scipy minimization doesn't seem readily
# compatible for this. We could either have the gradients
# stored somewhere or rewrite the L-BFGS algorithm.  I'm
# concerned that the L-BFGS-B algorithm gives slightly
# different results.
    
def g_r(wr,Y,weights,q,lambda_h,lambda_J,r):
    '''optimized version of g_r_slow'''
    (B,N) = Y.shape
    
    h_r = wr[:q]
    J_r = np.reshape(wr[q:],(q,q,N-1))
    
    indices = np.concatenate([np.arange(r),np.arange(r+1,N)])
    fval = 0
    
    for i in range(B):
        seq = Y[i,indices]
        logPot = h_r + np.sum(J_r[:,seq,np.arange(N-1)],axis=1)
        z = np.sum(np.exp(logPot))
        fval += weights[i]*(np.log(z) - logPot[Y[i,r]]) # 91-92
    
    # Add contributions from R_l2
    fval += lambda_h*np.sum(h_r**2) + lambda_J*np.sum(J_r**2)
    
    return fval

def g_r_grad(wr,Y,weights,q,lambda_h,lambda_J,r):
    '''optimized version of the gradient of g_r_slow'''
    (B,N) = Y.shape
    
    h_r = wr[:q]
    J_r = np.reshape(wr[q:],(q,q,N-1))
    
    indices = np.concatenate([np.arange(r),np.arange(r+1,N)])
    grad1 = np.zeros(q)
    grad2 = np.zeros((q,q,N-1))
            
    for i in range(B):
        seq = Y[i,indices]
        logPot = h_r + np.sum(J_r[:,seq,np.arange(N-1)],axis=1)
        z = np.sum(np.exp(logPot))
        
        nodeBel = np.exp(logPot - np.log(z))
        grad1[Y[i,r]] -= weights[i]
        grad1 += np.sum(weights[i]*nodeBel)
        
        for n in range(N-1):
            grad2[Y[i,r],Y[i,indices[n]],n] -= weights[i]
            grad2[:,Y[i,indices[n]],n] += weights[i]*nodeBel
            
    # Add contributions from R_l2
    grad1 += lambda_h*2*np.sum(h_r)
    grad2 += lambda_J*2*J_r
    grad = np.concatenate([grad1.flatten(),grad2.flatten()])
    
    return grad
            
'''
def g_r_slow(wr,Y,weights,N,q,lambda_h,lambda_J,r):
    q=22
    h_r = wr[:q]
    J_r = np.reshape(wr[q:],(q,q,N-1))
    (nInstances,nNodes) = Y.shape
    
    grad1 = np.zeros(q)
    grad2 = np.zeros((q,q,nNodes-1))
    logPot = np.zeros(q)
    nodeBel = np.zeros(q)
    
    fval = 0
    
    for i in range(nInstances):
        for s in range(q):
            logPot[s] = h_r[s]
        for n in range(nNodes):
            if (n!=r):
                y2 = int(Y[i,n])
                for s in range(q):
                    logPot[s] += J_r[s,y2,int(n-(1 if n>r else 0))]
                    
        z=0
        for s in range(q):
            z += np.exp(logPot[s])
        fval -= weights[i]*logPot[Y[i,r]]
        fval += weights[i]*np.log(z)
        
        # Gradient
        
        #for s in range(q):
        #    nodeBel[s] = np.exp(logPot[s] - np.log(z))
        #    
        #y1 = Y[i,r]
        #grad1[y1] -= weights[i]
        #
        #for s in range(q):
        #    grad1[s] += weights[i]*nodeBel[s]
        #
        #for n in range(nNodes):
        #    if (n!=r):
        #        y2 = Y[i,n]
        #        grad2[y1,y2,n-(1 if n>r else 0)] -= weights[i]
        #        for s in range(q):
        #            grad2[s,y2,n-(1 if n>r else 0)] += weights[i]*nodeBel[s]
                    
    # Add contributions from R_l2
    for s in range(q):
        fval += lambda_h*h_r[s]*h_r[s]
        #grad1[s] += lambda_h*2*h_r[s]
        
    for n in range(nNodes):
        if (n!=r):
            for s in range(q):
                for t in range(q):
                    fval += lambda_J*J_r[s,t,n-(1 if n>r else 0)]*J_r[s,t,n-(1 if n>r else 0)]
                    #grad2[s,t,n-(1 if n>r else 0)] += lambda_J*2*J_r[s,t,n-(1 if n>r else 0)]
                    
    return fval #, grad1, grad2'''

# I'm working on making it work with lbfgs
def minFunc(funObj,funGrad,x0,options):
    res = minimize(funObj,x0, jac=funGrad,method='L-BFGS-B', options=options)
    return res.x

# plmDCA_main('ex1.txt','ex1.pkl',0.1)

#if __name__ == "__main__":
#    main(sys.argv[1],sys.argv[2],float(sys.argv[3]))
