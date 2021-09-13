from distutils import extension
from distutils.core import setup
from Cython.Build import cythonize
from numpy import *
import time
from operator import itemgetter
from multiprocessing import Pool





def GENIE3(expr_data, gene_names=None, regulators=all, tree_method="RF", K="sqrt", ntrees=1000,
        ncores=1, nthreads=1, seed=None):

    time_start = time.time()
	
    
    if not isinstance(expr_data,ndarray) or not len(expr_data.shape) == 2:
        raise ValueError('expr_data must be an array in which each row corresponds to a condition/sample and each column corresponds to a gene')

    #gives number of columns in array which must then correspond to gene names 
    num_genes = expr_data.shape[1]
	

	
	
    if gene_names is not None:
        if not isinstance(gene_names,(list,tuple)):
            raise ValueError('input argument gene_names must be a list of gene names')
        elif len(gene_names) != num_genes:
            raise ValueError('input argument gene_names must be a list of length p, where p is the number of columns/genes in the expr_data')	

    if regulators is not 'all':
        if not isinstance(regulators,(list,tuple)):
            raise ValueError('input argument regulators must be a list of gene names')

        if gene_names is None:
            raise ValueError('the gene names must be specified (in input argument gene_names)')
        else:
            sIntersection = set(gene_names).intersection(set(regulators))
            if not sIntersection:
                raise ValueError('the genes must contain at least one candidate regulator')



    if not isinstance(ncores,int):
            raise ValueError('input argument nthreads must be a stricly positive integer')
    elif ncores <= 0:
            raise ValueError('input argument nthreads must be a stricly positive integer')

    if not isinstance(nthreads,int):
            raise ValueError('input argument nthreads must be a stricly positive integer')
    elif nthreads <= 0:
            raise ValueError('input argument nthreads must be a stricly positive integer')



    if tree_method is not 'RF' and tree_method is not 'ET':
        raise ValueError('input argument tree_method must be "RF" (Random Forests) or "ET" (Extra-Trees)')

    
    if  K is not 'sqrt' and K is not 'all' and not isinstance(K,int): 
        raise ValueError('input argument K must be "sqrt", "all" or a stricly positive integer')
        
    if isinstance(K,int) and K <= 0:
        raise ValueError('input argument K must be "sqrt", "all" or a stricly positive integer')



    if not isinstance(ntrees,int):
        raise ValueError('input argument ntrees must be a stricly positive integer')
    elif ntrees <= 0:
        raise ValueError('input argument ntrees must be a stricly positive integer')
	

    if seed is not None:
        set.seed(seed)

    
 #confirming params
    print('Tree method: ' + str(tree_method))
    print('K: ' + str(K))
    print('Number of trees: ' + str(ntrees))
    print('\n')






	
 
	
    # Get the indices of the candidate regulators
    if regulators == 'all':
            input_idx = list(range(ngenes))
    else:
            input_idx = [i for i, gene in enumerate(gene_names) if gene in regulators]


    # Learn an ensemble of trees for each target gene, and compute scores for candidate regulators
    
    #parallel computing / multithreading
    
    VIM = zeros((ngenes,ngenes))
    
    if nthreads > 1:
        print('running jobs on %d threads' % nthreads)

        input_data = list()
        for i in range(ngenes):
            input_data.append( [expr_data,i,input_idx,tree_method,K,ntrees] )

        pool = Pool(nthreads)
        alloutput = pool.map(wr_GENIE3_single, input_data)
    
        for (i,vi) in alloutput:
            VIM[i,:] = vi

    else:
        print('running single threaded jobs')
        for i in range(ngenes):
            print('Gene %d/%d...' % (i+1,ngenes))
            
            vi = GENIE3_single(expr_data,i,input_idx,tree_method,K,ntrees)
            VIM[i,:] = vi

   
    VIM = transpose(VIM)
 
    time_end = time.time()
    print("Elapsed time: %.2f seconds" % (time_end - time_start))
    
    return VIM


    #read/write thread
def wr_GENIE3_single(args):
     return([args[1], GENIE3_single(args[0], args[1], args[2], args[3], args[4], args[5])])







cdef extension from "GENIE3.h":
void BuildTreeEns(int *nbobj, int *nbatt, CORETABLE_TYPE *X, CORETABLE_TYPE *Y, int *nm, int *et, int *rf, int *rfk, int *nbterms, int *bs, int *fsrf, SCORE_TYPE *vimp)

def py_BuildTreeEns(nbobj, nbatt, X, Y, nm, et, rf, rfk, nbterms, bs, fsrf, vimp):
    BuildTreeEns(nbobj, nbatt, X, Y, nm, et, rf, rfk, nbterms, bs, fsrf, vimp)




def GENIE3_single(expr_data,output_idx,input_idx,tree_method,K,ntrees):
    
    cdef int ngenes = expr_data.shape[1]
    
    # Expression of target gene
    output = expr_data[:,output_idx]
    
    # Normalize output data
    output = output / std(output)
    
    # Remove target gene from candidate regulators 
    input_idx = input_idx[:]
    if output_idx in input_idx:
        input_idx.remove(output_idx)

    expr_data_input = expr_data[:,input_idx]


    if tree_method == 'RF':
            randomForest = 1
            extraTrees = 0
            bootstrap_sampling = 1
    else: 
            randomForest = 0
            extraTrees = 1
            bootstrap_sampling = 0

    # set mtry
    if isinstance(K,int):
 	    mtry = K
    elif K == "sqrt": 
        mtry = round(sqrt(input_idx))
    else:
        mtry = input_idx
              
    cdef int numObj = expr_data.shape(0)

 	# some default parameters 
    nmin = 1
    bootstrap_sampling = 1
    permutation_importance = 0

    py_BuildTreeEns(numObj, ngenes, expr_data, output, nmin, extraTrees, randomForest, mtry, ntrees, bootstrap_sampling, permutation_importance, ("double",ngenes) )





def get_link_list(VIM,gene_names=None,regulators='all',maxcount='all',file_name=None):
    
    
    # Check input arguments      
    if not isinstance(VIM,ndarray):
        raise ValueError('VIM must be a square array')
    elif VIM.shape[0] != VIM.shape[1]:
        raise ValueError('VIM must be a square array')
        
    ngenes = VIM.shape[0]
        
    if gene_names is not None:
        if not isinstance(gene_names,(list,tuple)):
            raise ValueError('input argument gene_names must be a list of gene names')
        elif len(gene_names) != ngenes:
            raise ValueError('input argument gene_names must be a list of length p, where p is the number of columns/genes in the expression data')
        
    if regulators is not 'all':
        if not isinstance(regulators,(list,tuple)):
            raise ValueError('input argument regulators must be a list of gene names')

        if gene_names is None:
            raise ValueError('the gene names must be specified (in input argument gene_names)')
        else:
            sIntersection = set(gene_names).intersection(set(regulators))
            if not sIntersection:
                raise ValueError('The genes must contain at least one candidate regulator')
        
    if maxcount is not 'all' and not isinstance(maxcount,int):
        raise ValueError('input argument maxcount must be "all" or a positive integer')
        
    if file_name is not None and not isinstance(file_name,str):
        raise ValueError('input argument file_name must be a string')
    
    

    # Get the indices of the candidate regulators
    if regulators == 'all':
        input_idx = range(ngenes)
    else:
        input_idx = [i for i, gene in enumerate(gene_names) if gene in regulators]
    
    # Get the non-ranked list of regulatory links
    vInter = [(i,j,score) for (i,j),score in ndenumerate(VIM) if i in input_idx and i!=j]
    
    # Rank the list according to the weights of the edges 
    vInter_sort = sorted(vInter,key=itemgetter(2),reverse=True)
    nInter = len(vInter_sort)
    
    
    # Random permutation of edges with score equal to 0 
    flag = 1
    i = 0
    while flag and i < nInter:
        (TF_idx,target_idx,score) = vInter_sort[i]
        if score == 0:
            flag = 0
        else:
            i += 1
            
    if not flag:
        items_perm = vInter_sort[i:]
        items_perm = random.permutation(items_perm)
        vInter_sort[i:] = items_perm
       
    # Write the ranked list of edges


    #writing ranked list to file 

    nToWrite = nInter
    if isinstance(maxcount,int) and maxcount >= 0 and maxcount < nInter:
        nToWrite = maxcount
        
    if file_name:
    
        outfile = open(file_name,'w')
    
        if gene_names is not None:
            for i in range(nToWrite):
                (TF_idx,target_idx,score) = vInter_sort[i]
                TF_idx = int(TF_idx)
                target_idx = int(target_idx)
                outfile.write('%s\t%s\t%.6f\n' % (gene_names[TF_idx],gene_names[target_idx],score))
        else:
            for i in range(nToWrite):
                (TF_idx,target_idx,score) = vInter_sort[i]
                TF_idx = int(TF_idx)
                target_idx = int(target_idx)
                outfile.write('G%d\tG%d\t%.6f\n' % (TF_idx+1,target_idx+1,score))
            
        
        outfile.close()
        
    else:
        
        if gene_names is not None:
            for i in range(nToWrite):
                (TF_idx,target_idx,score) = vInter_sort[i]
                TF_idx = int(TF_idx)
                target_idx = int(target_idx)
                print('%s\t%s\t%.6f' % (gene_names[TF_idx],gene_names[target_idx],score))
        else:
            for i in range(nToWrite):
                (TF_idx,target_idx,score) = vInter_sort[i]
                TF_idx = int(TF_idx)
                target_idx = int(target_idx)
                print('G%d\tG%d\t%.6f' % (TF_idx+1,target_idx+1,score))
# cdefs not included, cprofiling TOdo, still runs but slower, maybe benchmark both impl



