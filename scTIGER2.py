#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Nishi Gupta
07.22.2024
Definitions file for run_scTIGER2.0.py

Updates list:
03.13.2026 - Madison Dautle | Added a flag for the causal validation step in TCDF
03.13.2026 - Madison Dautle | Updated pseudotime methods 
"""
import numpy as np
import pandas as pd
import os, sys
import scanpy as sc
from scipy.stats import pearsonr
import networkx as nx
import time
from datetime import timedelta
import bambi as bmb
import arviz as az
from scipy import stats 
import math
import shutil
from sklearn.linear_model import LinearRegression

def checkGeneList(geneList, normCD, allOutput_df):
    for gene in geneList:
        if gene not in normCD.columns:
            print(gene +" not in matrix")
            geneList.remove(gene)
            allOutput_df = allOutput_df.drop(gene, axis=1)
    return geneList, allOutput_df


def normalize(df, zeroThresh):
    df.replace(0,np.nan, inplace=True)
    df=df.dropna(thresh=len(df)*zeroThresh, axis=1)
    norm_df = ((df+1) - df.min())/(df.max()-df.min())
    norm_df.replace(np.nan, 0, inplace = True)
    return norm_df




def make_adata(norm_df):
    adata = sc.AnnData(norm_df.copy(), dtype=np.float64)
    adata.X = np.nan_to_num(adata.X, nan=0.0, posinf=0.0, neginf=0.0)

    # basic preprocessing shared by several methods
    sc.tl.pca(adata)
    sc.pp.neighbors(adata)
    return adata

def rank_pseudotime(time_df, col):
    time_df = time_df.copy()
    time_df["Rank"] = time_df[col].rank(method="first")
    time_df = time_df.sort_values(by="Rank")
    return list(time_df.index)

def reorder_cells(df, cell_order_list):
    df = df.T
    df = df[cell_order_list]
    df = df.T
    return df
    
def pseudotime_paga(norm_df):
    adata = sc.AnnData(norm_df.copy(), dtype=np.float64)
    adata.X = np.nan_to_num(adata.X, nan=0.0, posinf=0.0, neginf=0.0)

    sc.tl.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.leiden(adata, resolution=1.0)
    sc.tl.paga(adata, groups='leiden')

    df_res = adata.obs[['leiden']].copy()
    df_res['pseudotime'] = df_res['leiden'].rank(method='first')
    return df_res[['pseudotime']]

def pseudotime_dpt(norm_df, root_cell=None):
    adata = make_adata(norm_df)

    # DPT uses diffusion maps
    sc.tl.diffmap(adata)

    # pick a root cell index
    if root_cell is None:
        root_cell = adata.obs_names[0]

    # Scanpy wants integer index in adata.uns['iroot']
    adata.uns["iroot"] = np.where(adata.obs_names == root_cell)[0][0]

    sc.tl.dpt(adata)

    df_res = adata.obs[["dpt_pseudotime"]].copy()
    df_res.columns = ["pseudotime"]
    return df_res
      
def pseudotime_palantir(norm_df, root_cell=None):
    import palantir

    adata = make_adata(norm_df)

    # Palantir commonly works from PCA / diffusion representation
    pca_projections = pd.DataFrame(
        adata.obsm["X_pca"],
        index=adata.obs_names
    )

    dm_res = palantir.utils.run_diffusion_maps(pca_projections)
    ms_data = palantir.utils.determine_multiscale_space(dm_res)

    if root_cell is None:
        root_cell = adata.obs_names[0]

    pr_res = palantir.core.run_palantir(ms_data, early_cell=root_cell)

    df_res = pd.DataFrame(index=adata.obs_names)
    df_res["pseudotime"] = pr_res.pseudotime
    return df_res
     
def pseudotime_scfates(norm_df, root_cell=None):
    import numpy as np
    import pandas as pd
    import scanpy as sc
    import scFates as scf

    adata = sc.AnnData(norm_df.copy(), dtype=np.float64)
    adata.X = np.nan_to_num(adata.X, nan=0.0, posinf=0.0, neginf=0.0)

    sc.tl.pca(adata)
    sc.pp.neighbors(adata)

    # build tree on PCA space
    scf.tl.tree(adata, method="ppt", Nodes=50, use_rep="X_pca")

    # REQUIRED: define a root before pseudotime
    # simplest option: use a cell index
    if root_cell is None:
        root_cell = adata.obs_names[0]

    root_idx = np.where(adata.obs_names == root_cell)[0][0]
    scf.tl.root(adata, root=root_idx)

    # now pseudotime can be computed
    scf.tl.pseudotime(adata)

    # scFates writes pseudotime to adata.obs["t"]
    df_res = adata.obs[["t"]].copy()
    df_res.columns = ["pseudotime"]
    return df_res
    
def pseudotime_slingshot(norm_df, root_cell=None, root_cluster=None):
    import numpy as np
    import pandas as pd
    import scanpy as sc
    from pyslingshot import Slingshot

    adata = sc.AnnData(norm_df.copy(), dtype=np.float64)
    adata.X = np.nan_to_num(adata.X, nan=0.0, posinf=0.0, neginf=0.0)

    sc.tl.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.leiden(adata, resolution=1.0)

    # reduced coordinates
    coords = adata.obsm["X_pca"][:, :2]

    # original cluster labels as strings
    clusters = adata.obs["leiden"].astype(str)

    # choose root cluster label
    if root_cluster is None:
        if root_cell is not None and root_cell in adata.obs_names:
            root_cluster = clusters.loc[root_cell]
        else:
            root_cluster = clusters.iloc[0]

    # one-hot encode clusters
    cluster_onehot_df = pd.get_dummies(clusters)
    cluster_onehot = cluster_onehot_df.to_numpy()

    # map root cluster label -> integer column index
    cluster_names = list(cluster_onehot_df.columns)
    if root_cluster not in cluster_names:
        raise ValueError(
            f"Root cluster {root_cluster} not found in cluster names {cluster_names}"
        )
    start_node = cluster_names.index(root_cluster)

    sling = Slingshot(
        coords,
        cluster_onehot,
        start_node=start_node
    )

    sling.fit()

    df_res = pd.DataFrame(index=adata.obs_names)

    if hasattr(sling, "pseudotime"):
        pt = sling.pseudotime
    elif hasattr(sling, "unified_pseudotime"):
        pt = sling.unified_pseudotime
    else:
        print("Available Slingshot attributes:", dir(sling))
        raise ValueError("Could not find pseudotime output in Slingshot object.")

    df_res["pseudotime"] = np.asarray(pt).ravel()
    return df_res
    
def pseudotime_monocle3(norm_df, root_cell=None):
    import numpy as np
    import pandas as pd
    import scanpy as sc
    import py_monocle

    adata = sc.AnnData(norm_df.copy(), dtype=np.float64)
    adata.X = np.nan_to_num(adata.X, nan=0.0, posinf=0.0, neginf=0.0)

    sc.tl.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=1.0)

    matrix = np.asarray(adata.obsm["X_umap"])

    if root_cell is None:
        root_idx = 0
    else:
        matches = np.where(adata.obs_names == root_cell)[0]
        if len(matches) == 0:
            raise ValueError(f"Root cell {root_cell} not found in adata.obs_names")
        root_idx = int(matches[0])

    clusters = adata.obs["leiden"].astype("category").cat.codes.to_numpy()

    pt = py_monocle.pseudotime(
        matrix,
        root_cells=root_idx,
        clusters=clusters
    )

    df_res = pd.DataFrame(index=adata.obs_names)
    df_res["pseudotime"] = np.asarray(pt).ravel()
    return df_res
  
def get_pseudotime(norm_df, method="PAGA", root_cell=None, root_cluster=None):
    method = method.lower()

    if method == "paga":
        return pseudotime_paga(norm_df)
    elif method == "dpt":
        return pseudotime_dpt(norm_df, root_cell=root_cell)
    elif method == "palantir":
        return pseudotime_palantir(norm_df, root_cell=root_cell)
    elif method == "scfates":
        return pseudotime_scfates(norm_df, root_cell=root_cell)
    elif method == "slingshot":
        return pseudotime_slingshot(norm_df, root_cell=root_cell, root_cluster=root_cluster)
    elif method == "monocle3":
        return pseudotime_monocle3(norm_df, root_cell=root_cell)
    else:
        raise ValueError(f"Unsupported method: {method}")




def pearson_corr(normCD):
    corrM=normCD.corr(method='pearson')  
    corrM=corrM.rename_axis('temp')
    result = corrM.stack().reset_index() 
    result.columns = ['fromGene','toGene','corr'] 
    return result


def predict_GRN(geneList, corrResult, numTopRanked, normCD, cudaUse, timeDelay, noValidation):
    tcdfInput = "./TCDF_Output/inputTemp.csv"
    tcdfOutput= "./TCDF_Output/outputTemp.csv"
    allInteractions = pd.DataFrame(columns = geneList)
    for j in range(len(geneList)):
        interactionsPP = pd.DataFrame(columns = ('gene.pair', 'timesteps'))
        df = corrResult[corrResult['toGene'] == geneList[j]] 
        absList = pd.DataFrame({'AbsVal': df['corr'].abs()})
        df = pd.concat([df, absList], axis=1) 
        topGenes = df.sort_values(['AbsVal'], ascending=[False])
        gm = topGenes['fromGene'].head(numTopRanked + 1)
        gl = gm.tolist()
        geneMatrix = normCD.filter(gl)
        geneMatrix = geneMatrix.replace(np.nan, 0)
        geneMatrix.to_csv(tcdfInput, index=False)
        os.environ['MKL_THREADING_LAYER'] = 'GNU'
        if cudaUse == True and noValidation == False:
            command = 'python ./utils/runTCDF.py --data ' + tcdfInput + ' --cuda > ./TCDF_Output/tcdfRunning.txt'
        elif cudaUse == True and noValidation == True:
            command = 'python ./utils/runTCDF.py --data ' + tcdfInput + ' --cuda --no_validation > ./TCDF_Output/tcdfRunning.txt'
        elif cudaUse == False and noValidation == True: 
            command = 'python ./utils/runTCDF.py --data ' + tcdfInput + ' --no_validation > ./TCDF_Output/tcdfRunning.txt'
        else:
             command = 'python ./utils/runTCDF.py --data ' + tcdfInput + ' > ./TCDF_Output/tcdfRunning.txt'
        os.system(command)
        TCDFRes = pd.read_csv(tcdfOutput)
        interactionsPP = pd.concat([interactionsPP, TCDFRes], ignore_index=True)
        if len(interactionsPP) > 0:
            interactionsPP[['fromGene','toGene']] = interactionsPP.iloc[:,0].str.split('>',expand=True)
            interactionsPP.insert(len(interactionsPP.T), "Direction", np.nan)
            for r in range(len(interactionsPP)):
                fromGene = normCD.loc[:, interactionsPP.iloc[r,2]]
                toGene = normCD.loc[:, interactionsPP.iloc[r,3]]
                fromGene.replace(np.inf, 0, inplace=True)
                toGene.replace(np.inf, 0, inplace=True)
                corr, _ = pearsonr(fromGene, toGene) 
                if corr > 0:
                    interactionsPP.iloc[r, 4] = "+"
                else: 
                    interactionsPP.iloc[r, 4] = "-" 
            interactionsPP["gene.pair"] = interactionsPP["gene.pair"].astype(str) +","+ interactionsPP["Direction"]
            interactionsPP.drop(columns=['fromGene', 'toGene', 'Direction'], inplace = True)
        if timeDelay ==0:
            interactionsPP=interactionsPP[interactionsPP['timesteps'] ==0]
            interactionsPP=interactionsPP.drop('timesteps', axis = 1)
        else:
            interactionsPP=interactionsPP[interactionsPP['timesteps'] <= timeDelay]
            interactionsPP['gene.pair'] = interactionsPP['gene.pair'] + ',' + interactionsPP['timesteps'].astype(str)
            interactionsPP=interactionsPP.drop('timesteps', axis = 1)

        t = 0
        while t < len(interactionsPP):
            item = interactionsPP.iloc[t, interactionsPP.columns.get_loc('gene.pair')]
            new_df = pd.DataFrame({geneList[j]: [item]})
            allInteractions = pd.concat([allInteractions, new_df], ignore_index=True)
            t+=1
        print("Gene:\t" + geneList[j])
    return allInteractions


def combine_permutations(allInteractions, perms, geneList):
    for m in range(len(geneList)):
        counts = allInteractions[geneList[m]].value_counts()  
        outputFileName = geneList[m] + '.csv' 
        counts.to_csv(outputFileName,index=True)  
        
        os.chdir("./Graphs")
        from matplotlib import pyplot as plt
        df = counts.to_frame()
        plt.hist(np.array(df.iloc[:,0]))
        plt.yscale('log')
        plt.title(df.columns[0])
        outputGraphName = geneList[m] + '.png'
        plt.savefig(outputGraphName, dpi=300, bbox_inches='tight') 
        plt.clf()
        
        os.chdir('..')
        
def calSigLvl(bg_df, gene, alpha_val): #Significance level calculated using negative binomial
    bgAppCts = bg_df[gene].value_counts().sort_index() 
    bgAppCts = bgAppCts.to_frame() 
    bgAppCts['Freq'] = bgAppCts.index 
  
    try: 
        eq1 = gene+" ~ Freq" 
        model = bmb.Model(eq1, bgAppCts, family = "negativebinomial") 
        negB = model.fit()
        summary = az.summary(negB)
        z_val = stats.norm.ppf(1-alpha_val)    
        sig_lvl = summary.loc['Intercept', 'mean'] + ((z_val*summary.loc['Intercept', 'sd'])/math.sqrt(len(bg_df)))
    except:
        print("Negative binomial fitting failed. Taking average to calculate significance value.")
        sig_lvl = bg_df.loc[:, gene].mean() 
    return sig_lvl

def calSigLvl_avg(bg_df, gene, alpha_val): #Significance level calculated using mean of frequencies
    bgAppCts = bg_df[gene].value_counts().sort_index()
    bgAppCts = bgAppCts.to_frame()
    bgAppCts['Freq'] = bgAppCts.index
    sig_lvl = bg_df.loc[:, gene].mean() 
    return sig_lvl

def removeBackground(actual_df, sig_val, gene): 
    actual_noNoise = actual_df[actual_df[gene] >= sig_val]
    return actual_noNoise




def graphGRN(actual_df, gene, timeDelay): 
    actual_df.columns = ['genePair', 'weight']
    try:
        if timeDelay == 0:
            actual_df[['genePair','arrowShape']] = actual_df.iloc[:,0].str.split(',',expand=True)
        else:
            actual_df[['genePair','arrowShape','timeDelay']] = actual_df.iloc[:,0].str.split(',',expand=True)
        actual_df[['source','target']] = actual_df.iloc[:,0].str.split('>',expand=True)
        if timeDelay ==0:
            actual_df = actual_df[['source', 'target', 'weight', 'arrowShape']]
        else:
            actual_df = actual_df[['source', 'target', 'weight', 'arrowShape', 'timeDelay']]
        actual_df.replace("+", "delta", inplace=True)
        actual_df.replace("-", "T", inplace=True)
    
        G = nx.from_pandas_edgelist(actual_df, edge_attr=True)
        fileout = str(gene) + '.graphml'
        nx.write_graphml(G, fileout)
            
            
    except:
        if timeDelay ==0:
            actual_df = pd.DataFrame(columns = ['source', 'target', 'weight', 'arrowShape'])
        else:
            actual_df = pd.DataFrame(columns = ['source', 'target', 'weight', 'arrowShape', 'timeDelay'])
    return actual_df 





def scTIGER(outDir, permutations, cellNumber, rawCase, rawCtrl, zeroThresh, geneList, allInteractions, topGenes, cudaUse, timeDelay, noValidation, pseudotime):
    os.mkdir(outDir)
    for i in range(permutations):
        start = time.time()
        case=rawCase.sample(n=cellNumber,axis='columns')  
        ctrl=rawCtrl.sample(n=cellNumber, axis='columns') 
        case=case.T
        ctrl=ctrl.T
        
        
        # 1. normalization + pseudotime
        case = normalize(case, zeroThresh)
        ctrl = normalize(ctrl, zeroThresh)
        time_case = get_pseudotime(case, method=pseudotime, root_cell=case.index[0])
        time_ctrl = get_pseudotime(ctrl, method=pseudotime, root_cell=ctrl.index[0]) 
        case_cellOrder = rank_pseudotime(time_case, "pseudotime")
        ctrl_cellOrder = rank_pseudotime(time_ctrl, "pseudotime")
        case = reorder_cells(case, case_cellOrder)
        ctrl = reorder_cells(ctrl, ctrl_cellOrder)


        # 2. difference matrix
        case.reset_index(drop=True, inplace=True) 
        ctrl.reset_index(drop=True, inplace=True)
        caseMctrl = case.subtract(ctrl)
        if len(caseMctrl.columns) >= 10000: 
            print("You have " + str(len(caseMctrl.columns)) + " genes after normalization and filtering. You may want to increase the zero threshold for faster performance.")
        else: 
            print(str(len(caseMctrl.columns)) + " genes after normalization and filtering")
        geneList, allInteractions = checkGeneList(geneList, caseMctrl, allInteractions)


        # 3. correlation 
        corr = pearson_corr(caseMctrl) 
        
        
        # 4. Detect gene interactions
        interactions = predict_GRN(geneList, corr, topGenes, caseMctrl, cudaUse, timeDelay, noValidation)
        allInteractions = pd.concat([allInteractions, interactions], ignore_index=True)
        del(interactions) 

        exec_time = (time.time()-start)
        print("Permutation: \t" + str(i+1))
        print("Run time:\t" + str(timedelta(seconds=exec_time)))


    #combine permutations  
    os.chdir(outDir)
    os.mkdir("Graphs")
    combine_permutations(allInteractions, permutations, geneList)
    
    hold = 1
    return hold




def scTIGER2(outDir, permutations, cellNumber, rawCase, zeroThresh, geneList, allInteractions, topGenes, cudaUse, timeDelay, noValidation, pseudotime):
     os.mkdir(outDir)
     for i in range(permutations):
         start = time.time()
         case=rawCase.sample(n=cellNumber,axis='columns')
         case=case.T
         
         #1. normalization + pseudotime
         case = normalize(case, zeroThresh)
         time_case = get_pseudotime(case, method=pseudotime, root_cell=case.index[0])
         case_cellOrder = rank_pseudotime(time_case, "pseudotime")
         case = reorder_cells(case, case_cellOrder)
          
         #2. correlation
         case.reset_index(drop=True, inplace=True) 
         if len(case.columns) >= 10000:
             print("You have " + str(len(case.columns)) + " genes after normalization and filtering. You may want to increase the zero threshold for faster performance.")
         else:
             print(str(len(case.columns)) + " genes after normalization and filtering.")
         geneList, allInteractions = checkGeneList(geneList, case, allInteractions)
         
         corr = pearson_corr(case) 
         
         #3. detect gene interactions
         interactions = predict_GRN(geneList, corr, topGenes, case, cudaUse, timeDelay, noValidation)
         allInteractions = pd.concat([allInteractions, interactions], ignore_index=True)
         del(interactions)
         
         exec_time = (time.time()-start)
         print("Permutation: \t" + str(i+1))
         print("Run time:\t" + str(timedelta(seconds=exec_time)))
         
     #combine permutations
     os.chdir(outDir)
     os.mkdir("Graphs")
     combine_permutations(allInteractions, permutations, geneList)
     
     hold = 1
     return hold
 
 
 
    
def GRAPH(outDir, geneList, permutations, timeDelay, alpha_val, hold):
     if hold == 1: 
         try:
             os.mkdir('./GRN_Visualization')
             os.mkdir('./sig_gene_networks')
         except:
             ask_user = input("Do you want to overwrite the folder? y/n")
             if ask_user == 'y':
                 shutil.rmtree('./GRN_Visualization')
                 shutil.rmtree('./sig_gene_networks')
                 os.mkdir('./GRN_Visualization')
                 os.mkdir('./sig_gene_networks')
         os.chdir('./GRN_Visualization')
     else:
         os.chdir(outDir) 
         try:
             os.mkdir('./GRN_Visualization')
             os.mkdir('./sig_gene_networks')
         except:
             ask_user = input("Do you want to overwrite the folder? y/n")
             if ask_user == 'y':
                 shutil.rmtree('./GRN_Visualization')
                 shutil.rmtree('./sig_gene_networks')
                 os.mkdir('./GRN_Visualization')
                 os.mkdir('./sig_gene_networks')
         os.chdir('./GRN_Visualization') 
         
     if timeDelay ==0: 
         completeInteractions = pd.DataFrame(columns = ['source', 'target', 'weight', 'arrowShape'])
     else:
         completeInteractions = pd.DataFrame(columns = ['source', 'target', 'weight', 'arrowShape', 'timeDelay'])

     for g in geneList: 
         actual = pd.read_csv('../' + g + '.csv')  
         actual.rename(columns={g: 'interactions', 'count':g}, inplace=True) 
         sig_ct = calSigLvl(actual, g, alpha_val) 
         sig_ct_avg = calSigLvl_avg(actual, g, alpha_val)
             
         actual = removeBackground(actual, sig_ct, g) 
         os.chdir('../sig_gene_networks')
         filen = './' + g + '_sig.csv'
         actual.to_csv(filen, index=False)
         os.chdir('../GRN_Visualization')
         actual = graphGRN(actual, g, timeDelay)
         completeInteractions = pd.concat([completeInteractions, actual], ignore_index=True)
         
     G = nx.from_pandas_edgelist(completeInteractions, edge_attr=True)
     fileout = 'OverallInteractions.graphml'
     nx.write_graphml(G, fileout)
     
     return sig_ct, sig_ct_avg
           
