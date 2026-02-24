# Functions are adopted from divDyn
# https://CRAN.R-project.org/package=divDyn
# Functions can still be optimzied

import numpy as np

def expandFAD(fadladTable):
    # Unpack taxon list into matrix
    divDynMatrix = np.full(shape = (int(np.sum([fadladTable[i, 2] - fadladTable[i, 1] for i in range(fadladTable.shape[0])])), 2), fill_value = 0)
    counts = 0
    for i in range(fadladTable.shape[0]):
        for k in np.arange(fadladTable[i,1], fadladTable[i,2]):
            divDynMatrix[counts,0] = i
            if k >= fadladTable[i,1] and k <= fadladTable[i,2]:
                divDynMatrix[counts,1] = k
            counts += 1
    return divDynMatrix

def counts(divDynMat):
    tax = divDynMat[:,0]
    nbin = divDynMat[:,1]

    nOccPoints = len(nbin)
    nBins = np.max(nbin) + 1
    nTax = len(np.unique(tax)) + 1
    
    TotalMatrix = np.zeros(shape = (nTax, nBins), dtype=bool)

    # Fill the matrix
    for i in range(nOccPoints):
        TotalMatrix[int(tax[i]), int(nbin[i])] = True
    
    # Initialize vectors
    t1 = np.zeros(nBins)
    t2u = np.zeros(nBins)
    t2d = np.zeros(nBins)
    t3 = np.zeros(nBins)
    tP = np.zeros(nBins)
    tGFu = np.zeros(nBins)
    tGFd = np.zeros(nBins)
    s1d = np.zeros(nBins)
    s2d = np.zeros(nBins)
    s3d = np.zeros(nBins)
    s1u = np.zeros(nBins)
    s2u = np.zeros(nBins)
    s3u = np.zeros(nBins)
    singleton = np.zeros(nBins)
    tThrough = np.zeros(nBins)
    tExtNoSing = np.zeros(nBins)
    tOriNoSing = np.zeros(nBins)
    divSIB = np.zeros(nBins)
    
    # Process taxa
    for j in range(nTax):
        FAD, LAD = None, None
        
        for i in range(nBins):
            if TotalMatrix[j, i]:
                if FAD is None:
                    FAD = i
                LAD = i
                divSIB[i] += 1
                
                if i > 0 and i < nBins - 1:
                    if not TotalMatrix[j, i - 1] and not TotalMatrix[j, i + 1]:
                        t1[i] += 1
                    if TotalMatrix[j, i - 1] and TotalMatrix[j, i + 1]:
                        t3[i] += 1
                    if TotalMatrix[j, i - 1] and not TotalMatrix[j, i] and TotalMatrix[j, i + 1]:
                        tP[i] += 1
                
                if i > 0 and TotalMatrix[j, i - 1]:
                    t2d[i] += 1
                if i < nBins - 1 and TotalMatrix[j, i + 1]:
                    t2u[i] += 1
                if i > 0 and i < nBins - 2 and TotalMatrix[j, i - 1] and not TotalMatrix[j, i + 1] and TotalMatrix[j, i + 2]:
                    tGFu[i] += 1
                if i > 1 and i < nBins - 1 and TotalMatrix[j, i - 2] and not TotalMatrix[j, i - 1] and TotalMatrix[j, i + 1]:
                    tGFd[i] += 1
                
                if i > 0 and i < nBins - 2:
                    if TotalMatrix[j, i - 1] and TotalMatrix[j, i] and not TotalMatrix[j, i + 1] and not TotalMatrix[j, i + 2]:
                        s1d[i] += 1
                    if TotalMatrix[j, i - 1] and not TotalMatrix[j, i] and TotalMatrix[j, i + 1] and not TotalMatrix[j, i + 2]:
                        s2d[i] += 1
                    if TotalMatrix[j, i - 1] and not TotalMatrix[j, i] and not TotalMatrix[j, i + 1] and TotalMatrix[j, i + 2]:
                        s3d[i] += 1
                
                if i > 1 and i < nBins - 1:
                    if not TotalMatrix[j, i - 2] and not TotalMatrix[j, i - 1] and TotalMatrix[j, i] and TotalMatrix[j, i + 1]:
                        s1u[i] += 1
                    if not TotalMatrix[j, i - 2] and TotalMatrix[j, i - 1] and not TotalMatrix[j, i] and TotalMatrix[j, i + 1]:
                        s2u[i] += 1
                    if TotalMatrix[j, i - 2] and not TotalMatrix[j, i - 1] and not TotalMatrix[j, i] and TotalMatrix[j, i + 1]:
                        s3u[i] += 1
        
        if FAD is not None and LAD is not None:
            if FAD == LAD:
                singleton[FAD] += 1
            else:
                tOriNoSing[FAD] += 1
                tExtNoSing[LAD] += 1
                tThrough[FAD + 1:LAD] += 1
    
    endMatrix = np.column_stack([
        t1, t2d, t2u, t3, tP, tGFd, tGFu, s1d, s2d, s3d, 
        s1u, s2u, s3u, singleton, tOriNoSing, tExtNoSing, tThrough, divSIB
    ])
    # t1:0
    # t2d:1
    # t2u:2
    # t3:3 
    # tP:4
    # tGFd:5
    # tGFu:6
    # s1d:7
    # s2d:8
    # s3d:9
    # s1u:10
    # s2u:11
    # s3u:12
    # singleton:13
    # tOriNoSing:14
    # tExtNoSing:15
    # tThrough:16
    # divSIB:17
    return endMatrix

# Two-for-three extinction
def getExt2f3(countsMat):
    # s1d, s2d, s3d
    sSubE = np.sort(countsMat[:,7:10], axis = 1)[:,1]
    # s1d, sSubE, t2d, tP
    ext2f3 = (countsMat[:,7] - sSubE) / (countsMat[:,1] + countsMat[:,4])
    ext2f3P = np.log(1 / (1 - ext2f3))
    return ext2f3P

# Two-for-three orgination
def getOrg2f3(countsMat):
    # s1u, s2u, s3u
    sSubU = np.sort(countsMat[:,10:13], axis = 1)[:,1]
    # s1u, sSubU, t2u, tP
    ori2f3 = (countsMat[:,10] - sSubU) / (countsMat[:,2] + countsMat[:,4])
    ori2f3P = np.log(1 / (1 - ori2f3))
    return ori2f3P

# Boundary-Crosser Diversity (BC)
def getBCdiversity(countsMat):
    # tThrough, tExt
    divBC = countsMat[:,16] + countsMat[:,15]
    return divBC

# Range Through Diversity (RT)
def getdivRT(countsMat):
    # tThrough, tExt, tOri, tSing
    divRT = (countsMat[:,16] + countsMat[:,15] + countsMat[:,14] + countsMat[:,13])
    return divRT

# Proportional Extinction
def getPropExt(countsMat):
    # tExt, tSing, divRT
    extP = (countsMat[:,15] + countsMat[:,13]) /  (countsMat[:,16] + countsMat[:,15] + countsMat[:,14] + countsMat[:,13])
    return extP

# Proportional Origination
def getPropOri(countsMat):
    # tOri, tSing, divRT
    orgP = (countsMat[:,14] + countsMat[:,13]) /  (countsMat[:,16] + countsMat[:,15] + countsMat[:,14] + countsMat[:,13])
    return orgP

# Per capita extinction Rate
def extPC(countsMat):
    # tExt, tSing, divRT
    extPC = -np.log((countsMat[:,16]) /  (countsMat[:,16] + countsMat[:,15]))
    return extPC

# Per capita origination Rate
def oriPC(countsMat):
    # tOri, tSing, divRT
    orgPC = -np.log((countsMat[:,16]) /  (countsMat[:,16] + countsMat[:,14]))
    return orgPC

def binnedToExpandedFad(input_binned_occurence):
    # Expand binned matrixes
    A = input_binned_occurence.copy()
    # Total number of unique species per time bin
    v_l = 0
    for i in range(len(A)):
        v_l += np.unique(A[i]).shape[0]

    # Get species list
    tmp = np.full(shape = (v_l, 2), fill_value=0, dtype=np.int64)
    index_t = 0
    for i in range(len(A)):
        tmp_A = np.unique(A[i])
        for j in range(len(tmp_A)):
            tmp[index_t, 0] = tmp_A[j]
            tmp[index_t, 1] = i
            index_t += 1
    # Assign new species ID
    _, new_vals = np.unique(tmp[:,0], return_inverse=True)
    tmp[:,0] = new_vals

    # Sort by species column then by bin column
    tmp2 = np.full_like(tmp, fill_value=0)
    species = np.unique(tmp[:,0])
    lt = 0
    for i in range(len(species)):
        # Get rows with species
        sp_tmp = tmp[tmp[:,0] == species[i], :]
        # Sort
        sp_tmp = sp_tmp[sp_tmp[:, 1].argsort()]
        new_species_counter = sp_tmp.shape[0]
        # Put into new array
        tmp2[lt:lt+new_species_counter,:] = sp_tmp
        lt += new_species_counter
    return tmp2