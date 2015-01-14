import mergertreeHDF5

skipfac = 1
snapnum = 135
treebase = "/n/hernquistfs1/Illustris/Runs/Illustris-2/trees/treedata/"

combineLookup = False 

#just read merger tree header
mt = mergertreeHDF5.merger_tree(treebase, skipfac, snapnum, filenum = 0, tree_start = 0, tree_num = 0)


if (combineLookup == True):
	#combine partial lookup tables
	mt.combineSubhaloLookup(snapnum)

	#save global lookup table
	mt.saveSubhaloLookup(treebase, snapnum)
else:
	#read global look table
	mt.loadSubhaloLookup(treebase, snapnum)

#lookup halo number 23123
[filenum, ntree, nhalo] = mt.lookupSubhalo(23123)

#read trees
mt = mergertreeHDF5.merger_tree(treebase, skipfac, snapnum, filenum = filenum)

#print it 
print filenum, ntree, nhalo, mt.trees[ntree]["SubhaloNumber"][nhalo]

