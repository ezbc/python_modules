import mergertreeHDF5

skipfac = 1
snapnum = 135
treebase = "/n/hernquistfs1/Illustris/Runs/Illustris-2/trees/treedata/"
mt = mergertreeHDF5.merger_tree(treebase, skipfac, snapnum)

for snapnum in range(0,snapnum+1):
	print snapnum

	x,y,  h = mt.getNumberOfMergers(snapnum)

	print h

