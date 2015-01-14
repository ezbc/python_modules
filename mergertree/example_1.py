import mergertreeHDF5

skipfac = 1
snapnum = 135
ntree = 0
nhalo = 0
property = "SnapNum"
treebase = "/n/hernquistfs1/Illustris/Runs/Illustris-2/trees/treedata/"

mt = mergertreeHDF5.merger_tree(treebase, skipfac, snapnum)


print "NtreesPerFile:"
print mt.NtreesPerFile
print "NumberOfOutputFiles:"
print mt.NumberOfOutputFiles
print "ParticleMass"
print mt.ParticleMass
print "TreeNHalos:"
print mt.TreeNHalos 
print "TotNsubhalos:"
print mt.TotNsubhalos 
print "Redshifts:"
print mt.Redshifts

AllProgenitors = mt.getAllProgenitors(ntree, nhalo)
print mt.trees[ntree][property][AllProgenitors]
	
FirstProgenitors = mt.getFirstProgenitors(ntree, nhalo)
print mt.trees[ntree][property][FirstProgenitors]
	
Progenitors = mt.getProgenitors(ntree, nhalo)
print mt.trees[ntree][property][Progenitors]

