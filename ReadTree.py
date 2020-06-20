#  © Shahram Talei @ 2020 The University of Alabama
#you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation; either version 3 of the License, or
#(at your option) any later version.
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.
# This code is distributed as is and there is no warranty or technical support
# This code reads the merger tree

import numpy as np
#import ytree
from struct import *

if __name__ == "__main__":
    file="/media/shahram/SD/Sample100Mpc/717/treedata/trees_264.0"
    MergerTreeFile=open(file,"rb")
    #fileContent=list(MergerTreeFile.readline())
    fileContent=MergerTreeFile.read()
    count=unpack("ii", fileContent[:8])
    print(count)
    Ntrees=count[0]
    totNhalos=count[1]
    print("There are %d trees and %d halos."%(Ntrees,totNhalos))
    TreeNhalos=unpack("i"*Ntrees, fileContent[8:4*Ntrees+8])
    print(TreeNhalos[0])
    structSize=104 #6*4+14*4+8+3*4+4
    structFormat="iiiiiiffffffffffffffqiiif"
    step=len(structFormat)
    print("step:%d"%step)
    TreeIndex=0
    offset=4*Ntrees+8
    start=0
    if TreeIndex==0:
        start=0
    else:
        for i in range(0,TreeIndex):
            start+=TreeNhalos[i]*structSize
    Position=start+offset
    L=structSize*TreeNhalos[TreeIndex]-Position
    print(L)
    #Tree=unpack(structFormat*TreeNhalos[TreeIndex], fileContent[4*Ntrees:structSi*TreeNhalos[TreeIndex]])
    Tree=unpack(structFormat*TreeNhalos[TreeIndex], fileContent[Position:structSize*TreeNhalos[TreeIndex]+Position])
    print(np.array(Tree).shape)
    #print(Tree[21])
    k=0
    Descendant=np.array([0]*TreeNhalos[TreeIndex])
    FirstProgenitor=np.array([0]*TreeNhalos[TreeIndex])
    NextProgenitor=np.array([0]*TreeNhalos[TreeIndex])
    FirstHaloInFOFgroup=np.array([0]*TreeNhalos[TreeIndex])
    NextHaloInFOFgroup=np.array([0]*TreeNhalos[TreeIndex])
    Len=np.array([0]*TreeNhalos[TreeIndex])
    Mvir=np.array([0.0]*TreeNhalos[TreeIndex])
    SnapNum=np.array([0]*TreeNhalos[TreeIndex])
    SubhaloIndex=np.array([0]*TreeNhalos[TreeIndex])
    SubHalofMass=np.array([0.0]*TreeNhalos[TreeIndex])
    #
    for i in range(0,len(Tree),step):
        Descendant[k]=Tree[i]
        FirstProgenitor[k]=Tree[i+1]
        NextProgenitor[k]=Tree[i+2]
        FirstHaloInFOFgroup[k]=Tree[i+3]
        NextHaloInFOFgroup[k]=Tree[i+4]
        Len[k]=Tree[i+5]
        Mvir[k]=Tree[i+7]*1.0e10
        SnapNum[k]=Tree[i+21]
        SubhaloIndex[k]=Tree[i+23]
        SubHalofMass[k]=Tree[i+24]
        k+=1
    #print((Mvir[0:100]))
    #print(set(SnapNum))
    print(SnapNum[0:280])
    #print(Descendant[0:100])
    MergerTreeFile.close()
    ##
    #Descendant=np.array(Descendant0)
    #SnapNum=np.array(SnapNum0)
    #just playing with tree
    ss=264
    print(SnapNum == ss)
    print("test the extracts:")
    S=SnapNum[SnapNum == ss]
    D=Descendant[SnapNum==ss]
    SubH=SubhaloIndex[SnapNum ==ss]
    FirstFOF=FirstHaloInFOFgroup[SnapNum ==ss]
    #print(len(D))
    #print(D)
    print("FirstProgenitor:")
    print(FirstProgenitor[SnapNum==ss])
    print("NextProgenitor:")
    print(NextProgenitor[SnapNum==ss])
    print("FirstHaloInFOF:")
    print(FirstHaloInFOFgroup[SnapNum==ss])
    print("NextHaloInFOFgroup:")
    print(NextHaloInFOFgroup[SnapNum==ss])
    print("SubhaloIndex:")
    print(SubhaloIndex[SnapNum==ss])
    print("Mvir:")
    print(Mvir[SnapNum==ss])
    ############################
    #now let's extract our branch of the merger tree
    index=TreeIndex
    Smax=np.max(SnapNum)
    Smin=np.min(SnapNum)
    print(Smax,Smin)
    LL=Smax-Smin
    Snap=[0]*LL
    HID=[0]*LL
    SubID=[0]*LL
    j=0
