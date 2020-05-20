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
    structSize=96 #6*4+12*4+8+3*4+4
    structFormat="iiiiiiffffffffffffqiiif"
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
    MergerTreeFile.close()
