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
import ytree

if __name__ == "__main__":
    file="/media/shahram/SD/Sample100Mpc/717/treedata/trees_264.0"
    a=ytree.load(file)
    print(a[0])
