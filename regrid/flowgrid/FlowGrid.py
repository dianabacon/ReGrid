import os, io
import numpy as np
from datetime import *
import getpass
import string
import copy
import math

from mpl_toolkits.mplot3d import axes3d
from matplotlib.ticker import MaxNLocator
import matplotlib.pyplot as plt

from pyevtk.hl import pointsToVTK
import vtk

# import pkg_resources  # part of setuptools
# version = pkg_resources.require("ReGrid")[0].version

f2m = 0.3048  # ft to m


class FlowGrid(object):
    def __init__(self):
        self.skip = 0
        self.Prop = {}
        self.outputDir = 'output'

    def __getitem__(self, key):
        return getattr(self, key)

    def setOutputDir(self, dir):
        self.outputDir = dir

    # Checks if path outputDir/dirName exists, creates if not
    def checkOutputDir(self, subDir):
        if not os.path.exists(self.outputDir):
            os.mkdir(self.outputDir)
        if not os.path.exists(self.outputDir + '/' + subDir):
            os.mkdir(self.outputDir + '/' + subDir)

    def exportVTK(self, fname):
        """ Saves the SUTRA grid as a VTK file, either a VTKStructuredGrid (.vts)
            or a VTKUnstructuredGrid (.vtu) depending on mesh type.
            fname = the filename it will be saved at, if no extension is given,
            .vts is appended
        """
        filename, ext = os.path.splitext(fname)
        if self.GridType == "vtkStructuredGrid":
            sWrite = vtk.vtkXMLStructuredGridWriter()
            sWrite.SetInputData(self.Grid)
            sWrite.SetFileName(filename + ".vts")
            sWrite.Write()
        elif self.GridType == "vtkUnstructuredGrid":
            sWrite = vtk.vtkXMLUnstructuredGridWriter()
            sWrite.SetInputData(self.Grid)
            sWrite.SetFileName(filename + ".vtu")
            sWrite.Write()
        else:
            print("Grid type is not recognized")

    def printCOORDS(self, f, p, fstr):
        MAXL = 132
        # if self.skip:
        #    self.skip -= 1
        #    return fstr
        for point in p:
            up = " %2.2f" % (point)
            if len(fstr) + len(up) > MAXL:
                f.write(fstr + "\n")
                fstr = " "
            fstr += up
        return fstr

    def printAC(self, f, p, N, fstr):
        MAXL = 132
        if N == 1:
            up = " %i" % (p)
        else:
            up = " %i*%i" % (N, p)
        if len(fstr) + len(up) > MAXL:
            f.write(fstr + "\n")
            fstr = " "
        fstr += up
        return fstr

    def printPROP(self, f, p, N, fstr):
        MAXL = 132
        if N == 1:
            # up = " %1.4e" %(p) # standard notation
            up = " %1.4e" % (p)  # scientific notation
        else:
            up = " %i*%1.4e" % (N, p)
            # up = " %i*%1.4e" %(N,p) # scientific notation
        if len(fstr) + len(up) > MAXL:
            f.write(fstr + "\n")
            fstr = " "
        fstr += up
        return fstr

    def exportTOUGH2(self, fname):
        """Saves the grid as a fixed format TOUGH(2) grid.
        """
        STR = "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789"
        self.ne, self.nn, self.nz = np.array(self.Grid.GetDimensions())  # - 1 #
        filename, ext = os.path.splitext(fname)
        if self.GridType == "vtkStructuredGrid":
            with io.open(filename, 'w', newline='\r\n') as f:
                f.write("ELEME")
                # debug
                f.write(
                    """
                    1        10        20        30        40        50        60        70        80
                    |--------|---------|---------|---------|---------|---------|---------|---------|
                    12345678901234567890123456789012345678901234567890123456789012345678901234567890
                    """)

                ii = 0
                for iy in range(self.nn):
                    for ix in range(self.ne):
                        # f.write(str(iy)+str(ix)+"\n")
                        # first base
                        b2 = ii // (len(STR) * len(STR))
                        b1 = (ii - len(STR) * b2) // len(STR)
                        b0 = ii % len(STR)

                        f.write(STR[b2] + STR[b1] + STR[b0] + "\t" + str(ii) + "\n")
                        ii += 1

    def exportECL(self, fname):
        """ Saves the grid as an ECLIPSE grid. For the purposes of ECLIPSE
        """

        # TODO add consistency of dimensions across the inputs
        self.ne, self.nn, self.nz = np.array(self.Grid.GetDimensions()) - 1  # ECLIPSE
        filename, ext = os.path.splitext(fname)
        if self.GridType == "vtkStructuredGrid":
            with io.open(filename + ".GRDECL", 'w', newline='\r\n') as f:
                f.write('-- Generated [\n')
                f.write('-- Format      : ECLIPSE keywords (grid geometry and properties) (ASCII)\n')
                # f.write('-- Exported by : Petrel 2013.7 (64-bit) Schlumberger\n'
                f.write('-- Exported by : ReGrid v.' + version + "\n")
                f.write('-- User name   : ' + getpass.getuser() + "\n")
                f.write('-- Date        : ' + datetime.now().strftime("%A, %B %d %Y %H:%M:%S") + "\n")
                f.write('-- Project     : ' + "ReGrid project\n")
                f.write('-- Grid        : ' + "Description\n")
                f.write('-- Generated ]\n\n')

                f.write('SPECGRID                               -- Generated : ReGrid\n')
                f.write('  %i %i %i 1 F /\n\n' % (self.ne, self.nn, self.nz))
                f.write('COORDSYS                               -- Generated : ReGrid\n')
                f.write('  1 4 /\n\n')  # what is this line?

                f.write('COORD                                  -- Generated : ReGrid\n')
                nz = self.nz
                fstr = str(" ")

                for iy in range(self.nn):
                    for ix in range(self.ne):
                        p0 = self.Grid.GetCell(ix, iy, 0).GetPoints().GetPoint(0)
                        fstr = self.printCOORDS(f, p0, fstr)
                        p1 = self.Grid.GetCell(ix, iy, nz - 1).GetPoints().GetPoint(4)
                        fstr = self.printCOORDS(f, p1, fstr)
                    # outside edge on far x
                    p2 = self.Grid.GetCell(ix, iy, 0).GetPoints().GetPoint(1)
                    fstr = self.printCOORDS(f, p2, fstr)
                    p3 = self.Grid.GetCell(ix, iy, nz - 1).GetPoints().GetPoint(5)
                    fstr = self.printCOORDS(f, p3, fstr)
                # outside edge on far y
                for ix in range(self.ne):
                    p8 = self.Grid.GetCell(ix, iy, 0).GetPoints().GetPoint(3)
                    fstr = self.printCOORDS(f, p8, fstr)
                    p9 = self.Grid.GetCell(ix, iy, nz - 1).GetPoints().GetPoint(7)
                    fstr = self.printCOORDS(f, p9, fstr)
                # outside edge on far northeast
                p14 = self.Grid.GetCell(ix, iy, 0).GetPoints().GetPoint(2)
                fstr = self.printCOORDS(f, p14, fstr)
                p15 = self.Grid.GetCell(ix, iy, nz - 1).GetPoints().GetPoint(6)
                fstr = self.printCOORDS(f, p15, fstr)
                f.write(fstr)
                fstr = " "
                f.write(" /")
                f.write("\n")
                f.write("\n")

                f.write('ZCORN                                  -- Generated : ReGrid\n')
                for iz in range(self.nz):
                    for iy in range(self.nn):
                        # front face
                        for ix in range(self.ne):
                            p0 = self.Grid.GetCell(ix, iy, iz).GetPoints().GetPoint(0)
                            p1 = self.Grid.GetCell(ix, iy, iz).GetPoints().GetPoint(1)
                            fstr = self.printCOORDS(f, [p0[2]], fstr)
                            fstr = self.printCOORDS(f, [p1[2]], fstr)
                        # back face
                        for ix in range(self.ne):
                            p0 = self.Grid.GetCell(ix, iy, iz).GetPoints().GetPoint(3)
                            p1 = self.Grid.GetCell(ix, iy, iz).GetPoints().GetPoint(2)
                            fstr = self.printCOORDS(f, [p0[2]], fstr)
                            fstr = self.printCOORDS(f, [p1[2]], fstr)
                    # bottom layer
                    for iy in range(self.nn):
                        # front face
                        for ix in range(self.ne):
                            p0 = self.Grid.GetCell(ix, iy, iz).GetPoints().GetPoint(4)
                            p1 = self.Grid.GetCell(ix, iy, iz).GetPoints().GetPoint(5)
                            fstr = self.printCOORDS(f, [p0[2]], fstr)
                            fstr = self.printCOORDS(f, [p1[2]], fstr)
                        # back face
                        for ix in range(self.ne):
                            p0 = self.Grid.GetCell(ix, iy, iz).GetPoints().GetPoint(7)
                            p1 = self.Grid.GetCell(ix, iy, iz).GetPoints().GetPoint(6)
                            fstr = self.printCOORDS(f, [p0[2]], fstr)
                            fstr = self.printCOORDS(f, [p1[2]], fstr)
                f.write(fstr)
                fstr = " "
                f.write(" /")
                f.write("\n")
                f.write("\n")
                f.write('ACTNUM                                 -- Generated : ReGrid\n')

                c = -999
                N = 0
                for iac in self.ActiveCells.flatten(order='F'):
                    if iac == c:
                        N += 1
                    else:
                        if c != -999:
                            fstr = self.printAC(f, c, N, fstr)
                        c = iac
                        N = 1
                fstr = self.printAC(f, c, N, fstr)
                f.write(fstr)
                f.write(" /")
                f.write("\n")
                f.write("\n")
        else:
            print("Only structured grids can be converted to ECLIPSE files")

    def exportECLPropertyFiles(self, fname):
        """ Convert any point data to cell data
        """

        # Convert point data to cell data for output
        # verifying if this is necessary or if ECLIPSE can use point attributes
        pointConvert = True
        if pointConvert:
            p2c = vtk.vtkPointDataToCellData()
            p2c.SetInputDataObject(self.Grid)
            p2c.PassPointDataOn()
            p2c.Update()
            self.Grid = p2c.GetOutput()

        filename, ext = os.path.splitext(fname)
        for ia in range(self.Grid.GetCellData().GetNumberOfArrays()):
            prop = self.Grid.GetCellData().GetArray(ia).GetName()
            print("exporting prop", prop)
            if self.GridType == "vtkStructuredGrid":
                with io.open(filename + "prop-" + prop.lower() + ".GRDECL", 'w', newline='\r\n') as f:
                    f.write('-- Generated [\n')
                    f.write('-- Format      : ECLIPSE keywords (grid properties) (ASCII)\n')
                    f.write('-- Exported by : ReGrid v.' + version + "\n")
                    f.write('-- User name   : ' + getpass.getuser() + "\n")
                    f.write('-- Date        : ' + datetime.now().strftime("%A, %B %d %Y %H:%M:%S") + "\n")
                    f.write('-- Project     : ' + "ReGrid project\n")
                    f.write('-- Grid        : ' + "Description\n")
                    f.write('-- Unit system : ' + "ECLIPSE-Field\n")
                    f.write('-- Generated ]\n\n')

                    f.write(prop.upper() + '                                 -- Generated : ReGrid\n')
                    f.write('-- Property name in Petrel : ' + prop + '\n')

                    c = -999.9999
                    N = 0
                    ii = 0
                    fstr = " "
                    for iz in range(self.nz):
                        for iy in range(self.nn):
                            for ix in range(self.ne):
                                # iac = round(self.Grid.GetCellData().GetArray(ia).GetTuple1(ii), 4)
                                iac = '{:0.4e}'.format(self.Grid.GetCellData().GetArray(ia).GetTuple1(ii))
                                print(iac)
                                ii += 1
                                if iac == c:
                                    N += 1
                                else:
                                    if c != -999.9999:
                                        fstr = self.printPROP(f, c, N, fstr)
                                    c = eval(iac)
                                    N = 1
                    fstr = self.printPROP(f, c, N, fstr)
                    f.write(fstr)
                    f.write(" /")
                    f.write("\n")


class GRDECL(FlowGrid):
    """ GRDECL processes Schlumberger ECLIPSE files
    """

    def __init__(self):
        super(GRDECL, self).__init__()
        nx, ny, nz = 0, 0, 0

    def loadNodes(self, fname):
        """
            Reads I, J(max), K
                  iterates through I, then decriments J, increments K
                  I = easting
                  J = northing
                  K = depth or elevation?
        """
        with open(fname, "r") as fp:

            # Read in the header
            for line in fp:
                item = line.split()
                if len(item) > 0:
                    if item[0] == "SPECGRID":
                        self.SPECGRID = np.array(fp.readline().split()[0:3], dtype=int)
                    if item[0] == "COORDSYS":
                        self.COORDSYS = fp.readline().split()
                    if item[0] == "COORD":
                        break

            # Read in the coordinates
            self.coords = []
            for line in fp:
                if line.split()[-1] != "/":
                    item = line.split()
                    for c in item:
                        if '*' in c:
                            cc = c.split('*')
                            for i in range(int(cc[0])):
                                self.coords.append(cc[-1])
                        else:
                            self.coords.append(c)
                else:
                    if len(line.split()) > 1:
                        item = line.split()
                        for i in range(len(item) - 1):
                            cc = item[i]
                            if '*' in cc:
                                ccc = cc.split('*')
                                for j in range(int(ccc[0])):
                                    self.coords.append(ccc[-1])
                            else:
                                self.coords.append(c)
                        break
                    else:
                        break

            # Read in ZCORN
            self.zcorn = []
            i = 0
            for line in fp:
                item = line.split()
                if len(item) > 0:
                    if item[0] == "ZCORN":
                        for line in fp:
                            if line.split():
                                if line.split()[-1] != "/":
                                    self.zcorn += line.split()
                                else:
                                    self.zcorn += line.split()[0:-1]
                                    break
                if len(self.zcorn) > 0:
                    break

            # Read in (in)active cells
            self.active = []
            for line in fp:
                item = line.split()
                if len(item) > 0:
                    if item[0] == "ACTNUM":
                        for line in fp:
                            if line.split():
                                if line.split()[-1] != "/":
                                    c = line.split()
                                    if '*' in c:
                                        cc = c.split('*')
                                        for i in range(float(cc[0])):
                                            self.active += cc[-1]
                                    else:
                                        self.active += c
                                else:
                                    self.active += line.split()[0:-1]
                                    break

        self.coords = np.array(self.coords, dtype=float)
        print(self.coords)

        # In Petrel...
        self.ne = self.SPECGRID[0]  # x  i
        self.nn = self.SPECGRID[1]  # y  j
        self.nz = self.SPECGRID[2]  # z  k

        # build grid
        self.buildGrid(plot=False)
        self.buildActiveCells(plot=False)
        self.buildZGrid(plot=False)
        # self.calculateVolumes(plot=False)
        #
        # Convert to VTK
        self.GridType = "vtkStructuredGrid"
        self.Grid = vtk.vtkStructuredGrid()
        self.Grid.SetDimensions(self.ne+1, self.nn+1, self.nz+1)
        vtk_points = vtk.vtkPoints()
        ve = 1.

        for iz in range(self.nz):
            if iz == 0:
                for iy in range(self.nn+1):
                    for ix in range(self.ne+1):
                        vtk_points.InsertNextPoint( self.X0[ix,iy], \
                                                    self.Y0[ix,iy], \
                                               ve * self.ZZT[iz][ix,iy] )
            for iy in range(self.nn+1):
                for ix in range(self.ne+1):
                    vtk_points.InsertNextPoint( self.X0[ix,iy], \
                                                self.Y0[ix,iy], \
                                           ve * self.ZZB[iz][ix,iy] )
        self.Grid.SetPoints(vtk_points)

        # Add in active cells
        ac = vtk.vtkIntArray()
        ac.SetName( "ActiveCells" )
        for iac in self.ActiveCells.flatten( order='F' ):
            ac.InsertNextTuple1( iac )
        self.Grid.GetCellData().AddArray(ac)

    def buildGrid(self, plot=False):
        """
        Topology of COORD mesh, only describes first layer..

                  8--------10-------12-------14
                 /|       /|       /|       /|
                / |      / |      / |      / |
               0--------2--------4--------6  |
               |  9-----|--11----|--13----|--15
               | /      | /      | /      | /
               |/       |/       |/       |/
               1--------3--------5--------7            7  -->   (2*(NE+1))
                                                      15  -->   (2*(NE+1)*(NN+1))
        """

        print("Constructing grid")
        # print("Grid dims", self.ne, self.nn, self.nz)
        # print("Num points", 2*(self.ne+1)*(self.nn+1)*3, len(self.coords))

        # number of edges
        self.ndx = self.ne + 1
        self.ndy = self.nn + 1
        self.ndz = self.nz + 1

        # extract the triplets
        self.points = {}
        self.points["e"] = self.coords[0::3]
        self.points["n"] = self.coords[1::3]
        self.points["z"] = self.coords[2::3]

        print('points e')
        print(self.points["e"])

        # Here are the coordinates
        self.X0 = np.reshape(self.points["e"][0::2] , (self.ndx,self.ndy), order="F")
        self.Y0 = np.reshape(self.points["n"][0::2] , (self.ndx,self.ndy), order="F")
        self.Z0 = np.reshape(self.points["z"][0::2] , (self.ndx,self.ndy), order="F")

        self.X1 = np.reshape(self.points["e"][1::2] , (self.ndx,self.ndy), order="F")
        self.Y1 = np.reshape(self.points["n"][1::2] , (self.ndx,self.ndy), order="F")
        self.Z1 = np.reshape(self.points["z"][1::2] , (self.ndx,self.ndy), order="F")
        #
        # # visualize
        # if plot:
        #     print("plotting")
        #     fig = plt.figure()
        #     ax = fig.add_subplot(111, projection='3d')
        #     ax.plot_wireframe(f2m*self.X0, f2m*self.Y0, f2m*self.Z0, rstride=1, cstride=1)
        #     ax.plot_wireframe(f2m*self.X1, f2m*self.Y1, f2m*self.Z1, rstride=1, cstride=1)
        #     plt.show()

    def buildZGrid(self, plot=False):
        """
            Petrel provides the ZCORN in a truly arcane ordering--it's awful--and really, the programmers
            deserve a special place in hell for doing this. The ordering is as follows, for a given plane:

             29    36  30   37 31    38 32    39 33    40 34    41 35    42
              _______  _______  ______  _______  _______  _______  _______
             /      / /      / /     / /      / /      / /      / /      /|
            /      / /      / /     / /      / /      / /      / /      / |
           00----01 02----03 04----05 06----07 08----09 10----11 12----13 /
            |  A  | |  B   | |   C  | |   D  | |   E  | |  F   | |   G  |/
           14----15 16----17 18----19 20----21 22----23 24----25 26----27


            This pattern is then repeated for each depth layer, it isn't that clear, but my ASCII art skills
            are already sufficiently challenged.

        """

        print("Constructing Z corners")

        # self.zcorn = np.array(self.zcorn, dtype=float)
        # temp = np.zeros( ((self.ne+1)*(self.nn+1)*self.nz) )
        temp = []
        count = 0
        for item in self.zcorn:

            if "*" in item:
                ct = (int)(item.split("*")[0])
                vl = (float)(item.split("*")[1])
                temp += np.tile(vl, ct).tolist()
                count += ct
            else:
                temp += [(float)(item)]
                count += 1

        # layers = np.resize(temp, (8, self.ne*self.nn*self.nz ))
        layers = np.resize(temp, (self.nz * 2, self.ne * self.nn * 4))
        """
        plt.plot(newtemp[0,:])                    # TOP     0    0
        plt.plot(newtemp[1,:])       # SAME --    # BOTTOM  0    1
        #plt.plot(newtemp[2,:])      # SAME --    # TOP     1    2

        plt.plot(newtemp[3,:])       # SAME --    # BOTTOM  1    3
        #plt.plot(newtemp[4,:])      # SAME --    # TOP     2    4

        plt.plot(newtemp[5,:])       # SAME --    # BOTTOM  2    5
        #plt.plot(newtemp[6,:])      # SAME --    # TOP     3    6
        plt.plot(newtemp[7,:])                    # BOTTOM  3    7
        """
        self.ZZT = {}  # zztop ha ha...two year's later this is still funny -TI
        self.ZZB = {}
        for ilay in range(self.nz):
            self.ZZT[ilay] = np.zeros((self.ndx, self.ndy))
            self.ZZB[ilay] = np.zeros((self.ndx, self.ndy))
            iis = 0
            # plt.plot(layers[ilay*2])
            for iin in range(self.nn):
                nears = {}
                fars = {}
                bnears = {}
                bfars = {}
                for iif in range(2):
                    # top
                    nears[iif] = layers[ilay * 2][iis:iis + 2 * self.ne][0::2].tolist()
                    fars[iif] = layers[ilay * 2][iis:iis + 2 * self.ne][1::2].tolist()
                    layers[ilay * 2][iis:iis + 2 * self.ne][0::2] *= 0.  # check
                    layers[ilay * 2][iis:iis + 2 * self.ne][1::2] *= 0.
                    nears[iif].append(fars[iif][-1])
                    fars[iif] = [nears[iif][0]] + fars[iif]
                    # bottom
                    bnears[iif] = layers[ilay * 2 + 1][iis:iis + 2 * self.ne][0::2].tolist()
                    bfars[iif] = layers[ilay * 2 + 1][iis:iis + 2 * self.ne][1::2].tolist()
                    layers[ilay * 2 + 1][iis:iis + 2 * self.ne][0::2] *= 0.
                    layers[ilay * 2 + 1][iis:iis + 2 * self.ne][1::2] *= 0.
                    bnears[iif].append(bfars[iif][-1])
                    bfars[iif] = [bnears[iif][0]] + bfars[iif]
                    #
                    iis += 2 * self.ne

                self.ZZT[ilay][:, iin] = nears[0]
                self.ZZB[ilay][:, iin] = bnears[0]
                # NaN mask for visualizing, but can be sort of a pain to deal with
                # imask = np.nonzero( 1-self.ActiveCells[:,iin,ilay] )
                # self.ZZT[ilay][:,iin][1::][imask] = np.nan
                # self.ZZB[ilay][:,iin][1::][imask] = np.nan
                # if self.ActiveCells[0,iin,ilay] == 0:
                # self.ZZT[ilay][:,iin][0]  = np.nan
                # self.ZZB[ilay][:,iin][0]  = np.nan
                if iin == self.nn - 1:
                    self.ZZT[ilay][:, iin + 1] = fars[1]
                    self.ZZB[ilay][:, iin + 1] = bfars[1]
                    # NaN mask
                    # self.ZZT[ilay][:,iin+1][1::][imask] = np.nan
                    # self.ZZB[ilay][:,iin+1][1::][imask] = np.nan
                    # if self.ActiveCells[0,iin,ilay] == 0:
                    #    self.ZZT[ilay][:,iin+1][0]  = np.nan
                    #    self.ZZB[ilay][:,iin+1][0]  = np.nan

        print("Layers ||", np.linalg.norm(layers), "||")
        # exit()

        # visualize
        if plot:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            # ax.plot_wireframe( self.X0, self.Y0, self.Z0, rstride=1, cstride=1)

            ax.plot_wireframe(self.X0, self.Y0, self.ZZT[0], rstride=1, cstride=1, color="blue")
            # ax.plot_wireframe( self.X0, self.Y0, self.ZZT[1], rstride=1, cstride=1, color="blue")
            # ax.plot_wireframe( self.X0, self.Y0, self.ZZT[2], rstride=1, cstride=1, color="blue")
            # ax.plot_wireframe( self.X0, self.Y0, self.ZZT[3], rstride=1, cstride=1, color="blue")

            # ax.plot_wireframe( self.X0, self.Y0, self.ZZB[3], rstride=1, cstride=1, color="green")

            plt.gca().set_xlim(np.min(self.X0), np.max(self.X0))
            plt.gca().set_ylim(np.max(self.Y0), np.min(self.Y0))
            # plt.gca().set_zlim( np.max(self.ZZB[3]),  np.min(self.ZZT[0]) )
            plt.gca().set_zlim(5000, 4000)
            plt.savefig("mesh.png")
            plt.show()

    def buildActiveCells(self, plot=False):

        print("Constructing active cells")
        self.ActiveCells = np.zeros((self.ne * self.nn * self.nz), dtype=int)

        count = 0
        for item in self.active:
            if "*" in item:
                ct = (int)(item.split("*")[0])
                vl = (int)(item.split("*")[1])
                self.ActiveCells[count:count + ct] = vl
                count += ct
            else:
                self.ActiveCells[count] = (int)(item)
                count += 1

        self.ActiveCells = np.reshape(self.ActiveCells, (self.ne, self.nn, self.nz), order="F")

        if plot:
            plt.pcolor(self.X0.T, self.Y0.T, self.ActiveCells[:, :, 0].T, edgecolors='w', linewidths=.1)
            plt.xlabel("easting")
            plt.ylabel("northing")
            plt.gca().set_xlim(np.min(self.X0), np.max(self.X0))
            plt.gca().set_ylim(np.max(self.Y0), np.min(self.Y0))
            plt.gca().xaxis.tick_top()
            plt.gca().xaxis.set_label_position("top")
            plt.show()

    def calculateVolumes(self, plot=False):
        # Iterate over cells, assert that we are dealing with parallelpiped, if so
        #             | u1    u2   u3 |
        #    A = det  | v1    v2   v3 |
        #             | w1    w2   w3 |
        # self.Volumes = 10000*np.random.normal(0,1, (self.ne, self.nn, self.nz) )
        self.Volumes = np.zeros((self.ne, self.nn, self.nz))
        for iiz in range(self.nz):
            for iie in range(self.ne):
                for iin in range(self.nn):

                    if self.ActiveCells[iie, iin, iiz]:

                        u = np.array((self.X0[iie, iin], self.Y0[iie, iin], self.ZZT[iiz][iie, iin])) - \
                            np.array((self.X0[iie + 1, iin], self.Y0[iie + 1, iin], self.ZZT[iiz][iie, iin]))

                        v = np.array((self.X0[iie, iin], self.Y0[iie, iin], self.ZZT[iiz][iie, iin])) - \
                            np.array((self.X0[iie, iin + 1], self.Y0[iie, iin + 1], self.ZZT[iiz][iie, iin]))

                        w = np.array((self.X0[iie, iin], self.Y0[iie, iin], self.ZZT[iiz][iie, iin])) - \
                            np.array((self.X0[iie, iin], self.Y0[iie, iin], self.ZZB[iiz][iie, iin]))
                        if np.any(u != u) or np.any(v != v) or np.any(w != w):
                            print("NAN!", iie, iin, iiz)
                            exit()
                        V = np.linalg.det(np.array((f2m * u, f2m * v, f2m * w)))
                        self.Volumes[iie, iin, iiz] = np.abs(V)  # in m^3

        vr = ((3. / (4. * np.pi)) * self.Volumes) ** (1. / 3.)  # virtual radius, taking into account porosity

        print("Total grid volume: " + str(np.sum(self.Volumes)) + " m^3")

    def readProperty(self, fname, attr_name):
        """ Reads a single property from a file
        """
        temp = []
        with open(fname, "r") as fp:
            for line in fp:
                item = line.split()
                if len(item) > 0:
                    if item[0] != "--":
                        tag = item[0]
                        break

            for line in fp:
                attribute = line.split()
                if attribute:
                    if attribute[0] != "--":
                        if attribute[-1] != "/":
                            for c in attribute:
                                if '*' in c:
                                    cc = c.split('*')
                                    for i in range(int(cc[0])):
                                        temp.append(cc[-1])
                                else:
                                    temp.append(c)
                        else:
                            attribute.pop()
                            for c in attribute:
                                if '*' in c:
                                    cc = c.split('*')
                                    for i in range(int(cc[0])):
                                        temp.append(cc[-1])
                                else:
                                    temp.append(c)
                            break

                        # #attribute = fp.readline().split()[-1]
                        # attribute = fp.readline().split()
                        # print(attribute)
                        # if attribute[0] != "--":
                        #     self.Prop[tag] = attribute
                        # print("loading", attribute)
                        # for line in fp:
                        #     if line.split():
                        #         if line.split()[0] != "--":
                        #             if line.split()[-1] != "/":
                        #                 temp += line.split()
                        #             else:
                        #                 temp += line.split()[0:-1]
                        #                 break
        print(temp)
        data = np.zeros((self.ne * self.nn * self.nz), dtype=float)
        count = 0
        for item in temp:
            if "*" in item:
                ct = (int)(item.split("*")[0])
                vl = (float)(item.split("*")[1])
                data[count:count + ct] = vl
                count += ct
            else:
                data[count] = (float)(item)
                count += 1

        data = np.reshape(data, (self.ne, self.nn, self.nz), order="F")

        # Add to VTK grid
        ac = vtk.vtkDoubleArray()
        ac.SetName(attr_name)
        for iac in data.flatten(order='F'):
            ac.InsertNextTuple1(iac)
        self.Grid.GetCellData().AddArray(ac)

        return data

    def readOutputProperty(self, fname, prop_strings, toVTK=True, toNumpy=True):
        """
        Reads per-cell properties from .PRT file for all timesteps

        Pack property keywords to read from .PRT into list of lists
        This saves significant time for large grids as only one pass is required through .PRT file
        Inner lists should contain keywords that enable line denoting property section to be uniquely identified

        :param prop_strings: [[ECL prop keyword, subkey1, subkey2, ...], ...]

        Order props as they appear in .PRT file
        """
        print('Reading output properties\n')
        prop = {}
        prop_idx = 0
        for p in prop_strings:
            prop[p[0]] = {}
        with open(fname, "r") as fp:
            t = 0
            II = []
            build = False
            data = np.zeros(self.ne * self.nn * self.nz)
            for line in fp:
                # Find prop keywords
                if not build:
                    if all(e in line for e in prop_strings[prop_idx]):
                        data = np.zeros(self.ne * self.nn * self.nz)
                        # print('Reading output property: ' + prop_strings[prop_idx][0])
                        # print('t = ' + str(t))
                        build = True
                # Read prop data
                else:
                    item = line.split()
                    if len(item) > 0:
                        if 'I=' in line:
                            II = line.split('I=')[1].split()
                            II = list(map(int, II))
                        elif '(*,' in item[0]:
                            idxs = line.split('(')[1].split(')')[0].replace(',', ' ').split()
                            J = int(idxs[1])
                            K = int(idxs[2])
                            vals = line.split(')')[1].split()
                            for c,I in enumerate(II):
                                if '-' in vals[c]:
                                    vals[c] = 0
                                idx = ((self.ne * self.nn) * (K - 1)) + (self.ne * (J - 1)) + (I - 1)
                                data[idx] = vals[c]
                        elif '--' in item[0]:
                            build = False
                            pname = prop_strings[prop_idx][0]
                            prop[pname] = copy.deepcopy(data)
                            if prop_idx < len(prop_strings) - 1:
                                prop_idx += 1
                            # All properties for current t have been read
                            else:
                                print('Exporting grid for t = ' + str(t))
                                ids = []
                                for pn in prop.keys():
                                    data = prop[pn]
                                    if toNumpy:
                                        if t == 0:
                                            self.checkOutputDir(pn)
                                        grid_data = np.reshape(data, (self.ne, self.nn, self.nz), order="F")
                                        np.save(self.outputDir + '/' + pn + '/' + pn + '_' + str(t), grid_data)
                                    if toVTK:
                                        if t == 0:
                                            self.checkOutputDir('vtk')
                                        ac = vtk.vtkDoubleArray()
                                        ac.SetName(pn)
                                        for iac in data:
                                            ac.InsertNextTuple1(iac)
                                        id = self.Grid.GetCellData().AddArray(ac)
                                        ids.append(id)
                                if toVTK:
                                    self.exportVTK(self.outputDir + '/vtk/' + fname.split('/')[-1].split('.')[0] + '_' + str(t))
                                    for id in ids:
                                        self.Grid.GetCellData().RemoveArray(id)
                                prop_idx = 0
                                t += 1

    def readWellOutput(self, fname, keys):
        wellOutput = {}
        keyOrder = {}
        readNames = False
        build = False
        skip = 0

        for key in keys:
            wellOutput[key] = {}

        with open(fname, "r") as fp:
            for line in fp:
                item = line.split()
                if len(item) > 0:

                    # read time series values
                    if build:
                        # I don't know why '1' denotes timestep end in .RSM file
                        if item[0] == '1':
                            build = False
                            continue
                        t = item[0]
                        for idx in keyOrder.keys():
                            if t in wellOutput[keyOrder[idx][0]]:
                                wellOutput[keyOrder[idx][0]][t][keyOrder[idx][1]] = float(item[idx])
                            else:
                                if len(keyOrder[idx]) > 1:
                                    wellOutput[keyOrder[idx][0]][t] = {}
                                    wellOutput[keyOrder[idx][0]][t][keyOrder[idx][1]] = float(item[idx])
                                else:
                                    wellOutput[keyOrder[idx][0]][t] = float(item[idx])
                        continue

                    # get names of wells
                    if readNames:
                        for idx in keyOrder.keys():
                            curr = keyOrder[idx]
                            if keyOrder[idx][0][0] == 'W':
                                keyOrder[idx].append(item[idx - skip])
                        next(fp)
                        next(fp)
                        readNames = False
                        build = True
                        skip = 0
                        continue

                    # find line that contains keys
                    if item[0] == 'TIME':
                        keyOrder = {}
                        # if time found, then keys might be on this same line
                        for j,key in enumerate(keys):
                            for i,var in enumerate(item):
                                if key == var:
                                    if i in keyOrder:
                                        keyOrder[i].append(key)
                                    else:
                                        keyOrder[i] = [key]
                                if j == 0 and var[0] != 'W':
                                    # then it is time related or a field variable
                                    skip += 1
                        if len(keyOrder.keys()) > 0:
                            next(fp)
                            readNames = True
        return wellOutput

    # Exports Numpy array of property (can be wells)
    # TODO: inherit from FlowGrid
    def exportProp(self, title, prop):
        if not os.path.exists(self.outputDir):
            os.mkdir(self.outputDir)
        np.save(self.outputDir + '/' + title, prop)




class SUTRA(FlowGrid):
    """ SUTRA is a USGS flow modelling code.
    """

    def __init__(self):
        super(SUTRA, self).__init__()
        nx, ny, nz = 0, 0, 0

    def loadNodes(self, fname, nx, ny, nz, ve=-1):
        """ Reads in the points of the grid, ususally in a file called nodewise
            fname = nodes file
            nx = number of cells in the easting(x) direction
            ny = number of cells in the northing (y) direction
            nz = number of cells in depth, positive up
            ve = vertical exaggeration, default is 1 (none)
            This method results in the generation of a VtkStructuredGrid
        """
        self.nx = nx
        self.ny = ny
        self.nz = nz

        self.ActiveCells = np.ones((self.nx * self.ny * self.nz), dtype=int)

        X = np.loadtxt(fname, comments="#")
        self.points = np.reshape(np.array((X[:, 2], X[:, 3], X[:, 4])).T, (nx, ny, nz, 3))

        # work directly with VTK structures
        self.GridType = "vtkStructuredGrid"
        self.Grid = vtk.vtkStructuredGrid()
        self.Grid.SetDimensions(nx, ny, nz)
        vtk_points = vtk.vtkPoints()
        for iz in range(nz):
            for iy in range(ny):
                for ix in range(nx):
                    vtk_points.InsertNextPoint(self.points[ix, iy, iz][0], \
                                               self.points[ix, iy, iz][1], \
                                               ve * self.points[ix, iy, iz][2])
        self.Grid.SetPoints(vtk_points)

    def loadNodesConnections(self, nodes, connections):
        """ In contrast to the above method, the points and connections can be loaded instead.
            For non-regular grids this is necessary. This method results in the generation
            of a vtkUnstructuredGrid.
            nodes = node file, often called nodewise
            connections = element connections, often called incident
        """
        X = np.loadtxt(nodes, comments="#")
        #                   x       y       z
        points = np.array((X[:, 2], X[:, 3], X[:, 4])).T

        self.GridType = "vtkUnstructuredGrid"
        self.Grid = vtk.vtkUnstructuredGrid()
        vtk_points = vtk.vtkPoints()

        for point in range(np.shape(points)[0]):
            vtk_points.InsertNextPoint(points[point, 0], points[point, 1], points[point, 2])
        self.Grid.SetPoints(vtk_points)

        # Read in the connections, the format is as follows
        #  nodeid    p0, p1, p2, p3, p4, p5, p6, p7, p8
        C = np.loadtxt(connections, comments="#", skiprows=2, dtype=int)
        for line in range(np.shape(C)[0]):
            idList = vtk.vtkIdList()
            for node in C[line, :][1:]:
                idList.InsertNextId(node - 1)
            self.Grid.InsertNextCell(vtk.VTK_HEXAHEDRON, idList)

    def readPermeability(self, fname, label=("$\kappa_x$", "$\kappa_y$", "$\kappa_z$")):
        """ Reads in SUTRA permeability data
        """
        k = np.loadtxt(fname, comments="#")
        nr, nc = np.shape(k)
        if self.GridType == "vtkStructuredGrid":
            # Sutra and VTK use opposite ordering
            k = np.reshape(k, (self.nx - 1, self.ny - 1, self.nz - 1, np.shape(k)[1]))
            k = np.reshape(k, (nr, nc), order='F')
        kx = vtk.vtkDoubleArray()
        kx.SetName(label[0])
        ky = vtk.vtkDoubleArray()
        ky.SetName(label[1])
        kz = vtk.vtkDoubleArray()
        kz.SetName(label[2])
        for ik, K in enumerate(k):
            kx.InsertNextTuple1(K[2])
            ky.InsertNextTuple1(K[3])
            kz.InsertNextTuple1(K[4])
        self.Grid.GetCellData().AddArray(kx)
        self.Grid.GetCellData().AddArray(ky)
        self.Grid.GetCellData().AddArray(kz)

    def readPorosity(self, fname, label="phi"):  # LaTeX tags work too: $\phi$
        phi = np.loadtxt(fname)
        nr, nc = np.shape(phi)
        if self.GridType == "vtkStructuredGrid":
            # Sutra and VTK use opposite ordering
            phi = np.reshape(phi, (self.nx, self.ny, self.nz, np.shape(phi)[1]))
            phi = np.reshape(phi, (nr, nc), order='F')
        vphi = vtk.vtkDoubleArray()
        vphi.SetName(label)
        for ik, K in enumerate(phi):
            vphi.InsertNextTuple1(K[5])
        self.Grid.GetPointData().AddArray(vphi)

    def readPressure(self, fname, ts=2, label="$P$"):
        nnodes = self.nx * self.ny * self.nz
        P = np.loadtxt(fname, comments="#")[ts * nnodes:(ts + 1) * nnodes, :]
        C = np.loadtxt(fname, comments="#")[ts * nnodes:(ts + 1) * nnodes, :]
        nr, nc = np.shape(P)
        if self.GridType == "vtkStructuredGrid":
            # Sutra and VTK use opposite ordering
            P = np.reshape(P, (self.nx, self.ny, self.nz, np.shape(P)[1]))
            P = np.reshape(P, (nr, nc), order='F')
            C = np.reshape(C, (self.nx, self.ny, self.nz, np.shape(C)[1]))
            C = np.reshape(C, (nr, nc), order='F')
        vP = vtk.vtkDoubleArray()
        vP.SetName(label)

        vC = vtk.vtkDoubleArray()
        vC.SetName("Concentration")

        for ik in range(nnodes):
            vP.InsertNextTuple1(P[ik, 3])
            vC.InsertNextTuple1(C[ik, 4])
            # vP.InsertNextTuple1( P[2*nnodes+ik, 3] )a
        self.Grid.GetPointData().AddArray(vP)
        self.Grid.GetPointData().AddArray(vC)


# ========================================
# ===========================================
# Class for converting CMG grid to VTK grid
# ===========================================
# ========================================
class CMG(FlowGrid):
    def __init__(self):
        super(CMG, self).__init__()
        self.outputDir = 'output'

    def setOutputDir(self, dir):
        self.outputDir = dir

    # Get xyz coords of corner points defining cell
    def getCellCoords(self, cellIdxs):
        coordList = []
        xyz = [0, 0, 0]
        cell = self.Grid.GetCell(cellIdxs[0], cellIdxs[1], cellIdxs[2])
        pointIds = cell.GetPointIds()
        numIds = pointIds.GetNumberOfIds()
        for n in range(numIds):
            p = pointIds.GetId(n)
            self.Grid.GetPoint(p, xyz)
            coordList.append(copy.deepcopy(xyz))
        return coordList

    # Compute center of cell given corner point coords
    def centroid(self, coords):
        return np.mean(coords, axis=0)

    # Builds a corner point grid from a CMG formatted input file (.dat)
    # CMG output files (.out) do not always contain complete CP grid information
    def buildCorner(self, fname):
        print('Building corner point grid')
        self.iWidths = []
        self.jWidths = []

        # Reads DI or DJ, where dir is 'I' or 'J'
        def readBlockSizes(dir):
            count = 0
            if dir == "I":
                widths = self.iWidths
                nb = 0
            else:
                widths = self.jWidths
                nb = 1
            for line in fp:
                item = line.split()
                if item[0] == "D" + dir or item[0] == "*D" + dir:
                    if item[1] == dir + "VAR" or item[1] == dir + "*VAR":
                        self.iWidths += item[2:]
                        count += len(item) - 2
                        break
                    # Haven't added support for this yet
                    elif item[1] == "CON":
                        pass
            for line in fp:
                item = line.split()
                for zz in item:
                    if "*" in zz:
                        item = zz.split("*")
                        for i in range(0, int(item[0])):
                            widths.append(item[1])
                            count += 1
                    else:
                        widths.append(zz)
                        count += 1
                # If true, all attributes have been read
                if count == self.size[nb]:
                    break

        with open(fname, "r") as fp:
            # Read header
            for line in fp:
                item = line.split()
                if len(item) > 0:
                    # Searches for line of format *GRID *CORNER I J K
                    if item[0] == "GRID" or item[0] == "*GRID":
                        self.gridType = item[1]
                        self.size = np.array(item[2:5], dtype=int)
                        break
            # Read DI, DJ
            readBlockSizes('I')
            readBlockSizes('J')

            for line in fp:
                item = line.split()
                if item[0] == "ZCORN" or item[0] == "*ZCORN":
                    break
            self.calcCoords(fp)

            # # Read COORD
            # for line in fp:
            #     item = line.split()
            #     if item[0] == '*COORD':
            #         break
            # for line in fp:
            #     item = line.split()
            #     if item[0][:1] == "**":
            #         continue

            # Read NULL
            for line in fp:
                item = line.split()
                # Assumes NULL keyword followed by ALL
                if item[0] == "NULL" or item[0] == "*NULL":
                    break
            self.buildActiveCells(fp)

        # Add in active cells
        ac = vtk.vtkIntArray()
        ac.SetName("ActiveCells")
        for iac in self.ActiveCells.flatten(order='F'):
            ac.InsertNextTuple1(iac)
        self.Grid.GetCellData().AddArray(ac)

    # Builds a cartesian grid from a CMG output file (.out)
    def buildCart(self, fname):
        print('Building cartesian grid')
        self.iWidths = []
        self.jWidths = []
        with open(fname, "r") as fp:

            # Read header
            for line in fp:
                item = line.split()
                if len(item) > 0:
                    # Searches for line of format *GRID *CART I J K
                    if item[0] == "GRID" or item[0] == "*GRID":
                        self.gridType = item[1]
                        self.size = np.array(item[2:5], dtype=int)
                        break

            kSpacing = 0
            depth = 0
            # Assumes DEPTH is the final keyword describing grid structure (move 'break' if not)
            for line in fp:
                item = line.split()
                # Read DI
                if item[0] == "DI" or item[0] == "*DI":
                    if item[1] == "CON" or item[1] == "*CON":
                        self.X = np.arange(0, float(item[2]) * (self.size[0] + 1), float(item[2]))
                # Read DJ
                elif item[0] == "DJ" or item[0] == "*DJ":
                    if item[1] == "CON" or item[1] == "*CON":
                        self.Y = np.arange(0, float(item[2]) * (self.size[1] + 1), float(item[2]))
                # Read DK
                elif item[0] == "DK" or item[0] == "*DK":
                    if item[1] == "CON" or item[1] == "*CON":
                        kSpacing = float(item[2])
                # Read DEPTH (assumes of form *DEPTH *TOP I J K depth)
                elif item[0] == "DEPTH" or item[0] == "*DEPTH":
                    depth = float(item[5])
                    self.Z = np.arange(depth, depth - (kSpacing * self.size[2]), -kSpacing)
                    break

        # Write cell vertex coordinates
        XX, YY, ZZ = ([] for el in range(3))
        for k in range(self.size[2] + 1):
            for j in range(self.size[1] + 1):
                XX.extend(self.X)
                YY.extend([self.Y[j]] * (self.size[0] + 1))
            ZZ.extend([self.Z[k]] * (self.size[0] + 1) * (self.size[1] + 1))

        # Convert to vtk grid
        self.GridType = "vtkStructuredGrid"
        self.Grid = vtk.vtkStructuredGrid()
        self.Grid.SetDimensions(self.size[0] + 1, self.size[1] + 1, self.size[2] + 1)
        vtk_points = vtk.vtkPoints()
        for point in range(len(XX)):
            vtk_points.InsertNextPoint(XX[point], YY[point], ZZ[point])
        self.Grid.SetPoints(vtk_points)

    # Helps buildCorner in constructing cell vertex coordinates
    def calcCoords(self, fp):
        # -Convert CMG DI,DJ cell width data to x,y coordinate lists
        # -Write coordinate lists following ZCORN ordering
        # -There will be many duplicated points at this step, since
        #  grid blocks may share corners

        def writeCorners():
            x1 = 0
            x2 = 0
            for i in range(self.size[0]):
                x2 = x1 + float(self.iWidths[i])
                X.extend([x1, x2])
                Y.extend([y, y])
                x1 = x2

        def writeCoords(coord):
            XX.append(X[coord])
            YY.append(Y[coord])
            ZZ.append(Z[coord])

        def buildLayer():
            for j in range(self.size[1]):
                for i in range(self.size[0]):
                    # write NW corner
                    if i == 0:
                        nwCoord = 2 * i + 4 * self.size[0] * j + const
                        writeCoords(nwCoord)
                    # write NE corner
                    neCoord = 2 * i + 4 * self.size[0] * j + const + 1
                    writeCoords(neCoord)
                if j == self.size[1] - 1:
                    for i in range(self.size[0]):
                        # write SW corner
                        if i == 0:
                            swCoord = 2 * i + 4 * self.size[0] * j + 2 * self.size[0] + const
                            writeCoords(swCoord)
                        # write SE corner
                        seCoord = 2 * i + 4 * self.size[0] * j + 2 * self.size[0] + const + 1
                        writeCoords(seCoord)

        X, Y, Z, XX, YY, ZZ = ([] for i in range(6))
        # Write corner x,y coordinates for top, then bottom of entire i,j cell layer
        for layer in range(2):
            for k in range(self.size[2]):
                y = 0
                for j in range(self.size[1]):
                    # NW-T and NE-T corners
                    writeCorners()
                    y += float(self.jWidths[j])
                    # SW-T and SE-T corners
                    writeCorners()

        # Write z-coordinates from ZCORN
        for line in fp:
            item = line.split()
            if item[0][0] != "*":
                for zz in item:
                    if "*" in zz:
                        item = zz.split("*")
                        for i in range(0, int(item[0])):
                            Z.append(float(item[1]))
                    else:
                        Z.append(float(zz))
            else:
                break

        self.GridType = "vtkStructuredGrid"
        self.Grid = vtk.vtkStructuredGrid()
        self.Grid.SetDimensions(self.size[0] + 1, self.size[1] + 1, self.size[2] + 1)
        const = 0
        for k in range(self.size[2]):
            buildLayer()
            if k == self.size[2] - 1:
                const += self.size[0] * self.size[1] * 4
                buildLayer()
                break
            else:
                const += self.size[0] * self.size[1] * 8

        vtk_points = vtk.vtkPoints()
        for point in range(len(XX)):
            vtk_points.InsertNextPoint(XX[point], YY[point], ZZ[point])
        self.Grid.SetPoints(vtk_points)

    def buildActiveCells(self, fp):
        self.ActiveCells = []
        count = 0
        for line in fp:
            item = line.split()
            for zz in item:
                if "*" in zz:
                    item = zz.split("*")
                    for i in range(0, int(item[0])):
                        self.ActiveCells.append(int(item[1]))
                        count += 1
                else:
                    self.ActiveCells.append(int(zz))
                    count += 1
            if count == self.size[0] * self.size[1] * self.size[2]:
                break
        self.ActiveCells = np.array(self.ActiveCells)
        self.ActiveCells = np.reshape(self.ActiveCells, (self.size[0], self.size[1], self.size[2]), order="F")

    # Reads property attr_name from .dat file
    def readProperty(self, fname, attr_name, add=True):
        typeVal = None
        val = 0
        with open(fname, "r") as fp:
            for line in fp:
                item = line.split()
                if len(item) > 0:
                    if item[0] == attr_name:
                        if len(item) >= 2:
                            if item[1] == "*CON":
                                val = float(item[2])
                                typeVal = '*CON'
                            elif item[1] == '*EQUALSI' or item[1] == 'EQUALSI':
                                attr_I = attr_name[:-1] + 'I'
                                # Change 'PERMJ' to be the keyword that identifies the end of attribute section
                                data = self.readProperty(fname, attr_I, add=False)
                                if len(item) == 4:
                                    op = item[2]
                                    if op == '*':
                                        data *= float(item[3])
                                    elif op == '/':
                                        data /= float(item[3])
                                    elif op == '+':
                                        data += float(item[3])
                                    elif op == '-':
                                        data -= float(item[3])
                            elif item[1] == 'ALL':
                                typeVal = 'ALL'
                        break

            if typeVal == 'ALL':
                data = []
                count = 0
                for line in fp:
                    item = line.split()
                    for attr in item:
                        if "*" in attr:
                            item = attr.split("*")
                            for i in range(0, int(item[0])):
                                data.append(float(item[1]))
                                count += 1
                        else:
                            data.append(float(attr))
                            count += 1
                    # If true, all values have been read
                    if count == self.size[0] * self.size[1] * self.size[2]:
                        data = np.array(data)
                        data = np.reshape(data, (self.size[0], self.size[1], self.size[2]), order="F")
                        break
            elif typeVal == '*CON':
                data = np.full((self.size[0], self.size[1], self.size[2]), val)

        if add:
            self.addToGrid(data, attr_name)
        return data

    # Reads prop from external file denoted by INCLUDE in .dat
    def readExternalProperty(self, fname, attr_title, mult=1):
        data = []
        count = 0
        with open(fname, "r") as fp:
            for line in fp:
                item = line.split()
                for attr in item:
                    if "*" in attr:
                        item = attr.split("*")
                        for i in range(0, int(item[0])):
                            data.append(float(item[1]) * mult)
                            count += 1
                    else:
                        data.append(float(attr) * mult)
                        count += 1
                # If true, all values have been read
                if count == self.size[0] * self.size[1] * self.size[2]:
                    data = np.array(data)
                    data = np.reshape(data, (self.size[0], self.size[1], self.size[2]), order="F")
                    break
        self.addToGrid(data, attr_title)

    # Add data to VTK grid
    def addToGrid(self, data, attr_name):
        ac = vtk.vtkDoubleArray()
        ac.SetName(attr_name)
        for iac in data.flatten(order='F'):
            ac.InsertNextTuple1(iac)
        self.Grid.GetCellData().AddArray(ac)

    # Populates entire K-layer with val (for reading .out property)
    def buildConstLayer(self, val):
        jKeys = np.arange(self.size[1])
        kLayer = dict((el, []) for el in jKeys)
        for j in range(self.size[1]):
            iRow = []
            for i in range(self.size[0]):
                iRow.append(val)
            kLayer[j] = iRow
        return kLayer

    # -This builds a dictionary of a desired attribute for every timestep present in the .out file
    # -attr_name is the desired attribute as it appears in the .out file
    # -If a cell property is empty, then this will set it to null
    # -attr_title is how the attr will appear in the vts file
    def readOutputProperty(self, fname, outputProps):
        # Set up dictionary for attr_name
        # Doing so will allow us to read attr_name values for every time step
        self.outputProps = {}
        for prop in outputProps:
            attr_name = prop[0].replace(" ", "").strip()
            attr_title = prop[1]
            # setattr(self, attr_title, {})
            self.outputProps[attr_title] = {}
            if not hasattr(self, 'times'):
                self.times = []
            layers = {}
            propIdxs = []
            I, J, K = (None,) * 3
            time = None
            build = False
            buildTimes = True

            with open(fname, "r") as fp:
                for line in fp:
                    item = line.split()
                    if len(item) > 0:
                        # Find current time step
                        if item[0] == 'Time':
                            time = item[2]
                        attr = line.replace(" ", "").strip()
                        # Locate attribute name
                        if attr == attr_name:
                            build = True
                            layers = {}
                            propIdxs = []
                            I = None
                            J = None
                            K = '1'

                        if build:
                            if item[0] == 'All':
                                kKeys = np.arange(self.size[2])
                                grid = dict((el, {}) for el in kKeys)
                                for k in range(self.size[2]):
                                    kLayer = self.buildConstLayer(item[3])
                                    grid[k] = kLayer
                                self.outputProps[attr_title][time] = grid
                                build = False
                                continue

                            if item[0] == 'Plane':
                                K = item[3]
                                if len(item) > 4:
                                    if item[4] == 'All':
                                        kLayer = self.buildConstLayer(item[7])
                                        layers[K] = kLayer
                                else:
                                    I = None
                                    layers[K] = {}

                            if item[0] == 'I':
                                if K == '1' and I is None:
                                    layers[K] = {}
                                J = None
                                propIdxs = []
                                I = item[2:]
                                prevDigit = False
                                for i in range(len(line)):
                                    if line[i].isdigit():
                                        if not prevDigit:
                                            propIdxs.append(i)
                                            prevDigit = True
                                    else:
                                        prevDigit = False

                            # Check if there are any missing values in J line
                            skipItem = []
                            if item[0] == 'J=':
                                JIdx = item[1]
                                J = item[2:]
                                if JIdx not in layers[K].keys():
                                    layers[K][JIdx] = []
                                for i in range(len(propIdxs)):
                                    if line[propIdxs[i]] == ' ':
                                        skipItem.append(i)
                                numSkips = 0
                                for i in range(len(I)):
                                    if i in skipItem:
                                        layers[K][JIdx].append('NULL')
                                        numSkips += 1
                                    else:
                                        layers[K][JIdx].append(J[i - numSkips])

                            # Put entire grid worth of property in dictionary for current time step
                            if I and J:
                                if int(I[-1]) == self.size[0] and int(JIdx) == self.size[1] and int(K) == self.size[2]:
                                    if build:
                                        self.outputProps[attr_title][time] = layers
                                        layers = {}
                                        build = False
                # Convert layers from dictionary to arrays
                timeSeriesData = {}
                if len(self.times) > 0:
                    buildTimes = False
                for t in self.outputProps[attr_title].keys():
                    if buildTimes:
                        self.times.append(t)
                    timeSeriesData[t] = []
                    for j in self.outputProps[attr_title][t].keys():
                        for k in self.outputProps[attr_title][t][j].keys():
                            timeSeriesData[t] += self.outputProps[attr_title][t][j][k]
                self.outputProps[attr_title] = timeSeriesData

    # Build and export VTK grids and output prop arrays for every timestep
    def exportGrid(self, fname_vtk, toVTK=True, toNumpy=True):
        print('Exporting grids')
        tID = 0
        for t in self.times:
            propIds = []
            for prop in self.outputProps.keys():
                data = np.array(self.outputProps[prop][t])
                # Save to numpy
                if toNumpy:
                    self.exportProp(data, prop, tID)
                # Save to VTK
                if toVTK:
                    ac = vtk.vtkDoubleArray()
                    ac.SetName(prop)

                    def isfloat(value):
                        try:
                            float(value)
                            return True
                        except ValueError:
                            return False

                    idx = 0
                    for iac in data:
                        # po
                        if not iac[0].isdigit():
                            iac = 0
                        else:
                            if isfloat(iac):
                                iac = float(iac)
                            else:
                                iac = 0
                        idx += 1
                        ac.InsertNextTuple1(iac)
                    id = self.Grid.GetCellData().AddArray(ac)
                    propIds.append(id)

            if tID == 0:
                self.checkOutputDir('vtk')
            self.exportVTK(self.outputDir + '/vtk/' + fname_vtk + str(tID))
            for id in propIds:
                self.Grid.GetCellData().RemoveArray(id)
            tID += 1

    # Given an array of well dictionaries, this creates well data layers and adds to grid
    def addWellLayers(self, wells):
        # Fetch all well constraint types for injectors and producers
        # Create data layers for each type
        all_wells = np.zeros(self.size[0] * self.size[1] * self.size[2])
        well_data_layers = {'ALL': all_wells, 'INJ': {}, 'PRO': {}}
        for well in wells:
            i = 0
            loc = well['LOC']
            idx = ((self.size[0] * self.size[1]) * (loc[2] - 1)) + (self.size[0] * (loc[1] - 1)) + (loc[0] - 1)
            if well['TYPE'] == 'INJ':
                for con_type in well['CON_TYPE']:
                    if con_type not in well_data_layers['INJ'].keys():
                        well_data_layers['INJ'][con_type] = np.zeros(self.size[0] * self.size[1] * self.size[2])
                    else:
                        well_data_layers['INJ'][con_type][idx] = well['CON_VAL'][i]
                    i += 1
                well_data_layers['ALL'][idx] = 1
            else:
                for con_type in well['CON_TYPE']:
                    if con_type not in well_data_layers['PRO'].keys():
                        well_data_layers['PRO'][con_type] = np.zeros(self.size[0] * self.size[1] * self.size[2])
                    i += 1
                well_data_layers['ALL'][idx] = -1

        # Add data layers to grid
        data_layer = well_data_layers['ALL']
        self.addToGrid(data_layer, 'WELLS ALL')
        for con_type in well_data_layers['INJ'].keys():
            data_layer = well_data_layers['INJ'][con_type]
            self.addToGrid(data_layer, 'WELLS INJ ' + con_type)
        for con_type in well_data_layers['PRO'].keys():
            data_layer = well_data_layers['PRO'][con_type]
            self.addToGrid(data_layer, 'WELLS PRO ' + con_type)

    # Contains useful methods for working with wells
    class Wells:
        def __init__(self, wells, times, outputDir):
            self.wells = wells
            self.times = times
            self.outputDir = outputDir

        # Return list of well names
        def names(self):
            return [well['NAME'] for well in self.wells]

        # Get a single well by name
        def well(self, name):
            for well in self.wells:
                if well['NAME'] == name:
                    return well
            raise KeyError('Well name not recognized')

        # Builds a dictionary of wells and their associated properties
        def getWells(self, fname):
            getLoc = False
            wells = []
            well = {'NAME': None, 'TYPE': None, 'OP_MODE': [], 'CON_TYPE': [], 'CON_VAL': [], 'LOC': None}
            with open(fname, "r") as fp:
                for line in fp:
                    item = line.split()
                    if len(item) > 0:
                        keyword = item[0]
                        if getLoc:
                            if item[0][0] == '*':
                                continue
                            well['LOC'] = (int(item[0]), int(item[1]), int(item[2]))
                            getLoc = False
                        elif keyword == 'WELL':
                            # Add the previous well to the grid
                            if well['NAME'] is not None:
                                wells.append(copy.deepcopy(well))
                                well = {'NAME': None, 'TYPE': None, 'OP_MODE': [], 'CON_TYPE': [], 'CON_VAL': [],
                                        'LOC': None}
                            well['NAME'] = item[1].strip("'")
                        elif keyword == 'INJECTOR':
                            well['TYPE'] = 'INJ'
                        elif keyword == 'PRODUCER':
                            well['TYPE'] = 'PRO'
                        elif keyword == 'OPERATE':
                            well['OP_MODE'].append(item[1])
                            well['CON_TYPE'].append(item[2])
                            well['CON_VAL'].append(item[3])
                        elif keyword == 'PERF':
                            getLoc = True
                wells.append(well)
            return self.Wells(wells, self.times, self.outputDir)

        # Helper fn for readOutput
        # In GEMFIELDSUMMARY blocks, not all wells may be listed
        # This finds which ones are and returns an ordered list of their names
        def getOutputOrdering(self, fname):
            t = 1
            orderDict = {}
            order = []
            readWells = False
            lastBlock = False
            addOrder = False
            with open(fname, "r") as fp:
                for line in fp:
                    if readWells:
                        if lastBlock:
                            line = line.split('++')[0]
                            addOrder = True
                            lastBlock = False
                        item = list(map(str.strip, line.split('+')))
                        item = [e.split() for e in list(filter(None, item))]
                        order.extend([w[1] for w in item])
                        readWells = False
                        if addOrder:
                            orderDict[t] = order
                            order = []
                            addOrder = False
                            t += 1
                    elif len(line.split()) > 0:
                        if 'No.' and 'Name' in line:
                            if '++' in line:
                                lastBlock = True
                            readWells = True
                            next(fp)
                            continue
            return orderDict

        # Reads well information for all timesteps from .out
        def readOutput(self, fname, keys, subkeys):
            k = 0
            sk = 0
            tID = 0
            build_keys = False
            build_subkeys = False
            cur_key = None
            cur_block = 1
            n_blocks = 0
            order = self.getOutputOrdering(fname)
            # Initialize output dictionary
            wellOutput = {}
            for i,key in enumerate(keys):
                wellOutput[key] = {}
                for subkey in subkeys[i]:
                    wellOutput[key][subkey] = {}
                    for t in range(len(self.times) - 1):
                        wellOutput[key][subkey][t+1] = []
            with open(fname, "r") as fp:
                for line in fp:
                    item = line.split()
                    if len(item) > 0:
                        # Find current time step
                        if item[0] == 'TIME:':
                            head = ''.join(item[2:])
                            if 'GEMFIELDSUMMARY' in head:
                                build_keys = True
                                tID += 1
                                n_wells = len(order[tID])
                                n_blocks = math.ceil(n_wells / 4)
                                continue
                        # Assume that keywords are ordered as they appear in .out
                        if build_keys:
                            # Current block has been read, move to next one
                            if k == len(keys):
                                cur_block += 1
                                k = 0
                            if keys[k] in line:
                                cur_key = keys[k]
                                build_keys = False
                                build_subkeys = True
                                continue
                        elif build_subkeys:
                            if sk == len(subkeys[k]):
                                build_subkeys = False
                                build_keys = True
                                k += 1
                                sk = 0
                                continue
                            for subkey in subkeys[k]:
                                if subkey in line:
                                    if cur_block == n_blocks:
                                        line = line.split('++')[0]
                                    item = list(map(str.strip, line.split('+')))
                                    item = list(filter(None, item))
                                    wellOutput[cur_key][subkey][tID].extend(item[1:])
                                    sk += 1
            # Attach well names to well outputs
            for key in wellOutput:
                for subkey in wellOutput[key]:
                    for t in range(1, len(self.times)):
                        wellOutput[key][subkey][t] = {k:v for k,v in zip(order[t], wellOutput[key][subkey][t])}
            return wellOutput

    # Exports Numpy array of grid property
    def exportProp(self, prop, title, tID):
        data = np.reshape(prop, (self.size[0], self.size[1], self.size[2]))
        self.checkOutputDir(title)
        np.save(self.outputDir + '/' + title + '/' + title + '_' + str(tID), data)

    # Exports Numpy array of well dictionaries
    def exportWells(self, title, wells):
        self.checkOutputDir(title)
        np.save(self.outputDir + '/' + title, wells)
