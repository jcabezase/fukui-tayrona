#-*- coding: utf-8 -*-

# -------------------------------------------------------------------------------- #
# File        : FUKUI-Tayrona.py                                                   #
# Author      : Javiera Cabezas-Escares, Nicolas F. Barrera, Prof. Carlos Cardenas #
# Company     : http://tayrona.ciencias.uchile.cl                                  #
# Date        : 01/09/2024                                                         #
# -------------------------------------------------------------------------------- #

import numpy as np
from numpy.linalg import norm
from scipy import stats
from scipy import ndimage
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import sys
import math
import os
from scipy.interpolate import interp1d
###########################################################
#              Fukui f^- for Interpolation                #
###########################################################

def Electro_fukui_interpolation(FILE1,FILE2,FILE3,FILE4):
    # ******** Change dn according to your needs***
    dn=np.array([-0.15,-0.1,-0.05,0.0])
    print("\delta N is: ", dn)
    #**********************************************
    #FILE1=input("write CHARGE-FILE-1 using quotation marks: ")
    #FILE2=input("CHARGE-FILE-2: ")
    #FILE3=input("CHARGE-FILE-3: ")
    #FILE4=input("CHARGE-FILE-4: ")

    #Open the file
    #FILE1=open(FILE1, "r")
    
    #FILE1='CHGCAR-015'
    #FILE2='CHGCAR-010'
    #FILE3='CHGCAR-005'
    #FILE4='CHGCAR_00'
    
    #
    # read number of atoms and grid points 
    print ("FILE1: ",FILE1)
    NATOMS = 0
    with open(FILE1) as fp:
        for i, line in enumerate(fp):
            if i == 6:
                line = map(int, line.split())
                NATOMS = sum(line)
                print('NATOMS = ',NATOMS)
            elif i == NATOMS + 9:
                GRID = list(map(int, line.split()))
                print('GRID = ', GRID)
                POINTS = np.prod(GRID)
                print('POINTS= ', POINTS)
            elif i > NATOMS + 11:
                break
    #print('GRID IS  typre : ',type(GRID))     
    #print('POINTS IS  typre : ',type(POINTS))
    #Read heading of CHGCAR FILE
    FUKUIFILE= open("FUKUI","w+")
    with open(FILE1) as fp:
        for i, line in enumerate(fp):
            if i < NATOMS + 10: 
                FUKUIFILE.write(line)
            elif i >= NATOMS + 10:
                FUKUIFILE.close()
                break
    # read data
    skipelines = NATOMS + 9
    maxrows = int(np.ceil(POINTS/5))
    #print ('lines to skipe = ', skipelines)
    #print ('max rows to read = ', maxrows)
    #CHG1 = np.genfromtxt(FILE1,delimiter=" ",skip_header=skipelines,max_rows=maxrows,invalid_raise=True)
    #*****CHARGE file 1*****
    CHG1 = pd.read_fwf(FILE1,colspecs='infer', skiprows=skipelines,nrows=maxrows)

    #print('Changing data typre from: ',type(CHG1)) 
    CHG1=CHG1.to_numpy()
    #print('to: ',type(CHG1))

    #Flatten np.Array
    CHG1=CHG1.flatten('C')
    #Leave NaN out of the np.Array
    CHG1=CHG1[0:POINTS]
    
    # Check data 
    print('CHG1',CHG1)
    #print('lenght of ',FILE1,' : ',len(CHG1)) 
    #print ('first entry in ',FILE1,' : ',CHG1[0])
    #print ('last entry in',FILE1,' : ',CHG1[-1])
    
    #*****CHARGE file 2*****

    CHG2 = pd.read_fwf(FILE2,colspecs='infer', skiprows=skipelines,nrows=maxrows)

    #print('Changing data typre from: ',type(CHG2)) 
    CHG2=CHG2.to_numpy()
    #print('to: ',type(CHG2)) 

    #Flatten np.Array
    CHG2=CHG2.flatten('C')
    #Leave NaN out of the np.Array
    CHG2=CHG2[0:POINTS]

    # Check data 
    print('CHG2',CHG2)
    #print('lenght of ',FILE2,' : ',len(CHG2)) 
    #print ('first entry in ',FILE2,' : ',CHG2[0])
    #print ('last entry in',FILE2,' : ',CHG2[-1])

    #*****CHARGE file 3*****

    CHG3 = pd.read_fwf(FILE3,colspecs='infer', skiprows=skipelines,nrows=maxrows)

    #print('Changing data typre from: ',type(CHG3)) 
    CHG3=CHG3.to_numpy()
    #print('to: ',type(CHG3)) 

    #Flatten np.Array
    CHG3=CHG3.flatten('C')
    #Leave NaN out of the np.Array
    CHG3=CHG3[0:POINTS]

    # Check data 
    print('CHG3',CHG3)
    #print('lenght of ',FILE3,' : ',len(CHG3)) 
    #print ('first entry in ',FILE3,' : ',CHG3[0])
    #print ('last entry in',FILE3,' : ',CHG3[-1])

    #*****CHARGE file 4*****

    CHG4 = pd.read_fwf(FILE4,colspecs='infer', skiprows=skipelines,nrows=maxrows)

    #print('Changing data typre from: ',type(CHG4)) 
    CHG4=CHG4.to_numpy()
    #print('to: ',type(CHG4))

    #Flatten np.Array
    CHG4=CHG4.flatten('C')
    #Leave NaN out of the np.Array
    CHG4=CHG4[0:POINTS]

    # Check data 
    print('CHG4',CHG4)
    #print('lenght of ',FILE4,' : ',len(CHG4)) 
    #print ('first entry in ',FILE4,' : ',CHG4[0])
    #print ('last entry in',FILE4,' : ',CHG4[-1])

    #create array for Fukui
    FUKUI = np.zeros(POINTS)
    RSQUARED = np.zeros(POINTS)
    for i in range(POINTS):
            density= np.array([CHG1[i],CHG2[i],CHG3[i],CHG4[i]]).astype(float)
            #print("Shape densiti [i]", density.shape) 
            slope, intercept, r_value, p_value, std_err = stats.linregress(dn, density)
            FUKUI[i] = slope
            RSQUARED[i] = r_value**2
    
    # Number of full rows and size of last row
    num_full_rows = POINTS // 5
    last_row_size = POINTS % 5
    
    print('Fukui size', FUKUI.size)
    print('Fukui', FUKUI) 
    # Create an appropriately sized array
    if last_row_size == 0:
        FUKUI_matrix = FUKUI.reshape(-1, 5)
    else:
        FUKUI_matrix = np.zeros((num_full_rows + 1, 5))
        FUKUI_matrix[:num_full_rows] = FUKUI[:num_full_rows * 5].reshape(-1, 5)
        FUKUI_matrix[-1, :last_row_size] = FUKUI[num_full_rows * 5:]

    # Save to file with blanks instead of zeros if there is a partial last row
    with open("FUKUI", "a") as FUKUIFILE:
        for row in FUKUI_matrix:
            # If the row has correct values ​​(all data), it is saved normally
            if np.all(row != 0):
                formatted_row = " ".join(f"{value: .12E}" for value in row)
            else:
                # Convert row to a text string with blanks for zero values
                formatted_row = " ".join(f"{value: .12E}" if value != 0 else "" for value in row) 
            FUKUIFILE.write(formatted_row + "\n")
    print('Dimensions of FUKUI written to a FILE =', FUKUI_matrix.shape)
    print("")

############################################################
#            Fukui f^+ for Interpolation                   #
############################################################

def Nucleo_fukui_interpolation(FILE1,FILE2,FILE3,FILE4):
    
    # ******** Change dn according to your needs***
    #dn=np.array([-0.15,-0.1,-0.05,0.0])
    dn=np.array([0.0,0.05,0.1,0.15])
    print("\delta N is: ", dn)
    #**********************************************

    print ("FILE1: ",FILE1)
    NATOMS = 0
    with open(FILE1) as fp:
        for i, line in enumerate(fp):
            if i == 6:
                line = map(int, line.split())
                NATOMS = sum(line)
                print('NATOMS = ',NATOMS)
            elif i == NATOMS + 9:
                GRID = list(map(int, line.split()))
                print('GRID = ', GRID)
                POINTS = np.prod(GRID)
                print('POINTS= ', POINTS)
            elif i > NATOMS + 11:
                break
    
    #print('GRID IS  type : ',type(GRID))
    #print('POINTS IS  type : ',type(POINTS))
    #Read heading of CHGCAR FILE
    FUKUIFILE= open("FUKUI","w+")
    with open(FILE1) as fp:
        for i, line in enumerate(fp):
            if i < NATOMS + 10:
                FUKUIFILE.write(line)
            elif i >= NATOMS + 10:
                FUKUIFILE.close()
                break
    
    # read data
    skipelines = NATOMS + 9
    maxrows = int(np.ceil(POINTS/5))
    #print('maxrows',maxrows)
    #print ('lines to skipe = ', skipelines)
    #print ('max rows to read = ', maxrows)
    ##CHG1 = np.genfromtxt(FILE1,delimiter=" ",skip_header=skipelines,max_rows=maxrows,invalid_raise=True)

    #*****CHARGE file 1*****
    CHG1 = pd.read_fwf(FILE1,colspecs='infer', skiprows=skipelines,nrows=maxrows)

    #print('Changing data typre from: ',type(CHG1)) 
    CHG1=CHG1.to_numpy()
    #print('to: ',type(CHG1)) 

    #Flatten np.Array
    CHG1=CHG1.flatten('C')
    #Leave NaN out of the np.Array
    CHG1=CHG1[0:POINTS]

    # Check data 
    print('CHG1',CHG1)
    #print('lenght of ',FILE1,' : ',len(CHG1)) 
    #print ('first entry in ',FILE1,' : ',CHG1[0])
    #print ('last entry in',FILE1,' : ',CHG1[-1])

    #*****CHARGE file 2*****

    CHG2 = pd.read_fwf(FILE2,colspecs='infer', skiprows=skipelines,nrows=maxrows)

    #print('Changing data typre from: ',type(CHG2)) 
    CHG2=CHG2.to_numpy()
    #print('to: ',type(CHG2)) 

    #Flatten np.Array
    CHG2=CHG2.flatten('C')
    #Leave NaN out of the np.Array
    CHG2=CHG2[0:POINTS]

    # Check data
    print('CHG2',CHG2)
    #print('lenght of ',FILE2,' : ',len(CHG2)) 
    #print ('first entry in ',FILE2,' : ',CHG2[0])
    #print ('last entry in',FILE2,' : ',CHG2[-1])

    #*****CHARGE file 3*****

    CHG3 = pd.read_fwf(FILE3,colspecs='infer', skiprows=skipelines,nrows=maxrows)

    #print('Changing data typre from: ',type(CHG3)) 
    CHG3=CHG3.to_numpy()
    #print('to: ',type(CHG3)) 

    #Flatten np.Array
    CHG3=CHG3.flatten('C')
    #Leave NaN out of the np.Array
    CHG3=CHG3[0:POINTS]

    # Check data
    print('CHG3',CHG3)
    #print('lenght of ',FILE3,' : ',len(CHG3)) 
    #print ('first entry in ',FILE3,' : ',CHG3[0])
    #print ('last entry in',FILE3,' : ',CHG3[-1])

    #*****CHARGE file 4*****

    CHG4 = pd.read_fwf(FILE4,colspecs='infer', skiprows=skipelines,nrows=maxrows)

    #print('Changing data typre from: ',type(CHG4)) 
    CHG4=CHG4.to_numpy()
    #print('to: ',type(CHG4)) 

    #Flatten np.Array
    CHG4=CHG4.flatten('C')
    #Leave NaN out of the np.Array
    CHG4=CHG4[0:POINTS]

    # Check data
    print('CHG4',CHG4)
    #print('lenght of ',FILE4,' : ',len(CHG4)) 
    #print ('first entry in ',FILE4,' : ',CHG4[0])
    #print ('last entry in',FILE4,' : ',CHG4[-1])

    # do the linear regression 

    #dn=np.array([-0.15,-0.1,-0.05,0.0])
    #print dn
    #print("dn shape= ", dn.shape)
    #i = 2000
    #density= np.array([CHG1[i],CHG2[i],CHG3[i],CHG4[i]])
    #print density
    #print density.shape
    
    #create array for Fukui
    FUKUI = np.zeros(POINTS)
    RSQUARED = np.zeros(POINTS)
    for i in range(POINTS):
        density= np.array([CHG1[i],CHG2[i],CHG3[i],CHG4[i]]).astype(float)
        #print("Shape densiti [i]", density.shape) 
        slope, intercept, r_value, p_value, std_err = stats.linregress(dn, density)
        FUKUI[i] = slope
        RSQUARED[i] = r_value**2
        #print("slope: %f    intercept: %f" % (slope, intercept))
        #print ("\rho_Neutral: %f ",CHG4[i])
        #print("R-squared: %f" % r_value**2) 
        #Plots for testing
        #plt.plot(dn, density, 'o', label='original data')
        #plt.plot(dn, intercept + slope*dn, 'r', label='fitted line')
        #plt.legend()
        #plt.show()

    # Number of full rows and size of last row
    num_full_rows = POINTS // 5
    last_row_size = POINTS % 5

    # Create an appropriately sized array
    if last_row_size == 0:
        FUKUI_matrix = FUKUI.reshape(-1, 5)
    else:
        FUKUI_matrix = np.zeros((num_full_rows + 1, 5))
        FUKUI_matrix[:num_full_rows] = FUKUI[:num_full_rows * 5].reshape(-1, 5)
        FUKUI_matrix[-1, :last_row_size] = FUKUI[num_full_rows * 5:]
   
    # Save to file with blanks instead of zeros if there is a partial last row
    with open("FUKUI", "a") as FUKUIFILE:
        for row in FUKUI_matrix:
            # If the row has correct values ​​(all data), it is saved normally
            if np.all(row != 0):
                formatted_row = " ".join(f"{value: .12E}" for value in row)
            else: 
                # Convert row to a text string with blanks for zero values
                formatted_row = " ".join(f"{value: .12E}" if value != 0 else "" for value in row)
            FUKUIFILE.write(formatted_row + "\n")

    print('Dimensions of FUKUI written to a FILE =', FUKUI_matrix.shape)
    print("")
############################################################
#      Fukui potential via electrodes' method              #
############################################################

def fukui_electrodes(FILE0,FILE1,Epsilon):

    def local_maxima_3D(data, order):
            """Detects local maxima in a 3D array

            Parameters
            ---------
            data : 3d ndarray
            order : int
                How many points on each side to use for the comparison
            Returns
            -------
            coords : ndarray
                coordinates of the local maxima
            values : ndarray
                values of the local maxima
            """
            size = 1 + 2 * order
            footprint = np.ones((size, size, size))
            footprint[order, order, order] = 0
            filtered = ndimage.maximum_filter(data, footprint=footprint, mode='wrap')
            mask_local_maxima = data > filtered
            coords = np.asarray(np.where(mask_local_maxima)).T
            values = data[mask_local_maxima]
            return coords, values

    #print ("FILE1: ",FILE1)
    NATOMS = 0
    with open(FILE1) as fp:
        for i, line in enumerate(fp):
            if i == 2:
                alattvec = list(map(float, line.split()))
                #print('a latt vector = ', alattvec)
            elif i == 3:
                blattvec = list(map(float, line.split()))
            elif i == 4:
                clattvec = list(map(float, line.split()))
            elif i == 6:
                line = map(int, line.split())
                NATOMS = sum(line)
            elif i == NATOMS + 9:
                GRID = list(map(int, line.split()))
                POINTS = np.prod(GRID)
            elif i > NATOMS + 11:
                break

    #Set Grid dimension and such

    NGX = GRID[0]
    NGY = GRID[1]
    NGZ = GRID[2]
    #print("type NGX:",type(NGX))
    
    print("This will take a few seconds.")
    print("")

    # Lattice matrix in direct space in a.u (1 bohr = 0.529177210903 A)

    LATTCURA = (1.0/0.529177210903)*np.dstack([alattvec,blattvec,clattvec])
    # LATTCURA has one more dimension that needed and I don't know much Python
    LATTCURA = LATTCURA[0]

    # Compute lattice matrix in reciprocal space ... 2Pi has been  omitted to be consistent with VASP 
    # verbosity

    LATTCURB = np.zeros((3, 3))
    LATTCURB[0] = np.cross(LATTCURA[1],LATTCURA[2])
    LATTCURB[1] = np.cross(LATTCURA[2],LATTCURA[0])
    LATTCURB[2] = np.cross(LATTCURA[0],LATTCURA[1])

    VOL = np.abs(np.linalg.det(np.dstack([alattvec,blattvec,clattvec])))
    omega = np.dot(np.cross(LATTCURA[1],LATTCURA[2]),LATTCURA[0])

    LATTCURB = LATTCURB/omega

    # construc frequency arrays LPCTX, LPCTY, and LPCTZ. Standard form is used. Zero freq is the
    # first element, then it ingreases until max, then it jumps to the more negative and decreses. 
    # This is so to fit FFT libraries and because of hemriticity

    LPCTX =  np.zeros(NGX)
    LPCTY =  np.zeros(NGY)
    LPCTZ =  np.zeros(NGZ)

    LPCTX = [NX if NX < int(NGX/2)+1 else NX-NGX for NX in range(NGX)]
    LPCTY = [NY if NY < int(NGY/2)+1 else NY-NGY for NY in range(NGY)]
    LPCTZ = [NZ if NZ < int(NGZ/2)+1 else NZ-NGZ for NZ in range(NGZ)]

    GSQU =  np.zeros((NGX,NGY,NGZ))
    for N3 in range(NGZ):
        for N2 in range(NGY):
            for N1 in range(NGX):
                GX = LPCTX[N1]*LATTCURB[0,0] + LPCTY[N2]*LATTCURB[0, 1] + LPCTZ[N3]*LATTCURB[0, 2]
                GY = LPCTX[N1]*LATTCURB[1,0] + LPCTY[N2]*LATTCURB[1, 1] + LPCTZ[N3]*LATTCURB[1, 2]
                GZ = LPCTX[N1]*LATTCURB[2,0] + LPCTY[N2]*LATTCURB[2, 1] + LPCTZ[N3]*LATTCURB[2, 2]
                GSQU[N1,N2,N3] = (GX*GX + GY*GY + GZ*GZ)

    # make zero the charge to avoid divergency of the FT of 1/r = 4pi/G^2
    q = 0.0
    GSQU[0,0,0]=1.0
    GSQU = 1.0/(GSQU*(2.0*np.pi*2.0*np.pi))
    GSQU[0,0,0] = q

    #Read heading of CHGCAR FILE
    skipelines = NATOMS + 9
    if POINTS % 5 == 0:
        maxrows = int(POINTS/5)
        #print('maxrows',maxrows)
    elif POINTS % 5 != 0:
        maxrows = int(np.ceil(POINTS/5))
        #print('maxrows',maxrows)

    skipelines = skipelines + 1
    
    #*****read data file *****
    
    if (NGX * NGY * NGZ) % 5 != 0:
         with open(FILE0, 'r') as file:
            last_line = file.readlines()[maxrows + skipelines-1]
            #print(last_line)
            num_columns_last_line = len(last_line.split())
            #print('The number of columns in the row = ', maxrows, ' es=', num_columns_last_line)
            df0 = pd.read_table(FILE0,sep='\s+',skiprows=skipelines,names=range(5),nrows=maxrows)
            #print('df0:',df0)
    else:
        with open(FILE0, 'r') as file:
            df0 = pd.read_table(FILE0,sep='\s+',skiprows=skipelines,names=range(5),nrows=maxrows)
            #print('df0:',df0)

    if (NGX * NGY * NGZ) % 5 != 0:
        with open(FILE1, 'r') as file:
            last_line = file.readlines()[maxrows + skipelines-1]
            #print(last_line)
            num_columns_last_line = len(last_line.split())
            #print('The number of columns in the row =', maxrows, ' es=', num_columns_last_line)
            df = pd.read_table(FILE1,sep='\s+',skiprows=skipelines,names=range(5),nrows=maxrows)
            #print('df:',df)
    else:
        print("No review necessary")
        with open(FILE1, 'r') as file:
            df = pd.read_table(FILE1,sep='\s+',skiprows=skipelines,names=range(5),nrows=maxrows)
            #print('df:',df)

    #Find number of missing values, displayed as nan
    #remove two for "missing values" in dimension row

    def missing(data):
        nan_count=data.isna().sum().sum()
        #print('number of NAN counts:',nan_count)
        #Convert to 1D numpy array, removing excess nans/dimensions
        if nan_count==0:
            #CHG=data.to_numpy().flatten('C')
            CHG=data.to_numpy()
            CHG=CHG.flatten('C')
            CHG=CHG[0:POINTS]

            return CHG

        else:
            CHG=data.to_numpy()
            CHG=CHG.flatten('C')
            CHG=CHG[0:POINTS]
            
            hay_nan = any(np.isnan(valor) for valor in CHG)

            if hay_nan:
                print("Input error")
            else:
                print("")

            return CHG

    #print('revision df:')
    CHG = missing(df)
    #print('revision df0:')
    CHG00 = missing(df0)

    #print('CHG00',CHG00)
    #print('CHG', CHG)

    # WARNING: in VASP the electrosatic potential is sensed with negative probe. I prefer to use the
    # standard definition in quantum chemistry: the probe is positive

    # reshape array charge array
    def reshape_xyz(NGX,NGY,NGZ, CHG):
        GHGtem =  np.zeros((NGX,NGY,NGZ))
        for N3 in range(NGX):
            for N2 in range(NGY):
                for N1 in range(NGZ):
                    GHGtem[N3,N2,N1] = CHG[N1*NGX*NGY + N2*NGX + N3]
        return GHGtem

    CHGtem = reshape_xyz(NGX, NGY, NGZ, CHG)
    CHGtem00 = reshape_xyz(NGX, NGY, NGZ, CHG00)

    CHG = CHGtem/omega
    CHG00 = CHGtem00/omega

    #### prom_planar CHGCAR ###
    PLANAR_CHG00 = []
    for nz in range(NGZ):
        PLANAR_CHG00.append(np.sum(CHG00[:,:,nz])/(NGY*NGX))

    def valor_position(lista, numero, porcentaje=0.4):
        # Compute the elements amount corresponding to the percentage
        ultimos_elementos = int(len(lista) * porcentaje)
        posicion, valor_cercano = min(enumerate(lista[-ultimos_elementos:]), key=lambda x: abs(x[1] - numero))
        return valor_cercano, posicion + len(lista) - ultimos_elementos

    # Set number
    cota = 0.0001
    cota_ls = 0.02
    
    # Find the closest value and its position
    valor_cercano, pos0 = valor_position(PLANAR_CHG00, cota)
    val_ls, pos_ls = valor_position(PLANAR_CHG00, cota_ls)

    #print("El valor de ls es", pos_ls)

    #plt.plot(PLANAR_CHG00)
    #plt.scatter(pos0, valor_cercano, color='red', marker='o')
    #plt.scatter(pos_ls, val_ls, color='blue',marker='o')
    #plt.show()
    
    # compute the electrostatic potentential in reciprocal space and then FT back to real space
    CHGG = np.fft.fftn(CHG, norm='ortho')
    #print('shape FT(CHG) = ',CHGG.shape)
    #print('shape GSQU = ',GSQU.shape)
    #print("CHGG[1,1,1]=",CHGG[1,1,1])
    LOCPOTG = 4*np.pi*np.multiply(CHGG,GSQU)
    #LOCPOTG = CHGG 
    #print("LOCPOTG[1,1,1]",LOCPOTG[1,1,1])
    #print('shape LOCPOTG = ',LOCPOTG.shape)
    #print('type LOCPOTG = ',type(LOCPOTG))
    LOCPOT = np.fft.ifftn(LOCPOTG,norm='ortho')
    #print('shape IFT(LOCPOTG) = ',LOCPOT.shape)
    #print('type IFT(LOCPOTG) = ',type(LOCPOT))
    #print(LOCPOT[1,1,1])
    imagsum=np.sum(LOCPOT.imag)
    realsum=np.sum(LOCPOT.real)

    LOCPOT = 27.2114*LOCPOT.real ##### -27.2114 convenio quim

    ### CORRECTION

    ##Atoms positions

    def vector_atoms_position(nombre_archivo, N_ATOMS):
        fin = 8 + N_ATOMS
        with open(nombre_archivo, 'r') as archivo:
            lineas = archivo.readlines()
            lineas_interesantes = lineas[8 : fin]
            vectores = [list(map(float, linea.strip().split())) for linea in lineas_interesantes]
            return vectores
    
    atom_position= vector_atoms_position(FILE1, NATOMS)
    
    atom_positionx = [i[0] for i in atom_position]
    atom_positiony = [i[1] for i in atom_position]
    atom_positionz = [i[2] for i in atom_position]

    ##Lengths LZ and LS
    LZ = clattvec[2] #length of Z-axis. 

    max_atomz = max(atom_positionz)
    min_atomz = min(atom_positionz)

    LS = ((LZ/NGZ)*pos_ls - LZ/2)*2

    ######## Correction function ####### as function 
    def func_correction(omega_vol, eps, LS, LZ, NGZ):
            q = 1.6022 *1e-19
            eps_0 = 8.854*1e-22 
            modulo = (-q/(2.0*eps_0*omega_vol)) 
            correction = []
            recorrido = np.linspace(-(LZ*0.5),(LZ*0.5) , NGZ)
            for z in recorrido:
                if (-LS*0.5) < z < (LS*0.5):
                    func_dentro = (1/eps)*(z*z + LS*LS*0.25*(eps -1) - eps*LZ*LZ*0.25)
                    correction.append(func_dentro*modulo)
                elif (abs(z)>= (LS*0.5)):
                    func_fuera = z*z - (LZ*LZ*0.25)
                    correction.append(func_fuera*modulo)
                else:
                    print("Error func_correction", i)
            return correction
    
    salida = [item for sublist in func_correction(VOL, Epsilon, LS, LZ, NGZ) for item in np.array(sublist).flatten() ]
    
    LOCPOT_cor = LOCPOT.copy()

    for i in range(NGZ):
        LOCPOT_cor[:,:,i] = LOCPOT_cor[:,:,i] + salida[i] 
        #print (100,'correccion=', salida[100],'locpot original',LOCPOT[:,:,100],'salida=',  LOCPOT_cor[:,:,100], )

    #print("*****PLANAR****")
    PLANAR1 = []
    PLANAR2 = []
    for nz in range(NGZ):
        PLANAR1.append(np.sum(LOCPOT[:,:,nz])/(NGY*NGX))
        PLANAR2.append(np.sum(LOCPOT_cor[:,:,nz])/(NGY*NGX))
        #print (nz, PLANAR2[nz])

    #nz_end = len(PLANAR2)
    #print ('el largo de la lista es',nz_end)
    #rango_cero_inf = nz_end - pos0 -1
    #cota0 = PLANAR2[rango_cero_inf] 
    #print ('seran cero los valores antes de la posicion', rango_cero_inf, 'con un valor', cota0)

    #cota1 = PLANAR1[pos0]
    cota2 = PLANAR2[pos0]
    #print("La densidad es cercana a ", cota ,"en x=", pos0)
    #PLANAR1f = [elemento - cota1 for elemento in PLANAR1]
    #PLANAR2f = [elemento - cota2 for elemento in PLANAR2] #aqui debe haber un par de ifs 
    PLANAR2f_2 = [elemento - cota2 if 0 <= (elemento - cota2) else 0 for elemento in PLANAR2] #new correction

    correccion_final = [a - b for a, b in zip(PLANAR2f_2, PLANAR2)]
    #print(correccion_final)

    #print('LOCPOT_cor.copy()', LOCPOT_cor.copy())
    LOCPOT_cor2 = LOCPOT_cor.copy()

    for i in range(NGZ):
        LOCPOT_cor2[:,:,i] = LOCPOT_cor2[:,:,i] + correccion_final[i]

    PLANAR3= []
    for nz in range(NGZ):
        PLANAR3.append(np.sum(LOCPOT_cor2[:,:,nz])/(NGY*NGX))
      
    PLANAR1 = [-x for x in PLANAR1]
    PLANAR3 = [-x for x in PLANAR3]

    z_axis = np.linspace(0, LZ, NGZ)
    plt.plot(z_axis,PLANAR3,label='Corrected',linewidth=2)
    plt.plot(z_axis,PLANAR1,label='No correction',linewidth=2)
    #plt.plot(PLANAR2, label = 'corregido')
    #plt.axvline(pos0,color='red', linestyle='--', linewidth=1)
    plt.title('Planar Average Fukui Potential',fontsize=18)
    plt.ylabel(r'$v_{f}(r)$ (eV)',fontsize=12.5)
    plt.xlabel(r'Z-direction ($\AA$)',fontsize=12.5)
    plt.tick_params(axis='both', which='major', labelsize=10)
    plt.legend()
    plt.savefig('PA_vFuk.svg')

    # Supongamos que LOCPOT_cor2 es tu array original
    #LOCPOT = LOCPOT_cor2.flatten()

    # Número total de elementos
    #POINTS = LOCPOT.size
    #print('LOCPOT size',LOCPOT.size)

    # Número de filas completas y tamaño de la última fila
    num_full_rows = POINTS // 5
    last_row_size = POINTS % 5

    FILE= open("FUKUIPOT.LOCPOT","w+")
    with open(FILE1) as fp:
        for i, line in enumerate(fp):
            if i < NATOMS + 10:
                FILE.write(line)
            elif i >= NATOMS + 10:
                FILE.close()
                break
    ### write de locpot in a suitiable form to be written in vasp format (x is de fast index)
    LOCPOTtem =  np.zeros(NGX*NGY*NGZ)
    for N1 in range(NGZ):
        for N2 in range(NGY):
            for N3 in range(NGX):
                LOCPOTtem[N1*NGX*NGY + N2*NGX + N3] = LOCPOT_cor2[N3,N2,N1]

    print("*****LOCPOTtem****")
    print(LOCPOTtem)
    num_elements = LOCPOT_cor2.size
    print(f"El número de elementos en LOCPOT_cor2 es: {num_elements}")

    # Crear una matriz de tamaño adecuado
    if last_row_size == 0:
        LOCPOT = LOCPOTtem*(-1) 
        FUKUIPOT = np.resize(LOCPOT, (maxrows, 5))
        print ('Dimensions of FUKUI POTENTIAL  written to a file = ',FUKUIPOT.shape)
        with open("FUKUIPOT.LOCPOT", "a") as FILE:
            for row in FUKUIPOT:
                formatted_row = " " + " ".join(f"{value: .12E}" for value in row)
                FILE.write(formatted_row + "\n")
        '''        
        FILE= open("FUKUIPOT.LOCPOT","a")
        np.savetxt(FILE, FUKUIPOT, fmt='%.12E', delimiter=" ")
        FILE.close()
        '''
    else:
        LOCPOT = LOCPOTtem*(-1)
        print('LOCPOT size',LOCPOT.size)
        LOCPOT_matrix = np.zeros((num_full_rows + 1, 5))
        LOCPOT_matrix[:num_full_rows] = LOCPOT[:num_full_rows * 5].reshape(-1, 5)
        LOCPOT_matrix[-1, :last_row_size] = LOCPOT[num_full_rows * 5:]      
        print("VER ESTO")
        print('LOCPOT_matrix', LOCPOT_matrix)
        with open("FUKUIPOT.LOCPOT", "a") as FILE:
            for row in LOCPOT_matrix:
                # If the row has correct values <200b><200b>(all data), it is saved normally
                if np.all(row != 0):
                    formatted_row = " ".join(f"{value: .12E}" for value in row)
                else:
                    # Convert row to a text string with blanks for zero values
                    formatted_row = " ".join(f"{value: .12E}" if value != 0 else "" for value in row)
                FILE.write(formatted_row + "\n")
        print('Dimensions of FUKUI written to a FILE =', LOCPOT_matrix.shape)
        print("")
        '''
        with open("FUKUIPOT.LOCPOT", "a") as FILE:
            for row in LOCPOT_matrix:
                formatted_row = " " + " ".join(f"{value: .12E}" for value in row)
                #formatted_row = " ".join(f"{value: .12E}" if value != 0 else " " * 12 for value in row)
                FILE.write(formatted_row + "\n")
        '''

    '''
    # Guardar en el archivo con espacios en blanco en lugar de ceros si hay una última fila parcial
    with open("FUKUIPOT.LOCPOT", "w") as FILE:  
        for row in LOCPOT_matrix:
            formatted_row = " ".join(f"{value: .12E}" if value != 0 else " " * 12 for value in row)
            FILE.write(formatted_row + "\n")

    print('Dimensions of LOCPOT written to a FILE =', LOCPOT_matrix.shape)
    '''

    '''
    ### WRITE FUKUI POTENTIAL CHGCAR FILE 

    FILE= open("FUKUIPOT.LOCPOT","w+")
    with open(FILE1) as fp:
        for i, line in enumerate(fp):
            if i < NATOMS + 10:
                FILE.write(line)
            elif i >= NATOMS + 10:
                FILE.close()
                break
    
    LOCPOT = LOCPOT_cor2.flatten()
    #print('LOCPOT', LOCPOT)
    # Number of completed rows and size of last rows
    num_full_rows = POINTS // 5
    last_row_size = POINTS % 5

    # Create an appropiate matrix size
    if last_row_size == 0:
        LOCPOT_matrix = LOCPOT.reshape(-1, 5)
    else:
        LOCPOT_matrix = np.zeros((num_full_rows + 1, 5))
        LOCPOT_matrix[:num_full_rows] = LOCPOT[:num_full_rows * 5].reshape(-1, 5)
        LOCPOT_matrix[-1, :last_row_size] = LOCPOT[num_full_rows * 5:]

    # Guardar en el archivo con espacios en blanco en lugar de ceros si hay una última fila parcial
    with open("FUKUIPOT.LOCPOT", "a") as FILE:
        for row in LOCPOT_matrix:
            # If the row has a corrected values, they save normally
            if np.all(row != 0):
                formatted_row = " ".join(f"{value: .12E}" for value in row)
            else:
                # Convert the row a text chain with white spaces to zero values
                formatted_row = " ".join(f"{value: .12E}" if value != 0 else "" for value in row)
            FILE.write(formatted_row + "\n")
            #CHGCARSUMFILE.close()
    print('Dimensions of LOCPOT written to a FILE =', LOCPOT_matrix.shape)
    
    '''
############################################################
#             Fukui potential via SCPC method              #
############################################################

def vfukui_SCPC(FILE0,FILE1,FILE2,c):
    #print ("FILE1: ",FILE1)
    NATOMS = 0
    with open(FILE1) as fp:
        for i, line in enumerate(fp):
            if i == 2:
                alattvec = list(map(float, line.split()))
                #print('a latt vector = ', alattvec)
            elif i == 3:
                blattvec = list(map(float, line.split()))
            elif i == 4:
                clattvec = list(map(float, line.split()))
            elif i == 6:
                line = map(int, line.split())
                NATOMS = sum(line)
            elif i == NATOMS + 9:
                GRID = list(map(int, line.split()))
                POINTS = np.prod(GRID)
            elif i > NATOMS + 11:
                break

    #Set Grid dimension and such
    NGX = GRID[0]
    NGY = GRID[1]
    NGZ = GRID[2]
    #print("type NGX:",type(NGX))
    
    print("This will take a few seconds.")
    print("")


    # Lattice matrix in direct space in a.u (1 bohr = 0.529177210903 A)

    LATTCURA = (1.0/0.529177210903)*np.dstack([alattvec,blattvec,clattvec])
    LATTCURA = LATTCURA[0]

    # verbosity
    LATTCURB = np.zeros((3, 3))
    LATTCURB[0] = np.cross(LATTCURA[1],LATTCURA[2])
    LATTCURB[1] = np.cross(LATTCURA[2],LATTCURA[0])
    LATTCURB[2] = np.cross(LATTCURA[0],LATTCURA[1])

    VOL = np.abs(np.linalg.det(np.dstack([alattvec,blattvec,clattvec])))
    omega = np.dot(np.cross(LATTCURA[1],LATTCURA[2]),LATTCURA[0])

    LATTCURB = LATTCURB/omega

    # construc frequency arrays LPCTX, LPCTY, and LPCTZ. Standard form is used. Zero freq is the
    # first element, then it ingreases until max, then it jumps to the more negative and decreses.
    # This is so to fit FFT libraries and because of hermiticity

    LPCTX =  np.zeros(NGX)
    LPCTY =  np.zeros(NGY)
    LPCTZ =  np.zeros(NGZ)

    LPCTX = [NX if NX < int(NGX/2)+1 else NX-NGX for NX in range(NGX)]
    LPCTY = [NY if NY < int(NGY/2)+1 else NY-NGY for NY in range(NGY)]
    LPCTZ = [NZ if NZ < int(NGZ/2)+1 else NZ-NGZ for NZ in range(NGZ)]

    GSQU =  np.zeros((NGX,NGY,NGZ))
    for N3 in range(NGZ):
        for N2 in range(NGY):
            for N1 in range(NGX):
                GX = LPCTX[N1]*LATTCURB[0,0] + LPCTY[N2]*LATTCURB[0, 1] + LPCTZ[N3]*LATTCURB[0, 2]
                GY = LPCTX[N1]*LATTCURB[1,0] + LPCTY[N2]*LATTCURB[1, 1] + LPCTZ[N3]*LATTCURB[1, 2]
                GZ = LPCTX[N1]*LATTCURB[2,0] + LPCTY[N2]*LATTCURB[2, 1] + LPCTZ[N3]*LATTCURB[2, 2]
                GSQU[N1,N2,N3] = (GX*GX + GY*GY + GZ*GZ)
   
    # make zero the charge to avoid divergency of the FT of 1/r = 4pi/G^2
    q = 0.0
    GSQU[0,0,0]=1.0
    GSQU = 1.0/(GSQU*(2.0*np.pi*2.0*np.pi))
    GSQU[0,0,0] = q
    
    #Read heading of CHGCAR FILE
    skipelines = NATOMS + 9
    if POINTS % 5 == 0:
        maxrows = int(POINTS/5)
        #print('maxrows',maxrows)
    elif POINTS % 5 != 0:
        maxrows = int(np.ceil(POINTS/5))
        #print('maxrows',maxrows)

    skipelines = skipelines + 1

    #*****read data file *****
    if (NGX * NGY * NGZ) % 5 != 0:
        with open(FILE0, 'r') as file:
            last_line = file.readlines()[maxrows + skipelines-1]
            #print(last_line)
            num_columns_last_line = len(last_line.split())
            #print('The number of columns in the row = ', maxrows, ' es=', num_columns_last_line)
            df0 = pd.read_table(FILE0,sep='\s+',skiprows=skipelines,names=range(5),nrows=maxrows)
            #print('df0:',df0)

    else:
        with open(FILE0, 'r') as file:
            df0 = pd.read_table(FILE0,sep='\s+',skiprows=skipelines,names=range(5),nrows=maxrows)
            #print('df0:',df0)
    
    if (NGX * NGY * NGZ) % 5 != 0:
        with open(FILE1, 'r') as file:
            last_line = file.readlines()[maxrows + skipelines-1]
            #print(last_line)
            num_columns_last_line = len(last_line.split())
            #print('The number of columns in the row =', maxrows, ' es=', num_columns_last_line)
            df = pd.read_table(FILE1,sep='\s+',skiprows=skipelines,names=range(5),nrows=maxrows)
            #print('df:',df)
    else:
        print("No review necessary")
        with open(FILE1, 'r') as file:
            df = pd.read_table(FILE1,sep='\s+',skiprows=skipelines,names=range(5),nrows=maxrows)
            #print('df:',df)

    #Find number of missing values, displayed as nan
    #remove two for "missing values" in dimension row

    def missing(data):
        nan_count=data.isna().sum().sum()
        #print('number of NAN counts:',nan_count)
        #Convert to 1D numpy array, removing excess nans/dimensions
        if nan_count==0:
            #CHG=data.to_numpy().flatten('C')
            CHG=data.to_numpy()
            CHG=CHG.flatten('C')
            CHG=CHG[0:POINTS]

            return CHG
        else:
            CHG=data.to_numpy()
            CHG=CHG.flatten('C')
            CHG=CHG[0:POINTS]
            hay_nan = any(np.isnan(valor) for valor in CHG)

            if hay_nan:
                print("Input error")
            else:
                print("")

            return CHG

    #print('revision df:')
    CHG = missing(df)
    #print('revision df0:')
    CHG00 = missing(df0)
    #print('CHG00',CHG00)
    #print('CHG', CHG)
    # WARNING: in VASP the electrosatic potential is sensed with negative probe. I prefer to use the
    # standard definition in quantum chemistry: the probe is positive
       
    # reshape array charge array

    def reshape_xyz(NGX,NGY,NGZ, CHG):
        GHGtem =  np.zeros((NGX,NGY,NGZ))
        for N3 in range(NGX):
            for N2 in range(NGY):
                for N1 in range(NGZ):
                    GHGtem[N3,N2,N1] = CHG[N1*NGX*NGY + N2*NGX + N3]
        return GHGtem

    CHGtem = reshape_xyz(NGX, NGY, NGZ, CHG)
    CHGtem00 = reshape_xyz(NGX, NGY, NGZ, CHG00)

    CHG = CHGtem/omega
    CHG00 = CHGtem00/omega

    #### prom_planar CHGCAR ###

    PLANAR_CHG00 = []
    for nz in range(NGZ):
        PLANAR_CHG00.append(np.sum(CHG00[:,:,nz])/(NGY*NGX))

    def valor_position(lista, numero, porcentaje=0.4):
        # Compute the elements amount corresponding to the percentage
        ultimos_elementos = int(len(lista) * porcentaje)
        posicion, valor_cercano = min(enumerate(lista[-ultimos_elementos:]), key=lambda x: abs(x[1] - numero))
        return valor_cercano, posicion + len(lista) - ultimos_elementos
    
    # Set number
    cota = 0.0001
    cota_ls = 0.02

    # Find the closest value and its position
    valor_cercano, pos0 = valor_position(PLANAR_CHG00, cota)
    val_ls, pos_ls = valor_position(PLANAR_CHG00, cota_ls)

    #print("El valor de ls es", pos_ls)

    #plt.plot(PLANAR_CHG00)
    #plt.scatter(pos0, valor_cercano, color='red', marker='o')
    #plt.scatter(pos_ls, val_ls, color='blue',marker='o')
    #plt.show()

    # compute the electrostatic potentential in reciprocal space and then FT back to real space
    CHGG = np.fft.fftn(CHG, norm='ortho')
    #print('shape FT(CHG) = ',CHGG.shape)
    #print('shape GSQU = ',GSQU.shape)
    #print("CHGG[1,1,1]=",CHGG[1,1,1])
    LOCPOTG = 4*np.pi*np.multiply(CHGG,GSQU)
    #LOCPOTG = CHGG
    #print("LOCPOTG[1,1,1]",LOCPOTG[1,1,1])
    #print('shape LOCPOTG = ',LOCPOTG.shape)
    #print('type LOCPOTG = ',type(LOCPOTG))
    LOCPOT = np.fft.ifftn(LOCPOTG,norm='ortho')
    #print('shape IFT(LOCPOTG) = ',LOCPOT.shape)
    #print('type IFT(LOCPOTG) = ',type(LOCPOT))
    #print(LOCPOT[1,1,1])
    imagsum=np.sum(LOCPOT.imag)
    realsum=np.sum(LOCPOT.real)
    
    LOCPOT = 27.2114*LOCPOT.real 
    
    ### CORRECTION
    ##Lengths LZ and LS
    LZ = clattvec[2] #length of Z-axis.
    z_axis = np.linspace(0, LZ, NGZ)

    PLANAR1 = []
    for nz in range(NGZ):
        PLANAR1.append(np.sum(LOCPOT[:,:,nz])/(NGY*NGX))
    
    value1 = np.array(PLANAR1)

    # Read SCPC correction file (FILE2)
    with open(FILE2, 'r') as f:
        lines = f.readlines()

    # Find the last occurrence of the marker '#z-vcor.dat', ignoring leading spaces
    data_start = None
    for i in range(len(lines)-1, -1, -1):
        if lines[i].strip().startswith('#z-vcor.dat'):
            data_start = i + 1
            break

    if data_start is None:
        raise ValueError("El archivo FILE2 no contiene el marcador '#z-vcor.dat'")

    data_lines = lines[data_start:]
    
    data_lines = lines[data_start:]
    z2, value2 = [], []
    for line in data_lines:
        if line.strip():
            cols = line.split()
            z2.append(float(cols[0]))
            value2.append(float(cols[1]))
    z2 = np.array(z2)
   
    '''
    # Read SCPC correction file (FILE2)
    with open(FILE2, 'r') as f:
        lines = f.readlines()
        data_start = None
        for i in range(len(lines)-1, -1, -1):
            if lines[i].startswith('#z-vcor.dat'):
                data_start = i + 1
                break
    if data_start is None:
        raise ValueError("El archivo FILE2 no contiene el marcador '#z-vcor.dat'")

    # Leer los datos relevantes en FILE2
    data_lines = lines[data_start:]
    z2, value2 = [], []
    for line in data_lines: 
        if line.strip():
            cols = line.split()
            z2.append(float(cols[0]))
            value2.append(float(cols[1]))
    z2 = np.array(z2)
    value2 = np.array(value2)
    '''
    # Interpolar los valores del archivo SCPC utilizando los puntos del eje z calculados del promedio planar
    interp2 = interp1d(z2, value2, kind='cubic', fill_value="extrapolate")
    value2_interp = interp2(z_axis)
    
    print('value2_interp=',value2_interp)

    LOCPOT_cor = LOCPOT.copy()

    for i in range(NGZ):
        LOCPOT_cor[:,:,i] = LOCPOT_cor[:,:,i] - value2_interp[i]*(c)

    PLANAR2 = []
    for nz in range(NGZ):
        PLANAR2.append(np.sum(LOCPOT_cor[:,:,nz])/(NGY*NGX))

    cota2 = PLANAR2[pos0]

    PLANAR2f_2 = [elemento - cota2 if 0 <= (elemento - cota2) else 0 for elemento in PLANAR2]

    correccion_final = [a - b for a, b in zip(PLANAR2f_2, PLANAR2)]

    #print('LOCPOT_cor.copy()', LOCPOT_cor.copy())
    LOCPOT_cor2 = LOCPOT_cor.copy()

    for i in range(NGZ):
                LOCPOT_cor2[:,:,i] = LOCPOT_cor2[:,:,i] + correccion_final[i]

    PLANAR3= []
    for nz in range(NGZ):
        PLANAR3.append(np.sum(LOCPOT_cor2[:,:,nz])/(NGY*NGX))

    PLANAR1 = [-x for x in PLANAR1]
    PLANAR3 = [-x for x in PLANAR3]
    value2_interp = [-x for x in value2_interp]

    z_axis = np.linspace(0, LZ, NGZ)
    plt.plot(z_axis,PLANAR3,label='Corrected',linewidth=2)
    plt.plot(z_axis,PLANAR1,label='No correction',linewidth=2)
    #plt.plot(z_axis,value2_interp,label='Correction',linewidth=2)
    #plt.plot(PLANAR2, label = 'corregido')
    #plt.axvline(pos0,color='red', linestyle='--', linewidth=1)
    plt.title('Planar Average Fukui Potential',fontsize=18)
    plt.ylabel(r'$v_{f}(r)$ (eV)',fontsize=12.5)
    plt.xlabel(r'Z-direction ($\AA$)',fontsize=12.5)
    plt.tick_params(axis='both', which='major', labelsize=10)
    plt.legend()
    plt.savefig('PA_vFuk.svg')
    
    # Supongamos que LOCPOT_cor2 es tu array original
    #LOCPOT = LOCPOT_cor2.flatten()

    # Número total de elementos
    #POINTS = LOCPOT.size
    #print('LOCPOT size',LOCPOT.size)

    # Número de filas completas y tamaño de la última fila
    num_full_rows = POINTS // 5
    last_row_size = POINTS % 5

    FILE= open("FUKUIPOT.LOCPOT","w+")
    with open(FILE1) as fp:
        for i, line in enumerate(fp):
            if i < NATOMS + 10:
                FILE.write(line)
            elif i >= NATOMS + 10:
                FILE.close()
                break
    
    ### write de locpot in a suitiable form to be written in vasp format (x is de fast index)
    LOCPOTtem =  np.zeros(NGX*NGY*NGZ)
    for N1 in range(NGZ):
        for N2 in range(NGY):
            for N3 in range(NGX):
                LOCPOTtem[N1*NGX*NGY + N2*NGX + N3] = LOCPOT_cor2[N3,N2,N1]

    print("*****LOCPOTtem****")
    print(LOCPOTtem)
    num_elements = LOCPOT_cor2.size
    print(f"El número de elementos en LOCPOT_cor2 es: {num_elements}")

    # Crear una matriz de tamaño adecuado
    if last_row_size == 0:
        LOCPOT = LOCPOTtem*(-1)
        FUKUIPOT = np.resize(LOCPOT, (maxrows, 5))
        print ('Dimensions of FUKUI POTENTIAL  written to a file = ',FUKUIPOT.shape)
        with open("FUKUIPOT.LOCPOT", "a") as FILE:
            for row in FUKUIPOT:
                formatted_row = " " + " ".join(f"{value: .12E}" for value in row)
                FILE.write(formatted_row + "\n")

    else:
        LOCPOT = LOCPOTtem*(-1)
        print('LOCPOT size',LOCPOT.size)
        LOCPOT_matrix = np.zeros((num_full_rows + 1, 5))
        LOCPOT_matrix[:num_full_rows] = LOCPOT[:num_full_rows * 5].reshape(-1, 5)
        LOCPOT_matrix[-1, :last_row_size] = LOCPOT[num_full_rows * 5:]
        print("VER ESTO")
        print('LOCPOT_matrix', LOCPOT_matrix)
        with open("FUKUIPOT.LOCPOT", "a") as FILE:
            for row in LOCPOT_matrix:
                # If the row has correct values <200b><200b>(all data), it is saved normally
                if np.all(row != 0):
                    formatted_row = " ".join(f"{value: .12E}" for value in row)
                else:
                    # Convert row to a text string with blanks for zero values
                    formatted_row = " ".join(f"{value: .12E}" if value != 0 else "" for value in row)
                
                FILE.write(formatted_row + "\n")

        print('Dimensions of FUKUI written to a FILE =', LOCPOT_matrix.shape)
        print("")


############################################################
#                    Add                                   #
############################################################

def Add(FILE1,FILE2,c1,c2,c3):
    # read number of atoms and grid points
    print ("FILE1: ",FILE1)
    NATOMS = 0
    with open(FILE1) as fp:
        for i, line in enumerate(fp):
            if i == 6:
                line = map(int, line.split())
                NATOMS = sum(line)
                print('NATOMS = ',NATOMS)
            elif i == NATOMS + 9:
                GRID = list(map(int, line.split()))
                #print('GRID = ', GRID) 
                POINTS = np.prod(GRID) 
                print('POINTS= ', POINTS)
            elif i > NATOMS + 11:
    
                break
    print ("FILE2: ",FILE2)
    NATOMS2 = 0
    with open(FILE2) as fp:
        for i, line in enumerate(fp):
            if i == 6:
                line = map(int, line.split())
                NATOMS2 = sum(line)
                #print('NATOMS2 = ',NATOMS2)
            elif i == NATOMS2 + 9:
                GRID2 = list(map(int, line.split()))
                POINTS2 = np.prod(GRID2)
                print('POINTS2= ', POINTS2)
            elif i > NATOMS2 + 11:
                break
                
    #print('GRID IS  typre : ',type(GRID)) 
    #print('POINTS IS  typre : ',type(POINTS)) 
    #Read heading of CHGCAR FILE
    CHGCARSUMFILE= open("CHGCARSUM","w+")
    with open(FILE1) as fp:
        for i, line in enumerate(fp):
            if i < NATOMS + 10:
                CHGCARSUMFILE.write(line)
            elif i >= NATOMS + 10:
                CHGCARSUMFILE.close()
                break
    # read data
    skipelines = NATOMS + 9
    maxrows = int(np.ceil(POINTS/5))
    #print ('lines to skipe = ', skipelines)
    #print ('max rows to read = ', maxrows)
    
    #*****CHARGE file 1*****
    CHG1 = pd.read_fwf(FILE1,colspecs='infer', skiprows=skipelines,nrows=maxrows)
    #print('Changing data typre from: ',type(CHG1))
    CHG1=CHG1.to_numpy()
    #print('to: ',type(CHG1))
    
    #Flatten np.Array
    CHG1=CHG1.flatten('C')
    #Leave NaN out of the np.Array
    CHG1=CHG1[0:POINTS]

    # Check data
    print('CHG1',CHG1)
    #print('lenght of ',FILE1,' : ',len(CHG1))
    #print ('first entry in ',FILE1,' : ',CHG1[0])
    #print ('last entry in',FILE1,' : ',CHG1[-1])

    #*****CHARGE file 2*****
    CHG2 = pd.read_fwf(FILE2,colspecs='infer', skiprows=skipelines,nrows=maxrows)
    #print('Changing data typre from: ',type(CHG2))
    CHG2=CHG2.to_numpy()
    #print('to: ',type(CHG2))
    
    #Flatten np.Array
    CHG2=CHG2.flatten('C')
    #Leave NaN out of the np.Array
    CHG2=CHG2[0:POINTS]

    # Check data
    print('CHG2',CHG2)
    #print('lenght of ',FILE2,' : ',len(CHG2))
    #print ('first entry in ',FILE2,' : ',CHG2[0])
    #print ('last entry in',FILE2,' : ',CHG2[-1])
   

    # *****SUM OF CHARGE FILES*****
    
    if POINTS != POINTS2:
        print("The files have different grids")
        return

    CHGSUM = c1 * CHG1 + c2 * CHG2 + c3
    
    #print('CHGSUM', CHGSUM)
    # Number of full arrows and size of last arrow 
    num_full_rows = POINTS // 5
    last_row_size = POINTS % 5

    # Create an appropiate matrix size 
    if last_row_size == 0:
        CHGSUM_matrix = CHGSUM.reshape(-1, 5)
    else:
        CHGSUM_matrix = np.zeros((num_full_rows + 1, 5))
        CHGSUM_matrix[:num_full_rows] = CHGSUM[:num_full_rows * 5].reshape(-1, 5)
        CHGSUM_matrix[-1, :last_row_size] = CHGSUM[num_full_rows * 5:]

    # Save the file with withe spaces instead of zeros if there is a partial last arrow
    with open("CHGCARSUM", "a") as CHGCARSUMFILE:
        for row in CHGSUM_matrix:
            # Si la fila tiene valores correctos (todos los datos), se guarda normalmente
            if np.all(row != 0):
                formatted_row = " ".join(f"{value: .12E}" for value in row)
            else:
                # Convert the arrow to a text chain with white spaces for zero values 
                formatted_row = " ".join(f"{value: .12E}" if value != 0 else "" for value in row)
            CHGCARSUMFILE.write(formatted_row + "\n")
            #CHGCARSUMFILE.close()
    print('Dimensions of CHGARSUM written to a FILE =', CHGSUM_matrix.shape)


############################################################
#           Electron Density Planar Average                #
############################################################
def EDPaverage(FILE1):
    print ("FILE1: ",FILE1)
    NATOMS = 0
    with open(FILE1) as fp:
        for i, line in enumerate(fp):
            if i == 2:
                alattvec = list(map(float, line.split()))
                #print('a latt vector = ', alattvec)
            elif i == 3:
                blattvec = list(map(float, line.split()))
            elif i == 4:
                clattvec = list(map(float, line.split()))
            elif i == 6:
                line = map(int, line.split())
                NATOMS = sum(line) 
                print('NATOMS = ',NATOMS)
            elif i == NATOMS + 9:
                GRID = list(map(int, line.split()))
                #print('GRID = ', GRID)
                POINTS = np.prod(GRID)
                print('POINTS= ', POINTS)
            elif i > NATOMS + 11:
                break

    #Set Grid dimension and such
    NGX = GRID[0]
    NGY = GRID[1]
    NGZ = GRID[2]
    #print("type NGX:",type(NGX))

    # Lattice matrix in direct space in a.u (1 bohr = 0.529177210903 A)

    LATTCURA = (1.0/0.529177210903)*np.dstack([alattvec,blattvec,clattvec])
    # LATTCURA has one more dimension that needed and I don't know much Python
    LATTCURA = LATTCURA[0]

    #print("LATTCURA: ",LATTCURA)

    VOL = np.abs(np.linalg.det(np.dstack([alattvec,blattvec,clattvec])))
    omega = np.dot(np.cross(LATTCURA[1],LATTCURA[2]),LATTCURA[0])
    #print("omega: ",omega)
    LZ = clattvec[2]

    # read data
    skipelines = NATOMS + 9
    maxrows = int(np.ceil(POINTS/5))
    #print ('lines to skipe = ', skipelines)
    #print ('max rows to read = ', maxrows)

    #*****CHARGE file 1*****
    CHG1 = pd.read_fwf(FILE1,colspecs='infer', skiprows=skipelines,nrows=maxrows)
    #print('Changing data typre from: ',type(CHG1))
    CHG1=CHG1.to_numpy()
    #print('to: ',type(CHG1))
    #Flatten np.Array
    CHG1=CHG1.flatten('C')
    #Leave NaN out of the np.Array
    CHG1=CHG1[0:POINTS]
    # Check data
    print('CHG1',CHG1)
    #print('lenght of ',FILE1,' : ',len(CHG1))
    #print ('first entry in ',FILE1,' : ',CHG1[0])
    #print ('last entry in',FILE1,' : ',CHG1[-1])
    
    # reshape array charge array
    GHGtem =  np.zeros((NGX,NGY,NGZ))
    for N3 in range(NGX):
            for N2 in range(NGY):
                        for N1 in range(NGZ):
                            GHGtem[N3,N2,N1] = CHG1[N1*NGX*NGY + N2*NGX + N3]
    #print("DENSITY IN z the fast variable")
    #print(GHGtem[:,0,0])

    CHG = GHGtem/omega
    
    z_axis = np.linspace(0, LZ, NGZ)
    
    # Planar average
    
    AVEGPOTZcorr = []
    for i in range(NGZ):
        AVEGPOTZcorr.append(np.sum(CHG[:, :, i]) / (NGX * NGY))
    #with open('result.dat', 'w') as f:
    #    f.writelines([f'{value}\n' for value in AVEGPOTZcorr])
                
    # Write results to file with two columns
    data_to_save = np.column_stack((z_axis, AVEGPOTZcorr))

    #Format for the columns
    fmt = '%20.6f %20.5e'
 
    # Create header with column names
    header = f'{"#z_axis".rjust(20)}{"Planar_Avg".rjust(20)}'

    # Save with specified formatting
    np.savetxt('PA_VED.dat', data_to_save, fmt=fmt, header=header, comments='')
    plt.clf()
    plt.plot(z_axis,AVEGPOTZcorr,linewidth=2)
    plt.title('Planar Average Electron Density',fontsize=18)
    plt.ylabel(r'$\rho(r) \ (a_{0}^{-3})$',fontsize=12.5)
    plt.xlabel(r'Z-direction ($\AA$)',fontsize=12.5)
    plt.tick_params(axis='both', which='major', labelsize=10)
    #plt.legend()
    plt.tight_layout()
    plt.savefig('PA_VED.svg')


############################################################
#           Electrostatic Potential Planar Average         # 
############################################################

def PoPaverage(FILE1):
    print ("FILE1: ",FILE1)
    NATOMS = 0
    with open(FILE1) as fp:
        for i, line in enumerate(fp):
            if i == 2:
                alattvec = list(map(float, line.split()))
                #print('a latt vector = ', alattvec)
            elif i == 3:
                 blattvec = list(map(float, line.split()))
            elif i == 4:
                 clattvec = list(map(float, line.split()))
            elif i == 6:
                line = map(int, line.split())
                NATOMS = sum(line)
                print('NATOMS = ',NATOMS)
            elif i == NATOMS + 9:
                 GRID = list(map(int, line.split()))
                 #print('GRID = ', GRID)
                 POINTS = np.prod(GRID)
                 print('POINTS= ', POINTS)
            elif i > NATOMS + 11:
                break

    #Set Grid dimension and such
    NGX = GRID[0]
    NGY = GRID[1]
    NGZ = GRID[2]
    #print("type NGX:",type(NGX))

    # Lattice matrix in direct space in a.u (1 bohr = 0.529177210903 A)
    LATTCURA = (1.0/0.529177210903)*np.dstack([alattvec,blattvec,clattvec])
    # LATTCURA has one more dimension that needed and I don't know much Python
    LATTCURA = LATTCURA[0]
    #print("LATTCURA: ",LATTCURA)
    VOL = np.abs(np.linalg.det(np.dstack([alattvec,blattvec,clattvec])))
    omega = np.dot(np.cross(LATTCURA[1],LATTCURA[2]),LATTCURA[0])
    #print("omega: ",omega)
    LZ = clattvec[2]
       
    # read data
    skipelines = NATOMS + 9
    if POINTS % 5 == 0:
        maxrows = int(POINTS / 5)
    elif POINTS % 5 != 0:
        maxrows = int(np.ceil(POINTS / 5)) + 1

    skipelines = skipelines + 1

    if (NGX * NGY * NGZ) % 5 != 0:
        with open(FILE1, 'r') as file:
            last_line = file.readlines()[maxrows + skipelines - 2]
            num_columns_last_line = len(last_line.split())
            #print('La cantidad de columnas de la fila =', maxrows, ' es=', num_columns_last_line)
    else:
        print("No review necessary.")
    
    df = pd.read_table(FILE1, sep='\s+', skiprows=skipelines, names=range(5), nrows=maxrows-1)

    CHG = df.to_numpy().flatten('C')[0:POINTS]

    # Reshape charge array
    GHGtem =  np.zeros((NGX,NGY,NGZ))
    for N3 in range(NGX):
        for N2 in range(NGY):
            for N1 in range(NGZ):
                index = N1 * NGX * NGY + N2 * NGX + N3
                index %= len(CHG)
                GHGtem[N3, N2, N1] = CHG[index]

    #print(GHGtem[:,0,0])
    CHG = GHGtem*(1)

    z_axis = np.linspace(0, LZ, NGZ)
 
    AVEGPOTZcorr = []
    for i in range(NGZ):
        AVEGPOTZcorr.append(np.sum(CHG[:, :, i]) / (NGX * NGY))

    # Write results to file with two columns
    data_to_save = np.column_stack((z_axis, AVEGPOTZcorr))

    #Format for the columns
    fmt = '%20.6f %20.5e'

    # Create header with column names
    header = f'{"#z_axis".rjust(20)}{"Planar_Avg".rjust(20)}'

    # Save with specified formatting
    np.savetxt('PA_Pot.dat', data_to_save, fmt=fmt, header=header, comments='')

    # Plot
    plt.clf()
    plt.plot(z_axis,AVEGPOTZcorr,linewidth=2)
    plt.title('Planar Average Electrostatic Potential',fontsize=18)
    plt.ylabel(r'V(r) (eV)',fontsize=12.5)
    plt.xlabel(r'Z-direction ($\AA$)',fontsize=12.5)
    plt.tick_params(axis='both', which='major', labelsize=10)
    #plt.legend()
    plt.tight_layout()
    plt.savefig('PA_Pot.svg')


############################################################
#                        XYZ-Value CHG                     #
############################################################
def XYZCHG_value(FILE1):    
    print ("FILE1: ",FILE1)
    NATOMS = 0
    with open(FILE1) as fp:
        for i, line in enumerate(fp):
            if i == 2:
                alattvec = list(map(float, line.split()))
                print('a latt vector = ', alattvec)
            elif i == 3:
                blattvec = list(map(float, line.split()))
            elif i == 4:
                clattvec = list(map(float, line.split()))
            elif i == 6:
                line = map(int, line.split())
                NATOMS = sum(line)
                print('NATOMS = ',NATOMS)
            elif i == NATOMS + 9:
                GRID = list(map(int, line.split()))
                print('GRID = ', GRID)
                POINTS = np.prod(GRID)
                print('POINTS= ', POINTS)
            elif i > NATOMS + 11:
                 break

    #Set Grid dimension and such
    NGX = GRID[0]
    NGY = GRID[1]
    NGZ = GRID[2]
    print("type NGX:",type(NGX))

    # Lattice matrix in direct space in a.u (1 bohr = 0.529177210903 A)
    LATTCURA = (1.0/0.529177210903)*np.dstack([alattvec,blattvec,clattvec])
    # LATTCURA has one more dimension that needed and I don't know much Python
    LATTCURA = LATTCURA[0]
    #print("LATTCURA: ",LATTCURA)
    LATTCURADO = (1.0) * np.dstack([alattvec, blattvec, clattvec])
    LATTCURADO = LATTCURADO[0]
    VOL = np.abs(np.linalg.det(np.dstack([alattvec,blattvec,clattvec])))
    omega = np.dot(np.cross(LATTCURA[1],LATTCURA[2]),LATTCURA[0])
    #print("omega: ",omega)
    LZ = clattvec[2]
    # read data
    skipelines = NATOMS + 9
    maxrows = int(np.ceil(POINTS/5))
    #print ('lines to skipe = ', skipelines)
    #print ('max rows to read = ', maxrows)

    #*****CHARGE file 1*****
    CHG1 = pd.read_fwf(FILE1,colspecs='infer', skiprows=skipelines,nrows=maxrows)
    #print('Changing data typre from: ',type(CHG1))
    CHG1=CHG1.to_numpy()
    #print('to: ',type(CHG1))
    #Flatten np.Array
    CHG1=CHG1.flatten('C')
    #Leave NaN out of the np.Array
    CHG1=CHG1[0:POINTS]
    # Check data
    #print('CHG1',CHG1)
    #print('lenght of ',FILE1,' : ',len(CHG1))
    #print ('first entry in ',FILE1,' : ',CHG1[0])
    #print ('last entry in',FILE1,' : ',CHG1[-1])


    # reshape array charge array
    GHGtem =  np.zeros((NGX,NGY,NGZ))
    for N3 in range(NGX):
        for N2 in range(NGY):
            for N1 in range(NGZ):
                GHGtem[N3,N2,N1] = CHG1[N1*NGX*NGY + N2*NGX + N3]
    #print(GHGtem[:,0,0])
    CHG = GHGtem/omega
   
    # Start adding
    suma = 0.0

    with open('xyz_value.dat', 'w') as xyz_file:
        xyz_file.write(f"{NGX * NGY * NGZ}\n")
        xyz_file.write("Generated XYZ file\n")

        for N3 in range(NGZ):
            for N2 in range(NGY):
                for N1 in range(NGX):
                    x = (N1 * LATTCURADO[0, 0] + N2 * LATTCURADO[1, 0] + N3 * LATTCURADO[2, 0]) / NGX
                    y = (N1 * LATTCURADO[0, 1] + N2 * LATTCURADO[1, 1] + N3 * LATTCURADO[2, 1]) / NGY
                    z = (N1 * LATTCURADO[0, 2] + N2 * LATTCURADO[1, 2] + N3 * LATTCURADO[2, 2]) / NGZ

                    # Compute the density value
                    valor_densidad = CHG[N1, N2, N3]
                    # Add to the total
                    suma += valor_densidad

                    xyz_file.write(f"{x:>12.4f}{y:>12.4f}{z:>12.4f}{CHG[N1, N2, N3]:>20.8e}\n")

    # Print the sum
    print(f"Summing up all value in grid file: {suma}")
    Elemento_volumen = omega/(NGX*NGY*NGZ)
    print(f"Differential element: {Elemento_volumen}")
    Integral_rho = suma*Elemento_volumen
    print(f"After multiplied by differential element: {Integral_rho}")


############################################################
#                        XYZ-Value LOC                     #
############################################################
def XYZLOC_value(FILE1,c1,c2):
    print ("FILE1: ",FILE1)
    NATOMS = 0
    with open(FILE1) as fp:
                for i, line in enumerate(fp):
                    if i == 2:
                        alattvec = list(map(float, line.split()))
                        print('a latt vector = ', alattvec)
                    elif i == 3:
                        blattvec = list(map(float, line.split()))
                    elif i == 4:
                        clattvec = list(map(float, line.split())) 
                    elif i == 6:
                        line = map(int, line.split())
                        NATOMS = sum(line)
                        print('NATOMS = ',NATOMS)
                    elif i == NATOMS + 9:
                        GRID = list(map(int, line.split()))
                        print('GRID = ', GRID)
                        POINTS = np.prod(GRID)
                        print('POINTS= ', POINTS)
                    elif i > NATOMS + 11:
                        break

    #Set Grid dimension and such
    NGX = GRID[0]
    NGY = GRID[1]
    NGZ = GRID[2]
    #print("type NGX:",type(NGX))

    # Lattice matrix in direct space in a.u (1 bohr = 0.529177210903 A)
    LATTCURA = (1.0/0.529177210903)*np.dstack([alattvec,blattvec,clattvec])
    # LATTCURA has one more dimension that needed and I don't know much Python
    LATTCURA = LATTCURA[0]
    #print("LATTCURA: ",LATTCURA)
    LATTCURADO = (1.0) * np.dstack([alattvec, blattvec, clattvec])
    LATTCURADO = LATTCURADO[0]
    VOL = np.abs(np.linalg.det(np.dstack([alattvec,blattvec,clattvec])))
    omega = np.dot(np.cross(LATTCURA[1],LATTCURA[2]),LATTCURA[0])
    #print("omega: ",omega)
    LZ = clattvec[2]
    
    skipelines = NATOMS + 9
    if POINTS % 5 == 0:
        maxrows = int(POINTS / 5)
    elif POINTS % 5 != 0:
        maxrows = int(np.ceil(POINTS / 5)) + 1
    
    skipelines = skipelines + 1

    #print ('lines to skipe = ', skipelines)
    #print ('max rows to read = ', maxrows)


    if (NGX * NGY * NGZ) % 5 != 0:
        with open(FILE1, 'r') as file:
            last_line = file.readlines()[maxrows + skipelines - 2]
            num_columns_last_line = len(last_line.split())
            #print('La cantidad de columnas de la fila =', maxrows, ' es=', num_columns_last_line)
    else:
        print("No review necessary.")

    df = pd.read_table(FILE1, sep='\s+', skiprows=skipelines, names=range(5), nrows=maxrows - 1)

    CHG = df.to_numpy().flatten('C')[0:POINTS]

    # Reshape charge array
    GHGtem =  np.zeros((NGX,NGY,NGZ))
    for N3 in range(NGX):
        for N2 in range(NGY):
            for N1 in range(NGZ):
                index = N1 * NGX * NGY + N2 * NGX + N3
                index %= len(CHG)
                GHGtem[N3, N2, N1] = CHG[index]
    
    #print(GHGtem[:,0,0])
    CHG = GHGtem*(c1)
    
    AVEGPOTZcorr = []
    for i in range(NGZ):
        AVEGPOTZcorr.append(np.sum(CHG[:, :, i]) / (NGX * NGY))
    first_value_AVEGPOTZcorr = AVEGPOTZcorr[0]
   
    #print('primer valor:',first_value_AVEGPOTZcorr)
    
    if c2 == 1:
        CHG = CHG - first_value_AVEGPOTZcorr
    else:
        CHG = CHG*(1)

    #CHG = CHG - first_value_AVEGPOTZcorr

    AVEGPOTZ = []
    for i in range(NGZ):
        AVEGPOTZ.append(np.sum(CHG[:, :, i]) / (NGX * NGY))

    plt.plot(AVEGPOTZ)
    plt.xlabel('Z-direction')
    plt.ylabel('Planar Average of LOCPOT')
    plt.title('Electrostatic Potential')
    plt.savefig('PA_MEP.svg', format='svg', bbox_inches='tight')

    # Start adding
    suma = 0.0
    with open('xyz_value.dat', 'w') as xyz_file:
        xyz_file.write(f"{NGX * NGY * NGZ}\n")
        xyz_file.write("Generated XYZ file\n")

        for N3 in range(NGZ):
            for N2 in range(NGY):
                for N1 in range(NGX):
                    x = (N1 * LATTCURADO[0, 0] + N2 * LATTCURADO[1, 0] + N3 * LATTCURADO[2, 0]) / NGX
                    y = (N1 * LATTCURADO[0, 1] + N2 * LATTCURADO[1, 1] + N3 * LATTCURADO[2, 1]) / NGY
                    z = (N1 * LATTCURADO[0, 2] + N2 * LATTCURADO[1, 2] + N3 * LATTCURADO[2, 2]) / NGZ
                    # Compute the density value
                    valor_densidad = CHG[N1, N2, N3]
                    # Add to the total
                    suma += valor_densidad
                    xyz_file.write(f"{x:>12.4f}{y:>12.4f}{z:>12.4f}{CHG[N1, N2, N3]:>20.8e}\n")

    # Print the sum
    print(f"Summing up all value in grid file: {suma}")
    Elemento_volumen = omega/(NGX*NGY*NGZ)
    print(f"Differential element: {Elemento_volumen}")
    Integral_rho = suma*Elemento_volumen
    print(f"After multiplied by differential element: {Integral_rho}")


############################################################
#             Perturbative Perspectives                    #
############################################################

def Perturbative_point(FILE0,FILE1,q,N):
    # read number of atoms and grid points
    print ("FILE0: ",FILE0)
    NATOMS = 0
    with open(FILE0) as fp:
        for i, line in enumerate(fp):
            if i == 6:
                line = map(int, line.split())
                NATOMS = sum(line)
                print('NATOMS = ',NATOMS)
            elif i == NATOMS + 9:
                GRID = list(map(int, line.split()))
                #print('GRID = ', GRID)
                POINTS = np.prod(GRID)
                print('POINTS= ', POINTS)
            elif i > NATOMS + 11:
                 break
 
    #print ("FILE0: ",FILE0)
    with open(FILE0) as fp:
        for i, line in enumerate(fp):
            if i == 2:
                alattvec = list(map(float, line.split()))
                #print('a latt vector = ', alattvec)
            elif i == 3:
                blattvec = list(map(float, line.split()))
            elif i == 4:
                clattvec = list(map(float, line.split()))
                break

    #Set Grid dimension and such
    NGX = GRID[0]
    NGY = GRID[1]
    NGZ = GRID[2]
    #print("type NGX:",type(NGX))
                        
    CHGCARSUMFILE= open("MODELPOT.LOCPOT","w+")
    with open(FILE1) as fp:
        for i, line in enumerate(fp):
            if i < NATOMS + 10:
                CHGCARSUMFILE.write(line)
            elif i >= NATOMS + 10:
                CHGCARSUMFILE.close()
                break

    skipelines = NATOMS + 9
    maxrows = int(np.ceil(POINTS/5))

    LZ = clattvec[2]

    #*****CHARGE file 1*****
    CHG1 = pd.read_fwf(FILE0,colspecs='infer', skiprows=skipelines,nrows=maxrows)
    
    CHG1 = CHG1.to_numpy().flatten('C')[:POINTS]
    CHG1 = CHG1.astype(np.float64)
    
    # Check data
    print('CHG1',CHG1)
    '''
    #print('Changing data typre from: ',type(CHG1))
    CHG1=CHG1.to_numpy()
    CHG1=CHG1.flatten('C')
    #Leave NaN out of the np.Array
    CHG1=CHG1[0:POINTS]
    # Check data
    print('CHG1',CHG1)
    '''
    #*****CHARGE file 2*****
    CHG2 = pd.read_fwf(FILE1,colspecs='infer', skiprows=skipelines,nrows=maxrows)
    
    CHG2 = CHG2.to_numpy().flatten('C')[:POINTS]
    CHG2 = CHG2.astype(np.float64)
    
    # Check data
    print('CHG2',CHG2)

    '''
    #print('Changing data typre from: ',type(CHG2))
    CHG2=CHG2.to_numpy()
    #print('to: ',type(CHG2))

    #Flatten np.Array
    CHG2=CHG2.flatten('C')    
    #Leave NaN out of the np.Array
    CHG2=CHG2[0:POINTS]
    # Check data
    print('CHG2',CHG2)
    '''
    print("")
    print("Just a few seconds.")
    ###############################
    def reshape_xyz(NGX,NGY,NGZ, CHG):
        GHGtem =  np.zeros((NGX,NGY,NGZ))
        for N3 in range(NGX):
            for N2 in range(NGY):
                for N1 in range(NGZ):
                    GHGtem[N3,N2,N1] = CHG[N1*NGX*NGY + N2*NGX + N3]
        return GHGtem

    CHGtem1 = reshape_xyz(NGX, NGY, NGZ, CHG1)
    CHGtem2 = reshape_xyz(NGX, NGY, NGZ, CHG2)

    ###############################
    # *****SUM OF CHARGE FILES*****

    CHGSUM = q*CHGtem1*(-1) + N*q*CHGtem2*(-1)
    ##############################

    AVEGPOTZcorr = []
    for i in range(NGZ):
        AVEGPOTZcorr.append(np.sum(CHGSUM[:, :, i]) / (NGX * NGY))
    first_value_AVEGPOTZcorr = AVEGPOTZcorr[0]
    #print('AVEGPOTZcorr:', AVEGPOTZcorr)
    #print('primer valor:',first_value_AVEGPOTZcorr)
    
    CHG = CHGSUM - first_value_AVEGPOTZcorr
    
    AVEGPOTZchg = []
    for i in range(NGZ):
        AVEGPOTZchg.append(np.sum(CHG[:, :, i]) / (NGX * NGY))
    first_value_AVEGPOTZchg = AVEGPOTZchg[0]
    #print('primer valor:',first_value_AVEGPOTZchg)
    #print('AVEGPOTZchg:',  AVEGPOTZchg)
    
    z_axis = np.linspace(0, LZ, NGZ)
    
    '''
    # Write results to file with two columns
    #data_to_save = np.column_stack((z_axis, AVEGPOTZcorr))
    #Format for the columns
    #fmt = '%20.6f %20.5e'
    # Create header with column names
    #header = f'{"#z_axis".rjust(20)}{"Planar_Avg".rjust(20)}'
    # Save with specified formatting
    #np.savetxt('PA_VED.dat', data_to_save, fmt=fmt, header=header, comments='')
    
    # Plot
    plt.clf()
    plt.plot(z_axis,AVEGPOTZcorr,label='No corregido',linewidth=2)
    plt.plot(z_axis,AVEGPOTZchg,label='Corregido',linewidth=2)
    plt.title('Planar Average Electrostatic Potential',fontsize=18)
    plt.ylabel(r'V(r) (eV)',fontsize=12.5)
    plt.xlabel(r'Z-direction ($\AA$)',fontsize=12.5)
    plt.tick_params(axis='both', which='major', labelsize=10)
    plt.legend()
    plt.tight_layout()
    plt.savefig('PA_Model_U.svg')
    '''
    # Number of full arrows and size of last arrow
    num_full_rows = POINTS // 5
    last_row_size = POINTS % 5

    CHGCARSUMFILE= open("MODELPOT.LOCPOT","w+")
    with open(FILE1) as fp:
        for i, line in enumerate(fp):
            if i < NATOMS + 10:
                CHGCARSUMFILE.write(line)
            elif i >= NATOMS + 10:
                CHGCARSUMFILE.close()
                break
    ### write de locpot in a suitiable form to be written in vasp format (x is de fast index)

    LOCPOTtem =  np.zeros(NGX*NGY*NGZ)
    for N1 in range(NGZ):
        for N2 in range(NGY):
            for N3 in range(NGX):
                LOCPOTtem[N1*NGX*NGY + N2*NGX + N3] = CHG[N3,N2,N1]

    #print("*****LOCPOTtem****")
    #print(LOCPOTtem)
    num_elements = CHG.size
    #print(f"El número de elementos en LOCPOT_cor2 es: {num_elements}")

    # Crear una matriz de tamaño adecuado
    if last_row_size == 0:
        LOCPOT = LOCPOTtem
        FUKUIPOT = np.resize(LOCPOT, (maxrows, 5))
        print ('Dimensions of Interaction Model  written to a file = ',FUKUIPOT.shape)
        with open("MODELPOT.LOCPOT", "a") as FILE:
            for row in FUKUIPOT:
                formatted_row = " " + " ".join(f"{value: .12E}" for value in row)
                FILE.write(formatted_row + "\n")
    else:
         LOCPOT = LOCPOTtem
         #print('LOCPOT size',LOCPOT.size)
         LOCPOT_matrix = np.zeros((num_full_rows + 1, 5))
         LOCPOT_matrix[:num_full_rows] = LOCPOT[:num_full_rows * 5].reshape(-1, 5)
         LOCPOT_matrix[-1, :last_row_size] = LOCPOT[num_full_rows * 5:]
         #print("VER ESTO")
         #print('LOCPOT_matrix', LOCPOT_matrix)

         with open("MODELPOT.LOCPOT", "a") as FILE:
             for row in LOCPOT_matrix:
                 # If the row has correct values <200b><200b>(all data), it is saved normally
                 if np.all(row != 0):
                     formatted_row = " ".join(f"{value: .12E}" for value in row)
                 else:
                     # Convert row to a text string with blanks for zero values
                     formatted_row = " ".join(f"{value: .12E}" if value != 0 else "" for value in row)
                 FILE.write(formatted_row + "\n")

         print('Dimensions of Interaction Model written to a FILE =', LOCPOT_matrix.shape)
         print("")

    '''
    # Crear una matriz de tamaño adecuado
    if last_row_size == 0:
        LOCPOT = CHG
        FUKUIPOT = np.resize(LOCPOT, (maxrows, 5))
        print ('Dimensions of FUKUI POTENTIAL  written to a file = ',FUKUIPOT.shape)
        
        with open("MODELINT.LOCPOT", "a") as FILE:
            for row in FUKUIPOT:
                formatted_row = " " + " ".join(f"{value: .12E}" for value in row)
                FILE.write(formatted_row + "\n")
    else:
        LOCPOT = CHG
        print('VER ESTO')
        print('LOCPOT',LOCPOT)
        LOCPOT_matrix = np.zeros((num_full_rows + 1, 5))
        LOCPOT_matrix[:num_full_rows] = LOCPOT[:num_full_rows * 5].reshape(-1, 5)
        LOCPOT_matrix[-1, :last_row_size] = LOCPOT[num_full_rows * 5:]
    
        with open("MODELINT.LOCPOT", "a") as FILE:
            for row in LOCPOT_matrix:
                formatted_row = " " + " ".join(f"{value: .12E}" for value in row)
                #formatted_row = " ".join(f"{value: .12E}" if value != 0 else " " * 12 for value in row)
                FILE.write(formatted_row + "\n")
     '''
##################################################################################################

    
############################################################
#                 Main Menu                                #
############################################################

def main_menu():
        while True:

            #******Print Header*****
            print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
            print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@     @@@@@@@@@@@@@@@@@@@@@@@@@@@     @@@@")
            print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ @@@@")
            print("@@@     @@@@@@@@@@@@@@@@@@@     @@@ @@@@  @@@@@@ @@@@@@@@@@@@@@@@@@@@@@ @@@@")
            print("@@@ @@@@@@@@@@@       @@@@@ @@@@@@@ @@@@ @ @@@@@ @@@@@ @@@@@@@@@@ @@@@@ @@@@")
            print("@@@ @   @@@@@@@@@@@@@@@@@@@ @  @@@@ @@@@ @@ @@@@ @@@@@@ @@@@@@@@ @@@@@@ @@@@")
            print("@@@ @@@@@@@@@@@       @@@@@ @@@@@@@ @@@@ @@@@ @@ @@@@@@@ @@@@@@ @@@@@@@ @@@@")
            print("@@@     @@@@@@@@@@@@@@@@@@@     @@@ @@@@ @@@@@   @@@ @@@@@ @@@ @@@@@@@@ @@@@")
            print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ @@@@@@@@@@@@@@@ @@@@@@@ @ @@@@@@@@@ @@@@")
            print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@     @@@@@@@@@@@@@@@@@@@@@@@@@@@     @@@@")
            print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
            print("")
            print("Fukui-Tayrona -- A useful tool for Conceptual DFT in Solid-State.")
            print("Version 1.0, release date: September-2024")
            print("Developers: Javiera Cabezas-Escares, Nicolas F.Barrera, Prof. Carlos Cardenas -- ")
            print("            TheoChemPhys Group, University of Chile.")
            print("            http://tayrona.ciencias.uchile.cl/")
            print("")
            print("")
            print("")
            print("               ***** Main Menu *****               ")
            print("1 -- Fukui Function via Interpolation method")
            print("2 -- Fukui Potential via Electrodes' method")
            print("3 -- Fukui Potential via SCPC method")
            print("4 -- Process Grid Data")
            print("5 -- Perturbative Perspectives")
            print("6 -- Exit")
            print("")
            option = input("Choose an option: ")

            if option == "1":
                print("")
                print("Chose option 1: Fukui Function via Interpolation.")
                print("")
                print("11 Electrophilic Fukui function f^-(r).")
                print("12 Nucleophilic  Fukui function f^+(r).")

                option1 = input("Choose an option: ")
                if option1 == "11":
                    print("Name CHGCAR files with \u03B4N: -0.15, -0.10, -0.05, and 0.0")
                    FILE1 = input("Enter name of file 1: ")
                    FILE2 = input("Enter name of file 1: ")
                    FILE3 = input("Enter name of file 1: ")
                    FILE4 = input("Enter name of file 1: ")
                    print("")
                    print("This will take a few seconds.")
                    print("")

                    Electro_fukui_interpolation(FILE1, FILE2, FILE3, FILE4)

                    continue
    
                if option1 == "12":
                    print("Name CHGCAR files with \u03B4N: 0.0, +0.05, +010, and +0.15")
                    FILE1 = input("Enter name of file 1: ")
                    FILE2 = input("Enter name of file 1: ")
                    FILE3 = input("Enter name of file 1: ")
                    FILE4 = input("Enter name of file 1: ")
                    print("")
                    print("This will take a few seconds.")
                    print("")

                    Nucleo_fukui_interpolation(FILE1, FILE2, FILE3, FILE4)
                    
                    continue

                else:
                    print("Invalid option. Please select an option from the menu.")

            elif option == "2":
                print("You selected option 2: Fukui Potential via Electrodes' method")
                print("")
                print("Name CHGCAR file of charge density of the neutral slab.")
                FILE0 = input("Enter file name: ")
                print("")
                print("Name CHGCAR file of Fukui function.")
                FILE1 = input("Enter file name: ")
                print("")
                print("Dielectric constant value.")
                Epsilon = float(input("Value: ")) 
                fukui_electrodes(FILE0, FILE1, Epsilon)
                
                continue
            

            elif option == "3":
                print("You selected option4: ")
                print("")
                print("")
                print("31 -- Electrophilic Fukui function v_f^-(r).")
                print("32 -- Nucleophilic Fukui function v_f^+(r).")
                print("")

                option3 = input("Choose an option: ")
                if option3 == "31":
                    print("Name CHGCAR file of charge density of the neutral slab.")
                    FILE0 = input("Enter file name: ")
                    print("")
                    print("Name CHGCAR file of Fukui function.")
                    FILE1 = input("Enter file name: ")
                    print("")
                    print("Name SCPC correcton file.")
                    FILE2 = input("Enter file name: ")
                    c = -1
                    print("")

                    vfukui_SCPC(FILE0,FILE1,FILE2,c)

                    continue
                
                if option3 == "32":
                    print("Name CHGCAR file of charge density of the neutral slab.")
                    FILE0 = input("Enter file name: ")
                    print("")
                    print("Name CHGCAR file of Fukui function.")
                    FILE1 = input("Enter file name: ")
                    print("")
                    print("Name SCPC correcton file.")
                    FILE2 = input("Enter file name: ")
                    c = 1
                    print("")

                    vfukui_SCPC(FILE0,FILE1,FILE2.c)
                    
                    continue

                else:
                    print("Invalid option. Please select an option from the menu.")

            elif option == "4":
                print("You selected option 3: Process Grid Data")
                print("")
                print("")
                print("41 -- Add:            Add two CHGCAR or LOCPOT to produce a new one.")
                print("42 -- Subtract:       Subtracts two CHGCAR or LOCPOT to produce a new one.")
                print("43 -- Scale:          Scale a CHGCAR or POTCAR to produce a new one.")
                print("44 -- Add a constant: Adds a constant to a CHGCAR or LOCPOT file. ")
                print("45 -- Planar Average: Compute the planar average of CHGCAR or LOCPOT.")
                print("46 -- Value-xyz:      Convert a CHGCAR or LOCPOT into a simple list: X,Y,Z,Value format.")
                print("")
                print("")

                option4 = input("Choose an option: ")
                if option4 == "41":
                    print("Name CHGCAR or LOCPOT files to add: ")
                    FILE1 = input("Enter name of file 1: ")
                    FILE2 = input("Enter name of file 2: ")
                    print("Enter the constants c1 and c2 according to FILE1*c1 + FILE2*c2: ")
                    c1 = float(input("c1: "))
                    c2 = float(input("c2: "))
                    c3 = 0
                    print("")
                    print("This will take a few seconds.")
                    print("")
                    
                    Add(FILE1,FILE2,c1,c2,c3)
                   
                    # Change the file name
                    old_filename = "CHGCARSUM"
                    new_filename = "CHGCAR_SUM"
                    os.rename(old_filename, new_filename)

                    continue

                if option4 == "42":
                    print("Name CHGCAR or LOCPOT files to subtract: ")
                    FILE1 = input("Enter name of file 1: ")
                    FILE2 = input("Enter name of file 2: ")
                    c1 = 1
                    c2 = -1
                    c3 = 0
                    print("")
                    print("This will take a few seconds.")
                    print("")

                    Add(FILE1,FILE2,c1,c2,c3)
                    
                    # Change the file name
                    old_filename = "CHGCARSUM"
                    new_filename = "CHGCAR_DIFF"
                    os.rename(old_filename, new_filename)

                    continue

                if option4 == "43":
                    print("Name CHGCAR or LOCPOT files to scale: ")
                    FILE1 = input("Enter name of file 1: ")
                    FILE2 = FILE1
                    print("")
                    print("Enter the scale factor c:")
                    c1 = float(input("c: "))
                    c2 = 0
                    c3 = 0

                    Add(FILE1,FILE2,c1,c2,c3)

                    # Change the file name
                    old_filename = "CHGCARSUM"
                    new_filename = "CHGCAR_SCALE"
                    os.rename(old_filename, new_filename)
                    
                    continue
                
                if option4 == "44":
                    print("Name CHGCAR or LOCPOT files to add it a constant to:")
                    FILE1 = input("Enter name of file 1: ")
                    FILE2 =  FILE1
                    print("")
                    print("Enter the constant c:")
                    c3 = float(input("c: "))
                    c1 = 1
                    c2 = 0

                    Add(FILE1,FILE2,c1,c2,c3)

                    # Change the file name
                    old_filename = "CHGCARSUM"
                    new_filename = "CHGCAR_constant"
                    os.rename(old_filename, new_filename)

                    continue

                if option4 == "45":
                   print("")
                   print("451 -- Planar Average for a CHGCAR file.") 
                   print("452 -- Planar Average for a LOCPOT file.")
                   print("")
                   option45 = input("Choose an option: ")
                   
                   if option45 == "451":
                       print("")
                       FILE1 = input("Enter name of CHGCAR file: ")
                       EDPaverage(FILE1)

                       continue

                   if option45 == "452":
                       print("")
                       FILE1 = input("Enter name of LOCPOT file: ")
                       PoPaverage(FILE1)

                       continue

                if option4 == "46":
                    print("")
                    print("461 -- XYZ-value for a CHGCAR file.")
                    print("462 -- XYZ-value for a LOCPOT file.")
                    print("")
                    option46 = input("Choose an option:")

                    if option46 == "461":
                        print("")
                        FILE1 = input("Enter name of CHGCAR file: ")
                        XYZCHG_value(FILE1)

                        continue
                    
                    if option46 == "462":
                        print("")
                        FILE1 = input("Enter name of LOCPOT file: ")
                        print("")
                        print("Do you want to align potential to zero in the vacuum?")
                        option462 = input("yes or no: ")

                        if option462 == "yes" or option462 == "YES" or option462 == "y":
                            print("")
                            c2 = 1
                            print("Do you want to write the data in the chemical convention (i.e. take the negative)?")
                            option462y = input("yes or no: ")
                            if option462y == "yes" or option462y == "YES" or option462y == "y":
                                print("")
                                c1 = -1
                                XYZLOC_value(FILE1,c1,c2)

                                continue

                            elif option462y == "no" or option462y == "NO" or option462y == "n":
                                print("")
                                c1 = 1
                                XYZLOC_value(FILE1,c1,c2)

                                continue

                        if option462 == "no" or option462 == "NO" or option462 == "n":
                            print("")
                            c2 = 2
                            print("Do you want to write the data in the chemical convention (i.e. take the negative)?")
                            option462n = input("yes or no: ")
                            if option462n == "yes" or option462 == "YES" or option462 == "y":
                                print("")
                                c1 = -1
                                XYZLOC_value(FILE1,c1,c2)
                            
                            elif option462n == "no" or option462n == "NO" or option462n == "n":
                                 print("")
                                 c1 = 1
                                 XYZLOC_value(FILE1,c1,c2)
                                 continue
        
            elif option == "5":
                print("")
                print("You selected Perturbative Perspectives.")
                mu_pm = "\u03BC\u207A\u207B"
                #Delta_U = "\u0394U = " + mu_pm + "\u0394N + q\u03A6(r) - q\u0394Nvf\u207A\u207B(r)"
                #print(Delta_U)
                Delta_U = "\u0394U(r) = q\u03A6(r) - q\u0394Nvf\u207A\u207B(r)"
                print(Delta_U)
                print("Name of LOCPOT file with Electrostatic potential \u03A6(r).")
                FILE0 = input("Enter name of LOCPOT: ")
                print("")
                print("Name of LOCPOT file with Fukui potential vf\u207A\u207B(r).")
                FILE1 = input("Enter name of LOCPOT: ")
                print("")
                print("Enter the change in the number of electrons \u0394N: ")
                N = float(input("\u0394N: "))
                print("")
                print("Enter the charge q of active site: ")
                q = float(input("q: "))
                print("")
                Perturbative_point(FILE0,FILE1,q,N)

                continue


            elif option == "6":
                print("You selected option 6: Goodbye!")
                sys.exit()

            else:
                print("Invalid option. Please select an option from the menu.")

if __name__ == "__main__":
        main_menu()
