# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 21:33:12 2016

@author: mali
"""
from enum import Enum

#Region ExtraRayData Filters
class ExtraRayData(Enum):
    Source=1                 #1.0 if filter condition met, 0.0 otherwise.
    Surface =2               #1.0 if filter condition met, 0.0 otherwise.
    After_Surface=3          #1.0 if filter condition met, 0.0 otherwise.
    Element=4                #1.0 if filter condition met, 0.0 otherwise.
    After_Element=5          #1.0 if filter condition met, 0.0 otherwise.
    Property_Zone=6          #1.0 if filter condition met, 0.0 otherwise.
    After_Property_Zone=7    #1.0 if filter condition met, 0.0 otherwise.
    Hit_Number=8            #Number of times ray hit the receiver, as a real.
    Ray_Magnitude=9         #Ray magnitude of ray. -- RayDataMagnitude with LT.GetReceiverRayData is faster
    Wavelength=10             #Wavelength of ray (in nm).  -- RayDataWavelength with LT.GetReceiverRayData is faster
    Incident_Angle=11         #Angle of incidence of ray on receiver (in degrees).
    Exit_Angle=12             #Angle of exit of ray from receiver (in degrees).
    Path_Transmittance=13     #Path transmitance of ray.
    Volume_Interface=14       #1.0 if filter condition met, 0.0 otherwise.
    After_Volume_Interface=15 #1.0 if filter condition met, 0.0 otherwise.
    Optical_Path_Length=16    #Cumulative optical path length of ray at receiver.
    Optical_Property=17       #1.0 if filter condition met, 0.0 otherwise.
    After_Optical_Property=18 #1.0 if filter condition met, 0.0 otherwise.
    Polarization=19           #Value of the PolType
    RayPathIndex=20           #Base 0 Ray Path Index. -- RayPathIndex with LT.GetReceiverRayData is faster
    Filter_Group=21           #1.0 if filter condition met, 0.0 otherwise.
#End extra ray data filters
import clr
import math
import numpy as np
import System
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#This is the local LT pointer in this module.
#This must be set to a valid instance of LTCOM64.LTAPIx, after importing this module into another module
lt0=None 
ltu=None

def cart2sph(x,y,z):
    """Returns the azimuth, elevation, and magnitude from x, y, z"""
    azimuth = np.arctan2(y,x)
    elevation = np.arctan2(z,np.sqrt(x**2 + y**2))
    r = np.sqrt(x**2 + y**2 + z**2)
    return azimuth, elevation, r

pass

def unit_vector(vector):
    """Returns the unit vector of the vector."""
    return vector / np.linalg.norm(vector)

pass

def angle_between_vectors(v1, v2):
    """ Returns the angle between vectors 'v1' and 'v2' in radians::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            
            1.5707963267948966
            
            >>> angle_between((1, 0, 0), (1, 0, 0))
            
            0.0
            
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            
            3.141592653589793
    """
    try:
        v1_u = unit_vector(v1)
        v2_u = unit_vector(v2)
        ang=np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))
        return ang        
    except Exception as e:
        print('angle_between_vectors: ' + str(e))
pass

def GetKeyAndDataNameFromString(DataAccessString=''):
    """This is a utility function used for functions in this module

    """
    if DataAccessString=='':
        return ''
    else:
        DataAccessString=DataAccessString.strip()
        das=DataAccessString.split('.')
        dataname=das[len(das)-1]
        fulldataname=dataname
        #Strings for grids can contain [][] 
        tstr=dataname.split('[')
        dataname=tstr[0]
        key=DataAccessString[:len(DataAccessString)-len(fulldataname)-1]
        #print(key,dataname)
        return key,dataname

def GetLTDbItem(DataAccessString,returnstatus=-1, dbname=''):
    """Get a database item value
    
    Parameters
    ----------
    DataAccessString: String
        This is the string you get via Copy Data Access Name in LightTools
    returnstatus: Integer
        This is optional. Default is -1, no return status. Pass 0 to request the return status
    dbname: String
        This is an option to pass the datakey, dataname combination. When empty, this is not used
    
    Returns
    -------
    Requested data item
        Usually a floating point number or a string
    Status of the command execution (optional)
        An integer that matches LTReturncodeENUMs
    
    Examples
    --------
    Call without a return status
        L = GetLTDbItem('Solid[1].Primitive[1].Length')
    Call with the return status
        L,Stat = GetLTDbItem('Solid[1].Primitive[1].Length',0)
    Call with data key, data name combo
        L = GetLTDbItem('Solid[1].Primitive[1]', dbname='Length')
    """
    try:
        if dbname == '':
            key,dataname=GetKeyAndDataNameFromString(DataAccessString)
        else:
            dataname=dbname
            key=DataAccessString
        pass
        #print(key,dataname)
        [ltdata,stat]=lt0.DbGet(key,dataname,0)
        if returnstatus==0:
            return ltdata,stat
        else:
            return ltdata
    except Exception as e:
        print('GetLTDbItem: ' + str(e))
pass
##Local test
#Call without a return status
#L = GetLTDbItem('Solid[1].Primitive[1].Length')
#Call with the return status
#L,Stat = GetLTDbItem('Solid[1].Primitive[1].Length,0)
pass

def SetLTDbItem(DataAccessString,datavalue,dbname=''):
    """Set a database item value
    
    Parameters
    ----------
    DataAccessString: String
        This is the string you get via Copy Data Access Name in LightTools
    datavalue: string or numeric
        This is the new data value assigned to the database item
    dbname: String
        This is an option to pass the datakey, dataname combination. When empty, this is not used

    
    Returns
    -------
    Status of the command execution (optional)
        An integer that matches LTReturncodeENUMs
    
    Examples
    --------
        stat=SetLTDbItem('solid[1].primitive[1].radius',5)
    """
    try:
        if dbname == '':
            key,dataname=GetKeyAndDataNameFromString(DataAccessString)
        else:
            dataname=dbname
            key=DataAccessString
        pass
        stat=lt0.DbSet(key,dataname,datavalue)
        return stat
    except Exception as e:
        print('SetLTDbItem: ' + str(e))
pass
##Local test
#Call without a return status
#L = GetLTDbItem('Solid[1].Primitive[1].Length)
#Call with the return status
#L,Stat = GetLTDbItem('Solid[1].Primitive[1].Length,0)
pass

def GetLTGridItem(DataAccessString,RowIndex=1,ColIndex=-1,returnstatus=-1,dbname=''):
    """Get a data value from a 1D or 2D grid item value
    
    Parameters
    ----------
    DataAccessString: String
        This is the string you get via Copy Data Access Name in LightTools
        Note: Content in square brackets at the end (e.g. [4] or [4][6]) is ignored)
    RowIndex: Integer
        Index of the row in the grid, starting 1
    ColIndex: Integer
        Index of the column in the grid, starting 1
    returnstatus: Integer
        This is optional. Default is -1, no return status. Pass 0 to request the return status
    dbname: String
        This is an option to pass the datakey, dataname combination. When empty, this is not used

        
    Returns
    -------
    Requested data item
        Usually a floating point number or a string
    Status of the command execution (optional)
        An integer that matches LTReturncodeENUMs
    
    Examples
    --------
    #Call without a return status
        L = GetLTDbItem('Solid[1].Primitive[1].Length')
    #Call with the return status
        L,Stat = GetLTDbItem('Solid[1].Primitive[1].Length',0)
    """
    try:
        if dbname == '':
            key,dataname=GetKeyAndDataNameFromString(DataAccessString)
        else:
            dataname=dbname
            key=DataAccessString
        pass
        if ColIndex != -1:
            #2D grid
            [ltdata,stat]=lt0.DbGet(key,dataname,0,ColIndex,RowIndex)
        else:
            #1D grid    
            [ltdata,stat]=lt0.DbGet(key,dataname,0,RowIndex)
        
        if returnstatus==-1:
            return ltdata
        else:
            return ltdata,stat
    except Exception as e:
        print('GetLTGridItem: ' + str(e))
pass
##Local test
#Call without a return status
#L = GetLTDbItem('Solid[1].Primitive[1].Length)
#Call with the return status
#L,Stat = GetLTDbItem('Solid[1].Primitive[1].Length,0)
pass

def SetLTGridItem(DataAccessString,datavalue,RowIndex=1,ColIndex=-1,dbname=''):
    """Get a data value from a 1D or 2D grid item value
    
    Parameters
    ----------
    DataAccessString: String
        This is the string you get via Copy Data Access Name in LightTools
        Note: Content in square brackets at the end (e.g. [4] or [4][6]) is ignored)
    datavalue: string or numeric
        This is the new data value assigned to the database item
    RowIndex: Integer
        Index of the row in the grid, starting 1
    ColIndex: Integer
        Index of the column in the grid, starting 1
    dbname: String
        This is an option to pass the datakey, dataname combination. When empty, this is not used

        
    Returns
    -------
    Status of the command execution (optional)
        An integer that matches LTReturncodeENUMs
    
    Examples
    --------
    #Call without a return status
        L = GetLTDbItem('Solid[1].Primitive[1].Length')
    #Call with the return status
        L,Stat = GetLTDbItem('Solid[1].Primitive[1].Length',0)
    """
    try:
        if dbname == '':
            key,dataname=GetKeyAndDataNameFromString(DataAccessString)
        else:
            dataname=dbname
            key=DataAccessString
        pass
        if ColIndex != -1:
            #2D grid
            stat=lt0.DbSet(key,dataname,datavalue,ColIndex,RowIndex)
        else:
            #1D grid    
            stat=lt0.DbSet(key,dataname,datavalue,RowIndex)
        
        return stat
    except Exception as e:
        print('SetLTGridItem: ' + str(e))
pass
##Local test
#Call without a return status
#L = GetLTDbItem('Solid[1].Primitive[1].Length)
#Call with the return status
#L,Stat = GetLTDbItem('Solid[1].Primitive[1].Length,0)
pass



def GetLTMeshParams(DataAccessKey,CellValueType,paramsonly=False):
    """Get the data from a receiver mesh.
    
    Parameters
    ----------
    DataAccessKey: String
        Data access string for the receiver mesh, typically obtaned by 'Copy Data Access Name'
    CellValueType: String
        Data type to retrieve - e.g. 'CellValue', 'Flux', etc.
    paramsonly: Boolean
        If this is true then the mesh data array is not returned/retrieved, but other parameters such as xdim, ydim will be returned
        
    Returns
    -------
    X_Dimension
        Number of bins in X dimension
    Y_Dimension
        Number of bins in Y dimension
    Min_X_Bound
        Minimum X bound for the mesh
    Max_X_Bound
        Maximum X bound for the mesh
    Min_Y_Bound
        Minimum Y bound for the mesh
    Max_Y_Bound
        Maximum Y bound for the mesh
    Mesh_Data_Array
        An array of data, based on the cell value type requested
    
    Examples
    --------
    meshkey="receiver[1].Mesh[1]"    
    xdim,ydim,minx,maxx,miny,maxy,md=GetLTMeshParams(meshkey,"CellValue")
    xdim,ydim,minx,maxx,miny,maxy=GetLTMeshParams(meshkey,paramsonly=True)
    """
    try: 
        if DataAccessKey=='':
            print('Invalid Data Access String')
            return None
        key=DataAccessKey
        #Check the mesh key
        if str(lt0.DbGet(key,"Name"))=="":
            return None
        else:
            XDim=int(lt0.DbGet(key,"X_Dimension"))
            YDim=int(lt0.DbGet(key,"Y_Dimension"))
            MinX=lt0.DbGet(key,"Min_X_Bound")
            MaxX=lt0.DbGet(key,"Max_X_Bound")
            MinY=lt0.DbGet(key,"Min_Y_Bound")
            MaxY=lt0.DbGet(key,"Max_Y_Bound")
            if paramsonly==True:
                return XDim,YDim,MinX,MaxX,MinY,MaxY
            else:
                dblArray=System.Array.CreateInstance(System.Double,XDim,YDim)
                [Stat,mData]=lt0.GetMeshData(key,dblArray,CellValueType)
                MeshData=np.ones((XDim,YDim))
                print(XDim,YDim)
                for i in range(0,XDim):
                    for j in range(0,YDim):
                        MeshData[i,j]=mData[i,j]
                        #print(mData[i,j])
                MeshData=np.rot90(MeshData)
                return XDim,YDim,MinX,MaxX,MinY,MaxY,MeshData
            pass
    except Exception as e:
        print('GetLTMeshParams: ' + str(e))
pass

def PlotRaster(DataAccessString,CellValueType,colormap='jet',xlabel='X',ylabel='Y',zlabel='Value',title='',plottype='2D',plotsize=(5,5),returndata=False):
    """Creates a 2D or a 3D plot for a given receiver mesh. Optionally, data can be returned
    
    Parameters
    ----------
    DataAccessKey: String
        Data access string for the receiver mesh, typically obtaned by 'Copy Data Access Name'
    CellValueType: String
        Data type to retrieve - e.g. 'CellValue', 'Flux', etc.
    Chart_Parameters: Miscellaneous types
        See the Python/matplotlib documentation for charting for more details
    
    Returns
    -------
    X_Dimension
        Number of bins in X dimension
    Y_Dimension
        Number of bins in Y dimension
    Min_X_Bound
        Minimum X bound for the mesh
    Max_X_Bound
        Maximum X bound for the mesh
    Min_Y_Bound
        Minimum Y bound for the mesh
    Max_Y_Bound
        Maximum Y bound for the mesh
    Mesh_Data_Array
        An array of data, based on the cell value type requested
    
    Examples
    --------
    meshkey="receiver[1].Mesh[1]"    
    xdim,ydim,minx,maxx,miny,maxy,md=GetLTMeshParams(meshkey,"CellValue")
    
    """
    
    try:
        #This is the default width for the plot
        figw=plotsize[0]
        figh=plotsize[1]
        ar=1
        #for data with a single column, we need to add at least one more column
        xdim,ydim,minx,maxx,miny,maxy,d=GetLTMeshParams(DataAccessString,CellValueType)
        if d.shape[1]<2:
            nrows=d.shape[0]
            d2=np.zeros((nrows,2))
            d2[:,0]=d[:,0]
            d2[:,1]=d[:,0]
            d=d2
        #check the x and y scales to maintain the aspect
        w=abs(maxx-minx)
        h=abs(maxy-miny)
        if w>h:
            ar=h/w
            figh=figh*ar
        else:
            ar=w/h
            figw=figw*ar
        if plottype[0]=='2':
            cellx=np.linspace(minx,maxx,xdim+1)
            celly=np.linspace(miny,maxy,ydim+1)
            X,Y=np.meshgrid(cellx,celly)
            #plt.figure(figsize=(int(figw),int(figh)))
            plt.figure(figsize=plotsize)
            plt.pcolormesh(X,Y,np.flipud(d),cmap=colormap)
            plt.axis('equal')
            #plt.axis('tight')
            plt.xlim(minx,maxx)
            plt.ylim(miny,maxy)
            plt.xlabel(xlabel)
            plt.ylabel(ylabel)
            plt.colorbar()
            if title != '':
                plt.title(title)
        else: 
            cellx=np.linspace(minx,maxx,xdim)
            celly=np.linspace(miny,maxy,ydim)
            X,Y=np.meshgrid(cellx,celly)
            fig=plt.figure(figsize=(int(figw),int(figh)))
            ax = fig.add_subplot(111, projection='3d')
            ax.plot_surface(X, Y, d,cmap=colormap)
            #ax.plot_trisurf(X, Y, d,cmap=colormap,linewidth=0.2)
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
            ax.set_zlabel(zlabel)
            if title != '':
                ax.set_title(title)
#            for angle in range(0, 360):
#                ax.view_init(30, angle)
#                plt.draw()
        if returndata==True:
            return xdim,ydim,minx,maxx,miny,maxy,d
        else:
            return 0
   
    except Exception as e:
        print('PlotRaster: ' + str(e))
pass

def PlotSpectralDistribution(DataAccessKey,returndata=True):
    """Plots the spectral distribution on a receiver and returns the data
    
    Parameters
    ----------
    DataAccessKey: String
        This is the string you get via Copy Data Access Name in LightTools for the spectral distribution
    returndata: Boolean
        This is optional. Default is True, meaning the data will be returned. 
    
    Returns
    -------
    Row Count
        Number of wavelength/power pairs in the spectral distribution
    Data array
        Array containing the wavelength and the relative power
    
    Examples
    --------
    numrows,spd=PlotSpectralDistribution('receiver[2].spectral_distribution[1]')
    plt.plot(spd[:,0],spd[:,1])
        
    """
    
    try:
        rowcount=int(lt0.DbGet(DataAccessKey,'Count'))
        sd=np.zeros((rowcount,2))
        for i in range(1,rowcount+1):
            sd[i-1,0],stat=GetLTGridItem(DataAccessKey + '.Wavelength_At',RowIndex=i,returnstatus=0)
            sd[i-1,1],stat=GetLTGridItem(DataAccessKey + '.Power_At',RowIndex=i,returnstatus=0)
        plt.plot(sd[:,0],sd[:,1])
        plt.xlabel('Wavelength (nm)')
        plt.ylabel('Relative Power')
        ax = plt.gca()
        ax.grid(True)
        if returndata==True:
            return rowcount,sd
        else:
            return 0
            
    except Exception as e:
        print('PlotSpectralDistribution: ' + str(e))

    
def PlotTrueColorRster(DataAccessKey,xlabel='X',ylabel='Y',title='',plotsize=(5,5),returndata=False):
    """Creates an RGB image from R,G,B data on the CIE data mesh
    
    Parameters
    ----------
    DataAccessKey: String
        This is the string you get via Copy Data Access Name in LightTools for the CIE mesh (spatial or angular)
    returndata: Boolean
        This is optional. Default is False (no data is returned) 
    
    Returns
    -------
    Data arrays, R, G, and B. Note, for 'imshow', data must be of type 'numpy.uint8'
    
    Examples
    --------
    r,g,b=PlotTrueColorRster('receiver[1].mesh[2]',returndata=True)
        
    """
    
    try:
        xdim,ydim,minx,maxx,miny,maxy,d=GetLTMeshParams(DataAccessKey,'Red_Value_UI')
        R=d
        xdim,ydim,minx,maxx,miny,maxy,d=GetLTMeshParams(DataAccessKey,'Green_Value_UI')
        G=d
        xdim,ydim,minx,maxx,miny,maxy,d=GetLTMeshParams(DataAccessKey,'Blue_Value_UI')
        B=d
        im=np.zeros((xdim,ydim,3),dtype=np.uint8)
        im[:,:,0]=R
        im[:,:,1]=G
        im[:,:,2]=B
        plt.figure(figsize=plotsize)
        plt.imshow(im,origin='lower', extent=[minx, maxx, miny, maxy],interpolation='nearest')
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title=title
        if returndata==True:
            return R,G,B
        else:
            return 0
            
    except Exception as e:
         print('PlotTrueColorRster: ' + str(e))
pass

def GetViewImage(viewname):
    """Capture a screenshot of a given view (3D, 2D or Chart)
    
    Parameters
    ----------
    viewname: String
        This is the name of the view to capture. Can use the short name '3D' for the 3D View. Chart views have specific names.
        
    Returns
    -------
    File Name
        Path and the file name of the image saved in the working directory
    Image
        Image (returned from misc.imread() function)
        
    
    Examples
    --------
    viewname='3d'
    im,imname=GetViewImage(viewname)
    plt.imshow(im)
    
    """
    
    import System.Windows.Forms
    from System.Windows.Forms import Clipboard
    import tempfile
    from scipy import misc
    #import LTUtilities as ltu
    try:
        workdirstr=ltu.checkpyWorkDir()
        viewname=viewname.upper()
        if viewname[:2]=='3D':
            cmdstr=chr(92) + 'V3D'
        else:
            vname=lt0.ViewGet('View[' + viewname + ']','Name')
            vname=vname.replace(' ','_')
            cmdstr= chr(92) +'V' + 'Chart_' + vname
        lt0.Cmd(cmdstr)
        
        Clipboard.Clear()
        lt0.Cmd('CopyToClipboard')
        mybmap=Clipboard.GetImage()
        tempfname = next(tempfile._get_candidate_names())
        #print(tempfname)
        fname=workdirstr + '/' + tempfname + '.png'
        #print(fname)
        mybmap.Save(fname,System.Drawing.Imaging.ImageFormat.Png)
        del(mybmap)
        im=misc.imread(fname)
        return im,fname
    except Exception as e:
         print('GetViewImage: ' + str(e))
    
pass


def GetLTReceiverRays(DataAccessKey,descriptors=['raydatax','raydatay','raydataz'],simdata='Forward_Sim_Function[1]',usepassfilters=False):
    """Retrieves the ray data on the receiver
    
    Parameters
    ----------
    DataAccessKey : String
        Data access key for the receiver, typically obtaned by 'Copy Data Access Name' - e.g. 'receiver[1]'
    descriptors: String
        ray data items to retrieve. This should be a list containing all the desired data items.  
        These are the available descriptors (as of 10/2016):
        "RAYDATAX", "RAYDATAY", "RAYDATAZ", "RAYDATAENTERINGL", "RAYDATAENTERINGM", "RAYDATAENTERINGN",
        "RAYDATAMAGNITUDE", "RAYDATAWAVELENGTH", "RAYDATAABSORPTION",
        "RAYDATAORDINALNUMBER" ----> If you trace N rays, the ordinal number is between 1 and N,
        "RAYDATAPASSFILTERS" ---->  1, if the ray passes all filters; otherwise, 0,
        "RAYPATHINDEX" ---->  Base 1 index of the ray pathes (e.g., 1 would correspond to PathIndex=1 in the LightTools ray path),
        "STOKES0", "STOKES1", "STOKES2", "STOKES3", "POLDIRX", "POLDIRY", "POLDIRZ", "RAYTRANSMITTANCE"  
    simdata: String
        The type of simulation to use; the default is forward simulation (Forward_Sim_Function)
    usepassfilters: Boolean
        If False, all the ray data is returned. If True, only rays that pass receiver filter criteria will be returned
    
    Returns
    -------
    Number of Rays (N)
        Number of rays retrieved
    Number of Data Descriptors (M)
        Number of data descriptors received in the list
    Ray Data
        This is an array (N rows, M columns) with requested ray data
    
    Examples
    --------
    reckey="receiver[1]" --- assume forward simulation (default)
    descriptors=['RayDataX','RayDataY','RayDataZ','RayDataWavelength']
    N,M,RayData = GetReceiverRays(reckey,descriptors)
    
    
    """
    try:
        if len(descriptors)<1:
            return -1
        
        key=DataAccessKey + '.' + simdata
        numrays=int(lt0.DbGet(key,'NumberOfSamples'))
        if usepassfilters==True:
            descriptors.append('raydatapassfilters')
        pass
        numdes=len(descriptors)
        dblArray=System.Array.CreateInstance(System.Double,numrays,numdes)
        strArray=System.Array.CreateInstance(System.String,numdes,1)
        for i in range(0,numdes):
            strArray[i,0]=descriptors[i]
        pass
        stat,d1,rdata=lt0.GetReceiverRayData(key,strArray,dblArray)
        #If usepassfilters=True, we don't know how many qualify
        #Resizing arrays dynamically is inefficient, so first check the number of rays that qualify
        numpassedrays=0
        if usepassfilters==True:
            for i in range(0,numrays):
                passflag=int(rdata[i,numdes-1])
                if passflag>0:
                    numpassedrays += 1
                pass
            pass
            raydata=np.zeros((numpassedrays,numdes-1)) #We added passfilters flag at the end, so ignore it
            print('Qualified number of rays: ' + str(numpassedrays))
            numpassedrays=0
            for i in range(0,numrays):
                passflag=int(rdata[i,numdes-1])
                if passflag>0:
                    for j in range(0,numdes-1):
                        raydata[numpassedrays,j]=rdata[i,j]
                    pass
                    numpassedrays += 1
                pass
            pass
            #We need to remove the appended passfilters flag
            descriptors.remove('raydatapassfilters')
            return numpassedrays,numdes-1,raydata
        else:
            raydata=np.zeros((numrays,numdes))
            for i in range(0,numrays):
                for j in range(0,numdes):
                    raydata[i,j]=rdata[i,j]
                pass
            pass
            return numrays,numdes,raydata
    except Exception as e:
         print('GetLTReceiverRays: ' + str(e))
pass

def GetLTReceiverRays_Extra(DataAccessKey,FilterIndex,simdata='Forward_Sim_Function[1]',StartRay=-1,EndRay=-1):
    
    """Retrieves the extra ray data items, such as Optical Path length, that is not available with LTAPIx.GetReceiverRayData() function.
    This is significantly slower than the GetReceiverRays(), when there is a large number of rays on the receiver!
    
    Parameters
    ----------
    DataAccessString: String
        Data access key for the receiver, typically obtaned by 'Copy Data Access Name' - e.g. 'receiver[1]'
    FilterIndex: Integer
        This is the type of data to retrieve. Use the ENUM values at the top of this module to find the type of data. Value must be within [1, 21].
        For some data types, a receiver filter must be present. Refer to the documentation for details
    StartRay: Integer
        Default (-1) will start from ray number 1
    EndRay: Integer
        Default (-1) will use all rays, from the start ray.
        
    Returns
    -------
    Number of rays retrieved.
    
    Requested data as an array of floats
    
    Examples
    --------
    reckey="receiver[1]"
    N,exdata=GetReceiverRays_Extra(reckey,ExtraRayData.Optical_Path_Length.value)
    
    """
    try:  
        if FilterIndex>21:
            return -1
        
        key=DataAccessKey + '.' + simdata
        numrays=int(lt0.DbGet(key,'NumberOfSamples'))
        if numrays<1:
            return -1
        
        if StartRay>numrays:
            print('Start ray must be < number of saved rays on receiver')
            return -1
            
        if StartRay>EndRay:
            EndRay=StartRay+1
            print('End ray must be > Start Ray')
        
        if (StartRay==-1 and EndRay==-1):
            ExData=np.zeros((numrays))
            
            for i in range(1,numrays+1):
                ExData[i-1]=GetLTGridItem(key + '.ExtraRayData',i,FilterIndex)
                if i>numrays:
                    break
            return numrays,ExData
        else:
            if StartRay<0:
                StartRay=1
            if EndRay<0:
                EndRay=numrays
            n=0
            for i in range(StartRay,EndRay+1):
                ExData[n]=GetLTGridItem(key + '.ExtraRayData',i,FilterIndex)
                if i>numrays:
                    break
                n += 1
            return n,ExData
                
    #plt.hist(ExData,bins=21)
    except Exception as e:
         print('GetLTReceiverRays_Extra: ' + str(e))
pass

def GetRayPathData(DataAccessKey,simdata='forward_sim_function[1]', usevisibleonly=False):
    """Get the ray path data from a given receiver
    
    Parameters
    ----------
    DataAccessString: String
        Data access key for the ray paths, obtaned by 'Copy Data Access Name' - e.g. 'receiver[1]'
    usevisiblecolumns: Boolean
        Pass True in order to retrieve data only for the visible paths
        
    Returns
    -------
    The following data items from the grid are returned
        RayPathVisibleAt, RayPathPowerAt, RayPathNumRaysAt, RayPathStringAt
        
    Examples
    --------
    visibleat,powerat,raysat,stringat=GetRayPathData('receiver[1].forward_sim_function)
    
    """
    try:
        key=DataAccessKey + '.' + simdata
        numpaths=int(GetLTDbItem(key + '.NumberOfRayPaths'))
        if numpaths<1:
            return -1
        
        visibleat=[]
        powerat=[]
        numraysat=[]
        stringat=[]
        vakey=key + '.RayPathVisibleAt'
        pakey=key + '.RayPathPowerAt'
        rakey=key + '.RayPathNumRaysAt'
        sakey=key + '.RayPathStringAt'
    
        for i in range(1,numpaths+1):
            va=GetLTGridItem(vakey,i)
            if usevisibleonly==True:
                va=va.lower()
                if va[0]=='y':
                    visibleat.append(va)
                    pa=GetLTGridItem(pakey,i)
                    ra=GetLTGridItem(rakey,i)
                    sa=GetLTGridItem(sakey,i)
                    visibleat.append(va)
                    powerat.append(pa)
                    numraysat.append(ra)
                    stringat.append(sa)
            else:
                visibleat.append(va)
                pa=GetLTGridItem(pakey,i)
                ra=GetLTGridItem(rakey,i)
                sa=GetLTGridItem(sakey,i)
                powerat.append(pa)
                numraysat.append(ra)
                stringat.append(sa)
    
        return visibleat,powerat,numraysat,stringat
    except Exception as e:
        print('GetRayPathData: ' + str(e))
pass

def MakeRayFileUsingRayOrdinal(DataAccessKey,descriptors=['raydataX','raydataY','raydataZ','raydataenteringL','raydataenteringM','raydataenteringN','raydataMAGNITUDE'],simdata='Forward_Sim_Function[1]',DataAccessKey_Ordinal=''):
    """Creates a ray data file from receiver rays, based on ray ordinal numers from another receiver
    
    Parameters
    ----------
    DataAccessKey : String
        Data access key for the receiver, typically obtaned by 'Copy Data Access Name' - e.g. 'receiver[1]'
    descriptors: String
        ray data items to retrieve. This should be a list containing all the desired data items.  
        These are the available descriptors (as of 10/2016):
        "RAYDATAX", "RAYDATAY", "RAYDATAZ", "RAYDATAENTERINGL", "RAYDATAENTERINGM", "RAYDATAENTERINGN",
        "RAYDATAMAGNITUDE", "RAYDATAWAVELENGTH", "RAYDATAABSORPTION",
        "RAYDATAORDINALNUMBER" ----> If you trace N rays, the ordinal number is between 1 and N,
        "RAYDATAPASSFILTERS" ---->  1, if the ray passes all filters; otherwise, 0,
        "RAYPATHINDEX" ---->  Base 1 index of the ray pathes (e.g., 1 would correspond to PathIndex=1 in the LightTools ray path),
        "STOKES0", "STOKES1", "STOKES2", "STOKES3", "POLDIRX", "POLDIRY", "POLDIRZ", "RAYTRANSMITTANCE"  
    simdata: String
        The type of simulation to use; the default is forward simulation (Forward_Sim_Function)
    DataAccessKey_Ordinal: String
        This is the data access key for the receiver used to retrieve ray ordinal numbers
    
    Returns
    -------
    Number of Rays (N)
        Number of rays retrieved
    File Name
        This is the name of the new ray data file
        
    Examples
    --------
    reckey1='receiver[1]' --- this is the receiver for ray data, assume forward simulation (default)
    
    reckey2='receiver[2]'
    
    descriptors=['RayDataX','RayDataY','RayDataZ','RayDataWavelength']
    
    N,FName = MakeRayFileUsingRayOrdinal(reckey1,descriptors,DataAccessKey_Ordinal=reckey2)
    

    """
    
    try:
        import tempfile
        descriptors.append('RayDataOrdinalNumber')
        n1,m1,rd1=GetLTReceiverRays(DataAccessKey,descriptors,simdata)
        des=['raydataordinalnumber'] #This is all we need from this receiver
        n2,m2,rd2=GetLTReceiverRays(DataAccessKey_Ordinal,des,simdata,usepassfilters=True)
        tempfname = next(tempfile._get_candidate_names())
        #print(tempfname)
        workdirstr=ltu.checkpyWorkDir()
        fname=workdirstr + '/' + tempfname + '.txt'
        print('Using temporary file name: ' + fname)
        rf=open(fname,'w')
        rf.write(ltu.GetRayFileHeader())
        N=0
        for i in range(0,n1):
            tstr=''
            if rd1[i,m1-1] in rd2:
                N += 1
                for j in range(0,m1-1): #we added the ordinal, so don't write that
                    tstr=tstr + ' ' + str(rd1[i,j])
                pass
                rf.write(tstr + '\n\r')
            pass
        pass
        rf.write(ltu.GetRayFileFooter())
        rf.close()
        return N,fname
    except Exception as e:
        print('GetReceiverRaysByOrdinal: ' + str(e))
pass

def GetAOIDistributionFromReceiverRays(DataAccessKey,simdata='forward_sim_function[1]',usepassfilters=True,referencevector=(0,0,1)):
    """Returns the angle of incidence relative to a direction vector. Default is surface normal [0,0,1]
    
    Parameters
    ----------
    DataAccessKey : String
        Data access key for the receiver, typically obtaned by 'Copy Data Access Name' - e.g. 'receiver[1]'
    simdata: String
        The type of simulation to use; the default is forward simulation (Forward_Sim_Function)
    usepassfilters: Boolean
        This flag specifies whether to use any filters enabled on the receiver for data    
    referencevector: 3-element array like 0,0,1
        The default is the surface normal - 0,0,1
        
    Returns
    -------
    An array with computed incident angles
        
    Examples
    --------
    aoi=ltd.GetAOIDistributionFromReceiverRays('receiver[1]',referencevector=(0,0,1))
    plt.hist(aoi*180/math.pi,bins=30,range=(0,90))
    
    """
    
    try:
        des=['raydataEnteringL','raydataEnteringM','raydataEnteringN']
        n1,m1,rd1=GetLTReceiverRays(DataAccessKey,descriptors=des,simdata=simdata,usepassfilters=usepassfilters)
        ang=np.zeros(n1)
        #v2=(0,0,1) #surface normal
        v2=referencevector
        for i in range(0,n1):
            v1=(rd1[i,0],rd1[i,1],rd1[i,2])
            ang[i]=angle_between_vectors(v1,v2)
        
        return ang
    except Exception as e:
        print('GetAOIDistributionFromReceiverRays: ' + str(e))
pass

def ComputeLumDataForMesh(ReceiverDataAccessKey,carea,flux,datatype):
    """
    
    """
    try:
        print('Analysis type: ' + str(datatype))
        datatype=int(datatype)
        if carea.any() != 0:
            E=flux/carea #For intensity, CellArea is the solid angle
            if (datatype==0 or datatype==1): #Illum, Intensity
                return E
            elif datatype==2: #Spatial Lum
                meterkey=ReceiverDataAccessKey + '.SPATIAL_LUM_METER[1].CentralPSA'
                psa=GetLTDbItem(meterkey)
                if psa != 0:
                    E=E/psa
                else:
                    print('Central PSA must be > 0')
                    return None
                pass
            elif datatype==3: #Angular Lum
                aperturekey=ReceiverDataAccessKey + '.ANGULAR_LUM_METER[1].HalfSize'
                aperturesize=2*float(GetLTDbItem(aperturekey))
                if aperturesize != 0:
                    E=E/(aperturesize*aperturesize)
                else:
                    print('Aperture size must be > 0')
                    return None
                pass
            else:
                print('Unexpected data type.')
                return None
            pass
            return E
        else:
            print('Cell Area must be > 0')
            return None
    except Exception as e:
        print('ComputeLumDataForMesh: ' + str(e))
pass

def RebinReceiverRayData(ReceiverDataAccessKey,MeshDataAccessKey,simdata='forward_sim_function[1]',useraypathindex=False):
    """Only illuminance data is supported!
    """
    try:
        #Smoothing can produce different results in the UI
        #These calculations do not take smoothing into account
        smoothing=GetLTDbItem(MeshDataAccessKey + '.Do_Noise_Reduction')
        if smoothing != '':
            if smoothing.upper() == 'YES':
                print('Warning: mesh smoothing is ON.')
        else:
            print('Invalid mesh key.')
            return None
        pass
    
        #We need to determine the data type using the mesh key
        analysistype=-1
        datatype=['ILLUMINANCE_MESH', 'INTENSITY_MESH', 'SPATIAL_LUMINANCE_MESH', 'ANGULAR_LUMINANCE_MESH']
        filtertype=GetLTDbItem(MeshDataAccessKey + '.FilterName')
        filtertype=filtertype.upper()
        for i in range(0,4):
            if datatype[i] == filtertype:
                analysistype=i
                break
            pass
        pass
        if analysistype != 0:
            #Only illuminance is supported for now
            print('Data types other than illuminance are currently unsupported.')
            return None
            #We need to calculate some extra data
            #For luminance, filter rays inside the cone angle
            conekey=ReceiverDataAccessKey + '.SPATIAL_LUM_METER[1].CentralConeAngle'
            conehalfangle=GetLTDbItem(conekey)
            rayaoi=GetAOIDistributionFromReceiverRays(ReceiverDataAccessKey,simdata)
            raysincone=np.where(rayaoi<=conehalfangle) #These are the indices
            
            #For intensity, calculate angles
        pass
        
        xdim,ydim,minx,maxx,miny,maxy,carea=GetLTMeshParams(MeshDataAccessKey,'CellSurfaceArea')
        simrays=int(GetLTDbItem('SIMULATIONS[1].CurProgress'))       
        if useraypathindex==True:
            des=['raydatax','raydatay','raydatamagnitude','raypathindex']
            n,m,rd=GetLTReceiverRays(ReceiverDataAccessKey,simdata=simdata,descriptors=des)
            x=rd[:,0]
            y=rd[:,1]
            mag=rd[:,2]
            rpi=rd[:,3]
            if analysistype == 2 or analysistype == 3:
                #we need to clip rays
                x=x[raysincone]
                y=y[raysincone]
                mag=mag[raysincone]
            #A 3D array to save data for each ray path
            numpaths=int(GetLTDbItem(ReceiverDataAccessKey + '.' + simdata + '.NumberOfRayPaths'))
            pdata=np.zeros((ydim,xdim,numpaths)) #Rows/columns need to match, after flip
            for i in range(1,numpaths+1):
                xp=x[rpi==float(i)]
                yp=y[rpi==float(i)]
                magp=mag[rpi==float(i)]
                #sumpowerp, xedges, yedges = np.histogram2d(xp,yp,bins=(xdim,ydim),range=([sumpowerp=np.rot90(sumpowerp) #Same orientation
                sumpower, xedges, yedges = np.histogram2d(xp,yp,bins=(xdim,ydim),range=([minx,maxx],[miny,maxy]),weights=magp)
                sumpower=np.rot90(sumpower) #Same orientation
                fluxp=sumpower/simrays
                Ep=ComputeLumDataForMesh(ReceiverDataAccessKey,carea,fluxp,analysistype)
                pdata[:,:,i-1]=Ep
            pass
            return pdata
        else:
            des=['raydatax','raydatay','raydatamagnitude']
            n,m,rd=GetLTReceiverRays(ReceiverDataAccessKey,simdata=simdata,descriptors=des)
            x=rd[:,0]
            y=rd[:,1]
            mag=rd[:,2]
            sumpower, xedges, yedges = np.histogram2d(x,y,bins=(xdim,ydim),range=([minx,maxx],[miny,maxy]),weights=mag)
            sumpower=np.rot90(sumpower) #Same orientation
            flux=sumpower/simrays
            sumpower
            flux
            E=ComputeLumDataForMesh(ReceiverDataAccessKey,carea,flux,analysistype)
            return E
        pass
    except Exception as e:
        print('RebinReceiverRayData: ' + str(e))
pass

def GetAzimuthElevationFromRayData(DataAccessKey,simdata='forward_sim_function[1]',datarange=[0,360,0,90]):
    """Calculates and returns the azimuth, elevation from raya data. Mostly suitable for surface receivers
    
    Parameters
    ----------
    DataAccessKey: String
        This is the data access key for the receiver
    
    simdata: String
        Simulation data to use. The default is forward simulation
    
    datarange: List
        A list containing the min/max for Azimuth and Elevation. The full range will return for the whole hemisphere
    
    Returns
    -------
    Three arrays, containing Azimuth, Elevation, and Ray Data Magnitude
    
    Examples
    --------
    Get the ray data for azimuth 45,90 and elevation 45,60
        azm,elv,rdm=GetAzimuthElevationFromRayData('receiver[1]',datarange=[45,90,45,60])
    Get the ray data for the hemisphere
        azm,elv,rdm=GetAzimuthElevationFromRayData('receiver[1]') 
    """
    
    try:
        des=['raydataEnteringL','raydataEnteringM','raydataEnteringN','raydatamagnitude']
        n1,m1,rd0=GetLTReceiverRays(DataAccessKey,descriptors=des,simdata=simdata)
        azm=np.zeros(n1)
        elv=np.zeros(n1)
        for i in range(0,n1):
            azm[i],elv[i],r=cart2sph(rd0[i,0],rd0[i,1],rd0[i,2])
        
        azm=azm*180/math.pi
        elv=elv*180/math.pi
        azm[azm<0]=360-abs(azm[azm<0])
        elv=90-elv #We want elevation defined wrt the surface normal
        rdm=rd0[:,3]
        
        #if the datarange covers a hemisphere, then we're set
        if datarange==[0,360,0,90]:
            return azm,elv,rdm
        else:
            #If not, get a 'patch' of angles
            #First, isolate azimuth
            k=np.where((azm>=datarange[0]) & (azm<=datarange[1]))
            az2=azm[k]
            elv2=elv[k]
            rdm2=rdm[k]
            #From filtered azimuth,  isolate elevation
            j=np.where((elv2>=datarange[2]) & (elv2<=datarange[3]))
            az2=az2[j]
            elv2=elv2[j]
            rdm2=rdm2[j]
            return az2,elv2,rdm2
        pass
    except Exception as e:
        print('GetAzimuthElevationFromRayData: ' + str(e))
pass

def PlotCIE1931ForMeshData(DataAccessKey,returndata=False):
    #Plot CIE xy chart
    try:
        #First, construct the CIE diagram
        xyz=np.loadtxt('xyz31.txt',delimiter=',')
        xb=xyz[:,1]
        yb=xyz[:,2]
        zb=xyz[:,3]
        xlocus=xb/(xb+yb+zb)
        ylocus=yb/(xb+yb+zb)
        x0=np.zeros(2)
        y0=np.zeros(2)
        x0[0]=xlocus[0]
        x0[1]=xlocus[len(xlocus)-1]
        y0[0]=ylocus[0]
        y0[1]=ylocus[len(ylocus)-1]
        A=np.concatenate([xlocus,x0])
        B=np.concatenate([ylocus,y0])
        
        #Get the LT mesh data
        xdim,ydim,minx,maxx,miny,maxy,fc=GetLTMeshParams(DataAccessKey,"First_CIE_Coord_Value_UI")
        xdim,ydim,minx,maxx,miny,maxy,sc=GetLTMeshParams(DataAccessKey,"Second_CIE_Coord_Value_UI")
        
        plt.figure(figsize=(8,8))
        plt.plot(A,B)
        plt.hold
        plt.grid(which='major')
        plt.scatter(fc,sc,marker='+')
        axes = plt.gca()
        axes.set_xlim([0,0.8])
        axes.set_ylim([0,0.9])
        plt.axis('equal')
        #Return data based on the request
        if returndata==True:
            return A,B,fc,sc
        else:
            return 0
        pass
    except Exception as e:
        print('PlotCIE1931ForMeshData: ' + str(e))
pass



def GetVarsFromUDV(udvname='1'):
    """Get all variable values from a user variable collection [UDVG]
    """
    try:
        basekey='OPT_USERVARIABLECOLLECTION[' + udvname + ']'
        das=basekey + '.PARAMETER[Number_Of_Variables].ValueStrUI'
        numvar=int(GetLTDbItem(das))
        #print('Number of variables: ' + str(numvar))
        VarVals=np.zeros(numvar)
        for i in range(1,numvar+1):
            vardas=basekey + '.OPT_USERVARIABLE[' + str(i) + '].CurrentValue'
            #print(vardas)     
            currv=float(GetLTDbItem(vardas))       
            VarVals[i-1]=currv
            #print('Variable value is: ' + str(VarVals[i-1]))
        pass
        return VarVals
    except Exception as e:
        print('GetVarsFromUDV: ' + str(e))
pass
        
def ReadStokesFromASCIIRayFile(rayfile):
    #read the ray data file
    try:
        waveinfo=False
        startreaddata=False
        poldata=[]
        xyz=[]
        fname=rayfile
        rf=open(fname,'r')
        for line in rf:
            if line != '':
                line=line.upper()
                if 'LT_ENDOFDATA' in line:
                    startreaddata=True
                pass
                if 'LT_COLOR_INFO:' in line:
                    print(line)
                    if 'WAVELENGTH' in line:
                        waveinfo=True
                pass
                if 'LT_STARTOFDATA' in line:
                    startreaddata=True
                pass
                if startreaddata==True and ('LT' in line)==False:
                    line=line.replace(chr(9),' ')
                    s=line.split(' ')
                    s=filter(None,s)
                    s=list(s)
                    xyzstr=s[0:3]
                    nd_xyz=[float(i) for i in xyzstr]
                    xyz.append(nd_xyz)
                    if waveinfo==True:
                        s=s[8:12]
                    else:
                        s=s[7:11]
                    nd=[float(i) for i in s]
                    poldata.append(nd)
                pass
        pass
        rf.close
        poldata=np.asarray(poldata)
        xyz=np.asarray(xyz)
        #xdim,ydim,minx,maxx,miny,maxy=GetLTMeshParams(meshkey,'cellvalue',paramsonly=True)
        x=xyz[:,0]
        y=xyz[:,1]
        S0=poldata[:,0]
        S1=poldata[:,1]
        S2=poldata[:,2]
        S3=poldata[:,3]
        #LT ray data can contain numbers like -xxE-17
        thresholdvalue=1e-3
        S0[np.where(np.abs(S0)<thresholdvalue)]=0
        S1[np.where(np.abs(S1)<thresholdvalue)]=0
        S2[np.where(np.abs(S2)<thresholdvalue)]=0
        S3[np.where(np.abs(S3)<thresholdvalue)]=0
        return x,y,S0,S1,S2,S3
    except Exception as e:
         print('ReadStokesFromASCIIRayFile: ' + str(e))
pass