# -*- coding: utf-8 -*-
"""
Created on Thu Oct 13 16:41:49 2016

@author: mali
"""

import os
import clr
import inspect
import System

lt0=None

def ConvertArrayToNET(numpyarray,numrows,numcols):
    """Convert a numpy array to .NET array
    
    Parameters
    ----------
    numpyarray: 1D or 2D numpy array
    
    numrows: Number of rows
    
    numcols: Number of columns. If numcols==0, 1D array is assumed
    
    Returns
    -------
        Converted array in .NET format
        
    Examples
    --------
        myNET2darray=ConvertArrayToNET(mynumpyarray,4,5)
        myNET1darray=ConvertArrayToNET(mynumpyarray,4,0)
    
    """
    try:
        if numcols==0:
            #1D array
            dblArray=System.Array.CreateInstance(System.Double,numrows)
            for i in range(0,numrows):
                dblArray[i]=numpyarray[i]
            pass
            return dblArray
        else:
            dblArray=System.Array.CreateInstance(System.Double,numrows,numcols)
            for i in range(0,numrows):
                for j in range(0,numcols):
                    dblArray[i,j]=numpyarray[i,j]
                pass
            pass
            return dblArray
        pass
    except Exception as e:
        print('ConvertArrayToNET: ' + str(e))
    pass
pass

def GetRayFileFooter():
    """Returns the footer string for a ray file, without color info
    """   
    return 'lt_endofdata'    
pass


def GetRayFileHeader():
    """Returns the header string for a ray file, without color info
    """   
    hstr=''
    hstr=hstr + 'lt_rdf_version: 2.0' + '\n\r'
    hstr=hstr + 'dataname: untitledRayData' + '\n\r' 
    hstr=hstr + 'lt_datatype: radiant_power' + '\n\r' 
    hstr=hstr + 'lt_radiant_flux:  1' + '\n\r' 
    hstr=hstr + 'lt_far_field_data: NO' + '\n\r' 
    hstr=hstr + 'lt_color_info: none' + '\n\r' 
    hstr=hstr + 'lt_length_units: millimeters' + '\n\r' 
    hstr=hstr + 'lt_data_origin:  0  0  0' + '\n\r' 
    hstr=hstr + 'lt_startofdata' + '\n\r'
    return hstr
    
pass

def HideFileDialog():
    """Disables the open/save file dialog for operations like opening/importing a file
    """
    lt0.SetOption('ShowFileDialogBox',0)
pass

def SetViewFocus(viewname = '3D'):
    """Set the focus to a specific view. Default is the 3D View
    
    Parameters
    ----------
    viewname: String
        Name of the view that needs the focus
    
    Returns
    -------
        None
        
    Examples
    --------
        SetViewFocus('3D')
    
    """
    cmdstr=chr(92) + 'V' +viewname
    lt0.Cmd(cmdstr)

def InitJS(lthome='',LTPIDx=0):
    """Connect to a specific running LightTools session using the Process ID
    
    Parameters
    ----------
    lthome: LightTools installation folder
    
    LTPIDx: Integer/Long, obtained either from Windows Task Manager or LightTools Console Window
        
    Returns
    -------
        None    
    
    Examples
    --------
    Connect to the LightTools session with process ID 3096
        InitLT('C:\\Program Files\\Optical Research Associates\\LightTools 8.4.0', 3096)
    """
    try:
        ltcom64dll=lthome + '\\Utilities.NET\\LTCOM64.dll'
        if os.path.exists(ltcom64dll)==False:
            print('Specified path ' + ltcom64dll + ' is invalid for LTCOM64.dll.')
            return None
        else:
            clr.AddReference(ltcom64dll)
            from LTCOM64 import JSNET2
            js=JSNET2()
            js.LTPID=LTPIDx
            js.UpdateLTPointer
            print('LT pointer updated.')
            return js
    except Exception as e:
        print('InitJS(): ' + str(e))
    
def InitLT(lthome='',LTPIDx=0):
    """Connect to a specific running LightTools session using the Process ID
    
    Parameters
    ----------
    lthome: LightTools installation folder
    
    LTPIDx: Integer/Long, obtained either from Windows Task Manager or LightTools Console Window
        
    Returns
    -------
        None    
    
    Examples
    --------
    Connect to the LightTools session with process ID 3096
        InitLT('C:\\Program Files\\Optical Research Associates\\LightTools 8.4.0', 3096)
    """
    try:
        ltcom64dll=lthome + '\\Utilities.NET\\LTCOM64.dll'
        if os.path.exists(ltcom64dll)==False:
            print('Specified path ' + ltcom64dll + ' is invalid for LTCOM64.dll.')
            return None
        else:
            clr.AddReference(ltcom64dll)
            from LTCOM64 import LTAPIx
            lt=LTAPIx()
            lt.LTPID=LTPIDx
            lt.UpdateLTPointer
            print('LT pointer updated.')
            return lt
    except Exception as e:
        print('InitLT(): ' + str(e))
##Local test
#n=9484
#lthome='C:\\Program Files\\Optical Research Associates\\LightTools 8.4.0'
#InitLT(n)
#lt[0].Message('Hello to ' + str(n))


def PrintModulePath(modulehandle):
    """Print the full path for a given module
    
    Parameters
    ----------
    modulehandle: handle you create for a given module
    
    Returns
    -------
        Path string where the module file is located   
    
    Examples
    --------
        import LTData as ltd
        import LTUtilities as ltu
        ltu.PrintModulePath(ltd)
    """
    print(os.path.dirname(inspect.getfile(modulehandle)))
pass
   
def checkpyWorkDir(workdir=''):
    workdirstr='pyWorkDir'
    ltuser=str(lt0.DbGet('Lens_Manager[1]','LTUserDir'))
    ltuser=ltuser.replace(chr(92),'/')
    print(ltuser)
    if workdir != '':
        workdirstr=workdir
    else:
        workdirstr=ltuser + workdirstr
    if os.path.exists(workdirstr)==False:
        os.makedirs(workdirstr)
    workdirstr=workdirstr
    return workdirstr
pass

def TurnDbUpdateOff():
    lt0.SetOption('DbUpdate',0)
pass

def TurnDbUpdateOn():
    lt0.SetOption('DbUpdate',1)
pass

def GetLTUserDir():
    return str(lt0.DbGet('Lens_Manager[1]','LTUserDir'))