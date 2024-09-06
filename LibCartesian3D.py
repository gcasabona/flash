import os
import math
from yt.mods import *
import matplotlib.pyplot as plt
import numpy as np

@derived_field(name='CylR')
def CylR(field,data):
	return (data['x']**2.+data['y']**2.)**0.5

@derived_field(name='CylPhi')
def CylPhi(field,data):
	return np.arctan2(data['y'],data['x'])

@derived_field(name='CylZ')
def CylZ(field,data):
	return data['z']

@derived_field(name='SphR')
def SphR(field,data):
	return (data['x']**2.+data['y']**2.+data['z']**2.)**0.5

@derived_field(name='CylVelR')
def CylVelR(field,data):
	return data['velx']*np.cos(data['CylPhi'])+data['vely']*np.sin(data['CylPhi'])

@derived_field(name='CylVelPhi')
def CylVelPhi(field,data):
	return -data['velx']*np.sin(data['CylPhi'])+data['vely']*np.cos(data['CylPhi'])

@derived_field(name='CylVelZ')
def CylVelZ(field,data):
	return data['velz']

@derived_field(name='CylMagR')
def CylMagR(field,data):
        return data['magx']*np.cos(data['CylPhi'])+data['magy']*np.sin(data['CylPhi'])

@derived_field(name='CylMagPhi')
def CylMagPhi(field,data):
        return -data['magx']*np.sin(data['CylPhi'])+data['magy']*np.cos(data['CylPhi'])

@derived_field(name='CylMagZ')
def CylMagZ(field,data):
        return data['magz']

@derived_field(name='Beta')
def Beta(field,data):
	return data['pres']/data['magp']

@derived_field(name='OneOverBeta')
def OneOverBeta(field,data):
	return data['magp']/data['pres']	

@derived_field(name='CellVolume')
def CellVolume(field, data):
    if data['dx'].size == 1:
        try:
            return data['dx'] * data['dy'] * data['dz'] * \
                np.ones(data.ActiveDimensions, dtype='float64')
        except AttributeError:
            return data['dx'] * data['dy'] * data['dz']
    return data["dx"] * data["dy"] * data["dz"]

@derived_field(name='CellMass')
def CellMass(field,data):
	return data['CellVolume']*data['dens']

@derived_field(name='RadialMassFlux')
def RadialMassFlux(field,data):
	return data['CellMass']*data['CylVelR']

@derived_field(name="MaxwellStress")  
def MaxwellStress(field,data):
	return -data["CylMagR"]*data["CylMagPhi"] #rahul: CylMagPhi or CylMagZ

@derived_field(name="absMaxwellStress")
def absMaxwellStress(field,data):
        return abs(data["MaxwellStress"])

@derived_field(name='AveragedCylVelPhi')
def AveragedCylVelPhi(field,data):
    VelPhiFloor=0.
    dr = ((data['dx']**2.+data['dy']**2.)**0.5).min()
    dz = data['dz'].min()
    poolR=np.ceil(data['CylR']/dr)
    poolZ=np.ceil(data['CylZ']/dz)
    poolR[poolR==0.]=0.5
    poolZ[poolZ==0.]=0.5

    maxR=(data['CylR']*(abs(data['CylVelPhi'])>VelPhiFloor)).max()
    maxZ=(data['CylZ']*(abs(data['CylVelPhi'])>VelPhiFloor)).max()
    minZ=(data['CylZ']*(abs(data['CylVelPhi'])>VelPhiFloor)).min()

    poolR=poolR*(data['CylR']<=maxR)
    poolZ=poolZ*(data['CylZ']<=maxZ)*(data['CylZ']>=minZ)

    loopR=np.unique(poolR[poolR!=0])
    loopZ=np.unique(poolZ[poolZ!=0])

    PhiAverage=data['CylVelPhi']*0.
    if(loopR.size*loopZ.size!=0.):
        for i in loopR:
            filterR=(poolR==i)
            for j in loopZ:
                filterZ=(poolZ==j)
                filterTotal=filterR*filterZ
                if(filterTotal.sum()!=0.):
                    PhiAverage=PhiAverage+filterTotal*((data['CylVelPhi']*filterTotal).sum()/filterTotal.sum())
    return PhiAverage

@derived_field(name='CylR4_Omega2')
def CylR_sqr_Omega(field,data):
    return data['CylR']**4.*data['Omega']**2.

@derived_field(name='diff_x_CylR4_Omega2',validators=[ValidateSpatial(ghost_zones=1,fields=["CylR4_Omega2"])])
def diff_x_CylRSqr_Omega(field,data):
    if data.pf["HydroMethod"] == 2:
        sl_left = slice(None,-2,None)
        sl_right = slice(1,-1,None)
        div_fac = 1.0
    else:
        sl_left = slice(None,-2,None)
        sl_right = slice(2,None,None)
        div_fac = 2.0
    dx = div_fac * data['dx'].flat[0]
    fx = (data["CylR4_Omega2"][sl_right,1:-1,1:-1]-data["CylR4_Omega2"][sl_left,1:-1,1:-1])/dx
    new_field = np.zeros(data["CylR4_Omega2"].shape, dtype='float64')
    new_field[1:-1,1:-1,1:-1] = fx
    return new_field

@derived_field(name='diff_y_CylR4_Omega2',validators=[ValidateSpatial(ghost_zones=1,fields=["CylR4_Omega2"])])
def diff_y_CylRSqr_Omega(field,data):
    if data.pf["HydroMethod"] == 2:
        sl_left = slice(None,-2,None)
        sl_right = slice(1,-1,None)
        div_fac = 1.0
    else:
        sl_left = slice(None,-2,None)
        sl_right = slice(2,None,None)
        div_fac = 2.0
    dy = div_fac * data['dy'].flat[0]
    fy = (data["CylR4_Omega2"][1:-1,sl_right,1:-1]-data["CylR4_Omega2"][1:-1,sl_left,1:-1])/dy
    new_field = np.zeros(data["CylR4_Omega2"].shape, dtype='float64')
    new_field[1:-1,1:-1,1:-1] = fy
    return new_field

@derived_field(name='diff_r_CylR4_Omega2')
def diff_r_CylRSqr_Omega(field,data):
    return (data['x']*data['diff_x_CylR4_Omega2']+data['y']*data['diff_y_CylR4_Omega2'])/data['CylR']

@derived_field(name='EpicyclicFreqSqr')
def EpicyclicFreqSqr(field,data):
    return data['diff_r_CylR4_Omega2']/data['CylR']**3.

@derived_field(name='absEpicyclicFreqSqr')
def absEpicyclicFreqSqr(field,data):
    return abs(data['EpicyclicFreqSqr'])

@derived_field(name='absEpicyclicFreqSqr_Sqrt')
def absEpicyclicFreqSqr_Sqrt(field,data):
    return abs(data['EpicyclicFreqSqr'])**0.5

@derived_field(name='DeltaCylVelPhi')
def DeltaCylVelPhi(field,data):
        return (data['CylVelPhi']-data['AveragedCylVelPhi'])

@derived_field(name='ReynoldsStress')
def ReynoldsStress(field,data):
	return data['dens']*data['CylVelR']*data['DeltaCylVelPhi']

@derived_field(name='absReynoldsStress')
def absReynoldsStress(field, data):
        return abs(data['ReynoldsStress'])

@derived_field(name='Omega')
def Omega(field,data):
	return data['CylVelPhi']/data['CylR']
	
@derived_field(name='absOmega')
def absOmega(field,data):
	return abs(data['Omega'])

@derived_field(name='AlfvenSpeed')
def AlfvenSpeed(field,data):
	return (2*data['magp']/data['dens'])**0.5

@derived_field(name='Cs')
def Cs(field,data):
	return (data['gamc']*data['pres']/data['dens'])**0.5

@derived_field(name='Lambdac')
def Lambdac(field,data):
	return 4*math.pi*(16/15)**0.5*data['AlfvenSpeed']/data['absOmega']

@derived_field(name='FastSpeed')
def FastSpeed(field,data):
	return 0.5*((data['Cs']**2+data['AlfvenSpeed']**2)+((data['Cs']**2+data['AlfvenSpeed']**2)**2-4*data['Cs']**2*data['AlfvenSpeed']**2)**0.5)

@derived_field(name='dt')
def dt(field,data):
	dt_x=data['dx']/(abs(data['velx'])+data['FastSpeed']**0.5)
	dt_y=data['dy']/(abs(data['vely'])+data['FastSpeed']**0.5)
	dt_z=data['dz']/(abs(data['velz'])+data['FastSpeed']**0.5)
	return np.minimum(dt_x,dt_y,dt_z)

@derived_field(name='LambdacOverdx')
def LambdacOverdx(field,data):
	return data['Lambdac']/data['dx'] 

@derived_field(name='RotationalEnergy')
def RotationalEnergy(field,data):
	return 0.5*data['CylVelPhi']**2

@derived_field(name='MagneticEnergy')
def MagneticEnergy(field,data):
	return 0.5*(data['magx']**2+data['magy']**2+data['magz']**2)

@derived_field(name='EnergyFluxX')
def EnergyFluxR(field,data):
	return data['dens']*data['velx']**2*data['AlfvenSpeed']

@derived_field(name='EnergyFluxY')
def EnergyFluxR(field,data):
	return data['dens']*data['vely']**2*data['AlfvenSpeed']

@derived_field(name='EnergyFluxZ')
def EnergyFluxR(field,data):
	return data['dens']*data['velz']**2*data['AlfvenSpeed']	

@derived_field(name='absBindingEnergy')
def BindingEnergy(field,data):
	return abs(0.5*data['gpot'])

@derived_field(name='Ekin')
def Ekin(field,data):
        return 0.5*(data['velx']**2+data['vely']**2+data['velz']**2)

@derived_field(name='Epot')
def Epot(field,data):
        return 0.5*data['gpot']

@derived_field(name='RatioEpotEkin')
def RatioEpotEner(field,data):
	        return abs(data['Epot'])/data['Ekin']
	
@derived_field(name='MergerArea')
def MergerArea(field,data):
	return (data['eint']>=data['RotationalEnergy'])*(data['dens']>1.e4) #make changes here

@derived_field(name='MergerArea1')
def MergerArea(field,data):
	return (data['SphR']<=.7e9) # make changes here/ previously .5e9 as defined by Ji

@derived_field(name='OtherArea1')
def MergerArea(field,data):
	return (data['SphR']>.5e9)  # make changes here/ 

@derived_field(name='DiskArea')
def DiskArea(field,data):
	return (data['eint']<data['RotationalEnergy'])

@derived_field(name='DiskArea1')    # make changes here/
def DiskArea1(field,data):
        return ((data['CylR'] > 5.0e8) * (data['CylR'] < 3.0e9) * (data['CylZ'] < 7.5e9) * (data['CylZ'] > -7.5e9))

@derived_field(name='MergerAreaDensity')
def MergerAreaDensity(field,data):
	return data['MergerArea']*data['dens']

@derived_field(name='MergerAreaDensity1')
def MergerAreaDensity(field,data):
	return data['MergerArea1']*data['dens']

@derived_field(name='DiskAreaDensity')
def DiskAreaDensity(field,data):
	return data['DiskArea']*data['dens']

@derived_field(name='MergerAreaAngularMomentum')
def MergerAreaAngularMomentum(field,data):
	return data['MergerAreaDensity']*data['CylVelPhi']*data['CylR']

@derived_field(name='DiskAreaAngularMomentum')
def MergerAreaDiskMomentum(field,data):
	return data['DiskAreaDensity']*data['CylVelPhi']*data['CylR']	

@derived_field(name='TotalEnergyMHD')  # Volume averged, not cellmass averaged
def TotalEnergyMHD(field,data):
	return data['dens']*(0.5*(data['velx']**2+data['vely']**2+data['velz']**2)+data['eint']+data['gpot'])+data['magp']

@derived_field(name='AngularMomentum')
def AngularMomentum(field,data):
	return data['CylVelPhi']*data['CylR']*data['dens']

@derived_field(name='absDivB')
def absDivB(field,data):
	return abs(data['divb'])

@derived_field(name='KineticEnergy') # total kinetic energy per cell
def KineticEnergy(field,data):
	return 0.5 * data['dens'] *data['CellVolume'] * (data['velx']**2+data['vely']**2+data['velz']**2)

@derived_field(name='GravX',validators=[ValidateSpatial(ghost_zones=1,fields=["gpot"])])
def GravX(field,data):
    if data.pf["HydroMethod"] == 2:
        sl_left = slice(None,-2,None)
        sl_right = slice(1,-1,None)
        div_fac = 1.0
    else:
        sl_left = slice(None,-2,None)
        sl_right = slice(2,None,None)
        div_fac = 2.0
    dx = div_fac * data['dx'].flat[0]
    fx = (data["gpot"][sl_right,1:-1,1:-1]-data["gpot"][sl_left,1:-1,1:-1])/dx
    new_field = np.zeros(data["gpot"].shape, dtype='float64')
    new_field[1:-1,1:-1,1:-1] = -fx
    return new_field

@derived_field(name='GravY',validators=[ValidateSpatial(ghost_zones=1,fields=["gpot"])])
def GravY(field,data):
    if data.pf["HydroMethod"] == 2:
        sl_left = slice(None,-2,None)
        sl_right = slice(1,-1,None)
        div_fac = 1.0
    else:
        sl_left = slice(None,-2,None)
        sl_right = slice(2,None,None)
        div_fac = 2.0
    dy = div_fac * data['dy'].flat[0]
    fy = (data["gpot"][1:-1,sl_right,1:-1]-data["gpot"][1:-1,sl_left,1:-1])/dy
    new_field = np.zeros(data["gpot"].shape, dtype='float64')
    new_field[1:-1,1:-1,1:-1] = -fy
    return new_field

@derived_field(name='GravCylR')
def GravCylR(field,data):
    return data['GravX']*np.cos(data['CylPhi'])+data['GravY']*np.sin(data['CylPhi'])

@derived_field(name='GravCylPhi')
def GravCylPhi(field,data):
    return -data['GravX']*np.sin(data['CylPhi'])+data['GravY']*np.cos(data['CylPhi'])

@derived_field(name='GravStress')
def GravStress(field,data):
    G=6.67e-8
    return data['GravCylR']*data['GravCylPhi']/(4.*np.pi*G)

#===================================================

def getServerDir():
        # getting export directory
        exportFolder = os.getenv('SERVER_DIR')
        print "Export Folder = ",exportFolder
        if exportFolder is None:
                print "Please define $SERVER_DIR to be the folder where you want to export your plots to"
                exit()
        return exportFolder

def use_subplot(nHorizontal, nVertical):
	F = plt.gcf()
	DPI = F.get_dpi()
	DefaultSize = F.get_size_inches()
	F.set_size_inches( DefaultSize[0]*nHorizontal, DefaultSize[1]*nVertical )

