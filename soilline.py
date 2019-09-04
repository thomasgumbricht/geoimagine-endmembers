'''
Created on 12 juni 2018

@author: thomasgumbricht
'''
import os
import sys
import csv

import array as arr
import numpy as np


#import landsat_osgeo_v50 as mj_osgeo

from collections import OrderedDict
from operator import itemgetter
from math import sqrt,atan,cos,sin
from xml.dom import minidom

class SoilLine:
    
    """Layer is the parentid class for all spatial layers."""
    def __init__(self, process): 
        """The constructor expects an instance of the composition class."""
        # LayerCommon.__init__(self)
        self.params = self.proccess.params
        
    def CandidateSoilLine(self):
        vegKernelSize = soilKernelSize = 3            
        soilsdmax = vegsdmax = self.params.process.kernelsd
        pbiMax = self.params.process.pbisoilmax
        pbiMin = self.params.process.pbisoilmin
        pviSoilMean = self.params.process.pvisoilmean
        pviSoilRange = pviIniSoilRange = self.params.process.pvisoilrange
        pviMin = self.params.process.pvivegmin
        searchFraction = self.params.searchpercent*0.01
        nlins = self.params.layerInD['pvi'].lins
        ncols = self.params.layerInD['pvi'].cols
        PVI = self.layerInD['pvi'].BAND
        PBI = self.layerInD['pbi'].BAND
        PVInull = self.layerInD['pvi'].comp.cellnull
        PBInull = self.layerInD['pbi'].comp.cellnull
        adjoinBandD = OrderedDict()
        for key in self.layerInD:
            if key in ['pvi','pbi']:
                continue
            adjoinBandD[key] = self.layerInD[key].BAND
        
        soilD = {'searchpercent' : self.params.searchpercent, 'vegindex': 'pvi'}
        vegD = {'searchpercent' : self.params.searchpercent, 'vegindex': 'pvi'}
        soilD['soilindex'] = 'pbi'
        soilD['kernelsd'] = vegD['kernelsd'] = soilsdmax
        soilD['sisoilmin'] = int(pbiMin)
        soilD['sisoilmax'] = int(pbiMax)
        soilD['visoilmin'] = pviSoilMean - pviSoilRange
        soilD['visoilmax'] = pviSoilMean + pviSoilRange
        soilD['fraction'] = vegD['fraction'] = 0
        soilD['pixels'] = vegD['pixels'] = 0
        soilD['jump'] = vegD['jump'] = 0
        vegD['vivegmin'] = pviMin
        headL = ['pvi','pbi']
        for key in adjoinBandD:
            headL.append(key)
    
        sceneD = {'source':self.layerInD['pvi'].comp.source,'product':self.layerInD['pvi'].comp.product, 'acqdatestr':self.layerInD['pvi'].acqdate.acqdatestr, \
                   'path':self.layerInD['pvi'].location.path, 'row':self.layerInD['pvi'].location.row, 'srcmethod':'pvipbidefault','bounds':self.process.bounds, \
                   'srcdata': self.layerInD['rl'].comp.folder}
        print ('    Checking number of valid cells...'),
        nonValid = len(np.where(PVI==PVInull)[0])
        nValid = nlins*ncols-nonValid
        print ('valid cells in region', nValid)
        if nValid < 500/searchFraction:
            self._WriteSearchResults(self.process,headL,[[]],[[]],sceneD, soilD, vegD)
            print ('SKIPPING - too few clear pixels')
            return sceneD, False, False
        #validfraction = nValid/(nlins*ncols)  
        allSearch = max(500,int(nValid*searchFraction))
        print ( '    Scene search:', allSearch ) 
        _bareSoil = list(); _denseVeg = list()
        kernelSize = max([vegKernelSize, soilKernelSize])
        kernelhalf = kernelSize-1
        #Alway aim for 1500
        linjump = coljump = max(1,int(round(sqrt(allSearch/1500))))
        nSearch = min(1500, allSearch / (linjump*coljump))
        s = nSearch+1 # s is the number of soil points identified, if too large the parameters are changed
        v = nSearch/2+1 
        m = 0 
        print ('    Searching for spectral end-members, targetting for:', int(nSearch) )
        _prevBareSoil = False
        _prevDenseVeg = False
        veg = soil = False
        while (s > nSearch or v > nSearch/2):
            m += 1
            if m > 1:
                #save the previous results in an independent vector
                _prevBareSoil = _bareSoil[:]
                _prevDenseVeg = _denseVeg[:] 
            if not soil:
                _bareSoil = []
            if not veg:
                _denseVeg = []
            soilPVImin = pviSoilMean - pviSoilRange
            soilPVImax = pviSoilMean + pviSoilRange    
            printstr = '    Nr of lines: %(l)d, skipping to every %(j)d line' %{'l':nlins,'j':linjump} 
            print (printstr)
            printstr = '    Nr of columns: %(c)d, skipping to every %(j)d column' %{'c':ncols,'j':coljump} 
            print (printstr)
            for lin in range(kernelhalf,nlins-kernelhalf):
                if (lin%linjump !=0): #TGTODO only check every yth line, otherwise it takes too long
                    continue            
                ln = lin
                for col in range (kernelhalf,ncols-kernelhalf):
                    if (col%coljump != 0): #only check every xth column, otherwise it takes too long
                        continue
                    cl = col
                    if PBI.item((ln, cl)) == PBInull:
                        continue
                    #Soil search
                    PVIsd = np.std(PVI[ln-1:ln+2,cl-1:cl+2])
                    if PVIsd >= max(soilsdmax,vegsdmax):
                        continue
                    if (pbiMin <= PBI.item((ln, cl)) <= pbiMax and soilPVImin <= PVI.item((ln, cl)) <= soilPVImax and not soil): 
                        soiljump = linjump
                        if PVIsd < soilsdmax:
                            bareSoil = list()
                            bareSoil.append(PVI.item((ln, cl)))
                            bareSoil.append(PBI.item((ln, cl)))
                            for key in adjoinBandD:
                                bareSoil.append(adjoinBandD[key].item((ln, cl)))  
                            #append as row to the global list
                            _bareSoil.append(bareSoil)
                    #vegetation search
                    if (PVI.item((ln, cl)) > pviMin and not veg):
                        vegjump = linjump
                        if PVIsd < vegsdmax:
                            denseVeg = list() 
                            denseVeg.append(PVI.item((ln, cl)))
                            denseVeg.append(PBI.item((ln, cl)))
                            for key in adjoinBandD:
                                denseVeg.append(adjoinBandD[key].item((ln, cl)))                               
                            _denseVeg.append(denseVeg)
            s = len(_bareSoil)  
            v = len(_denseVeg)
            print ('    iteration result:', m)
            printstr = '        soil-pixels: %(d1)d; veg-pixels:%(d2)d' %{'d1':s, 'd2': v}
            print (printstr)
            if m > 5:
                break
            if not veg:
                if (v < nSearch/2):
                    if _prevDenseVeg and len(_prevDenseVeg) > len(_denseVeg):
                        _denseVeg = _prevDenseVeg[:] 
                    pviMin *= 0.9; vegsdmax *= 1.1
                else:
                    veg = True
            if not soil:   
                if s <  nSearch:
                    if _prevBareSoil and len(_prevBareSoil) > len(_bareSoil):
                        _bareSoil = _prevBareSoil[:]
                    soilsdmax *= 1.1; 
                    pviSoilRange = min(pviIniSoilRange*1.33,pviSoilRange+25)
                    if linjump > 1:
                        linjump -= 1
                    if coljump > 1:
                        coljump -= 1 
                    nSearch = allSearch / (linjump*coljump)
                else:
                    soil = True
            if veg and soil:
                break
            print ('        Changing search criteria')
        printstr = '    Target pixels:%(t)d, Found soil-pixels: %(d1)d ; veg-pixels:%(d2)d' %{'t': nSearch,'d1':len(_bareSoil), 'd2': len(_denseVeg)}  
        print (printstr) 
        #Sort _baresoil after PBI
        _bareSoil.sort(key=itemgetter(1))
        #Sort _denseveg after PVI
        _denseVeg.sort(key=itemgetter(0)) 
        printstr = '    Final soil-pixels: %(d1)d; veg-pixels:%(d2)d' %{'d1':len(_bareSoil), 'd2': len(_denseVeg)}  
        print (printstr)
        print ('    Saving final search criteria to db')
        soilfraction = s*coljump*linjump/nValid
        vegfraction = v*coljump*linjump/nValid
    
        soilD['jump'] = soiljump
        vegD['jump'] = vegjump
        soilD['soilindex'] = 'pbi'
        soilD['kernelsd'] = soilsdmax
        soilD['sisoilmin'] = int(pbiMin)
        soilD['sisoilmax'] = int(pbiMax)
        soilD['visoilmin'] = int(soilPVImin)
        soilD['visoilmax'] = int(soilPVImax)
        soilD['fraction'] = soilfraction
        soilD['pixels'] = s
        vegD['kernelsd'] = vegsdmax
        vegD['vivegmin'] = int(pviMin)
        vegD['fraction'] = vegfraction
        vegD['pixels'] = v
        #first line
        headL = ['pvi','pbi']
        for key in adjoinBandD:
            headL.append(key)
        with open(self.layerOutD['soil'].FPN, "wb") as f:
            writer = csv.writer(f)
            writer.writerow(headL)
            writer.writerows(_bareSoil)   
        with open(self.layerOutD['veg'].FPN, "wb") as f:
            writer = csv.writer(f)
            writer.writerow(headL)
            writer.writerows(_denseVeg)
        soilxmlFPN = self.layerOutD['soil'].FPN.replace('.csv','.xml')
        vegxmlFPN = self.layerOutD['veg'].FPN.replace('.csv','.xml')
        layer = self.layerInD['pvi']
        self._WriteSearchxml(sceneD, soilD, soilxmlFPN)
        self._WriteSearchxml(sceneD, vegD, vegxmlFPN) 
        return sceneD,soilD,vegD 
    
    def LoopVIsoil(self,allvalues,retrievefraction,viSoilMin,viSoilMax):
        BSLArray = []
        minArray = []
        for values in allvalues:      
            values.sort(key=itemgetter(2)) #lowest nir
            sp = int(len(values)*retrievefraction)
            for i in range(sp):
                if (viSoilMin <= values[i][0] <= viSoilMax):
                    BSLArray.append(values[i])
            for i in range(sp):
                if (viSoilMin <= values[i][0] <= viSoilMax):
                    minArray.append(values[0])
                    break
        return BSLArray, minArray
    
    def WriteEMxml(self, emD, emsampleD, emparamD, trimD, emxmlFPN):
           
        doc = minidom.Document()
        root = doc.createElement('endmembers')
        doc.appendChild(root)
        for ptD in emparamD:
            paramtag = doc.createElement(ptD)
            for p in emparamD[ptD]:
                val = '%s' %(emparamD[ptD][p])
                paramtag.setAttribute(p, val)
            root.appendChild(paramtag)
            if ptD == 'soilline' and emparamD[ptD]['trimming'] == 'Y':
                trimtag = doc.createElement('trimming')
                paramtag.appendChild(trimtag)
                for band in trimD:
                    bandtag = doc.createElement('band')
                    bandtag.setAttribute('id', band)
                    for item in trimD[band]:
                        val = '%s' %(trimD[band][item])
                        bandtag.setAttribute(item, val)
                    trimtag.appendChild(bandtag)        
        for em in emD:
            emtag = doc.createElement('endmember')
            emtag.setAttribute('id', em)
            val = '%s' %(emsampleD[em])
            emtag.setAttribute('n', val)
            root.appendChild(emtag)
            for band in emD[em]:
                bandtag = doc.createElement('band')
                bandtag.setAttribute('id', band)
                for item in emD[em][band]:
                    val = '%s' %(emD[em][band][item])
                    bandtag.setAttribute(item, val)
                emtag.appendChild(bandtag)   
        xml_str = doc.toprettyxml(indent="  ")
        with open(emxmlFPN, "w") as f:
            f.write(xml_str)
    
    def _WriteSearchxml(self, sceneD, soilD, soilxmlFPN):  
        doc = minidom.Document()
        root = doc.createElement('scenesearch')
        doc.appendChild(root)
        scenetag = doc.createElement('scene')
        for p in sceneD:
            val = '%s' %(sceneD[p])
            scenetag.setAttribute(p, val)
        root.appendChild(scenetag)
        searchtag = doc.createElement('search')
        for p in soilD:
            val = '%s' %(soilD[p])
            searchtag.setAttribute(p, val)
        root.appendChild(searchtag)
        xml_str = doc.toprettyxml(indent="  ")
        with open(soilxmlFPN, "w") as f:
            f.write(xml_str)
            
    def WriteSoilLinexml(self,soillineD, slxmlFPN, soillineparamD, trimD):
        doc = minidom.Document()
        root = doc.createElement('soillines')
        doc.appendChild(root)
        paramtag = doc.createElement('parameters')
        for p in soillineparamD:
            val = '%s' %(soillineparamD[p])
            paramtag.setAttribute(p, val)
        if soillineparamD['trimming'] == 'Y':
            trimtag = doc.createElement('trimming')
            paramtag.appendChild(trimtag)
            for band in trimD:
                bandtag = doc.createElement('band')
                bandtag.setAttribute('id', band)
                for item in trimD[band]:
                    val = '%s' %(trimD[band][item])
                    bandtag.setAttribute(item, val)
                trimtag.appendChild(bandtag) 
        root.appendChild(paramtag)    
        for band in soillineD:
            bandtag = doc.createElement('soilline')
            bandtag.setAttribute('id', band)
            root.appendChild(bandtag)
            for item in soillineD[band]:
                val = '%s' %(soillineD[band][item])
                bandtag.setAttribute(item, val)
            root.appendChild(bandtag)
        xml_str = doc.toprettyxml(indent="  ")
        with open(slxmlFPN, "w") as f:
            f.write(xml_str)  
            
    def ReadsoillineXML(self, slxmlFPN):

        dom = minidom.parse(slxmlFPN)   
        params = ['defsamples', 'inputsamples', 'seek', 'slsamples', 'totsamples', 'trimming']
        paramtag = dom.getElementsByTagName('parameters')
        soillineparamD = {}
        soillineD = {}
        trimD = {}
        for param in params:
            if param == 'trimming':
                soillineparamD[param] = str(paramtag[0].getAttribute(param))
            else:
                soillineparamD[param] = int(paramtag[0].getAttribute(param))
        trimtag = paramtag[0].getElementsByTagName('trimming')
        bandtags = trimtag[0].getElementsByTagName('band')  
        for bandtag in bandtags:
            band = bandtag.getAttribute('id')
            mini = int(bandtag.getAttribute('min'))
            maxi = int(bandtag.getAttribute('max'))
            trimD[band] = {'min':mini, 'max':maxi}
        sltags = dom.getElementsByTagName('soilline')
        for sltag in sltags:
            band = sltag.getAttribute('id')
            r2 = float(sltag.getAttribute('r2'))
            slope = float(sltag.getAttribute('slope'))
            intercept = int(sltag.getAttribute('intercept'))
            stderr = float(sltag.getAttribute('stderr'))
            soillineD[band] = {'r2':r2, 'intercept':intercept,'slope':slope,'stderr':stderr}
        return soillineD,soillineparamD,trimD  
    
    def WriteSoilSpectraXML(self,soilxmlFPN,soilL):  
        doc = minidom.Document()
        root = doc.createElement('soilspectra')
        doc.appendChild(root)
        for band in soilL[0]:
            bandtag = doc.createElement('band')
            bandtag.setAttribute('id', band)
            root.appendChild(bandtag)
            root.appendChild(bandtag)
        for spectraD in soilL: 
            spectratag = doc.createElement('spectra')   
            for item in spectraD:  
                '''
                specratag.setAttribute('id', band)
                root.appendChild(bandtag)
                for item in soillineD[band]:
                '''
                val = '%s' %(spectraD[item])
                spectratag.setAttribute(item, val)
            root.appendChild(spectratag)
        xml_str = doc.toprettyxml(indent="  ")
        with open(soilxmlFPN, "w") as f:
            f.write(xml_str)
            
    def ReadSoilSpectraXML(self,soilxmlFPN):
        dom = minidom.parse(soilxmlFPN)   
        bandL = []
        spectraL = []
        bandtags = dom.getElementsByTagName('band')
        for bandtag in bandtags:
            bandL.append( bandtag.getAttribute('id') ) 
        spectratags = dom.getElementsByTagName('spectra')
        for spectratag in spectratags:
            spectraD = {}
            for band in bandL:
                spectraD[band] = spectratag.getAttribute(band)
            spectraL.append(spectraD)
        return spectraL

   
    def WriteSoillinePlotCSV(self, band,slLayersOutD,soillineD,emD):
        #Soil lines are always written for endmembers darksoil, brightsoil
        dark = emD['darksoil'][band]['mean']
        bright = emD['brightsoil'][band]['mean']
        #and then in the regression
        nirdark = dark*soillineD[band]['slope']+soillineD[band]['intercept']
        nirbright = bright*soillineD[band]['slope']+soillineD[band]['intercept']
        sl = [[dark,nirdark],[bright,nirbright]]
        #open file and write the point pairs
        with open(slLayersOutD[band], "wb") as f:
            writer = csv.writer(f)
            writer.writerows(sl)