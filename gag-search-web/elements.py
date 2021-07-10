import numpy as np
import pandas as pd
from bokeh.io import curdoc
from bokeh.layouts import row, column
from bokeh.models import Button,ColumnDataSource,DataTable,Paragraph,Select,Slider,TableColumn
from bokeh.models.widgets import TextInput, Div
from bokeh.plotting import figure


class DbChooser:
    def __init__(self, glycanMZDBsToChoose, defaultElement, fullWidth=800):
        """
        Selector to choose from predefined glycan databases
        """
        self.gList = glycanMZDBsToChoose
        self.selector = Select(title='Choose the glycan m/z database to match against:',
                             value=self.gList[defaultElement], options=self.gList,
                             width=fullWidth)
        self.selector.margin = (0,0,30,0)

class PpmPlot:
    def __init__(self, mzPPMdict, fullWidth):
        self.cds = ColumnDataSource(mzPPMdict)
        self.mzPPMdict = mzPPMdict
        self.fullWidth = fullWidth
        tooltips = [('Glycan:', '@composition'),
                    ('Exp. m/z:', '@mz'),
                    ('PPM error:', '@err_ppm')]
        self.plot = figure(plot_height=250, plot_width=fullWidth,
                           toolbar_location='right', tools='pan,wheel_zoom,box_zoom,reset,save',
                           tooltips=tooltips)
        self.plot.circle(x='mz', y='err_ppm', source=self.cds, size=6, color='navy', alpha=0.5)
        self.plot.xaxis.axis_label = 'Experimental m/z'
        self.plot.yaxis.axis_label = 'Relative mass error, ppm'
        self.plot.xaxis.axis_label_text_font_size = '14px'
        self.plot.yaxis.axis_label_text_font_size = '14px'
        self.plot.xaxis.major_label_text_font_size = '14px'
        self.plot.yaxis.major_label_text_font_size = '14px'
        self.paragr = Paragraph(text='', width=fullWidth, height=30, margin=(20,0,0,0))
        self.calculate_stats()
        
    def calculate_stats(self):
        nMasses = len(self.cds.data['err_ppm'])
        if nMasses > 1:
            errArr = np.array(self.cds.data['err_ppm'])
            meanErr = np.nanmean(errArr)
            meanAbsErr = np.nanmean(np.abs(errArr))
            stDevErr = np.nanstd(errArr)
            firstString = f'PPM error array has {nMasses} values.\n'
        elif nMasses == 1:
            errArr = np.array(self.cds.data['err_ppm'])
            meanErr = errArr[0]
            meanAbsErr = errArr[0]
            stDevErr = 0
            firstString = f'PPM error array has {nMasses} value.\n'
        else:
            meanErr = 0
            meanAbsErr = 0
            stDevErr = 0
            firstString = f'PPM error array has no values.\n'
        self.paragr.text = (firstString + 
                            f'Mean error: {meanErr:.2f} ppm, mean absolute error: {meanAbsErr:.2f} ppm, st. dev: {stDevErr:.2f} ppm')

class RawFileChooser:
    def __init__(self, rawFilesToChoose, paramObj, fullWidth=800):
        """paramObj should be an instance of SearchParamCollector"""
        self.paramObj = paramObj
        self.rawList = rawFilesToChoose
        self.selector = Select(title='Choose the raw files one by one. You can select more than one file. Select a file again to remove:',
                             value=self.rawList[0], options=self.rawList,
                             width=int(fullWidth/2))
        self.selector.margin = (0,0,30,0)
        self.rawFilesToProcess = self.paramObj.returnRawFiles()
        self.paragr = Paragraph(text=('Selected raw files:\n' + str(self.rawFilesToProcess)),
                                width=fullWidth, height=50, margin=(20,0,0,0))

    def on_rawFileUpdate(self, attr, old, newRawFile):
        """Has attr and old arguments for compliance with Bokeh on_change format"""
        self.paramObj.passRawFile(newRawFile)
        self.rawFilesToProcess = self.paramObj.returnRawFiles()
        self.paragr.text = ('Selected raw files:\n' + str(self.rawFilesToProcess))
        return True

class RelAbPlot:
    def __init__(self, relAbDict, fname, fullWidth):
        self.cds = ColumnDataSource(relAbDict)
        self.fname = fname
        self.fullWidth = fullWidth
        tooltips = [('Glycan:', '@composition'),
                    ('Abund. %:', '@percent_abund'),
                    ('Precursors:', '@num_precursors'),
                    ('Avg. RT, min:', '@avg_rt_min')]
        self.plot = figure(x_range=self.cds.data['composition'], plot_height=400, plot_width=fullWidth,
                           title='Twenty Most Abundant Glycans (as % of Total Abundance)',
                           toolbar_location='right', tools='xpan,xwheel_zoom,box_zoom,reset,save',
                           tooltips=tooltips)
        self.plot.vbar(x='composition', top='percent_abund', width=0.9,
                       color='navy', alpha=0.5, source=self.cds)
        self.plot.xgrid.grid_line_color = None
        self.plot.y_range.start = 0
        self.plot.yaxis.axis_label = 'Abundance as % of Total'
        self.plot.xaxis.major_label_orientation = 0.25*np.pi
        self.plot.xaxis.major_label_text_font_size = '14px'
        self.plot.yaxis.major_label_text_font_size = '14px'

        self.div = Div(text=('Showing results for file: ' + self.fname),
                      width=fullWidth, height=50)
        self.div.style = {'background-color': '#00ace6',
                    'color': 'white',
                    #'margin': '5px',
                    'width': (str(fullWidth) + 'px'),
                    'padding-top': '10px',
                    'padding-bottom': '10px',
                    'padding-left': '2px',
                    'box-sizing': 'border-box',
                    'font-size': '22px'}

class ResultTable:
    def __init__(self, cds, fullWidth=800):

        self.colList = [TableColumn(field='composition', title='Glycan Composition'),
                      TableColumn(field='sum_abundance', title='Combined Precursor Intensity'),
                      TableColumn(field='num_precursors', title='Identified Precursors'),
                      TableColumn(field='avg_rt_min', title='Avg Retention Time, min')]

        self.tab = DataTable(source=cds, columns=self.colList,
                             width=fullWidth, height=300)


class SearchParamCollector:
   def __init__(self):
      self.gagDBdict = {}
      self.gagDBPath = None
      self.gagDbName = None
      self.minAbsInt = 0
      self.minSN = 0
      self.searchTolPPM = 10
      self.rawFilesToProcess = []
      
   def clearRawFiles(self):
       self.rawFilesToProcess = []
       return self.rawFilesToProcess
   
   def passDbDictionary(self, dbDict):
      """Passes the dictionary that will allow to match the database name to it's path"""
      self.gagDBdict = dbDict
      return True

   def passGagDbName(self, attr, old, gagDbName):
      """Has attr and old arguments for compliance with Bokeh on_change format"""
      self.gagDbName = gagDbName
      try:
         self.gagDBPath = self.gagDBdict[gagDbName]
         print('Selected the database', gagDbName)
      except:
         print('Could not find the database named', gagDbName)
      return self.gagDBPath

   def passRawFile(self, fname):
      """Takes names of raw files 1 at a time
         If the name is absent from the list, adds it
         If the name is already on the list, it removes the name from the list."""
      n = 3
      if len(fname) > n:
         if fname in self.rawFilesToProcess:
            try:
               self.rawFilesToProcess.remove(fname)
               print(f'Removed file {fname} from selection')
            except:
               print('Could not remove file', fname)
         else:
            self.rawFilesToProcess.append(fname)
            print(f'Added file {fname} to selection')
      else:
         print(f'File name is too short, should be at least {n} characters')
      self.rawFilesToProcess.sort()
      return self.rawFilesToProcess
   
   def printParams(self):
      print(self.gagDBPath,'\n',self.minAbsInt,'\n',self.minSN,'\n',
            self.searchTolPPM,'\n',self.rawFilesToProcess)
      return True

   def setMinAbsInt(self, attr, old, absint):
      """Has attr and old arguments for compliance with Bokeh on_change format"""
      self.minAbsInt = absint
      print('Set min absolute intensity to',self.minAbsInt)
      return self.minAbsInt

   def setMinSN(self, attr, old,  sn):
      """Has attr and old arguments for compliance with Bokeh on_change format"""
      self.minSN = sn
      print('Set min S/N to',self.minSN)
      return self.minSN
 
   def setPPMtol(self, attr, old,  tol):
      """Has attr and old arguments for compliance with Bokeh on_change format"""
      self.searchTolPPM = tol
      print('Set ppm tolerance to',self.searchTolPPM)
      return self.searchTolPPM

   def returnGagDbName(self):
      return self.gagDbName

   def returnGagDbPath(self):
      return self.gagDBPath

   def returnRawFiles(self):
      return self.rawFilesToProcess

   def returnTol_MinInt_MinSN(self):
      """Returns a tuple
         (ppm tolerance, min abs intensity, min S/N)"""
      return (self.searchTolPPM, self.minAbsInt, self.minSN)

class TitleHeader:
   def __init__(self, fullWidth=800):
      self.head = Div(text=('GAG QUANTITATION DASHBOARD'),
                      width=fullWidth, height=100)
      titleStyle = {'background-color': '#00ace6',
                    'color': 'white',
                    #'margin': '5px',
                    'width': (str(fullWidth) + 'px'),
                    'padding-top': '20px',
                    'padding-bottom': '20px',
                    'padding-left': '2px',
                    'box-sizing': 'border-box',
                    'font-size': '30px'}
      self.head.style = titleStyle

def df_to_sorted_perc_dict(dfShort, dfAllFeat, outputLength=20, abundCName='abundance'):
    #dfAllFeat is the LC-MS feature table with all LC-MS features
    
    if len(dfShort.index) > 0:
        #dfShort is the short result table
        dfRes = dfShort.copy()
        dfRes.sort_values(by='sum_abundance', ascending=False, inplace=True)
        sortedDict = {c:list(dfRes[c]) for c in dfRes.columns}
        totalAbund = np.nansum(dfAllFeat[abundCName])
        sortedDict['percent_abund'] = list(dfRes['sum_abundance'] * 100 / totalAbund)
        if len(sortedDict['composition']) <= outputLength:
            toAdd = outputLength - len(sortedDict['composition'])
            sortedDict['composition'] += [('_'+str(x)) for x in range (0,toAdd)]
            sortedDict['sum_abundance'] += ([0] * toAdd)
            sortedDict['num_precursors'] += ([0] * toAdd)
            sortedDict['avg_rt_min'] += ([0] * toAdd)
            sortedDict['percent_abund'] += ([0] * toAdd)
        else:
            sortedDict = {k:sortedDict[k][:outputLength] for k in sortedDict.keys()}
    else:
        sortedDict = {}
        sortedDict['composition'] = [('_'+str(x)) for x in range (0,outputLength)]
        sortedDict['sum_abundance'] = [0] * outputLength
        sortedDict['num_precursors'] = [0] * outputLength
        sortedDict['avg_rt_min'] = [0] * outputLength
        sortedDict['percent_abund'] = [0] * outputLength
    
    return sortedDict

def res_dict_to_mzppm_dict(resdict):
    if len(resdict.keys()) > 0:
        mzPPMdict = {'composition':[],'mz':[],'err_ppm':[]}
        for g in resdict.keys():
            for resultLine in resdict[g]:
                for ion in resultLine[0]:
                    mzPPMdict['composition'].append(g)
                    mzPPMdict['mz'].append(ion[0])
                    mzPPMdict['err_ppm'].append(ion[2])
    else:
        mzPPMdict = {'composition':['_',],'mz':[0,],'err_ppm':[0,]}
    df = pd.DataFrame(mzPPMdict)
    df.to_excel('ppm_errs.xlsx')
    return mzPPMdict
