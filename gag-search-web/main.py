from bokeh.io import curdoc
from bokeh.layouts import row, column
from bokeh.models import Button,ColumnDataSource,DataTable,Paragraph,Select,Slider,TableColumn
from bokeh.models.widgets import TextInput, Div
from bokeh.plotting import figure
import configparser
import pandas as pd
import os

from gag_search import retrieve_all_fnames, search_main_single_sql
from elements import DbChooser, SearchParamCollector, RawFileChooser
from elements import RelAbPlot, ResultTable, PpmPlot, TitleHeader
from elements import res_dict_to_mzppm_dict, df_to_sorted_perc_dict

CONFIG_FNAME = 'config.ini'

DEFAULT_DF = pd.DataFrame((('',0,0,0),('',0,0,0)),
                          columns=('composition','sum_abundance','num_precursors','avg_rt_min'))

class SearchLauncher:
    def __init__(self, paramObj, savePath, dbPath, fileTabName, featureTabName,
                 cds, rawFilesWidget, abPlotWidget, tabWidget, errPlotWidget):
        """
        paramObj is an instance of SearchParamCollector from elements
        tabWidget is an instance of ResultTable from elements
        cds is an instance of ColumnDataSource from bokeh
        procStatus is a Div object from bokeh
        """
        self.paramObj = paramObj
        self.savePath = savePath
        self.dbPath = dbPath
        self.fileTabName = fileTabName
        self.featureTabName = featureTabName
        self.rawFilesWidget = rawFilesWidget
        self.abPlotWidget = abPlotWidget
        self.tabWidget = tabWidget
        self.errPlotWidget = errPlotWidget
        self.cds = cds
        self.resultDfList = []

    def start_search(self):
        tol, minAbsInt, minSN = self.paramObj.returnTol_MinInt_MinSN()
        flist = self.paramObj.returnRawFiles()
        gagDB = self.paramObj.returnGagDbPath()
        if len(flist) == 0:
            print('No files to process, did not launch the search')
        else:
            print(f'Starting search with tol {tol} ppm, min S/N {minSN} and min abs intensity {minAbsInt}')
            print(f'Matching against glycan database {self.paramObj.returnGagDbName()}')
                        
            #Calculate the search results
            self.resultDfList = search_main_single_sql(flist, gagDB, self.savePath,
                                                       self.dbPath, self.fileTabName, self.featureTabName,
                                                       ppmTol=tol, minSN=minSN, minAbsInt=minAbsInt)
            #Clear the list of raw files that has just been searched
            clearedFList = self.paramObj.clearRawFiles()
            self.rawFilesWidget.paragr.text = ('Selected raw files:\n' + str(clearedFList))
                        
            #Each element of the resultDfList has format:
            # (fname, short_out_df, long_out_df, res_dict, df_filtered)
            
            #Just take the first filename and the first short dataframe from the results to display
            displayedRawFname = self.resultDfList[0][0]
            df = self.resultDfList[0][1]
            newDataDict = {c:list(df[c]) for c in df.columns}
            self.cds.data = newDataDict
            self.tabWidget.tab.update()
                        
            #Update the mass error scatter plot
            resDict = self.resultDfList[0][3]
            newErrDict = res_dict_to_mzppm_dict(resDict)
            self.errPlotWidget.cds.data = newErrDict
            self.errPlotWidget.calculate_stats()
            self.errPlotWidget.plot.update()
                        
            #Update the relative abundance plot
            dfFiltered = self.resultDfList[0][4]
            dictRelAb = df_to_sorted_perc_dict(df, dfFiltered, outputLength=20, abundCName='abundance')
            self.abPlotWidget.cds.data = dictRelAb
            self.abPlotWidget.plot.x_range.factors = dictRelAb['composition']
            self.abPlotWidget.div.text = ('Showing results for file: ' + displayedRawFname)
                        
            return True
        
config = configparser.ConfigParser()
config.read(
    os.path.join(
        os.path.dirname(__file__), CONFIG_FNAME
    )
)
print('Parsed config sections:', config.sections() )
full_width = config.getint('params', 'full_width')
gag_mz_libs = dict(config['database_paths'])
exp_path = config['export_path']['path']
spectra_db_path = config['spectra_db']['path']
file_table_name = config['spectra_db']['file_table_name']
feature_table_name = config['spectra_db']['feature_table_name']

cds = ColumnDataSource(DEFAULT_DF)

#Intitialize the object that collects the search parameters
params = SearchParamCollector()
params.passDbDictionary(gag_mz_libs)

#Construct the Bokeh dashboard
titleHead = TitleHeader(full_width)

textPath = Paragraph(text=('Path for exporting the reports: ' + exp_path), width=full_width, height=50)

availableDBs = list(gag_mz_libs.keys())
defaultElement = 0
gagDBSel = DbChooser(availableDBs, defaultElement, full_width)
params.passGagDbName(None, None, availableDBs[defaultElement])
gagDBSel.selector.on_change('value', params.passGagDbName)

ppmSlider = Slider(start=1, end=50, value=10, step=1,
                   title='Choose relative mass tolerance in ppm')
ppmSlider.on_change('value', params.setPPMtol)
snSlider = Slider(start=0, end=50, value=0, step=1,
                   title='Choose minimum precursor S/N value to be considered')
snSlider.on_change('value', params.setMinSN)
intSlider = Slider(start=0, end=20000, value=0, step=1000,
                   title='Choose minimum absolute precursor intensity to be considered')
intSlider.on_change('value', params.setMinAbsInt)

rawFilesToChoose = retrieve_all_fnames(dbPath = spectra_db_path, fileTable = file_table_name)
rawFiles = RawFileChooser(rawFilesToChoose, params, full_width)
rawFiles.selector.on_change('value', rawFiles.on_rawFileUpdate)

#Pass empty dataframe to intitiate the plot. It will get non-zero values after the search
dictRelAb = df_to_sorted_perc_dict(pd.DataFrame(), pd.DataFrame(), outputLength=20, abundCName='abundance')
relAbunds = RelAbPlot(dictRelAb, '', full_width)
results = ResultTable(cds, full_width)

#Pass empty dictionary to intitiate the plot. It will get non-zero values after the search
mzPPMdict = res_dict_to_mzppm_dict(dict())
errorPlot = PpmPlot(mzPPMdict, full_width)

#Initiate an instance of the Launcher class
launch = SearchLauncher(params, exp_path, spectra_db_path, file_table_name, feature_table_name,
                        cds, rawFiles, relAbunds, results, errorPlot)
runButton = Button(label='Start Search', button_type='success', width=int(full_width/4))
runButton.on_click(launch.start_search)


curdoc().add_root(column(titleHead.head, textPath, gagDBSel.selector,
                         ppmSlider, snSlider, intSlider,
                         rawFiles.paragr, rawFiles.selector, runButton,
                         relAbunds.div, relAbunds.plot, results.tab,
                         errorPlot.paragr, errorPlot.plot))

curdoc().title = 'GAG Quan Dashboard'
