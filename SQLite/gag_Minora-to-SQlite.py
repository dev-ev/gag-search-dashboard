import datetime
import integrationutils as ius
import numpy as np
import os
import pathlib
import pandas as pd
import sqlite3
import sys
import warnings


''' Direct questions and concerns regarding this script to Egor Vorontsov
    egor.vorontsov@gu.se
'''

def find_abund_col(df):
    clist = []
    for c in df.columns:
        if 'bundances (per File)' in c:
            clist.append(c)
    if len(clist) == 1:
        colname = clist[0]
        print('Found abundance column', colname)
    elif len(clist) == 0:
        colname = None
        print('Could not find abundance column in the Feature table')
    elif len(clist) > 1:
        colname = None
        print('Found more than one potential abundance column in the Feature table')
    return colname

def fname_to_instrument(infname):
    '''Decides on the instrument name based on the file name'''
    infname = infname.lower()
    instrument_name = None
    if ('lumos_' in infname) and ('faims' in infname):
        instrument_name = 'Lumos_FAIMS'
    elif ('lumos_' in infname) and ('faims' not in infname):
        instrument_name = 'Lumos'
    elif ('fusion_' in infname) and ('faims' in infname):
        instrument_name = 'Fusion_FAIMS'
    elif ('fusion_' in infname) and ('faims' not in infname):
        instrument_name = 'Fusion'
    elif ('qehf_' in infname):
        instrument_name = 'QEHF'
    elif ('qe_' in infname):
        instrument_name = 'QE'
    elif ('elite_' in infname):
        instrument_name = 'Elite'

    return instrument_name

def get_console_arg():
    try:
        assert (len(sys.argv) == 2), "Script requires only one parameter"
        
        print('Input file',sys.argv[1])
        injs = sys.argv[1]
        print('Current work dir', os.getcwd())
        return injs
    except:
        print('Could not pass json file to the script')
        return None

def parse_cons_features(in_dfs):
    """ If the Consensus Features table is present, returns the tuple with peak properties
        (for features with charge 2 or higher):
        (mean peak width in s, st dev of peak width in s,
        numpy array of peak widths in s, numpy array of feature intensity)
        If the Consensus Features table is absent, returns
        (None, None, numpy array of zeroes, None)
    """
    if 'Consensus Features' in in_dfs.keys():
        try:
            df = in_dfs['Consensus Features']
            #df1 = df[df['Charge'] > 1]
            diff_array = np.subtract(np.array(df1['Right RT [min]']),np.array(df1['Left RT [min]']))
            diff_array = 60*diff_array
            try:
                for c in df1.columns:
                    if 'Abundances (per File)' in c:
                        abundcol = c
                int_array = np.array(df1[abundcol])
            except:
                warnings.warn('Could not find Consensus Feature intinsities', UserWarning)
                int_array = None

            if int_array is None:
                sum_abund = None
            else:
                sum_abund = np.sum(int_array)
                print('Sum of all intensities is',sum_abund)
            
            return (round(np.nanmean(diff_array),1),
                    round(np.nanstd(diff_array),1),
                    diff_array,
                    sum_abund, int_array)
        except:
            warnings.warn('Could not process Consensus Features table', UserWarning)
            return (None, None, np.zeros(100,dtype=int), None, None)
    else:
        print('Could not find the Consensus Features table')
        return (None, None, np.zeros(100,dtype=int), None, None)

def parse_table_input_file(in_dfs):
    """ 
    Returns a tuple
    (instrument name list, creation date list, instr from name list, instr from metadata list)
    """
    try:
        # Select the Input Files df from the dictionary
        df = in_dfs['Input Files']
        # Select the rows with .RAW files (df in principle can contain many .RAW files, .MSF etc)
        df1 = df[df['File Name'].str.contains('.raw')]
        # Create a list of base filenames for .RAW files
        shortnames = []
        for longname in df1['File Name']:  
            shortnames.append(pathlib.Path(longname).name)
        filedates = list(df1['Creation Date'])
        instr_from_metadata = list(df1['Instrument Name'])
        instr_from_fnames = []
        for n in shortnames:
            instr_from_fnames.append(fname_to_instrument(n))
        return (shortnames,filedates,instr_from_fnames,instr_from_metadata)
    except:
        warnings.warn('Could not process Input Files table', UserWarning)
        return None

def print_all_rows(conn, table):
    cur = conn.cursor()
    sql_command = "SELECT * FROM " + table
    cur.execute(sql_command)

    rows = cur.fetchall()
    for row in rows:
        print(row)
    cur.close()
    return True


def return_last_rows(conn, table, index_col, num, colnames):
    cur = conn.cursor()
    sql_command = ('''SELECT * FROM (
                        SELECT * FROM ''' + table + " ORDER BY " +
                   index_col + ''' DESC LIMIT ''' + str(num) + ''')
                    ORDER BY ''' + index_col + " ASC;")
    cur.execute(sql_command)

    rows = cur.fetchall()

    list_of_tuples = []
    for row in rows:
        list_of_tuples.append(row)

    cur.close()
    df = pd.DataFrame(list_of_tuples, columns=colnames)
    print('Fetched', num, 'latest results from the QC database')
    
    return df

def return_latest_psm_is(df, id_col, file_col, instr_col, psm_col):
    ''' Extracts info on PSM number, search ID and Instrument from the last row in DB
    '''
    last_row = df.iloc[-1]
    search_id = last_row[id_col]
    instr = last_row[instr_col]
    psm = last_row[psm_col]
    psm_string = str(psm) + ' PSMs in file ' + str(last_row[file_col])
    print('String to put on the graph', psm_string)
    return (search_id, instr, psm, psm_string)

def return_peptide_number(in_dfs):
    ''' Returns the number of peptides based on the Peptide table
        '''
    try:
        # Select the Peptide Groups df from the dictionary
        df = in_dfs['Peptide Groups']
        return (len(df.index))
    except:
        warnings.warn('Could not process Peptide Groups table', UserWarning)
        return None

def return_protein_number(in_dfs):
    """Returns the number of Master proteins based on the Proteins table"""
    try:
        # Select the Proteins df from the dictionary
        df = in_dfs['Proteins']
        # Select the Master proteins
        df1 = df[df['Master'] == 'IsMasterProtein']
        return (len(df1.index))
    except:
        warnings.warn('Could not process Proteins table', UserWarning)
        return None

def testing_load_example_files():
    df_dict = {}
    df_dict['Proteins'] = pd.read_csv('TargetProtein.txt',delimiter='\t')
    df_dict['Peptide Groups'] = pd.read_csv('TargetPeptideGroup.txt',delimiter='\t')
    df_dict['PSMs'] = pd.read_csv('TargetPeptideSpectrumMatch.txt',delimiter='\t')
    df_dict['MS/MS Spectrum Info'] = pd.read_csv('MSnSpectrumInfo.txt',delimiter='\t')
    df_dict['Input Files'] = pd.read_csv('WorkflowInputFile.txt',delimiter='\t')
    
    return df_dict

def write_dataFrame(conn, df, featureTabName, rawFileID, cDct):
    
    writing_successful = False
        
    cur = conn.cursor()
    sql_command = ("INSERT INTO " + featureTabName +
                   " VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)")    
    try:
        for r in df.iterrows():
            d = r[1]
            #Order of the columns in the feature DB:
            #feature_id, raw_file_id,
            #feature_mz, charge,
            #max_sn, apex_rt_min,
            #left_rt_min, right_rt_min,
            #isotopes, abundance,
            #comment
            dataTuple = (None, rawFileID,
                         d[cDct['mz']], int(d[cDct['ch']]),
                         d[cDct['SN']], d[cDct['apexRT']],
                         d[cDct['leftRT']], d[cDct['rightRT']],
                         d[cDct['iso']], int(d[cDct['abund']]),
                         None)
            cur.execute(sql_command, dataTuple)
            
        conn.commit()
        writing_successful = True
        print(f'Successfuly saved the feature into the database {featureTabName}' +
              f' with file ID {rawFileID}')
    except:
        print(f'Could not write the results to the database into table {featureTabName}')

    print()    
    cur.close()

    return writing_successful

def write_row(conn, sql_command, data_tuple, tableName):
    
    writing_successful = False
    cur = conn.cursor()
    try:
        cur.execute(sql_command, data_tuple)
        conn.commit()
        writing_successful = True
        print('Successfuly saved the file info into the database:',tableName)
        print(data_tuple)
    except:
        print('Could not write the results to the database into table',tableName)
        print()
        print()
        
    cur.close()

    return writing_successful

if __name__ == '__main__':

    renameDict = {'Avg. m/z [Da]': 'feature_mz',
                  'Charge': 'charge',
                  'Max. Apex S/N': 'max_sn',
                  'Avg. Apex RT [min]': 'apex_rt_min',
                  'Left RT [min]': 'left_rt_min',
                  'Right RT [min]': 'right_rt_min',
                  '# Isotopes': 'isotopes',
                  'Abundances (per File): F50': 'abundance'}

    colNamesDict = {'mz': 'Avg. m/z [Da]', 'ch': 'Charge',
                    'SN': 'Max. Apex S/N', 'apexRT': 'Avg. Apex RT [min]',
                    'leftRT': 'Left RT [min]', 'rightRT': 'Right RT [min]',
                    'iso': '# Isotopes'}
    
    #Modify the path according to your situation
    dbase = 'D:/Data/gag-search-web/SQLite/gagFeaturesDB_v210710.db'
    fileTabName = 'raw_files'
    featureTabName = 'lcms_features'


    #-----------------------------------------@@@@-----------------------------------------@@@@
    

    input_json = get_console_arg()

    # Read the parameters from json files and load the PD tables into dataframes
    input_params = ius.NodeArgs(input_json)
    tableObj = ius.InputTables(input_params.return_all_table_properties())
    in_dfs = tableObj.return_all_tables()

    # The function below returns lists, in case there are multiple files
    # We expect one file per workflow, so choose the 0-th element of each list
    shortnames,filedates,_,instrFromMetadata = parse_table_input_file(in_dfs)
    print('Found raw file info:', shortnames)
    rawFN = shortnames[0]

    df = in_dfs['Consensus Features']
    #print('Columns in the Feature DF:\n', df.columns)
    abund_cname = find_abund_col(df)
    colNamesDict['abund'] = abund_cname
    
    conn = sqlite3.connect(dbase)
    conn.execute("PRAGMA foreign_keys = ON")

    # Let's write the raw file info into the appropriate table
    file_tuple = (None, rawFN, filedates[0], instrFromMetadata[0],
                  None, None, None, None)
    sql_command = ("INSERT INTO " + fileTabName +
                   " VALUES (?, ?, ?, ?, ?, ?, ?, ?)")


    writing_successful = write_row(conn, sql_command, file_tuple, fileTabName)
    
    if writing_successful is True:
        print(f'Ok, wrote the file name {rawFN} into the file info table')
        cur = conn.cursor()
        sql_command = ("SELECT raw_file_id, filename FROM " + fileTabName +
                       " WHERE filename = '" + rawFN + "'")
        cur.execute(sql_command)

        rows = cur.fetchall()
        if len(rows) == 0:
            print('No files found with the file name ', rawFN)
        elif len(rows) == 1:
            print('OK, found file record in the database, ID, file name',
                  rows[0])
            fileID = rows[0][0]
        else:
            print('Spmething is wrong, retrieved more than one files record',
                  rows)
        
        print('File ID is', fileID)
        
        
        cur.close()

        sql_command = ("INSERT INTO " + featureTabName +
                       " VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)")


        data_successful = write_dataFrame(conn, df,
                                          featureTabName, fileID,
                                          colNamesDict)

    conn.close()

    
    



    
