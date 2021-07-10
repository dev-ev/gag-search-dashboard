import glob
from pathlib import Path
import pandas as pd
import numpy as np
from sklearn.cluster import DBSCAN
from sklearn.cluster import MeanShift
import collections
import datetime
import json
import pickle
import sqlite3

COL_NAMES = {
    'mz_colname': 'feature_mz', 'theo_mz_colname': 'theo_mz',
    'charge_colname': 'charge', 'sn_colname': 'max_sn',
    'isotope_colname': 'isotopes', 'rt_colname': 'apex_rt_min',
    'label_colname': 'cluster', 'abund_colname': 'abundance',
    'sql_fileID': 'raw_file_id', 'rawFile_colname': 'filename'
}

SQL_COLLIST = [
    'feature_id','raw_file_id','feature_mz','charge',
    'max_sn','apex_rt_min','left_rt_min','right_rt_min',
    'isotopes','abundance','comment'
]

DBA_MASS = 129.15175
PROTON_MASS = 1.00728
SO3_MASS = 79.95682

N_CORES = 6

BANDWIDTH = 0.075
DBSCAN_EPS = 0.06
MIN_MATCH_PERCENT = 1


DB_FORMAT = '*.json'
INPUT_FORMAT = '*.txt'
XLSX_SUFF = '_id_abund_%.fppm_%s_%s.xlsx'

def filter_feature_df(df,min_charge=1,min_mz=250,min_isotopes=3):

    df_out = df[
        (df[COL_NAMES['charge_colname']] >= min_charge) &
        (df[COL_NAMES['mz_colname']] >= min_mz) &
        (df[COL_NAMES['isotope_colname']] >= min_isotopes)
    ]
    
    return df_out.copy()

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

def find_longest_result(input_reslist):
    output_reslist = []
    max_l = max([len(i[1]) for i in input_reslist])
    output_reslist = [x for x in input_reslist if (len(x[1]) == max_l)]
    return output_reslist

def filter_matching_MZs(inDf, allMZsDict, colsDict, relTol=1E-5):
    """
    Takes the dataframe with all the features,
    and the disctionary with theorethical m/zs:
        {charges: [lsit of m/z at that charge], }
    returns the dataframe with the rowas that match the theoretchical ions by m/z and charge
    """

    #Dividing the tolerance by 2 since we use the SUM matrix for division
    #when we calculate ther relative mass difference: np.divide(product_matrix,sum_matrix)
    #instead of using the average values that would have been sum_matrix/2.
    #Thus the relative delta is 2 times lower than expected after that operation
    relTol = relTol/2 
    matched_labels = []
    for ch in allMZsDict.keys():
        exp_data_by_charge = inDf[inDf[colsDict['charge_colname']] == ch].copy()
        if len(exp_data_by_charge.index) > 0:
            expmzs = np.array(exp_data_by_charge[colsDict['mz_colname']])
            theo_mzs = np.array(allMZsDict[ch])
            product_matrix = np.subtract.outer(theo_mzs,expmzs)
            product_matrix = np.absolute(product_matrix)
            sum_matrix = np.add.outer(theo_mzs,expmzs)
            rel_delta_matrix = np.divide(product_matrix,sum_matrix)
            rel_delta_matrix = pd.DataFrame(rel_delta_matrix)
            bool_matrix = rel_delta_matrix.le(relTol)
            matching_mz_indices = np.where(bool_matrix == True)
            
            if matching_mz_indices[1].shape[0] > 0:
                exp_indices = list(matching_mz_indices[1])
                ##print('Matched masses', len(exp_indices), 'for charge', ch)
                matchedDf = exp_data_by_charge.iloc[exp_indices].copy()
                ##print(matchedDf.index)
                matched_labels += list(matchedDf.index)
                
    matched_labels = list(set(matched_labels))

    matchedDf = inDf.loc[matched_labels,:].copy()
    matchedDf.sort_values(by=[colsDict['rt_colname']],inplace=True,ascending=True)
    
    return matchedDf

def get_lib_MZs_by_charge(list_of_list, identicalMZsThreshold = 5E-7):
    """
    Takes the dtabase list of list (imported from json).
    Retruns a dictionary:
        charge: [all theorethical m/z in the DB with that charge]
    """

    all_masses_inDB = {}
    for g in list_of_list:
        all_masses_for_glycan = g[1] + g[2]
        for m in all_masses_for_glycan:
            ch = int(m[1])
            if ch not in all_masses_inDB.keys():
                all_masses_inDB[ch] = []
            all_masses_inDB[ch].append(m[0])
    for arr in all_masses_inDB.values():
        arr.sort()

    # Exclude the thorethical m/z that are too close to each other
    # Use the relative threshold provided as a prameter
    for k in all_masses_inDB.keys():
        rangeToCheck = list(range(len(all_masses_inDB[k])-1))
        filteredMZs = [all_masses_inDB[k][x] for x in rangeToCheck if not
                       (np.isclose(all_masses_inDB[k][x], all_masses_inDB[k][x+1],
                                   rtol=identicalMZsThreshold))]
        if not (np.isclose(all_masses_inDB[k][-1],filteredMZs[-1],rtol=identicalMZsThreshold)):
            filteredMZs.append(all_masses_inDB[k][-1])
        all_masses_inDB[k] = filteredMZs
    
    return all_masses_inDB

def get_result_DFs(resDict, minInt = 0):
    """
    Takes the result dictionary from the search function.
    Returns two dataframes with the results:
    1) short dataframe, with only one signal for each glycan. The reported glycan is chosen as follows:
        - the minimum absolute abundance threshold is applied first
        - the group with the most coeluting precursors is chosen from the filtered results
            - if there is 1 result at that number of precursors, it is reported
            - if there is more than 1 group with the given max number of precursors, the most abundant one is chosen
    2) long dataframe with all the matches for a given glycan composition
    """

    short_out_dict = {'composition':[], 'sum_abundance':[],
                      'num_precursors':[], 'avg_rt_min':[]}
    long_out_dict = {'composition':[], 'sum_abundance':[],
                      'num_precursors':[], 'avg_rt_min':[]}
    
    for k in resDict.keys():
        results = resDict[k]
        # Add all the result lines to the long result dictionary
        for result_line in results:
            long_out_dict['composition'].append(k)
            long_out_dict['sum_abundance'].append(int(result_line[3]))
            long_out_dict['num_precursors'].append(len(result_line[0]))
            long_out_dict['avg_rt_min'].append(np.mean(result_line[1]))

        # Filter the identifications by the minumum absolute abundance
        filteredResults = [x for x in results if x[3] >= minInt]

        if len(filteredResults) > 0:
            # Select the group with the highest numer of precursors
            all_lengths = [len(x[0]) for x in filteredResults]
            most_precursors = [x for x in filteredResults if len(x[0]) == max(all_lengths)]
            if len(most_precursors) == 1:
                to_report = most_precursors[0]
            else:
                ints = [x[3] for x in most_precursors]
                maxint = [x for x in most_precursors if x[3] == max(ints)]
                if len(maxint) > 1:
                    print(('There were two IDs with the same intensity and number of precursors for ' +
                           k + ', reporting only the first one'))
                to_report = maxint[0]
            short_out_dict['composition'].append(k)
            short_out_dict['sum_abundance'].append(int(to_report[3]))
            short_out_dict['num_precursors'].append(len(to_report[0]))
            short_out_dict['avg_rt_min'].append(np.mean(to_report[1]))            

    shortDf = pd.DataFrame(short_out_dict, columns=['composition','sum_abundance',
                                                    'num_precursors','avg_rt_min'])
    shortDf.sort_values(by=['composition'], inplace=True)
    longDf = pd.DataFrame(long_out_dict, columns=['composition','sum_abundance',
                                                  'num_precursors','avg_rt_min'])
    longDf.sort_values(by=['composition'], inplace=True)
         
    return shortDf, longDf

def max_retTimeSpan(clusteredDf):
    spans = []
    for c in set(clusteredDf[COL_NAMES['label_colname']]):
        clDf = clusteredDf[clusteredDf[COL_NAMES['label_colname']] == c]
        rtArray = clDf[COL_NAMES['rt_colname']]
        span = max(rtArray) - min(rtArray)
        spans.append(span)
    maxSpan = max(spans)
    medianSpan = np.median(spans)
    avgSpan = np.mean(spans)
    print(f'Maximum span for a cluster was {maxSpan:.2f}, median span {medianSpan:.2f} and mean span {avgSpan:.2f} (min)')
    return maxSpan, medianSpan, avgSpan

def retrieve_all_fnames(dbPath, fileTable, colDict=COL_NAMES):
    allFiles = []

    conn = sqlite3.connect(dbPath)
    conn.execute("PRAGMA foreign_keys = ON")
    cur = conn.cursor()

    sql_command = ("SELECT " + colDict['rawFile_colname'] +
                   " FROM " + fileTable + ";")
    #print('SQL command to execute:\n', sql_command)
    
    cur.execute(sql_command)
    rows = cur.fetchall()
    #print(rows)
    
    allFiles = [x[0] for x in rows]
    allFiles = list(set(allFiles))
    allFiles.sort()
    
    conn.commit()
    conn.close()
    
    return allFiles

def retrieveDF_from_SQLite(rawFName, dbPath, fileTabName, featureTabName,
                           sqlColList, colDict):

    conn = sqlite3.connect(dbPath)
    conn.execute("PRAGMA foreign_keys = ON")
    cur = conn.cursor()

    sql_command = ("SELECT * FROM " + featureTabName +
                   " WHERE " + colDict['sql_fileID'] +
                   " IN (SELECT " + colDict['sql_fileID'] +
                   " FROM " + fileTabName +
                        " WHERE " + colDict['rawFile_colname'] +
                   " = '" + rawFName + "');")
    #print('SQL command to execute:\n', sql_command)
    
    cur.execute(sql_command)
    rows = cur.fetchall()
    print(f'Retrieved {len(rows)} rows from the database for file {rawFName}')
    #print(rows[:2])
    #print(rows[-2:])

    df = pd.DataFrame(rows, columns=sqlColList)
    print(df)

    conn.commit()
    conn.close()
    
    return df

def return_ppm_diff_list(exp_mz_list, theo_mz_list):
    exp_arr = np.array(exp_mz_list)
    theo_arr = np.array(theo_mz_list)
    ppm_arr = 1000_000 * 2 * ((exp_arr - theo_arr) / (exp_arr + theo_arr))
    
    return list(ppm_arr)

def rt_clustering_dbscan(df, eps, min_neighb=1,metric='euclidean'):
    all_rt_df = df[[COL_NAMES['rt_colname']]].copy()
    print(all_rt_df.head(10))
    X = np.array(all_rt_df)
    model_dbscan = DBSCAN(eps=eps, min_samples=min_neighb)
    model_dbscan.fit(X)
    print('Found clusters', len(np.unique(model_dbscan.labels_)))
##    print(collections.Counter(model_dbscan.labels_))
    all_rt_df[COL_NAMES['label_colname']] = model_dbscan.labels_
    df[COL_NAMES['label_colname']] = model_dbscan.labels_

    return model_dbscan.labels_

def rt_clustering_meanshift(df, bandwidth):
    all_rt_df = df[[COL_NAMES['rt_colname']]].copy()
    print(all_rt_df.head(10))
    X = np.array(all_rt_df)
    ms = MeanShift(bandwidth=bandwidth, cluster_all=True)
    ms.fit(X)
    print('Found clusters', len(np.unique(ms.labels_)))
##    print(collections.Counter(ms.labels_))
    all_rt_df[COL_NAMES['label_colname']] = ms.labels_
    df[COL_NAMES['label_colname']] = ms.labels_

    return ms.labels_


def search_gags_ml_v2(df, sorted_gag_list, mass_tol, min_match_perc,
                      col_names):
    """The improved version relying on numpy and the comprehensive library"""
    
    if len(df.index) > 0:
        res_list = []
        sorted_label_list = sorted(list(set(df[col_names['label_colname']])))
        for lbl in sorted_label_list:
            ldf = df[df[col_names['label_colname']] == lbl]
                
            results_for_cluster = search_gags_v2(ldf, sorted_gag_list, col_names,
                                                 mass_tolerance=mass_tol,
                                                 min_match_percent=min_match_perc)
            res_list += results_for_cluster
        return res_list
    else:
        return []

def search_gags_v2(df, gaglist, col_names,
                   mass_tolerance=1e-5, min_match_percent=0):
    """
    The improved version relying on numpy and the comprehensive library
    Gaglist should be in a particular format from the translate_lib function,
    it should be also sorted in the order of the descending number of signals (seniority)
    Returns a list with the search output:
    [[saccharide name, [(mz1,mz1_theo,ppm diff,charge1),(mz2,mz2_theo,ppm diffcharge2)...],
    [RT1, RT2...], [abund1, abund2...], sum abund], ]
    """

    min_match_percent = min_match_percent/100
    #Dividing the tolerance by 2 since we use the SUM matrix for division
    #when we calculate the relative mass difference: np.divide(product_matrix,sum_matrix)
    #instead of using the average values that would have been sum_matrix/2.
    #Thus the relative delta is 2 times lower than expected after that operation
    mass_tolerance = mass_tolerance/2
    
    charges_exp = set(df[col_names['charge_colname']].to_numpy())
    exp_data_by_charge = {}
    for ch in charges_exp:
        exp_data_by_charge[ch] = df[df[col_names['charge_colname']] == ch].copy()
    
    search_results = []

    #The library items are expected to be filtered by the descending seniority,
    #so that the glycans with more peaks are searched first
    for seniority, unique, non_unique in gaglist:
        
        min_match_num = min_match_percent*(len(unique)+len(non_unique))

        identified_glycans = {}

        #Searching for matches to the unique ions of glycans first
        charges_to_survey = charges_exp.intersection(unique[0])
        for ch in charges_to_survey:

            subdf = exp_data_by_charge[ch]
            expmzs = np.array(subdf[col_names['mz_colname']])

            theo_mzs = np.array(unique[1][ch]['mz'])

            product_matrix = np.subtract.outer(theo_mzs,expmzs)

            product_matrix = np.absolute(product_matrix)

            sum_matrix = np.add.outer(theo_mzs,expmzs)

            rel_delta_matrix = np.divide(product_matrix,sum_matrix)
            rel_delta_matrix = pd.DataFrame(rel_delta_matrix)

            bool_matrix = rel_delta_matrix.le(mass_tolerance)

            matching_mz_indices = np.where(bool_matrix == True)

            if matching_mz_indices[0].shape[0] > 0:


                for theo_index, exp_index in zip(matching_mz_indices[0],
                                                 matching_mz_indices[1]):
                    exp_rt = subdf[col_names['rt_colname']].iloc[exp_index]
                    min_allowed_rt = unique[1][ch]['min_RT'].iloc[theo_index]
                    max_allowed_rt = unique[1][ch]['max_RT'].iloc[theo_index]
                    if (exp_rt<=max_allowed_rt) and (exp_rt>=min_allowed_rt):
                        gl = unique[1][ch]['glycan'].iloc[theo_index]
                        if gl not in identified_glycans.keys():
                            identified_glycans[gl] = {col_names['mz_colname']:[],
                                                      col_names['theo_mz_colname']:[],
                                                      col_names['charge_colname']:[],
                                                      'indices_exp_mz':[],
                                                      col_names['rt_colname']:[],
                                                      col_names['abund_colname']:[]}

                        identified_glycans[gl][col_names['mz_colname']].append(expmzs[exp_index])
                        identified_glycans[gl][col_names['theo_mz_colname']].append(theo_mzs[theo_index])
                        identified_glycans[gl][col_names['charge_colname']].append(ch)
                        identified_glycans[gl]['indices_exp_mz'].append(exp_index)
                        identified_glycans[gl][col_names['rt_colname']].append(exp_rt)
                        identified_glycans[gl][col_names['abund_colname']].append(subdf[col_names['abund_colname']].iloc[exp_index])

        #Searching for the non-unique signals, for the glycans that have already got unique matches
        for g in identified_glycans.keys():
            #Finding what charge states to survey
            nu_ch_to_survey = charges_exp.intersection(non_unique[g][0])
            for ch in nu_ch_to_survey:
                subdf = exp_data_by_charge[ch]
                expmzs = np.array(subdf[col_names['mz_colname']])

                theo_mzs = non_unique[g][1][ch]

                product_matrix = np.subtract.outer(theo_mzs,expmzs)

                product_matrix = np.absolute(product_matrix)

                sum_matrix = np.add.outer(theo_mzs,expmzs)

                rel_delta_matrix = np.divide(product_matrix,sum_matrix)
                rel_delta_matrix = pd.DataFrame(rel_delta_matrix)

                bool_matrix = rel_delta_matrix.le(mass_tolerance)

                matching_mz_indices = np.where(bool_matrix == True)
                if matching_mz_indices[0].shape[0] > 0:
                    for theo_index, exp_index in zip(matching_mz_indices[0],
                                                     matching_mz_indices[1]):
                        identified_glycans[g][col_names['mz_colname']].append(expmzs[exp_index])
                        identified_glycans[g][col_names['theo_mz_colname']].append(theo_mzs[theo_index])
                        identified_glycans[g][col_names['charge_colname']].append(ch)
                        identified_glycans[g]['indices_exp_mz'].append(exp_index)
                        identified_glycans[g][col_names['rt_colname']].append(subdf[col_names['rt_colname']].iloc[exp_index])
                        identified_glycans[g][col_names['abund_colname']].append(subdf[col_names['abund_colname']].iloc[exp_index])

        #Checking if the matched glycans have enough ions.
        #If so, calculate the ppm differences between the exp and theo m/zs
        #Then add the glycan IDs to the result and marks the indices that have been matched for removal from the experimental data
        matched_indices = {}
        for g in identified_glycans.keys():
            if len(identified_glycans[g][col_names['mz_colname']]) >= min_match_num:
                #Calculate the ppm mass difference between exp and theo m/zs
                exp_mz_list = identified_glycans[g][col_names['mz_colname']]
                theo_mz_list = identified_glycans[g][col_names['theo_mz_colname']]
                identified_glycans[g]['ppm_diff'] = return_ppm_diff_list(exp_mz_list, theo_mz_list)
                mzch = []
                for idx, exp_mz in enumerate(exp_mz_list):
                    ch = identified_glycans[g][col_names['charge_colname']][idx]
                    theo_mz = round(theo_mz_list[idx], 5)
                    ppm_diff = round(identified_glycans[g]['ppm_diff'][idx], 2)
                    mzch.append((exp_mz, theo_mz, ppm_diff, ch))

                search_results.append([g, mzch, identified_glycans[g][col_names['rt_colname']],
                                       identified_glycans[g][col_names['abund_colname']],
                                       int(np.sum(identified_glycans[g][col_names['abund_colname']]))])
                
                
                for ch, expind in zip(identified_glycans[g][col_names['charge_colname']],
                                      identified_glycans[g]['indices_exp_mz']):
                    if ch not in matched_indices.keys():
                        matched_indices[ch] = []
                    matched_indices[ch].append(expind)
                        
        #Removes the matched experimental data from the dataframes            
        for ch in matched_indices.keys():
            if len(matched_indices[ch]) > 0:
                data_before = exp_data_by_charge[ch].copy()
                inds = set(range(len(data_before.index)))

                for expind in matched_indices[ch]:
                    inds.discard(expind)

                exp_data_by_charge[ch] = data_before.iloc[list(inds)].copy()
          
    return search_results

def translate_lib(inlist):
    """Translates the GAG library from the simple list format
        into the comprehensive format for searching.
        The main list:
        [[seniority, unique_precursors, non_unique_precursors],...]
        where:
            seniority - number of precursor from the GAG list,
                items are sorted in descending order of seniority
            unique_precursors list:
                [set_of_charges,dict_unique]
                where:
                    set_of_charges - set of all charges for unique ions at given seniority
                    dict_unique - {charge : df, ...}
                    where:
                        charge - guess what, it is the charge
                        df - dataframe for a given charge with cols <glycan,m/z,min_RT,max_RT>
            non_unique_precursors dict:
                {glycan : [set_of_charges,dict_non_unique],...}
                where:
                    set_of_charges - set of all charges for non-unique ions for this glycan
                    dict_non_unique - {charge : list_m/z, ...}"""

    all_seniorities = list(reversed(sorted(set([x[0] for x in inlist]))))
    out_list = []
    for s in all_seniorities:
        charge_set_unique = set()
        dict_unique = {}
        dict_nonunique = {}
        for g in inlist:
            if s == g[0]:
                charge_set_nonunique = set()
                inner_dict_nonunique = {}
                for unique_ion in g[1]:
                    charge_set_unique.add(unique_ion[1])
                    if unique_ion[1] not in dict_unique.keys():
                        dict_unique[unique_ion[1]] = {'glycan':[],'mz':[],
                                                      'min_RT':[],'max_RT':[]}
                    dict_unique[unique_ion[1]]['glycan'].append(g[4])
                    dict_unique[unique_ion[1]]['mz'].append(unique_ion[0])
                    dict_unique[unique_ion[1]]['min_RT'].append(g[3][0])
                    dict_unique[unique_ion[1]]['max_RT'].append(g[3][1])
                for non_unique_ion in g[2]:
                    charge_set_nonunique.add(non_unique_ion[1])
                    if non_unique_ion[1] not in inner_dict_nonunique.keys():
                        inner_dict_nonunique[non_unique_ion[1]] = []
                    inner_dict_nonunique[non_unique_ion[1]].append(non_unique_ion[0])
                    
                dict_nonunique[g[4]] = [charge_set_nonunique,
                                        inner_dict_nonunique]

        for ch in dict_unique.keys():
            df = pd.DataFrame(dict_unique[ch],columns=['glycan','mz','min_RT','max_RT'])
            dict_unique[ch] = df

        for g in dict_nonunique.keys():
            for ch in dict_nonunique[g][1].keys():
                dict_nonunique[g][1][ch] = np.array(dict_nonunique[g][1][ch])
                    
        out_list.append([s,
                         [charge_set_unique,dict_unique],
                         dict_nonunique])

    return out_list

def search_main_single_sql(flist, gagDBFName, saveDirPath,
                           dbPath, fileTabName, featureTabName,
                           ppmTol=10, minSN=0, minAbsInt=0):
    
    start = datetime.datetime.now()
    #One version script the first file in the folder with the extension as specified in DB_FORMAT
    #You are better have only one file with that extension in the folder in order to avoid confusion
    #For example, you can use the extension .json
##    fname = glob.glob(DB_FORMAT)[0]

    massTol = ppmTol/1000000
    
    print(f'Using the database {gagDBFName}')
    
    with open(gagDBFName,'r') as f:
        gag_list = json.load(f)
    print(f'Got {len(gag_list)} glycans in the database')

    gag_lib = translate_lib(gag_list)

    baseGagListFName = Path(gagDBFName).name

    theoMZsByCharge = get_lib_MZs_by_charge(gag_list,
                                            identicalMZsThreshold=(massTol/10))

    finish = datetime.datetime.now()
    print('Database reading and interpretation took', (finish-start),'s')
    print('Start processing files', flist)

    # The list to save all the 
    resultList = []
    for fname_data in flist:

        # Set the full paths of the report files
        print('Saving reports to\n', saveDirPath)
        out_suffix = XLSX_SUFF % ((ppmTol),
                                  datetime.datetime.now().strftime("%y%m%d"),
                                  datetime.datetime.now().strftime("%H%M%S"))
        textReportFile = saveDirPath + '/' + fname_data[:-4]+ out_suffix[:-4] + 'rep'
        shortTableFile = saveDirPath + '/' + fname_data[:-4] + '_SHORT' + out_suffix
        longTableFile = saveDirPath + '/' + fname_data[:-4] + '_LONG' + out_suffix

        print('Processing file %s' % fname_data)
        df = retrieveDF_from_SQLite(fname_data, dbPath = dbPath,
                                    fileTabName = fileTabName,
                                    featureTabName = featureTabName,
                                    sqlColList = SQL_COLLIST,
                                    colDict = COL_NAMES)

        print('Number of processing threads %s\nTolerance %.f ppm'
              % (N_CORES, ppmTol))
        print('Minimum matching percent for identification %.f\nRT width %f'
              % (MIN_MATCH_PERCENT, DBSCAN_EPS))
        print('Minimum S/N %.f and absolute feature intensity %.f'
              % (minSN, minAbsInt))

        print('Loaded table', df.shape)

        #Filter the table on m/z, N of isotopes, minimum abundance and minimum S/N
        dfFiltNotMatched = filter_feature_df(df,min_mz=240,min_isotopes=2)
        dfFiltNotMatched = dfFiltNotMatched[dfFiltNotMatched[COL_NAMES['abund_colname']] >= minAbsInt]
        dfFiltNotMatched = dfFiltNotMatched[dfFiltNotMatched[COL_NAMES['sn_colname']] >= minSN]
        dfFiltNotMatched.sort_values(by=[COL_NAMES['rt_colname']], inplace=True)
        print('Filtered table before matching', dfFiltNotMatched.shape)

        df_filtered = filter_matching_MZs(dfFiltNotMatched, allMZsDict=theoMZsByCharge,
                                          colsDict=COL_NAMES, relTol=massTol)
        print('Filtered table after matching', df_filtered.shape)
        print(df_filtered.head())

        charges = df_filtered[COL_NAMES['charge_colname']].to_numpy()
        print('Detected charges:', set(charges))

        all_rts = df_filtered[COL_NAMES['rt_colname']].to_numpy()
        print('RT array has length', len(all_rts))
        print('First 10 RTs are', all_rts[:10])

##        label_list = rt_clustering_meanshift(df_filtered, BANDWIDTH)
        label_list = rt_clustering_dbscan(df_filtered, DBSCAN_EPS)
        maxSpan, medianSpan, avgSpan = max_retTimeSpan(df_filtered)
##        df_filtered.to_excel(fname_data[:-5] + '_clust_RTs_DBSCAN.xlsx')

##        search_gags_v2(df_filtered, gag_lib, COL_NAMES,
##                       mass_tolerance=MASS_TOL, min_match_percent=MIN_MATCH_PERCENT)

        start = datetime.datetime.now()

        reslist = search_gags_ml_v2(df_filtered, gag_lib, massTol,
                                    MIN_MATCH_PERCENT,COL_NAMES)

        # Transform the list of lists into one result dictionary, where the keys are the glycans
        res_dict = {}
        if len(reslist) > 0:
            for r in reslist:
                if r[0] not in res_dict.keys():
                    res_dict[r[0]] = []
                res_dict[r[0]].append(r[1:])

        finish = datetime.datetime.now()
        print('Searching took', (finish-start),'s')
        
        
        with open(textReportFile,'w') as report_fh:
            report_fh.write('Using the database %s\nProcessing file %s\n'
                            % (baseGagListFName,fname_data))
            report_fh.write('Tolerance %.f ppm\n'
                            % (ppmTol,))
            report_fh.write('Minimum matching percent for identification %.f\nRT width %f\n'
                            % (MIN_MATCH_PERCENT, DBSCAN_EPS))
            report_fh.write('Minimum absolute peak intensity %.f and minimum S/N %.f\n'
                            % (minAbsInt, minSN))
            report_fh.write('Detected charges:' + str(set(charges)) + '\n')
            report_fh.write('Searching took ' + str((finish-start)) + ' s\n\n')
            
            for k in sorted(list(res_dict.keys())):
                report_fh.write('\nComposition %s found in:\n' % k)
                
                results = res_dict[k]
                for result_line in results:

                    report_fh.write('The following signals were found as a cluster\nm/z and z values '+
                                    str(result_line[0]) + '\n')
                    report_fh.write('retention times '+str(result_line[1])+'\n')
                    report_fh.write('abundances '+str(result_line[2])+'\n')
                    report_fh.write('sum abundance '+str(int(result_line[3]))+'\n\n')
                    
            print('\n\n\n\n\n')

        short_out_df, long_out_df = get_result_DFs(res_dict, minInt = minAbsInt)

        short_out_df.to_excel(shortTableFile)
        long_out_df.to_excel(longTableFile)

        resultList.append((fname_data, short_out_df, long_out_df, res_dict, df_filtered))

    return resultList
        
    


