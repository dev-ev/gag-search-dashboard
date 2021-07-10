import sqlite3

def create_feature_table(sql_cursor, featureTabName, fileTabName):

    table_command = (''' CREATE TABLE IF NOT EXISTS ''' + featureTabName +
                     ''' (
                                    feature_id INTEGER PRIMARY KEY,
                                    raw_file_id INTEGER NOT NULL,
                                    feature_mz REAL NOT NULL,
                                    charge INTEGER,
                                    max_sn REAL,
                                    apex_rt_min REAL,
                                    left_rt_min REAL,
                                    right_rt_min REAL,
                                    isotopes INTEGER,
                                    abundance INTEGER,
                                    comment TEXT,
                                    FOREIGN KEY (raw_file_id)
                                        REFERENCES ''' + fileTabName + ''' (raw_file_id)
                                            ON DELETE RESTRICT
                                )''')

    sql_cursor.execute(table_command)

    return True

def create_file_table(sql_cursor, fileTabName):

    table_command = (''' CREATE TABLE IF NOT EXISTS ''' + fileTabName +
                     ''' (
                                    raw_file_id INTEGER PRIMARY KEY,
                                    filename TEXT NOT NULL UNIQUE,
                                    file_date TEXT,
                                    instrument TEXT,
                                    instr_method TEXT,
                                    method_length_min REAL,
                                    comment_from_sequence TEXT,
                                    comment TEXT
                                )''')

    sql_cursor.execute(table_command)

    return True

if __name__ == '__main__':

    renameDict = {'Avg. m/z [Da]': 'feature_mz',
                  'Charge': 'charge',
                  'Max. Apex S/N': 'max_sn',
                  'Avg. Apex RT [min]': 'apex_rt_min',
                  'Left RT [min]': 'left_rt_min',
                  'Right RT [min]': 'right_rt_min',
                  '# Isotopes': 'isotopes',
                  'Abundances (per File): F50': 'abundance'}

    db_file = 'gagFeaturesDB_v210710.db'
    fileTabName = 'raw_files'
    featureTabName = 'lcms_features'

    conn = sqlite3.connect(db_file)
    conn.execute("PRAGMA foreign_keys = ON")
    cur = conn.cursor()

    create_file_table(cur, fileTabName)
    create_feature_table(cur, featureTabName, fileTabName)
    
    conn.commit()
    conn.close()

    
