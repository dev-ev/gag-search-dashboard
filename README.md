# gag-search-dashboard
Script for quantification of glucosaminoglycans (GAGs) based on Proteome Discoverer output, with web-interface powered by Bokeh

The scripts work well in Python 3.7 and 3.8. The script for Proteome Discoverer intgration is only tested in Windows, but the search and the web-interface have been working both in Windows 10 and Ubuntu 20.04.

Developed for precursor quantification on the chondroitinase-digested glucosaminoglycanes. LC-MS peaks are quantified using the Minora Feature Detector node in Proteome Discoverer (PD). Why PD? Because it was available, and I was familiar with it's performance and capabilities. Unfortunately, PD is a commercial software, so I do not combine the processing into one integrated workflow. On the contrary, the processing is separated into several steps, with only the first step being dependent on PD. The steps are as follows:
1. Create an SQLite database for storing the quantified chromatographic information using the *gagFeatureDB_create* script. I used a single SQLite database for all the GAG files regardless of the project. The format allows for a simple and fast access to each of the pre-processed files. The 
```
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
```
2. LC-MS raw files are processed in PD, Minora Feature Detection node quantifies the LC-MS features. The list of the quantified features is then saved into the SQLite database via the Scripting Node with the *gag_Minora-to-SQlite* script. 
