# gag-search-dashboard
## Script for quantification of glucosaminoglycans (GAGs) based on Proteome Discoverer output, with web-interface powered by Bokeh

The scripts work well in Python 3.7 and 3.8. The script for Proteome Discoverer intgration has only been tested in Windows, but the search and the web-interface have been working both in Windows 10 and Ubuntu 20.04.

Developed for precursor quantification on the chondroitinase-digested glucosaminoglycanes. LC-MS peaks are quantified using the Minora Feature Detector node in Proteome Discoverer (PD). Why PD? Because it was available, and I was familiar with it's performance and capabilities. Unfortunately, PD is a commercial software, so I do not combine the processing into one integrated workflow. On the contrary, the processing is separated into several steps, with only the first step being dependent on PD. The steps are as follows:

1. Create an SQLite database for storing the quantified chromatographic information using the *gagFeatureDB_create* script. I used a single SQLite database for all the GAG files regardless of the project. The format allows for a simple and fast access to each of the pre-processed files. The database contains information about the intensity, *m/z*, charge, number of detected isotopes and retention times for the LC-MS peaks for all the quantified LC-MS files.
2. LC-MS raw files are processed in PD, Minora Feature Detection node quantifies the LC-MS features. The list of the quantified features is then saved into the SQLite database via the Scripting Node with the *gag_Minora-to-SQlite* script. An example of a database containing the quantification results for a few files can be found in this repository (*gagFeaturesDB_v210710.db*).
3. When all the interesting LC-MS files have been pre-processed, we can open the terminal, change the directory to the folder with out *Bokeh* app and launch the web-interface: 
```
bokeh serve --show gag-search-web
```

4. Results
