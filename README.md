# gag-search-dashboard
## Script for quantification of glucosaminoglycans (GAGs) based on Proteome Discoverer output, with web-interface powered by Bokeh

The scripts work well in Python 3.7 and 3.8. The script for Proteome Discoverer integration has only been tested in Windows, but the search and the web-interface have been working both in Windows 10 and Ubuntu 20.04.

### Processing steps

The script have been developed for precursor quantification on the chondroitinase-digested glucosaminoglycanes. This version is specifically tailored for the GAG workflow with di-n-butylamine in the mobile phases, as described in [the publication by Persson *et al.*](https://www.nature.com/articles/s41598-020-60526-0).<br>
LC-MS peaks are quantified using the Minora Feature Detector node in Proteome Discoverer (PD) in this workflow. Why PD? Because it was available, and I was familiar with it's performance and capability. Unfortunately, PD is a proprietary commercial software, so I do not combine the workflow into one integrated distribution. On the contrary, the processing is separated into several steps, with only the first step being dependent on PD. The steps are as follows:

1. Create the library of theoretical GAG forms that will be matched against the experimental data. I decided to start with the human readable tab-separated files that are defined by glycobiologists according to their expectations:
```
dp2S0_int	0	1	0	379.1115	0	1000
dp2S0NS1_int	1	1	0	417.0577	0	1000
...
```
Where the columns are as follows:
```
      1) unique name of the glycan
      2) number of sulfates
      3) number of hexuronic acid residues
      4) number of NeuAc residues
      5) neutral monoisotopic mass of the molecule
      6) minimal retention time in minutes
      7) maximal retention time in minutes
```

I added the retention time columns for a potential future use, if there would be an extensive library of glycans with retention times. I haven't really used this option down the line.<br>
The script *gag-lib-gen* is then applied to create a json library from each tab file. For each of the molecules in the tab file, the script creates the chemically intact forms, as well as sulfate loss forms (prevalent for sulfated sugars) and amine adduct forms (di-n-butylamine was used as an additive to the mobile phases). The charge states are decided based on the detectable *m/z* values and on empirical rules, depending on the number of sulfate and carboxyl groups. Examples of the tab files and json libraries can be located in the folder *lib-generator*.<br>

2. Create an SQLite database for storing the quantified chromatographic information using the *gagFeatureDB_create* script. I used a single SQLite database for all the GAG files regardless of the project. The format allows for a simple and fast access to each of the pre-processed files. The database contains information about the intensity, *m/z*, charge, number of detected isotopes and retention times for the LC-MS peaks for all the quantified LC-MS files.<br>
3. LC-MS raw files are processed in PD, Minora Feature Detection node quantifies the LC-MS features. The list of the quantified features is then saved into the SQLite database via the Scripting Node with the *gag_Minora-to-SQlite* script. An example of a database containing the quantification results for a few files can be found in this repository (*gagFeaturesDB_v210710.db*).<br>
4. Provide the paths to the json libraries, LC-MS peak database and to the output folder by editing the *config.ini* file. By default, the configuration file is in the main application folder alongside the Python code.
5. When the configuration has been set up, we can open the terminal, change the directory to the folder with out *Bokeh* app and launch the web-interface: 
```
bokeh serve --show gag-search-web
```
The dashboard will open in the default web browser. If it does not, check the IP address and the port that are currently used by bokeh in the console.<br>

<img src="https://github.com/dev-ev/gag-search-dashboard/blob/main/images/gag_search_overview1.png" alt="drawing" width="400"/>

The dropdown menu allows to select the GAG databases

6. Results

### Search algorithm

