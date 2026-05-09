### Metprofiler
Metprofiler is a python pipeline built to automate the analysis of DNA methylation from this [dataset](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE59685)

```
Disclaimer: 'This code was handwritten. LLMs were used for consultation on programming logic, documentation lookup, and API/function usage'

```

## Structure
```
[main.py](./bin/main.py) <── [DMA.py](./bin/DMA.py)
^
├── samplesheet.csv <── [samplesheet.workout.py](./bin/samplesheet_workout.py)
└── betas.csv
```

## Overarching logic
1) [samplesheet.workout.py](./bin/samplesheet_workout.py)
- reads the xml file provided by the authors and 
2) [DMA.py](./bin/DMA.py)
- Betas have a 3-rows header indicating the corresponding GSMs. For standardisation, the header is cut off and the beta values Samplesheet instead is used to subset the relevant barcodes.
- Due to the size of the table, the pipeline has the option to downsample the input beta table to make it more manageble during development. Default = 'Off'
- After loading the input files, the code asserts for missing Sample_Barcodes in the column names of the beta table
