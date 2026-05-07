### Metprofiler

Disclaimer: 'This code was handwritten. LLMs were used for consultation on programming logic, documentation lookup, and API/function usage'

## Structure
```
main.py <── DMA.py
^
├── samplesheet.csv
└── betas.csv
```

## Overarching logic
- Betas have a 3-rows header indicating the corresponding GSMs. For standardisation, the header is cut off and the beta values Samplesheet instead is used to subset the relevant barcodes.
- Due to the size of the table, the pipeline has the option to downsample the input beta table to make it more manageble during development. Default = 'Off'
- After loading the input files, the code asserts for missing Sample_Barcodes in the column names of the beta table
