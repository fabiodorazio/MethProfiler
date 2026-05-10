import xml.etree.ElementTree as ET
import pandas as pd

import utils

def clean_classify_samplesheet(samplesheet, control_threshold=1):
    '''
    removes samples that are labelled Exclude
    classifies  
    '''
    # remove rows to exclude
    cleaned_samplesheet = samplesheet[samplesheet['disease_status'] != 'Exclude']
    cleaned_samplesheet = cleaned_samplesheet[cleaned_samplesheet['braak_stage'] != 'Exclude']

    # convert braak_stage to numeric after removing strings
    cleaned_samplesheet['braak_stage'] = pd.to_numeric(cleaned_samplesheet['braak_stage'], errors='coerce')

    # binary classification of disease
    cleaned_samplesheet['Condition'] = cleaned_samplesheet['braak_stage'].apply(lambda x: 'Disease' if x >= control_threshold else 'Control')

    return cleaned_samplesheet


xml_file = "/Users/fdorazio/Desktop/Projects/Metprofiler/assets/GSE59685_family.xml"

# Parse XML
tree = ET.parse(xml_file)
root = tree.getroot()

# GEO MINiML namespace
ns = {"geo": "http://www.ncbi.nlm.nih.gov/geo/info/MINiML"}

rows = []

for sample in root.findall("geo:Sample", ns):
    sample_id = sample.attrib.get("iid")

    channel = sample.find("geo:Channel", ns)

    tissue_source = None
    braak_stage = None
    disease_status = None
    sex = None
    age = None

    if channel is not None:
        # Prefer Characteristics tag="source tissue"
        for char in channel.findall("geo:Characteristics", ns):
            tag = char.attrib.get("tag", "").lower()
            value = char.text.strip() if char.text else None

            # source tissue
            if tag == "source tissue":
                tissue_source = value

            # braak stage
            elif tag == "braak.stage":
                braak_stage = value
            # Alzheimer disease status
            elif tag == "ad.disease.status":
                disease_status = value
            # Sex
            elif tag == "sex":
                sex = value
            # age: using age_brain
            elif tag == "age.brain":
                age = value

        # Fallback to <Source>
        if tissue_source is None:
            source = channel.find("geo:Source", ns)
            if source is not None and source.text:
                tissue_source = source.text.strip()
        

    rows.append({
        "sampleId": sample_id,
        "tissue_source": tissue_source,
        "braak_stage": braak_stage,
        "disease_status":disease_status,
        "Sex": sex,
        "age": age
    })

df = pd.DataFrame(rows)

# apply filtering and binarization
df_final = clean_classify_samplesheet(df, 3)
print(df_final)

# assets/ must exist
df_final.to_csv("/Users/fdorazio/Desktop/Projects/Metprofiler/assets/samplesheet_processed.csv")
