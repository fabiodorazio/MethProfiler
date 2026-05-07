import xml.etree.ElementTree as ET
import pandas as pd

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

    if channel is not None:
        # Prefer Characteristics tag="source tissue"
        for char in channel.findall("geo:Characteristics", ns):
            if char.attrib.get("tag", "").lower() == "source tissue":
                tissue_source = char.text.strip() if char.text else None
                break

        # Fallback to <Source>
        if tissue_source is None:
            source = channel.find("geo:Source", ns)
            if source is not None and source.text:
                tissue_source = source.text.strip()

    rows.append({
        "sampleId": sample_id,
        "tissue_source": tissue_source
    })

df = pd.DataFrame(rows)
df.to_csv("/Users/fdorazio/Desktop/Projects/Metprofiler/assets/samplesheet_processed.csv")
print(df.head())
