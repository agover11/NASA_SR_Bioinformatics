#!/bin/bash

S3_BASE="s3://nasa-osdr/OSD-404/version-6/scRNA_Seq"
OUTDIR="./GLDS404_STARsolo"
mkdir -p "$OUTDIR"

SAMPLES=(
  "FLT_LAR_YNG_FY9" "FLT_LAR_YNG_FY11" "FLT_LAR_YNG_FY14" "FLT_LAR_YNG_FY15"
  "FLT_LAR_OLD_FO4" "FLT_LAR_OLD_FO16" "FLT_LAR_OLD_FO19" "FLT_LAR_OLD_FO20"
  "GC_LAR_YNG_GY4" "GC_LAR_YNG_GY6" "GC_LAR_YNG_GY9" "GC_LAR_YNG_GY10"
  "GC_LAR_OLD_GO14" "GC_LAR_OLD_GO16" "GC_LAR_OLD_GO19" "GC_LAR_OLD_GO20"
)

for SAMPLE in "${SAMPLES[@]}"; do
  for TYPE in GeneFull Velocyto; do
    FNAME="GLDS-404_scRNA_Seq_RRRM2_PBMC_${SAMPLE}_${TYPE}.zip"
    echo "Downloading $FNAME..."
    aws s3 cp "${S3_BASE}/${FNAME}" "${OUTDIR}/${FNAME}" --no-sign-request
    echo "Unzipping $FNAME..."
    unzip -q "${OUTDIR}/${FNAME}" -d "${OUTDIR}/${SAMPLE}_${TYPE}"
    rm "${OUTDIR}/${FNAME}"
  done
done

echo "Done!"
