# Single_Sample_vcf_PRS_and_Calibration
This project provides a set of R and Python scripts for calculating Polygenic Risk Scores (PRS) from single-sample imputed VCF files and performing ancestry-adjusted PRS calibration.



    ######PRS Calculation: Use either PRS.R (R) or PRS_Python.py (Python) to compute PRS.

    Input files:

        Per-chromosome VCF files (one sample each).

        A PRS weight file, containing SNP IDs and effect sizes.

    Ensure:

        SNP IDs in the weight file match those in the VCF files.

        Quality control (QC) steps are performed, including:

            Checking SNP missingness.

            Removing problematic variants.

    Output: Raw PRS value for the individual.
    
    #####PRS Calibration (Ancestry Adjustment)

    Ancestry-adjusted PRS is calculated to account for population structure. The adjustment method is described in our previous publication:

    DOI: 10.1038/s41467-023-41819-0
