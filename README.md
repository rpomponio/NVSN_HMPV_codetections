# HMPV co-detections with other respiratory viruses, 2016-2026

Project lead: Anna Wang-Erickson
Proposal date: 3/26/2026

# Proposed Timeline

* Q3 2026: Clinical data anlysis complete
* Q4 2026: Manuscript drafted
* Q2 2027: Manuscript publication

# Data Source (from LG at CDC)

Please see the ShareFile link below for the updated dataset and DD
(with new additions highlighted), that include data through end of
April 2026 to account for most of the data lag for testing and chart
review (though these will still be incomplete for non-closed out data). 

<https://centersfordiseasecontrol.sharefile.com/d-sf825f47e261f4ca8aaf2df2ff352a1fb>

Most of the testing result variables are calculated variables, which
combine results from multiple fields. We do not have the same calculated
variables for CT value, so I have included the individual CT results. For
example, the c_rsv_result variable combines the variables trsv (untyped),
trsva (RSV-A), and trsvb (RSV-B) and I have provided the CT value variables
for each of these. Let me know if you have any questions or if I can assist
with interpretation. For the most part, CT values are available for all
pathogens for the same four sites (Houston, Pittsburgh, Rochester, and
Vanderbilt). The trsv and HCoV CT variables also have KC, I believe, but
I would double check.

The co-detections variables include all viral and bacterial co-detections.
Specifically, SARS-CoV-2, influenza, HCoV, PIV, RV/EV,  RSV, adenovirus,
mycoplasma, legionella, bordetella pertussis, bocavirus, and chlamydophila
pneumoniae. If you are interested in a variable for just certain co-detections,
I can create that for you.

# Directory Structure

```
./
│   .gitignore
│   HMPV Co-detections.Rproj
│   HMPV Codetections Request Form_2026_04_06.docx
│   ...
├───Archive
│   ...
└───Data
        nvsn_PITT_ANNA_HMPV_DD_6.9.26.xlsx
        Pitt_Anna_HMPV_JUN26.csv
```