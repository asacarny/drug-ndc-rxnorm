# Opioids, Benzodiazepines, and Other Drug NDCs from RxNorm
This repository contains code to create a list of NDCs (National Drug Codes) using [RxNorm](https://www.nlm.nih.gov/research/umls/rxnorm/index.html), the public database of drug nomenclature. I wrote this code to to compile NDCs for opioids and other related controlled substances, but it can work with any drug tracked in RxNorm.

One goal of this code is to replace the CDC's now-discontinued file of opioid NDCs. Unfortunately, the CDC stopped supplying this list with the publication of its 2022 clinical practice guidelines. [The discontinuation message](https://archive.cdc.gov/www_cdc_gov/opioids/data-resources/index.html) referred users to RxNorm but did not provide any code or guidance on how to build a comparable dataset. Since RxNorm is tricky to work with, I hope this code will save people time if they are trying to take the CDC's suggestion.

I would love to post the NDC list data online, but I believe the [RxNorm license agreement](https://uts.nlm.nih.gov/uts/assets/LicenseAgreement.pdf) prohibits it. So at least for now, I'm only distributing the code.

# Installation

The installation process should be straightforward. You will need a recent version of RStudio as well as a small set of R packages installed (`tidyverse`, `jsonlite`, `curl`, `pbapply`, and `haven`).

Directions:

1. Download this repository from Github either by clicking Code -> Download ZIP on this page, or by issuing the command from your terminal `git clone https://github.com/asacarny/drug-ndc-rxnorm.git`.
2. Navigate to the repository and open `drug-ndc-rxnorm.Rproj` in RStudio.
4. Open the main R script, `drug-ndc-rxnorm.R` and run it.

The code will produce two files, `ndc-list-final.Rdata` and `ndc-list-final.dta`, the NDC list in R and Stata format, respectively.

# Drug List

The code compiles a list of drug NDCs related to my own research, but it may meet the needs of other users. I have aimed to capture the major opioids, gabapentinoids, benzodiazepines, and stimulants filled in outpatient pharmacies and tracked in databases like PDMPs and pharmacy claims. Buprenorphine is classified separately as an OUD drug while methadone is retained as an opioid -- you may want to reclassify these drugs depending on what is tracked in your data. (In our data, we only observe outpatient pharmacy dispenses for methadone, where it is used for pain treatment.)

To bring more drugs into your NDC list, add rows to the `drug_list` tibble at the top of the code. The `drug` column is the drug's ingredient. To determine the right name of the ingredient to put here, you can search [RxNav](https://mor.nlm.nih.gov/RxNav/) and review what's listed in the `Ingredient` box.

# Equivalency Calculations

The `drug_list` tibble contains a column, `equiv`, that indicates the conversion from 1mg of the drug to a base drug in the same class. The base drugs are morphine (opioids), gabapentin (gabapentinoids), and diazepam (benzodiazepines). This field is key for users hoping to calculate e.g. morphine milligram equivalents (MME) or morphine equivalent daily dose (MEDD) and analogous fields for gabapentinoids or benzodiazepines.

To calculate MME or MEDD (and their analogs for other classes), I run: 

```
mme=equiv*strength_per_unit*quantity_dispensed
medd=mme/days_supply
```

where `quantity_dispensed` and `days_supply` come from claims or PDMP data.

For opioids, I pull `equiv` from the conversion factors in the 2020 CDC Opioid NDC list. Conversions for fentanyl are specific to its dose form and are done later in the code, not in `drug_list`. Note that the 2022 CDC clinical practice guidelines have different conversion factors than the ones used here for some opioids. I expect to update the conversion factors in the future.

For benzodiazepines, there is no generally recognized set of conversion factors. The values in this code are the average of the conversion factors from three textbooks. If the textbook gave a range, we took its midpoint. For clobazepam and halazepam, which were not mentioned in the textbooks, we searched the academic literature. See `equivalency-table.xlsx`. Note that we developed these conversions several years ago and did so without extensive subject matter expertise.

To convert pregabalin to gabapentin, we used the equivalence ratio mentioned in [an academic study](https://www.sciencedirect.com/science/article/pii/S2451865421001526).

# Notes

* I am really not sure how complete the RxNorm NDC list is for older NDCs!
* Products can contain multiple ingredients, creating a challenge for how to represent them in this data. In the current drug list, this is only an issue for Adderall and its generics, which contain amphetamine *and* dextroamphetamine. We give these records an `ingredient` of amphetamine. Note that `strength_per_unit` for these NDCs only refers to amphetamine, it ignores the dextroamphetamine in the product.
* Products can also contain just one ingredient, but multiple "precise ingredients". For instance, a product could contain amphetamine aspartate and amphetamine sulfate. In this situation, we add together the `strength_per_unit` of the precise ingredients and keep their names in the `precise_ingredient` field separated by slashes.
* In RxNorm, certain sprays seem to have strengths given per actuation, but PDMP or claims data might give quantity dispensed in terms of mL. For products that are packs, it's pretty unclear what strength means in the RxNorm data. In these cases, I recode `equiv` and `strength_per_unit` to missing.
* For products delivered intravenuously, I recode `equiv` to missing because the conversion factors may not be correct for these delivery routes.

# Change log

June 27, 2025

* Initial release