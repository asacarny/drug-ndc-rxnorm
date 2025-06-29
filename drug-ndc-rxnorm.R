library(tidyverse)
library(jsonlite)
library(curl)
library(pbapply)
library(haven)

source("rxnorm-utility.R")

drug_list <- tribble(
  ~class,   ~drug,            ~equiv,
  "oud"   , "buprenorphine",  NA,
  "opioid",	"butorphanol",    7,
  "opioid",	"codeine",        0.15,
  "opioid",	"dihydrocodeine", 0.25,
  "opioid",	"fentanyl",       NA,
  "opioid",	"hydrocodone",    1,
  "opioid",	"benzhydrocodone",1.2,
  "opioid",	"hydromorphone",  4,
  "opioid",	"levomethadyl",   8,
  "opioid",	"levorphanol",    11,
  "opioid",	"meperidine",     0.1,
  "opioid",	"methadone",      3,
  "opioid",	"morphine",       1,
  "opioid",	"opium",          1,
  "opioid",	"oxycodone",      1.5,
  "opioid",	"oxymorphone",    3,
  "opioid",	"pentazocine",    0.37,
  "opioid",	"propoxyphene",   0.23,
  "opioid",	"tapentadol",     0.4,
  "opioid",	"tramadol",       0.1,
  "gaba"  , "gabapentin",     1,
  "gaba"  , "pregabalin",     6, # https://www.sciencedirect.com/science/article/pii/S2451865421001526
  "benzo" , "alprazolam",     15,
  "benzo" , "chlordiazepoxide", 0.4,
  "benzo" , "clobazam",       0.5,
  "benzo" , "clonazepam",     14.16666667,
  "benzo" , "clorazepate",    0.888888889,
  "benzo" , "diazepam",       1,
  "benzo" , "estazolam",      8.333333333,
  "benzo" , "flunitrazepam",  10,
  "benzo" , "flurazepam",     0.5,
  "benzo" , "halazepam",      0.25,
  "benzo" , "lorazepam",      6.666666667,
  "benzo" , "midazolam",      0.666666667,
  "benzo" , "oxazepam",       0.555555556,
  "benzo" , "prazepam",       0.5,
  "benzo" , "quazepam",       0.611111111,
  "benzo" , "temazepam",      0.555555556,
  "benzo" , "triazolam",      33.33333333,
  "stim"  , "methylphenidate",    NA,
  "stim"  , "dexmethylphenidate", NA,
  "stim"  , "amphetamine",        NA,
  "stim"  , "dextroamphetamine",  NA,
  "stim"  , "lisdexamfetamine",   NA,
  "stim"  , "methamphetamine",    NA
)

print("> Fetching rxcui of ingredients.")
ingredients <-
  drug_list$drug |>
  pblapply(findByString) |>
  bind_rows() |>
  inner_join(
    drug_list,
    join_by(drug),
    unmatched="error",
    relationship="one-to-one"
  )

# get the products (SCDs, GPCKs, SBDs, and BPCKs) of all the INs
print("> Fetching products of each ingredient.")
product <- ingredients$rxcui |>
  pblapply(getRelated,tty=c("SCD","GPCK","SBD","BPCK"),src_name="ingredient") |>
  bind_rows() |>
  rename(rxcui_product=rxcui_scd_gpck_sbd_bpck)

# Adderall and its generic equivalents come under ingredient amphetamine and
# and ingredient dextroamphetamine because the drugs contain both.
# These are duplicate records in the product table which cause problems
# later. I'll handle by removing the dex observations

### WARNING WARNING ###
# This means the strength_per_unit for adderall and its generics only
# counts amphetamine

duplicates <- product |>
  group_by(rxcui_product) |>
  filter(n()>1) |>
  ungroup()

if (nrow(duplicates)>0) {
  print("duplicate product rows (these should all be adderall and generics):")
  print(duplicates |> nrow())
  
  if (duplicates |> group_by(rxcui_product) |> filter(n()>2) |> ungroup() |> nrow()) {
    stop("there were duplicated products with more than 2 rows")
  }
  
  n_non_add <- duplicates |>
    filter(!str_detect(name,"amphetamine aspartate [0-9.]+ MG / amphetamine sulfate [0-9.]+ MG / dextroamphetamine saccharate [0-9.]+ MG / dextroamphetamine sulfate [0-9.]+ MG")) |>
    nrow()
  if (n_non_add > 0) {
    stop("there were duplicate products other than adderall and its generics")
  }
  
  # bring in ingredient name
  duplicates <- duplicates |>
    inner_join(
      ingredients |> select(drug,rxcui_ingredient=rxcui),
      join_by(rxcui_ingredient),
      unmatched=c("error","drop"),
      relationship="many-to-one"
    )
  
  # confirm each duplicated drug has amphetamine and dex row
  n_with_both <- duplicates |>
    group_by(rxcui_product) |>
    filter(any(drug=="amphetamine") & any(drug=="dextroamphetamine")) |>
    ungroup() |>
    nrow()
  if (n_with_both != nrow(duplicates)) {
    stop("each duplicated product did not have both an amphetamine and dextroamphetamine row")
  }
  
  # just the dex records
  dex_records <- duplicates |>
    filter(drug=="dextroamphetamine") |>
    select(rxcui_product,rxcui_ingredient)
  print("dex records to remove from product list, keeping amphetamine version:")
  print(dex_records |> nrow())
  
  # remove them from the product list
  print("removing")
  product <- product |>
    anti_join(
      dex_records,
      join_by(rxcui_product,rxcui_ingredient)
    )
  
}

# get details on the products
print("> Fetching details of each product.")
product_detail <-
#  lapply(X=scd_gpck$rxcui_scd_gpck,FUN=getHistoryStatus,rxcui_base=rxcui_ingredient) |>
  pbmapply(
    FUN=getHistoryStatus,
    rxcui_target=product$rxcui_product,
    rxcui_base=product$rxcui_ingredient,
    SIMPLIFY=FALSE
  ) |>
  bind_rows() |>
  rename(rxcui_product=rxcui_target)

# dose forms to flag as IV drugs - we'll set equiv to NA for them
dose_forms_iv_list <- c(
  "Cartridge",
  "Injectable Solution",
  "Injectable Suspension",
  "Injection",
  "Prefilled Syringe"
)

# add flag as IV
product_detail <- product_detail |>
  mutate(dose_form_iv = dose_form %in% dose_forms_iv_list)

# find all NDCs of the products

# get NDC list for every product
print("Fetching NDCs of each product.")
ndc_hist <-
  product_detail |>
#  filter(dose_form_iv==FALSE) |>
  pull(rxcui_product) |>
  pblapply(getHistNDC) |>
  bind_rows() |>
  rename(rxcui_product=rxcui)

# now merge together with product details
ndc_hist_detail <-
  ndc_hist |>
  inner_join(
    product_detail,
    join_by(rxcui_product),
    unmatched=c("error","drop"),
    relationship="many-to-one"
  )

# it turns out multiple product rxcui's often point to the same ndc
# in practice this tends to be harmless bc the key parameters for us are
# constant across the products

# let's get the data down to the unique NDC level
print("NDCs with a unique row in ndc_hist_detail")
ndc_unique_row <- 
  ndc_hist_detail |>
  group_by(ndc) |>
  filter(n()==1) |>
  ungroup()
print(ndc_unique_row |> nrow())

print("NDCs with multiple rows in ndc_hist_detail")
ndc_multiple_row <- 
  ndc_hist_detail |>
  group_by(ndc) |>
  filter(n()>1) |>
  ungroup()
print(ndc_multiple_row |> group_by(ndc) |> summarize() |> nrow())

# most recent row for each NDC
ndc_multiple_recent <-
  ndc_multiple_row |>
  group_by(ndc) |>
  filter(end_date==max(end_date)) |>
  ungroup()

# NDCs that are now unique (good to use)
ndc_multiple_recent_unique <-
  ndc_multiple_recent |>
  group_by(ndc) |>
  filter(n()==1) |>
  ungroup()

# NDCs that are still not unique (need more work)
ndc_multiple_recent_multiple <-
  ndc_multiple_recent |>
  group_by(ndc) |>
  filter(n()>1) |>
  ungroup()

print("NDCs with multiple rows even after limiting to most recent end dates")
print(ndc_multiple_recent_multiple |> group_by(ndc) |> summarize() |> nrow())

# pick the first row with the combo of (ndc,drug,uom,strength,doseform)
ndc_multiple_recent_multiple_fixed <-
  ndc_multiple_recent_multiple |>
  group_by(ndc,drug,uom,strength_per_unit,dose_form) |>
  arrange(pick(!c(ndc,drug,uom,strength_per_unit,dose_form))) |>
  slice(1) |>
  ungroup()

# this is a hack to deal with some old morphine records
# the only non-unique NDC rows should be these records...
if (
    ndc_multiple_recent_multiple_fixed |>
      group_by(ndc) |>
      filter( !( n()==1 | (n()==2 & drug=="morphine" & dose_form %in% c("Injection","Cartridge")) ) ) |>
      nrow()
) {
  stop("NDCs weren't unique outside of known issue with morphine")
}

# for the morphine records, just pick the injection ones over the cartridge ones
ndc_multiple_recent_multiple_fixed <-
  ndc_multiple_recent_multiple_fixed |>
  group_by(ndc) |>
  filter( n()==1 | (n()==2 & drug=="morphine" & dose_form=="Injection") ) |>
  ungroup()


if (ndc_multiple_recent_multiple_fixed |> group_by(ndc) |> filter(n()>1) |> nrow()) {
  stop("NDCs weren't unique")
}

ndc_list_final <-
  bind_rows( # bring together distinct NDCs
    ndc_unique_row,
    ndc_multiple_recent_unique,
    ndc_multiple_recent_multiple_fixed
  ) |>
  arrange(drug,active_ingredient,ndc) |>
  inner_join( # add rxcui of ingredient
    product |> select(starts_with("rxcui")),
    join_by(rxcui_product),
    unmatched=c("error","drop"),
    relationship("many-to-one")
  ) |>
  inner_join( # add drug class and equivalency
    ingredients |> select(!drug),
    join_by(rxcui_ingredient==rxcui),
    unmatched=c("error","drop"),
    relationship("many-to-one")
  ) |>
  mutate( # update equivalencies for fentanyl
    equiv = case_when(
      drug=="fentanyl" & 
        dose_form %in% c("Buccal Tablet","Sublingual Tablet","Oral Lozenge") &
        str_detect(uom,"^MG")
        ~ 0.13*1000,
      drug=="fentanyl" &
        dose_form %in% c("Mucosal Spray","Pack") & # Packs are actually mucosal sprays
        str_detect(uom,"^MG")
        ~ 0.18*1000,
      drug=="fentanyl" & 
        dose_form=="Metered Dose Nasal Spray" &
        str_detect(uom,"^MG")
        ~ 0.16*1000,
      drug=="fentanyl" &
        dose_form=="Transdermal System" &
        str_detect(uom,"^MG")
        ~ 7.2*1000,
      .default = equiv
    )
  )

if (ndc_list_final |> filter(drug=="fentanyl" & is.na(equiv) & !dose_form_iv) |> nrow()) {
  stop("Didn't fill in equivalency for all (non-IV) fentanyl")
}    

if (ndc_list_final |> group_by(ndc) |> filter(n()>1) |> nrow()) {
  stop("NDCs weren't unique")
}

# equivalency for IV drugs is not always clear. let's not even try
# set equiv to NA for IV drugs.
ndc_list_final <- ndc_list_final |>
  mutate(equiv = replace(equiv,dose_form_iv,NA))

#### WARNING WARNING ###
# RXNORM SEEMS TO CLASSIFY SPRAYS BY MG PER ACTUATION
# THUS THE DOSE STRENGTH IS MG PER ACTUATION, NOT MG PER ML
# HOWEVER THE QUANTITY DISPENSED FOR SUCH SPRAYS TENDS TO BE IN ML

print("NDCs in which dose strength given in MG/ACTUAT (setting strength per unit & equiv to NA)")
print(ndc_list_final|>filter(uom=="MG/ACTUAT")|>nrow())

ndc_list_final <- ndc_list_final |>
  mutate(
    equiv = replace(equiv,uom=="MG/ACTUAT",NA),
    strength_per_unit = replace(strength_per_unit,uom=="MG/ACTUAT",NA)
  )

#### WARNING WARNING ####
# I WONDER IF WE ARE NOT HANDLING PACKS CORRECTLY
# WHAT WOULD QUANTITY DISPENSED BE FOR A PACK?
# SEE NDC 52427081599
# let's set strength and equiv to NA for packs
ndc_list_final <- ndc_list_final |>
  mutate(
    equiv = replace(equiv,dose_form=="Pack",NA),
    strength_per_unit = replace(strength_per_unit,dose_form=="Pack",NA)
  )

ndc_list_final |>
  write_dta("ndc-list-final.dta")

ndc_list_final |>
  save(file="ndc-list-final.Rdata")
