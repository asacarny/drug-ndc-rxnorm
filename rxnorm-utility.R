findByString <- function(str_name) {
  url <- paste0(
    "https://rxnav.nlm.nih.gov/REST/rxcui.json?name=",
    str_name,
    "&search=0"
  )
  #print(paste0("Fetching ",url))
  
  result <- fromJSON(url)
  
  if (length(result$idGroup)==0) {
    stop(paste0("Couldn't find name ",str_name))
  }
  if (length(result$idGroup)>1) {
    stop(paste0("Too many rxcui's for ",str_name))
  }
  
  tibble(drug=str_name, rxcui=as.integer(result$idGroup$rxnormId))
}

classMembers <- function(classid) {
  url <- paste0(
    "https://rxnav.nlm.nih.gov/REST/rxclass/classMembers?classId=",
    classid,
    "&relaSource=DAILYMED"
  )
  
  fromJSON(url)$drugMemberGroup$drugMember$minConcept |>
    mutate(class=classid)
}

getRelated <- function(rxcui,tty,src_name) {
  url <- paste0(
    "https://rxnav.nlm.nih.gov/REST/rxcui/",
    rxcui,
    "/related?tty=",
    paste(tty,sep="+",collapse="+")  
  )
  #print(paste0("Fetching ",url))
  
  rxcui_src_name <- paste("rxcui",src_name,sep="_")
  rxcui_tty_name <- paste(c("rxcui",tolower(tty)),sep="_",collapse="_")
  
  fromJSON(url)$relatedGroup$conceptGroup |>
    apply(1,function(v) {v$conceptProperties}) |>
    bind_rows() |>
    rename({{rxcui_tty_name}} := rxcui) |>
    mutate({{rxcui_src_name}} := rxcui, .before=1) |>
    mutate(across(starts_with("rxcui"),as.integer)) |>
    #mutate(class=class) |>
    select(-language,-suppress,-umlscui)
  
}

getHistoryStatus <- function(rxcui_target,rxcui_base) {
  #print(paste0(
  #  "Working on rxcui target ",
  #  as.character(rxcui_target),
  #  " rxcui base ",
  #  as.character(rxcui_base)
  #))
  
  url <- paste0(
    "https://rxnav.nlm.nih.gov/REST/rxcui/",
    rxcui_target,
    "/historystatus"
  )  
  
  json <- fromJSON(url)
  #str(json)
  
  if(!is.list(json$rxcuiStatusHistory)) {
    stop("Did JSON request fail?")
  }
  else if (json$rxcuiStatusHistory$metaData$status=="UNKNOWN") {
    stop("Was rxcui valid?")
  }
  
  name_target <- json$rxcuiStatusHistory$attributes$name
  branded_target <- json$rxcuiStatusHistory$attributes$isBranded
  
  ingredient_strength <-
    json$rxcuiStatusHistory$definitionalFeatures$ingredientAndStrength |>
    filter(baseRxcui==rxcui_base) |> # down to just the ingredient we care about
    mutate(
      uom = paste0( # unit of measure is EACH, or numeratorunit/denominatorunit
        numeratorUnit,
        if_else(denominatorUnit!="EACH",paste0("/",denominatorUnit),"")
      ),
      strength_per_unit = if_else( # if denominatorValue isn't 1, set strength unknown
        as.double(denominatorValue)==1,
        as.double(numeratorValue),NA_real_
      )
    ) |>
    rename(
      drug=baseName,
      active_ingredient=activeIngredientName
    ) |>
    select(drug,active_ingredient,uom,strength_per_unit)
  
  if (nrow(ingredient_strength)==0) {
    stop("Couldn't find active ingredient")
  }
  else if (nrow(ingredient_strength)>1) {
    # some rxcuis have multiple ingredient records with the ingredient
    # we care about. we simplify them by adding them up, assuming
    # the units are the same
    cat("\nToo many active ingredients, simplifying")
    print(ingredient_strength)
    ingredient_strength <- ingredient_strength |>
      group_by(drug,uom) |> # they have to be same drug & uom
      summarize(
        active_ingredient=paste(active_ingredient,collapse="/"),
        strength_per_unit=sum(strength_per_unit),
        .groups="drop"
      ) |>
      relocate(drug,active_ingredient,uom,strength_per_unit)
    print(ingredient_strength)
    
    # this should get us to one ingredient record, otherwise stop... 
    if (nrow(ingredient_strength)>1) {
      stop("Still couldn't simplify down to 1 common ingredient")
    }
  }
  
  if (nrow(json$rxcuiStatusHistory$definitionalFeatures$doseFormConcept)>1) {
    stop("Too many dose forms")
  }
  dose_form <- json$rxcuiStatusHistory$definitionalFeatures$doseFormConcept$doseFormName
  
  # now build the one-row tibble to be returned
  ingredient_strength |>
    mutate(
      rxcui_target=rxcui_target,
      full_name=name_target,
      .before=1
    ) |>
    mutate(
      #      branded=case_when(
      #        branded_target=="YES" ~ TRUE,
      #        branded_target=="NO"  ~ FALSE,
      #        .default=NA
      #      ),
      dose_form=dose_form
    )
}

getHistNDC <- function (rxcui) {
  #print(paste0("Working on rxcui ",as.character(rxcui)))
  
  url <- paste0(
    "https://rxnav.nlm.nih.gov/REST/rxcui/",
    rxcui,
    "/allhistoricalndcs.json"
  )
  
  json <- fromJSON(url)
  
  if (length(json)==0) {
    return (NULL)
  }
  
#  json$historicalNdcConcept$historicalNdcTime$ndcTime[[1]] |>
  json$historicalNdcConcept$historicalNdcTime$ndcTime |> # this is a list of data frames
    # which could be 2 elements long, with first element status=INDIRECT and 
    # second status=DIRECT. INDIRECT are NDCs associated with rxcui that were
    # remapped to this one, which I think we want to keep.
    bind_rows() |>
    mutate(rxcui=rxcui,.before=1) |>
    mutate(start_date=ym(startDate),end_date=rollforward(ym(endDate))) |>
    select(!c(startDate,endDate)) |>
    unnest(ndc) # ndc was stored as a list-column, switch to atomic vector
}

getStrength <- function(rxcui) {
  url <- paste0(
    "https://rxnav.nlm.nih.gov/REST/rxcui/",
    rxcui,
    "/property.json?propName=STRENGTH"
  )
  
  fromJSON(url)$propConceptGroup$propConcept |>
    select(propValue) |>
    #rename(strength = propValue) |>
    mutate(
      rxcui = rxcui,
      strength_val = as.double(word(propValue,1)),
      strength_unit = word(propValue,2),
    ) |>
    select(-propValue)
}

