# --------------------- OBJECTIVE: generate synthetic data for code publication -------------------------------------


library(truncnorm)

generateData <- function(){
  
  writeLines("Generating synthetic data")
  
  ids <- 1:3710
  
  # generates a vector with random integers that represent the number of times a participant has follow-up
  reps <- round(rtruncnorm(length(ids), a = 0.5, b = Inf, mean = 18, sd = 12)) 
  reps_zphi <- round(rtruncnorm(length(ids), a = -Inf, b = Inf, mean = 1.9, sd = 0.9)) # fewer fups in ZPHI
  
  
  df <- NULL
  df_zphi <- NULL
  for (i in ids){                       # create data for each participant separately
    
    if (i %% 100 == 0){cat("-")}
    
    # tried to sample from distributions with roughly correct means, but most likely made waay too many 
    # assumptions about normal distributions ...
    
    n <- reps[i]                        # number of follow-ups
    id <- rep(as.character(i), n)       # id
    risk <- rtruncnorm(1, a = 0, b = 2, mean = 1, sd = 0.5) # multiply with nsP/nsCAI/STI/syphilis odds to add some covariance
    age <- rep(round(rtruncnorm(1, a = 16, b = 100, mean = 47, sd = 13)), n)
    fupdate <- as.character(sample(seq(as.Date('2000-01-01'), as.Date('2020-06-30'), by="day"), n))
    nsp <- rbinom(n, 1, 0.24*risk)           # sampled nsP trajectory
    nscai <- rbinom(n, 1, 0.14*risk)         # sampled nsCAI trajectory
    syph <- rbinom(n, 1, 0.02*risk)          # sampled syphilis events
    std <- rbinom(n, 1, 0.05*risk)           # sampled STI events
    
    df.i <- cbind(id, age, fupdate, nsp, nscai, syph, std)
    df <- rbind(df, df.i)
    
    # make zphi block
    n_zphi <- reps_zphi[i] 
    
    if (n_zphi > 0){           # number of fups allowed to be zero because ZPHI is subcohort (in our analyses)
      
      id_zphi <- rep(as.character(i), n_zphi)
      screening_date <- as.character(sample(seq(as.Date('2015-06-01'), as.Date('2016-06-16'), by="day"), n_zphi))
      number_sexpartners <- round(rtruncnorm(n_zphi, a = 0, b = 100, mean = 4*risk, sd = 8))
      
      df.i_zphi <- cbind(id_zphi, screening_date, number_sexpartners)
      df_zphi <- rbind(df_zphi, df.i_zphi)
    
      }
  }
  
  
  df <- data.frame(df)
  df$fupdate <- as.Date(df$fupdate)
  df$timesfup <- 1
  
  colnames(df_zphi)[colnames(df_zphi) == "id_zphi"] <- "id"
  df_zphi <- data.frame(df_zphi)
  df_zphi$screening_date <- as.Date(df_zphi$screening_date)
  
  
  write_csv(df, "datagen"); writeLines(paste("\nCSV created: datagen"))
  write_csv(df_zphi, "datagen_zphi"); writeLines(paste("CSV created: datagen_zphi"))
  
}



dataPrepSHCS <- function(studystart, studyend, endend){
  
  writeLines("Preparing SHCS dataset (data_shcs_ ...)")
  
  studystart <- as.Date(studystart,'%Y%m%d')
  studyend <- as.Date(studyend, '%Y%m%d')
  endend <- as.Date(endend, '%Y%m%d')
  
  if (file.exists("datagen") == F){generateData()}
  df <- read.csv("datagen")
  
  
  # ---- SEX ---------------------------------------------------------------------------------------------------------
  
  df$fupdate <- as.Date(df$fupdate)
  
  # SEX: get last answer on nsp and nsCAI before studyend
  nsp_before <- findLastValue(df, 'nsp', 'fupdate', studystart, studyend); df <- merge(df, nsp_before, by = 'id', all = T)
  nscai_before <- findLastValue(df, 'nscai', 'fupdate', studystart, studyend); df <- merge(df, nscai_before, by = 'id', all = T)
  
  # SEX: get first answer on nsp and nsCAI after studyend
  if (studyend != endend){
    nsp_after <- findFirstValue(df, 'nsp', 'fupdate', as.Date(studyend)+1, as.Date(endend)); df <- merge(df, nsp_after, by = 'id', all = T)
    nscai_after <- findFirstValue(df, 'nscai', 'fupdate', as.Date(studyend)+1, as.Date(endend)); df <- merge(df, nscai_after, by = 'id', all = T) 
  }
  
  # SEX: get 'ever' for nsp and nsCAI between studystart and studyend
  ever_nsp <- getEver(df, 'nsp', 'fupdate', studystart, studyend); df <- merge(df, ever_nsp, by = 'id', all = T)
  ever_nscai <- getEver(df, 'nscai', 'fupdate', studystart, studyend); df <- merge(df, ever_nscai, by = 'id', all = T)
  
  
  
  # ---- STDs -------------------------------------------------------------------------------------------------------------
  
  if (studyend <= as.Date('2017-05-01')){
    
    std_after <- getEver(df, 'std', 'fupdate', as.Date(studyend)+1, as.Date(endend))
    df <- merge(df, std_after, by = 'id', all = T); colnames(df)[colnames(df) == 'ever_std'] <- 'std_after'
    num_std_after <- getNumEvents(df, 'std', 'fupdate', as.Date(studyend)+1, as.Date(endend))
    df <- merge(df, num_std_after, by = 'id', all = T); colnames(df)[colnames(df) == 'num_std'] <- 'num_std_after'
    
  }
  
  # ---- SYPH -------------------------------------------------------------------------------------------------------------
  
  
  # SYPH: get ever syph event
  syph_before <- getEver(df, 'syph', 'fupdate', studystart, studyend)
  df <- merge(df, syph_before, by = 'id', all = T); colnames(df)[colnames(df) == 'ever_syph'] <- 'syph_before'
  
  # SYPH: get number of syph events before cutoff
  num_syph_before <- getNumEvents(df, 'syph', 'fupdate', studystart, studyend)
  df <- merge(df, num_syph_before, by = 'id', all = T); colnames(df)[colnames(df) == 'num_syph'] <- 'num_syph_before'
  
  # SYPH: get syph and number of syph after cutoff
  if (studyend != endend){
    
    syph_after <- getEver(df, 'syph', 'fupdate', as.Date(studyend)+1, as.Date(endend))
    df <- merge(df, syph_after, by = 'id', all = T); colnames(df)[colnames(df) == 'ever_syph'] <- 'syph_after'
    num_syph_after <- getNumEvents(df, 'syph', 'fupdate', as.Date(studyend)+1, as.Date(endend))
    df <- merge(df, num_syph_after, by = 'id', all = T); colnames(df)[colnames(df) == 'num_syph'] <- 'num_syph_after'
    
  }
  
  df <- unique(select(df, -c(fupdate, nscai, nsp, syph, std, timesfup)))
  name_csv <- paste('data_shcs', format.Date(studystart, '%Y%m%d'), format.Date(studyend, '%Y%m%d'), format.Date(endend, '%Y%m%d'), sep ='_')
  write_csv(df, name_csv); writeLines(paste("CSV created:", name_csv))
  
}


dataPrepZPHI <- function(studystart, studyend, endend){
  
  writeLines("Preparing ZPHI dataset (data_zphi_ ...)")
  
  startdate <- studystart
  enddate <-  studyend
  
  if (file.exists("datagen") == F){generateData()}
  sex <- read.csv("datagen")
  zphi <- read.csv("datagen_zphi")
  
  sex$fupdate <- as.Date(sex$fupdate)
  sex <- sex[sex$fupdate <= as.Date(enddate, '%Y%m%d'),] 
  
  num.entries <- NULL
  c <- NULL
  for (j in unique(zphi$id)){
    a <- merge(zphi[which(zphi$id == j), c('id', 'screening_date')], sex[c('id', 'age', 'fupdate', 'nscai', 'nsp')], by='id', all = F)
    a <- a[!is.na(a$screening_date),]
    a$screening_date <- as.character(a$screening_date)
    num.entries <- c(num.entries, nrow(zphi[which(zphi$id == j), c('id', 'screening_date')]))
    
    #print(nrow(zphi[which(zphi$id == j), c('id', 'screening_date')]))
    if (nrow(a) != 0){
      b <- NULL
      for (k in unique(a$screening_date)){
        b <- rbind(b, a[tapply(1:nrow(a), a$id, function(ii) ii[which.min(abs(as.Date(a$fupdate) - as.Date(k)))]),])
      }
      b$screening_date <- unique(a$screening_date)
      b$mean_number_sexpartners <- rep(mean(zphi[zphi$id == j, ]$number_sexpartners), nrow(b))
      b$mean_number_sexpartners <- rep(mean(zphi[zphi$id == j, ]$number_sexpartners), nrow(b))
      c <- rbind(c, b)
    }
    if (j %% 100 == 0){cat("-")}
  }
  
  # summary of number of entries per participant
  #num.entries <- data.frame(num.entries); summary(num.entries$num.entries)
  
  c$screening_date <- as.Date(c$screening_date)
  zphi$screening_date <- as.Date(zphi$screening_date)
  db <- merge(c, zphi, by = c('id', 'screening_date'))
  
  
  # write CSV
  name_csv <- paste('data_zphi', studystart, studyend, endend, sep ='_')
  write_csv(db, name_csv); writeLines(paste("\nCSV created:", name_csv))
  
}


