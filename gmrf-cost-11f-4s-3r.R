#################################################################################
##### Benjamin M. Wagman      3/19/2019            
#####  11f, 3r, 4s, for 8 years of obs (Dec 2006- Nov 2014) 
#####  Imports: S_q, S_tot, S_rad, balance.txt, Q.Rdata. 
#####           ...without these, it wont run. 
#####
#####  FIELDS: 1) PRECIP RATE                             
#####          2) SEA LEVEL PRESSURE                               
#####          3) TREFHT               
#####          4) WIND SPEED @ 10m                             
#####          5) SWCF                             
#####          6) LWCF                             
#####          7) CALIPSO LOW FRACTION             
#####          8) CALIPSO MED FRACTION             
#####          9) CALIPSO HGH FRACTION             
#####          10) RH @ 300 hPa                           
#####          11) U  @ 300 hPa                   
#####          ...and Radiative balance scalar     
#####                                              
#####  Calls: process_obs from process_obs.R
#####         cost.nfields from docosts.R    
#####
#####  Outputs: AVG.COV.MAT aka "Sq" and obs.climate.Rdata (if these files don't already exist)
#####          results (dir) 
#####             -unscaled_cost (dir)
#####             -scaled_cost   (dir)
#####             -cost_for_MECS (dir)

#################################################################################

# Summary of script: 
        
  #   Get processed observational monthly data
  #     If the obs have not been processed (i.e. obs. climatologies and covariances have not been computed)
  #       Calculate seasonal obs. climatologies at gridpoints
  #       Calculate S, the spatial mean of the gridpoint variance and covariance, and invert it. (Only needs to be run once)
  #           Remove linear trend at each gridpoint using a linear model
  #           Calculate the variance of the detrended data (the residuals of the linear model) at each grid point
  #           Compute the spatial mean of the variances and covariances between variables at each gridpoint. 
  #           Save S as AVG.COV.MAT 
  #     else
  #       Load processed observational climatologies and S from Rdata
  #   Load Q, S, 
  #   For each region,
  #     and each season,
  #       Load model climatologies 
  #       Load obs climatologies
  #       Call cost.nfields with S=1, alpha=1 (Traditional cost) 
  #       Call cost.nfields with S containing off-diagonals, alpha=1 (Fields cost)
  #       Call cost.nfields with S containing off-diagonals, alpha < 1 (Field-and-space cost)

#################################################################################
# Load libraries, Q operator, weights, and constants.  
#################################################################################
library(ncdf4)  #Install these packages if needed. 
library(lattice)
library(akima)
library(fields)
library(R.matlab)

source('Rcode/process_obs.R') 
source('Rcode/docosts.R') # Does the actual model-data difference in function cost.nfields

setwd("./") 
        
#DECLARE NYEARS of obs data
nyears <- 8

# Load externally-generated operators and weights. 
print("loading Q operator...takes a little time", quote=FALSE)
load("external_weights/Q.Rdata") 
# Q[1] is the precision matrix for tropics
# Q[2] is the precision matrix for south
# Q[3] is the precision matrix for north (although S and N are interchangeable because same grid.)

load("external_weights/alpha_solutions.Rdata")
# optimal alpha value for GMRF for Tropics and extratropics.
# User should generate their own if using different grid size or regions. 
# See equation 5 in:
# Nosedal-Sanchez, A., C. Jackson, and G. Huerta, 2016. A new test statistic for climate models that includes field and spatial
# dependencies using Gaussian Markov random fields. Geoscientific Model Development, 9(7):2407–2414.
alpha_optimal<-c( tropic_alpha$root, poleward_alpha$root, poleward_alpha$root)

balance<-unlist( read.table("model_seasonal_climo/balance.txt")) # extract the annual mean radiative balance at TOA for this model. Scalar. 

fieldnames<-c("precip_obs" , "psl_obs" , "trefht_obs" , "speed_obs" , "swcf_obs" , "lwcf_obs", "cll_obs" , "clm_obs" , "clh_obs" , "rh300_obs" , "u300_obs")
nfields = length(fieldnames) 
seasonnames<- c("DJF", "MAM", "JJA","SON")

# Divide globe into 3 regions by latitude 
## Region: Tropics.
#Lats -30.63 to 30.63
tropics<-seq(64,129,by=1)  # lat index for tropics
## Region: Southern midlats
#Lats  -64.55 to -30.63
south<-seq(28,63,by=1)
## Region: Northern midlats
# Lats 30.63 to 64.55
north<-seq(130,165,by=1)
region<-list(tropics, south, north) # start with tropics, then south, then north. 
regionnames<-c("tropics" , "south_extratropics" , "north_extratropics")

#################################################################
# Process the obs, or skip it if already done. 
# If you've already processed the obs and want to do it again anyway, delete the file 'processed_obs/AVG.COV.MAT'

if (!file.exists('processed_obs/AVG.COV.MAT')){
  
  print("File processed_obs/AVG.COV.MAT not found; calling fn process_obs to produce it ", quote=FALSE)
  
  ## Monthly obs files
  # All time coordinates standardized to match ERA_interim
  # Obs span Dec 2006 through Nov 2014
  gpcp_prect_obs     <-nc_open("obs_monthly/PRECIP_RATE.nc") #
  era_psl_obs        <-nc_open("obs_monthly/MSL_ERA.nc")
  era_trefht_obs     <-nc_open("obs_monthly/T2M_ERA.nc")
  era_SI10_obs       <-nc_open("obs_monthly/SI10_ERA.nc")
  ceres_swcf_obs     <-nc_open("obs_monthly/SWCF.nc")
  ceres_lwcf_obs     <-nc_open("obs_monthly/LWCF.nc")
  calipso            <-nc_open("obs_monthly/calipso.nc")
  era_rh_obs         <-nc_open("obs_monthly/RH_300_ERA.nc")
  era_u300_obs       <-nc_open("obs_monthly/U_300_ERA.nc")
  
  # Obs data 
  precip_obs         <-ncvar_get(gpcp_prect_obs,"PRECIP_GPCP")      # obs field 1 PRECIP
  psl_obs            <-ncvar_get(era_psl_obs,"MSL_ERA")             # obs field 2 PSL
  trefht_obs         <-ncvar_get(era_trefht_obs,"T2M_ERA")          # obs field 3 TREFHT
  speed_obs          <-ncvar_get(era_SI10_obs,"SI10_ERA")           # obs field 4 10 m wind speed
  swcf_obs           <-ncvar_get(ceres_swcf_obs)                    # obs field 5 SWCF  
  lwcf_obs           <-ncvar_get(ceres_lwcf_obs)                    # obs field 6 LWCF
  cll_obs            <-ncvar_get(calipso,"CLL")                     # obs field 7 CALIPSO LOW CLOUD FRACTION
  clm_obs            <-ncvar_get(calipso,"CLM")                     # obs field 8 CALIPSO MED CLOUD FRACTION
  clh_obs            <-ncvar_get(calipso,"CLH")                     # obs field 9 CALIPSO HGH CLOUD FRACTION
  rh300_obs          <-ncvar_get(era_rh_obs)                        # obs field 10 RH @ 300 hPa
  u300_obs           <-ncvar_get(era_u300_obs)                      # obs field 11 U  @ 300 hPa 
    
  obs_list = list( precip_obs, psl_obs, trefht_obs, speed_obs, swcf_obs, lwcf_obs, cll_obs, clm_obs, clh_obs, rh300_obs, u300_obs) 
    
  process_obs( regionnames, list( tropics, south, north), fieldnames, obs_list ) # Fieldnames and obs_list must be ordered consistently. 
  
  print("finished processing observations. wrote processed_obs/AVG.COV.MAT and obs.climate.Rdata", quote=FALSE)
  
}  # close process_obs loop
  
#################################################################################
# Compare Obs to Model
#################################################################################
## Load Q, S, and model output. 

load('processed_obs/AVG.COV.MAT')   
S<-AVG.COV.MAT

load('processed_obs/obs.climate.Rdata')

# Load GMRF precision-matrix makier, Q.
# Q is a function of grid size only. To make your own 'Q', see:
# Nosedal-Sanchez, A., C. Jackson, and G. Huerta, 2016. A new test statistic for climate models that includes field and spatial
# dependencies using Gaussian Markov random fields. Geoscientific Model Development, 9(7):2407–2414.

# Load experiments
model<-("./") #local data
seas <- list("DJF", "MAM", "JJA", "SON") # Ordering or seasons consistenc with S

#WITHIN THIS LOOP...
# SELECT A REGION
  # SELECT A SEASON 
      # READ MODEL OUTPUT
      # READ CORRESPONDING OBS  
      # COMPUTE EACH TYPE OF COST

trad.cost.vec <- fields.cost.vec  <- fs.cost.vec  <- matrix(0,nrow=length(seas),ncol=3) # a row for each season and a column for each region (sum over fields)
trad.cost.matrix <- fields.cost.matrix <- fs.cost.matrix <-  array(0,dim=c(4, 3, nfields, nfields))  # "full cost matrix" (no summing over fields)

dimnames(trad.cost.matrix)[[1]] <-dimnames(fields.cost.matrix)[[1]]<- dimnames(fs.cost.matrix)[[1]] <- seasonnames
dimnames(trad.cost.matrix)[[2]] <- dimnames(fields.cost.matrix)[[2]] <- dimnames(fs.cost.matrix)[[2]] <-  c("TROP", "S", "N")
dimnames(trad.cost.matrix)[[3]] <- dimnames(fields.cost.matrix)[[3]]<- dimnames(fs.cost.matrix)[[3]]<- fieldnames
dimnames(trad.cost.matrix)[[4]]<-dimnames(fields.cost.matrix)[[4]]<-dimnames(fs.cost.matrix)[[4]]<-fieldnames

files=0 # path to model file names
files_p=0
for (i in 1:length(seas)){
  files[i]   <-Sys.glob( paste( "model_seasonal_climo/", seas[i], ".nc" , sep="")) # Character string of seasonal model output file
  files_p[i] <-Sys.glob( paste( "model_seasonal_climo/", seas[i], "-p.nc" , sep="")) # Character string of seasonal model output file on 300 hPa surface
}

print("Begin computing cost", quote=FALSE)
for (r in 1:length(region)){ #loop over regions
  print(paste ( "Region: ", regionnames[r], sep=""), quote=FALSE)  
  for (i in 1:length(seas)){ #loop over seasons 
    print(paste("     Season: ", seasonnames[i], sep=""), quote=FALSE)
    file<-files[i] # open the model output file
    file_p<-files_p[i]
    # Get model climo for region and season
    
    model_arg<-vector("list", length = length(fieldnames))                      
    obs_arg  <-vector("list", length = length(fieldnames))                     
    
    experiment<-nc_open(file)
    precc_exp<-ncvar_get(experiment,"PRECC")[ ,region[[r]]] 
    precl_exp<-ncvar_get(experiment,"PRECL")[ ,region[[r]]]
    model_arg[[1]]<-(precc_exp + precl_exp) * (1000*60*60*24)                          # Field 1
    model_arg[[2]]<-ncvar_get(experiment,"PSL")[ ,region[[r]]] # Pa.                   # Field 2 
    model_arg[[3]]<-ncvar_get(experiment,"TREFHT")[ ,region[[r]]]                      # Field 3
    model_arg[[4]]<-ncvar_get(experiment,"U10")[ ,region[[r]]] #Wind speed (not zonal) # Field 4
    model_arg[[5]]<-ncvar_get(experiment,"SWCF")[ ,region[[r]]]                        # Field 5                       
    model_arg[[6]]<-ncvar_get(experiment,"LWCF")[ ,region[[r]]]                        # Field 6
    model_arg[[7]]<-ncvar_get(experiment,"CLDLOW_CAL")[ ,region[[r]]] * 0.01           # Field 7
    model_arg[[8]]<-ncvar_get(experiment,"CLDMED_CAL")[ ,region[[r]]] * 0.01           # Field 8
    model_arg[[9]]<-ncvar_get(experiment,"CLDHGH_CAL")[ ,region[[r]]] * 0.01           # Field 9
    experiment <- nc_open(file_p) # Pressure_level data
    model_arg[[10]]<-ncvar_get(experiment,"RELHUMonP")[ ,region[[r]]]                  # Field 10
    model_arg[[11]]<-ncvar_get(experiment,"UonP")[ ,region[[r]]]                       # Field 11

    # Get obs climo for region and season. 
    obs_file<-paste(seas[i],"_list", sep="") 
    obs_data<-get(obs_file)  
    obs<-obs_data[[r]] 
    obs_arg = list(obs[ , , 1], obs[ , , 2], obs[ , , 3], obs[ , , 4], obs[ , , 5] , obs[ , , 6], 
                   obs[ , , 7], obs[ , , 8], obs[ , , 9], obs[ , , 10],obs[ , , 11]  )
    
    rm(obs) # will change size so must be removed 

    #COMPUTE COST FOR CURRENT SEASON AND REGION

    ## TRADITIONAL COST: S = diag(S) and alpha =1.
    #Set alpha, and format S.
    alpha<-1
    trad.S<-S[ , , i,r]                                                    #initialize traditional S as the same as S
    trad.S[]<-0                                                            #then zero-it out. 
    diag(trad.S)<-diag(S[ , , i,r])                                        #The traditional S only includes variance 
    inv.trad.S<-solve(trad.S)
    # For each seson and region, cost.matrix is nfields x nfields. 
    trad.cost.matrix[i,r, , ]<-cost.nfields(obs_arg,model_arg,Q[[r]],inv.trad.S,alpha)
    trad.cost.vec[i,r]<-sum(trad.cost.matrix[i,r, , ]) # scalar version of cost 
    print("          Traditional cost = ...", quote=FALSE)
    print( format(round(trad.cost.matrix[i,r, , ], 2), nsmall = 2), quote=FALSE)
    
    
    ## FIELDS COST: S = full(S) and alpha =1.
    #set alpha, and format S.
    alpha<-1
    fields.S<-S[ , ,i,r]
    inv.fields.S<-solve(fields.S)
    fields.cost.matrix[i,r, , ]<-cost.nfields(obs_arg,model_arg,Q[[r]],inv.fields.S,alpha)
    fields.cost.vec[i,r]<-sum(fields.cost.matrix[i,r, , ])
    print("          Fields cost = ...", quote=FALSE)
    print( format(round(fields.cost.matrix[i,r, , ], 2), nsmall = 2), quote=FALSE)
    
    ## FIELDS-SPACE COST: S = full(S) and alpha = optimal alpha for regional grid. 
    ## set alpha, and format S.
    alpha <- alpha_optimal[r]
    fs.S<-S[ ,  ,i,r]
    inv.fs.S<-solve(fs.S)
    fs.cost.matrix[i,r, , ]<-cost.nfields(obs_arg,model_arg,Q[[r]],inv.fs.S,alpha)
    fs.cost.vec[i,r]<-sum(fs.cost.matrix[i,r, , ])
    print("          fields and space cost = ...", quote=FALSE)
    print( format(round(fs.cost.matrix[i,r, , ], 2), nsmall = 2), quote=FALSE)
    
  } #close the season loop. 

  
} #close the region loop. 

trad.cost.scalar=sum(trad.cost.vec)
fields.cost.scalar=sum(fields.cost.vec)
fs.cost.scalar=sum(fs.cost.vec)
#################################################################################
# Compute Radiative balance and begin to write results. 
#################################################################################
# Compute the annual mean radiative balance at TOA from seasonal climatology files of model output. 
# Targeting + 0.75 W/m^2, meaning FSNT - FLNT = 0.75 W/m^2
# Edited to target 0.9 W/m^2 based on ocean heat uptake estimates:
# Trenberth et al., (2016). Insights into Earth’s energy imbalance from multiple sources

dir.create('results',showWarnings = FALSE) # cost output before multipyling by S_tot and Sq

dir.create('results/unscaled_cost',showWarnings = FALSE) # cost output before multipyling by S_tot and Sq
dir.create('results/scaled_cost',showWarnings = FALSE)   # cost output after  multipyling by S_tot and Sq
dir.create('results/cost_for_MECS',showWarnings = FALSE) # cost output after  multipyling by Sq but NOT multiplying by S_tot. 

rad_cost_unsc<- (balance - 0.90)^2/(2 * 0.5^2) # Target = + 0.75 W/m^2. Standard deviation = 1/2 W/m^2, but this is arbitrary and gets renormalized by S_rad.  
write(rad_cost_unsc,file="results/unscaled_cost/rad_cost.dat_unsc",append = FALSE)

rownames(trad.cost.vec)<-rownames(fields.cost.vec)<-rownames(fs.cost.vec)<-seasonnames
colnames(trad.cost.vec)<-colnames(fields.cost.vec)<-colnames(fs.cost.vec)<-c("TROP" , "S" , "N" )

#################################################################################
# Save 'unscaled' costs in all forms, before scaling by Sq and Stot
#################################################################################

#write scalar cost
write(trad.cost.scalar,file="results/unscaled_cost/result_trad.dat_unsc")
write(fields.cost.scalar,file="results/unscaled_cost/result_fields.dat_unsc")
write(fs.cost.scalar,file="results/unscaled_cost/result_fs_scalar.dat_unsc")

# save season by region costs
save(trad.cost.vec,fields.cost.vec,fs.cost.vec,file="results/unscaled_cost/cost_vecs_unsc.Rdata")

#Save text versions.  4x3 (DJF,MAM,JJA,SON) x (TROP,S,N)
write.table(trad.cost.vec, file="results/unscaled_cost/trad.cost.vec_unsc.txt", row.names=FALSE, col.names=FALSE)
write.table(fields.cost.vec, file="results/unscaled_cost/fields.cost.vec_unsc.txt", row.names=FALSE, col.names=FALSE)
write.table(fs.cost.vec, file="results/unscaled_cost/fs.cost.vec_unsc.txt", row.names=FALSE, col.names=FALSE)

save(trad.cost.matrix,fields.cost.matrix,fs.cost.matrix,file="results/unscaled_cost/cost_mats_unsc.Rdata")

#################################################################################
# Scale costs by Sq and Stot and save
#################################################################################
# Import the S_q and S_tot, calculated on matlab dof.m using output from the perfect model experiments. 
S_weights = readMat('external_weights/S_and_q_coeffs.mat')

# S_weights$S_q is 3x3x4 (cost type x region x season)
# Immediately reshape S_q such that it is season x region to match cost vecs. 
S_q_tmp     <- S_weights$S.q[,,] 
S_q         <-array(0,dim=c(3,4,3))
S_q[1, , ]  <- t(S_q_tmp[1,,])
S_q[2, , ]  <- t(S_q_tmp[2,,])
S_q[3, , ]  <- t(S_q_tmp[3,,])

S_tot <- S_weights$S.tot
S_rad <- S_weights$S.rad

# Scaled cost vector  = S_tot * Sq(season x region) * cost(season x region) 
trad.cost.vec.s   = S_tot[1] *  S_q[1,,] * trad.cost.vec 
fields.cost.vec.s = S_tot[2] *  S_q[2,,] * fields.cost.vec 
fs.cost.vec.s     = S_tot[3] *  S_q[3,,] * fs.cost.vec 

#Save text versions.  4x3 (DJF,MAM,JJA,SON) x (TROP,S,N)
write.table(trad.cost.vec.s, file="results/scaled_cost/trad.cost.vec.s.txt", row.names=FALSE, col.names=FALSE)
write.table(fields.cost.vec.s, file="results/scaled_cost/fields.cost.vec.s.txt", row.names=FALSE, col.names=FALSE)
write.table(fs.cost.vec.s, file="results/scaled_cost/fs.cost.vec.s.txt", row.names=FALSE, col.names=FALSE)

# Scaled cost matrix = S_tot * Sq (season x region) * cost (season x region x field x field)
trad.cost.matrix.s  <- fields.cost.matrix.s<- fs.cost.matrix.s    <- array(0,dim=c(length(seas),length(region),nfields,nfields))
dimnames(trad.cost.matrix.s)[[1]] <- dimnames(fields.cost.matrix.s)[[1]] <- dimnames(fs.cost.matrix.s)[[1]] <- seasonnames
dimnames(trad.cost.matrix.s)[[2]] <- dimnames(fields.cost.matrix.s)[[2]] <- dimnames(fs.cost.matrix.s)[[2]] <- c("TROP", "S", "N")
dimnames(trad.cost.matrix.s)[[3]] <- dimnames(fields.cost.matrix.s)[[3]] <- dimnames(fs.cost.matrix.s)[[3]] <- fieldnames
dimnames(trad.cost.matrix.s)[[4]] <- dimnames(fields.cost.matrix.s)[[4]] <- dimnames(fs.cost.matrix.s)[[4]] <- fieldnames

for (season in 1:4){
  for (region in 1:3){
    trad.cost.matrix.s[season,region, , ]   = drop(S_q[1,season,region]) * trad.cost.matrix[season, region, , ]
    fields.cost.matrix.s[season,region, , ] = drop(S_q[2,season,region]) * fields.cost.matrix[season, region, , ]
    fs.cost.matrix.s[season,region, , ]     = drop(S_q[3,season,region]) * fs.cost.matrix[season, region, , ]    
  }  
}
save(trad.cost.matrix.s,fields.cost.matrix.s,fs.cost.matrix.s,file="results/scaled_cost/cost_mats.Rdata")

# The ensemble control system generates its own S_tot; it does not generate an S_q. 
# Therefore, for the ensemble control system, write cost scaled by Sq but not S_tot. 
trad.cost.scalar.s    =  sum( S_q[1,,] * trad.cost.vec )    + (rad_cost_unsc * S_rad)
fields.cost.scalar.s  =  sum( S_q[2,,] * fields.cost.vec )  + (rad_cost_unsc * S_rad)
fs.cost.scalar.s      =  sum( S_q[3,,] * fs.cost.vec )      + (rad_cost_unsc * S_rad)

write(trad.cost.scalar.s,file="results/cost_for_MECS/result_trad.dat")
write(fields.cost.scalar.s,file="results/cost_for_MECS/result_fields.dat")
write(fs.cost.scalar.s,file="results/cost_for_MECS/result_fs_scalar.dat")
write(fs.cost.scalar.s,file="results/cost_for_MECS/result.dat")
#################################################################################
print("finished. Check results dir", quote=FALSE)
# End of cost. 
#################################################################################