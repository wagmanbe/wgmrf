###########################################################
### FUNCTION FOR PROCESSING OBS
### Benjamin M. Wagman

# Inputs: Arg1: A list of regions and their lat indicies
# Inputs: Arg2: A list of obs data by field
# Returns: nothing
# Saves: AVG.COV.MAT, AVG.COR.MAT, obs.climate.Rdata
###########################################################

process_obs<-function(regionnames, region, fieldnames, fields_list){

  nyears=8  
  nlon = 288
  nfields = length(fieldnames)

  print("Generating obs climos and covariances. This takes a few minutes", quote=FALSE)
  print("...but it only needs to be done once.", quote=FALSE)
  print("Processing the obs produces obs.climate.Rdata and AVG.COV.MAT", quote=FALSE)

  ## Step 1. Get index for each month.
  # Overlapping time frame severely limited by Calipso. Obs netcdf files start in December. 
  dec<-seq(1, 1 + (nyears-1) * 12, by=12)   
  jan<-dec+1
  feb<-jan+1
  mar<-jan+2
  apr<-jan+3
  may<-jan+4
  jun<-jan+5
  jul<-jan+6
  aug<-jan+7
  sep<-jan+8
  oct<-jan+9
  nov<-jan+10
  
  # Initialization.  
  s_djf<- s_mam<- s_jja<- s_son <- array(0,dim=c(nfields,nfields,length(region)))
    
  #initiate the climate lists (lon x lat x nfields), for each region.  
  DJF_list<- MAM_list<- JJA_list<- SON_list<- list(array(0,dim=c(nlon,length(region[1]),nfields)),
                                                   array(0,dim=c(nlon,length(region[2]),nfields)),
                                                   array(0,dim=c(nlon,length(region[3]),nfields)))
    
  AVG.COV.MAT<- S<-  AVG.COR.MAT<- array(0,dim=c(nfields,nfields,4,length(region))) #one for each season and region
 
  for (r in 1:length(region)){ #For each region...
    print(regionnames[r], quote=FALSE)
    
    #################################################################################
    # Calculate S, the spatial mean of the gridpoint variance and covariance, and invert it. 
    #################################################################################
    
    #Calculate a seasonal mean for EACH year at each grid point (big_fields_xxx)
    print( paste("     calculating seasonal mean for each year at each gridpoint for region",regionnames[r]), quote=FALSE)
    big_fields_djf<- big_fields_mam<- big_fields_jja<- big_fields_son<- array(0,dim=c(nlon*length(region[[r]]), nyears, nfields)) #array of npoints, nyears, nfields
        
    #Eval works by operating on the first ROW, then the second, and preserving the order. 
    for (l in 1:nfields){ #loop over all fields 
      field <- as.name(fieldnames[l]) # select the field
      #print(fieldnames[l])
      year<-1
      indice<-seq(1,nyears,by=1)
      
      for (i in indice ){
        avg.dec<-matrix(eval(field)[ , region[[r]], dec[i]])  #select the winter data from that field for a particular year. 
        avg.jan<-matrix(eval(field)[ , region[[r]], jan[i]])   
        avg.feb<-matrix(eval(field)[ , region[[r]], feb[i]])  
        
        avg.mar<-matrix(eval(field)[ , region[[r]], mar[i]])  #select the spring data from that field
        avg.apr<-matrix(eval(field)[ , region[[r]], apr[i]])   
        avg.may<-matrix(eval(field)[ , region[[r]], may[i]])  
        
        avg.jun<-matrix(eval(field)[ , region[[r]], jun[i]])  #select the summer data from that field
        avg.jul<-matrix(eval(field)[ , region[[r]], jul[i]])   
        avg.aug<-matrix(eval(field)[ , region[[r]], aug[i]])  
        
        avg.sep<-matrix(eval(field)[ , region[[r]], sep[i]])  #select the fall data from that field
        avg.oct<-matrix(eval(field)[ , region[[r]], oct[i]])   
        avg.nov<-matrix(eval(field)[ , region[[r]], nov[i]])  
        
        big_fields_djf[ , year, l ]<-(avg.dec + avg.jan + avg.feb)/3    # 1-year mean winter at each lat and lon
        big_fields_mam[ , year, l ]<-(avg.mar + avg.apr + avg.may)/3    # 1-year mean 
        big_fields_jja[ , year, l ]<-(avg.jun + avg.jul + avg.aug)/3    # 1-year mean 
        big_fields_son[ , year, l ]<-(avg.sep + avg.oct + avg.nov)/3    # 1-year mean 
        
        if (year < nyears) {year<- year + 1}
      }
    }# close fields loop. 
    
    #REBUILD THE LAT X LON GRID FROM THE VECTOR IN BIG_FIELDS TO MAKE THE CLIMATOLOGY at each point
    DJF_climate<-array(0,dim=c(nlon, length(region[[r]]), nfields))
    MAM_climate<-array(0,dim=c(nlon, length(region[[r]]), nfields))
    JJA_climate<-array(0,dim=c(nlon, length(region[[r]]), nfields))
    SON_climate<-array(0,dim=c(nlon, length(region[[r]]), nfields))
    
    j<-seq(from=1, to=length(avg.dec), by = nlon) # 1: 19008 by 288
    
    for (l in 1:nfields){ #loop over all fields 
      for (i in 1:length(j)){           # i goes from 1:length(region[[r]]) and is the lat index. use to select all 288 lons for that lat. 
        
        field <- as.name(fieldnames[l]) # select the field
        
        time_mean_djf<-rowMeans(big_fields_djf[ , ,l]) # average that field over all years. 
        time_mean_mam<-rowMeans(big_fields_mam[ , ,l]) # 
        time_mean_jja<-rowMeans(big_fields_jja[ , ,l]) # 
        time_mean_son<-rowMeans(big_fields_son[ , ,l]) # 
        
        index<-seq(from= j[i], to= j[i]+287,by=1) #select the lats for a particular lon. 
        
        lats_djf<-c(time_mean_djf[index])
        lats_mam<-c(time_mean_mam[index])
        lats_jja<-c(time_mean_jja[index])
        lats_son<-c(time_mean_son[index])
        
        DJF_climate[ ,i , l]<-lats_djf # put data from vector back on Lat x Lon grid
        MAM_climate[ ,i , l]<-lats_mam #
        JJA_climate[ ,i , l]<-lats_jja #
        SON_climate[ ,i , l]<-lats_son #
        
        
      }
    } #close fields loop
    
    DJF_list[[r]]<-DJF_climate[ , , ] # Climatologies (Region x Lat x Lon x Variable)
    MAM_list[[r]]<-MAM_climate[ , , ]
    JJA_list[[r]]<-JJA_climate[ , , ]
    SON_list[[r]]<-SON_climate[ , , ]
    dir.create('processed_obs',showWarnings = FALSE) 
    save(DJF_list,MAM_list,JJA_list,SON_list,file="processed_obs/obs.climate.Rdata")
    print("     saved Lat x Lon x Variable climatology for each season in obs.climate.Rdata", quote=FALSE)
    
    # Remove linear trend at each gridpoint 
    #print( "     removing linear trend from climatology at each gridpoint")
    
    resid1_djf<-matrix(0, nyears, nlon*length(region[[r]]))   # initializing matrix of residuals 
    resid1_mam<-matrix(0, nyears, nlon*length(region[[r]]))   # initializing matrix of residuals 
    resid1_jja<-matrix(0, nyears, nlon*length(region[[r]]))   # initializing matrix of residuals 
    resid1_son<-matrix(0, nyears, nlon*length(region[[r]]))   # initializing matrix of residuals 
    
    resid2_djf<-matrix(0, nyears, nlon*length(region[[r]]) )   
    resid2_mam<-matrix(0, nyears, nlon*length(region[[r]]) )   
    resid2_jja<-matrix(0, nyears, nlon*length(region[[r]]) )   
    resid2_son<-matrix(0, nyears, nlon*length(region[[r]]) )   
    
    cv_djf<-matrix(0,nlon*length(region[[r]])) #initializing the vector of covariances at each point 
    cv_mam<-matrix(0,nlon*length(region[[r]])) #initializing the vector of covariances at each point 
    cv_jja<-matrix(0,nlon*length(region[[r]])) #initializing the vector of covariances at each point 
    cv_son<-matrix(0,nlon*length(region[[r]])) #initializing the vector of covariances at each point  
    
    x<-seq(1,nyears,by=1)  
    
    #calculate the variance and covariance between fields by using the residals of the linear models at each gridpoint
    print("     at each gridpoint, compute variance in each field for each season over nyears", quote=FALSE)
    print("     and similarly the covariance with other fields at the same gridpoint (expensive)", quote=FALSE)
    
    for (field1 in 1:nfields){    #loop over fields
      #print(fieldnames[field1])
      for (field2 in 1:nfields){
        if (field1 >= field2 ){    #Only calculate half of the covariance matrix--later it will be reflected over diag
          
          model1_djf<-lm(t(big_fields_djf[ , ,field1])~x) 
          model1_mam<-lm(t(big_fields_mam[ , ,field1])~x) 
          model1_jja<-lm(t(big_fields_jja[ , ,field1])~x) 
          model1_son<-lm(t(big_fields_son[ , ,field1])~x) 
          
          resid1_djf[,]<-model1_djf$res
          resid1_mam[,]<-model1_mam$res
          resid1_jja[,]<-model1_jja$res
          resid1_son[,]<-model1_son$res
          
          if (field1 == field2){
            cv_djf[,1]<-apply (resid1_djf[,], 2 , var)   # Variance of the vectors of yearly residuals at the same gridpoint
            cv_mam[,1]<-apply (resid1_mam[,], 2 , var) 
            cv_jja[,1]<-apply (resid1_jja[,], 2 , var) 
            cv_son[,1]<-apply (resid1_son[,], 2 , var) 
          }
          else      {
            model2_djf<-lm(t(big_fields_djf[ , , field2])~x) 
            model2_mam<-lm(t(big_fields_mam[ , , field2])~x)
            model2_jja<-lm(t(big_fields_jja[ , , field2])~x)
            model2_son<-lm(t(big_fields_son[ , , field2])~x)
            
            resid2_djf[,]<- model2_djf$res
            resid2_mam[,]<- model2_mam$res
            resid2_jja[,]<- model2_jja$res
            resid2_son[,]<- model2_son$res
            
            for (i in 1:(nlon*length(region[[r]]))){
              
              cv_djf[i,1]<-var(resid1_djf[,i],resid2_djf[,i]) # Covariance of the vectors of yearly residuals at the same gridpoint. 
              cv_mam[i,1]<-var(resid1_mam[,i],resid2_mam[,i])
              cv_jja[i,1]<-var(resid1_jja[,i],resid2_jja[,i])
              cv_son[i,1]<-var(resid1_son[,i],resid2_son[,i])
              
            }
          }
          # Step 2.5 
          # Compute the spatial mean of the each gridpoint's variance and field-covariance 
          s_djf[field1,field2,r]<-apply(cv_djf, 2, FUN=mean) #the spatial mean of the variance or the covariance of the gridpoint residuals. 
          s_mam[field1,field2,r]<-apply(cv_mam, 2, FUN=mean) 
          s_jja[field1,field2,r]<-apply(cv_jja, 2, FUN=mean) 
          s_son[field1,field2,r]<-apply(cv_son, 2, FUN=mean) 
          
        }   
      } 
    } #close fields loop. 
    
    
    #s
    s_djf[ , ,r]<-(t(s_djf[ , ,r])+s_djf[ , ,r])#flip it over the diagonal. 
    s_mam[ , ,r]<-(t(s_mam[ , ,r])+s_mam[ , ,r])
    s_jja[ , ,r]<-(t(s_jja[ , ,r])+s_jja[ , ,r])
    s_son[ , ,r]<-(t(s_son[ , ,r])+s_son[ , ,r])
    
    diag(s_djf[ , ,r])<-diag(s_djf[ , ,r])/2
    diag(s_mam[ , ,r])<-diag(s_mam[ , ,r])/2
    diag(s_jja[ , ,r])<-diag(s_jja[ , ,r])/2
    diag(s_son[ , ,r])<-diag(s_son[ , ,r])/2
    
    s_djf[ , ,r]  
    s_mam[ , ,r]
    s_jja[ , ,r]
    s_son[ , ,r]
    
    
  }# Close the regions loop 
  
  ## Average covariance matrix (formatted). Dim is nfields x nfields x season x region. 
  
  for (r in 1:3){
    
    AVG.COV.MAT[ , , 1, r]<-s_djf[ , ,r]
    AVG.COV.MAT[ , , 2, r]<-s_mam[ , ,r]
    AVG.COV.MAT[ , , 3, r]<-s_jja[ , ,r]
    AVG.COV.MAT[ , , 4, r]<-s_son[ , ,r]
    
    AVG.COR.MAT[ , , 1, r]<-cov2cor(AVG.COV.MAT[ , , 1, r]) # Compute the correlation matrix to see if it makes sense. It will not be used for any math. 
    AVG.COR.MAT[ , , 2, r]<-cov2cor(AVG.COV.MAT[ , , 2, r])
    AVG.COR.MAT[ , , 3, r]<-cov2cor(AVG.COV.MAT[ , , 3, r])
    AVG.COR.MAT[ , , 4, r]<-cov2cor(AVG.COV.MAT[ , , 4, r])
  }
  
  ## Average correlation matrix. 
  for (season in 1:4){
    for (r in 1:3){
      AVG.COR.MAT[ , , season, r]<-cov2cor(AVG.COV.MAT[ , , season, r]) #Dim is nfields x nfields x season x region. 
    }
  }
  
  ## Finally, S = Average covariance matrix
  rownames(AVG.COV.MAT)<-colnames(AVG.COV.MAT)<- rownames(AVG.COR.MAT)<- colnames(AVG.COR.MAT)<- 
    c("precip_obs" , "psl_obs" , "trefht_obs" , "speed_obs" , "swcf_obs" , "lwcf_obs" , "cll_obs" , "clm_obs" , "clh_obs" , "rh300_obs" , "u300_obs") 
    
  # "S" is AVG.COV.MAT
  save(AVG.COV.MAT, file = "processed_obs/AVG.COV.MAT")
  save(AVG.COR.MAT, file = "processed_obs/AVG.COR.MAT")
  

}
