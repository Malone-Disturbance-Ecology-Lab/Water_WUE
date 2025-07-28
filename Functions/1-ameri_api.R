
library("amerifluxr")
#######################################################################################################################################

# download ameriflux data directly from ameriflux site

amf_download_base(user_id = "ammara_786",  # add your user id
                  user_email = "talib@wisc.edu",  # add your email
                  site_id = c("US-EvM","US-KS3", "US-KS4", "US-Skr", "US-TaS", 
                              "US-Dmg", "US-EDN","US-EKH","US-EKP","US-Myb","US-Srr","US-Tw1",    
                              "US-StJ",
                              "US-LA1","US-LA2","US-LA3",
                              "US-PHM",
                              "US-HB1","US-HB2","US-HB3",
                              "US-A03","US-A10","US-Atq","US-NGB"
                  ),
                  data_product = "BASE-BADM",
                  data_policy = "CCBY4.0",
                  agree_policy = TRUE,
                  intended_use = "model",
                  intended_use_text = "salinity effect on WUE, mutisynthesis",
                  verbose = TRUE,
                  out_dir = "\\\\corellia.environment.yale.edu\\MaloneLab\\Research\\WUE_CUE\\ameri_data")

# Call the function
amf_download_base("ammara_786", "talib@wisc.edu", c("US-KS3", "US-KS4"))


################################################################################



#download ameriflux data directly in your folder
site <- amf_site_info()

# if code does not run for all sites, you can select fewer sites


amf_download_base <- function(user_id, user_email, site_id) {
  # Ensure you are calling the correct function for downloading
  amerifluxr::amf_download_base(
    user_id = user_id,
    user_email = user_email,
    site_id = site_id,
    data_product = "BASE-BADM",
    data_policy = "CCBY4.0",
    agree_policy = TRUE,
    intended_use = "model",
    intended_use_text = "salinity effect on WUE, mutisynthesis",
    verbose = TRUE,
    out_dir = "\\\\corellia.environment.yale.edu\\MaloneLab\\Research\\WUE_CUE\\ameri_data"
  )
}






############################################################################################

site_id = c("US-EvM","US-KS3", "US-KS4", "US-Skr", "US-TaS", 
            "US-Dmg", "US-EDN","US-EKH","US-EKP","US-Myb","US-Srr","US-Tw1",    
            "US-StJ",
            "US-LA1","US-LA2","US-LA3",
            "US-PHM",
            "US-HB1","US-HB2","US-HB3",
            "US-A03","US-A10","US-Atq","US-NGB"
)

# Florida sites 

# manually download  warnings
"US-KS4", "US-Skr", "US-Atq"
"US-Hpy","US-Mrm" # not on ameriflux site 
