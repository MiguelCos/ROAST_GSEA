# Function to get the normalized data from MSstats/MSstatsTMT; format it ----
# to ROAST and generate the design file 

MSstatsToRoast <- function(data,
                           MSstatsAnnotation,
                           MSstatsType, # "MSstats" or "MSsstatsTMT"
                           Paired = FALSE, # if paired design, set to TRUE
                           Conditions = c("Control", "Treatment") # a character vector of two elements corresponding to the same experimental conditions of your setting in MSstats 
                           ) {
      
} 