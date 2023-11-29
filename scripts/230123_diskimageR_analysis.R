#load diskImageR
library(diskImageR)

#Run the ImageJ analysis component, save the output. "newProject" should be changed to something of your choice (and then the same name used throughout); note that the quotation marks are required.
#To use a pop-up box interface:
IJMacro("fortyeighthourdataset")

#Plot the result of ImageJ analysis (averaged among 72 lines draft outward from the center of the diffusion disk). Type ?plotRaw for additional parameter options.
plotRaw("fortyeighthourdataset")

#Use maximum likelihood to fit a bilogistic and single logistic model to the data from each photograph. "clearHalo" is used to specify a picture that has a clear halo; this is used to standardize all photographs and will be most effective when photographs are taken with equal lighting without shadows.
maxLik("fortyeighthourdataset", clearHalo=1, RAD="all")

createDataframe("fortyeighthourdataset", clearHalo = 1)

addType("fortyeighthourdataset", typeName="rep", typePlace = 2)

addType("fortyeighthourdataset", typePlace = 3, typeName="time")

library(readr)
library(here)

#strip off the "R", we want the replicate to be numeric
fortyeighthourdataset_df<-read.csv(here("parameter_files", "fortyeighthourdataset_df.csv")) 
fortyeighthourdataset_df$rep <- sub("R", " ", fortyeighthourdataset_df$rep)

#strip off the replicate letter, just keep replicate number
fortyeighthourdataset_df$rep <- sub('A', '', fortyeighthourdataset_df$rep)
fortyeighthourdataset_df$rep <- sub('B', '', fortyeighthourdataset_df$rep)
fortyeighthourdataset_df$rep <- sub('C', '', fortyeighthourdataset_df$rep)

library(here)
write.csv(fortyeighthourdataset_df, here("fortyeighthourdataset_df_edited.csv"))
