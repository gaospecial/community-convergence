# create local KEGG db
library(createKEGGdb)
createKEGGdb::create_kegg_db(species = c("eco","ppu"))
install.packages("KEGG.db_1.0.tar.gz", type = "source")
