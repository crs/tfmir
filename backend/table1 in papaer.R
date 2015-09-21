cat="mirna-gene"

unique(dbs[dbs$category ==cat ,"source"])

df=dbs[dbs$category ==cat & dbs$source %in% c("Starbase") , ]

df=dbs
print( "## no of edges>>>"  )
dim(df)
print("==================")

nodes=unique( c(as.character(df$regulator), as.character(df$target)  ))
mirnas=nodes[grep("hsa-",nodes)]
print( " no of mirnas is >>>>>")
print(length(unique(mirnas)))
print("==================")

genes=setdiff(nodes,mirnas)
print( " no of genes is >>>>>")
length(unique(genes))



