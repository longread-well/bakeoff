args = commandArgs( trailingOnly = TRUE )
input = args[1]
output = args[2]

# read in data in jellyfish dump column format
cat( sprintf( "Reading kmers from %s...\n", input ))
X = read.table( input, head = F, sep = " ", as.is = T )
# convert counts to proportions
cat( sprintf( "%d kmers read.  Converting counts...\n", nrow(X)) )
total = sum( X[,2] )
X[,2] = sprintf( "%.8f", X[,2] / total )

cat( sprintf( "Writing results to %s...\n", output ))
cat( sprintf( "%d\n", nrow(X)), file = output )
write.table( X, file = output, append = T, col.names = F, row.names = F, sep = "\t", quote = F )
cat( sprintf( "Done!  Thanks for using make_mhap_kmer_file.R!\n" ))
