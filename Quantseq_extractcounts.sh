#!/bin/bash

# Create a file with Gene IDs
awk 'NR>4 {print $1}' $(ls *ReadsPerGene.out.tab | head -1) > Gene_IDs.txt

# Initialize the header with "Gene_ID"
echo -n "Gene_ID" > header.txt

# Start the merged counts table with gene IDs
cp Gene_IDs.txt Merged_CountsTable.txt

# Extract counts and append to the merged table
for file in *ReadsPerGene.out.tab; do
    # Extract sample/library name from filename
    sample_name=$(basename "$file" "_ReadsPerGene.out.tab")
    
    # Add sample name to header
    echo -ne "\t$sample_name" >> header.txt
    
    # Extract counts from column 3
    awk 'NR>4 {print $3}' $file > temp_counts.txt
    
    # Merge with the existing table
    paste Merged_CountsTable.txt temp_counts.txt > temp.txt && mv temp.txt Merged_CountsTable.txt
done

# Insert a newline at the end of the header
echo "" >> header.txt

# Prepend header to the merged table
cat header.txt Merged_CountsTable.txt > temp.txt && mv temp.txt Merged_CountsTable.txt

# Clean up temporary files
rm temp_counts.txt Gene_IDs.txt header.txt

echo "Merged table created: Merged_CountsTable.txt"