for file in Data/chr1_388*.txt; do
    cut -f 1,2,3,4,5,6,7,8,9,10 "$file" > "${file%.txt}_short.txt"
done 
