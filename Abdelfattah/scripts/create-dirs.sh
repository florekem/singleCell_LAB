# example: 
# from file names in a directory:
# GSM5518596_rGBM-01-A_barcodes.tsv.gz
#                     _features.tsv.gz
#                     _matrix.mtx.gz
# creates a dir: GSM5518596_rGBM-01-A
# and puts those files inside.
# https://mywiki.wooledge.org/BashFAQ/100
# https://askubuntu.com/questions/618663/creating-a-folder-based-on-a-portion-of-a-file-name

for file in ./*; do
    dir=${file%_*} # remove everything after last _ 
    mkdir -p "./$dir" &&
    mv -iv "$file" "./$dir/${file##*_}"
done

# to print dirnames with newline 
# printf "%s\n" * > file_names.txt