#!/usr/bin/env bash

# uses codedown: https://github.com/earldouglas/codedown

MDFILES=`grep -o "[a-zA-Z0-9_\/]*\.md" mkdocs.yml | tr '\n' ' '`

echo """
####################################
#                                  #
#       Auto-extracted code        #
# DO NOT MANUALLY CHANGE THIS FILE #
#                                  #
####################################
""" > scripts/auto_extracted_code.R

for file in `echo $MDFILES`
do
  printf "\n## Code found in: $file\n" >>  scripts/auto_extracted_code.R
  cat docs/$file | sed 's/^ *//' | codedown R >> scripts/auto_extracted_code.R
done
