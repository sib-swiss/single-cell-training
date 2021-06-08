#!/usr/bin/env bash

# uses codedown: https://github.com/earldouglas/codedown

cd ~/Documents/repositories/single-cell-training

MDFILES=`grep -o "[a-zA-Z0-9_\/]*\.md" mkdocs.yml | tr '\n' ' '`

echo "# Auto-extracted code from markdown files" > scripts/auto_extracted_code.R

for file in `echo $MDFILES`
do
  printf "\n## Code found in: $file\n" >>  scripts/auto_extracted_code.R
  cat docs/$file | codedown R >> scripts/auto_extracted_code.R
done
