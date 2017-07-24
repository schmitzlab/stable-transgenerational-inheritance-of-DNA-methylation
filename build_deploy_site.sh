#!/bin/bash

cd _doc
# build site
hugo --theme hugo-material-docs
cd ..
# add changes to git
git add -A

# Commit changes.
msg="rebuilding site `date`"
if [ $# -eq 1 ]
  then msg="$1"
fi
git commit -m "$msg"

# Push source
git push origin gh-pages
