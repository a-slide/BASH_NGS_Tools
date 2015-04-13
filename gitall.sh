#! /bin/bash

echo -e "Add all modifications\n"
git add --all
git status
echo -e "Commit all modifications\n"
git commit --all -m "$1" 
echo -e "Pushing to Github\n"
git push

exit 0
