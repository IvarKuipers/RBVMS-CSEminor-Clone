# This file is part of the RBVMS application. For more information and source code
# availability visit https://idoakkerman.github.io/
#
# RBVMS is free software; you can redistribute it and/or modify it under the
# terms of the BSD-3 license.

# Get build info
PWD=`pwd`
TEMP=$PWD/.tmp.buildInfo.hpp


if [ $# -le 0 ]; then
  echo "Illegal number of inputs"
  echo "Provide at least the RBVMS source directory"
  exit -1
fi

# RBVMS directory
RBVMS=$1
shift

# Check if other repos need to be included
MFEM=false
HYPRE=false
for var in $@; do
  if [ ${var^^} == HYPRE ]; then
     HYPRE=true
  fi
  if [ ${var^^} == MFEM ]; then
     MFEM=true
  fi
done

# Info per repo
repoInfo () {
  echo -n "    dir     "| tee -a $TEMP;pwd | tee -a $TEMP
  echo -n "    repo    "| tee -a $TEMP;git remote get-url origin | tee -a $TEMP
  echo -n "    branch  "| tee -a $TEMP;git rev-parse --abbrev-ref HEAD | tee -a $TEMP
  echo -n "    commit  "| tee -a $TEMP;git rev-parse HEAD | tee -a $TEMP
  echo -n "    changes "| tee -a $TEMP;git status -sb | tee -a $TEMP
}

# Create new hpp file from git info
echo "// Autogenerated using buildInfo.sh" >  $TEMP
echo "#include <sstream>"                  >> $TEMP
echo 'std::istringstream buildInfo(R"('    >> $TEMP
echo "------------------------------------"   | tee -a $TEMP
echo "  Build info:"                          | tee -a $TEMP
echo -n "    host    "| tee -a $TEMP;hostname | tee -a $TEMP
echo -n "    date    "| tee -a $TEMP;date +%F | tee -a $TEMP

echo "  RBVMS git info:"                      | tee -a $TEMP
cd $RBVMS
repoInfo

if $HYPRE; then
  echo "  HYPRE git info:"| tee -a $TEMP
  cd $RBVMS/../hypre
  repoInfo
fi

if $MFEM; then
  echo "  MFEM git info:"| tee -a $TEMP
  cd $RBVMS/../mfem
  repoInfo
fi

echo "------------------------------------"| tee -a $TEMP
echo ')");'>>$TEMP

# Replace hpp if necessary
cd $RBVMS
if test -f buildInfo.hpp; then
  cmp -s $TEMP buildInfo.hpp || mv $TEMP buildInfo.hpp
else
  mv $TEMP buildInfo.hpp
fi
ls -ltrh buildInfo.hpp

# Cleanup tmp file
rm -f $TEMP



