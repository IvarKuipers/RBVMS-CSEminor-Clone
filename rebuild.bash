git fetch origin
git pull
cd build
make clean
cmake ..
make -j 16