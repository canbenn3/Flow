# Hydrological Flow Simulation
This utility is designed to model heavy rainfall over a period of time and track the movement and aggregation of water over a GIS dataset representing an area's elevation data.

### **1. Prerequisites**

You must install the GDAL development headers and the binary utilities.

Ubuntu / Debian / WSL
```Bash
sudo apt-get update
sudo apt-get install libgdal-dev gdal-bin g++
```
MacOS (using Homebrew)
```Bash
brew install gdal
```
### **2. Compilation**

To compile, you must link the GDAL library using the -lgdal flag.

```Bash
g++ -O3 <filename>.cpp -o <output> -lgdal
```

## Get started with Simple testing files

**1. Compiling the test script**
```Bash
g++ -03 generate_test.cpp -0 generate_test -lgdal
```

**2. Directory Setup**

The script expects a `test_data` folder to exist in the same directory as the executable.

```Bash
mkdir -p test_data
```


**3. Running the Generator script**

Run the executable to produce a few dummy `.tif` files modeling extremely simple geological features:

```Bash
./generate_test 5
```
**4. Verifying the Output**

You can verify that the GeoTIFFs were created correctly and inspect their numerical metadata using the gdalinfo command-line tool:

```Bash
gdalinfo test_data/hill.tif
```
To see the raw floating-point values printed to your terminal:

```Bash
gdalinfo -stats test_data/hill.tif
```
Pro-Tip for CHPC Users
If you are moving this workflow to the University of Utah's Notchpeak or Kingspeak clusters, skip the apt-get step and simply run:

```Bash
module load gdal
g++ -O3 generate_test.cpp -o generate_test -lgdal
```

## Elevation Map generation
**1. Compile code**
```Bash
g++ -O3 parse.cpp -o parse -I/usr/include/gdal -lgdal
```

**2. Directory setup**
```Bash
mkdir output
```

**2. Run code**
```Bash
./parse test_data/bowl.tif
```

You should now see a `.bmp` file inside the `output` directory.