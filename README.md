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

# Using OpenTopography to get Digital Elevation Models

We can use [OpenTopography](https://opentopography.org/) to download a digital elevation model (DEM). To get a DEM for a specific area, first click on the DATA tab in the nav bar. This will take you to a new page containing a map showing the different data sources available. For this project, we'll deselect all data sources (found in the top right corner of the map) except for the OpenTopography source.

To fetch the elevation data, zoom in on the target area (e.g. Logan) and click the `SELECT A REGION` button in the top left corner of the map. You can then click and drag a box over the area to download. After selecting the area, there will be some results listed underneath the map. You might see tabs such as `Global & Regional DEM`, `High Resolution Topography`, and `Community Contributed`. To get your data, click on the `High Resolution Topography`, which will display another tab selection allowing you to select which datasource you want to collect. Select `OpenTopography` to see the specific datasets that OpenTopography has associated with your selected area, and click `Get Raster Data` on the far right side of the table row of the dataset to be redirected to the page where you can submit the job to get your data.

On the page you have just been redirected to, scroll down past the map again to select options for the format of your dataset. Make sure that the Data Output format is `GeoTiff`. Finally, you can scroll to the bottom of the page and click `Submit`. This will take you to a new page showing the status of jobs you've submitted (These datasets are in point cloud form, so the server needs to manually process the point cloud data and format it to the `.tif` format).

It may take some time for the job to complete, but when it's finished you can download a compressed file containing your processed `.tif` file.

> Note:
> 
> You'll likely want to create an account on OpenTopography's website, as this will extend the size of data you're allowed to download within a single job. If you associate your account with a University email, it will also grant you access to more precise datasets (e.g. the 1m<sup>2</sup> dataset).
>
> With your DEM download, you can also get an image detailing elevation and other features, such as slope, through the options on the job submission page