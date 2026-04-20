# Hydrological Flow Simulation
This utility is designed to model heavy rainfall over a period of time and track the movement and aggregation of water over a GIS dataset representing an area's elevation data. It is meant to primarily be run on the University of Utah's CHPC cluster, and all instructions to run the files will be given with that in mind.

## Background

The function we used to calculate water runoff is based off the paper [Two-dimensional kinematic wave model of overland-flow](https://www.sciencedirect.com/science/article/abs/pii/S0022169403005146?via%3Dihub) published by the Journal of Hydrology. The overall concept of the function used is that given a grid of elevation data, we calculate the gradient formed by each group of 4 cells. This results in another set of coordinates (with a width of 1 cell less than the total elevation grid size). Using this grid of gradients, we calculate the amount of water that flows to the cell most aligned with this gradient (for example, if a group of 4 cells has a gradient pointing to 45 degrees from the "x-axis", the water would flow to the upper right cell).

As this methodology doesn't account for the pooling of water, we modified it slightly. Each timestep of our process includes the following:
1. Adding rain water to each cell.
2. Compute the total surface level elevation (including both terrain elevation and water level).
3. Compute the grid geometry (that allows us to predict the flow of water).
4. Calculate where the water flows for the next time step.

The algorithm mentioned above computes grid geometry from terrain. We modified it to calculate the grid geometry from the combination of elevation and water level since pooled water will change where water can flow.

## Output
Each version of our function (serial, OpenMP, etc.) produces a bitmap file that represents the *water level* at each point of the given grid. A space that is completely black represents no water, and the cells with the maximum amount of blue represent an area with over 12 inches of water.

## Running the Code

View [this link](https://usu.box.com/s/03uzjd1dsny32k0cef6m2b8bgbqxw4uz) to get the same data that we've used for our testing.

We use CMake to compile all the code on our cluster. To best match our results, run these commands on the Kingspeak cluster. Run this to load all of the prerequisite models:
```bash
$ module load cuda/11.8.0 openmpi/4.1.6-gpu gdal cmake gcc && export OMPI_MCA_opal_cuda_support=1 && ulimit -l unlimited
```

To compile each implementation, run the following from the project directory:
```bash
$ cmake .
$ make
```

No that all of the binary files are compiled, you can run the scripts defined in `./launch_scripts` to produce the desired output through `$ sbatch ./launch_scripts/<script>`.

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