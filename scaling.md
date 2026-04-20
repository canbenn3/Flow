# Serial VS OpenMP

We used four different files of various sizes to measure how well our serial and OpenMP implementations handle various amounts of data.

## SERIAL IMPLEMENTATION

| FILENAME      | TIME(s)  | NUM ITER |
|---------------|----------|----------|
| Hospital      | 12.0065  | 100      |
| USU           | 110.758  | 25       |
| Hyrum         | 298.107  | 25       |
| cache_valley  | 1430     | 25       |

## OpenMP

All OpenMP studies performed 100 iterations over the map with a rainfall amount of 2.5 inches. 

### WEAK SCALING STUDY- OpenMP

| FILENAME      | FILE SIZE | NUM THREADS | TIME(s)     |
|---------------|-----------|-------------|----------|
| cache_valley  | 612 MB    | 16          | 522.446  |
| Hospital      | 1.8 MB    | 2           | 7.7037   |
| Hyrum         | 174 MB    | 8           | 222.365  |
| USU           | 70.8 MB   | 4           | 151.582  |

---

### STRONG SCALING STUDY- OpenMP

### `cache_valley.tif` (612 MB)

| NUM THREADS | TIME(s) |
|-------------|---------|
| 8           | 867.336 |
| 12          | 786.261 |
| 16          | 540.338 |

### `USU.tif` (70.8 MB)

| NUM THREADS | TIME     |
|-------------|----------|
| 4           | 152.501  |
| 8           | 118.736  |
| 12          | 57.6354  |
| 16          | 58.2229  |

### `Hyrum.tif` (174 MB)

| NUM THREADS | TIME    |
|-------------|---------|
| 4           | 422.956 |
| 8           | 214.55  |
| 12          | 232.105 |
| 16          | 132.658 |

### `Hospital.tif` (1.8 MB)

| NUM THREADS | TIME     |
|-------------|----------|
| 4           | 3.93257  |
| 8           | 2.07409  |
| 12          | 2.26485  |
| 16          | 1.23746  |


## Analysis

As expected, using OpenMP helped us see a vast improvement in the time required to perform the same task on a given file. As serial performance for large files was extremely slow, we capped the larger files at a quarter of the iterations we had OpenMP perform, and still saw faster times than the serial implementation. Here are a few specific speedup times (for serial implementations of 25 iterations, we multiply by 4 to match the time it would have taken to perform the same number of iterations performed with the OpenMP version):

| FILENAME      | SERIAL TIME | NUM THREADS | PARALLEL TIME | SPEEDUP |
|---------------|-------------|-------------|---------------|---------|
| Hospital      | 12.0065     | 2           | 7.7037        | 1.558   |
| USU           | 443.032     | 4           | 151.582       | 2.922   |
| Hyrum         | 1192.428    | 8           | 222.365       | 5.362   |
| cache_valley  | 5720        | 16          | 522.446       | 10.94   |

We never reached the full desired speedup (desired speed up = num threads), but the speedup was still significant. 


---

# SINGLE GPU IMPLEMENTATION

We tested different block sizes for one of our elevation data files:

| FILENAME           | BLOCKSIZE | TIME(s)     |
|--------------------|-----------|-------------|
| Hyrum.tif          | 32        | 69.511      |
|                    | 16        | 67.1297     |
|                    | 64        | 69.429      |
|                    | 8         | 67.3319     |

Results were fairly consistent across different block sizes

---

# Distributed CPU vs. Distributed GPU

## Distributed GPU

### STRONG SCALING STUDY - Distributed GPU

Increasing from 2 to 3 GPU nodes decreased simulation time across all input
sizes except for the smallest file (`Hospital.tif`), which is understandable as
the communication overhead associated with using another node is significant
compared to the rest of the computation at this scale.

### `cache_valley.tif` (612 MB)

| NUM GPU NODES | TIME(s) |
|---------------|---------|
| 2             | 227.647 |
| 3             | 58.690  |
| 4             | 100.765 |

### `USU.tif` (70.8 MB)

| NUM GPU NODES | TIME(s)  |
|---------------|----------|
| 2             | 21.578   |
| 3             | 7.842    |
| 4             | 11.799   |

### `Hyrum.tif` (174 MB)

| NUM GPU NODES | TIME(s) |
|---------------|---------|
| 2             | 57.477  |
| 3             | 16.402  |
| 4             | 26.954  |

### `Hospital.tif` (1.8 MB)

| NUM GPU NODES | TIME(s)  |
|---------------|----------|
| 2             | 2.670    |
| 3             | 4.849    |
| 4             | 3.420    |

---

# VALIDATION

1000 iterations over Hospital with 2.5 inches to `.bmp` file

Compare results using `cmp file1.bmp file2.bmp`
- Serial and OpenMP results are identical

GPU implementations were compared using the ImageMagick compare utility (run
`module load imagemagick` on CHPC):

```
magick compare -metric RMSE file1.bmp file2.bmp NULL:
```

All GPU results reported a normalized error of less than 0.0003, which can be
attributed to floating point error.
