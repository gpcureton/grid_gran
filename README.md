# Gridding and Granulation using Ctypes

This is a module gridding satellite swath data, and granulating global gridded data. It contains the following files...

```
grid_gran
    ├── griddingAndGranulation.c
    ├── __init__.py
    ├── LandWaterMask.py
    ├── PrecipWater.py
    ├── README.md
    └── snapToGrid.py
```

- `README.md` : This file.
- `__init__.py`: Pythion init file which loads the required modules.
- `snapToGrid.py`: Contains the `SnapToGrid` class, which contains a variety of methods which do gridding with a variety of methods.
- `griddingAndGranulation.c` : Contains various `C` functions for doing gridding and granulation.
- `LandWaterMask.py`: Contains example use of the granulation of a discrete dataset using nearest-neighbour.
- `PrecipWater.py`: Contains example use of the granulation of a continuous dataset using bilinear interpolation.


## Compilation instructions for the `C` dynamic library

- Statically linked

```
gcc -c -fPIC griddingAndGranulation.c
ar rcs libgriddingAndGranulation.a griddingAndGranulation.o
gcc -c -fPIC main.c
gcc -static main.o -L/usr/lib64 -L. -lgriddingAndGranulation -o main_static
```

- Dynamically linked, default optimisation

```
gcc -c -fPIC -O0 griddingAndGranulation.c -o griddingAndGranulation.o
gcc -shared -O0 -Wl,-soname,libgriddingAndGranulation.so.1 -lm -o libgriddingAndGranulation.so.1.0.1 griddingAndGranulation.o
```

- Dynamically linked, more optimisation

```
gcc -c -fPIC -O3 griddingAndGranulation.c -o griddingAndGranulation.o
gcc -shared -O3 -Wl,-soname,libgriddingAndGranulation.so.1 -lm -o libgriddingAndGranulation.so.1.0.1 griddingAndGranulation.o
```
