## r-plot

### Introduction

A simple plot library based on `bash`, `R` mainly `ggplot2` to make
the plot easily and flexibly, Developed by Dechun Lin.
This lists the basic information for using [`r-plot`](https://github.com/lindechun/r-plot).

### r-plot Installation Notes

#### Required dependencies

* A UNIX based operating system.

* R v3.4.1 or higher installed.

#### Installation

Download r-plot from GitHub:

```
git clone https://github.com/lindechun/r-plot
```

And then, install required modules.

```
$ cd /your/path/to/r-plot/bin/
$ Rscript dependencies.r
$ R
> webshot::install_phantomjs()
> q()
```
* you will need to add `/your/path/to/r-plot/bin/` to your PATH environment variable.

### Currently supported plots

Please type `r-plot` in command line directly to have the latest list.

```
#### Dots
point3d

#### Histogram
histogram

#### Areaplot
areaplot

#### Pie
pie
```

### Basic test data set

See `/your/path/to/r-plot/data/` for test data set.

### Usage

```r-plot [pie/point3d/histogram/areaplot] [options]```

See `/your/path/to/r-plot/example/` for example of 4 plot classes.

### Basic layouts and themes

* Width, Height, Resolution, type of output pictures

  Default, width and height is automatic set-up, one can give `number` to `-w` to change it, give `number` to `-u` to change height. Give `number` to `-r` to alter resolution instead of using `300` as default.

  8 picture formats are supported, `eps/ps`, `tex` (pictex), `pdf`, `jpeg`, `tiff`, `bmp`, `svg` and `wmf`, with `png` as default. Give any mentioned string to `-E` to change output format.

* Legend position

  Defult, the legend is posited at the right of pictures. One can give `top` to `-P` to put the legend above pictures. Other accepted strings to `-p` is `bottom`,`left`,`right`, or `c(0.008,0.8)`. The two element numerical vactor indicats the reltaive position of legend in pictures. 0.008 means position relative y-axis and 0.8 means position relative to x-axis. Specially, `c(1,1)` put legends at top-right.

* Title, xlab, ylab of picture

  One can set title, xlab, ylab with `-t`, `-x`, `-y`.

* themes for output picture
    One can give theme\_classic2, theme\_classic3, theme\_cin to `-T` to set theme.
