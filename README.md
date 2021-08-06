# data-structures-algorithms

Codes for Master's Data Structures and Algorithms (FGV, 2021).

## Voronoi image

Voronoi diagram implementation in C++ and an application to images.

Original             |  Voronoi
:-------------------------:|:-------------------------:
![](https://raw.githubusercontent.com/lucasresck/data-structures-algorithms/main/images/cover_1.jpg)  |  ![](https://raw.githubusercontent.com/lucasresck/data-structures-algorithms/main/images/cover_2.png)

### Instructions

First, clone the repository:

```bash
git clone https://github.com/lucasresck/data-structures-algorithms.git
cd data-structures-algorithms/voronoi_image/
```

Run the Python script to sample 10000 pixels from the image at `../images/image_1.jpg`:

```bash
python sample_pixels.py --image_path ../images/image_1.jpg --n_pixels 10000
```

It will generate the file `points.txt`, which contains the sampled pixels from the image.

Compile and execute the C++ file:

```bash
g++ voronoi.cpp -o voronoi
./voronoi
```

It will generate the file `cells.txt`, which contains the cells' vertices from the Voronoi diagram.

Double click `visualization.html` to open the webpage in your browser and upload `points.txt` and `cells.txt` files, in this order.
The result will be the Voronoi diagram constructed from the sampled pixels.

