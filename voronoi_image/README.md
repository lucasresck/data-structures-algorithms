# Voronoi image

Voronoi diagram implementation in C++ and an application to images.

Original image             |  Voronoi image
:-------------------------:|:-------------------------:
![](https://raw.githubusercontent.com/lucasresck/data-structures-algorithms/main/images/cover_1.jpg)  |  ![](https://raw.githubusercontent.com/lucasresck/data-structures-algorithms/main/images/cover_2.png)

This was an evaluated assignment for the first part of this course, that covered Computational Geometry problems.

## Usage

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

## Details

The Voronoi diagram implementation is Fortune's, based on the book [Computational Geometry](https://link.springer.com/book/10.1007/978-3-662-03427-9).
The C++ file receives a list of n points, constructs a doubly-connected edge list (DCEL) that represents the Voronoi diagram, and saves the cells' vertices as a file.

The algorithm, as original designed, having the status structure as a balanced binary search tree, has complexity O(n log n). However, due to time constraints, I implemented a linked list.
This should not be a problem, given the fact the main structure was the DCEL and the complexity went to O(n²).
Another limitation was the bounding box, which is not (at least yet) present in this implementation.
The main effect is the empty cells on the image border.

## References

De Berg, Mark, et al. "Computational geometry." Computational geometry. Springer, Berlin, Heidelberg, 1997.

[Matt Brubeck](https://github.com/mbrubeck/)'s "[Fortune's Algorithm in C++](https://www.cs.hmc.edu/~mbrubeck/voronoi.html)".
