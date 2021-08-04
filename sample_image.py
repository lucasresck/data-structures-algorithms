from PIL import Image
import argparse
import numpy as np

parser = argparse.ArgumentParser(description='Read an image.')
parser.add_argument('--image_path', default='image.jpg', type=str, help='Path of the image.')
parser.add_argument('--n_pixels', default=100, type=int, help='Number of random pixels.')
args = parser.parse_args()

im = Image.open(args.image_path)
pix = im.load()

np.random.seed(42)
random_pixels = set()
while len(random_pixels) < args.n_pixels:
    x = np.random.randint(0, im.size[0])
    y = np.random.randint(0, im.size[1])
    color = pix[x, y]
    random_pixels.add((color, x, y))

with open('points.txt', 'w') as f:
    for pixel in random_pixels:
        f.write(str(pixel[0][0]) + ' ' + str(pixel[0][1]) + ' ' + str(pixel[0][2]) + ' ')
        f.write(str(pixel[1]) + ' ' + str(pixel[2]) + '\n')
