import sys
import glob
from PIL import Image

if __name__ == "__main__":
    if sys.argv > 1:
        fp_in = f"../data/{argv[1]}/frames/*.png"
        fp_out = f"../data/{argv[1]}/gifs/{argv[1]}.gif"
    else:
        fp_in = "../data/debug/frames/*.png"
        fp_out = "../data/debug/gifs/image.gif"
    print(f"Reading Frames from: {fp_in}.")
    print(f"Writing GIF to: {fp_out}.")

    prefix_len = len(fp_in[:-5])
    img, *imgs = [Image.open(f) for f in sorted(glob.glob(fp_in), key=lambda name: int(name[prefix_len:-4]))]
    img.save(fp=fp_out, format='GIF', append_images=imgs, save_all=True, duration=100, loop=0)