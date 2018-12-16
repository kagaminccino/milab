from PIL import Image, ImageDraw, ImageFont
import glob
import re


def extract_time(filename):
    s = re.search(r'\d+', filename).group()
    return int(s) * 0.001


def add_text(filename, time):
    img = Image.open(filename).convert('RGB')
    draw = ImageDraw.Draw(img)
    fnt = ImageFont.truetype(
        '/usr/share/fonts/truetype/freefont/FreeMono.ttf', 40)
    timestamp = "time: %2.3f [s]" % time
    draw.text((10, 10), timestamp, font=fnt, fill=(0, 0, 0))
    return img
    img.save('cap00000.png')


src = glob.glob('*.png')
src.sort()
for i in src:
    time = extract_time(i)
    img = add_text(i, time)
    img.save(i)
# img.show()
